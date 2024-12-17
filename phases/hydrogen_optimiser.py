from Bio.PDB import Select, PDBIO, MMCIFIO, MMCIFParser, Superimposer, NeighborSearch, PDBParser
from os import system, path
from math import dist
from rdkit import Chem
import tqdm


class AtomSelector(Select):
    """
    Support class for Biopython.
    After initialization, a set with all full ids of the atoms to be written into the substructure must be stored in self.full_ids.
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)


class HydrogenOptimiser:
    """
    This class optimises hydrogen positions. It uses GFN-FF force field method
    from software xtb to optimise hydrogens. Cutoff approach is employed to achieve more
    calculation speed.
    """
    def __init__(self,
                 input_mmCIF_file: str,
                 logger,
                 output_mmCIF_file: str,
                 data_dir: str,
                 delete_auxiliary_files: bool):
        """
        :param input_mmCIF_file: PDB file containing the structure which should be prepared
        :param logger: loger of workflow to unify outputs
        :param data_dir: directory where the results will be stored
        :param output_mmCIF_file: mmCIF file in which prepared structure will be stored
        :param delete_auxiliary_files: auxiliary files created during the preraparation will be deleted
        """
        self.logger = logger
        self.logger.print("\nHYDROGEN OPTIMISER")
        self.logger.print("Hydrogen optimiser initialization... ", end="")
        self.input_mmCIF_file = input_mmCIF_file
        self.output_mmCIF_file = output_mmCIF_file
        self.data_dir = data_dir
        system(f"mkdir {self.data_dir}")
        self.delete_auxiliary_files = delete_auxiliary_files
        self.logger.print("ok")

    def optimise(self):

        # load structure by Biopython
        self.logger.print("Loading structure... ", end="")
        structure = MMCIFParser(QUIET=True).get_structure(structure_id="structure",
                                                          filename=self.input_mmCIF_file)
        io = PDBIO()
        io.set_structure(structure)
        self.io = io
        self.structure = io.structure[0]
        self.selector = AtomSelector()
        self.logger.print("ok")

        self.logger.print("Optimisation of hydrogens... ", end="", silence=True)
        for atom in tqdm.tqdm(list(self.structure.get_atoms()),
                              desc="Hydrogen optimisation",
                              unit="atoms",
                              smoothing=0,
                              delay=0.1,
                              mininterval=0.4,
                              maxinterval=0.4):
            if atom.element == "H":
                continue
            self.optimise_atom(atom)

        # write logs
        for residue in self.structure.get_residues():
            non_optimised_hydrogens = []
            for atom in residue.get_atoms():
                try:
                    if not atom.optimised:
                        non_optimised_hydrogens.append(atom.name)
                except AttributeError:
                    continue
            if non_optimised_hydrogens:
                warning = f"Optimisation of hydrogen(s) {' '.join(non_optimised_hydrogens)} failed."
                self.logger.add_warning(chain=residue.get_parent().id,
                                        resnum=residue.id[1],
                                        resname=residue.resname,
                                        warning=warning)
        self.logger.print("ok", silence=True)

        self.logger.print("Writing structure with optimised hydrogens to file... ", end="")
        self.io = MMCIFIO()
        self.io.set_structure(self.structure)
        self.io.save(f"{self.data_dir}/{self.output_mmCIF_file}")
        if self.delete_auxiliary_files:
            system(f"for au_file in {self.data_dir}/sub_* ; do rm -fr $au_file ; done &")
        self.logger.print("ok")

    def optimise_atom(self,
                      central_atom):

        # creation of substructure
        self.kdtree = NeighborSearch(list(self.structure.get_atoms()))
        substructure_data_dir = f"{self.data_dir}/sub_{central_atom.serial_number}"
        system(f"mkdir {substructure_data_dir}")
        bonded_hydrogens = [atom for atom in self.kdtree.search(center=central_atom.coord,
                                                                radius=3,
                                                                level="A") if atom.element == "H"]
        if not bonded_hydrogens:
            return
        bonded_hydrogens_full_ids = (set(atom.full_id for atom in bonded_hydrogens))

        # create and save min_radius and max_radius substructures by biopython
        atoms_in_6A = self.kdtree.search(center=central_atom.coord,
                                         radius=6,
                                         level="A")
        atoms_in_12A = self.kdtree.search(center=central_atom.coord,
                                          radius=12,
                                          level="A")
        self.selector.full_ids = set([atom.full_id for atom in atoms_in_6A])
        self.io.save(file=f"{substructure_data_dir}/atoms_in_6A.pdb",
                     select=self.selector,
                     preserve_atom_numbering=True)
        self.selector.full_ids = set([atom.full_id for atom in atoms_in_12A])
        self.io.save(file=f"{substructure_data_dir}/atoms_in_12A.pdb",
                     select=self.selector,
                     preserve_atom_numbering=True)

        # load substructures by RDKit to determine bonds
        mol_min_radius = Chem.MolFromPDBFile(molFileName=f"{substructure_data_dir}/atoms_in_6A.pdb",
                                             removeHs=False,
                                             sanitize=False)
        mol_min_radius_conformer = mol_min_radius.GetConformer()
        mol_max_radius = Chem.MolFromPDBFile(molFileName=f"{substructure_data_dir}/atoms_in_12A.pdb",
                                             removeHs=False,
                                             sanitize=False)
        mol_max_radius_conformer = mol_max_radius.GetConformer()

        # dictionaries allow quick and precise matching of atoms from mol_min_radius and mol_max_radius
        mol_min_radius_coord_dict = {}
        for i, mol_min_radius_atom in enumerate(mol_min_radius.GetAtoms()):
            coord = mol_min_radius_conformer.GetAtomPosition(i)
            mol_min_radius_coord_dict[(coord.x, coord.y, coord.z)] = mol_min_radius_atom
        mol_max_radius_coord_dict = {}
        for i, mol_max_radius_atom in enumerate(mol_max_radius.GetAtoms()):
            coord = mol_max_radius_conformer.GetAtomPosition(i)
            mol_max_radius_coord_dict[(coord.x, coord.y, coord.z)] = mol_max_radius_atom

        # find atoms from mol_min_radius with broken bonds
        atoms_with_broken_bonds = []
        for mol_min_radius_atom in mol_min_radius.GetAtoms():
            coord = mol_min_radius_conformer.GetAtomPosition(mol_min_radius_atom.GetIdx())
            mol_max_radius_atom = mol_max_radius_coord_dict[(coord.x, coord.y, coord.z)]
            if len(mol_min_radius_atom.GetNeighbors()) != len(mol_max_radius_atom.GetNeighbors()):
                atoms_with_broken_bonds.append(mol_max_radius_atom)

        # create a substructure that will have only C-C bonds broken
        carbons_with_broken_bonds_coord = []  # hydrogens will be added only to these carbons
        substructure_coord_dict = mol_min_radius_coord_dict
        while atoms_with_broken_bonds:
            atom_with_broken_bonds = atoms_with_broken_bonds.pop(0)
            bonded_atoms = atom_with_broken_bonds.GetNeighbors()
            for bonded_atom in bonded_atoms:
                coord = mol_max_radius_conformer.GetAtomPosition(bonded_atom.GetIdx())
                if (coord.x, coord.y, coord.z) in substructure_coord_dict:
                    continue
                else:
                    if atom_with_broken_bonds.GetSymbol() == "C" and bonded_atom.GetSymbol() == "C":
                        carbons_with_broken_bonds_coord.append(
                            mol_max_radius_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                        continue
                    else:
                        atoms_with_broken_bonds.append(bonded_atom)
                        substructure_coord_dict[(coord.x, coord.y, coord.z)] = bonded_atom

        # create substructure in Biopython library
        # we prefer to use kdtree because serial_id may be discontinuous in some pdbs files
        # for example, a structure with PDB code 107d and its serial numbers 218 and 445
        substructure_atoms = [self.kdtree.search(center=coord,
                                                 radius=0.1,
                                                 level="A")[0] for coord in substructure_coord_dict.keys()]
        self.selector.full_ids = set([atom.full_id for atom in substructure_atoms])
        self.io.save(file=f"{substructure_data_dir}/substructure.pdb",
                     select=self.selector,
                     preserve_atom_numbering=True)
        substructure = PDBParser(QUIET=True).get_structure(id="structure",
                                                           file=f"{substructure_data_dir}/substructure.pdb")[0]

        # define constrained atoms
        constrained_atom_indices = []
        for atom_index, atom in enumerate(substructure.get_atoms(), start=1):
            if atom.full_id in bonded_hydrogens_full_ids:
                continue
            constrained_atom_indices.append(str(atom_index))

        # add hydrogens to broken C-C bonds by openbabel
        system(f"cd {substructure_data_dir} ; obabel -iPDB -oPDB substructure.pdb -h > readded_hydrogens_substructure.pdb 2>/dev/null")
        with open(f"{substructure_data_dir}/readded_hydrogens_substructure.pdb") as readded_hydrogens_substructure_file:
            atom_lines = [line for line in readded_hydrogens_substructure_file.readlines() if line[:4] in ["ATOM", "HETA"]]
            original_atoms_lines = atom_lines[:len(substructure_atoms)]
            added_hydrogens_lines = atom_lines[len(substructure_atoms):]
        with open(f"{substructure_data_dir}/repaired_substructure.pdb", "w") as repaired_substructure_file:
            added_hydrogen_indices = []
            added_hydrogen_indices_counter = len(substructure_atoms) + 1
            repaired_substructure_file.write("".join(original_atoms_lines))
            for added_hydrogen_line in added_hydrogens_lines:
                added_hydrogen_coord = (float(added_hydrogen_line[30:38]),
                                        float(added_hydrogen_line[38:46]),
                                        float(added_hydrogen_line[46:54]))
                if any([dist(added_hydrogen_coord, carbon_coord) < 1.3 for carbon_coord in carbons_with_broken_bonds_coord]):
                    repaired_substructure_file.write(added_hydrogen_line)
                    added_hydrogen_indices.append(str(added_hydrogen_indices_counter)) # added hydrogens should be also constrained
                    added_hydrogen_indices_counter += 1

        # optimise substructure by xtb
        xtb_settings_template = """$constrain
        atoms: xxx
        force constant=10
        $end
        $opt
        engine=rf
        $end
        """
        substructure_settings = xtb_settings_template.replace("xxx", ", ".join(constrained_atom_indices + added_hydrogen_indices))
        with open(f"{substructure_data_dir}/xtb_settings.inp", "w") as xtb_settings_file:
            xtb_settings_file.write(substructure_settings)
        run_xtb = (f"cd {substructure_data_dir} ;"
                   f"ulimit -s unlimited ;"
                   f"export OMP_NUM_THREADS=1,1 ;"
                   f"export OMP_MAX_ACTIVE_LEVELS=1 ;"
                   f"export MKL_NUM_THREADS=1 ;"
                   f"xtb repaired_substructure.pdb --gfnff --input xtb_settings.inp --opt --gbsa water --verbose > xtb_output.txt 2>&1")
        # second try by L-ANCOPT
        if not path.isfile(f"{substructure_data_dir}/xtbopt.pdb"):
            substructure_settings = open(f"{substructure_data_dir}/xtb_settings.inp", "r").read().replace("rf","lbfgs")
            with open(f"{substructure_data_dir}/xtb_settings.inp", "w") as xtb_settings_file:
                xtb_settings_file.write(substructure_settings)
            system(run_xtb)

        if path.isfile(f"{substructure_data_dir}/xtbopt.pdb"):
            optimised_substructure = PDBParser(QUIET=True).get_structure(id="structure",
                                                                         file=f"{substructure_data_dir}/xtbopt.pdb")[0]
            original_constrained_atoms = []
            optimised_constrained_atoms = []
            for atom in substructure.get_atoms():
                if atom.full_id in bonded_hydrogens_full_ids:
                    continue
                original_constrained_atoms.append(atom)
                optimised_constrained_atoms.append(optimised_substructure[atom.get_parent().get_parent().id][atom.get_parent().id][atom.id])
            sup = Superimposer()
            sup.set_atoms(original_constrained_atoms, optimised_constrained_atoms)
            sup.apply(optimised_substructure.get_atoms())
            for atom in optimised_substructure.get_atoms():
                if atom.full_id in bonded_hydrogens_full_ids:
                    self.structure[atom.get_parent().get_parent().id][atom.get_parent().id][atom.name].coord = atom.coord
        else:
            for hydrogen in bonded_hydrogens:
                hydrogen.optimised = False
