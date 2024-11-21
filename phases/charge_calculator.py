import argparse
from math import dist
from os import path, system

import gemmi
from Bio import PDB
from rdkit import Chem



class AtomSelector(PDB.Select):
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)



def load_arguments():
    print("\nParsing arguments... ", end="")
    parser = argparse.ArgumentParser()
    parser.add_argument("--PDB_file",
                        help="PDB file with protein structure.",
                        type=str,
                        required=True)
    parser.add_argument("--data_dir",
                        help="Directory for saving results.",
                        type=str,
                        required=True)
    parser.add_argument("--delete_auxiliary_files",
                        help="Auxiliary calculation files can be large. With this argument, "
                             "the auxiliary files will be continuously deleted during the calculation.",
                        action="store_true")
    args = parser.parse_args()
    if not path.isfile(args.PDB_file):
        exit(f"\nERROR! File {args.PDB_file} does not exist!\n")
    if path.exists(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument!\n")
    print("ok")
    #TODO argument pro estimation, kontrola, zda existuje
    return args


class ChargeCalculator:
    """
    TODO
    """
    def __init__(self,
                 PDB_file: str,
                 charges_estimation: str,
                 data_dir: str,
                 delete_auxiliary_files: bool):
        """
        :param PDB_file: PDB file containing the structure for which the partial atomic charges are to be calculated
        :param charges_estimation TODO
        :param data_dir: directory where the results will be stored
        :param delete_auxiliary_files: auxiliary files created by the calculation taking up a significant amount of space will be deleted
        """
        self.PDB_file = PDB_file
        self.charges_estimation = charges_estimation
        self.data_dir = data_dir
        self.delete_auxiliary_files = delete_auxiliary_files
        system(f"mkdir {self.data_dir}")

    def calculate_charges(self):
        print("Calculating of charges... ", end="")

        # load structure by Biopython
        structure = PDB.PDBParser(QUIET=True).get_structure("structure", self.PDB_file)[0]
        structure_atoms = list(structure.get_atoms())
        structure_num_of_atoms = len(structure_atoms)
        selector = AtomSelector()
        io = PDB.PDBIO()
        io.set_structure(structure)
        kdtree = PDB.NeighborSearch(structure_atoms)

        # load partial atomic charges estimation
        if self.charges_estimation:
            with open(self.charges_estimation, "r") as estimated_charges_file:
                charge_estimations = [float(x) for x in estimated_charges_file.readlines()[1].split()]
        else:
            charge_estimations = [0 for _ in range(structure_num_of_atoms)]

        # creating charge attributes to make them easy to work with in Biopython library
        for atom, atomic_charge_estimation in zip(structure_atoms, charge_estimations):
            atom.charge_estimation = atomic_charge_estimation
            atom.cm5_charge = None

        # definition of radii limiting the substructure
        # all atoms that are closer to the calculated atom than min_radius are included in the substructure
        # atoms more distant from the calculated atom than max_radius are never included in the substructure
        min_radius = 6
        max_radius = 12

        # calculate the charges for each atom using the cutoff approach.
        for calculated_atom_i, calculated_atom in enumerate(structure_atoms):

            print(calculated_atom_i)

            # To speed up the calculation, the charges of the hydrogen and oxygen atoms bound to one atom
            # are calculated together with the nearest other heavy atoms
            if calculated_atom.element == "H":
                continue
            elif calculated_atom.element == "O":
                if len(kdtree.search(calculated_atom.coord, 1.5, level="A")) <= 2:
                    continue

            substructure_data_dir = f"{self.data_dir}/sub_{calculated_atom_i}"
            system(f"mkdir {substructure_data_dir}")

            # save substructures by biopython
            atoms_in_min_radius = kdtree.search(calculated_atom.coord, min_radius, level="A")
            selector.full_ids = set([atom.full_id for atom in atoms_in_min_radius])
            io.save(f"{substructure_data_dir}/atoms_in_{min_radius}_angstroms.pdb", selector)
            atoms_in_max_radius = kdtree.search(calculated_atom.coord, max_radius, level="A")
            selector.full_ids = set([atom.full_id for atom in atoms_in_max_radius])
            io.save(f"{substructure_data_dir}/atoms_in_{max_radius}_angstroms.pdb", selector)

            # load substructures by RDKit to determine bonds
            mol_min_radius = Chem.MolFromPDBFile(f"{substructure_data_dir}/atoms_in_{min_radius}_angstroms.pdb", removeHs=False, sanitize=False)
            mol_min_radius_conformer = mol_min_radius.GetConformer()
            mol_max_radius = Chem.MolFromPDBFile(f"{substructure_data_dir}/atoms_in_{max_radius}_angstroms.pdb", removeHs=False, sanitize=False)
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
                            carbons_with_broken_bonds_coord.append(mol_max_radius_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                            continue
                        else:
                            atoms_with_broken_bonds.append(bonded_atom)
                            substructure_coord_dict[(coord.x, coord.y, coord.z)] = bonded_atom

            # create substructure in Biopython library
            substructure_atoms = []
            for coord in substructure_coord_dict.keys():
                substructure_atoms.append(kdtree.search(coord, 0.1, level="A")[0])
            substructure_atoms.sort(key=lambda x: x.serial_number)
            selector.full_ids = set([atom.full_id for atom in substructure_atoms])
            io.save(f"{substructure_data_dir}/substructure.pdb", selector)

            # add hydrogens to broken C-C bonds by openbabel
            system(f"cd {substructure_data_dir} ; obabel -iPDB -oPDB substructure.pdb -h > readded_hydrogens_substructure.pdb 2>/dev/null")
            with open(f"{substructure_data_dir}/readded_hydrogens_substructure.pdb") as readded_hydrogens_substructure_file:
                atom_lines = [line for line in readded_hydrogens_substructure_file.readlines() if line[:4] in ["ATOM", "HETA"]]
                original_atoms_lines = atom_lines[:len(substructure_atoms)]
                added_hydrogens_lines = atom_lines[len(substructure_atoms):]
            with open(f"{substructure_data_dir}/repaired_substructure.pdb", "w") as repaired_substructure_file:
                repaired_substructure_file.write("".join(original_atoms_lines))
                for added_hydrogen_line in added_hydrogens_lines:
                    added_hydrogen_coord = (float(added_hydrogen_line[30:38]), float(added_hydrogen_line[38:46]), float(added_hydrogen_line[46:54]))
                    if any([dist(added_hydrogen_coord, carbon_coord) < 1.3 for carbon_coord in carbons_with_broken_bonds_coord]):
                        repaired_substructure_file.write(added_hydrogen_line)

            # calculate charges for substructure
            substructure_charge = round(sum([atom.charge_estimation for atom in substructure_atoms]))
            system(f"cd {substructure_data_dir} ; "
                   f"xtb repaired_substructure.pdb --gfn 1 --gbsa water --acc 1000 --chrg {substructure_charge}   > xtb_output.txt 2> xtb_error_output.txt ")

            # read calculated charges from xtb output file
            xtb_output_file_lines = open(f"{substructure_data_dir}/xtb_output.txt").readlines()
            try:
                cm5_charges_headline = "  Mulliken/CM5 charges         n(s)   n(p)   n(d)\n"
                charge_headline_index = xtb_output_file_lines.index(cm5_charges_headline)
            except ValueError:  # charge calculation failed
                continue

            # define for which atoms we have calculated the charges
            calculated_atoms = set([calculated_atom])
            for near_atom in kdtree.search(calculated_atom.coord, 1.5, level="A"):
                if near_atom.element == "H":
                    calculated_atoms.add(near_atom)
                elif near_atom.element == "O":
                    if len(kdtree.search(near_atom.coord, 1.5, level="A")) <= 2:
                        calculated_atoms.add(near_atom)
            calculated_atoms_full_ids = set([calculated_atom.full_id for calculated_atom in calculated_atoms])

            # write the charges into the Biopython structure
            for substructure_index, substructure_atom in enumerate(substructure_atoms):
                if substructure_atom.full_id in calculated_atoms_full_ids:
                    charge = float(xtb_output_file_lines[charge_headline_index + substructure_index + 1].split()[3])
                    substructure_atom.cm5_charge = charge

            if self.delete_auxiliary_files:
                system(f"rm -r {substructure_data_dir}")

        self.cm5_charges = [atom.cm5_charge for atom in structure_atoms]
        print("ok")

    def write_charges_to_files(self):
        print("Writing charges to file... ", end="")
        with open(f"{self.data_dir}/charges.txt", "w") as charges_file:
            charges_string = f"{path.basename(self.PDB_file)[:-4]}\n" + " ".join([str(x) for x in self.cm5_charges])
            charges_file.write(charges_string)

        exit()
        # write charges to mmCIF file
        structure = gemmi.cif.read_file(self.mmCIF_file)
        block = structure.sole_block()
        block.find_mmcif_category('_chem_comp.').erase()
        sb_ncbr_partial_atomic_charges_meta_prefix = "_sb_ncbr_partial_atomic_charges_meta."
        sb_ncbr_partial_atomic_charges_meta_attributes = ["id",
                                                          "type",
                                                          "method"]
        metadata_loop = block.init_loop(sb_ncbr_partial_atomic_charges_meta_prefix,
                                        sb_ncbr_partial_atomic_charges_meta_attributes)
        metadata_loop.add_row(['1',
                               "'QM'",
                               "'GFN1-xTB/CM5 (with cutoff)'"])
        sb_ncbr_partial_atomic_charges_prefix = "_sb_ncbr_partial_atomic_charges."
        sb_ncbr_partial_atomic_charges_attributes = ["type_id",
                                                     "atom_id",
                                                     "charge"]
        charges_loop = block.init_loop(sb_ncbr_partial_atomic_charges_prefix,
                                       sb_ncbr_partial_atomic_charges_attributes)
        for atomId, charge in enumerate(self.cm5_charges):
            if isinstance(charge, float):
                charge = f"{charge: .4f}"
            else:
                charge = "?"
            charges_loop.add_row(["1",
                                  f"{atomId + 1}",
                                  f"{charge}"])
        block.write_file(f"{self.data_dir}/with_charges.cif")
        print("ok")


if __name__ == "__main__":
    args = load_arguments()
    calculator = ChargeCalculator(PDB_file=args.PDB_file,
                                  data_dir=args.data_dir,
                                  delete_auxiliary_files=args.delete_auxiliary_files)
    calculator.calculate_charges()
    calculator.write_charges_to_files()
