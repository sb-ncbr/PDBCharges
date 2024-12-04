import argparse
from math import dist
from os import path, system

import gemmi
import tqdm
from Bio import PDB
from rdkit import Chem


class AtomSelector(PDB.Select):
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)


def load_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_mmCIF_file",
                        help="mmCIF file with protein structure.",
                        type=str,
                        required=True)
    parser.add_argument("--output_mmCIF_file",
                        help="mmCIF file to store structure with charges. "
                             "mmCIF file will be save into directory defined by --data_dir argument.",
                        type=str,
                        required=True)
    parser.add_argument("--data_dir",
                        help="Directory for saving results.",
                        type=str,
                        required=True)
    parser.add_argument("--charges_estimation",
                        help="File with estimation of partial atomic charges.",
                        type=str,
                        required=False)
    parser.add_argument("--delete_auxiliary_files",
                        help="Auxiliary calculation files can be large. With this argument, "
                             "the auxiliary files will be continuously deleted during the calculation.",
                        action="store_true")
    return parser.parse_args()


class ChargeCalculator:
    """
    This class calculated partial atomic charges for proteins. Specifically, it uses GFN1 semiempirical QM method
    from software xtb to calculate reproduction of PBE0/TZVP/CM5 charges. Cutoff approach is employed to achieve more
    calculation speed.

    Calculated charges are stored directly into mmCIF file and also in to txt file
    in a user-defined data directory in mmCIF format.
    """
    def __init__(self,
                 input_mmCIF_file: str,
                 charges_estimation: str,
                 output_mmCIF_file: str,
                 data_dir: str,
                 delete_auxiliary_files: bool):
        """
        :param input_mmCIF_file: mmCIF file containing the structure for which the partial atomic charges are to be calculated
        :param charges_estimation: txt file with estimation of partial atomic charges for more accurate results
        :param output_mmCIF_file: mmCIF file in which calculated partial atomic charges will be stored
        :param data_dir: directory where the results will be stored
        :param delete_auxiliary_files: auxiliary files created by the calculation taking up a significant amount of space will be deleted
        """

        print("\nCHARGE CALCULATOR")
        print("Charge calculator initialization... ", end="")

        if not path.isfile(input_mmCIF_file):
            exit(f"\nERROR! File {input_mmCIF_file} does not exist!\n")
        if path.exists(data_dir):
            exit(f"\nError! Directory with name {data_dir} exists. "
                 f"Remove existed directory or change --data_dir argument!\n")
        if charges_estimation and not path.isfile(charges_estimation):
            exit(f"\nERROR! File {charges_estimation} does not exist!\n")

        self.output_mmCIF_file = output_mmCIF_file
        self.charges_estimation = charges_estimation
        self.data_dir = data_dir
        self.delete_auxiliary_files = delete_auxiliary_files

        system(f"mkdir {self.data_dir}")
        system(f"cp {input_mmCIF_file} {self.data_dir}/{self.output_mmCIF_file}")

        print("ok")


    def calculate_charges(self):

        # load structure by Biopython
        print("Loading structure... ", end="")
        structure = PDB.MMCIFParser(QUIET=True).get_structure(structure_id="structure",
                                                              filename=f"{self.data_dir}/{self.output_mmCIF_file}")[0]
        structure_atoms = list(structure.get_atoms())
        selector = AtomSelector()
        io = PDB.PDBIO()
        io.set_structure(structure)
        kdtree = PDB.NeighborSearch(structure_atoms)
        print("ok")

        # load partial atomic charges estimation
        print("Loading patial atomic charges estimation... ", end="")
        if self.charges_estimation:
            charge_estimations = [float(x) for x in open(self.charges_estimation, "r").read().split()]
        else:
            charge_estimations = [0 for _ in range(len(structure_atoms))]
        total_charge = round(sum(charge_estimations))
        # creating charge attributes to make them easy to work with in Biopython library
        for atom, atomic_charge_estimation in zip(structure_atoms, charge_estimations):
            atom.charge_estimation = atomic_charge_estimation
            atom.cm5_charge = None
        print("ok")

        # calculate the charges for each atom using the cutoff approach.
        for calculated_atom_i, calculated_atom in enumerate(tqdm.tqdm(structure_atoms,
                                                                      desc="Charges calculation",
                                                                      unit="atoms",
                                                                      smoothing=0,
                                                                      delay=0.5,
                                                                      mininterval=0.5,
                                                                      maxinterval=0.5),
                                                            start=1):

            # To speed up the calculation, the charges of the hydrogen and oxygen atoms bound to one atom
            # are calculated together with the nearest other heavy atoms
            if calculated_atom.element == "H":
                continue
            elif calculated_atom.element == "O":
                if len(kdtree.search(center=calculated_atom.coord,
                                     radius=1.5,
                                     level="A")) <= 2:
                    continue

            substructure_data_dir = f"{self.data_dir}/sub_{calculated_atom_i}"
            system(f"mkdir {substructure_data_dir}")

            # definition of radii limiting the substructure
            # all atoms that are closer to the calculated atom than min_radius are included in the substructure
            # atoms more distant from the calculated atom than max_radius are never included in the substructure
            min_radius = 6
            max_radius = 12

            # xtb calculation may not converge
            # in this case we try the calculation four times with min_radius and max_radius increased
            calculation_converged = False
            while not calculation_converged:

                # create and save min_radius and max_radius substructures by biopython
                atoms_in_min_radius = kdtree.search(center=calculated_atom.coord,
                                                    radius=min_radius,
                                                    level="A")
                selector.full_ids = set([atom.full_id for atom in atoms_in_min_radius])
                io.save(file=f"{substructure_data_dir}/atoms_in_{min_radius}_angstroms.pdb",
                        select=selector)
                atoms_in_max_radius = kdtree.search(center=calculated_atom.coord,
                                                    radius=max_radius,
                                                    level="A")
                selector.full_ids = set([atom.full_id for atom in atoms_in_max_radius])
                io.save(file=f"{substructure_data_dir}/atoms_in_{max_radius}_angstroms.pdb",
                        select=selector)

                # load substructures by RDKit to determine bonds
                mol_min_radius = Chem.MolFromPDBFile(molFileName=f"{substructure_data_dir}/atoms_in_{min_radius}_angstroms.pdb",
                                                     removeHs=False,
                                                     sanitize=False)
                mol_min_radius_conformer = mol_min_radius.GetConformer()
                mol_max_radius = Chem.MolFromPDBFile(molFileName=f"{substructure_data_dir}/atoms_in_{max_radius}_angstroms.pdb",
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
                                carbons_with_broken_bonds_coord.append(mol_max_radius_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                                continue
                            else:
                                atoms_with_broken_bonds.append(bonded_atom)
                                substructure_coord_dict[(coord.x, coord.y, coord.z)] = bonded_atom

                # create substructure in Biopython library
                # we prefer to use kdtree because serial_id may be discontinuous in some pdbs files
                # for example, a structure with PDB code 107d and its serial numbers 218 and 445
                substructure_atoms = [kdtree.search(center=coord,
                                                    radius=0.1,
                                                    level="A")[0] for coord in substructure_coord_dict.keys()]
                selector.full_ids = set([atom.full_id for atom in substructure_atoms])
                io.save(file=f"{substructure_data_dir}/substructure.pdb",
                        select=selector)

                # add hydrogens to broken C-C bonds by openbabel
                system(f"cd {substructure_data_dir} ; obabel -iPDB -oPDB substructure.pdb -h > readded_hydrogens_substructure.pdb 2>/dev/null")
                with open(f"{substructure_data_dir}/readded_hydrogens_substructure.pdb") as readded_hydrogens_substructure_file:
                    atom_lines = [line for line in readded_hydrogens_substructure_file.readlines() if line[:4] in ["ATOM", "HETA"]]
                    original_atoms_lines = atom_lines[:len(substructure_atoms)]
                    added_hydrogens_lines = atom_lines[len(substructure_atoms):]
                with open(f"{substructure_data_dir}/repaired_substructure.pdb", "w") as repaired_substructure_file:
                    repaired_substructure_file.write("".join(original_atoms_lines))
                    for added_hydrogen_line in added_hydrogens_lines:
                        added_hydrogen_coord = (float(added_hydrogen_line[30:38]),
                                                float(added_hydrogen_line[38:46]),
                                                float(added_hydrogen_line[46:54]))
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
                    calculation_converged = True
                except ValueError:  # charge calculation failed
                    min_radius += 1
                    max_radius += 1
                if max_radius > 15:
                    break

            if calculation_converged:

                # define for which atoms we have calculated the charges
                calculated_atoms = set([calculated_atom])
                for near_atom in kdtree.search(center=calculated_atom.coord,
                                               radius=1.5,
                                               level="A"):
                    if near_atom.element == "H":
                        calculated_atoms.add(near_atom)
                    elif near_atom.element == "O":
                        if len(kdtree.search(center=near_atom.coord,
                                             radius=1.5,
                                             level="A")) <= 2:
                            calculated_atoms.add(near_atom)
                calculated_atoms_full_ids = set([calculated_atom.full_id for calculated_atom in calculated_atoms])

                # write the charges into the Biopython structure
                # we prefer to use kdtree because serial_id may be discontinuous in some pdbs files
                # for example, a structure with PDB code 107d and its serial numbers 218 and 445
                calculated_substructure = PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                                  file=f"{substructure_data_dir}/substructure.pdb")[0]
                for calculated_substructure_atom_i, substructure_atom in enumerate(calculated_substructure.get_atoms()):
                    if substructure_atom.full_id in calculated_atoms_full_ids:
                        charge = float(xtb_output_file_lines[charge_headline_index + calculated_substructure_atom_i + 1][19:28])
                        kdtree.search(center=substructure_atom.coord,
                                      radius=0.1,
                                      level="A")[0].cm5_charge = charge

            if self.delete_auxiliary_files:
                system(f"rm -r {substructure_data_dir}")

        # create final array of charges
        cm5_charges = [atom.cm5_charge for atom in structure_atoms]
        cm5_charges_sum = sum([charge for charge in cm5_charges if charge]) # filter None values from failed xtb calculations
        correction = (cm5_charges_sum - total_charge) / len(cm5_charges)
        self.cm5_charges = [round(charge - correction, 5) if charge else None for charge in cm5_charges]

        # todo log residues with some non calculated atoms

    def write_charges_to_files(self):
        print("Writing charges to files... ", end="")
        with open(f"{self.data_dir}/charges.txt", "w") as charges_file:
            charges_string = " ".join([str(x) for x in self.cm5_charges])
            charges_file.write(charges_string)

        # write charges to mmCIF file
        structure = gemmi.cif.read_file(f"{self.data_dir}/{self.output_mmCIF_file}")
        block = structure.sole_block()
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
        block.write_file(f"{self.data_dir}/{self.output_mmCIF_file}")
        print("ok\n")


if __name__ == "__main__":
    args = load_arguments()
    calculator = ChargeCalculator(input_mmCIF_file=args.input_mmCIF_file,
                                  charges_estimation=args.charges_estimation,
                                  output_mmCIF_file=args.output_mmCIF_file,
                                  data_dir=args.data_dir,
                                  delete_auxiliary_files=args.delete_auxiliary_files)
    calculator.calculate_charges()
    calculator.write_charges_to_files()
