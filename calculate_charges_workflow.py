import argparse
from os import path, system

from phases.charge_calculator import ChargeCalculator
from phases.structure_preparer import StructurePreparer


def load_arguments():
    print("\nParsing arguments... ",
          end="")
    parser = argparse.ArgumentParser()
    parser.add_argument("--mmCIF_file",
                        help="mmCIF file with protein structure.",
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
    if not path.isfile(args.mmCIF_file):
        exit(f"\nERROR! File {args.mmCIF_file} does not exist!\n")
    if path.exists(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument!\n")
    print("ok")
    return args


if __name__ == "__main__":
    args = load_arguments()

    # check number of chains
    # Number of chains is greater than 52. Moleculekit (pdb2pqr) is not able to add hydrogens to structure.

    # prepare directories to store data
    system(f"mkdir {args.data_dir}; "
           f"mkdir {args.data_dir}/input_mmCIF; "
           f"cp {args.mmCIF_file} {args.data_dir}/input_mmCIF")

    structure_preparer = StructurePreparer(mmCIF_file=args.mmCIF_file,
                                           data_dir=f"{args.data_dir}/structure_preparer",
                                           delete_auxiliary_files=args.delete_auxiliary_files)
    structure_preparer.fix_structure()  # fix structure by PDBFixer
    structure_preparer.remove_hydrogens()
    structure_preparer.add_hydrogens_by_hydride()  # add hydrogens to the residues that the moleculekit can't process
    structure_preparer.add_hydrogens_by_moleculekit()  # add hydrogens to rest of structure

    calculator = ChargeCalculator(mmCIF_file=f"{args.data_dir}/structure_preparer/combined.cif",
                                  data_dir=f"{args.data_dir}/charges_calculator",
                                  delete_auxiliary_files=args.delete_auxiliary_files)
    calculator.calculate_charges()
    calculator.write_charges_to_files()
