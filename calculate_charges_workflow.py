import argparse
from os import path, system

from phases.charge_calculator import ChargeCalculator
from phases.structure_preparer import StructurePreparer


def load_arguments():
    print("\nParsing arguments... ",
          end="")
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
    return args


if __name__ == "__main__":
    args = load_arguments()

    # prepare directories to store data
    system(f"mkdir {args.data_dir}; "
           f"mkdir {args.data_dir}/input_PDB; "
           f"cp {args.PDB_file} {args.data_dir}/input_PDB")

    # prepare structure for main calculation of partial atomic charges
    structure_preparer_input = args.PDB_file
    structure_preparer_data_directory = f"{args.data_dir}/structure_preparer"
    structure_preparer_output = f"{path.basename(args.PDB_file)[:-4]}_prepared.cif"
    structure_preparer = StructurePreparer(input_PDB_file=structure_preparer_input,
                                           data_dir=structure_preparer_data_directory,
                                           output_mmCIF_file=structure_preparer_output,
                                           delete_auxiliary_files=args.delete_auxiliary_files,
                                           save_charges_estimation=True)
    structure_preparer.fix_structure()
    structure_preparer.remove_hydrogens()
    structure_preparer.add_hydrogens_by_hydride()
    structure_preparer.add_hydrogens_by_moleculekit()

    # calculate parcial atomic charges
    charge_calculator_input = f"{structure_preparer_data_directory}/{structure_preparer_output}"
    charge_calculator_data_directory = f"{args.data_dir}/charge_calculator"
    charge_calculator_output = f"{path.basename(args.PDB_file)[:-4]}.cif"
    charges_estimation = f"{structure_preparer_data_directory}/estimated_charges.txt"
    charge_calculator = ChargeCalculator(input_mmCIF_file=charge_calculator_input,
                                         charges_estimation=charges_estimation,
                                         output_mmCIF_file=charge_calculator_output,
                                         data_dir=charge_calculator_data_directory,
                                         delete_auxiliary_files=args.delete_auxiliary_files)
    charge_calculator.calculate_charges()
    charge_calculator.write_charges_to_files()
