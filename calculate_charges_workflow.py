#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from os import path, system, listdir

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
    parser.add_argument("--CCD_file",
                        help="SDF file with Chemical Component Dictionary.",
                        type=str,
                        required=True)
    parser.add_argument("--delete_auxiliary_files",
                        help="Auxiliary calculation files can be large. With this argument, "
                             "the auxiliary files will be continuously deleted during the calculation.",
                        action="store_true")

    args = parser.parse_args()
    if not path.isfile(args.PDB_file):
        exit(f"\nERROR! File {args.PDB_file} does not exist!\n")
    if path.exists(args.data_dir) and listdir(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists and is not empty. "
             f"Remove existed directory or change --data_dir argument!\n")
    print("ok")
    return args

class Logger:
    """This class handles output and warnings for residues and stores them in a defined files."""
    def __init__(self,
                 output_file: str,
                 warning_file: str):
        self.output_file = output_file
        self.warning_file = warning_file

    def print(self,
              text: str,
              end: str="\n",
              silence=False):
        if not silence:
            print(text, end=end)
        with open(self.output_file, "a") as output_file:
            output_file.write(f"{text}{end}")

    def add_warning(self,
                    chain: str,
                    resname: str,
                    resnum: str,
                    warning: str):
        with open(self.warning_file, "a") as log_file:
            log_file.write(f"{chain};{resnum};{resname};{warning}\n")

if __name__ == "__main__":
    args = load_arguments()

    # prepare directories to store data
    results_directory = f"{args.data_dir}/results_{path.basename(args.PDB_file)[:-4].lower()}"
    if not path.exists(args.data_dir):
        system(f"mkdir {args.data_dir}")
    system(f"mkdir {args.data_dir}/input_PDB; "
           f"mkdir {results_directory}; "
           f"cp {args.PDB_file} {args.data_dir}/input_PDB")

    residual_warnings_file = f"{results_directory}/residual_warnings.csv"
    logger = Logger(output_file=f"{results_directory}/output.txt",
                    warning_file=residual_warnings_file)

    # prepare structure for main calculation of partial atomic charges
    structure_preparer_input = args.PDB_file
    structure_preparer_data_directory = f"{args.data_dir}/structure_preparer"
    structure_preparer_output = f"{path.basename(args.PDB_file)[:-4]}_prepared.cif"
    structure_preparer = StructurePreparer(input_PDB_file=structure_preparer_input,
                                           CCD_file=args.CCD_file,
                                           logger=logger,
                                           data_dir=structure_preparer_data_directory,
                                           output_mmCIF_file=structure_preparer_output,
                                           delete_auxiliary_files=args.delete_auxiliary_files,
                                           save_charges_estimation=True)
    structure_preparer.fix_structure()
    structure_preparer.remove_hydrogens()
    structure_preparer.add_hydrogens_by_hydride()
    structure_preparer.add_hydrogens_by_moleculekit()

    # calculate partial atomic charges
    charge_calculator_input = f"{structure_preparer_data_directory}/{structure_preparer_output}"
    charge_calculator_data_directory = f"{args.data_dir}/charge_calculator"
    charge_calculator_output = f"{path.basename(args.PDB_file)[:-4]}.cif"
    charges_estimation = f"{structure_preparer_data_directory}/estimated_charges.txt"
    charge_calculator = ChargeCalculator(input_mmCIF_file=charge_calculator_input,
                                         charges_estimation=charges_estimation,
                                         logger=logger,
                                         output_mmCIF_file=charge_calculator_output,
                                         data_dir=charge_calculator_data_directory,
                                         delete_auxiliary_files=args.delete_auxiliary_files)
    charge_calculator.calculate_charges()
    charge_calculator.write_charges_to_files()

    system(f"cp {charge_calculator_data_directory}/{charge_calculator_output} {results_directory}")

    # merge and sort warnings for residues
    if path.exists(residual_warnings_file):
        residual_warnings = defaultdict(list)
        with open(residual_warnings_file, "r") as f:
            for chain, res_num, res_name, warning in list(csv.reader(f, delimiter=";")):
                residual_warnings[(chain, int(res_num), res_name)].append(warning)
        warnings = sorted([(chain, res_num, res_name, " ".join(warning)) for (chain, res_num, res_name), warning in residual_warnings.items()])
        with open(residual_warnings_file, 'w') as f:
            csv_writer = csv.writer(f, delimiter=";", quoting=csv.QUOTE_NONNUMERIC)
            for line in warnings:
                csv_writer.writerow(line)

# todo openmm náboje
    # todo optimalizace vodíků
