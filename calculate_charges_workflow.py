#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict
from os import path, system, listdir

from phases.charge_calculator import ChargeCalculator
from phases.structure_preparer import StructurePreparer
from phases.hydrogen_optimiser import HydrogenOptimiser


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
        self.warnings = defaultdict(list)

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
        self.warnings[(chain, int(resnum), resname)].append(warning)

    def write_warnings(self):
        json_warnings = []
        for (chain_id, residue_id, residue_name), warnings in sorted(self.warnings.items()):
            json_warnings.append({"chain_id": chain_id,
                                  "residue_id": residue_id,
                                  "residue_name": residue_name,
                                  "warning": " ".join(warnings)})
        with open(self.warning_file, 'w') as warning_file:
            warning_file.write(json.dumps(json_warnings, indent=4))

if __name__ == "__main__":
    args = load_arguments()

    # prepare directories to store data
    results_directory = f"{args.data_dir}/results_{path.basename(args.PDB_file)[:-4].lower()}"
    if not path.exists(args.data_dir):
        system(f"mkdir {args.data_dir}")
    system(f"mkdir {args.data_dir}/input_PDB; "
           f"mkdir {results_directory}; "
           f"cp {args.PDB_file} {args.data_dir}/input_PDB")

    residual_warnings_file = f"{results_directory}/residual_warnings.json"
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

    # optimize added hydrogens
    hydrogen_optimiser_input = f"{structure_preparer_data_directory}/{structure_preparer_output}"
    hydrogen_optimiser_data_directory = f"{args.data_dir}/hydrogen_optimiser"
    hydrogen_optimiser_output = f"{path.basename(args.PDB_file)[:-4]}_optimisedH.cif"
    hydrogen_optimiser = HydrogenOptimiser(input_mmCIF_file=hydrogen_optimiser_input,
                                           logger=logger,
                                           output_mmCIF_file=hydrogen_optimiser_output,
                                           data_dir=hydrogen_optimiser_data_directory,
                                           delete_auxiliary_files=args.delete_auxiliary_files)
    hydrogen_optimiser.optimise()

    # calculate partial atomic charges
    charge_calculator_input = f"{hydrogen_optimiser_data_directory}/{hydrogen_optimiser_output}"
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

    logger.write_warnings()
