import argparse
from math import dist
from os import path, system
from time import time


from biotite.structure import BondType, BondList, residue_iter, get_residue_starts
from biotite.structure.io import pdbx as mmCIF

import gemmi
import hydride
import numpy as np
from Bio import PDB
from biotite.structure.residues import get_residue_starts_for
from dimorphite_dl import DimorphiteDL
from moleculekit import molecule as molit
from moleculekit.tools.preparation import systemPrepare
from numba.cuda import atomic
from openmm.app import PDBxFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import rdFMCS


def load_arguments():
    print("\nParsing arguments... ", end="")
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


def fix_structure(input_mmCIF_file: str,
                  data_dir: str):
    """
    mmCIF file is fixed by tool PDBFixer.
    https://github.com/openmm/pdbfixer

    PDBFixer solves common problems in protein structure files.
    It selects the first model, alternative locations, fills in missing heavy atoms, etc.
    """

    print("Fixing structure... ", end="")
    fixer = PDBFixer(filename=input_mmCIF_file)
    set_of_resnames = set([x.name for x in fixer.topology.residues()])

    # download templates for heteroresidues
    fixer_available_resnames = set(fixer.templates.keys())
    for resname in set_of_resnames:
        if resname not in fixer_available_resnames:
            try:
                fixer.downloadTemplate(resname)
            except:
                exit(f"ERROR! Old version of pdbfixer installed or heteroresiduum {resname} does not exist!")

    # add heavy atoms
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    # with open(f"{self.data_dir}/logs/added_heavy_atoms.txt", "w") as added_heavy_atoms_file: # log it TODO
    #     for residue, residue_list in fixer.missingAtoms.items():
    #         for atom in residue_list:
    #             added_heavy_atoms_file.write(f"{residue} {atom}\n")
    fixer.addMissingAtoms()

    # write fixed structure to file
    PDBxFile.writeFile(fixer.topology, fixer.positions, open(f"{data_dir}/fixed.cif", 'w'), keepIds=True)

    # PDBFixer removes values from a struct_conn block
    # and therefore only atom_site block is overwritten from PDBFixer to fixed.cif.
    biotite_mmCIF_file = mmCIF.CIFFile.read(input_mmCIF_file)
    fixed_mmCIF_file = mmCIF.CIFFile.read(f"{data_dir}/fixed.cif")
    biotite_mmCIF_file.block["atom_site"] = fixed_mmCIF_file.block["atom_site"]
    biotite_mmCIF_file.write(f"{data_dir}/fixed.cif")
    # these three lines can be probably removed, after biotite 1.0.2 will be released
    mmcif_string = open(f"{data_dir}/fixed.cif").read()
    repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
    open(f"{data_dir}/fixed.cif", "w").write(repaired_mmcif_string)
    print("ok")


class StructurePreparer:
    def __init__(self,
                 input_mmCIF_file: str,
                 data_dir: str,
                 pH: float=7.2):
        self.mmCIF_file = mmCIF.CIFFile.read(input_mmCIF_file)
        self.data_dir = data_dir
        self.pH = pH

        self.remove_hydrogens()
        self.protonate_by_hydride() # add hydrogens to the residues that the moleculekit can't process
        self.protonate_by_moleculekit() # add hydrogens to rest of structure



    def remove_hydrogens(self):
        print("Removing hydrogens ... ", end="")
        protein = mmCIF.get_structure(self.mmCIF_file,
                                              model=1,
                                              extra_fields=["b_factor", "occupancy", "charge"],
                                              include_bonds=True)
        protein_without_hydrogens = protein[protein.element != "H"]
        mmCIF.set_structure(self.mmCIF_file, protein_without_hydrogens)
        print("ok")


    def protonate_by_hydride(self):
        """
        This function is based on the biotite, hydride, RDKit and dimorphite_dl libraries.
        https://github.com/biotite-dev/biotite
        https://github.com/biotite-dev/hydride
        https://github.com/rdkit/rdkit
        https://github.com/durrantlab/dimorphite_dl

        The heteroresidues are protonated by the hydride library, which is built on top of the biotite library.
        Prior to the actual protonation, the formal charges are loaded from CCD dictionary and extended by the dimorphite_dl.
        Dimorphite_dl formal charges are mapped to CCD molecule by RDKit library.
        """

        print("Adding hydrogens to heteroresidues... ", end="")
        # pdb2pqr is part of moleculekit
        residues_processed_by_pdb2pqr = set(
            ['004', '03Y', '0A1', '0AF', '0BN', '1MH', '2AS', '2GX', '2ML', '2MR', '4IN', '4PH', '4PQ', '5JP', 'AA4',
             'ABA', 'AHP', 'ALA', 'ALC', 'ALN', 'ALY', 'APD', 'ARG', 'ASN', 'ASP', 'BB8', 'BCS', 'BTK', 'CCS', 'CGU',
             'CSA', 'CSO', 'CSP', 'CSS', 'CYS', 'D4P', 'DA2', 'DAB', 'DAH', 'DPP', 'ESC', 'FGL', 'GHG', 'GLN', 'GLU',
             'GLY', 'GME', 'GNC', 'HHK', 'HIS', 'HLU', 'HLX', 'HOX', 'HPE', 'HQA', 'HTR', 'HYP', 'I2M', 'IGL', 'IIL',
             'ILE', 'IML', 'KYN', 'LEU', 'LME', 'LMQ', 'LYS', 'LYZ', 'M3L', 'ME0', 'MEA', 'MEN', 'MEQ', 'MET', 'MLE',
             'MLY', 'MLZ', 'MME', 'MMO', 'MVA', 'NAL', 'NCY', 'NLE', 'NVA', 'NZC', 'OCY', 'OMX', 'ONL', 'ORM', 'P1L',
             'PCA', 'PHE', 'PRK', 'PRO', 'PTR', 'SEP', 'SER', 'THR', 'TPO', 'TRO', 'TRP', 'TY2', 'TYQ', 'TYR', 'VAL',
             'WAT', 'YCM', 'YNM', 'RA', 'RC', 'RG', 'DT', 'RU', 'ASH', 'CYM', 'CYX', 'GLH', 'HSE', 'HSD', 'HSP', 'HID',
             'HIE', 'HIP', 'AR0', 'LYN', 'TYM', 'C004', 'C03Y', 'C0A1', 'C0AF', 'C0BN', 'C1MH', 'C2AS', 'C2GX', 'C2ML',
             'C2MR', 'C4IN', 'C4PH', 'C4PQ', 'C5JP', 'CAA4', 'CABA', 'CAHP', 'CALA', 'CALC', 'CALN', 'CALY', 'CAPD',
             'CARG', 'CASN', 'CASP', 'CBB8', 'CBCS', 'CBTK', 'CCCS', 'CCGU', 'CCSA', 'CCSO', 'CCSP', 'CCSS', 'CCYS',
             'CD4P', 'CDA2', 'CDAB', 'CDAH', 'CDPP', 'CESC', 'CFGL', 'CGHG', 'CGLN', 'CGLU', 'CGLY', 'CGME', 'CGNC',
             'CHHK', 'CHIS', 'CHLU', 'CHLX', 'CHOX', 'CHPE', 'CHQA', 'CHTR', 'CHYP', 'CI2M', 'CIGL', 'CIIL', 'CILE',
             'CIML', 'CKYN', 'CLEU', 'CLME', 'CLMQ', 'CLYS', 'CLYZ', 'CM3L', 'CME0', 'CMEA', 'CMEN', 'CMEQ', 'CMET',
             'CMLE', 'CMLY', 'CMLZ', 'CMME', 'CMMO', 'CMVA', 'CNAL', 'CNCY', 'CNLE', 'CNVA', 'CNZC', 'COCY', 'COMX',
             'CONL', 'CORM', 'CP1L', 'CPCA', 'CPHE', 'CPRK', 'CPRO', 'CPTR', 'CSEP', 'CSER', 'CTHR', 'CTPO', 'CTRO',
             'CTRP', 'CTY2', 'CTYQ', 'CTYR', 'CVAL', 'CWAT', 'CYCM', 'CYNM', 'CASH', 'CCYM', 'CCYX', 'CGLH', 'CHSE',
             'CHSD', 'CHSP', 'CHID', 'CHIE', 'CHIP', 'CAR0', 'CLYN', 'CTYM', 'NEUTRAL-C004', 'NEUTRAL-C03Y',
             'NEUTRAL-C0A1', 'NEUTRAL-C0AF', 'NEUTRAL-C0BN', 'NEUTRAL-C1MH', 'NEUTRAL-C2AS', 'NEUTRAL-C2GX',
             'NEUTRAL-C2ML', 'NEUTRAL-C2MR', 'NEUTRAL-C4IN', 'NEUTRAL-C4PH', 'NEUTRAL-C4PQ', 'NEUTRAL-C5JP',
             'NEUTRAL-CAA4', 'NEUTRAL-CABA', 'NEUTRAL-CAHP', 'NEUTRAL-CALA', 'NEUTRAL-CALC', 'NEUTRAL-CALN',
             'NEUTRAL-CALY', 'NEUTRAL-CAPD', 'NEUTRAL-CARG', 'NEUTRAL-CASN', 'NEUTRAL-CASP', 'NEUTRAL-CBB8',
             'NEUTRAL-CBCS', 'NEUTRAL-CBTK', 'NEUTRAL-CCCS', 'NEUTRAL-CCGU', 'NEUTRAL-CCSA', 'NEUTRAL-CCSO',
             'NEUTRAL-CCSP', 'NEUTRAL-CCSS', 'NEUTRAL-CCYS', 'NEUTRAL-CD4P', 'NEUTRAL-CDA2', 'NEUTRAL-CDAB',
             'NEUTRAL-CDAH', 'NEUTRAL-CDPP', 'NEUTRAL-CESC', 'NEUTRAL-CFGL', 'NEUTRAL-CGHG', 'NEUTRAL-CGLN',
             'NEUTRAL-CGLU', 'NEUTRAL-CGLY', 'NEUTRAL-CGME', 'NEUTRAL-CGNC', 'NEUTRAL-CHHK', 'NEUTRAL-CHIS',
             'NEUTRAL-CHLU', 'NEUTRAL-CHLX', 'NEUTRAL-CHOX', 'NEUTRAL-CHPE', 'NEUTRAL-CHQA', 'NEUTRAL-CHTR',
             'NEUTRAL-CHYP', 'NEUTRAL-CI2M', 'NEUTRAL-CIGL', 'NEUTRAL-CIIL', 'NEUTRAL-CILE', 'NEUTRAL-CIML',
             'NEUTRAL-CKYN', 'NEUTRAL-CLEU', 'NEUTRAL-CLME', 'NEUTRAL-CLMQ', 'NEUTRAL-CLYS', 'NEUTRAL-CLYZ',
             'NEUTRAL-CM3L', 'NEUTRAL-CME0', 'NEUTRAL-CMEA', 'NEUTRAL-CMEN', 'NEUTRAL-CMEQ', 'NEUTRAL-CMET',
             'NEUTRAL-CMLE', 'NEUTRAL-CMLY', 'NEUTRAL-CMLZ', 'NEUTRAL-CMME', 'NEUTRAL-CMMO', 'NEUTRAL-CMVA',
             'NEUTRAL-CNAL', 'NEUTRAL-CNCY', 'NEUTRAL-CNLE', 'NEUTRAL-CNVA', 'NEUTRAL-CNZC', 'NEUTRAL-COCY',
             'NEUTRAL-COMX', 'NEUTRAL-CONL', 'NEUTRAL-CORM', 'NEUTRAL-CP1L', 'NEUTRAL-CPCA', 'NEUTRAL-CPHE',
             'NEUTRAL-CPRK', 'NEUTRAL-CPRO', 'NEUTRAL-CPTR', 'NEUTRAL-CSEP', 'NEUTRAL-CSER', 'NEUTRAL-CTHR',
             'NEUTRAL-CTPO', 'NEUTRAL-CTRO', 'NEUTRAL-CTRP', 'NEUTRAL-CTY2', 'NEUTRAL-CTYQ', 'NEUTRAL-CTYR',
             'NEUTRAL-CVAL', 'NEUTRAL-CWAT', 'NEUTRAL-CYCM', 'NEUTRAL-CYNM', 'NEUTRAL-CASH', 'NEUTRAL-CCYM',
             'NEUTRAL-CCYX', 'NEUTRAL-CGLH', 'NEUTRAL-CHSE', 'NEUTRAL-CHSD', 'NEUTRAL-CHSP', 'NEUTRAL-CHID',
             'NEUTRAL-CHIE', 'NEUTRAL-CHIP', 'NEUTRAL-CAR0', 'NEUTRAL-CLYN', 'NEUTRAL-CTYM', 'N004', 'N03Y', 'N0A1',
             'N0AF', 'N0BN', 'N1MH', 'N2AS', 'N2GX', 'N2ML', 'N2MR', 'N4IN', 'N4PH', 'N4PQ', 'N5JP', 'NAA4', 'NABA',
             'NAHP', 'NALA', 'NALC', 'NALN', 'NALY', 'NAPD', 'NARG', 'NASN', 'NASP', 'NBB8', 'NBCS', 'NBTK', 'NCCS',
             'NCGU', 'NCSA', 'NCSO', 'NCSP', 'NCSS', 'NCYS', 'ND4P', 'NDA2', 'NDAB', 'NDAH', 'NDPP', 'NESC', 'NFGL',
             'NGHG', 'NGLN', 'NGLU', 'NGLY', 'NGME', 'NGNC', 'NHHK', 'NHIS', 'NHLU', 'NHLX', 'NHOX', 'NHPE', 'NHQA',
             'NHTR', 'NHYP', 'NI2M', 'NIGL', 'NIIL', 'NILE', 'NIML', 'NKYN', 'NLEU', 'NLME', 'NLMQ', 'NLYS', 'NLYZ',
             'NM3L', 'NME0', 'NMEA', 'NMEN', 'NMEQ', 'NMET', 'NMLE', 'NMLY', 'NMLZ', 'NMME', 'NMMO', 'NMVA', 'NNAL',
             'NNCY', 'NNLE', 'NNVA', 'NNZC', 'NOCY', 'NOMX', 'NONL', 'NORM', 'NP1L', 'NPCA', 'NPHE', 'NPRK', 'NPRO',
             'NPTR', 'NSEP', 'NSER', 'NTHR', 'NTPO', 'NTRO', 'NTRP', 'NTY2', 'NTYQ', 'NTYR', 'NVAL', 'NWAT', 'NYCM',
             'NYNM', 'NASH', 'NCYM', 'NCYX', 'NGLH', 'NHSE', 'NHSD', 'NHSP', 'NHID', 'NHIE', 'NHIP', 'NAR0', 'NLYN',
             'NTYM', 'NEUTRAL-N004', 'NEUTRAL-N03Y', 'NEUTRAL-N0A1', 'NEUTRAL-N0AF', 'NEUTRAL-N0BN', 'NEUTRAL-N1MH',
             'NEUTRAL-N2AS', 'NEUTRAL-N2GX', 'NEUTRAL-N2ML', 'NEUTRAL-N2MR', 'NEUTRAL-N4IN', 'NEUTRAL-N4PH',
             'NEUTRAL-N4PQ', 'NEUTRAL-N5JP', 'NEUTRAL-NAA4', 'NEUTRAL-NABA', 'NEUTRAL-NAHP', 'NEUTRAL-NALA',
             'NEUTRAL-NALC', 'NEUTRAL-NALN', 'NEUTRAL-NALY', 'NEUTRAL-NAPD', 'NEUTRAL-NARG', 'NEUTRAL-NASN',
             'NEUTRAL-NASP', 'NEUTRAL-NBB8', 'NEUTRAL-NBCS', 'NEUTRAL-NBTK', 'NEUTRAL-NCCS', 'NEUTRAL-NCGU',
             'NEUTRAL-NCSA', 'NEUTRAL-NCSO', 'NEUTRAL-NCSP', 'NEUTRAL-NCSS', 'NEUTRAL-NCYS', 'NEUTRAL-ND4P',
             'NEUTRAL-NDA2', 'NEUTRAL-NDAB', 'NEUTRAL-NDAH', 'NEUTRAL-NDPP', 'NEUTRAL-NESC', 'NEUTRAL-NFGL',
             'NEUTRAL-NGHG', 'NEUTRAL-NGLN', 'NEUTRAL-NGLU', 'NEUTRAL-NGLY', 'NEUTRAL-NGME', 'NEUTRAL-NGNC',
             'NEUTRAL-NHHK', 'NEUTRAL-NHIS', 'NEUTRAL-NHLU', 'NEUTRAL-NHLX', 'NEUTRAL-NHOX', 'NEUTRAL-NHPE',
             'NEUTRAL-NHQA', 'NEUTRAL-NHTR', 'NEUTRAL-NHYP', 'NEUTRAL-NI2M', 'NEUTRAL-NIGL', 'NEUTRAL-NIIL',
             'NEUTRAL-NILE', 'NEUTRAL-NIML', 'NEUTRAL-NKYN', 'NEUTRAL-NLEU', 'NEUTRAL-NLME', 'NEUTRAL-NLMQ',
             'NEUTRAL-NLYS', 'NEUTRAL-NLYZ', 'NEUTRAL-NM3L', 'NEUTRAL-NME0', 'NEUTRAL-NMEA', 'NEUTRAL-NMEN',
             'NEUTRAL-NMEQ', 'NEUTRAL-NMET', 'NEUTRAL-NMLE', 'NEUTRAL-NMLY', 'NEUTRAL-NMLZ', 'NEUTRAL-NMME',
             'NEUTRAL-NMMO', 'NEUTRAL-NMVA', 'NEUTRAL-NNAL', 'NEUTRAL-NNCY', 'NEUTRAL-NNLE', 'NEUTRAL-NNVA',
             'NEUTRAL-NNZC', 'NEUTRAL-NOCY', 'NEUTRAL-NOMX', 'NEUTRAL-NONL', 'NEUTRAL-NORM', 'NEUTRAL-NP1L',
             'NEUTRAL-NPCA', 'NEUTRAL-NPHE', 'NEUTRAL-NPRK', 'NEUTRAL-NPRO', 'NEUTRAL-NPTR', 'NEUTRAL-NSEP',
             'NEUTRAL-NSER', 'NEUTRAL-NTHR', 'NEUTRAL-NTPO', 'NEUTRAL-NTRO', 'NEUTRAL-NTRP', 'NEUTRAL-NTY2',
             'NEUTRAL-NTYQ', 'NEUTRAL-NTYR', 'NEUTRAL-NVAL', 'NEUTRAL-NWAT', 'NEUTRAL-NYCM', 'NEUTRAL-NYNM',
             'NEUTRAL-NASH', 'NEUTRAL-NCYM', 'NEUTRAL-NCYX', 'NEUTRAL-NGLH', 'NEUTRAL-NHSE', 'NEUTRAL-NHSD',
             'NEUTRAL-NHSP', 'NEUTRAL-NHID', 'NEUTRAL-NHIE', 'NEUTRAL-NHIP', 'NEUTRAL-NAR0', 'NEUTRAL-NLYN',
             'NEUTRAL-NTYM', 'HOH', 'DA', 'DA3', 'DA5', 'RA3', 'RA5', 'DC', 'DC3', 'DC5', 'RC3', 'RC5', 'DG', 'DG3',
             'DG5', 'RG3', 'RG5', 'DT3', 'RU3', 'RU5'])
        # shortcuts for RNA, also protonated by pdb2pqr, defined in RNA_MAPPING
        residues_processed_by_pdb2pqr.update(["A", "C", "G", "U"])

        protein = mmCIF.get_structure(self.mmCIF_file,
                                              model=1,
                                              extra_fields=["b_factor", "occupancy", "charge"],
                                              include_bonds=True)


        protein.set_annotation("hydride_mask", [False for _ in range(len(protein))])

        residues_processed_by_hydride = [(res, res_start_i) for res, res_start_i in zip(residue_iter(protein), get_residue_starts(protein)) if res.res_name[0] not in residues_processed_by_pdb2pqr]
        # todo!! je to tak správně? titruje je propka?

        if residues_processed_by_hydride:
            # load formal charges for ligand from CCD. Add other formal charges by dimorphite
            dimorphite = DimorphiteDL(min_ph=self.pH,
                                      max_ph=self.pH,
                                      max_variants=1,
                                      label_states=False,
                                      pka_precision=0.001)
            residues_processed_by_hydride_names = set([res.res_name[0] for res, _ in residues_processed_by_hydride])
            residues_processed_by_hydride_formal_charges = {}
            for CCD_mol_sdf in open("CCD/Components-pub.sdf", "r").read().split("$$$$\n"):
                mol_name = CCD_mol_sdf.partition("\n")[0]
                if mol_name in residues_processed_by_hydride_names:
                    supplier = Chem.SDMolSupplier()
                    supplier.SetData(CCD_mol_sdf)
                    CCD_mol = next(supplier)
                    if CCD_mol is None:
                        residues_processed_by_hydride_formal_charges[mol_name] = (None,
                                                                                  "The structure cannot be loaded by RDKit.")
                    else:
                        CCD_mol = Chem.RemoveAllHs(CCD_mol)
                        CCD_mol_formal_charges = [atom.GetFormalCharge() for atom in CCD_mol.GetAtoms()]
                        CCD_mol_smiles = Chem.MolToSmiles(CCD_mol)
                        dimorphite_smiles = dimorphite.protonate(CCD_mol_smiles)[0]
                        dimorphite_mol = Chem.MolFromSmiles(dimorphite_smiles)
                        dimorphite_mol = Chem.RemoveAllHs(dimorphite_mol)
                        params = rdFMCS.MCSParameters()
                        params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                        params.BondTyper = rdFMCS.BondCompare.CompareOrder
                        params.BondCompareParameters.RingMatchesRingOnly = True
                        params.BondCompareParameters.CompleteRingsOnly = True
                        params.AtomCompareParameters.MatchFormalCharge = False
                        params.Timeout = 60
                        res = rdFMCS.FindMCS([CCD_mol, dimorphite_mol], params)
                        atom_indices_map = [x[1] for x in
                                            sorted(list(zip(dimorphite_mol.GetSubstructMatch(res.queryMol),
                                                            CCD_mol.GetSubstructMatch(res.queryMol))))]
                        if len(dimorphite_mol.GetAtoms()) != len(atom_indices_map):
                            residues_processed_by_hydride_formal_charges[mol_name] = (CCD_mol_formal_charges,
                                                                                      "Atom mapping of structures from CCD and dimorphite failed. The formal charges are taken from the CCD.")
                        else:
                            dimorphite_formal_charges = [0 for _ in dimorphite_mol.GetAtoms()]
                            for atom, a_i in zip(dimorphite_mol.GetAtoms(), atom_indices_map):
                                dimorphite_formal_charges[a_i] = atom.GetFormalCharge()
                            combined_charges = []
                            for CCD_atom, CCD_formal_charge, dimorphite_formal_charge in zip(CCD_mol.GetAtoms(),
                                                                                             CCD_mol_formal_charges,
                                                                                             dimorphite_formal_charges):
                                if CCD_formal_charge != 0:
                                    combined_charges.append(CCD_formal_charge)
                                elif dimorphite_formal_charge != 0:
                                    bonded_atoms = list(CCD_atom.GetNeighbors())
                                    bonded_atoms_over_two_bonds = []
                                    for bonded_atom in bonded_atoms:
                                        bonded_atoms_over_two_bonds.extend(bonded_atom.GetNeighbors())
                                    if all([atom.GetFormalCharge() == 0 for atom in
                                            bonded_atoms + bonded_atoms_over_two_bonds]):  # CCD_atom is already in bonded_atoms_over_two_bonds
                                        combined_charges.append(dimorphite_formal_charge)
                                    else:
                                        combined_charges.append(0)
                                else:
                                    combined_charges.append(0)
                            residues_processed_by_hydride_formal_charges[mol_name] = (combined_charges,
                                                                                      "The formal charges are taken from CCD and dimorphite.")


            for res, res_start_i in residues_processed_by_hydride:
                try:
                    formal_charges, log = residues_processed_by_hydride_formal_charges[res.res_name[0]]
                except KeyError:
                    # zalogovat, že ligand není v CCD, necháváme neutrální
                    continue
                # logovat podle logu! 1) není načetnutelný, nešlo dimorphite, šlo dimorphite

                if len(res) != len(formal_charges):
                    # zalogovat, že ligand má jiný počet atomů než v CCD, necháváme neutrální
                    continue

                # kontrolovat pořadí!!! TODO
                # toto souvisí i s tím, jak budeme nakládat s ligandy, které budou mít menší počet atomů než v CCD, domyslet!



                for atom_i, (atom, formal_charge) in enumerate(zip(res, formal_charges), start=res_start_i):
                    protein.charge[atom_i] = formal_charge
                    protein.hydride_mask[atom_i] = True
                # zalogovat, že ligandu jsou přiřazeny náboje buď podle CCD a nebo podle dimorphite_dl


        bond_array = protein.bonds.as_array()
        unknown_order_mask = bond_array[:, 2] == BondType.ANY
        if unknown_order_mask.any():
            print("For some bonds the bond order is unknown, hence single bonds are assumed")
            bond_array[unknown_order_mask, 2] = BondType.SINGLE
        protein.bonds = BondList(protein.bonds.get_atom_count(), bond_array)


        protein_with_hydrogens, _ = hydride.add_hydrogen(protein, mask=protein.hydride_mask)
        mmCIF.set_structure(self.mmCIF_file, protein_with_hydrogens)
        bond_array = protein_with_hydrogens.bonds.as_array()
        residue_starts_1, residue_starts_2 = (
            get_residue_starts_for(protein_with_hydrogens, bond_array[:, :2].flatten()).reshape(-1, 2).T)
        interresidual_bonds = bond_array[residue_starts_1 != residue_starts_2]
        ptnr1_auth_asym_id = []
        ptnr2_auth_asym_id = []
        ptnr1_auth_seq_id = []
        ptnr2_auth_seq_id = []
        protein_auth_asym_id = self.mmCIF_file.block["atom_site"]["auth_asym_id"].as_array()
        protein_auth_seq_id = self.mmCIF_file.block["atom_site"]["auth_seq_id"].as_array()
        for a1, a2, _ in interresidual_bonds:
            ptnr1_auth_asym_id.append(protein_auth_asym_id[a1])
            ptnr2_auth_asym_id.append(protein_auth_asym_id[a2])
            ptnr1_auth_seq_id.append(protein_auth_seq_id[a1])
            ptnr2_auth_seq_id.append(protein_auth_seq_id[a2])
        self.mmCIF_file.block["struct_conn"]["ptnr1_auth_asym_id"] = ptnr1_auth_asym_id
        self.mmCIF_file.block["struct_conn"]["ptnr2_auth_asym_id"] = ptnr2_auth_asym_id
        self.mmCIF_file.block["struct_conn"]["ptnr1_auth_seq_id"] = ptnr1_auth_seq_id
        self.mmCIF_file.block["struct_conn"]["ptnr2_auth_seq_id"] = ptnr2_auth_seq_id

        self.label_asym_ids = sorted(set(self.mmCIF_file.block["atom_site"]["label_asym_id"].as_array()))
        self.auth_asym_ids = sorted(set(self.mmCIF_file.block["atom_site"]["auth_asym_id"].as_array()))
        self.one_symbol_chains = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        if len(self.label_asym_ids) > 52 or len(self.auth_asym_ids) > 52:
            exit(
                "Number of chains is greater than 52. Moleculekit (pdb2pqr) is not able to add hydrogens to structure.")

        presented_label_asym_id = {key: self.one_symbol_chains[index] for index, key in enumerate(self.label_asym_ids)}
        presented_auth_asym_id = {key: self.one_symbol_chains[index] for index, key in enumerate(self.auth_asym_ids)}
        self.mmCIF_file.block["atom_site"]["label_asym_id"] = [presented_label_asym_id[key] for key in
                                                                       self.mmCIF_file.block["atom_site"][
                                                                           "label_asym_id"].as_array()]
        self.mmCIF_file.block["atom_site"]["auth_asym_id"] = [presented_auth_asym_id[key] for key in
                                                                      self.mmCIF_file.block["atom_site"][
                                                                          "auth_asym_id"].as_array()]

        self.mmCIF_file.write(f"{self.data_dir}/protonated_ligands.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/protonated_ligands.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/protonated_ligands.cif", "w").write(repaired_mmcif_string)
        print("ok")

    def protonate_by_moleculekit(self):
        print("Adding hydrogens to protein... ", end="")
        molecule = molit.Molecule(f"{self.data_dir}/protonated_ligands.cif")
        prepared_molecule = systemPrepare(molecule,
                                          pH=self.pH,
                                          ignore_ns_errors=True,
                                          _molkit_ff=False)
        prepared_molecule.write(f"{self.data_dir}/protonated_protein.cif", )
        pdb2pqr_charges = np.nan_to_num(prepared_molecule.charge)
        if all(chg == 0 for chg in pdb2pqr_charges):
            exit("ERROR! Moleculekit is not modified!")
            # Modification of moleculekit:
            # .../moleculekit/tools/preparation.py
            #  line 827 ("ffcharge", "charge")
            # https://github.com/Acellera/mdoleculekit/issues/136

        # print(sum(pdb2pqr_charges)) # jakto že je to takto divné? todo dopsat warning!
        self.atomic_charges = pdb2pqr_charges + prepared_molecule.formalcharge
        print("ok")


class SelectAtoms(PDB.Select):
    def accept_atom(self, atom):
        if atom.full_id in self.full_ids:
            return 1
        else:
            return 0


class ChargesCalculator:
    def __init__(self,
                 mmCIF_file: str,
                 atomic_charge_estimations: np.array,
                 data_dir: str,
                 delete_auxiliary_files: bool):
        print("Loading structure... ", end="")
        self.mmCIF_file = mmCIF_file
        self.data_dir = data_dir
        self.atomic_charge_estimations = atomic_charge_estimations
        self.delete_auxiliary_files = delete_auxiliary_files
        print("ok")


    def calculate_charges(self):
        structure = PDB.MMCIFParser(QUIET=True).get_structure("structure", self.mmCIF_file)[0]
        kdtree = PDB.NeighborSearch(list(structure.get_atoms()))

        for atom, atomic_charge_estimation in zip(structure.get_atoms(), self.atomic_charge_estimations):
            atom.atomic_charge_estimation = atomic_charge_estimation
            atom.cm5_charge = None

        selector = SelectAtoms()
        io = PDB.PDBIO()
        io.set_structure(structure)
        radius = 6
        for i, atom in enumerate(structure.get_atoms()):
            if atom.element == "H":
                continue
            if atom.element == "O":
                if len(kdtree.search(atom.coord, 1.5, level="A")) <= 2:
                    continue

            t = time()
            substructure_data_dir = f"{self.data_dir}/sub_{i}"
            system(f"mkdir {substructure_data_dir}")

            atoms_up_to_radius = kdtree.search(atom.coord, radius, level="A")
            selector.full_ids = set([atom.full_id for atom in atoms_up_to_radius])
            io.save(f"{substructure_data_dir}/atoms_up_to_{radius}_angstroms.pdb", selector)

            atoms_up_to_12A = kdtree.search(atom.coord, 12, level="A")
            selector.full_ids = set([atom.full_id for atom in atoms_up_to_12A])
            io.save(f"{substructure_data_dir}/atoms_up_to_12_angstroms.pdb", selector)

            calculated_atoms = kdtree.search(atom.coord, 1.5, level="A") # ať se počítají jen centrální atom, dvouvazné kyslíky a vodíky
            calculated_atoms_full_ids = set([calculated_atom.full_id[1:] for calculated_atom in calculated_atoms])

            mol_6A = Chem.MolFromPDBFile(f"{substructure_data_dir}/atoms_up_to_{radius}_angstroms.pdb", removeHs=False, sanitize=False)
            mol_6A_conformer = mol_6A.GetConformer()
            mol_12A = Chem.MolFromPDBFile(f"{substructure_data_dir}/atoms_up_to_12_angstroms.pdb", removeHs=False, sanitize=False)
            mol_12A_conformer = mol_12A.GetConformer()


            mol_6A_coord_dict = {}
            for aatom in mol_6A.GetAtoms():
                position = mol_6A_conformer.GetAtomPosition(aatom.GetIdx())
                mol_6A_coord_dict[(position.x, position.y, position.z)] = aatom
            mol_12A_coord_dict = {}
            for aatom in mol_12A.GetAtoms():
                position = mol_12A_conformer.GetAtomPosition(aatom.GetIdx())
                mol_12A_coord_dict[(position.x, position.y, position.z)] = aatom

            atoms_with_broken_bonds = []
            for aatom in mol_6A.GetAtoms():
                position = mol_6A_conformer.GetAtomPosition(aatom.GetIdx())
                mol_12A_atom = mol_12A_coord_dict[(position.x, position.y, position.z)]
                if len(aatom.GetNeighbors()) != len(mol_12A_atom.GetNeighbors()):
                    atoms_with_broken_bonds.append(mol_12A_atom)

            carbons_with_broken_bonds_positions = []
            while atoms_with_broken_bonds:
                atom_with_broken_bonds = atoms_with_broken_bonds.pop(0)
                bonded_atoms = atom_with_broken_bonds.GetNeighbors()
                for ba in bonded_atoms:
                    position = mol_12A_conformer.GetAtomPosition(ba.GetIdx())
                    if (position.x, position.y, position.z) in mol_6A_coord_dict:
                        continue
                    else:
                        if atom_with_broken_bonds.GetSymbol() == "C" and ba.GetSymbol() == "C":
                            carbons_with_broken_bonds_positions.append(mol_12A_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                            continue
                        else:
                            atoms_with_broken_bonds.append(ba)
                            mol_6A_coord_dict[(position.x, position.y, position.z)] = ba

            substructure_atoms = []
            for atom in atoms_up_to_12A:
                if tuple(round(float(x),3) for x in atom.coord) in mol_6A_coord_dict:
                    substructure_atoms.append(atom)


            selector.full_ids = set([atom.full_id for atom in substructure_atoms])
            io.save(f"{substructure_data_dir}/substructure.pdb", selector)

            substructure_charge = round(sum([atom.atomic_charge_estimation for atom in substructure_atoms]))
            # log, pokud to není celé číslo!

            system(
             f"cd {substructure_data_dir} ; obabel -iPDB -oPDB substructure.pdb -h > reprotonated_substructure.pdb 2>/dev/null")
            with open(f"{substructure_data_dir}/reprotonated_substructure.pdb") as reprotonated_substructure_file:
                atom_lines = [line for line in reprotonated_substructure_file.readlines() if line[:4] in ["ATOM", "HETA"]]
                original_atoms = atom_lines[:len(substructure_atoms)]
                added_atoms = atom_lines[len(substructure_atoms):]


            with open(f"{substructure_data_dir}/repaired_substructure.pdb", "w") as repaired_substructure_file:
                repaired_substructure_file.write("".join(original_atoms))
                for added_atom in added_atoms:

                    if any([dist([float(added_atom[30:38]), float(added_atom[38:46]), float(added_atom[46:54])], (p.x, p.y, p.z)) < 1.3 for p in carbons_with_broken_bonds_positions]):
                        repaired_substructure_file.write(added_atom)

            system(f"cd {substructure_data_dir} ; "
                   f"xtb repaired_substructure.pdb --gfn 1 --gbsa water --acc 1000 --chrg {substructure_charge}   > xtb_output.txt 2> xtb_error_output.txt ")


            xtb_output_file_lines = open(f"{substructure_data_dir}/xtb_output.txt").readlines()
            charge_headline_index = xtb_output_file_lines.index(
                "  Mulliken/CM5 charges         n(s)   n(p)   n(d)\n")
            for substructure_index, substructure_atom in enumerate(PDB.PDBParser().get_structure("substructure",
                                                                                 f"{substructure_data_dir}/substructure.pdb").get_atoms()):
                if substructure_atom.full_id[1:] in calculated_atoms_full_ids: # todo, nejspíše přepisujeme přesněji vypočítané náboje za horši!
                    charge = float(xtb_output_file_lines[charge_headline_index + substructure_index + 1].split()[3])
                    structure[substructure_atom.full_id[2]][substructure_atom.get_parent().id][substructure_atom.id].cm5_charge = charge

            print(i, time() - t)
            if self.delete_auxiliary_files:
                system(f"rm -r {substructure_data_dir}")

        from os import path
        charges = [atom.cm5_charge for atom in structure.get_atoms()]
        open(f"{self.data_dir}/charges.txt", "w").write(f"{path.basename(self.mmCIF_file)[:-4]}\n" + " ".join([str(x) for x in charges]))


        input_file = f"{self.data_dir}/protonated_protein.cif"
        structure = gemmi.cif.read_file(input_file)
        block = structure.sole_block()
        block.find_mmcif_category('_chem_comp.').erase() # remove pesky _chem_comp category >:(
        sb_ncbr_partial_atomic_charges_meta_prefix = "_sb_ncbr_partial_atomic_charges_meta."
        sb_ncbr_partial_atomic_charges_meta_attributes = ["id",
                                                  "type",
                                                  "method"]
        metadata_loop = block.init_loop(sb_ncbr_partial_atomic_charges_meta_prefix,
                                        sb_ncbr_partial_atomic_charges_meta_attributes)
        metadata_loop.add_row(['1',
                               "'empirical'",
                               "'SQE+qp/Schindler 2021 (PUB_pept)'"])
        sb_ncbr_partial_atomic_charges_prefix = "_sb_ncbr_partial_atomic_charges."
        sb_ncbr_partial_atomic_charges_attributes = ["type_id",
                                             "atom_id",
                                             "charge"]
        charges_loop = block.init_loop(sb_ncbr_partial_atomic_charges_prefix,
                                       sb_ncbr_partial_atomic_charges_attributes)
        for atomId, charge in enumerate(charges):
            charges_loop.add_row(["1",
                                  f"{atomId + 1}",
                                  f"{charge: .4f}"])
        block.write_file(f"{self.data_dir}/final.cif")

        if self.delete_auxiliary_files:
            system(f"rm {self.mmCIF_file};"
                   f"cd {self.data_dir};"
                   f"rm without_hydrogens.cif protonated_protein.cif protonated_ligands.cif fixed.cif")



if __name__ == "__main__":
    args = load_arguments()

    system(f"mkdir {args.data_dir}; "
           f"cp {args.mmCIF_file} {args.data_dir}")

    fix_structure(input_mmCIF_file=f"{args.data_dir}/{path.basename(args.mmCIF_file)}",
                  data_dir=args.data_dir)

    prepared_structure = StructurePreparer(input_mmCIF_file=f"{args.data_dir}/fixed.cif",
                                           data_dir=args.data_dir)

    calculator = ChargesCalculator(mmCIF_file=f"{args.data_dir}/protonated_protein.cif",
                                   atomic_charge_estimations=prepared_structure.atomic_charges,
                                   data_dir=args.data_dir,
                                   delete_auxiliary_files=args.delete_auxiliary_files)
    calculator.calculate_charges()

# dotazy na chlapy
# stahovat vždy jedno pdb a nebo stáhnout celou PDB a ?
# jak cesty k pdb2pqr, hydride, xtb?
# výsledky potřebujeme někam uložit aby si to pak webovka mohla tahat
# uděláme testovací run? Třeba 1000 struktur?
# pořešit, zda to jde na vícekrát?




# možná nepůjde stáhnout všechno
# určitě log pro všechny rezidua
# nechat si vždycky verzi, se kteoru se pracovalo
