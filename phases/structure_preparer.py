from os import system,path
import argparse

from biotite.structure import BondType, BondList
from biotite.structure.io import pdbx as biotite_mmCIF
from .atom_selector import AtomSelector
from .chain_changer import change_chain_names

import hydride
from Bio import PDB

import numpy as np
from biotite.structure.residues import get_residue_starts_for
from dimorphite_dl import DimorphiteDL
from moleculekit import molecule as molit
from moleculekit.tools.preparation import systemPrepare
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


class StructurePreparer:
    def __init__(self,
                 mmCIF_file: str,
                 data_dir: str,
                 delete_auxiliary_files: bool,
                 save_charges_estimation: bool = True):
        self.mmCIF_file = mmCIF_file
        self.data_dir = data_dir
        self.delete_auxiliary_files = delete_auxiliary_files
        self.biotite_mmCIF_file = biotite_mmCIF.CIFFile.read(self.mmCIF_file)
        self.save_charges_estimation = save_charges_estimation
        system(f"mkdir {self.data_dir}")
        self.pH = 7.2

    def fix_structure(self):
        """
        mmCIF file is fixed by tool PDBFixer.
        https://github.com/openmm/pdbfixer

        PDBFixer solves common problems in protein structure files.
        It selects the first model, alternative locations, fills in missing heavy atoms, etc.
        """

        print("Fixing structure... ", end="")
        fixer = PDBFixer(filename=self.mmCIF_file)
        set_of_resnames = set([x.name for x in fixer.topology.residues()])


        # napsat komentář, proč to neděláme!

        # # download templates for heteroresidues
        # fixer_available_resnames = set(fixer.templates.keys())
        # for resname in set_of_resnames:
        #     if resname not in fixer_available_resnames:
        #         try:
        #             fixer.downloadTemplate(resname)
        #         except:
        #             exit(f"ERROR! Old version of pdbfixer installed or heteroresiduum {resname} does not exist!")

        # add heavy atoms
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        # with open(f"{self.data_dir}/logs/added_heavy_atoms.txt", "w") as added_heavy_atoms_file: # log it TODO
        #     for residue, residue_list in fixer.missingAtoms.items():
        #         for atom in residue_list:
        #             added_heavy_atoms_file.write(f"{residue} {atom}\n")
        fixer.addMissingAtoms()

        # write fixed structure to file
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/pdbfixer.cif", 'w'), keepIds=True)

        # PDBFixer removes values from a struct_conn block
        # and therefore only atom_site block is overwritten from PDBFixer to pdbfixer.cif.
        fixed_mmCIF_file = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/pdbfixer.cif")
        self.biotite_mmCIF_file.block["atom_site"] = fixed_mmCIF_file.block["atom_site"]
        self.biotite_mmCIF_file.write(f"{self.data_dir}/pdbfixer.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/pdbfixer.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/pdbfixer.cif", "w").write(repaired_mmcif_string)
        print("ok")

    def remove_hydrogens(self):
        print("Removing hydrogens ... ", end="")
        protein = biotite_mmCIF.get_structure(self.biotite_mmCIF_file,
                                      model=1,
                                      extra_fields=["b_factor", "occupancy", "charge"],
                                      include_bonds=True)
        protein_without_hydrogens = protein[protein.element != "H"]
        biotite_mmCIF.set_structure(self.biotite_mmCIF_file, protein_without_hydrogens)
        self.biotite_mmCIF_file.write(f"{self.data_dir}/without_hydrogens.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/without_hydrogens.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/without_hydrogens.cif", "w").write(repaired_mmcif_string)


        print("ok")


    def add_hydrogens_by_hydride(self):
        """
        This function is based on the biotite, hydride, RDKit and dimorphite_dl libraries.
        https://github.com/biotite-dev/biotite
        https://github.com/biotite-dev/hydride
        https://github.com/rdkit/rdkit
        https://github.com/durrantlab/dimorphite_dl

        Hydrogens are added to the heteroresidues by the hydride library, which is built on top of the biotite library.
        Prior to the adding of hydrogens, the formal charges are loaded from CCD dictionary and extended by the dimorphite_dl.
        Dimorphite_dl formal charges are mapped to CCD molecule by RDKit library.
        """

        print("Adding hydrogens by hydride... ", end="")
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
        # shortcuts for RNA, also processed by pdb2pqr, defined in RNA_MAPPING
        residues_processed_by_pdb2pqr.update(["A", "C", "G", "U"])

        self.label_asym_ids_map, self.auth_asym_ids_map = change_chain_names(f"{self.data_dir}/without_hydrogens.cif",
                                                                             f"{self.data_dir}/without_hydrogens_m.cif")
        asdf = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/without_hydrogens_m.cif")
        self.biotite_mmCIF_file.block["atom_site"] = asdf.block["atom_site"]




        structure = PDB.MMCIFParser(QUIET=True).get_structure("structure", f"{self.data_dir}/without_hydrogens_m.cif")[0]

        structure_atoms = list(structure.get_atoms())
        structure_atoms.sort(key=lambda x: x.serial_number)
        selector = AtomSelector()
        io = PDB.PDBIO()
        io.set_structure(structure)
        kdtree = PDB.NeighborSearch(structure_atoms)
        for atom in structure_atoms:
            atom.charge_estimation = 0
            atom.charged_by_dimorphite = False
            atom.hydride_mask = False

        residues_processed_by_hydride = [res for res in structure.get_residues() if res.resname not in residues_processed_by_pdb2pqr]

        residues_protonated_by_hydride = set() # špatně, popsat rozdíl mezi processed a protonated

        # todo!! je to tak správně? titruje je propka?


        protein = biotite_mmCIF.get_structure(self.biotite_mmCIF_file,
                                              model=1,
                                              extra_fields=["b_factor", "occupancy", "charge"],
                                              include_bonds=True)

        # for biopython_atom, biotite_atom in zip(structure_atoms, protein):
        #     print(biopython_atom.full_id, biotite_atom)
        # exit()

        if residues_processed_by_hydride:
            # load formal charges for ligand from CCD. Add other formal charges by dimorphite
            dimorphite = DimorphiteDL(min_ph=self.pH,
                                      max_ph=self.pH,
                                      max_variants=1,
                                      label_states=False,
                                      pka_precision=0.001)
            residues_processed_by_hydride_names = set([res.resname for res in residues_processed_by_hydride])
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
                        CCD_mol_smiles = Chem.MolToSmiles(CCD_mol)
                        dimorphite_smiles = dimorphite.protonate(CCD_mol_smiles)[0]
                        dimorphite_mol = Chem.MolFromSmiles(dimorphite_smiles)
                        dimorphite_mol = Chem.RemoveAllHs(dimorphite_mol)
                        params = rdFMCS.MCSParameters()
                        params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                        params.BondTyper = rdFMCS.BondCompare.CompareOrder # fakt? Jakto že to funguje? Testnout ješte jednou!
                        params.BondCompareParameters.RingMatchesRingOnly = True
                        params.BondCompareParameters.CompleteRingsOnly = True
                        params.AtomCompareParameters.MatchFormalCharge = False
                        params.Timeout = 60
                        res = rdFMCS.FindMCS([CCD_mol, dimorphite_mol], params)
                        atom_indices_map = [x[1] for x in
                                            sorted(zip(dimorphite_mol.GetSubstructMatch(res.queryMol),
                                                            CCD_mol.GetSubstructMatch(res.queryMol)))]
                        if len(dimorphite_mol.GetAtoms()) != len(atom_indices_map):
                            residues_processed_by_hydride_formal_charges[mol_name] = (CCD_mol,
                                                                                      "Atom mapping of structures from CCD and dimorphite failed. The formal charges are taken from the CCD.")
                        else:
                            dimorphite_formal_charges = [0 for _ in dimorphite_mol.GetAtoms()] # přepsat na generátor
                            for atom, a_i in zip(dimorphite_mol.GetAtoms(), atom_indices_map):
                                dimorphite_formal_charges[a_i] = atom.GetFormalCharge()

                            for CCD_atom in CCD_mol.GetAtoms():
                                CCD_atom.SetProp("ChargedByDimorphite", "0")

                            for CCD_atom, dimorphite_formal_charge in zip(CCD_mol.GetAtoms(),
                                                                          dimorphite_formal_charges):

                                if CCD_atom.GetFormalCharge() == 0 and dimorphite_formal_charge != 0:
                                    bonded_atoms = list(CCD_atom.GetNeighbors())
                                    bonded_atoms_over_two_bonds = []
                                    for bonded_atom in bonded_atoms:
                                        bonded_atoms_over_two_bonds.extend(bonded_atom.GetNeighbors())
                                    if all([atom.GetFormalCharge() == 0 for atom in
                                            bonded_atoms + bonded_atoms_over_two_bonds]):  # CCD_atom is already in bonded_atoms_over_two_bonds
                                        CCD_atom.SetProp("ChargedByDimorphite", "1")
                                        CCD_atom.SetFormalCharge(dimorphite_formal_charge)

                            residues_processed_by_hydride_formal_charges[mol_name] = (CCD_mol,
                                                                                      "The formal charges are taken from CCD and dimorphite.")

            for res in residues_processed_by_hydride:
                residues_protonated_by_hydride.add(res)

                try:
                    CCD_mol, log = residues_processed_by_hydride_formal_charges[res.resname]
                except KeyError:
                    # zalogovat, že ligand není v CCD, necháváme neutrální
                    print("DIVNY ERROR!!!")
                    continue



                # logovat podle logu! 1) není načetnutelný, nešlo dimorphite, šlo dimorphite
                # logovat zda sedí počet atomů!
                # if len(res) != len(formal_charges):
                #     # zalogovat, že ligand má jiný počet atomů než v CCD, necháváme neutrální, protonuje hydride


                # ulož reziduum do substruktury
                res_atoms = [atom for atom in res.get_atoms()]
                res_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in res_atoms])
                io.save(f"{self.data_dir}/{res.resname}_{res.id[1]}.pdb", selector)
                rdkit_mol = Chem.MolFromPDBFile(f"{self.data_dir}/{res.resname}_{res.id[1]}.pdb",
                                                removeHs=False,
                                                sanitize=False)
                rdkit_mol = Chem.RemoveAllHs(rdkit_mol)

                params = rdFMCS.MCSParameters()
                params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                params.BondTyper = rdFMCS.BondCompare.CompareAny  # fakt?
                # params.BondCompareParameters.RingMatchesRingOnly = True # fakt?
                # params.BondCompareParameters.CompleteRingsOnly = True # fakt?
                params.AtomCompareParameters.MatchFormalCharge = False
                params.Timeout = 60
                resss = rdFMCS.FindMCS([rdkit_mol, CCD_mol], params)
                atom_indices_map = {x[0]:x[1] for x in
                                    sorted(zip(rdkit_mol.GetSubstructMatch(resss.queryMol), CCD_mol.GetSubstructMatch(resss.queryMol)))}

                #if len(map) != len(res): log it!, continue

                for atom_i, atom in enumerate(res.get_atoms()):
                    atom.charge_estimation = CCD_mol.GetAtoms()[atom_indices_map[atom_i]].GetFormalCharge()
                    atom.charged_by_dimorphite = bool(int(CCD_mol.GetAtoms()[atom_indices_map[atom_i]].GetProp("ChargedByDimorphite")))
                # zalogovat, že ligandu jsou přiřazeny náboje buď podle CCD a nebo podle dimorphite_dl



                from math import dist
                res_center = res.center_of_mass(geometric=True)
                res_radius = max([dist(res_center, atom.coord) for atom in res.get_atoms()])
                substructure_atoms = kdtree.search(res_center, res_radius + 5, level="A")
                substructure_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in substructure_atoms])
                io.save(f"{self.data_dir}/{res.resname}_{res.id[1]}_substructure.pdb", selector, preserve_atom_numbering=True) # možná změnit res.id[1]

                rdkit_mol = Chem.MolFromPDBFile(f"{self.data_dir}/{res.resname}_{res.id[1]}_substructure.pdb",
                                                removeHs=False,
                                                sanitize=False)

                # přemístit mimo cyklus
                biotite_bonds_set = set([frozenset((a1, a2)) for a1, a2 in protein.bonds.as_array()[:, :2]]) # vyřazujeme třetí sloupec s typem vazby
                rdkit_biotite_bonds_converter = {Chem.BondType.SINGLE: BondType.SINGLE,
                                                 Chem.BondType.DOUBLE: BondType.DOUBLE}

                for bond in rdkit_mol.GetBonds():
                    ba1_serial_number = bond.GetBeginAtom().GetPDBResidueInfo().GetSerialNumber()
                    ba2_serial_number = bond.GetEndAtom().GetPDBResidueInfo().GetSerialNumber()
                    ba1_index = ba1_serial_number - 1
                    ba2_index = ba2_serial_number - 1
                    ba1 = structure_atoms[ba1_index]
                    ba2 = structure_atoms[ba2_index]
                    ba1_res = ba1.get_parent()
                    ba2_res = ba2.get_parent()

                    if ba1_res != ba2_res and res in (ba1_res, ba2_res): # inter-residual bond

                        # nastav nulový náboj pro struktury v peptidové vazbě
                        # toto platí jak pro CCD tak pro dimorphite
                        if set([ba1.element, ba2.element]) == {"N", "C"}:
                            carbon = [atom for atom in [ba1, ba2] if atom.element == "C"][0]
                            if "O" in [atom.element for atom in kdtree.search(carbon.coord, 1.3, level="A")]:
                                ba1.charge_estimation = 0
                                ba2.charge_estimation = 0

                        # nastav nulový náboj pro dimorphite
                        if ba1.charged_by_dimorphite:
                            ba1.charge_estimation = 0
                        if ba2.charged_by_dimorphite:
                            ba2.charge_estimation = 0




                        if frozenset((ba1_index, ba1_index)) not in biotite_bonds_set:
                            biotite_bond_type = rdkit_biotite_bonds_converter.get(bond.GetBondType(), BondType.ANY)
                            protein.bonds.add_bond(ba1_index, ba2_index, biotite_bond_type)
                            print(f"added bond!!! {ba1_index} {ba2_index} {biotite_bond_type}")
                            # log it!

                        residues_protonated_by_hydride.update({ba1_res, ba2_res})

        for res in residues_protonated_by_hydride:
            for atom in res.get_atoms():
                atom.hydride_mask = True


        protein.set_annotation("hydride_mask", [atom.hydride_mask for atom in structure_atoms])


        hydride_estimated_charges = hydride.estimate_amino_acid_charges(protein, ph=self.pH)
        for atom, hydride_estimated_charge in zip(structure_atoms, hydride_estimated_charges):
            if atom.hydride_mask and atom.charge_estimation == 0:
                atom.charge_estimation = hydride_estimated_charge






        protein.charge = [atom.charge_estimation for atom in structure_atoms]


        bond_array = protein.bonds.as_array()
        unknown_order_mask = bond_array[:, 2] == BondType.ANY
        if unknown_order_mask.any():
            print("For some bonds the bond order is unknown, hence single bonds are assumed")
            bond_array[unknown_order_mask, 2] = BondType.SINGLE
        protein.bonds = BondList(protein.bonds.get_atom_count(), bond_array)


        protein_with_hydrogens, _ = hydride.add_hydrogen(protein, mask=protein.hydride_mask)

        self.hydride_charges = protein_with_hydrogens.charge
        self.hydride_mask = protein_with_hydrogens.hydride_mask



        biotite_mmCIF.set_structure(self.biotite_mmCIF_file, protein_with_hydrogens)
        bond_array = protein_with_hydrogens.bonds.as_array()
        residue_starts_1, residue_starts_2 = (get_residue_starts_for(protein_with_hydrogens, bond_array[:, :2].flatten()).reshape(-1, 2).T)
        interresidual_bonds = bond_array[residue_starts_1 != residue_starts_2]
        ptnr1_auth_asym_id = []
        ptnr2_auth_asym_id = []
        ptnr1_auth_seq_id = []
        ptnr2_auth_seq_id = []
        protein_auth_asym_id = self.biotite_mmCIF_file.block["atom_site"]["auth_asym_id"].as_array()
        protein_auth_seq_id = self.biotite_mmCIF_file.block["atom_site"]["auth_seq_id"].as_array()
        for a1, a2, _ in interresidual_bonds:
            ptnr1_auth_asym_id.append(protein_auth_asym_id[a1])
            ptnr2_auth_asym_id.append(protein_auth_asym_id[a2])
            ptnr1_auth_seq_id.append(protein_auth_seq_id[a1])
            ptnr2_auth_seq_id.append(protein_auth_seq_id[a2])
        self.biotite_mmCIF_file.block["struct_conn"]["ptnr1_auth_asym_id"] = ptnr1_auth_asym_id
        self.biotite_mmCIF_file.block["struct_conn"]["ptnr2_auth_asym_id"] = ptnr2_auth_asym_id
        self.biotite_mmCIF_file.block["struct_conn"]["ptnr1_auth_seq_id"] = ptnr1_auth_seq_id
        self.biotite_mmCIF_file.block["struct_conn"]["ptnr2_auth_seq_id"] = ptnr2_auth_seq_id

        self.biotite_mmCIF_file.write(f"{self.data_dir}/hydride.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/hydride.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/hydride.cif", "w").write(repaired_mmcif_string)
        print("ok")


    def add_hydrogens_by_moleculekit(self):
        print("Adding hydrogens by moleculekit... ", end="")

        molecule = molit.Molecule(f"{self.data_dir}/hydride.cif")


        prepared_molecule = systemPrepare(molecule,
                                          pH=self.pH,
                                          hold_nonpeptidic_bonds=False,
                                          ignore_ns_errors=True,
                                          _molkit_ff=False) # todo
        prepared_molecule.write(f"{self.data_dir}/moleculekit_changed_chains.cif")

        moleculekit_mmCIF_file = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/moleculekit_changed_chains.cif")
        hydride_mmCIF_file = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/hydride.cif")
        self.biotite_mmCIF_file.block["atom_site"] = moleculekit_mmCIF_file.block["atom_site"]
        self.biotite_mmCIF_file.block["struct_conn"] = hydride_mmCIF_file.block["struct_conn"]
        self.biotite_mmCIF_file.write(f"{self.data_dir}/moleculekit_changed_chains.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/moleculekit_changed_chains.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/moleculekit_changed_chains.cif", "w").write(repaired_mmcif_string)


        pdb2pqr_charges = np.nan_to_num(prepared_molecule.charge)


        hydride_structure = PDB.MMCIFParser(QUIET=True).get_structure("structure", f"{self.data_dir}/hydride.cif")[0]
        for atom, hydride_charge, hydride_mask in zip(hydride_structure.get_atoms(), self.hydride_charges, self.hydride_mask):
            atom.charge_estimation = hydride_charge
            atom.hydride_mask = hydride_mask



        # structure which combine hydrogens and charges from hydride and moleculekit
        combined_structure = PDB.MMCIFParser(QUIET=True).get_structure("structure", f"{self.data_dir}/moleculekit_changed_chains.cif")[0]
        for atom, moleculekit_charge in zip(combined_structure.get_atoms(), pdb2pqr_charges):
            atom.charge_estimation = moleculekit_charge

        for h_chain, c_chain in zip(sorted(hydride_structure), sorted(combined_structure)):
            for res_i, (h_res, c_res) in enumerate(zip(h_chain, c_chain)):
                if h_res.full_id != c_res.full_id:
                    exit(f"ERROR {h_res.full_id} {c_res.full_id}")
                if any([atom.hydride_mask for atom in h_res]):
                    combined_structure[c_chain.id].detach_child(c_res.id)
                    combined_structure[c_chain.id].insert(res_i, h_res)
                    # log
                else:
                    pass
                    # also log




        io = PDB.MMCIFIO()
        io.set_structure(combined_structure)
        io.save(f"{self.data_dir}/combined.cif", preserve_atom_numbering=True)
        exit()


        change_chain_names(f"{self.data_dir}/combined.cif",
                           f"{self.data_dir}/combined2.cif",
                           self.label_asym_ids_map,
                           self.auth_asym_ids_map)
        exit()

        if all(chg == 0 for chg in pdb2pqr_charges):
            print("Warning! Moleculekit is probably not modified!")
            # Modification of moleculekit:
            # .../moleculekit/tools/preparation.py
            #  line 827 ("ffcharge", "charge")
            # https://github.com/Acellera/mdoleculekit/issues/136

        # print(sum(pdb2pqr_charges)) # jakto že je to takto divné? todo dopsat warning!

        # moleculekit_mmCIF_file = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/moleculekit.cif")
        # self.biotite_mmCIF_file.block["atom_site"] = moleculekit_mmCIF_file.block["atom_site"]
        # if self.save_charges_estimation:
        #     self.biotite_mmCIF_file.block["atom_site"]["charge_estimation"] = [round(charge, 4) for charge in pdb2pqr_charges + prepared_molecule.formalcharge]
        # self.biotite_mmCIF_file.write(f"{self.data_dir}/moleculekit.cif")
        # # these three lines can be probably removed, after biotite 1.0.2 will be released
        # mmcif_string = open(f"{self.data_dir}/moleculekit.cif").read()
        # repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        # open(f"{self.data_dir}/moleculekit.cif", "w").write(repaired_mmcif_string)
        # print("ok")

        moleculekit_mmCIF_file = biotite_mmCIF.CIFFile.read(f"{self.data_dir}/combined.cif")
        self.biotite_mmCIF_file.block["atom_site"] = moleculekit_mmCIF_file.block["atom_site"]
        if self.save_charges_estimation:
            self.biotite_mmCIF_file.block["atom_site"]["charge_estimation"] = [round(charge, 4) for charge in pdb2pqr_charges + prepared_molecule.formalcharge]
        self.biotite_mmCIF_file.write(f"{self.data_dir}/combined.cif")
        # these three lines can be probably removed, after biotite 1.0.2 will be released
        mmcif_string = open(f"{self.data_dir}/combined.cif").read()
        repaired_mmcif_string = mmcif_string.replace("\n# ", "\n# \n")
        open(f"{self.data_dir}/combined.cif", "w").write(repaired_mmcif_string)
        print("ok")