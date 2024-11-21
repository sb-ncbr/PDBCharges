import argparse
from math import dist
from os import system, path

import hydride
import numpy as np
from Bio import PDB as biopython_PDB
from biotite.structure import BondType, BondList
from biotite.structure import io as biotite
from dimorphite_dl import DimorphiteDL
from moleculekit import molecule as moleculekit_PDB
from moleculekit.tools.preparation import systemPrepare as moleculekit_system_prepare
from openmm.app import PDBFile as openmm_PDB
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import rdFMCS


class AtomSelector(biopython_PDB.Select):
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
    return args


class StructurePreparer:
    """
    This class prepares the protein for further research. Specifically, it fixes common problems encountered in PDB files,
    adding hydrogens for specific pH and estimating partial atomic charges.

    The prepared structure is stored in a user-defined data directory as prepared.pdb.
    """

    def __init__(self,
                 PDB_file: str,
                 data_dir: str,
                 delete_auxiliary_files: bool,
                 save_charges_estimation: bool = True):
        """
        :param PDB_file: PDB file containing the structure which should be prepared
        :param data_dir: directory where the results will be stored
        :param delete_auxiliary_files: auxiliary files created during the preraparation will be deleted
        :param save_charges_estimation: TODO
        """
        self.PDB_file = PDB_file
        self.data_dir = data_dir
        self.delete_auxiliary_files = delete_auxiliary_files
        self.save_charges_estimation = save_charges_estimation
        self.pH = 7.2
        system(f"mkdir {self.data_dir}")

    def fix_structure(self):
        """
        PDB file is fixed by tool PDBFixer.
        https://github.com/openmm/pdbfixer

        PDBFixer solves common problems in protein structure files.
        It selects the first model, alternative locations, fills in missing heavy atoms, etc.
        """

        print("Fixing structure... ", end="")
        fixer = PDBFixer(filename=self.PDB_file)

        # download templates for heteroresidues
        set_of_resnames = set([x.name for x in fixer.topology.residues()])
        fixer_available_resnames = set(fixer.templates.keys()) # musí to být set?
        for resname in set_of_resnames:
            if resname not in fixer_available_resnames:
                try:
                    fixer.downloadTemplate(resname)
                    # log downloaded
                except:
                    pass
                    # log heteroresiduum {resname} does not exist!

        # add heavy atoms
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        # with open(f"{self.data_dir}/logs/added_heavy_atoms.txt", "w") as added_heavy_atoms_file: # log it TODO
        #     for residue, residue_list in fixer.missingAtoms.items():
        #         for atom in residue_list:
        #             added_heavy_atoms_file.write(f"{residue} {atom}\n")
        fixer.addMissingAtoms()

        # write fixed structure to file
        openmm_PDB.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/pdbfixer.pdb", 'w'), keepIds=True)
        print("ok")

    def remove_hydrogens(self):
        print("Removing hydrogens ... ", end="")
        protein = biotite.load_structure(file_path=f"{self.data_dir}/pdbfixer.pdb",
                                         model=1,
                                         include_bonds=True)
        protein_without_hydrogens = protein[protein.element != "H"]
        biotite.save_structure(file_path=f"{self.data_dir}/without_hydrogens.pdb",
                               array=protein_without_hydrogens)
        print("ok")

    def add_hydrogens_by_hydride(self):
        """
        This function is based on the Biotite, Hydride, Biopython and Dimorphite-DL libraries.
        https://github.com/biotite-dev/biotite
        https://github.com/biotite-dev/hydride
        https://github.com/biopython/biopython
        https://github.com/rdkit/rdkit
        https://github.com/durrantlab/dimorphite_dl

        In this function, hydrogens are added to all residues to which the MoleculeKit library cannot add hydrogens.
        These are hetero-residues and standard residues that are covalently bound to hetero-residues.

        Hydrogens are added by the hydride library.
        Before adding hydrogens with the hydride tool,
        the formal charges are loaded from CCD dictionary and extended by the dimorphite_dl.
        Dimorphite_dl formal charges are mapped to CCD molecule by RDKit library.
        The RDKit library is also used to search for bonds between hetero-residues and standard residues.
        """

        def get_molecules_from_CCD(molecule_names: list):
            dimorphite = DimorphiteDL(min_ph=self.pH,
                                      max_ph=self.pH,
                                      max_variants=1,
                                      label_states=False,
                                      pka_precision=0.001)
            molecules = {}
            for CCD_mol_sdf in open("CCD/Components-pub.sdf", "r").read().split("$$$$\n"):
                mol_name = CCD_mol_sdf.partition("\n")[0]
                if mol_name in molecule_names:
                    supplier = Chem.SDMolSupplier()
                    supplier.SetData(CCD_mol_sdf)
                    CCD_mol = next(supplier)
                    if CCD_mol is None:
                        molecules[mol_name] = (None,
                                               "The molecule cannot be loaded by RDKit.")
                    else:
                        CCD_mol = Chem.RemoveAllHs(CCD_mol)
                        CCD_mol_smiles = Chem.MolToSmiles(CCD_mol)

                        # add charges to structure by Dimorphite-DL
                        dimorphite_smiles = dimorphite.protonate(CCD_mol_smiles)[0]
                        dimorphite_mol = Chem.MolFromSmiles(dimorphite_smiles)
                        dimorphite_mol = Chem.RemoveAllHs(dimorphite_mol)

                        # Map original mol and mol processed by Dimorphite-DL
                        params = rdFMCS.MCSParameters()
                        params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                        params.BondTyper = rdFMCS.BondCompare.CompareOrder
                        params.BondCompareParameters.RingMatchesRingOnly = True
                        params.BondCompareParameters.CompleteRingsOnly = True
                        params.AtomCompareParameters.MatchFormalCharge = False
                        params.Timeout = 60
                        res = rdFMCS.FindMCS([CCD_mol, dimorphite_mol], params)
                        atom_indices_map = [x[1] for x in sorted(zip(dimorphite_mol.GetSubstructMatch(res.queryMol),
                                                                     CCD_mol.GetSubstructMatch(res.queryMol)))]

                        if len(dimorphite_mol.GetAtoms()) != len(atom_indices_map):
                            molecules[mol_name] = (CCD_mol,
                                                   "Atom mapping of structures from CCD and dimorphite failed. The formal charges are taken from the CCD.")
                        else:
                            dimorphite_formal_charges = [atom.GetFormalCharge() for _, atom in sorted(zip(atom_indices_map,
                                                                                                          dimorphite_mol.GetAtoms()))]
                            for CCD_atom, dimorphite_formal_charge in zip(CCD_mol.GetAtoms(),
                                                                          dimorphite_formal_charges):
                                CCD_atom.SetProp("ChargedByDimorphite", "0")
                                if CCD_atom.GetFormalCharge() == 0 and dimorphite_formal_charge != 0:
                                    bonded_atoms = list(CCD_atom.GetNeighbors())
                                    bonded_atoms_over_two_bonds = []
                                    for bonded_atom in bonded_atoms:
                                        bonded_atoms_over_two_bonds.extend(bonded_atom.GetNeighbors())
                                    # We consider Dimorphite-DL charge only if no atoms across two bonds are charged from CCD charge
                                    if all([atom.GetFormalCharge() == 0 for atom in
                                            bonded_atoms + bonded_atoms_over_two_bonds]):  # CCD_atom is already in bonded_atoms_over_two_bonds
                                        CCD_atom.SetProp("ChargedByDimorphite", "1")
                                        CCD_atom.SetFormalCharge(dimorphite_formal_charge)
                            molecules[mol_name] = (CCD_mol,
                                                   "The formal charges are taken from CCD and dimorphite.")
            return molecules



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

        # load structure by Biopython
        structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                      file=f"{self.data_dir}/without_hydrogens.pdb")[0]
        structure_atoms = list(structure.get_atoms())
        structure_atoms.sort(key=lambda x: x.serial_number)
        selector = AtomSelector()
        io = biopython_PDB.PDBIO()
        io.set_structure(structure)
        kdtree = biopython_PDB.NeighborSearch(structure_atoms)
        for atom in structure_atoms:
            atom.charge_estimation = 0
            atom.charged_by_dimorphite = False
            atom.hydride_mask = False

        # load structure by Biotite
        protein = biotite.load_structure(file_path=f"{self.data_dir}/without_hydrogens.pdb",
                                         model=1,
                                         extra_fields=["charge"],
                                         include_bonds=True)
        biotite_bonds_set = set([frozenset((a1, a2)) for a1, a2 in
                                 protein.bonds.as_array()[:, :2]])  # we exclude the third column with the bond type
        rdkit_biotite_bonds_converter = {Chem.BondType.SINGLE: BondType.SINGLE,
                                         Chem.BondType.DOUBLE: BondType.DOUBLE}

        residues_processed_by_hydride = [res for res in structure.get_residues() if
                                         res.resname not in residues_processed_by_pdb2pqr]

        if residues_processed_by_hydride:
            # load formal charges for ligand from CCD. Add other formal charges by dimorphite
            residues_processed_by_hydride_formal_charges = get_molecules_from_CCD(set(res.resname for res in residues_processed_by_hydride))

            for res in residues_processed_by_hydride:
                res.hydride_mask = True

                # NH2 nabité

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
                params.BondTyper = rdFMCS.BondCompare.CompareAny
                params.AtomCompareParameters.MatchFormalCharge = False
                params.Timeout = 60
                resss = rdFMCS.FindMCS([rdkit_mol, CCD_mol], params)
                atom_indices_map = {x[0]: x[1] for x in
                                    sorted(zip(rdkit_mol.GetSubstructMatch(resss.queryMol),
                                               CCD_mol.GetSubstructMatch(resss.queryMol)))}

                # if len(map) != len(res): log it!, continue

                for atom_i, atom in enumerate(res.get_atoms()):
                    atom.charge_estimation = CCD_mol.GetAtoms()[atom_indices_map[atom_i]].GetFormalCharge()
                    atom.charged_by_dimorphite = bool(int(CCD_mol.GetAtoms()[atom_indices_map[atom_i]].GetProp("ChargedByDimorphite")))
                    # tady udělat kontrolu zda neprohodit náboje na kyslíku na COO
                    # pouze když je to z dimorphite!

                # zalogovat, že ligandu jsou přiřazeny náboje buď podle CCD a nebo podle dimorphite_dl



                # find interrezidual covalent bonds
                res_center = res.center_of_mass(geometric=True)
                res_radius = max([dist(res_center, atom.coord) for atom in res.get_atoms()])
                substructure_atoms = kdtree.search(res_center, res_radius + 5, level="A")
                substructure_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in substructure_atoms])
                io.save(f"{self.data_dir}/{res.resname}_{res.id[1]}_substructure.pdb", selector,
                        preserve_atom_numbering=True)  # možná změnit res.id[1]

                rdkit_mol = Chem.MolFromPDBFile(f"{self.data_dir}/{res.resname}_{res.id[1]}_substructure.pdb",
                                                removeHs=False,
                                                sanitize=False)



                for bond in rdkit_mol.GetBonds():
                    ba1_serial_number = bond.GetBeginAtom().GetPDBResidueInfo().GetSerialNumber()
                    ba2_serial_number = bond.GetEndAtom().GetPDBResidueInfo().GetSerialNumber()
                    ba1_index = ba1_serial_number - 1
                    ba2_index = ba2_serial_number - 1
                    ba1 = structure_atoms[ba1_index]
                    ba2 = structure_atoms[ba2_index]
                    ba1_res = ba1.get_parent()
                    ba2_res = ba2.get_parent()

                    if ba1_res != ba2_res and res in (ba1_res, ba2_res):  # inter-residual bond

                        # nastav nulový náboj pro struktury v peptidové vazbě
                        # toto platí jak pro CCD tak pro dimorphite
                        if set([ba1.element, ba2.element]) == {"N", "C"}:
                            carbon = [atom for atom in [ba1, ba2] if atom.element == "C"][0]

                            # zkontroluj, neměl by to být count?
                            if "O" in [atom.element for atom in kdtree.search(carbon.coord, 1.3, level="A")]:
                                ba1.charge_estimation = 0
                                ba2.charge_estimation = 0
                            # taky O atom by měl být neutrální!

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
                        ba1_res.hydride_mask = True
                        ba2_res.hydride_mask = True

        # final definition which atoms should be processed by hydride
        for res in structure.get_residues():
            if hasattr(res, "hydride_mask"):
                for atom in res.get_atoms():
                    atom.hydride_mask = True

        # set bonds with unknown order as single
        bond_array = protein.bonds.as_array()
        unknown_order_mask = bond_array[:, 2] == BondType.ANY
        if unknown_order_mask.any():
            bond_array[unknown_order_mask, 2] = BondType.SINGLE
        protein.bonds = BondList(protein.bonds.get_atom_count(), bond_array)

        # estimation of charges for standard aminoacids processed by hydride
        hydride_estimated_charges = hydride.estimate_amino_acid_charges(protein, ph=self.pH)
        for atom, hydride_estimated_charge in zip(structure_atoms, hydride_estimated_charges):
            if atom.hydride_mask and atom.charge_estimation == 0:
                atom.charge_estimation = hydride_estimated_charge

        # adding of hydrogens
        protein.charge = [atom.charge_estimation for atom in structure_atoms]
        protein.set_annotation("hydride_mask", [atom.hydride_mask for atom in structure_atoms])
        protein_with_hydrogens, _ = hydride.add_hydrogen(protein, mask=protein.hydride_mask)
        self.hydride_charges = protein_with_hydrogens.charge
        self.hydride_mask = protein_with_hydrogens.hydride_mask
        biotite.save_structure(file_path=f"{self.data_dir}/hydride.pdb",
                               array=protein_with_hydrogens)
        print("ok")

    def add_hydrogens_by_moleculekit(self):
        """
        Hydrogens are added to the protein structure by a library of molecular libraries based on the pdb2pqr tool.
        Moleculekit take into account any non-protein, non-nucleic molecules for the pKa calculation and hydrogen addition.
        However, the Moleculekit is not universal and is only able to work with a limited set of residues.
        """

        print("Adding hydrogens by moleculekit... ", end="")

        molecule = moleculekit_PDB.Molecule(f"{self.data_dir}/hydride.pdb")

        prepared_molecule = moleculekit_system_prepare(molecule,
                                          pH=self.pH,
                                          hold_nonpeptidic_bonds=False,
                                          ignore_ns_errors=True,
                                          _molkit_ff=False)  # todo

        prepared_molecule.write(f"{self.data_dir}/moleculekit.pdb")

        pdb2pqr_charges = np.nan_to_num(prepared_molecule.charge)

        hydride_structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                              file=f"{self.data_dir}/hydride.pdb")[0]
        for atom, hydride_charge, hydride_mask in zip(hydride_structure.get_atoms(), self.hydride_charges,
                                                      self.hydride_mask):
            atom.charge_estimation = hydride_charge
            atom.hydride_mask = hydride_mask

        # structure which combine hydrogens and charges from hydride and moleculekit
        combined_structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                               file=f"{self.data_dir}/moleculekit.pdb")[
            0]
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

        io = biopython_PDB.PDBIO()
        io.set_structure(combined_structure)
        io.save(f"{self.data_dir}/prepared.pdb", preserve_atom_numbering=True)

        if all(chg == 0 for chg in pdb2pqr_charges):
            print("Warning! Moleculekit is probably not modified!")
            # Modification of moleculekit:
            # .../moleculekit/tools/preparation.py
            #  line 827 ("ffcharge", "charge")
            # https://github.com/Acellera/mdoleculekit/issues/136

        # print(sum(pdb2pqr_charges)) # jakto že je to takto divné? todo dopsat warning!



        if self.save_charges_estimation:
            with open(f"{self.data_dir}/estimated_charges.txt", "w") as charges_file:
                charges_string = f"{path.basename(self.PDB_file)[:-4]}\n" + " ".join(
                    [str(round(charge, 4)) for charge in pdb2pqr_charges + prepared_molecule.formalcharge])
                charges_file.write(charges_string)

        print("ok")

        # todo delete auxiliary files!

if __name__ == "__main__":
    args = load_arguments()
    structure_preparer = StructurePreparer(PDB_file=args.PDB_file,
                                           data_dir=args.data_dir,
                                           delete_auxiliary_files=args.delete_auxiliary_files)
    structure_preparer.fix_structure()
    structure_preparer.remove_hydrogens()
    structure_preparer.add_hydrogens_by_hydride()  # add hydrogens to the residues that the moleculekit can't process
    structure_preparer.add_hydrogens_by_moleculekit()  # add hydrogens to rest of structure