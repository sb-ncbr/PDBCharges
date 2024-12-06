import logging
import sys
from math import dist
from os import system, path

import hydride
import numpy as np
from Bio import PDB as biopython_PDB
from biotite.structure import BondType, BondList
from biotite.structure import io as biotite
from dimorphite_dl import DimorphiteDL
from moleculekit import molecule as moleculekit_PDB
from moleculekit.tools.preparation import systemPrepare as moleculekit_system_prepare, logger
from openmm.app import PDBFile as openmm_PDB
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import rdFMCS


class AtomSelector(biopython_PDB.Select):
    """
    Support class for Biopython.
    After initialization, a set with all full ids of the atoms to be written into the substructure must be stored in self.full_ids.
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)


class StructurePreparer:
    """
    This class prepares the protein for further research. Specifically, it fixes common problems encountered in PDB files,
    adding hydrogens for specific pH and estimating partial atomic charges.

    The prepared structure is stored in a user-defined data directory in mmCIF format.
    """

    def __init__(self,
                 input_PDB_file: str,
                 CCD_file: str,
                 logger,
                 data_dir: str,
                 output_mmCIF_file: str,
                 delete_auxiliary_files: bool,
                 save_charges_estimation: bool = False):
        """
        :param input_PDB_file: PDB file containing the structure which should be prepared
        :param CCD_file: SDF file with Chemical Component Dictionary
        :param data_dir: directory where the results will be stored
        :param output_mmCIF_file: mmCIF file in which prepared structure will be stored
        :param delete_auxiliary_files: auxiliary files created during the preraparation will be deleted
        :param save_charges_estimation: save estimation of partial atomic charges from pdb2pqr, Dimorphite-DL and CCD
        """
        self.logger = logger
        self.logger.print("\nSTRUCTURE PREPARER")
        self.logger.print("Structure preparer initialization... ", end="")
        self.input_PDB_file = input_PDB_file
        self.CCD_file = CCD_file
        self.output_mmCIF_file = output_mmCIF_file
        self.data_dir = data_dir
        system(f"mkdir {self.data_dir}")
        self.delete_auxiliary_files = delete_auxiliary_files
        self.save_charges_estimation = save_charges_estimation
        self.pH = 7.2
        self.logger.print("ok")


    def _get_molecules_from_CCD(self,
                                molecule_names: list):
        """
        The function retrieves molecules and their formal charges from the CCD dictionary
        and adds additional formal charges to them using the Dimorphite-DL library.
        """

        dimorphite = DimorphiteDL(min_ph=self.pH,
                                  max_ph=self.pH,
                                  max_variants=1,
                                  label_states=False,
                                  pka_precision=0.001)
        molecules = {}
        for CCD_mol_sdf in open(self.CCD_file, "r").read().split("$$$$\n"):
            mol_name = CCD_mol_sdf.partition("\n")[0]
            if mol_name in molecule_names: # we process only molecules defined in molecule names
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
                    MCS_results = rdFMCS.FindMCS([CCD_mol, dimorphite_mol], params)
                    atom_indices_map = [x[1] for x in sorted(zip(dimorphite_mol.GetSubstructMatch(MCS_results.queryMol),
                                                                 CCD_mol.GetSubstructMatch(MCS_results.queryMol)))]

                    if len(dimorphite_mol.GetAtoms()) != len(atom_indices_map):
                        molecules[mol_name] = (CCD_mol,
                                               "Atom mapping of structures from CCD and Dimorphite-DL failed. The formal charges are taken from the CCD.")
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
                                # CCD_atom is already in bonded_atoms_over_two_bonds
                                if all([atom.GetFormalCharge() == 0 for atom in bonded_atoms + bonded_atoms_over_two_bonds]):
                                    CCD_atom.SetProp("ChargedByDimorphite", "1")
                                    CCD_atom.SetFormalCharge(dimorphite_formal_charge)
                        molecules[mol_name] = (CCD_mol,
                                               "The formal charges are taken from CCD and Dimorphite-DL.")
        return molecules

    def fix_structure(self):
        """
        PDB file is fixed by tool PDBFixer.
        https://github.com/openmm/pdbfixer

        PDBFixer solves common problems in protein structure files.
        It selects the first model, alternative locations, fills in missing heavy atoms, etc.
        """

        self.logger.print("Fixing structure... ", end="")
        fixer = PDBFixer(filename=self.input_PDB_file)

        # download templates for heteroresidues
        set_of_resnames = set([x.name for x in fixer.topology.residues()])
        for resname in set_of_resnames:
            if resname not in fixer.templates.keys():
                try:
                    fixer.downloadTemplate(resname)
                except:
                    pass
                    # log heteroresiduum {resname} does not exist or cannot be downloaded! todo

        # add heavy atoms
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        for residue, missing_atoms in fixer.missingAtoms.items():
            warning = f"atom(s) {' '.join(atom.name for atom in missing_atoms)} added by pdbfixer"
            self.logger.add_warning(chain=residue.chain.id,
                                    resnum=residue.id,
                                    resname=residue.name,
                                    warning=warning)
        fixer.addMissingAtoms()
        openmm_PDB.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/pdbfixer.pdb", 'w'), keepIds=True)
        self.logger.print("ok")

    def remove_hydrogens(self):
        self.logger.print("Removing hydrogens ... ", end="")
        protein = biotite.load_structure(file_path=f"{self.data_dir}/pdbfixer.pdb",
                                         model=1,
                                         include_bonds=True)
        protein_without_hydrogens = protein[protein.element != "H"]
        biotite.save_structure(file_path=f"{self.data_dir}/without_hydrogens.pdb",
                               array=protein_without_hydrogens)
        self.logger.print("ok")

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
        the formal charges are loaded from CCD dictionary and extended by the Dimorphite-DL.
        Dimorphite-DL formal charges are mapped to CCD molecule by RDKit library.
        The RDKit library is also used to search for bonds between hetero-residues and standard residues.
        """

        self.logger.print("Adding hydrogens by hydride... ", end="")
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
             'DG5', 'RG3', 'RG5', 'DT3', 'RU3', 'RU5']) # todo
        # shortcuts for RNA, also processed by pdb2pqr, defined in RNA_MAPPING
        residues_processed_by_pdb2pqr.update(["A", "C", "G", "U"])

        # load structure by Biopython
        # we can use serial numbers because biotite removes any inconsistencies during hydrogen removing
        # for example, a structure with PDB code 107d and its serial numbers 218 and 445
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
            # load formal charges for ligand from CCD. Add other formal charges by Dimorphite-DL
            residues_processed_by_hydride_formal_charges = self._get_molecules_from_CCD(set(res.resname for res in residues_processed_by_hydride))

            for res in residues_processed_by_hydride:
                res.hydride_mask = True
                try:
                    CCD_mol, log = residues_processed_by_hydride_formal_charges[res.resname]
                    # logovat podle logu! 1) je v ccd, ale není načetnutelný 2) nešlo dimorphite
                except KeyError:
                    # zalogovat, že reziduum není v CCD, necháváme neutrální
                    continue

                # map charges from CCD and Dimorphite-DL to residuum from structure
                res_atoms = [atom for atom in res.get_atoms()]
                res_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in res_atoms])
                io.save(file=f"{self.data_dir}/{res.resname}_{res.id[1]}.pdb",
                        select=selector)
                rdkit_mol = Chem.MolFromPDBFile(molFileName=f"{self.data_dir}/{res.resname}_{res.id[1]}.pdb",
                                                removeHs=False,
                                                sanitize=False)
                rdkit_mol = Chem.RemoveAllHs(rdkit_mol)
                params = rdFMCS.MCSParameters()
                params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                params.BondTyper = rdFMCS.BondCompare.CompareAny
                params.AtomCompareParameters.MatchFormalCharge = False
                params.Timeout = 60
                MCS_results = rdFMCS.FindMCS([rdkit_mol, CCD_mol], params)
                atom_indices_map = {x[0]: x[1] for x in
                                    sorted(zip(rdkit_mol.GetSubstructMatch(MCS_results.queryMol),
                                               CCD_mol.GetSubstructMatch(MCS_results.queryMol)))}
                if len(atom_indices_map) <= len(res) - 1:
                    print(f"Warning! {res}")
                # může se lišit o jeden kyslík, pak pravděpodobně v peptidové vazbě -> list?
                # pokud se liší o více, tak warning, něco je špatně!
                CCD_mol_atoms = CCD_mol.GetAtoms()
                for atom_i, atom in enumerate(res.get_atoms()):
                    try:
                        CCD_atom = CCD_mol_atoms[atom_indices_map[atom_i]]
                        atom.charge_estimation = CCD_atom.GetFormalCharge()
                        atom.charged_by_dimorphite = bool(int(CCD_atom.GetProp("ChargedByDimorphite")))
                    except KeyError: # Mapping for atom failed. It is already logged by previous "if len(atom_indices_map) <= len(res) - 1 statement"
                        continue

                # because of mapping without bond orders there can be negative charge at double-bond oxygen
                for atom in res.get_atoms():
                    if atom.charge_estimation == -1 and atom.element == "O":
                        bonded_atom_indices, bond_types = protein.bonds.get_bonds(atom.serial_number - 1)
                        if 2 in bond_types: # oxygen is bonded by double bond
                            if len(bonded_atom_indices) > 1:
                                # log
                                continue
                            bonded_atom_2_indices, bond_2_types = protein.bonds.get_bonds(bonded_atom_indices[0]) # bonded atoms over two bonds
                            for bonded_atom_2_index, bond_type in zip(bonded_atom_2_indices, bond_2_types):
                                if bonded_atom_2_index == atom.serial_number - 1:
                                    continue
                                elif protein.element[bonded_atom_2_index] == "O" and bond_type == 1:
                                    right_oxygen = structure_atoms[bonded_atom_2_index]
                                    atom.charge_estimation, right_oxygen.charge_estimation = right_oxygen.charge_estimation, atom.charge_estimation
                                    atom.charged_by_dimorphite, right_oxygen.charged_by_dimorphite = right_oxygen.charged_by_dimorphite, atom.charged_by_dimorphite
                                    break
                            else:
                                # log
                                pass

                # find interrezidual covalent bonds and modify charge estimation for specific cases
                res_center = res.center_of_mass(geometric=True)
                res_radius = max([dist(res_center, atom.coord) for atom in res.get_atoms()])
                substructure_atoms = kdtree.search(center=res_center,
                                                   radius=res_radius + 5,
                                                   level="A")
                substructure_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in substructure_atoms])
                residuum_file = f"{self.data_dir}/{res.get_parent().id}_{'_'.join([str(id_part) for id_part in res.id if id_part != ' '])}_substructure.pdb"
                io.save(file=residuum_file,
                        select=selector,
                        preserve_atom_numbering=True)
                rdkit_mol = Chem.MolFromPDBFile(molFileName=residuum_file,
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
                        # set zero charge for interresidual peptide bonds
                        # this is true for both CCD and Dimorphite-DL formal charges
                        if set([ba1.element, ba2.element]) == {"N", "C"}:
                            carbon = [atom for atom in [ba1, ba2] if atom.element == "C"][0]
                            bonded_oxygens_to_carbon = [atom for atom in kdtree.search(center=carbon.coord,
                                                                                       radius=1.3,
                                                                                       level="A") if atom.element == "O"]
                            if len(bonded_oxygens_to_carbon) == 1:
                                ba1.charge_estimation = 0
                                ba2.charge_estimation = 0
                                bonded_oxygens_to_carbon[0].charge_estimation = 0
                        # setting formal charge from Dimorphite-DL to zero for all interresidual bonds
                        if ba1.charged_by_dimorphite:
                            ba1.charge_estimation = 0
                        if ba2.charged_by_dimorphite:
                            ba2.charge_estimation = 0
                        # if the bond was detected by RDKit and not Biotite, create it
                        if frozenset((ba1_index, ba1_index)) not in biotite_bonds_set:
                            biotite_bond_type = rdkit_biotite_bonds_converter.get(bond.GetBondType(), BondType.ANY)
                            protein.bonds.add_bond(ba1_index, ba2_index, biotite_bond_type)
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
        self.logger.print("ok")

    def add_hydrogens_by_moleculekit(self):
        """
        Hydrogens are added to the protein structure by a library of molecular libraries based on the pdb2pqr tool.
        Moleculekit take into account any non-protein, non-nucleic molecules for the pKa calculation and hydrogen addition.
        However, the Moleculekit is not universal and is only able to work with a limited set of residues.
        """

        self.logger.print("Adding hydrogens by moleculekit... ", end="")
        molecule = moleculekit_PDB.Molecule(f"{self.data_dir}/hydride.pdb")
        original_stdout = sys.stdout # redirect moleculekit output to files
        sys.stdout = open(f"{self.data_dir}/moleculekit_chains_report.txt", 'w')
        logger.propagate = False
        file_handler = logging.FileHandler(f"{self.data_dir}/moleculekit_report.txt")
        logger.addHandler(file_handler)
        prepared_molecule, details = moleculekit_system_prepare(molecule,
                                                                pH=self.pH,
                                                                hold_nonpeptidic_bonds=False,
                                                                ignore_ns_errors=True,
                                                                _molkit_ff=False,
                                                                return_details=True)
        sys.stdout = original_stdout
        prepared_molecule.write(f"{self.data_dir}/moleculekit.pdb")

        # combine structures from hydride and moleculekit
        pdb2pqr_charges = np.nan_to_num(prepared_molecule.charge)
        hydride_structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                              file=f"{self.data_dir}/hydride.pdb")[0]
        for atom, hydride_charge, hydride_mask in zip(hydride_structure.get_atoms(), self.hydride_charges,
                                                      self.hydride_mask):
            atom.charge_estimation = hydride_charge
            atom.hydride_mask = hydride_mask
        # structure which combine hydrogens and charges from hydride and moleculekit
        combined_structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                               file=f"{self.data_dir}/moleculekit.pdb")[0]
        for atom, moleculekit_charge in zip(combined_structure.get_atoms(), pdb2pqr_charges):
            atom.charge_estimation = moleculekit_charge
        for h_chain, c_chain in zip(sorted(hydride_structure), sorted(combined_structure)):
            for res_i, (h_res, c_res) in enumerate(zip(h_chain, c_chain)):
                c_res.resname = h_res.resname # rename amber residue names from moleculekit back
                if any([atom.hydride_mask for atom in h_res]):
                    combined_structure[c_chain.id].detach_child(c_res.id)
                    combined_structure[c_chain.id].insert(res_i, h_res)
        for i, atom in enumerate(combined_structure.get_atoms(),
                                 start=1):
            atom.serial_number = i
        io = biopython_PDB.PDBIO()
        io.set_structure(combined_structure)
        io.save(file=f"{self.data_dir}/combined.pdb",
                preserve_atom_numbering=True)

        # export pdb to mmcif, because biopython is not so good in such exporting
        protein = biotite.load_structure(f"{self.data_dir}/combined.pdb",
                                         model=1,
                                         extra_fields=["b_factor", "occupancy"],
                                         include_bonds=True)
        biotite.save_structure(f"{self.data_dir}/{self.output_mmCIF_file}", protein)

        if self.save_charges_estimation:
            with open(f"{self.data_dir}/estimated_charges.txt", "w") as charges_file:
                charges_string = " ".join([str(round(charge, 4)) for charge in pdb2pqr_charges + prepared_molecule.formalcharge])
                charges_file.write(charges_string)

        if self.delete_auxiliary_files:
            system(f"cd {self.data_dir} ; rm *.txt *.pdb")
        self.logger.print("ok")
