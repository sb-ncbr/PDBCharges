import logging
import sys
from math import dist
from os import system

import hydride
import numpy as np
from Bio import PDB as biopython_PDB
from biotite.structure import BondType, BondList
from biotite.structure import io as biotite
from dimorphite_dl import DimorphiteDL
from moleculekit import molecule as moleculekit_PDB
from moleculekit.tools.preparation import systemPrepare as moleculekit_system_prepare, logger
from collections import defaultdict
from openmm.app import PDBFile as openmm_PDB, ForceField
from openmm import NonbondedForce
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

class NucleicSelector(biopython_PDB.Select):
    """
    Support class for Biopython.
    """
    def accept_residue(self, residue):
        return int(biopython_PDB.Polypeptide.is_nucleic(residue))

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
        :param logger: loger of workflow to unify outputs
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
                if CCD_mol is None or mol_name in ["UNX", "UNL"]:
                    molecules[mol_name] = (None,
                                           "The molecule cannot be loaded by RDKit and therefore the residue is left neutral.")
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
                                               "Mapping of formal charges from Dimorphite-DL to CCD failed and therefore formal charges are taken from CCD only.")
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
                                               None)
        return molecules

    def fix_structure(self):
        """
        PDB file is fixed by tool PDBFixer.
        https://github.com/openmm/pdbfixer

        PDBFixer solves common problems in protein structure files.
        It selects the first model, alternative locations, fills in missing heavy atoms, etc.
        """

        self.logger.print("Fixing structure... ", end="")

        # load structure by PDBFixer
        fixer = PDBFixer(filename=self.input_PDB_file)

        # download templates for heteroresidues
        for residue in fixer.topology.residues():
            if residue.name not in fixer.templates.keys():
                try:
                    fixer.downloadTemplate(residue.name)
                except:
                    warning = f"PDBFixer could not download the template for this residue."
                    self.logger.add_warning(chain=residue.chain.id,
                                            resnum=residue.id,
                                            resname=residue.name,
                                            warning=warning)

        # add heavy atoms
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        for residue, missing_atoms in fixer.missingAtoms.items():
            warning = f"Atom(s) {' '.join(atom.name for atom in missing_atoms)} were added by PDBFixer."
            self.logger.add_warning(chain=residue.chain.id,
                                    resnum=residue.id,
                                    resname=residue.name,
                                    warning=warning)
        fixer.addMissingAtoms()
        openmm_PDB.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/pdbfixer.pdb", 'w'), keepIds=True)

        # remove atoms with duplicit names
        # duplicit atoms could be directly in PDB or can be produces as bug in PDBFixer (2ki4)
        structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                      file=f"{self.data_dir}/pdbfixer.pdb")[0]
        io = biopython_PDB.PDBIO()
        io.set_structure(structure)
        io.save(file=f"{self.data_dir}/duplicate_atoms_removed.pdb")

        self.logger.print("ok")

    def remove_hydrogens(self):
        self.logger.print("Removing hydrogens... ", end="")
        protein = biotite.load_structure(file_path=f"{self.data_dir}/duplicate_atoms_removed.pdb",
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
        residues_processed_by_pdb2pqr = {'ALA', 'AR0', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX', 'DA', 'DA3',
                                         'DA5', 'DC', 'DC3', 'DC5', 'DG', 'DG3', 'DG5', 'DT', 'DT3', 'GLH', 'GLN',
                                         'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'HIS', 'HOH', 'HSD', 'HSE', 'HSP', 'ILE',
                                         'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'RA', 'RA3', 'RA5', 'RC', 'RC3',
                                         'RC5', 'RG', 'RG3', 'RG5', 'RU', 'RU3', 'RU5', 'SER', 'THR', 'TRP', 'TYM',
                                         'TYR', 'VAL', 'WAT', "A", "C", "G", "U"}

        # load structure by Biopython
        # we can use serial numbers because biotite removes any inconsistencies during hydrogen removing
        # for example, a structure with PDB code 107d and its serial numbers 218 and 445
        structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                      file=f"{self.data_dir}/without_hydrogens.pdb")[0]
        # Biopython works with atoms hierarchically through chain, residue, atom
        # however, the order of atoms in the file may be different
        # we create a list that is sorted by serial_number and work with it
        structure_atoms = sorted(structure.get_atoms(),
                                 key=lambda x: x.serial_number)
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

        residues_processed_by_hydride = [res for res in structure.get_residues() if res.resname not in residues_processed_by_pdb2pqr]
        # first residues from DNA and RNA chains should be processed by hydride because moleculekit removes atom
        for chain in structure.get_chains():
            first_residue = list(chain.get_residues())[0]
            if biopython_PDB.Polypeptide.is_nucleic(first_residue) and first_residue not in residues_processed_by_hydride:
                residues_processed_by_hydride.append(first_residue)
                for atom in first_residue.get_atoms():
                    atom.hydride_mask = True

        if residues_processed_by_hydride:
            # load formal charges for ligand from CCD. Add other formal charges by Dimorphite-DL
            residues_processed_by_hydride_formal_charges = self._get_molecules_from_CCD(set(res.resname for res in residues_processed_by_hydride))

            for residue in residues_processed_by_hydride:

                # skip unknown residues
                if residue.resname in ["UNX", "UNL"]:
                    continue

                # get residue from CCD (protonated also by Dimorphite-DL)
                residue.hydride_mask = True
                try:
                    CCD_mol, warning = residues_processed_by_hydride_formal_charges[residue.resname]
                    if warning:
                        self.logger.add_warning(chain=residue.get_parent().id,
                                                resnum=residue.id[1],
                                                resname=residue.resname,
                                                warning=warning)
                except KeyError:
                    warning = f"Name {residue.resname} was not found in the CCD and therefore the residue is left neutral."
                    self.logger.add_warning(chain=residue.get_parent().id,
                                            resnum=residue.id[1],
                                            resname=residue.resname,
                                            warning=warning)
                else:
                    # map charges from CCD and Dimorphite-DL to residuum from structure
                    res_atoms = sorted(residue.get_atoms(),
                                       key=lambda x: x.serial_number)
                    selector.full_ids = set([atom.full_id for atom in res_atoms])
                    io.save(file=f"{self.data_dir}/{residue.resname}_{residue.id[1]}.pdb",
                            select=selector)
                    rdkit_mol = Chem.MolFromPDBFile(molFileName=f"{self.data_dir}/{residue.resname}_{residue.id[1]}.pdb",
                                                    removeHs=False,
                                                    sanitize=False)
                    rdkit_mol = Chem.RemoveAllHs(mol=rdkit_mol,
                                                 sanitize=False)
                    params = rdFMCS.MCSParameters()
                    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
                    params.BondTyper = rdFMCS.BondCompare.CompareAny
                    params.AtomCompareParameters.MatchFormalCharge = False
                    params.Timeout = 60
                    MCS_results = rdFMCS.FindMCS([rdkit_mol, CCD_mol], params)
                    atom_indices_map = {x[0]: x[1] for x in
                                        sorted(zip(rdkit_mol.GetSubstructMatch(MCS_results.queryMol),
                                                   CCD_mol.GetSubstructMatch(MCS_results.queryMol)))}

                    if len(atom_indices_map) <= len(residue) - 1: # one oxygen can miss because of peptide bond
                        warning = "Mapping of formal charges from Dimorphite-DL and CCD to residue failed and therefore the residue is left neutral."
                        self.logger.add_warning(chain=residue.get_parent().id,
                                                resnum=residue.id[1],
                                                resname=residue.resname,
                                                warning=warning)
                    else:
                        CCD_mol_atoms = CCD_mol.GetAtoms()
                        for atom_i, atom in enumerate(residue.get_atoms()):
                            try:
                                CCD_atom = CCD_mol_atoms[atom_indices_map[atom_i]]
                                atom.charge_estimation = CCD_atom.GetFormalCharge()
                                atom.charged_by_dimorphite = bool(int(CCD_atom.GetProp("ChargedByDimorphite")))
                            except KeyError: # Mapping for atom failed. It is already logged by previous "if len(atom_indices_map) <= len(res) - 1 statement"
                                continue

                        # because of mapping without bond orders there can be negative charge at double-bond oxygen
                        for atom in residue.get_atoms():
                            if atom.charge_estimation == -1 and atom.element == "O":
                                bonded_atom_indices, bond_types = protein.bonds.get_bonds(atom.serial_number - 1)
                                if 2 in bond_types: # oxygen is bonded by double bond
                                    if len(bonded_atom_indices) > 1:
                                        warning = f"Oxygen {atom.name} with double bond has more neighbors than one."
                                        self.logger.add_warning(chain=residue.get_parent().id,
                                                                resnum=residue.id[1],
                                                                resname=residue.resname,
                                                                warning=warning)
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
                                        warning = f"Oxygen {atom.name} has probably wrong formal charge."
                                        self.logger.add_warning(chain=residue.get_parent().id,
                                                                resnum=residue.id[1],
                                                                resname=residue.resname,
                                                                warning=warning)

                # find interrezidual covalent bonds and modify charge estimation for specific cases
                residue_center = residue.center_of_mass(geometric=True)
                residue_radius = max([dist(residue_center, atom.coord) for atom in residue.get_atoms()])
                substructure_atoms = kdtree.search(center=residue_center,
                                                   radius=residue_radius + 5,
                                                   level="A")
                substructure_atoms.sort(key=lambda x: x.serial_number)
                selector.full_ids = set([atom.full_id for atom in substructure_atoms])
                residuum_file = f"{self.data_dir}/{residue.get_parent().id}_{'_'.join([str(id_part) for id_part in residue.id if id_part != ' '])}_substructure.pdb"
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
                    if ba1_res != ba2_res and residue in (ba1_res, ba2_res):  # inter-residual bond
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
                        if frozenset((ba1_index, ba2_index)) not in biotite_bonds_set:
                            biotite_bond_type = rdkit_biotite_bonds_converter.get(bond.GetBondType(), BondType.ANY)
                            protein.bonds.add_bond(ba1_index, ba2_index, biotite_bond_type)
                        ba1_res.hydride_mask = True
                        ba2_res.hydride_mask = True

                # nitrogens in DNA and RNA should be neutral
                if biopython_PDB.Polypeptide.is_nucleic(residue):
                    for atom in residue.get_atoms():
                        if atom.element == "N":
                            atom.charge_estimation = 0

        # final definition which atoms should be processed by hydride
        # hydrogens should by added to DNA and RNA by moleculekit because of charge consistency
        # (Dimorphite-DL charges nucleic acids differently then moleculekit)
        for residue in structure.get_residues():
            if hasattr(residue, "hydride_mask") and not biopython_PDB.Polypeptide.is_nucleic(residue):
                for atom in residue.get_atoms():
                    atom.hydride_mask = True

        # set bonds with unknown order as single
        bond_array = protein.bonds.as_array()
        unknown_order_mask = bond_array[:, 2] == BondType.ANY
        if unknown_order_mask.any():
            bond_array[unknown_order_mask, 2] = BondType.SINGLE
        protein.bonds = BondList(protein.bonds.get_atom_count(), bond_array)

        # remove metalic bonds
        metal_single_atom_ligands = {'0BE', '3CO', '3NI', '4MO', '4PU', '4TI', '6MO', 'AG', 'AL', 'AM', 'AR', 'ARS',
                                     'AU', 'AU3', 'BA', 'BR', 'BRO', 'BS3', 'CA', 'CD', 'CE', 'CF', 'CL', 'CLO', 'CO',
                                     'CR', 'CS', 'CU', 'CU1', 'CU3', 'DY', 'ER3', 'EU', 'EU3', 'F', 'FE', 'FE2', 'FLO',
                                     'GA', 'GD', 'GD3', 'HG', 'HO', 'HO3', 'IDO', 'IN', 'IOD', 'IR', 'IR3', 'K', 'KR',
                                     'LA', 'LI', 'LU', 'MG', 'MN', 'MN3', 'MO', 'NA', 'ND', 'NI', 'OS', 'OS4', 'PB',
                                     'PD', 'PR', 'PT', 'PT4', 'RB', 'RE', 'RH', 'RH3', 'RHF', 'RU', 'SB', 'SE', 'SM',
                                     'SR', 'TA0', 'TB', 'TE', 'TH', 'TL', 'U1', 'V', 'W', 'XE', 'Y1', 'YB', 'YB2', 'YT3',
                                     'ZCM', 'ZN', 'ZN2', 'ZR', 'ZTM'}
        interresidual_bonds = []
        bonds = protein.bonds.as_array()
        for a1_index, ch1, r1, a2_index, ch2, r2 in zip(bonds[:, 0],
                                                        protein.chain_id[bonds[:, 0]],
                                                        protein.res_id[bonds[:, 0]],
                                                        bonds[:, 1],
                                                        protein.chain_id[bonds[:, 1]],
                                                        protein.res_id[bonds[:,1]]):
            if (ch1, r1) != (ch2, r2):
                interresidual_bonds.append((a1_index, a2_index))
        for a1_index, a2_index in interresidual_bonds:
            if protein.res_name[a1_index] in metal_single_atom_ligands or protein.res_name[a2_index] in metal_single_atom_ligands:
                protein.bonds.remove_bond(a1_index, a2_index)

        # estimation of charges for standard aminoacids processed by hydride
        # the function hydride.estimate_amino_acid_charges has errors, and therefore
        # we control to avoid assigning charges to atoms involved in an interresidual bond
        hydride_estimated_charges = hydride.estimate_amino_acid_charges(protein, ph=self.pH)
        interresidual_bonds_atom_indices = set(atom_i for bond in interresidual_bonds for atom_i in bond)
        for i, (atom, hydride_estimated_charge) in enumerate(zip(structure_atoms, hydride_estimated_charges)):
            if atom.hydride_mask and atom.charge_estimation == 0 and i not in interresidual_bonds_atom_indices:
                atom.charge_estimation = hydride_estimated_charge

        # adding of hydrogens
        protein.charge = [atom.charge_estimation for atom in structure_atoms]
        protein.set_annotation("hydride_mask", [atom.hydride_mask for atom in structure_atoms])
        original_stderr = sys.stderr # redirect hydride output to file
        sys.stderr = open(f"{self.data_dir}/hydride.txt", 'w')
        protein_with_hydrogens, _ = hydride.add_hydrogen(protein, mask=protein.hydride_mask)
        sys.stderr = original_stderr
        self.hydride_charges = protein_with_hydrogens.charge
        self.hydride_mask = protein_with_hydrogens.hydride_mask

        # modify atom names with more symbols than 4
        residue_atom_names = defaultdict(set)
        for atom_name, chain_id, res_id in zip(protein_with_hydrogens.atom_name,
                                               protein_with_hydrogens.chain_id,
                                               protein_with_hydrogens.res_id):
            residue_atom_names[(chain_id, res_id)].add(atom_name)

        for i, (atom_name, chain_id, res_id) in enumerate(zip(protein_with_hydrogens.atom_name,
                                                              protein_with_hydrogens.chain_id,
                                                              protein_with_hydrogens.res_id)):
            if len(atom_name) > 4:
                for x in range(1,1000):
                    candidate_name = f"{atom_name[0]}{x}"
                    if candidate_name not in residue_atom_names[(chain_id, res_id)]:
                        protein_with_hydrogens.atom_name[i] = candidate_name
                        residue_atom_names[(chain_id, res_id)].add(candidate_name)
                        break
        biotite.save_structure(file_path=f"{self.data_dir}/hydride.pdb",
                               array=protein_with_hydrogens)

        # parse warnings from hydride
        hydride_residual_warnings = defaultdict(list)
        warning_lines = [line.strip() for line in open(f"{self.data_dir}/hydride.txt", 'r').readlines() if "Missing fragment for atom" in line]
        for warning_line in warning_lines:
            sl = warning_line.split()
            a_index = int(sl[-1])
            a_name = sl[-4]
            hydride_residual_warnings[(protein.chain_id[a_index],
                                       protein.res_id[a_index],
                                       protein.res_name[a_index])].append(a_name[1:-1]) # remove the quotes
        for (chain, resnum, resname), atoms in hydride_residual_warnings.items():
            warning = f"Hydrogens or formal charges can be added incorrectly because the hydride tool does not have fragment(s) for the atom(s) {' '.join(atoms)}."
            self.logger.add_warning(chain=chain,
                                    resnum=resnum,
                                    resname=resname,
                                    warning=warning)

        self.logger.print("ok")

    def add_hydrogens_by_moleculekit(self):
        """
        Hydrogens are added to the protein structure by a library of molecular libraries based on the pdb2pqr tool.
        Moleculekit take into account any non-protein, non-nucleic molecules for the pKa calculation and hydrogen addition.
        However, the Moleculekit is not universal and is only able to work with a limited set of residues.
        """

        self.logger.print("Adding hydrogens by moleculekit... ", end="")
        try:
            original_stdout = sys.stdout  # redirect moleculekit output to files
            sys.stdout = open(f"{self.data_dir}/moleculekit_chains_report.txt", 'w')
            logger.propagate = False
            file_handler = logging.FileHandler(f"{self.data_dir}/moleculekit_report.txt")
            logger.addHandler(file_handler)
            molecule = moleculekit_PDB.Molecule(f"{self.data_dir}/hydride.pdb")
            prepared_molecule, details = moleculekit_system_prepare(molecule,
                                                                    pH=self.pH,
                                                                    hold_nonpeptidic_bonds=False,
                                                                    ignore_ns_errors=True,
                                                                    _molkit_ff=False,
                                                                    return_details=True)
            prepared_molecule.write(f"{self.data_dir}/moleculekit.pdb")
            sys.stdout = original_stdout
        except:
            sys.stdout = original_stdout
            self.logger.print("\nERROR! The molecule is not processable by the moleculekit library.", end="\n")
            exit()
        self.logger.print("ok")

        self.logger.print("Writing prepared structure to file... ", end="")
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
        for h_chain, c_chain in zip(sorted(hydride_structure),
                                    sorted(combined_structure)):
            for res_i, (h_res, c_res) in enumerate(zip(sorted(h_chain),
                                                       sorted(c_chain))):
                c_res.resname = h_res.resname # rename amber residue names from moleculekit back
                if any([atom.hydride_mask for atom in h_res]):
                    combined_structure[c_chain.id].detach_child(c_res.id)
                    combined_structure[c_chain.id].insert(res_i, h_res)
        for i, atom in enumerate(combined_structure.get_atoms(),
                                 start=1):
            atom.serial_number = i
        io = biopython_PDB.PDBIO()
        io.set_structure(combined_structure)
        try:
            io.save(file=f"{self.data_dir}/combined.pdb",
                    preserve_atom_numbering=True)
        except:
            self.logger.print("\nERROR! The molecule has more then 99999 atoms.", end="\n")
            exit()

        # export pdb to mmcif, because biopython is not so good in such exporting
        protein = biotite.load_structure(f"{self.data_dir}/combined.pdb",
                                         model=1,
                                         extra_fields=["b_factor", "occupancy"])
        biotite.save_structure(f"{self.data_dir}/{self.output_mmCIF_file}", protein)

        if self.save_charges_estimation:
            if any(biopython_PDB.Polypeptide.is_nucleic(residue) for residue in combined_structure.get_residues()):
                # estimate charges also for DNA and RNA
                try:
                    io.save(file=f"{self.data_dir}/only_DNA_and_RNA.pdb",
                            select=NucleicSelector(),
                            preserve_atom_numbering=True)
                    pdb = openmm_PDB(f"{self.data_dir}/only_DNA_and_RNA.pdb")
                    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
                    ff_system = forcefield.createSystem(pdb.topology)
                    nonbonded = [f for f in ff_system.getForces() if isinstance(f, NonbondedForce)][0]
                    charges = [nonbonded.getParticleParameters(i)[0]._value for i in range(ff_system.getNumParticles())]
                    DNARNA_structure = biopython_PDB.PDBParser(QUIET=True).get_structure(id="structure",
                                                                                         file=f"{self.data_dir}/only_DNA_and_RNA.pdb")[0]
                    for atom, charge in zip(sorted(DNARNA_structure.get_atoms(),
                                                   key=lambda x: x.serial_number),
                                            charges):
                        combined_structure[atom.get_parent().get_parent().id][atom.get_parent().id][atom.id].charge_estimation = charge
                except:
                    kdtree = biopython_PDB.NeighborSearch(list(combined_structure.get_atoms()))
                    for residue in combined_structure.get_residues():
                        if biopython_PDB.Polypeptide.is_nucleic(residue):
                            phosphorus_atoms = [atom for atom in residue.get_atoms() if atom.element == "P"]
                            for phosphorus_atom in phosphorus_atoms:
                                neighbors_atoms = kdtree.search(center=phosphorus_atom.coord,
                                                                radius=1.8,
                                                                level="A")
                                neighbors_oxygens = [atom for atom in neighbors_atoms if atom.element == "O"]
                                if len(neighbors_oxygens) == 4 and sum(atom.charge_estimation for atom in neighbors_oxygens) == 0:
                                    phosphorus_atom.charge_estimation = -1

            with open(f"{self.data_dir}/estimated_charges.txt", "w") as charges_file:
                charges_string = " ".join([str(round(atom.charge_estimation, 4)) for atom in combined_structure.get_atoms()])
                charges_file.write(charges_string)

        if self.delete_auxiliary_files:
            system(f"cd {self.data_dir} ; rm *.pdb hydride.txt moleculekit_chains_report.txt moleculekit_report.txt")
        self.logger.print("ok")
