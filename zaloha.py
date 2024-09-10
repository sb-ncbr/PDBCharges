import argparse
from os import path, system
from Bio.PDB import Select, PDBIO, PDBParser, NeighborSearch
from rdkit.Chem import SDMolSupplier, RemoveHs, rdMolAlign, rdmolfiles
from collections import Counter
import biotite.structure as struc
import biotite.structure.io.pdb as biotitePDB
import hydride


def load_arguments():
    print("\nParsing arguments... ", end="")
    parser = argparse.ArgumentParser()
    parser.add_argument("--mmCIF_file", type=str, required=True,
                        help="mmCIF file with protein structure.")
    parser.add_argument("--data_dir", type=str, required=True,
                        help="Directory for saving results.")
    parser.add_argument("--pH", type=float, required=False, default=7.2,
                        help="PDB file with protein structure.")
    args = parser.parse_args()
    if not 0 <= args.pH <= 14:
        exit(f"\nERROR! The pH value must be between 0 and 14!\n")
    if not path.isfile(args.PDB_file):
        exit(f"\nERROR! File {args.PDB_file} does not exist!\n")
    if path.exists(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument!\n")
    print("ok")
    return args


class ChargesCalculator:
    def __init__(self,
                 PDB_file: str,
                 data_dir: str,
                 pH: float):
        self.data_dir = data_dir
        self.pH = pH

        # create data dir and copy PDB file into data dir
        system(f"mkdir {self.data_dir}; "
               f"cp {PDB_file} {self.data_dir}")

        self.PDB_file = f"{self.data_dir}/{path.basename(PDB_file)}"

    def protonate_ligands(self):
        """
        Ligands are compared with their eqvivalents in the CCD dictionary.
        If the element counts are identical, the ligand is left unchanged.
        If the number of heavy atoms is the same but the number of hydrogens is different,
            the formal charges from the CCD dictionary are read and the ligand is protonated.
        If the number of heavy atoms is different, the ligand is protonated with the hydride tool
            without loading the formal charges from the CCD dictionary.
        Ligands are protonated by the hydride tool (https://hydride.biotite-python.org/api.html).
        CCD dictionary is downloaded from http://ligand-expo.rcsb.org/ld-download.html.
        Waters are in this function ignored.
        """

        class HeteroresidueSelector(Select):
            def accept_residue(self, residue):
                if residue == self.heteroresidue:
                    return 1
                else:
                    return 0

        heteroresidue_selector = HeteroresidueSelector()

        # get all heteroresidues in structure (except water)
        structure_heteroresidues = [res for res in self.structure.get_residues() if res.id[0][0] == "H"]

        if structure_heteroresidues:
            # find ideal eqvivalents of heteroresidues in the CCD dictionary
            structure_heteroresidue_names = set([heteroresidue.id[0][2:] for heteroresidue in structure_heteroresidues])
            ideal_heteroresidues = {}
            for ideal_heteroresidue in open("CCD/Components-pub.sdf", "r").read().split("$$$$\n"):
                name = ideal_heteroresidue.partition("\n")[0]
                if name in structure_heteroresidue_names:
                    supplier = SDMolSupplier()
                    supplier.SetData(ideal_heteroresidue, removeHs=False)
                    rdkit_mol = next(supplier)
                    if rdkit_mol is not None:
                        ideal_heteroresidues[name] = rdkit_mol

            # compare heteroresidues with its eqvivalents in CCD
            for heteroresidue in structure_heteroresidues:
                if heteroresidue.resname not in ideal_heteroresidues.keys():
                    # zalogovat, že ligand buď není v CCD a nebo není načetnutelný
                    compare_num_of_Hs = False
                    any_atom_charged = False
                else:
                    ideal_heteroresidue = ideal_heteroresidues[heteroresidue.resname]
                    heteroresidue_elements_counter = Counter([atom.element for atom in heteroresidue.get_atoms()])
                    ideal_heteroresidue_elements_counter = Counter([atom.GetSymbol() for atom in ideal_heteroresidue.GetAtoms()])

                    if heteroresidue_elements_counter == ideal_heteroresidue_elements_counter:
                        continue
                    else:
                        del heteroresidue_elements_counter["H"]
                        ideal_heteroresidue_num_of_Hs = ideal_heteroresidue_elements_counter["H"]
                        del ideal_heteroresidue_elements_counter["H"]

                        if heteroresidue_elements_counter == ideal_heteroresidue_elements_counter: # compare element counters with only heavy atoms
                            compare_num_of_Hs = True
                            ideal_heteroresidue_without_Hs = RemoveHs(ideal_heteroresidue)
                            ideal_heteroresidue_without_Hs_formal_charges = [atom.GetFormalCharge() for atom in ideal_heteroresidue_without_Hs.GetAtoms()]
                            if any([formal_charge != 0 for formal_charge in ideal_heteroresidue_without_Hs_formal_charges]):
                                any_atom_charged = True
                            else:
                                any_atom_charged = False
                        else:
                            # zalogovat, že má jiný počet těžkých atomů než v CCD
                            compare_num_of_Hs = False
                            any_atom_charged = False

                heteroresidue_selector.heteroresidue = heteroresidue
                self.io.save(f"{self.data_dir}/ligand.pdb", heteroresidue_selector)

                biotite_pdb_file = biotitePDB.PDBFile.read(f"{self.data_dir}/ligand.pdb")
                model = biotitePDB.get_structure(biotite_pdb_file,
                                                 model=1,
                                                 extra_fields=["charge"],
                                                 include_bonds=True)
                bond_array = model.bonds.as_array()
                unknown_order_mask = bond_array[:, 2] == struc.BondType.ANY
                if unknown_order_mask.any():
                    bond_array[unknown_order_mask, 2] = struc.BondType.SINGLE
                    model.bonds = struc.BondList(model.array_length(), bond_array)
                model = model[model.element != "H"]
                if any_atom_charged:
                    rdkit_heteroresidue_without_Hs = rdmolfiles.MolFromPDBFile(f"{self.data_dir}/ligand.pdb")
                    align = rdMolAlign.GetCrippenO3A(ideal_heteroresidue_without_Hs, rdkit_heteroresidue_without_Hs)
                    for match in align.Matches():
                        ideal_heteroresidue_atom, heteroresidue_atom = match
                        model.charge[heteroresidue_atom] = ideal_heteroresidue_without_Hs_formal_charges[ideal_heteroresidue_atom]
                molecule, mask = hydride.add_hydrogen(model)
                pdb_file = biotitePDB.PDBFile()
                biotitePDB.set_structure(pdb_file, molecule)
                pdb_file.write(f"{self.data_dir}/protonated_ligand.pdb")

                # napsat funkci, která přejmenuje duplicitní atomy

                heteroresidue_chain = heteroresidue.full_id[2]
                self.structure[heteroresidue_chain].detach_child(heteroresidue.id)
                self.structure[heteroresidue_chain].add(list(PDBParser(QUIET=True).get_structure("ligand", f"{self.data_dir}/protonated_ligand.pdb").get_residues())[0])
                if compare_num_of_Hs:
                    if len([atom.element for atom in self.structure[heteroresidue_chain][heteroresidue.id].get_atoms() if atom.element == "H"]) != ideal_heteroresidue_num_of_Hs:
                        print(heteroresidue.id, len([atom.element for atom in self.structure[heteroresidue_chain][heteroresidue.id].get_atoms() if atom.element == "H"]), ideal_heteroresidue_num_of_Hs)
                        print("error")
                        print(ideal_heteroresidue_without_Hs_formal_charges)
                        # log, že není stejný počet vodíků!
                        self.io.set_structure(self.structure)
                        self.io.save("pokus.pdb")

                        pass

                        exit()
        self.io.set_structure(self.structure)
        self.io.save("pokus.pdb")








    def calculate_charges(self):

        print("Fixing structure by PDBFixer... ", end="")
        system(f"pdbfixer {self.PDB_file} --add-atoms=heavy --output {self.data_dir}/fixed.pdb")
        print("ok")

        self.structure = PDBParser(QUIET=True).get_structure("structure", f"{self.data_dir}/fixed.pdb")[0]
        self.io = PDBIO()
        self.io.set_structure(self.structure)

        print("Adding missing hydrogens to ligands... ", end="")
        self.protonate_ligands()
        print("ok")







        # postupovat budeme per ligand
        # pokud je počet vodíků stejný jako v CCD, tak vodíky necháme a opíšeme formální náboje
        # pokud je počet vodíků jiný, tak dosadíme formální náboje z CCD a doplníme vodíky pomocí hydride

        # sdf from string




        # keep hydrogens
        # protonate protein by pdb2pqr
        # system(f"pdb2pqr30 --titration-state-method propka "
        #        f"--with-ph {self.pH} --pdb-output {self.data_dir}/protein_and_waters_protonated.pdb {self.data_dir}/fixed.pdb "
        #        f"{self.data_dir}/protein_and_waters_protonated.pqr > {self.data_dir}/propka.log 2>&1 ")


        # sežen formální náboje
            # z pqr souboru - jak budeme přiřazovat?
            # CCD dictionary

        # oprav znovu pdbfixerem
        # bude nutné?



if __name__ == "__main__":
    args = load_arguments()
    ChargesCalculator(args.PDB_file, args.data_dir, args.pH).calculate_charges()


# dotazy na Adriana
# stahovat vždy jedno pdb a nebo stáhnout celou PDB a ?
# jak cesty k pdb2pqr, hydride, xtb?
# výsledky potřebujeme někam uložit aby si to pak webovka mohla tahat
# uděláme testovací run? Třeba 100 struktur?
# pořešit, zda to jde na vícekrát?


# dotazy na Radkou
# V jakých formátech chceme náboje poskytovat? Jako v AlphaCharges?
# Chceme počítat náboje i pro vody?
# jen první model pokud jich je více?
# protonace ligandů kolem kovů? podle délky vazby? vandervalls poloměr? Asi Karel

# belbínův test osobnosti



# možná nepůjde stáhnout všechno
# určitě log pro všechny rezidua
# nechat si vždycky verzi, se kteoru se pracovalo

#todo - volba smazání tmp souborů.







