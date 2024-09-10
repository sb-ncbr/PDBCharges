import argparse
from os import path, system
from Bio import PDB
from rdkit.Chem import SDMolSupplier, RemoveHs, rdMolAlign, rdmolfiles
from collections import Counter
import biotite.structure.io.pdbx as biotite_mmCIF
import biotite.structure as biotite_structure
from moleculekit import molecule as molit
from moleculekit.tools.preparation import systemPrepare
import math

import hydride
from pdbfixer import PDBFixer
from openmm.app import PDBxFile



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
    if not path.isfile(args.mmCIF_file):
        exit(f"\nERROR! File {args.mmCIF_file} does not exist!\n")
    if path.exists(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument!\n")
    print("ok")
    return args


class SelectResidues(PDB.Select):
    def accept_residue(self, residue):
        if residue.full_id in self.full_ids:
            return 1
        else:
            return 0

class SelectAtoms(PDB.Select):
    def accept_atom(self, atom):
        if atom.full_id in self.full_ids:
            return 1
        else:
            return 0

class ChargesCalculator:
    def __init__(self,
                 mmCIF_file: str,
                 data_dir: str,
                 pH: float):
        self.data_dir = data_dir
        self.pH = pH
        system(f"mkdir {self.data_dir}; "
               f"cp {mmCIF_file} {self.data_dir}")
        self.mmCIF_file = f"{self.data_dir}/{path.basename(mmCIF_file)}"


    def calculate_charges(self):
        # self.fix_structure() # asi budeme mazat
        self.protonate_ligands()
        self.protonate_protein()
        exit()


        # structure = PDB.MMCIFParser(QUIET=True).get_structure("structure", f"{self.data_dir}/protonated_protein.cif")[0]
        structure = PDB.MMCIFParser(QUIET=True).get_structure("pdb", self.mmCIF_file)[0] # todo change to protonated_protein.cif

        residues = list(structure.get_residues())
        kdtree = PDB.NeighborSearch(list(structure.get_atoms()))

        # nearest_residues = []
        # for residue in residues:
        #     center_of_mass = residue.center_of_mass(geometric=True)
        #     if residue.resname in amk_max_radius:
        #         radius = amk_max_radius[residue.resname]
        #     else:
        #         radius = max([math.dist(atom.coord, center_of_mass) for atom in residue.get_atoms()])
        #     nearest_residues.append(set(kdtree.search(center_of_mass, radius+8, level="R")))

        # from time import time
        #
        # for i, calculated_residue in enumerate(residues, start=1):
        #     t = time()
        #     substructure_data_dir = f"{self.data_dir}/sub_{i}"
        #     system(f"mkdir {substructure_data_dir}")
        #     center_of_mass = calculated_residue.center_of_mass(geometric=True)
        #     radius = max([math.dist(atom.coord, center_of_mass) for atom in calculated_residue.get_atoms()]) + 8
        #     nearest_atoms = kdtree.search(center_of_mass, radius, level="A")
        #     selector = SelectAtoms()
        #     selector.full_ids = set([atom.full_id for atom in nearest_atoms])
        #     io = PDB.PDBIO()
        #     io.set_structure(structure)
        #     io.save(f"{substructure_data_dir}/substructure.pdb", selector)
        #     system(f"cd {substructure_data_dir} ; "
        #            f"xtb substructure.pdb --gfn 1 --gbsa water --acc 1000 > xtb_output.txt")
        #     print(i, len(nearest_atoms), time() - t)

        from time import time
        charges = []
        for i, atom in enumerate(structure.get_atoms()):
            t = time()
            substructure_data_dir = f"{self.data_dir}/sub_{i}"
            system(f"mkdir {substructure_data_dir}")
            selector = SelectAtoms()
            nearest_atoms = kdtree.search(atom.coord, 5, level="A")
            selector.full_ids = set([atom.full_id for atom in nearest_atoms])
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(f"{substructure_data_dir}/substructure.pdb", selector)
            system(f"cd {substructure_data_dir} ; "
                   f"xtb substructure.pdb --gfn 1 --gbsa water --acc 1000 > xtb_output.txt")

            xtb_output_file_lines = open(f"{substructure_data_dir}/xtb_output.txt").readlines()
            charge_headline_index = xtb_output_file_lines.index("  Mulliken/CM5 charges         n(s)   n(p)   n(d)\n")
            for si, substructure_atom in enumerate(PDB.PDBParser().get_structure("substructure", f"{substructure_data_dir}/substructure.pdb").get_atoms()):
                if substructure_atom.full_id[1:] == atom.full_id[1:]:
                    atom_substructure_index = si
                    break

            charges.append(float(xtb_output_file_lines[charge_headline_index + atom_substructure_index + 1].split()[3]))

            print(i, len(nearest_atoms), time() - t)
        print(charges)
        original = [-0.79196,
         0.11357,
         0.16227,
         -0.55082,
         0.20112,
         -0.45673,
         -0.19932,
         0.35776,
         -0.43203,
         0.7385,
         -0.67032,
         -0.63281,
         0.48483,
         -0.78635,
         -0.10857,
         0.17081,
         0.1205,
         0.1303,
         0.15295,
         0.17801,
         0.17086,
         0.16596,
         0.52675,
         0.48009,
         0.14171,
         0.16793,
         0.16154,
         0.56678,
         0.54527,
         -0.55008,
         -0.55026,
         -0.46457,
         0.11796,
         0.1516,
         -0.55454,
         0.19825,
         -0.46017,
         -0.19304,
         0.3536,
         -0.4308,
         0.73208,
         -0.66081,
         -0.64521,
         0.47787,
         -0.81451,
         -0.11817,
         0.17731,
         0.13056,
         0.12056,
         0.14428,
         0.18531,
         0.17259,
         0.15488,
         0.51826,
         0.46949,
         0.14029,
         0.13905,
         0.15779,
         0.54309,
         -0.55629,
         -0.56637,
         -0.45234,
         0.12383,
         0.1579,
         -0.54426,
         0.1938,
         -0.47833,
         -0.21903,
         0.3576,
         -0.45594,
         0.72902,
         -0.6113,
         -0.62328,
         0.57579,
         -0.63892,
         -0.04918,
         -0.22961,
         0.14311,
         0.13753,
         0.11995,
         0.14247,
         0.17526,
         0.17007,
         0.13201,
         0.13025,
         0.12678,
         0.44715,
         0.12683,
         0.14951,
         0.18254,
         0.54609,
         -0.56043,
         -0.58905,
         -0.45623,
         0.12032,
         0.158,
         -0.54198,
         0.20012,
         -0.46968,
         -0.20133,
         0.35409,
         -0.43772,
         0.72887,
         -0.63624,
         -0.65477,
         0.57144,
         -0.62349,
         -0.04934,
         -0.22735,
         0.14839,
         0.13049,
         0.12388,
         0.13637,
         0.1794,
         0.19031,
         0.13826,
         0.12691,
         0.12203,
         0.48383,
         0.14323,
         0.14499,
         0.16811,
         0.5449,
         -0.55975,
         -0.56887,
         -0.48059,
         0.11334,
         0.15525,
         -0.56427,
         0.19374,
         -0.47969,
         -0.21369,
         0.34965,
         -0.45243,
         0.72703,
         -0.63782,
         -0.66407,
         0.57437,
         -0.62792,
         -0.05518,
         -0.22664,
         0.14729,
         0.13533,
         0.12325,
         0.11963,
         0.16609,
         0.17329,
         0.13059,
         0.12684,
         0.12562,
         0.48837,
         0.14358,
         0.14671,
         0.15975,
         0.54014,
         -0.56903,
         -0.56815,
         -0.46769,
         0.12665,
         0.15555,
         -0.56271,
         0.1957,
         -0.46853,
         -0.21728,
         0.3485,
         -0.45018,
         0.72727,
         -0.62138,
         -0.63988,
         0.56971,
         -0.62528,
         -0.06734,
         -0.23287,
         0.14871,
         0.12056,
         0.12605,
         0.11226,
         0.16062,
         0.16845,
         0.13568,
         0.12861,
         0.11804,
         0.47259,
         0.13452,
         0.15225,
         0.16619,
         0.5464,
         -0.56953,
         -0.56939,
         -0.48496,
         0.12225,
         0.15927,
         -0.55506,
         0.18777,
         -0.79932,
         -0.19324,
         0.35026,
         -0.41702,
         0.73778,
         -0.66052,
         -0.63284,
         0.47699,
         -0.79085,
         -0.11366,
         0.16053,
         0.12126,
         0.12258,
         0.13738,
         0.19055,
         0.17215,
         0.16523,
         0.51188,
         0.47938,
         0.13125,
         0.16224,
         0.17231,
         0.57473,
         -0.82437,
         0.10797,
         0.15364,
         -0.5719,
         0.19236,
         -0.47123,
         -0.21751,
         0.36243,
         -0.4249,
         0.39053,
         -0.48152,
         0.16175,
         0.58564,
         -0.60483,
         -0.65857,
         0.65038,
         -0.77276,
         -0.52523,
         0.42558,
         0.09944,
         0.12062,
         0.13753,
         0.15886,
         0.1744,
         0.52193,
         0.52162,
         0.49096,
         0.13168,
         0.15399,
         0.16149,
         0.56836,
         0.54276,
         -0.57306,
         -0.57437,
         -0.47495,
         0.1153,
         0.16355,
         -0.55057,
         0.19577,
         -0.45839,
         -0.20781,
         0.35941,
         -0.43407,
         0.38435,
         -0.53394,
         0.15406,
         0.4541,
         -0.77029,
         -0.5921,
         0.40739,
         -0.51655,
         0.41025,
         0.12772,
         0.12693,
         0.14532,
         0.16623,
         0.17439,
         0.50977,
         0.48204,
         0.16358,
         0.13283,
         0.1558,
         0.16647,
         0.54261,
         -0.56527,
         -0.56778,
         -0.47792,
         0.12008,
         0.16222,
         -0.55288,
         0.19998,
         -0.46599,
         -0.20369,
         0.36338,
         -0.44456,
         0.36293,
         -0.52582,
         0.14589,
         0.466,
         -0.75943,
         -0.60758,
         0.40716,
         -0.55392,
         0.41587,
         0.13368,
         0.13211,
         0.15462,
         0.16195,
         0.16722,
         0.49596,
         0.47873,
         0.15813,
         0.13792,
         0.16503,
         0.16625,
         0.54421,
         -0.57536,
         -0.55262,
         -0.45469,
         0.11912,
         0.15688,
         -0.56704,
         0.19554,
         -0.47461,
         -0.20648,
         0.35035,
         -0.43183,
         0.38441,
         -0.53471,
         0.1522,
         0.46527,
         -0.76896,
         -0.60475,
         0.39081,
         -0.53619,
         0.40283,
         0.13425,
         0.12279,
         0.11095,
         0.14836,
         0.17708,
         0.49777,
         0.47774,
         0.17284,
         0.14568,
         0.14796,
         0.15431,
         0.54378,
         -0.56531,
         -0.56659,
         -0.47891,
         0.11553,
         0.15421,
         -0.57626,
         0.19707,
         -0.48239,
         -0.20579,
         0.35244,
         -0.44397,
         0.36966,
         -0.53212,
         0.16963,
         0.47579,
         -0.80709,
         -0.56727,
         0.37728,
         -0.43705,
         0.36276,
         0.11894,
         0.13337,
         0.14248,
         0.14635,
         0.16199,
         0.50511,
         0.47517,
         0.19139,
         0.14514,
         0.16063,
         0.16225,
         0.54187,
         -0.57113,
         -0.55851,
         -0.46192,
         0.11539,
         0.15095,
         -0.56534,
         0.18794,
         -0.46439,
         -0.21025,
         0.3597,
         -0.43329,
         0.37885,
         -0.50238,
         0.13138,
         0.57,
         -0.62489,
         -0.67751,
         0.64305,
         -0.76247,
         -0.53683,
         0.40452,
         0.12609,
         0.11883,
         0.15488,
         0.16184,
         0.18773,
         0.51495,
         0.52056,
         0.48017,
         0.13692,
         0.14653,
         0.15624,
         0.54078,
         -0.56893,
         -0.56214,
         -0.48735,
         0.11062,
         0.15105,
         -0.57022,
         0.16977,
         -0.78661,
         -0.23499,
         0.36466,
         -0.43554,
         0.37358,
         -0.48724,
         0.15266,
         0.58791,
         -0.61163,
         -0.65919,
         0.6515,
         -0.76521,
         -0.52632,
         0.41836,
         0.12491,
         0.11076,
         0.13956,
         0.16295,
         0.17042,
         0.52448,
         0.52139,
         0.49705,
         0.12186,
         0.16051,
         0.16486,
         0.5409,
         -0.52795,
         0.14606,
         0.43436,
         -0.04293,
         0.21071,
         0.32251,
         -0.07222,
         0.23152,
         0.04394,
         0.03043,
         -0.07819,
         0.01949,
         -0.4041,
         0.57207,
         0.17859,
         -0.0366,
         -0.01671,
         -0.05155,
         0.27754,
         0.25138,
         0.28578,
         0.16477,
         -0.54591,
         -0.47571,
         0.0309,
         -0.51609,
         0.03798,
         -0.50816,
         0.04062,
         -0.54505,
         -0.70041,
         -0.54602,
         -0.21794,
         0.65671,
         -0.54092,
         -0.51931,
         0.0417,
         0.51897,
         0.19319,
         0.16593,
         0.16616,
         0.19975,
         0.1468,
         0.178,
         0.15668,
         0.16688,
         0.5311,
         0.13704,
         0.13421,
         0.12275,
         0.10912,
         0.12814,
         0.10819,
         0.10691,
         0.13279,
         0.11848,
         0.60139,
         0.14117,
         0.13814,
         0.14439,
         0.1179,
         0.13603,
         0.11579]
        from bokeh.plotting import figure, show
        f = figure()
        f.circle(original, charges, radius=0.1)
        show(f)
        import numpy as np
        charges = np.array(charges)
        original = np.array(original)
        rmsd = np.sqrt((((charges - original)**2).sum()).mean())
        print(rmsd)
        r = np.corrcoef(charges, original)
        print(r)

        exit()


# $constrain
#    elements: 6,7,8,16
#    force constant=1.0
# $end




    def fix_structure(self):
        print("Fixing structure by PDBFixer... ", end="")
        fixer = PDBFixer(filename=self.mmCIF_file)
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/fixed.cif", 'w'),keepIds=True)
        print("ok")


    def protonate_ligands(self):
        """
        všechny vodíky jsou odstraněny


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
        # todo co tady s tím vlastně děláme
        print("Adding hydrogens to heteroresidues... ", end="")
        biotite_mmCIF_file = biotite_mmCIF.PDBxFile.read(self.mmCIF_file)
        biotite_protein = biotite_mmCIF.get_structure(biotite_mmCIF_file,
                                                      model=1,
                                                      extra_fields=["charge"],
                                                      include_bonds=True)
        biotite_protein.bonds = biotite_structure.connect_via_residue_names(biotite_protein)

        biotite_protein = biotite_protein[biotite_protein.element != "H"] # remove all hydrogens
        bond_array = biotite_protein.bonds.as_array()
        unknown_order_mask = bond_array[:, 2] == biotite_structure.BondType.ANY
        if unknown_order_mask.any():
            bond_array[unknown_order_mask, 2] = biotite_structure.BondType.SINGLE
            biotite_protein.bonds = biotite_protein.BondList(biotite_protein.array_length(), bond_array)


        biotite_protein_with_hydrogens, _ = hydride.add_hydrogen(biotite_protein, mask=biotite_protein.hetero)
        biotite_mmCIF_file_with_hydrogens = biotite_mmCIF.PDBxFile()
        biotite_mmCIF.set_structure(biotite_mmCIF_file_with_hydrogens, biotite_protein_with_hydrogens)
        biotite_mmCIF_file_with_hydrogens.write(f"{self.data_dir}/protonated_ligands.cif")
        # todo vázané ligandy jsou protonované špatně kvůli nenačtení vazby - github issue
        # todo jména vodíků rezidua jsou duplicitní a proto jsou pdb fixerem v další iteraci smazány smazány - github issue

    def protonate_protein(self):

        print("Fixing structure by PDBFixer... ", end="")
        fixer = PDBFixer(filename=f"{self.data_dir}/protonated_ligands.cif")
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(f"{self.data_dir}/protonated_ligands_fixed.cif", 'w'))
        print("ok")

        print("Adding hydrogens to protein... ", end="")
        molecule = molit.Molecule(f"{self.data_dir}/protonated_ligands.cif")
        prepared_molecule = systemPrepare(molecule, pH=self.pH)
        prepared_molecule.write(f"{self.data_dir}/protonated_protein.cif")
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
    ChargesCalculator(args.mmCIF_file, args.data_dir, args.pH).calculate_charges()


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



# belbínův test osobnosti



# možná nepůjde stáhnout všechno
# určitě log pro všechny rezidua
# nechat si vždycky verzi, se kteoru se pracovalo

#todo - volba smazání tmp souborů.







