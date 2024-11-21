Requirements:

Python 3.11.6

openbabel 3.1.1       # conda install conda-forge::openbabel=3.1.1
xtb 6.6.1             # conda install conda-forge::xtb=6.6.1

pip install dimorphite_dl==1.3.2 hydride==1.2.3 biopython==1.84 numpy==1.26.4 biotite=1.0.1 gemmi=0.6.6 rdkit==2023.09.6 moleculekit==1.9.15 pdb2pqr==3.6.1 openmm==8.2.0 pdbfixer==1.10




je potřeba upravit moleculekit
v souboru miniconda3/lib/python3.11/site-packages/moleculekit/tools/preparation.py 
ve funkci _biomolecule_to_molecule do seznamu propmap dopsat ("ffcharge", "charge")
v definici funkce _pdb2pqr změnit defaultní nastavení parametru opt=False, 


Pořešit:

jak cesty k pdb2pqr, hydride, xtb?
výsledky potřebujeme někam uložit aby si to pak webovka mohla tahat
možná nepůjde stáhnout všechno
určitě log pro všechny rezidua
nechat si vždycky verzi, se kteoru se pracovalo
