Requirements:
Python 3.11.6
biotite 1.0.1         # conda install -c conda-forge biotite=1.0.1
gemmi 0.6.6           # conda install -c conda-forge gemmi=0.6.6
rdkit 2023.09.6       # conda install -c conda-forge rdkit=2023.09.6
openmm 8.1.2          # conda install -c conda-forge openmm=8.1.2
moleculekit 1.9.15    # conda install -c acellera moleculekit=1.9.15
openbabel 3.1.1       # conda install conda-forge::openbabel=3.1.1
pdb2pqr 3.6.1         # conda install conda-forge::pdb2pqr=3.6.1
numpy 1.26.4          # conda install numpy=1.26.4
xtb 6.6.1             # conda install conda-forge::xtb=6.6.1

dimorphite_dl 1.3.2   # pip install dimorphite_dl==1.3.2
hydride 1.2.2         # pip install hydride==1.2.2
biopython 1.84        # pip install biopython==1.84


knihovnu pdbfixer je potřeba nainstalovat přímo z gitu. 
git clone https://github.com/openmm/pdbfixer
cd pdbfixer
python setup.py install

je potřeba upravit moleculekit
v souboru miniconda3/lib/python3.11/site-packages/moleculekit/tools/preparation.py 
ve funkci _biomolecule_to_molecule do seznamu propmap dopsat ("ffcharge", "charge")


Pořešit:

jak cesty k pdb2pqr, hydride, xtb?
výsledky potřebujeme někam uložit aby si to pak webovka mohla tahat
možná nepůjde stáhnout všechno
určitě log pro všechny rezidua
nechat si vždycky verzi, se kteoru se pracovalo
