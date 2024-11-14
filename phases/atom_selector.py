from Bio import PDB

class AtomSelector(PDB.Select):
    """
    This class is used by the biopython library to create substructures and store them in a file.
    More at https://biopython.org/docs/1.76/api/Bio.PDB.PDBIO.html?highlight=accept_atom#Bio.PDB.PDBIO.Select.accept_atom
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)