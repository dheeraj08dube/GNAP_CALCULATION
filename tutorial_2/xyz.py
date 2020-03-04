from Bio import PDB

parser = PDB.PDBParser()
io = PDB.PDBIO()
struct = parser.get_structure('ref','ref.pdb')

for model in struct:
    for chain in model:
        for residue in chain:
            for atom in residue:
                x,y,z = atom.get_coord()
                print(x,y,z)
                
                
                
import parmed
pdb = parmed.load_file('ref.pdb')
A = pdb.get_box()
print(A[0][0],A[0][1],A[0][2])                
                
