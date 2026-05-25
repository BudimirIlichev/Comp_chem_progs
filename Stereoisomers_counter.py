from rdkit import Chem
mol=Chem.MolFromSmiles("enter your molecule's smiles code")
if mol is None:
    print("Invalid SMILES string")
else: 
 s=0
for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'C':
           neighbors = atom.GetNeighbors()
    subs=[n for n in neighbors if n.GetSymbol() != 'H']
    if len(subs)==4:
            print("Chiral molecule")
            s+=1
print(f"Number of chiral centers: {s}")
Chir_Isoms=2**s
print(f"Number of stereoisomers: {Chir_Isoms}")
