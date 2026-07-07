import numpy as np
import math
from rdkit import Chem
Rad_cutoff=2
RADIUS_of_ATOMS={
    'H':0.31,
    'C':0.76,
   'N':0.71,
   'O':0.66
}
file_xyz=input("filename:")
atoms=[]
with open(file_xyz,"r",encoding="utf-8") as file:
    lines=file.readlines()
    num_atoms=int(lines[0].strip())
    for line in lines[2:]:
        parts=line.split()
        if len(parts)==4:
            element=parts[0]
            x,y,z=float(parts[1]),float(parts[2]),float(parts[3])
            atoms.append({"element":element,"coords":(x,y,z)})
corr_matrix=np.zeros((len(atoms),len(atoms)))

for i in range(len(atoms)):
    for j in range (i+1,len(atoms)):
        x1,y1,z1=atoms[i]["coords"]
        x2,y2,z2=atoms[j]["coords"]
        Rij=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

        r_i=RADIUS_of_ATOMS.get(atoms[i]["element"],0.7)
        r_j=RADIUS_of_ATOMS.get(atoms[j]["element"],0.7)
        if 0.9<=Rij<=1.6:
            corr_matrix[i][j]=1
            corr_matrix[j][i]=1

for row in corr_matrix:
    print(row)
s=0
for idx,row in enumerate(corr_matrix):
    Chir=int(row.sum())
    if Chir==4:
        s=s+1
        chiral_isomers=2**s

print("Number of chiral isomers:", chiral_isomers)



