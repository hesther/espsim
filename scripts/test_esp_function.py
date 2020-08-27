from rdkit import Chem
from espsim import GetEspSim

#The following block of code reads in prealigned molecules:
with open('scripts/prbmol1.mol') as f:
    read_data = f.read()
print("#"*30)
print("Loading molecule A:")
print(read_data)
mol11=Chem.MolFromMolBlock(read_data,sanitize=False)
Chem.SanitizeMol(mol1,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
with open('scripts/refmol1.mol') as f:
    read_data = f.read()
print("#"*30)
print("Loading molecule B:")
print(read_data)
mol2=Chem.MolFromMolBlock(read_data,sanitize=False)
Chem.SanitizeMol(mol2,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
print("#"*30)


#Now here is the actual calculation of similarity, it is a single function call:
sim_esp=GetEspSim(mol1,mol2)
print("%15s %5.2f" % ("ESP similarity:",sim_esp))
