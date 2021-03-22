from espsim import GetEspSim
from espsim.helpers import readMolFile, readMol2File, readSdfFile

#The following block of code reads in prealigned molecules in mol format and calculates their ESP similarity:
mol1=readMolFile("scripts/prbmol1.mol")
mol2=readMolFile("scripts/refmol1.mol")
sim_esp=GetEspSim(mol1,mol2)
print("%15s %5.2f" % ("ESP similarity (mols read from mol file):",sim_esp))


#The following block of code reads in prealigned molecules in mol2 format with custom charges and calculates their ESP similarity:
mol1,charge1=readMol2File("scripts/prbmol1.mol2")
mol2,charge2=readMol2File("scripts/refmol1.mol2")
sim_esp=GetEspSim(mol1,mol2,prbCharge=charge1,refCharge=charge2)
print("%15s %5.2f" % ("ESP similarity (mols read from mol2 file):",sim_esp))


#The following block of code reads in prealigned molecules in sdf format with custom charges as a comma-separated list
#in the SDF file and calculates their ESP similarity:
mol1,charge1=readSdfFile("scripts/prbmol1.sdf")
mol2,charge2=readSdfFile("scripts/refmol1.sdf")
sim_esp=GetEspSim(mol1,mol2,prbCharge=charge1,refCharge=charge2)
print("%15s %5.2f" % ("ESP similarity (mols read from sdf file):",sim_esp))
