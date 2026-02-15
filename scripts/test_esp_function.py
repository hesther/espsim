from espsim import GetEspSim, readMolFile, readMol2File, readSdfFile
from importlib.resources import files, as_file

# Get paths to test files from package data
testfiles = files('espsim.testfiles')

#The following block of code reads in prealigned molecules in mol format and calculates their ESP similarity:
with as_file(testfiles.joinpath('prbmol1.mol')) as prbmol_path, \
     as_file(testfiles.joinpath('refmol1.mol')) as refmol_path:
    mol1=readMolFile(str(prbmol_path))
    mol2=readMolFile(str(refmol_path))
    sim_esp=GetEspSim(mol1,mol2)
    print("%15s %5.2f" % ("ESP similarity (mols read from mol file):",sim_esp))


#The following block of code reads in prealigned molecules in mol2 format with custom charges and calculates their ESP similarity:
with as_file(testfiles.joinpath('prbmol1.mol2')) as prbmol_path, \
     as_file(testfiles.joinpath('refmol1.mol2')) as refmol_path:
    mol1,charge1=readMol2File(str(prbmol_path))
    mol2,charge2=readMol2File(str(refmol_path))
    sim_esp=GetEspSim(mol1,mol2,prbCharge=charge1,refCharge=charge2)
    print("%15s %5.2f" % ("ESP similarity (mols read from mol2 file):",sim_esp))


#The following block of code reads in prealigned molecules in sdf format with custom charges as a comma-separated list
#in the SDF file and calculates their ESP similarity:
with as_file(testfiles.joinpath('prbmol1.sdf')) as prbmol_path, \
     as_file(testfiles.joinpath('refmol1.sdf')) as refmol_path:
    mol1,charge1=readSdfFile(str(prbmol_path))
    mol2,charge2=readSdfFile(str(refmol_path))
    sim_esp=GetEspSim(mol1,mol2,prbCharge=charge1,refCharge=charge2)
    print("%15s %5.2f" % ("ESP similarity (mols read from sdf file):",sim_esp))
