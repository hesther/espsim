error=False
try:
    import numpy
except:
    print("ERROR: Numpy could not be imported.")
    error=True
try:
    import scipy.spatial
except:
    print("ERROR: Scipy.spatial could not be imported.")
    error=True
try:
    import rdkit
except:
    print("ERROR: RDKit could not be imported.")
    error=True
try:
    import espsim
except:
    print("ERROR: epsim could not be imported.")
    error=True

if not error:
    print("Test passed, imports work fine.")
