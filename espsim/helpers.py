import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from argparse import Namespace
import os
import contextlib
import pkg_resources
import warnings
import logging

logging_level = logging.root.level

def check_hs(mol):
    """Checks whether a molecule misses explicit hydrogens and throws a warning.
    :param mol: RDKit molecule
    """
    if not mol.GetNumAtoms() == Chem.AddHs(mol).GetNumAtoms():
        warnings.warn('Input molecule missing hydrogens. Use Chem.AddHs() to make hydrogens explicit.')
        
def Renormalize(similarity,
                metric = "carbo",
                customrange = None,
):
    """
    Renormalizes a similarity metric to the range [0:1]
    :param similarity: Similarity score.
    :param mode: (optional) Mode of similarity score
    :param customrange: (optional) Custom range of similarity score, overrides mode parameter. Tuple or list of two values.
    :return: Renormalized similarity score
    """
    if customrange != None:
        similarity=(similarity-customrange[0])/(customrange[1]-customrange[0])
    elif metric=="carbo":
        similarity=(similarity+1)/2
    elif metric=="tanimoto":
        similarity=(similarity+1/3)/(4/3)
    else:
        raise ValueError("Unknown metric.")
    return similarity

def SimilarityMetric(intPrbPrb,
                     intRefRef,
                     intPrbRef,
                     metric = "carbo",
):
    """
    Calculates a similarity metrics from the overlap integrals of the electrostatic potentials
    a probe and reference molecule.
    :param intPrbPrb: Value of self-overlap integral of probe molecule.
    :param intRefRef: Value of self-overlap integral of reference molecule.
    :param intPrbRef: Value of overlap integral between probe and reference molecule.
    :param mode: (optional) Similarity metric.
    :return: Similarity score
    """
    if metric=='carbo':
        numerator=intPrbRef
        denominator=np.sqrt(intPrbPrb*intRefRef)
    elif metric == 'tanimoto':
        numerator=intPrbRef
        denominator=intPrbPrb+intRefRef-intPrbRef
    else:
        raise ValueError("Unknown metric.")

    if denominator!= 0:
        similarity=numerator/denominator
    else:
        raise ValueError("Denominator in similarity calculation equals zero.")
    return similarity

try:
    import psi4
    import resp

    def psi4Charges(xyz,
                    basisPsi4 = '3-21G',
                    methodPsi4 = 'scf',
                    gridPsi4 = 1,
    ):
        """
        Calculates RESP charges via Psi4.
        :param xyz: String of xyz file of an embedded molecule.
        :param basisPsi4: (optional) Basis set.
        :param methodPsi4: (optional) Method.
        :param gridPsi4: (optional) Integer grid point density for ESP evaluation.
        :return: Array of RESP partial charges.
        """
        mol = psi4.core.Molecule.from_string(xyz, dtype='xyz')
        mol.update_geometry()

        options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
                   'VDW_POINT_DENSITY'  : int(gridPsi4),
                   'RESP_A'             : 0.0005,
                   'RESP_B'             : 0.1,
                   'BASIS_ESP' : basisPsi4,
                   'METHOD_ESP' : methodPsi4,
        }

        charge = resp.resp([mol], options)[1]
        return charge
except ImportError:
    def psi4Charges(xyz,
                    basisPsi4 = '3-21G',
                    methodPsi4 = 'scf',
                    gridPsi4 = 1,
    ):
        """
        Mock implementation raising an ImportError if psi4 and resp cannot be imported.
        """
        raise ImportError("Failed to import Psi4 and RESP. Please install via 'conda install -c psi4 psi4 resp'")
finally:
    logging.getLogger().setLevel(logging_level)


def num_atoms_bonds(smiles):
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    return len(m.GetAtoms()), len(m.GetBonds())

def get_reorder_list(mol):
    """ Returns the order of atoms or bonds of the current molecule compared to a molecule in Chemprop
    :param mol: RDKit molecule object
    :return: List of indices, to be used to reorder molecules or atom/bond feature lists.
    """
    chemprop_mol=Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))

    if chemprop_mol.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError("Could not set up canonical molecule for reordering.")

    atom_order=list(chemprop_mol.GetSubstructMatch(mol))    

    if chemprop_mol.GetNumAtoms() != len(atom_order):
        raise ValueError("Could not reorder partial charges for molecule of length",chemprop_mol.GetNumAtoms(),"with subgraph match of only length",len(atom_order))
    
    return atom_order

try:
    from chemprop.train import make_predictions

    def mlCharges(mols):
        """
        Calc ML charges
        """

        #MUST use mols with hydrogens!
        smiles = [Chem.MolToSmiles(mol) for mol in mols]

        path = pkg_resources.resource_filename(__name__, 'QM_137k.pt')
        args = Namespace(batch_size=50, checkpoint_dir=None, checkpoint_path=path, checkpoint_paths=[path], cuda=False, features_generator=None, features_path=None, gpu=None, max_data_size=None, no_features_scaling=False, preds_path=None, test_path=None, use_compound_names=False)

        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                test_preds, test_smiles = make_predictions(args, smiles=smiles)
        n_atoms, n_bonds = zip(*[num_atoms_bonds(x) for x in smiles])
        partial_charge = test_preds[0]
        partial_charge = np.split(partial_charge.flatten(), np.cumsum(np.array(n_atoms)))[:-1]

        charges=[]
        for i,mol in enumerate(mols):
            try:
                reorder_list = get_reorder_list(mol)
                charges.append([partial_charge[i][reorder_list[x]] for x in range(mol.GetNumAtoms())])
            except ValueError:
                #Could not get prediction, default to Gasteiger
                print("Warning: could not obtain prediction, defaulting to Gasteiger charges for one molecule")
                AllChem.ComputeGasteigerCharges(mol)
                charges.append([a.GetDoubleProp('_GasteigerCharge') for a in mol.GetAtoms()])
            

        return charges

except ImportError:
    def mlCharges(mols):
        """
        Mock implementation raising an ImportError if the ML model cannot be imported.
        """
        raise ImportError("Failed to import the atomic Chemprop model. Please install via 'pip install git+https://github.com/hesther/chemprop-atom-bond.git'")


def readMolFile(f):
    """
    Reads a molecule and its coordinates from a mol file.
    :param f: Path to file.
    :return: RDKit molecule object.
    """
    try:
        mol=Chem.MolFromMolFile(f,removeHs=False)
    except:
        raise ValueError("File could not be read.")
    return mol

def readMol2File(f):
    """
    Reads a molecule and its coordinates and charges from a mol2 file.
    :param f: Path to file.
    :return: RDKit molecule object, list of partial charges
    """
    try:
        mol=Chem.MolFromMol2File(f,removeHs=False)
    except:
        raise ValueError("File could not be read. Supported mol2 format: Corina")
    charge=[atom.GetDoubleProp("_TriposPartialCharge") for atom in mol.GetAtoms()]
    return mol, charge

def readSdfFile(f,key='CHARGES'):
    """
    Read one molecule from an SDF file, and atomic features from the property block named according to the variable key.
    :param f: String of SDF file location.
    :param key: Name of the property block.
    :return: RDKit molecule, list of features (floats).
    """
    mol=Chem.SDMolSupplier(f,removeHs=False)[0]
    if key not in mol.GetPropsAsDict().keys(): 
        raise ValueError("Unknown property key supplied. Check choice of key and/or the supplied sdf file") 
    charge=list(mol.GetProp(key).split(","))
    if len(charge) != mol.GetNumAtoms(): 
        raise ValueError("List of partial charges must contain exactly N_atoms entries.")
    return mol,charge
