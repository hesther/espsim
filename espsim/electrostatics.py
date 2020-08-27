from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AlignMol, EmbedMolecule, EmbedMultipleConfs
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
import numpy as np
import scipy.spatial


def GetMolProps(mol,cid,charge):
    """
    Extracts the coordinates, van der Waals radii and charges from a given conformer cid of a molecule mol.
    :param mol: RDKit mol object.
    :param cid: Index of the conformer for 3D coordinates.
    :param charge: List or array of charge. If None, RDKit Gasteiger Charges are read, or if not available yet, calculated.
    :return: 2D array of coordinates, 1D array of charges.
    """
    
    coor=mol.GetConformer(cid).GetPositions()
    if charge == None:
        try:
            charge=np.array([np.float(a.GetProp('_GasteigerCharge')) for a in mol.GetAtoms()])
        except KeyError:
            AllChem.ComputeGasteigerCharges(mol)
            charge=np.array([np.float(a.GetProp('_GasteigerCharge')) for a in mol.GetAtoms()])
    else:
        charge=np.array(charge,dtype=np.float).flatten()

    return coor,charge
    
def GetEspSim(prbMol,refMol,prbCid=-1,refCid=-1,prbCharge=None,refCharge=None):
    """
    Calculates the similarity of the electrostatic potential around two previously aligned molecules.
    :param prbMol: RDKit mol object of the probe molecule.
    :param refMol: RDKit mol object of the reference molecule.
    :param prbCid: Index of the conformer of the probe molecule to be used for 3D coordinates.
    :param refCid: Index of the conformer of the reference molecule to be used for 3D coordinates.
    :param prbCharge: (optional) List or array of partial charges of the probe molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param refCharge: (optional) List or array of partial charges of the reference molecule. If not given, RDKit Gasteiger Charges are used as default.
    :return: Similarity score.
    """

    #Set up probe molecule properties:
    prbCoor,prbCharge=GetMolProps(prbMol,prbCid,prbCharge)
    refCoor,refCharge=GetMolProps(refMol,refCid,refCharge)

    similarity=GetIntegralsViaGaussians(prbCoor,refCoor,prbCharge,refCharge)

    return similarity

def GetIntegralsViaGaussians(prbCoor,refCoor,prbCharge,refCharge):
    """
    Calculates the integral of the overlap between the point charges prbCharge and refCharge at coordinates prbCoor and refCoor via fitting to Gaussian functions and analytic integration.
    :param prbCoor: 2D array of coordinates of the probe molecule. 
    :param refCoor: 2D array of coordinates of the reference molecule. 
    :param prbCharge: 1D array of partial charges of the probe molecule. 
    :param refCharge: 1D array of partial charges of the reference molecule.
    :return: Similarity of the overlap integrals.
    """

    distPrbPrb = scipy.spatial.distance.cdist(prbCoor,prbCoor)
    distPrbRef = scipy.spatial.distance.cdist(prbCoor,refCoor)
    distRefRef = scipy.spatial.distance.cdist(refCoor,refCoor)

    intPrbPrb=GaussInt(distPrbPrb,prbCharge,prbCharge)
    intPrbRef=GaussInt(distPrbRef,prbCharge,refCharge)
    intRefRef=GaussInt(distRefRef,refCharge,refCharge)

    similarity=intPrbRef/np.sqrt(intPrbPrb*intRefRef)
    return similarity
  
def GaussInt(dist,charge1,charge2):
    """Calculates the analytic Gaussian integrals.
    :param dist: Distance matrix.
    :param charge1: 1D array of partial charges of first molecule.
    :param charge2: 1D array of partial charges of second molecule.
    :return: Analytic overlap integral.
    """

    #These are precomputed coefficients:
    a=np.array([[ 15.90600036,   3.9534831 ,  17.61453176],[  3.9534831 ,   5.21580206,   1.91045387],[ 17.61453176,   1.91045387, 238.75820253]])
    b=np.array([[-0.02495   , -0.04539319, -0.00247124],[-0.04539319, -0.2513    , -0.00258662],[-0.00247124, -0.00258662, -0.0013    ]])
    
    intOverall=0
    for idx1 in range(charge1.shape[0]):
        for idx2 in range(charge2.shape[0]):
            intGauss=0
            for i in range(3):
                intGauss+=a[i,i]*np.exp(b[i,i]*dist[idx1,idx2]**2) #Diagonal (11, 22, 33)
                for j in range(i+1,3):
                    intGauss+=2*a[i,j]*np.exp(b[i,j]*dist[idx1,idx2]**2) #Off diagonal, 12, 13, 21, 23, 31, 32 = 2*(12,13,23)
            intOverall+=charge1[idx1]*charge2[idx2]*intGauss
    return intOverall


def ConstrainedEmbedMultipleConfs(mol, core, numConfs=10, useTethers=True, coreConfId=-1, randomseed=2342,
                     getForceField=UFFGetMoleculeForceField, **kwargs):
    """
    Function to obtain multiple constrained embeddings per molecule. This was taken as is from:
    from https://github.com/rdkit/rdkit/issues/3266
    :param mol: RDKit molecule object to be embedded.
    :param core: RDKit molecule object of the core used as constrained. Needs to hold at least one conformer coordinates.
    :param numCons: Number of conformations to create
    :param useTethers: (optional) boolean whether to pull embedded atoms to core coordinates, see rdkit.Chem.AllChem.ConstrainedEmbed
    :param coreConfId: (optional) id of the core conformation to use
    :param randomseed: (optional) seef for the random number generator
    :param getForceField: (optional) force field to use for the optimization of molecules
    :return: RDKit molecule object containing the embedded conformations.
    """

    match = mol.GetSubstructMatch(core)
    if not match:
        raise ValueError("molecule doesn't match the core")
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI

    cids = EmbedMultipleConfs(mol, numConfs=numConfs, coordMap=coordMap, randomSeed=randomseed, **kwargs)
    cids = list(cids)
    if len(cids) == 0:
        raise ValueError('Could not embed molecule.')

    algMap = [(j, i) for i, j in enumerate(match)]
    
    if not useTethers:
        # clean up the conformation
        for cid in cids:
            ff = getForceField(mol, confId=cid)
            for i, idxI in enumerate(match):
                for j in range(i + 1, len(match)):
                    idxJ = match[j]
                    d = coordMap[idxI].Distance(coordMap[idxJ])
                    ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.)
            ff.Initialize()
            n = 4
            more = ff.Minimize()
            while more and n:
                more = ff.Minimize()
                n -= 1
            # rotate the embedded conformation onto the core:
            rms = AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        for cid in cids:
            rms = AlignMol(mol, core, prbCid=cid, atomMap=algMap)
            ff = getForceField(mol, confId=cid)
            conf = core.GetConformer()
            for i in range(core.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
                ff.AddDistanceConstraint(pIdx, match[i], 0, 0, 100.)
            ff.Initialize()
            n = 4
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            while more and n:
                more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
                n -= 1
            # realign
            rms = AlignMol(mol, core, prbCid=cid, atomMap=algMap)
    return mol


def EmbedAlignConstrainedScore(prbMol,refMols,core,prbNumConfs=10,refNumConfs=10,prbCharge=None,refCharges=None):
    """Calculates a constrained alignment based on a common pattern in the input molecules. Caution: Will fail if the pattern does not match. 
    Calculates a shape and electrostatic potential similarity of the best alignment.

    :param prbMol: RDKit molecule for which shape and electrostatic similarities are calculated.
    :param refMol: RDKit molecule or list of RDKit molecules serving as references.
    :param core: Common pattern for the constrained embedding as embedded RDKit molecule
    :param prbNumConfs: Number of conformers to create for the probe molecule. A higher number creates better alignments but slows down the algorithm.
    :param refNumConfs: Number of conformers to create for each reference molecule. A higher number creates better alignments but slows down the algorithm.
    :param prbCharge: (optional) List or array of partial charges of the probe molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param refCharge: (optional) List of list or 2D array of partial charges of the reference molecules. If not given, RDKit Gasteiger Charges are used as default.
    :return: shape similarity and ESP similarity.
    """
    
    if type(refMols) != list:
        refMols=[refMols]

    if refCharges == None:
        refCharges=[None]*len(refMols)
        
    prbMol=ConstrainedEmbedMultipleConfs(prbMol, core, numConfs=prbNumConfs)
    for refMol in refMols:
        refMol=ConstrainedEmbedMultipleConfs(refMol, core, numConfs=refNumConfs)
        
    prbMatch = prbMol.GetSubstructMatch(core)
    allShapeDist = []
    allEspSim = []
    
    for idx,refMol in enumerate(refMols):
        shapeDist=1
        prbBestConf=0
        refBestConf=0
        refMatch = refMol.GetSubstructMatch(core)
        for i in range(refNumConfs):
            for j in range(prbNumConfs):
                AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=j,refCid=i)
                shape = AllChem.ShapeTanimotoDist(prbMol,refMol,confId1=j,confId2=i)
                if shape<shapeDist:
                    shapeDist=shape
                    prbBestConf=j
                    refBestConf=i
        espSim=GetEspSim(prbMol,refMol,prbBestConf,refBestConf,prbCharge,refCharges[idx])
        allShapeDist.append(1-shapeDist)
        allEspSim.append(espSim)

    return allShapeDist,allEspSim

