from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AlignMol, EmbedMolecule, EmbedMultipleConfs
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
import numpy as np
import scipy.spatial
from .helpers import Renormalize, SimilarityMetric, psi4Charges, mlCharges, check_hs

def GetMolProps(mol,
                cid,
                charge,
                partialCharges = "gasteiger",
                basisPsi4 = '3-21G',
                methodPsi4 = 'scf',
                gridPsi4 = 1,
):
    """
    Extracts the coordinates, van der Waals radii and charges from a given conformer cid of a molecule mol.
    :param mol: RDKit mol object.
    :param cid: Index of the conformer for 3D coordinates.
    :param charge: List or array of charge. If empty list, charges are calculated based on the parameter partialCharges.
    :param partialCharges: (optional) Partial charge distribution.
    :param basisPsi4: (optional) Basis set for Psi4 calculation.
    :param methodPsi4: (optional) Method for Psi4 calculation.
    :param gridPsi4: (optional) Integer grid point density for ESP evaluation for Psi4 calculation.
    :return: 2D array of coordinates, 1D array of charges.
    """
    
    coor=mol.GetConformer(cid).GetPositions()
    if len(charge) == 0:
        if partialCharges == "gasteiger":
            try:
                charge=np.array([a.GetDoubleProp('_GasteigerCharge') for a in mol.GetAtoms()])
            except KeyError:
                AllChem.ComputeGasteigerCharges(mol)
                charge=np.array([a.GetDoubleProp('_GasteigerCharge') for a in mol.GetAtoms()])
        elif partialCharges == "mmff":
            mp = AllChem.MMFFGetMoleculeProperties(mol)
            if mp:
                charge=np.array([mp.GetMMFFPartialCharge(i) for i in range(mol.GetNumAtoms())])
            else:
                print("MMFF charges not available for the input molecule, defaulting to Gasteiger charges.")
                AllChem.ComputeGasteigerCharges(mol)
                charge=np.array([a.GetDoubleProp('_GasteigerCharge') for a in mol.GetAtoms()])
        elif partialCharges == 'ml':
            charge=np.array(mlCharges([mol])[0])

        elif partialCharges == "resp":
            xyz=Chem.rdmolfiles.MolToXYZBlock(mol,confId=cid)
            charge=psi4Charges(xyz,basisPsi4,methodPsi4,gridPsi4)
        else:
            raise ValueError("Unknown partial charge distribution.")
        if charge.shape[0] != coor.shape[0]:
            raise ValueError("Error in partial charge calculation.")
    else:
        charge=np.array(charge,dtype=float).flatten()
        if charge.shape[0] != coor.shape[0]:
            raise ValueError("Dimensions of the supplied charges does not match dimensions of coordinates of molecule")

    return coor,charge
    
def GetShapeSim(prbMol,
              refMol,
              prbCid = -1,
              refCid = -1):
    """
    Calculates the similarity of the shape between two previously aligned molecules.
    :param prbMol: RDKit mol object of the probe molecule.
    :param refMol: RDKit mol object of the reference molecule.
    :param prbCid: Index of the conformer of the probe molecule to be used for 3D coordinates.
    :param refCid: Index of the conformer of the reference molecule to be used for 3D coordinates.
    :return: Shape score
    """

    return 1 - AllChem.ShapeTanimotoDist(prbMol,refMol,confId1=prbCid,confId2=refCid)

def GetEspSim(prbMol,
              refMol,
              prbCid = -1,
              refCid = -1,
              prbCharge = [],
              refCharge = [],
              metric = "carbo",
              integrate = "gauss",
              partialCharges = "gasteiger",
              renormalize = False,
              customrange = None,
              marginMC = 10,
              nMC = 1,
              basisPsi4 = '3-21G',
              methodPsi4 = 'scf',
              gridPsi4 = 1,
              nocheck=False,
              randomseed = 2342,
):
    """
    Calculates the similarity of the electrostatic potential around two previously aligned molecules.
    :param prbMol: RDKit mol object of the probe molecule.
    :param refMol: RDKit mol object of the reference molecule.
    :param prbCid: Index of the conformer of the probe molecule to be used for 3D coordinates.
    :param refCid: Index of the conformer of the reference molecule to be used for 3D coordinates.
    :param prbCharge: (optional) List or array of partial charges of the probe molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param refCharge: (optional) List or array of partial charges of the reference molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param metric:  (optional) Similarity metric.
    :param integrate: (optional) Integration method.
    :param partialCharges: (optional) Partial charge distribution.
    :param renormalize: (optional) Boolean whether to renormalize the similarity score to [0:1].
    :param customrange: (optional) Custom range to renormalize to, supply as tuple or list of two values (lower bound, upper bound).
    :param marginMC: (optional) Margin up to which to integrate (added to coordinates plus/minus their vdW radii) if MC integration is utilized.
    :param nMC: (optional) Number of grid points per 1 Angstrom**3 volume of integration vox if MC integration is utilized.
    :param basisPsi4: (optional) Basis set for Psi4 calculation.
    :param methodPsi4: (optional) Method for Psi4 calculation.
    :param gridPsi4: (optional) Integer grid point density for ESP evaluation for Psi4 calculation.
    :param nocheck: (optional) whether no checks on explicit hydrogens should be run. Speeds up the function, but use wisely.
    :param randomseed: (optional) seed for the random number generator. Only used with the `mc` integration method.
    :return: Similarity score.
    """

    #Check hydrogens
    if not nocheck:
        check_hs(prbMol)
        check_hs(refMol)
    
    #Set up probe molecule properties:
    prbCoor,prbCharge=GetMolProps(prbMol,prbCid,prbCharge,partialCharges,basisPsi4,methodPsi4,gridPsi4)
    refCoor,refCharge=GetMolProps(refMol,refCid,refCharge,partialCharges,basisPsi4,methodPsi4,gridPsi4)

    if integrate=='gauss':
        similarity=GetIntegralsViaGaussians(prbCoor,refCoor,prbCharge,refCharge,metric)
    elif integrate=='mc':
        prbVdw = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in prbMol.GetAtoms()]).reshape(-1,1)
        refVdw = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in refMol.GetAtoms()]).reshape(-1,1)
        similarity=GetIntegralsViaMC(prbCoor,refCoor,prbCharge,refCharge,prbVdw,refVdw,metric,marginMC,nMC, randomseed=randomseed)

    if renormalize:
        similarity=Renormalize(similarity,metric,customrange)   

    return similarity

def GetIntegralsViaGaussians(prbCoor,
                             refCoor,
                             prbCharge,
                             refCharge,
                             metric,
):
    """
    Calculates the integral of the overlap between the point charges prbCharge and refCharge at coordinates prbCoor and refCoor via fitting to Gaussian functions and analytic integration.
    :param prbCoor: 2D array of coordinates of the probe molecule. 
    :param refCoor: 2D array of coordinates of the reference molecule. 
    :param prbCharge: 1D array of partial charges of the probe molecule. 
    :param refCharge: 1D array of partial charges of the reference molecule.
    :param metric: Metric of similarity score.
    :return: Similarity of the overlap integrals.
    """

    distPrbPrb = scipy.spatial.distance.cdist(prbCoor,prbCoor)
    distPrbRef = scipy.spatial.distance.cdist(prbCoor,refCoor)
    distRefRef = scipy.spatial.distance.cdist(refCoor,refCoor)

    intPrbPrb=GaussInt(distPrbPrb,prbCharge,prbCharge)
    intPrbRef=GaussInt(distPrbRef,prbCharge,refCharge)
    intRefRef=GaussInt(distRefRef,refCharge,refCharge)

    similarity=SimilarityMetric(intPrbPrb,intRefRef,intPrbRef,metric)
    return similarity
  
def GaussInt(dist,
             charge1,
             charge2,
):
    """Calculates the analytic Gaussian integrals.
    :param dist: Distance matrix.
    :param charge1: 1D array of partial charges of first molecule.
    :param charge2: 1D array of partial charges of second molecule.
    :return: Analytic overlap integral.
    """

    #These are precomputed coefficients:
    a=np.array([[ 15.90600036,   3.9534831 ,  17.61453176],
                [  3.9534831 ,   5.21580206,   1.91045387],
                [ 17.61453176,   1.91045387, 238.75820253]])
    b=np.array([[-0.02495   , -0.04539319, -0.00247124],
                [-0.04539319, -0.2513    , -0.00258662],
                [-0.00247124, -0.00258662, -0.0013    ]])

    a_flat = a.flatten() 
    b_flat = b.flatten()
    dist = (dist**2).flatten()
    charges = (charge1[:,None]*charge2).flatten() #pairwise products of atomic charges, flattened 
    return ((a_flat[:,None] * np.exp(dist * b_flat[:,None])).sum(0) * charges).sum()

def GetIntegralsViaMC(prbCoor,
                      refCoor,
                      prbCharge,
                      refCharge,
                      prbVdw,
                      refVdw,
                      metric,
                      marginMC = 10,
                      nMC = 1,
                      randomseed = 2342,
):
    """
    Calculates the integral of the overlap between the point charges prbCharge and refCharge at coordinates prbCoor and refCoor via Monte Carlo numeric integration (up to 10 Angstrom away from .
    :param prbCoor: 2D array of coordinates of the probe molecule. 
    :param refCoor: 2D array of coordinates of the reference molecule. 
    :param prbCharge: 1D array of partial charges of the probe molecule. 
    :param refCharge: 1D array of partial charges of the reference molecule.
    :param metric: Metric of similarity score.
    :param marginMC: (optional) Margin up to which to integrate (added to coordinates plus/minus their vdW radii).
    :param nMC: (optional) Number of grid points per 1 Angstrom**3 volume of integration vox.
    :param randomseed: (optional) seed for the random number generator
    :return: Similarity of the overlap integrals.
    """
    np.random.seed(randomseed)
    margin=marginMC
    allCoor=np.concatenate((prbCoor,refCoor))
    allVdw=np.concatenate((prbVdw,refVdw))
    lenAll=allCoor.shape[0]
    lenPrb=prbCoor.shape[0]
    lenRef=refCoor.shape[0]
    minValues=np.min(allCoor-allVdw-np.array([[margin]]),axis=0)
    maxValues=np.min(allCoor+allVdw+np.array([[margin]]),axis=0)
    boxvolume=np.prod(maxValues-minValues)
    N=int(boxvolume*nMC)

    nInMargin=0
    intPrbPrb=0
    intPrbRef=0
    intRefRef=0
    
    for i in range(N):
        x=np.random.uniform(minValues[0],maxValues[0])
        y=np.random.uniform(minValues[1],maxValues[1])
        z=np.random.uniform(minValues[2],maxValues[2])

        distPrb=scipy.spatial.distance.cdist(np.array([[x,y,z]]),prbCoor)
        distRef=scipy.spatial.distance.cdist(np.array([[x,y,z]]),refCoor)
        distAll=np.concatenate((distPrb,distRef),axis=1)
        distMinVdw=distAll-allVdw
        minDist=np.min(distMinVdw)
        if minDist<=margin and minDist>0:
            nInMargin+=1
            fPrb=sum([prbCharge[ii]/distPrb[0,ii] for ii in range(lenPrb)])
            fRef=sum([refCharge[ii]/distRef[0,ii] for ii in range(lenRef)])
            intPrbPrb+=fPrb*fPrb
            intPrbRef+=fPrb*fRef
            intRefRef+=fRef*fRef
    factor=nInMargin/N*boxvolume/N
    intPrbPrb*=factor
    intPrbRef*=factor
    intRefRef*=factor

    similarity=SimilarityMetric(intPrbPrb,intRefRef,intPrbRef,metric)
        
    return similarity


def ConstrainedEmbedMultipleConfs(mol,
                                  core,
                                  numConfs = 10,
                                  useTethers=True,
                                  coreConfId = -1,
                                  randomSeed = 2342,
                                  getForceField = UFFGetMoleculeForceField,
                                  **kwargs,
):
    """
    Function to obtain multiple constrained embeddings per molecule. This was taken as is from:
    from https://github.com/rdkit/rdkit/issues/3266
    :param mol: RDKit molecule object to be embedded.
    :param core: RDKit molecule object of the core used as constrained. Needs to hold at least one conformer coordinates.
    :param numCons: Number of conformations to create
    :param useTethers: (optional) boolean whether to pull embedded atoms to core coordinates, see rdkit.Chem.AllChem.ConstrainedEmbed
    :param coreConfId: (optional) id of the core conformation to use
    :param randomSeed: (optional) seed for the random number generator
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

    cids = EmbedMultipleConfs(mol, numConfs=numConfs, coordMap=coordMap, randomSeed=randomSeed, **kwargs)
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


def EmbedAlignConstrainedScore(prbMol,
                               refMols,
                               core,
                               prbNumConfs = 10,
                               refNumConfs = 10,
                               prbCharge = [],
                               refCharges = [],
                               metric = "carbo",
                               integrate = "gauss",
                               partialCharges = "gasteiger",
                               renormalize = False,
                               customrange = None,
                               marginMC = 10,
                               nMC = 1,
                               basisPsi4 = '3-21G',
                               methodPsi4 = 'scf',
                               gridPsi4 = 1,
                               getBestESP = False,
                               randomseed = 2342):
    """Calculates a constrained alignment based on a common pattern in the input molecules. Caution: Will fail if the pattern does not match. 
    Calculates a shape and electrostatic potential similarity of the best alignment.

    :param prbMol: RDKit molecule for which shape and electrostatic similarities are calculated.
    :param refMol: RDKit molecule or list of RDKit molecules serving as references.
    :param core: Common pattern for the constrained embedding as embedded RDKit molecule
    :param prbNumConfs: Number of conformers to create for the probe molecule. A higher number creates better alignments but slows down the algorithm.
    :param refNumConfs: Number of conformers to create for each reference molecule. A higher number creates better alignments but slows down the algorithm.
    :param prbCharge: (optional) List or array of partial charges of the probe molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param refCharge: (optional) List of list or 2D array of partial charges of the reference molecules. If not given, RDKit Gasteiger Charges are used as default.
    :param metric:  (optional) Similarity metric.
    :param integrate: (optional) Integration method.
    :param partialCharges: (optional) Partial charge distribution.
    :param renormalize: (optional) Boolean whether to renormalize the similarity score to [0:1].
    :param customrange: (optional) Custom range to renormalize to, supply as tuple or list of two values (lower bound, upper bound).
    :param marginMC: (optional) Margin up to which to integrate (added to coordinates plus/minus their vdW radii) if MC integration is utilized.
    :param nMC: (optional) Number of grid points per 1 Angstrom**3 volume of integration vox if MC integration is utilized.
    :param basisPsi4: (optional) Basis set for Psi4 calculation.
    :param methodPsi4: (optional) Method for Psi4 calculation.
    :param gridPsi4: (optional) Integer grid point density for ESP evaluation for Psi4 calculation.
    :param getBestESP: (optional) Whether to select best alignment via ESP instead of shape.
    :param randomseed: (optional) seed for the random number generator
    :return: shape similarity and ESP similarity.
    """
    
    if type(refMols) != list:
        refMols=[refMols]

    if refCharges == []:
        refCharges=[[]]*len(refMols)
        
    prbMol=ConstrainedEmbedMultipleConfs(prbMol, core, numConfs=prbNumConfs, randomSeed=randomseed)
    for refMol in refMols:
        refMol=ConstrainedEmbedMultipleConfs(refMol, core, numConfs=refNumConfs, randomSeed=randomseed)
        
    prbMatch = prbMol.GetSubstructMatch(core)
    allShapeSim = []
    allEspSim = []

    if not getBestESP:
        for idx,refMol in enumerate(refMols):
            shapeSim=0
            prbBestConf=0
            refBestConf=0
            refMatch = refMol.GetSubstructMatch(core)
            for i in range(refNumConfs):
                for j in range(prbNumConfs):
                    AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=j,refCid=i)
                    shape = GetShapeSim(prbMol,refMol,j,i)
                    if shape>shapeSim:
                        shapeSim=shape
                        prbBestConf=j
                        refBestConf=i
            #Go back to best alignment
            AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=prbBestConf,refCid=refBestConf)
        
            espSim=GetEspSim(prbMol,refMol,prbBestConf,refBestConf,prbCharge,refCharges[idx],metric,integrate,
                             partialCharges,renormalize,customrange,marginMC,nMC,basisPsi4,methodPsi4,gridPsi4,
                             randomseed=randomseed)
            allShapeSim.append(shapeSim)
            allEspSim.append(espSim)
    else:
        for idx,refMol in enumerate(refMols):
            espSim=0
            shapeSim=0
            refMatch = refMol.GetSubstructMatch(core)
            for i in range(refNumConfs):
                for j in range(prbNumConfs):
                    AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=j,refCid=i)
                    score = GetEspSim(prbMol,refMol,j,i,prbCharge,refCharges[idx],metric,integrate,
                                      partialCharges,renormalize,customrange,marginMC,nMC,basisPsi4,methodPsi4,
                                      gridPsi4, randomseed=randomseed)
                    if score>espSim:
                        espSim=score
                    shape = GetShapeSim(prbMol,refMol,j,i)
                    if shape>shapeSim:
                        shapeSim=shape
            allShapeSim.append(shapeSim)
            allEspSim.append(espSim)
            
    return allShapeSim,allEspSim


def EmbedAlignScore(prbMol,
                    refMols,
                    prbNumConfs = 10,
                    refNumConfs = 10,
                    prbCharge = [],
                    refCharges = [],
                    metric = "carbo",
                    integrate = "gauss",
                    partialCharges = "gasteiger",
                    renormalize = False,
                    customrange = None,
                    marginMC = 10,
                    nMC = 1,
                    basisPsi4 = '3-21G',
                    methodPsi4 = 'scf',
                    gridPsi4 = 1,
                    getBestESP = False,
                    randomseed = 2342):
    """Calculates a general alignment in the input molecules.
    Calculates a shape and electrostatic potential similarity of the best alignment.

    :param prbMol: RDKit molecule for which shape and electrostatic similarities are calculated.
    :param refMol: RDKit molecule or list of RDKit molecules serving as references.
    :param prbNumConfs: Number of conformers to create for the probe molecule. A higher number creates better alignments but slows down the algorithm.
    :param refNumConfs: Number of conformers to create for each reference molecule. A higher number creates better alignments but slows down the algorithm.
    :param prbCharge: (optional) List or array of partial charges of the probe molecule. If not given, RDKit Gasteiger Charges are used as default.
    :param refCharge: (optional) List of list or 2D array of partial charges of the reference molecules. If not given, RDKit Gasteiger Charges are used as default.
    :param metric:  (optional) Similarity metric.
    :param integrate: (optional) Integration method.
    :param partialCharges: (optional) Partial charge distribution.
    :param renormalize: (optional) Boolean whether to renormalize the similarity score to [0:1].
    :param customrange: (optional) Custom range to renormalize to, supply as tuple or list of two values (lower bound, upper bound).
    :param marginMC: (optional) Margin up to which to integrate (added to coordinates plus/minus their vdW radii) if MC integration is utilized.
    :param nMC: (optional) Number of grid points per 1 Angstrom**3 volume of integration vox if MC integration is utilized.
    :param basisPsi4: (optional) Basis set for Psi4 calculation.
    :param methodPsi4: (optional) Method for Psi4 calculation.
    :param gridPsi4: (optional) Integer grid point density for ESP evaluation for Psi4 calculation.
    :param getBestESP: Whether to select best alignment via ESP instead of shape.
    :param randomseed: (optional) seed for the random number generator
    :return: shape similarity and ESP similarity.
    """
    
    if type(refMols) != list:
        refMols=[refMols]

    if refCharges == []:
        refCharges=[[]]*len(refMols)
        
    AllChem.EmbedMultipleConfs(prbMol, prbNumConfs, randomSeed=randomseed)
    for refMol in refMols:
        AllChem.EmbedMultipleConfs(refMol, refNumConfs, randomSeed=randomseed)

    prbCrippen = rdMolDescriptors._CalcCrippenContribs(prbMol)

    allShapeSim = []
    allEspSim = []

    if not getBestESP:
        for idx,refMol in enumerate(refMols):
            shapeSim=0
            prbBestConf=0
            refBestConf=0
            refCrippen = rdMolDescriptors._CalcCrippenContribs(refMol)
            for i in range(refNumConfs):
                for j in range(prbNumConfs):
                    alignment = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippen, refCrippen, j, i)
                    alignment.Align()
                    shape = GetShapeSim(prbMol,refMol,j,i)
                    if shape>shapeSim:
                        shapeSim=shape
                        prbBestConf=j
                        refBestConf=i
            #Go back to best alignment
            alignment = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippen, refCrippen, prbBestConf, refBestConf)
            alignment.Align()
        
            espSim=GetEspSim(prbMol,refMol,prbBestConf,refBestConf,prbCharge,refCharges[idx],metric,integrate,
                             partialCharges,renormalize,customrange,marginMC,nMC,basisPsi4,methodPsi4,gridPsi4,
                             randomseed=randomseed)
            allShapeSim.append(shapeSim)
            allEspSim.append(espSim)
    else:
         for idx,refMol in enumerate(refMols):
            espSim=0
            shapeSim=0
            prbBestConf=0
            refBestConf=0
            refCrippen = rdMolDescriptors._CalcCrippenContribs(refMol)
            for i in range(refNumConfs):
                for j in range(prbNumConfs):
                    alignment = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippen, refCrippen, j, i)
                    alignment.Align()
                    score = GetEspSim(prbMol,refMol,j,i,prbCharge,refCharges[idx],metric,integrate,
                                      partialCharges,renormalize,customrange,marginMC,nMC,basisPsi4,
                                      methodPsi4,gridPsi4,randomseed=randomseed)
                    if score>espSim:
                        espSim=score
                    shape = GetShapeSim(prbMol,refMol,j,i)
                    if shape>shapeSim:
                        shapeSim=shape
            allShapeSim.append(shapeSim)
            allEspSim.append(espSim)           

    return allShapeSim,allEspSim
