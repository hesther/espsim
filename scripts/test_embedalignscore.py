from espsim import EmbedAlignConstrainedScore
from rdkit import Chem
from rdkit.Chem import AllChem


#Load molecules (9 reference molecules, 3 probe molecules)
train_smiles=['C1=CC=C(C=C1)C(C(=O)O)O','CCC(C(=O)O)O','OC(C(O)=O)c1ccc(Cl)cc1','C1=CC(=CC=C1C(C(=O)O)O)O','COc1ccc(cc1)C(O)C(O)=O','OC(C(O)=O)c1ccc(cc1)[N+]([O-])=O','CCCC(C(=O)O)O','CCC(C)C(C(=O)O)O','CC(C(=O)O)O']
test_smiles=['c1c(C)cc(cc1)C(O)C(O)=O','Cc1ccc(cc1)C(O)C(O)=O','C(C(C(=O)O)O)O']
patt=Chem.MolFromSmiles("[H]OC([H])(C)C(=O)O[H]",sanitize=False)
train_mols=[Chem.AddHs(Chem.MolFromSmiles(x)) for x in train_smiles]
test_mols=[Chem.AddHs(Chem.MolFromSmiles(x)) for x in test_smiles]
print("Created reference molecules:")
for train_smi in train_smiles:
    print(train_smi)
print("#"*30)
print("Created probe molecules:")
for test_smi in test_smiles:
    print(test_smi)
print("#"*30)


#Embed first training molecule and optimize structure to extract the 3D pattern to get a 3D core molecule:
AllChem.EmbedMolecule(train_mols[0],AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(train_mols[0])
core = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(train_mols[0],patt),Chem.MolFromSmiles('*'))
core.UpdatePropertyCache()
train_mols[0].RemoveAllConformers()


#For each probe molecule, call EmbedAlignConstrainedScore() and print all scores.    
for j,test_mol in enumerate(test_mols):
    results_shape,results_electrostatics=EmbedAlignConstrainedScore(test_mol,train_mols,core)
    print("Results for probe molecule",test_smiles[j])
    print('%40s %8s %8s' % ("Reference","Shape","ESP"))
    for i in range(len(results_shape)):
        print('%40s %8.2f %8.2f' % (train_smiles[i],results_shape[i],results_electrostatics[i]))

