import sys
from morgan_24_06_24 import *
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit import DataStructs
from rdkit import Chem
import warnings
import os
#from bioster_replacement import *

warnings.filterwarnings("ignore")
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)  # Only show errors, not warnings.

def search_bioisosters_morgan(PDB_ligands_file_name, replacement_file_name, deep):
    
    ''' this function determines the mol files that correspond
       to the same molecules; it opens 2 sdf files, counts the 
       Morgan coefficient (or fingerprint) for each mol from the sdf,
       and then moves between two sdf files. This function uses
       the local 'morgan_types' library '''
    
    #inf = open("C:/Users/moise/Desktop/projects/lab44/ITOG/parser/replacements/Components-pub.sdf",'rb')
    Fingerprints_pdb = [] # list with original PDB fingerprints
    inf = open(PDB_ligands_file_name,'rb') # paste the path to the required folder from your local directory
    suppl = Chem.ForwardSDMolSupplier(inf)
    for mol in suppl: # iterating over the mol from the sdf file
        if mol != None:
            name = mol.GetProp('_Name') # get mol name
            try:
                morgan_mol = MorganTypes(mol) # create MorganTypes object 
                atom_weights = morgan_mol.calc_morgan_weight(k=deep) # calculate Morgan's weights for all atoms in mol 
                wmp = 0 # calculate Morgan coefficient (or fingerprint)
                for i in atom_weights.keys():
                    wmp += atom_weights[i] 
            except: 
                wmp = AllChem.GetMorganFingerprint(mol, 2)
            Fingerprints_pdb.append([wmp, name]) # add to list

    # repeat algorithm for the second sdf file 
    inf = open(replacement_file_name,'rb') # paste the path to the required folder from your local directory
    Fingerprints_bioster = [] # list with bioisosters fingerprints
    f1suppl = Chem.ForwardSDMolSupplier(inf)
    for mol in f1suppl:
        if mol != None:
            name = mol.GetProp('_Name')
            try:
                morgan_mol = MorganTypes(mol)
                atom_weights = morgan_mol.calc_morgan_weight(k=deep)
                wmp = 0
                for i in atom_weights.keys():
                    wmp += atom_weights[i]
            except: 
                wmp = AllChem.GetMorganFingerprint(mol, 2)
            Fingerprints_bioster.append([wmp, name])
            
    inf.close()

    # we compare the lists and find the same molecules
    # we write down the names of identical molecules in the list
    biosters = []
    for fp in Fingerprints_bioster:
        name = fp[1]
        fp1 = fp[0]

        for ffpp in Fingerprints_pdb:
            fp2 = ffpp[0]
            try:
                if fp1 == fp2 and [name, ffpp[1]] not in biosters:
                    biosters.append([name, ffpp[1]])
            except:
                # for reliability, check the identity using the RDKit function
                # this is necessary in order not to lose the bioisosters 
                if DataStructs.DiceSimilarity(fp1,fp2) == 1 and [name, ffpp[1]] not in biosters: 
                    biosters.append([name, ffpp[1]])
    
    return biosters


def search_bioisosters_fingerprints(PDB_ligands_filename, replacement_file_name, deep):

    """ This function not in use; this function determines the mol files that correspond
       to the same molecules """
    
    inf = open(PDB_ligands_filename,'rb') # paste the path to the required folder from your local directory
    Fingerprints_pdb = []
    f1suppl = Chem.ForwardSDMolSupplier(inf)
    for mol in f1suppl:
        if mol != None:
            name = mol.GetProp('_Name')
            fp = AllChem.GetMorganFingerprint(mol, deep)
            Fingerprints_pdb.append([fp, name])
    inf.close()
    
    
    inf = open(replacement_file_name,'rb') # paste the path to the required folder from your local directory
    Fingerprints_bioster = []
    f2suppl = Chem.ForwardSDMolSupplier(inf)
    for mol in f2suppl:
        if mol != None:
            name = mol.GetProp('_Name')
            fp = AllChem.GetMorganFingerprint(mol, deep)
            Fingerprints_bioster.append([fp, name])
    inf.close()
 
    #file = open(path, "r")
    #for line in file:
    #    h = line.strip().split()
    #    smi = h[0]
    #    name = h[1]
    #    m1 = Chem.MolFromSmiles(smi)
    #    ffp = AllChem.GetMorganFingerprint(m1,2)
    #    Fingerprints_bioster.append([ffp, name])

    biosters = []
    for fp in Fingerprints_bioster:
        name = fp[1]
        fp1 = fp[0]
        q += 1
        for ffpp in Fingerprints_pdb:
            fp2 = ffpp[0]

            if DataStructs.DiceSimilarity(fp1,fp2) == 1:
                biosters.append([name, ffpp[1]])
    return biosters


