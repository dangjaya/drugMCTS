# helpers

# import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

def index2SMILE(atomList,
                eos_token,
                max_len,
                smileVocab):
    smileLookup = {v: k for k, v in smileVocab.stoi.items() if k not in ['<unk>','<pad>','<sos>','<eos>'] }
    length      = atomList.index( eos_token ) if eos_token in atomList else max_len
    smile       = [smileLookup[tok] for tok in atomList[1:length] if tok in smileLookup ]
    return ''.join(smile)


def validSMILES(smile, verbose=False):
    m = Chem.MolFromSmiles(smile,sanitize=False)
    if m is None:
       if verbose:
          print('invalid SMILES')
       return 0
    else:
      try:
        Chem.SanitizeMol(m)
      except:
        if verbose:
           print('invalid chemistry')
        return 0
    return 1


def checkLipinski(smile,verbose=False):
    if verbose:
       print ("Checking lipinski for ", smile )
    molecule = Chem.MolFromSmiles(smile,sanitize=True)
    if molecule is None:
       return 0
    else:
        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)

        if molecular_weight <= 500 and \
           logp <= 5 and \
           h_bond_donor <= 5 and \
           h_bond_acceptors <= 10 and \
           rotatable_bonds <= 5:
            return 1
        else:
            return 0

def checkDrugLikeness(smile,verbose=False):
    if verbose:
       print ("Checking drug-likeness for ", smile )
    molecule = Chem.MolFromSmiles(smile,sanitize=True)
    if molecule is None:
       print('invalid SMILES')
       return 0
    else:
        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)

        #Drug Like Filter
        if molecular_weight < 400 and \
           num_of_rings > 0 and \
           rotatable_bonds < 5 and \
           h_bond_donor <= 5 and \
           h_bond_acceptors <= 10 and \
           logp < 5:
            return 1
        else:
            return 0

def computeQED(smile, verbose=False):
    if verbose:
        print ("Compute QED ", smile )
    molecule = Chem.MolFromSmiles(smile,sanitize=True)
    if molecule is None:
       return 0
    else:
        return Chem.QED.qed( molecule )

