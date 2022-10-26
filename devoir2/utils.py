import numpy as np

def nL1(v1, v2):
    """ 
    Norme L1 
    v1 : premier vecteur
    v2 : deuxième vecteur
    """
    return np.sum(np.abs(v1 - v2))/len(v1)

def nL2(v1, v2):
    """ 
    Norme L2 
    v1 : premier vecteur
    v2 : deuxième vecteur
    """
    return np.sqrt(np.sum((v1 - v2)**2))/len(v1)

def nLi(v1, v2):
    """ 
    Norme Linf 
    v1 : premier vecteur
    v2 : deuxième vecteur
    """
    return np.max(np.abs(v1 - v2))
