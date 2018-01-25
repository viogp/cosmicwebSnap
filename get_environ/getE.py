import numpy as np
from readVweb import readVweb

def sepregion(lamda, thres):
    indc = np.zeros(lamda.shape[0], dtype='int32') - 1
    idck = np.where(lamda[:, 2] >= thres)[0]
    indc[idck] = 3
    Nknots = len(idck)
    idck = np.where((lamda[:, 2] < thres) & (lamda[:, 1] >= thres))[0]
    indc[idck] = 2
    Nfilms = len(idck)
    idck = np.where((lamda[:, 1] < thres) & (lamda[:, 0] >= thres))[0]
    indc[idck] = 1
    Nwalls = len(idck)
    idck = np.where(lamda[:, 0] < thres)[0]
    indc[idck] = 0
    Nviods = len(idck)
    return(indc, Nknots, Nfilms, Nwalls, Nviods)


def identify(pos,vweb):
    N = np.round(vweb.size**(1./3.))
    Nn = np.int32(N)
    xyz = np.int32(pos/(500./N))
    ids = xyz==128
    xyz[ids] = 0
    xyz = xyz[:,0]+Nn*xyz[:,1]+Nn*Nn*xyz[:,2]
    return vweb[xyz]


def return_env(Pos, threshold=0.1, webf="../VPweb_data/VP_039_DM_MHD.000256.Vweb", Vweb=True):
    """
    Parameters:
    ----------
    Pos: must be 3D data array in shape of (N,3) in units of Mpc/h
    threshold: the \lambda_th for classifying environments, default 0.1
    webf: the file location for Vweb output, default: ../VPweb_data/VP_039_DM_MHD.000256.Vweb
    Vweb: If you want to use Vweb method, Set False for the Pweb method.
    
    Return:
    An 1D array indicates the position environments: 0 Void; 1 sheet; 2 filament; 3 Knots.
    """
    
    print("Loading file: ", webf)
    if Vweb:
        Evp = readVweb(webf, selected=[10,11,12])
    else:
        Evp = readVweb(webf, DWEB=True, selected=[22,23,24])

    print("Preparing data...")
    web = sepregion(Evp,threshold)
    

    return np.int32(identify(pos,web[0]))
