# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

# Functions

def add2sat(xgal,ygal,zgal,gtype,add=True):
    '''
    Add or remove the position of the central to the satellites
    '''
    xlast = xgal[0]
    ylast = ygal[0]
    zlast = zgal[0]

    for i in range(len(xgal)):
        if gtype[i]:
            if (add):
                xgal[i] = xgal[i] + xlast
                ygal[i] = ygal[i] + ylast
                zgal[i] = zgal[i] + zlast
            else:
                xgal[i] = xgal[i] - xlast
                ygal[i] = ygal[i] - ylast
                zgal[i] = zgal[i] - zlast
        else:
            xlast = xgal[0]
            ylast = ygal[0]
            zlast = zgal[0]


def correct_periodic(xgal,ygal,zgal,gtype,lbox,halfbox=False):
    '''
    Correct for peridic boundary conditions
    '''
    lbox2 = lbox/2.
    
    for arr in [xgal,ygal,zgal]:
        if halfbox:
            ind = np.where((arr > lbox2) & (gtype > 0))
            arr[ind] = arr[ind] - lbox

            ind = np.where((arr < -lbox2) & (gtype > 0))
            arr[ind] = arr[ind] + lbox
        else:
            ind = np.where((arr >= lbox) & (gtype > 0))
            arr[ind] = arr[ind] - lbox

            ind = np.where((arr < 0) & (gtype > 0))
            arr[ind] = arr[ind] + lbox


def mhhalo2sat(mhhalo,gtype):
    '''
    Assign the co 
    '''
##here
############################################

Testing = True

nvol = 64
sn_list = ['39','41']

if Testing:
    nvol = 2
    sn_list = ['39']

############################################

path = '/cosma6/data/dp004/dc-gonz3/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

############################################

nbin = 70
mmin, mmax = 9., 16.

for iz, sn in enumerate(sn_list):
    volume = 0. ; iifil = -1
    for ivol in range(nvol):
        gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (not os.path.isfile(gfile)):
            continue
        iifil += 1
        if Testing: 
            if iifil>2: break

        # Read the halo information
        f = h5py.File(gfile,'r')
        group = f['Output001']
        vol1 = f['Parameters/volume'][()] ; volume = volume + vol1
        tjm = group['Trees/jm'][:]
        tngals = group['Trees/ngals'][:] 

        ixgal   = group['xgal'][:]   # Mpc/h
        iygal   = group['ygal'][:]
        izgal   = group['zgal'][:]
        imhhalo = np.log10(group['mhhalo'][:])   # log10(M/Msun/h)
        igtype  = group['type'][:] # 0= Centrals; 1,2= Satellites
        ijm     = np.repeat(tjm,tngals) 
        iihhalo = group['ihhalo'][:] 
        
        f.close()

        # Array with subvolumes
        ivols = np.zeros(shape=len(ixgal),dtype=int) ; ivols.fill(ivol)

        # Sort haloes within one subvolume
        ind = np.lexsort((igtype,iihhalo,ijm)) # Sort by jm, then by ihhalo, then by gtype

        if (iifil==0):
            xgal   = ixgal[ind]
            ygal   = iygal[ind]
            zgal   = izgal[ind]
            mhhalo = imhhalo[ind]
            gtype  = igtype[ind]
            jm     = ijm[ind]
            ihhalo = iihhalo[ind]
            vols   = ivols
        else:
            xgal   = np.append(xgal,ixgal[ind])
            ygal   = np.append(ygal,iygal[ind])
            zgal   = np.append(zgal,izgal[ind])
            mhhalo = np.append(mhhalo,imhhalo[ind])
            gtype  = np.append(gtype,igtype[ind])
            jm     = np.append(jm,ijm[ind])
            ihhalo = np.append(ihhalo,iihhalo[ind])
            vols   = np.append(vols,ivols)

    lbox = pow(volume,1./3.)
    print('sn={}, Box side (Mpc/h) ={}'.format(sn,lbox))

    # Initialise arrays to contain shuffled properties
    s_xgal = np.copy(xgal)
    s_ygal = np.copy(ygal)
    s_zgal = np.copy(zgal)
    s_mhhalo = np.copy(mhhalo)

    # Find the bin each halo mass belongs to
    imass = ((np.copy(mhhalo)-mmin)/(mmax-mmin)*nbin).astype(int)

    # Change the position of the satellites to be relative to the halo
    add2sat(xgal,ygal,zgal,gtype,add=False)

    # Correct for periocid boundary conditions (for satellites)
    correct_periodic(xgal,ygal,zgal,gtype,lbox,halfbox=False) 

    # Shuffle the centrals in each halo mass bin
    for i in range(nbin):
        # Select all central galaxies in the halo mass bin
        ind = np.where((gtype == 0) & (i == imass))
        if (np.shape(ind)[1] < 1): continue

        # Get shuffled indexes
        ind_shuffle = np.copy(ind)
        np.random.shuffle(ind_shuffle)

        # Reassigned properties
        s_xgal[ind] = xgal[ind_shuffle]
        s_ygal[ind] = ygal[ind_shuffle]
        s_zgal[ind] = zgal[ind_shuffle]
        s_mhhalo[ind] = mhhalo[ind_shuffle]

    # Add the central position
    add2sat(s_xgal,s_ygal,s_zgal,gtype,add=True)

    # Change the mass of the host halo for the satellites
    mhhalo2sat(s_mhhalo,gtype)

