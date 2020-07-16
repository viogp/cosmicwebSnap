# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

# Functions

def add2sat(pos,gtype,mult=1):
    '''
    Add (mult=1) or remove (mult=-1) the position of the central to the satellites
    '''
    last = pos[0]
    for i in range(len(pos)):
        if gtype[i]:
            #here: gtype working as expected?

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

nbin, mmin, mmax = 70.,9.,16.

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

    # Find the bin each halo mass belongs to
    imass = ((np.copy(mhhalo)-mmin)/(mmax-mmin)*nbin).astype(int)
