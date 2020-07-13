# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

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
        if (os.path.isfile(gfile)):        
            iifil += 1
            if Testing: 
                if iifil>2: break

            # Read the halo information
            f = h5py.File(gfile,'r')
            group = f['Output001']
            vol1 = f['Parameters/volume'][()] ; volume = volume + vol1
            tjm = group['Trees/jm'][:]
            tngals = group['Trees/ngals'][:] 
            jm = np.repeat(tjm,tngals) ;print(gfile)
            ihhalo = group['ihhalo'][:] 
            haloid = 10e10*ivol + 10e6*jm + ihhalo
            print(haloid,np.shape(haloid),type(haloid)) ; sys.exit()
            ###here Is this a good id for the shuffling?
            if (iifil==0):
                xgal   = group['xgal'][:]   # Mpc/h
                ygal   = group['ygal'][:]
                zgal   = group['zgal'][:]
                mhhalo = group['mhhalo'][:]   # Msun/h
                gtype  = group['type'][:] # 0= Centrals; 1,2= Satellites
            else:
                xgal   = np.append(xgal,group['xgal'][:])   # Mpc/h
                ygal   = np.append(ygal,group['ygal'][:])
                zgal   = np.append(zgal,group['zgal'][:])
                mhhalo = np.append(mhhalo,group['mhhalo'][:])   # Msun/h
                gtype  = np.append(gtype,group['type'][:]) # 0= Centrals; 1,2= Satellites
            f.close()
    
            
    lbox = pow(volume,1./3.)
    print('sn={}, Box side (Mpc/h) ={}'.format(sn,lbox))
