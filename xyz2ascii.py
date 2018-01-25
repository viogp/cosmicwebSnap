# Output the x,y,z coordinates of galaxies selected in OII flux into 
# an ascii file, either in redshift or real-space

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

zspace = True # Set either real or redshift-space

nvol = 64

inleg = ['centrals,20 particles']
ntypes = len(inleg)

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/gp15newmg/iz41/'
#############################
line = 'OII3727' ; lline = '[OII]'

m20 = 1.0 #0.0 #1.87*10.**10
############################################
outpath = path+'ascii_files/'
outfile = outpath+'centrals1.dat'

ngal = 0
volume = 0. ; firstpass = True
for ivol in range(nvol):
    gfile = path+'ivol'+str(ivol)+'/galaxies.hdf5'
    if (os.path.isfile(gfile)):
        # Get some of the model constants
        f = h5py.File(gfile,'r')
        group = f['Parameters']
        vol1 = group['volume'].value ; volume = volume + vol1
        h0 = group['h0'].value ; lambda0 =group['lambda0'].value
        omega0 = group['omega0'].value ; omegab = group['omegab'].value

        group = f['Output001']
        zz     = group['redshift'].value
        set_cosmology(omega0=omega0,omegab=omegab,\
                          lambda0=lambda0,h0=h0,\
                          universe="Flat",\
                          include_radiation=False)
        xgal   = group['xgal'].value   # Mpc/h
        ygal   = group['ygal'].value
        zgal   = group['zgal'].value
        vxgal  = group['vxgal'].value*(1.+zz)/H(zz)  # km/s
        vygal  = group['vygal'].value*(1.+zz)/H(zz)
        vzgal  = group['vzgal'].value*(1.+zz)/H(zz)

        mhhalo = group['mhhalo'].value   # Msun/h
        gtype  = group['type'].value # 0= Centrals; 1,2= Satellites

        f.close()

        for index in range(ntypes):
            ind = np.where((gtype<1) & (mhhalo>m20))
            ngal = ngal + np.shape(ind)[1]
            mass = np.log10(mhhalo[ind])

            xzs = xgal[ind] 

            tofile = np.column_stack((xzs,ygal[ind],\
                                          zgal[ind],\
                                          vxgal[ind],\
                                          vygal[ind],\
                                          vzgal[ind],\
                                          mass))

            with open(outfile,'a') as f_handle:
                np.savetxt(f_handle, tofile, \
                               fmt=('%.5f %.5f %.5f %.5f %.5f %.5f %.5f'))

        # Continue loop over subvolumes
        firstpass = False

lbox = pow(volume,1./3.)
print zz,'Box side (Mpc/h) =',lbox
print 'Output :',outfile
