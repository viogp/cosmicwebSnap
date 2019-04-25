# Output the x,y,z coordinates of galaxies selected in OII flux into 
# an ascii file, either in redshift or real-space

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

nvol = 64

sval = 3.0

sn = '41' ; inleg = ['eBOSS']

#sn = '41' ; inleg = ['VVDS-DEEP']

#sn = '39' ; inleg = ['DEEP2']

#sn = '39' ; inleg = ['DESI']

#############################
path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

line = 'OII3727' ; lline = '[OII]'
############################################
ntypes = len(inleg) + 1

outpath = '/cosma5/data/durham/violeta/lines/cosmicweb/selections/'
for ileg in inleg:
    outm = outpath+model+'ascii_files/scut_'+ileg+'_s'+str(sval)+'_sn'+str(sn)+'.dat'
    print('Output: ',outm)
    outf = open(outm, 'w')
    outf.write('# xgal,ygal,zgal (Mpc/h), vxgal,vygal,vzgal (Km/s), log10(massh),log10(mass/Msun/h), log10(sfr/Msun/h/Gyr), lum,lum_ext (10^40 h^-2 erg/s) \n')
    outf.close() 

volume = 0.
for ivol in range(nvol):
    gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/galaxies.hdf5'
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
        tomag = band_corrected_distance_modulus(zz)

        xgal   = group['xgal'].value   # Mpc/h
        ygal   = group['ygal'].value
        zgal   = group['zgal'].value
        vxgal  = group['vxgal'].value*(1.+zz)/H(zz)  # km/s
        vygal  = group['vygal'].value*(1.+zz)/H(zz)
        vzgal  = group['vzgal'].value*(1.+zz)/H(zz)

        mhhalo = group['mhhalo'].value   # Msun/h
        gtype  = group['type'].value # 0= Centrals; 1,2= Satellites

        mdisk = group['mstars_disk'].value # Msun/h
        mbulge = group['mstars_bulge'].value
        mass1 = mdisk + mbulge

        sdisk = group['mstardot'].value # Msolar/h/Gyr
        sbulge = group['mstardot_burst'].value
        sfr1 = sdisk + sbulge

        f.close()

        gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/elgs.hdf5'
        if (not os.path.isfile(gfile)):        
            print('STOP {} not found'.format(gfile)) ; sys.exit()
        f = h5py.File(gfile,'r')
        group = f['Output001']

        lum = group['L_tot_'+line].value # 10^40 h^-2 erg/s
        lum_ext = group['L_tot_'+line+'_ext'].value 

        for index,ileg in enumerate(inleg):
            seln = ileg.split('_')

            if (seln[0] == 'DEEP2'):
                ielg = 0
            elif (seln[0] == 'VVDS-DEEP'):
                ielg = 1
            elif (seln[0] == 'eBOSS'):
                ielg = 2
            elif (seln[0] == 'DESI'):
                ielg = 3
            #elif (seln[0] == 'VVDS-WIDE'):
            #    ielg = 2

            # Mass selected
            filcut = outpath+model+'sfr_cuts_sn'+sn+'.dat'
            if (not os.path.isfile(filcut)):
                print('STOP: {} not found'.format(filcut)) ; sys.exit()
            ielgs, svals, dum, ngals = np.loadtxt(filcut,unpack=True)
            ind = np.where((ielgs == ielg) & (svals == sval))
            val = ngals[ind][0]

            ind = np.where(sfr1>10**val)

            massh = np.log10(mhhalo[ind])
            mass = np.log10(mass1[ind])
            sfr = np.log10(sfr1[ind])

            tofile = np.column_stack((xgal[ind],\
                                          ygal[ind],\
                                          zgal[ind],\
                                          vxgal[ind],\
                                          vygal[ind],\
                                          vzgal[ind],\
                                          massh,mass,sfr,\
                                          lum[ind],lum_ext[ind]))

            outm = outpath+model+'ascii_files/scut_'+ileg+'_s'+str(sval)+'_sn'+str(sn)+'.dat'
            with open(outm,'a') as outf:
                np.savetxt(outf, tofile, fmt ='%.5f')

        f.close()

lbox = pow(volume,1./3.)
print zz,'Box side (Mpc/h) =',lbox
