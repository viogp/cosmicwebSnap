# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    
import h5py
import numpy as np
from Cosmology import *
from Corrfunc.theory.xi import xi # For test plot
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec 
import mpl_style 
plt.style.use(mpl_style.style1)

# Functions

def add2sat(xgal,ygal,zgal,gtype,add=True):
    '''
    Add or remove the position of the central to the satellites
    '''
    xlast = xgal[0]
    ylast = ygal[0]
    zlast = zgal[0]

    for i in range(len(xgal)):
        if (gtype[i] > 0):
            if (add):
                xgal[i] = xgal[i] + xlast 
                ygal[i] = ygal[i] + ylast
                zgal[i] = zgal[i] + zlast
            else:
                xgal[i] = xgal[i] - xlast
                ygal[i] = ygal[i] - ylast
                zgal[i] = zgal[i] - zlast
        else:
            xlast = xgal[i]
            ylast = ygal[i]
            zlast = zgal[i]
        

def correct_periodic(xgal,ygal,zgal,gtype,lbox,halfbox=False):
    '''
    Correct satellites for periodic boundary conditions
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

def mhhalo2sat(mh,gtype):
    '''
    Assign the new halo mass to the satellite galaxies 
    '''
    mhlast = mh[0]

    for i in range(len(mh)):
        if (gtype[i] > 0):
            mh[i] = mhlast
        else:
            mhlast = mh[i]

############################################

Testing = False

nvol = 64
sn_list = ['39','41']

if Testing:
    nvol = 2
    sn_list = ['39']
    lbox_test = 500.

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
        ijm     = np.repeat(tjm,tngals) # tree index 
        iihhalo = group['ihhalo'][:]    # dhalo index within a tree
        
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
    if Testing: lbox=lbox_test
    print('sn={}, Box side (Mpc/h) ={}'.format(sn,lbox))

    # Initialise arrays to contain shuffled properties
    s_xgal = np.copy(xgal)
    s_ygal = np.copy(ygal)
    s_zgal = np.copy(zgal)
    s_mhhalo = np.copy(mhhalo)

    # Change the position of the satellites to be relative to the halo
    add2sat(s_xgal,s_ygal,s_zgal,gtype,add=False)

    # Correct for periocid boundary conditions (for satellites)
    correct_periodic(s_xgal,s_ygal,s_zgal,gtype,lbox,halfbox=True) 

    # Find the bin each halo mass belongs to
    imass = ((np.copy(s_mhhalo)-mmin)/(mmax-mmin)*nbin).astype(int)

    # Shuffle the centrals in each halo mass bin
    for i in range(nbin):
        # Select all central galaxies in the halo mass bin
        ind = np.where((gtype == 0) & (i == imass))[0]
        if (len(ind) < 1): continue

        # Get shuffled indexes
        ind_shuffle = np.copy(ind)
        np.random.shuffle(ind_shuffle)

        # Reassigned properties
        s_xgal[ind] = s_xgal[ind_shuffle]
        s_ygal[ind] = s_ygal[ind_shuffle]
        s_zgal[ind] = s_zgal[ind_shuffle]
        s_mhhalo[ind] = s_mhhalo[ind_shuffle]

    # Add the central position
    add2sat(s_xgal,s_ygal,s_zgal,gtype,add=True)

    # Change the mass of the host halo for the satellites
    mhhalo2sat(s_mhhalo,gtype)

    # Correct satellites again for periodic boundaries
    correct_periodic(s_xgal,s_ygal,s_zgal,gtype,lbox,halfbox=False) 

    # Save the shuffled information into a file
    outff = path+model+'iz'+sn+'/shuffled.hdf5'

    hf = h5py.File(outff,'w')
    hf.create_dataset('vols',data=vols)
    hf.create_dataset('jm',data=jm)
    hf.create_dataset('ihhalo',data=ihhalo)
    hf.create_dataset('gtype',data=gtype)
    hf.create_dataset('s_xgal',data=s_xgal)
    hf.create_dataset('s_ygal',data=s_ygal)
    hf.create_dataset('s_zgal',data=s_zgal)
    hf.create_dataset('s_mhhalo',data=s_mhhalo)
    hf.close()

    print('Output: {}'.format(outff))

    # Plot to test the shuffling
    nd = 3.16e-4   # Number density cut
    s_mhhalo = np.sort(s_mhhalo)[::-1] # Sort in reverse order
    mcut = s_mhhalo[int(nd*lbox**3)]
    print('Mcut for plot ={}'.format(mcut))

    mask = mhhalo>mcut
    x1 = xgal[mask]
    y1 = ygal[mask]
    z1 = zgal[mask]

    mask = s_mhhalo>mcut
    x2 = s_xgal[mask]
    y2 = s_ygal[mask]
    z2 = s_zgal[mask]

    rbins = np.logspace(-2.,2.,40)
    rbins_mid = rbins[:-1] + (rbins[1]-rbins[0])/2.

    nthreads = 4
    xi1_cf = xi(lbox, nthreads, rbins, x1, y1, z1)
    xi2_cf = xi(lbox, nthreads, rbins, x2, y2, z2)
    ratio = xi1_cf['xi']/xi2_cf['xi']
    #print('Clustering ratios = {}'.format(ratio))

    # Plot the clustering as a test
    fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
    gs = gridspec.GridSpec(2, 1)
    gs.update(hspace=0.0)
    ax = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])

    ax.plot(np.log10(rbins_mid), np.log10(xi1_cf['xi']),label='nd='+str(nd))
    ax.plot(np.log10(rbins_mid), np.log10(xi2_cf['xi']),label='Shuffle')
    ax.legend(loc='best',frameon=False)
                                                                                 
    ax2.plot([-3,3], [1,1], color='k',ls=':',linewidth=1)
    ax2.plot(np.log10(rbins_mid), ratio)

    ax.set_xlim([-1.,1.5]) ; ax.set_xlim([-2.,3.5])
    ax2.set_xlim([-1.,1.5]) ; ax2.set_ylim([0.8,1.2])
    ax.set_ylabel(r'$\rm log(\xi(r))$')
    ax2.set_ylabel(r'$\rm log(r/h^{-1}Mpc$')
    ax2.set_ylabel(r'$\rm Ratios$')

    plotff = path+model+'iz'+sn+'/shuffled.pdf'
    plt.savefig(plotff,bbox_inches='tight')
    print('Plot at {}'.format(plotff))
