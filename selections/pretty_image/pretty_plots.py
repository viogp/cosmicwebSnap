#!/bin/env python

import os.path, sys
import numpy as np
import h5py
import subprocess
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distinct_colours import get_distinct

testing = True 
#------------------------------------------
model = 'gp19/'

path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
plotdir =  path+'plots/'+model+'pretty/'
ndpath = path+'selections/'

surveys = ['DEEP2','DESI'] ; sn = 39 # 0.988
#surveys = ['VVDS-DEEP','eBOSS-SGC'] ; sn = 41 # 0.83

nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr']
marks = ['o','*'] ; cols = ['darkred','dodgerblue']

# Plot only a section of the box
vols = 512
zlow = 10. ; zup  = 20.


lmin = 0. ; lmax= 500.
#lmin = 420. ; lmax= 470.

verbose = True ; dm = False
#------------------------------------------
if testing:
    verbose = True
    vols = 3

    surveys = ['DEEP2'] ; nds = ['-4.2','-2.0'] ; sn = 39
    lmin = 0. ; lmax= 500.

    #surveys = ['eBOSS-SGC'] ; nds = ['-2.0'] ; sn = 41
    #lmin = 200. ; lmax= 300.   
#------------------------------------------

def density_image(x, y, min_count, max_count, shape, log=True, xyrange=None):    
    """
    Author: John Helly
    Contributions: Violeta Gonzalez-Perez

    Generate a floating point image of a 2D particle distribution

    x, y                 - particle coordinates
    min_count, max_count - range of counts that map to 0 and 1 in
                           the output image
    shape                - shape of the output array, (nx,ny)
    log                  - whether to log scale the counts
    range                - range of coordinates to plot

    Returns:

    2D image array with values in range 0-1
    """

    if xyrange is None:
        xyrange = ((np.amin(x), np.amax(x)), (np.amin(y), np.amax(y)))

    image,xedge,yedge = np.histogram2d(x, y, range=xyrange, bins=shape)
    image = np.asarray(image, dtype=float)

    if log:
        # Log scaling - negative values get set to zero in output
        ind = image > 0
        image[np.logical_not(ind)] = 0.0 
        lgmin = np.log10(min_count)
        lgmax = np.log10(max_count)
        image[ind] = (np.log10(image[ind]) - lgmin) / (lgmax-lgmin)
    else:
        # Linear scaling - clamp to range 0-1
        image = (image-min_count)-(max_count-min_count)
        image[image<0.0] = 0.0
        image[image>1.0] = 1.0

    return np.transpose(image)

def extract_particles(isnap, zlow, zup, lmin, lmax, cm):
    """
    Author: John Helly
    Contributions: Violeta Gonzalez-Perez
    """
    basename = "/cosma5/data/jch/MillGas/dm/500/snapdir_%03d/500_dm_%03d" % (isnap, isnap)

    shift = lmin

    # Read in particle data in a 10Mpc/h slice
    # First loop over snapshot files.
    x = [] ; y = []
    for ifile in range(vols):
        # Read coordinates from one file
        fname = "%s.%d.hdf5" % (basename,ifile)
        print "Reading: %s" % fname
        f = h5py.File(fname, "r")
        pos_x = f["PartType1/Coordinates"].value[:,0]
        pos_y = f["PartType1/Coordinates"].value[:,1]
        pos_z = f["PartType1/Coordinates"].value[:,2]
        f.close()

        # Discard particles not in z slice and outside lmin,lmac
        ind = np.where((pos_z>zlow) & (pos_z<zup)  & \
                        (pos_x>lmin) & (pos_x<lmax) & \
                        (pos_y>lmin) & (pos_y<lmax) )

        if (np.shape(ind)[1]>1):
            x = np.append(x,pos_x[ind]) 
            y = np.append(y,pos_y[ind])

    if (len(x) < 2):
        print('STOP: No dark matter particles selected')
        sys.exit()

    # Make image as 2D float array
    img = density_image(x-shift,y-shift,min_count=0.1, max_count=10000, 
                        shape=(1024,1024), log=True)    

    # Display
    return plt.imshow(img, extent=(lmin-shift, lmax-shift, \
                                       lmin-shift, lmax-shift), \
                          cmap=cm, vmin=0.0, vmax=1.0, \
                          interpolation="nearest", origin="lower")

def check_jump(infile,verbose=False):
    # Check the existance of the file
    if (not os.path.isfile(infile)):
        if verbose:
            print('WARNING: {} not found'.format(infile))
            return True
    # Jump files with only a one line header
    wcl_line = subprocess.check_output(["wc", "-l",infile])
    wcl = int(wcl_line.split()[0])
    if (wcl <= 1):
        if verbose:
            print('WARNING: {} has too few lines'.format(infile))
            return True
    return False


###############################
# DM
if dm:
    extract_particles(sn, zlow, zup, lmin, lmax, plt.cm.Greys)
    plt.xlabel("x(Mpc $h^{-1})$") ; plt.ylabel("y(Mpc $h^{-1})$")
    plt.xlim((0.,lmax-lmin)) ;plt.ylim((0.,lmax-lmin))
###############################

val = 1./8.

# Loop over different files
for survey in surveys:
    for nd in nds:
        doplot = False

        # ELG
        for ic,cut in enumerate(cuts):
            infile = ndpath+model+'ascii_files/'+cut+\
                     'cut_'+survey+'_nd'+nd+'_sn'+str(sn)+'.dat'
            jump = check_jump(infile,verbose=verbose)
            if not jump:
                #print('   File: {}'.format(infile))
                xgal,ygal,zgal,lo2 = np.loadtxt(
                    infile,usecols=(0,1,2,10),unpack=True)
                # xgal,ygal,zgal (Mpc/h)
                # lo2: L_tot_OII_ext (10^40 h^-2 erg/s)
                ind = np.where((xgal>lmin) & (xgal<lmax) & 
                               (ygal>lmin) & (ygal<lmax) & 
                               (zgal>zlow) & (zgal<zup))
                if (np.shape(ind)[1]<1):
                    print('No gal. in plot range for {}'.format(infile))
                else:
                    doplot = True
                    # Shifted positions
                    x = xgal[ind] - lmin
                    y = ygal[ind] - lmin
                    # Luminosity for the size of the marker
                    size = np.zeros(shape=len(x)) ; size.fill(1./val)
                    nlo2 = lo2[ind]
                    jnd = np.where(nlo2>0.)
                    size[jnd] = (np.log10(nlo2[jnd]) + 40.)*val 
                    # Plot
                    for i in range(len(size)):
                        plt.plot(x[i],y[i],marks[ic],
                                 c=cols[1],markeredgecolor=cols[1],
                                 markersize=1,alpha=0.8)
                                 #markersize=size[i],alpha=0.8)

        # All
        for ic,cut in enumerate(cuts):
            infile = ndpath+model+'ascii_files/'+cut+\
                     'cut_All_nd'+nd+'_sn'+str(sn)+'.dat'
            jump = check_jump(infile,verbose=verbose)
            if not jump:
                xgal,ygal,zgal,lo2 = np.loadtxt(
                    infile,usecols=(0,1,2,10),unpack=True)
                ind = np.where((xgal>lmin) & (xgal<lmax) & 
                               (ygal>lmin) & (ygal<lmax) & 
                               (zgal>zlow) & (zgal<zup))
                if (np.shape(ind)[1]<1):
                    print('No gal. in plot range for {}'.format(infile))
                else:
                    doplot = True
                    x = xgal[ind] - lmin
                    y = ygal[ind] - lmin
                    size = np.zeros(shape=len(x)) ; size.fill(1./val)
                    nlo2 = lo2[ind]
                    jnd = np.where(nlo2>0.)
                    size[jnd] = (np.log10(nlo2[jnd]) + 40.)*val 
                    # Plot
                    for i in range(len(size)):
                        plt.plot(x[i],y[i],marks[ic],
                                 c='none',markeredgecolor=cols[0],
                                 markersize=1,alpha=0.8)
                                 #markersize=size[i],alpha=0.8)

        if doplot:
            # Save figures
            plotfile = plotdir+'pretty_'+survey+\
                       '_nd'+nd+'_sn'+str(sn)+'.pdf'
            
            plt.savefig(plotfile)
            print '* Output: ',plotfile
