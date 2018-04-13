#! /usr/bin/env python

import os.path, sys
import h5py
from numpy import *
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt

#### John's scripts###########
def density_image(x, y, min_count, max_count, shape, log=True, range=None):    
    """
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
    if range is None:
        range = ((amin(x), amax(x)), (amin(y), amax(y)))
    image,xedge,yedge = histogram2d(x, y, range=range, bins=shape)
    image = asarray(image, dtype=float)
    if log:
        # Log scaling - negative values get set to zero in output
        ind = image > 0
        image[logical_not(ind)] = 0.0        
        lmin = log10(min_count)
        lmax = log10(max_count)
        image[ind] = (log10(image[ind]) - lmin) / (lmax-lmin)
    else:
        # Linear scaling - clamp to range 0-1
        image = (image-min_count)-(max_count-min_count)
        image[image<0.0] = 0.0
        image[image>1.0] = 1.0
        
    return image.transpose()


def plot_density(x, y, min_count, max_count, shape, range=None,
                 cm=None, log=True):
    """
    Calculate and display a density image

    x, y                 - particle coordinates
    min_count, max_count - range of counts that map to 0 and 1 in
                           the output image
    shape                - shape of the output array, (nx,ny)
    log                  - whether to log scale the counts
    range                - range of coordinates to plot
    cm                   - colour map to use

    Returns:

    Image object created by imshow()
    """
    if range is None:
        range = ((amin(x), amax(x)), (amin(y), amax(y)))
    image = density_image(x, y, min_count, max_count, shape, log, range)
    return plt.imshow(image, extent=(range[0][0], range[0][1],
                                     range[1][0], range[1][1]),
                      cmap=cm, vmin=0.0, vmax=1.0,
                      interpolation="nearest", origin="lower")


##############################
plotfile =  '/gpfs/data/violeta/lines/desi_hod_o2/plots/pretty_images/snap61.pdf'
simpath = '/gpfs/data/Galform/Merger_Trees/MillGas/dm/500/new/particle_lists/particle_list_061.'
snap = 61 #37
zlow = 10.
zup  = 20.

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp14/'
line = 'OII3727' ; lline = '[OII]'
##############################

# Look for OII emitters
bands = ['RK','m2','m2']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.]

inleg = ['DEEP2','VVDS-WIDE','VVDS-DEEP']
colors = ['g','k','r']

# Loop over volumes
firstpass = True ; firstelgs = True
for ivol in range(64): #(64):
    pfile = simpath+str(ivol)+'.hdf5'
    if (os.path.isfile(pfile)):
        #print pfile
        f = h5py.File(pfile,'r')
        group = f['particles']
        snapshot = group['snapshotNumber'].value
        inx = group['position'].value[:,0]
        iny = group['position'].value[:,1]
        inz = group['position'].value[:,2]
        f.close()

        ind = where((snapshot == snap) & (inz>zlow) & (inz<zup))
        if firstpass:
            x = inx[ind] ; y = iny[ind] ; z = inz[ind]
            firstpass = False
        else:
            x = append(x,inx[ind])
            y = append(y,iny[ind])
            z = append(z,inz[ind])
    else:
        print pfile,' not found'

    # OII emitters 
    pfile = path+model+'/iz'+str(snap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
    if (os.path.isfile(pfile)):
        #print pfile
        # Get some of the model constants
        f = h5py.File(pfile,'r')
        zz   = f['Output001/redshift'].value
        xgal = f['Output001/xgal'].value
        ygal = f['Output001/ygal'].value
        zgal = f['Output001/zgal'].value
        group = f['Parameters']
        h0 = group['h0'].value 
        omega0 = group['omega0'].value
        omegab = group['omegab'].value
        lambda0 =group['lambda0'].value
        f.close()

        set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
        tomag = band_corrected_distance_modulus(zz)
        
        efile = path+model+'/iz'+str(snap)+'/ivol'+str(ivol)+'/elgs.hdf5'
        if (os.path.isfile(efile)):
            f = h5py.File(efile,'r')            
            for index,ib in enumerate(bands):
                mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                icut = mcuts[index] ; fluxcut = fcuts[index]
                lcut = emission_line_luminosity(fluxcut,zz)
                
                # ELGs
                if (np.shape(mag)==np.shape(zgal) and np.shape(lum_ext)==np.shape(zgal)):
                    ind = np.where((mag<icut) & (lum_ext>lcut)  & \
                                       (zgal>zlow) & (zgal<zup))                
                    if firstelgs:
                        if (index ==0): 
                            xd2 = xgal[ind] ; yd2 = ygal[ind]
                            ld2 = lum_ext[ind]
                        elif (index ==1): 
                            xvw = xgal[ind] ; yvw = ygal[ind] 
                            lvw = lum_ext[ind]
                        else: 
                            xvd = xgal[ind] ; yvd = ygal[ind] 
                            lvd = lum_ext[ind]
                            firstelgs = False
                    else:
                        if (index ==0): 
                            xd2 = append(xd2,xgal[ind])
                            yd2 = append(yd2,ygal[ind])
                            ld2 = append(ld2,lum_ext[ind])
                        elif (index ==1): 
                            xvw = append(xvw,xgal[ind])
                            yvw = append(yvw,ygal[ind])
                            lvw = append(lvw,lum_ext[ind])
                        else: 
                            xvd = append(xvd,xgal[ind])
                            yvd = append(yvd,ygal[ind])
                            lvd = append(lvd,lum_ext[ind])
                else:
                    print 'Array mismatch: ',efile
            f.close()
    else:
        print pfile,' not found'

lmin = 0. ; lmax= 200.
plot_density(x, y, 1, 500, (1000,1000), range=None, cm=plt.cm.Greys, log=True)
plt.xlabel("x(Mpc $h^{-1})$") ; plt.ylabel("y(Mpc $h^{-1})$")
plt.xlim((lmin,lmax)) ;plt.ylim((lmin,lmax))

for i in range(len(xd2)):
    plt.plot(xd2[i],yd2[i],"o",c=colors[0],markeredgecolor=colors[0],markersize=ld2[i]*2./45.)

for i in range(len(xvw)):
    plt.plot(xvw[i],yvw[i],"^",c=colors[1],markeredgecolor=colors[1],markersize=lvw[i]*2./45.)

for i in range(len(xvw)):
    plt.plot(xvd[i],yvd[i],"x",c=colors[2],markeredgecolor=colors[2],markersize=lvd[i]*2./45.)

# Save figures
plt.savefig(plotfile)
print 'Output: ',plotfile
