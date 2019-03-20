#! /usr/bin/env python

import numpy as np
import h5py, os.path, sys
from Cosmology import *
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from distinct_colours import get_distinct,pault_cmap
import matplotlib.gridspec as gridspec
from scipy import ndimage
from stats import *
import mpl_style
plt.style.use(mpl_style.style1) ; ptmap=pault_cmap(1)

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
nvol = 64
sn = '61'

gaea = '/cosma5/data/durham/violeta/gaea/NebCat-H17_FIRE-H16_z0.0_5.dat'

plotdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/modelplots/gaea.'
models = ['gp19','gp19.font','gp19.starvation','gp17']  
inleg = ['This work','10% stripping ','Starvation','GP18','GAEA']

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5

ntot  = np.zeros(shape=(len(inleg),len(mbins)))
npas = np.zeros(shape=(len(inleg),len(mbins)))
nsat = np.zeros(shape=(len(inleg),len(mbins)))
nsp = np.zeros(shape=(len(inleg),len(mbins)))

###########################################
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


# Loop over redshifts
for index,model in enumerate(models):
    volume = 0. ; first = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+sn+\
            '/ivol'+str(ivol)+'/galaxies.hdf5'
        #print gfile
        if(os.path.isfile(gfile)):
            # Read the relevant information from galaxies.hdf5
            f = h5py.File(gfile,'r') #; print gfile
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value


            group = f['Output001']
            zz   = group['redshift'].value
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0,h0=h0, \
                              universe="Flat",include_radiation=False)
            slim = 1./tHubble(zz) 
            if first:
                print 'zz, slim, 0.3*slim=',zz, slim, 0.3*slim
                first = False

            mdisk = group['mstars_disk'].value
            mbulge = group['mstars_bulge'].value
            mass1 = mdisk + mbulge
            sdisk = group['mstardot'].value # Msolar/h/Gyr
            sbulge = group['mstardot_burst'].value
            sfr1 = sdisk + sbulge            

            sat = group['type'].value
            f.close()

            efile = path+model+'/iz'+sn+\
                '/ivol'+str(ivol)+'/elgs.hdf5'
            if(os.path.isfile(efile)):
                ff = h5py.File(efile,'r') 
                bot= ff['Output001/BoT'].value
                ff.close()
            else:
                print 'STOP: No ',efile ; sys.exit()


            # All
            ind = np.where((mass1>0.) & (sfr1>0.))
            ssfr = np.zeros(shape=(len(sfr1)))
            ssfr[ind] = sfr1[ind]/mass1[ind]
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            ntot[index,:] = ntot[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<0.3*slim))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            npas[index,:] = npas[index,:] +H

            # Satellites
            ind = np.where((mass1>0.) & (sfr1>0.) & (sat>0))
            ssfr = np.zeros(shape=(len(sfr1)))
            ssfr[ind] = sfr1[ind]/mass1[ind]
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nsat[index,:] = nsat[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<0.3*slim) & (sat>0))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nsp[index,:] = nsp[index,:] +H

        else:
            print 'NOT found:',gfile
            
    if (volume>0.):
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
        for i in range(len(mbins)):
            if ntot[index,i]>0.:
                npas[index,i] = npas[index,i]/ntot[index,i]
                nsat[index,i] = nsat[index,i]/ntot[index,i]
                nsp[index,i] = nsp[index,i]/ntot[index,i]
            else:
                npas[index,i] = -999.
                nsat[index,i] = -999.
                nsp[index,i] = -999.

############
# Read GAEA
if(not os.path.isfile(gaea)):
    print ('STOP: GAEA file not found ',gaea) ; sys.exit()

zz = 0.05
# Cosmology
h0_g = 0.73 ; lambda0_g = 0.75
omega0_g = 0.25 ; omegab_g = 0.045
set_cosmology(omega0=omega0_g,omegab=omegab_g, \
                  lambda0=lambda0_g,h0=h0_g, \
                  universe="Flat",include_radiation=False)
slim = 1./tHubble(zz) # SFR cut

# Initialize arrays
lmass, ssfr, sat = [np.array([]) for i in range(3)]

# Read file
ff = open(gaea, 'r') 
for iline, line in enumerate(ff):
    sat1 = float(line.split()[4])
    lmass1 =float(line.split()[6]) + np.log10(h0_g)
    ssfr1 = float(line.split()[7])*(10.**9.)/(10.**float(line.split()[6]))

    if (lmass1 > 0.):
        lmass = np.append(lmass, lmass1)
        ssfr = np.append(ssfr, ssfr1)
        sat =  np.append(sat, sat1)

    ##Testing-----------
    #if (iline > 10000):
    #    break
    ##------------------

index = len(inleg)-1
# All
H, bins_edges = np.histogram(lmass,bins=np.append(mbins,mmax))
ntot[index,:] = ntot[index,:] +H

ind = np.where(ssfr<0.3*slim) ; mass = lmass[ind]
H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
npas[index,:] = npas[index,:] +H

# Satellites
ind = np.where(sat>0) ; mass = lmass[ind]
H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
nsat[index,:] = nsat[index,:] +H

ind = np.where((ssfr<0.3*slim) & (sat>0)) ; mass = lmass[ind]
H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
nsp[index,:] = nsp[index,:] +H

# Normalize
volume = pow(500.,3.)/5.            
for i in range(len(mbins)):
    if ntot[index,i]>0.:
            npas[index,i] = npas[index,i]/ntot[index,i]
            nsp[index,i] = nsp[index,i]/ntot[index,i]
            nsat[index,i] = nsat[index,i]/ntot[index,i]
    else:
        npas[index,i] = -999.
        nsp[index,i] = -999.
        nsat[index,i] = -999.

############

# Figure http://matplotlib.org/users/gridspec.html
fig = plt.figure(figsize=(8.5,9.))
ax = plt.subplot()
cols = get_distinct(len(inleg)+1) 
colors = cols ; g = ['grey'] #; colors.extend(g) ; colors.extend(g)
colors[len(colors)-1] = 'grey' ; colors.extend(g)

xtit = "$log(\\rm{M_{*}/M_{\odot} h^{-1}})$"
ytit="Passive fraction"
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ax.set_autoscale_on(False) ;  ax.minorticks_on()
ax.set_xlim(mmin,12) ; ax.set_ylim(0.,1.) 

# Observations
col = colors[len(colors)-1]

fobs = 'Obs_Data/gilbank10.txt'
#log(M_*)   passive_fraction error
lm, p, erro = np.loadtxt(fobs, unpack=True)
xo = lm + np.log10(0.7)
ax.errorbar(xo,p, yerr=erro, color=col, ecolor=col,\
                label ='Gilbank+10', fmt = 'o')

fobs = 'Obs_Data/bauer13.txt'
#log(M_*)   passive_fraction error
lm, p = np.loadtxt(fobs, unpack=True)
xo = lm + np.log10(0.7)
erro = lm*0.
ax.errorbar(xo,p, yerr=erro, color=col, ecolor=col,\
                label ='Bauer+10', fmt = '^')

lsty = ['-','--',':','-','-'] 
lwdt = [3.,1.5,1.5,1.5,2.] 

# Models
for im in range(len(inleg)):
    py = npas[im,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    ax.plot(x,y,color=colors[im],label=inleg[im],\
                linestyle=lsty[im],linewidth=lwdt[im])

    py = nsp[im,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    ax.plot(x,y,color=colors[im],\
                linestyle='-.',linewidth=lwdt[im])

## Legend
leg = ax.legend(loc=2,fontsize='small')
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
plotfile = plotdir + 'passive_sn'+sn+'.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile
