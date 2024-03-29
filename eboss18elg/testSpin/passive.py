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

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/MillGas/'
nvol = 3 #64

plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/spin.acce'
models = ['gp17','gp17.spin','gp17.spin.acce0.1'] ; inleg = models

#plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/spin.'
#models = ['gp17','gp17.spin']
#inleg = ['GP18','GP18.BH'] 

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5

ntot  = np.zeros(shape=(len(models),len(mbins)))
npas1 = np.zeros(shape=(len(models),len(mbins)))
npas2 = np.zeros(shape=(len(models),len(mbins)))

nsat  = np.zeros(shape=(len(models),len(mbins)))
nsp1 = np.zeros(shape=(len(models),len(mbins)))
nsp2 = np.zeros(shape=(len(models),len(mbins)))

ndis  = np.zeros(shape=(len(models),len(mbins)))
nsph  = np.zeros(shape=(len(models),len(mbins)))

nbhm  = np.zeros(shape=(len(models),len(mbins)))
nbhd  = np.zeros(shape=(len(models),len(mbins)))
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
        gfile = path+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'
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
            mbh = group['M_SMBH'].value
            f.close()

            efile = path+model+'/iz61/ivol'+str(ivol)+'/elgs.hdf5'
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
            npas1[index,:] = npas1[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<slim))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            npas2[index,:] = npas2[index,:] +H

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
            nsp1[index,:] = nsp1[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<slim) & (sat>0))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nsp2[index,:] = nsp2[index,:] +H

            # Disks and bulges
            ind = np.where((mass1>0.) & (sfr1>0.))
            ssfr = np.zeros(shape=(len(sfr1)))
            ssfr[ind] = sfr1[ind]/mass1[ind]

            ind = np.where((mass1>0.) & (ssfr<0.3*slim) & (bot<0.5))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            ndis[index,:] = ndis[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<0.3*slim) & (bot>=0.5))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nsph[index,:] = nsph[index,:] +H

            # Cut in BH mass
            ind = np.where((mass1>0.) & (sfr1>0.))
            ssfr = np.zeros(shape=(len(sfr1)))
            ssfr[ind] = sfr1[ind]/mass1[ind]

            val = 1000.

            ind = np.where((mass1>0.) & (ssfr<0.3*slim) & (mbh>val))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nbhm[index,:] = nbhm[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<0.3*slim) &\
                               (mbh>val) & (bot<0.5))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            nbhd[index,:] = nbhd[index,:] +H

        else:
            print 'NOT found:',gfile
            
    if (volume>0.):
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
        for i in range(len(mbins)):
            if ntot[index,i]>0.:
                npas1[index,i] = npas1[index,i]/ntot[index,i]
                npas2[index,i] = npas2[index,i]/ntot[index,i]

                nsp1[index,i] = nsp1[index,i]/ntot[index,i]
                nsp2[index,i] = nsp2[index,i]/ntot[index,i]

                ndis[index,i] = ndis[index,i]/ntot[index,i]
                nsph[index,i] = nsph[index,i]/ntot[index,i]

                nbhm[index,i] = nbhm[index,i]/ntot[index,i]
                nbhd[index,i] = nbhd[index,i]/ntot[index,i]
            else:
                npas1[index,i] = -999.
                npas2[index,i] = -999.

                nsp1[index,i] = -999.
                nsp2[index,i] = -999.

                ndis[index,i] = -999.
                nsph[index,i] = -999.

                nbhm[index,i] = -999.
                nbhd[index,i] = -999.


# Figure http://matplotlib.org/users/gridspec.html
fig = plt.figure(figsize=(8.5,9.))
ax = plt.subplot()
cols = get_distinct(len(models)) 
colors = cols ; g = ['grey'] ; colors.extend(g) ; colors.extend(g)
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

# Models
for ii in range(len(models)):
    py = npas1[ii,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    ax.plot(x[ind],y[ind],color=colors[ii],label=inleg[ii])

    #py = npas2[ii,:] ; ind = np.where(py>0.)
    #y = py[ind]
    #ax.plot(x[ind],y[ind],color=colors[ii],\
    #            linestyle=':')

    #py = nsp1[ii,:] ; ind = np.where(py>0.)
    #x = mhist[ind] ; y = py[ind]
    #print ii,np.shape(x),np.shape(y)
    #ax.plot(x[ind],y[ind],color=colors[ii],linestyle='--')

    #py = npas2[ii,:] ; ind = np.where(py>0.)
    #y = py[ind]
    #ax.plot(x[ind],y[ind],color=colors[ii],\
    #            linestyle='.-')

    py = ndis[ii,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    print ii,np.shape(x),np.shape(y)
    ax.plot(x[ind],y[ind],color=colors[ii],linestyle='--')

    py = nsph[ii,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    print ii,np.shape(x),np.shape(y)
    ax.plot(x[ind],y[ind],color=colors[ii],linestyle=':')

    #py = nbhm[ii,:] ; ind = np.where(py>0.)
    #x = mhist[ind] ; y = py[ind]
    #print ii,np.shape(x),np.shape(y)
    #ax.plot(x[ind],y[ind],color=colors[ii],linestyle='-.')
    #
    #py = nbhd[ii,:] ; ind = np.where(py>0.)
    #x = mhist[ind] ; y = py[ind]
    #print ii,np.shape(x),np.shape(y)
    #ax.plot(x[ind],y[ind],color=colors[ii],linestyle=':')

## Legend
leg = ax.legend(loc=2,fontsize='small')
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
plotfile = plotdir + 'passive.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile

