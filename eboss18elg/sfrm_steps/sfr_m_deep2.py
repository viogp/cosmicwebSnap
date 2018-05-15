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
model = 'gp18'
nvol = 2 #64

#snap = '41' ; survey = 'eBOSS'
snap = '39' ; survey = 'DEEP2'  

plotfile = '/gpfs/data/violeta/lines/cosmicweb/plots/'+model+\
    '/sfr_m_'+survey+'.pdf'

inleg = ['All','F$_{[OII]}>2.7\cdot10^{-17}$',\
             '$r<24.1$','DEEP2  cuts']

ntypes = len(inleg)

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5
# Initialize SFR
smin = 3.
smax = 12.
ds = 0.1
sbins = np.arange(smin,smax,ds)
shist = sbins + ds*0.5

nlevel = 10

gsmf  = np.zeros(shape=(ntypes,len(mbins)))
sfrf = np.zeros(shape=(ntypes,len(sbins)))
smf   = np.zeros(shape=(ntypes,len(sbins),len(mbins)))

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

# Loop over volumes
volume = 0. ; first = True
for ivol in range(nvol):
    gfile = path+model+'/iz'+snap+\
        '/ivol'+str(ivol)+'/galaxies.hdf5'
    if(not os.path.isfile(gfile)):
        print 'STOP: NOT ',gfile ; sys.exit()

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

    mdisk = group['mstars_disk'].value
    mbulge = group['mstars_bulge'].value
    mass1 = mdisk + mbulge

    if (model=='gp14'):
        mass1 = mass1/0.81

    sdisk = group['mstardot'].value # Msolar/h/Gyr
    sbulge = group['mstardot_burst'].value
    sfr1 = sdisk + sbulge
    #if (model=='gp14'):
    #    sfr1 = sfr1/1.26

    f.close()

    for index in range(ntypes):
        sel0 = (mass1>10**mmin) & (sfr1>10**smin)

        if (index == 0):
            sel = sel0

        elif (index == 1):
            if (survey == 'DEEP2'):
                fluxcut = 2.7*10.**-17
            elif (survey == 'VVDS-DEEP'):
                fluxcut = 1.9*10.**-17.
            elif (survey == 'VVDS-WIDE'):
                fluxcut = 3.5*10.**-17.
            elif (survey == 'eBOSS'):
                fluxcut = 10.**-16. #erg/s/cm^2
            elif (survey == 'DESI'):
                fluxcut = 8.*10.**-17. #erg/s/cm^2

            lcut = emission_line_luminosity(fluxcut,zz)
            sel = sel0 & (lum_ext>lcut)
############################
        elif (index == 2):
            if (survey == 'DEEP2'):
                fluxcut = 2.7*10.**-17
                mcut = 24.1
                band = 'DEIMOS-R'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'VVDS-DEEP'):
                fluxcut = 1.9*10.**-17.
                mcut = 24.
                band = 'MegaCam-i-atmos'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'VVDS-WIDE'):
                fluxcut = 3.5*10.**-17.
                mcut = 22.5
                band = 'MegaCam-i-atmos'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'eBOSS'):
                fluxcut = 10.**-16. #erg/s/cm^2

                g = group['mag_DES-g_o_tot_ext'].value + tomag 
                r = group['mag_DES-r_o_tot_ext'].value + tomag 
                z = group['mag_DES-z_o_tot_ext'].value + tomag 
                rz = r-z ; gr = g-r

                sel0 = (g>21.825) & (g<22.825) & \
                    (gr>-0.068*rz + 0.457) & \
                    (gr< 0.112*rz + 0.773) & \
                    (rz> 0.218*gr + 0.571) & \
                    (rz<-0.555*gr + 1.901)

            elif (survey == 'DESI'):
                fluxcut = 8.*10.**-17. #erg/s/cm^2

                g = group['mag_DES-g_o_tot_ext'].value + tomag 
                r = group['mag_DES-r_o_tot_ext'].value + tomag 
                z = group['mag_DES-z_o_tot_ext'].value + tomag 
                rz = r-z ; gr = g-r

                sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                    (rz>0.9*gr+0.12) & (rz<1.345-0.85*gr)


            sel = sel0 & (lum_ext>lcut)
########################
        elif:
            if (survey == 'DEEP2'):
                fluxcut = 2.7*10.**-17
                mcut = 24.1
                band = 'DEIMOS-R'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'VVDS-DEEP'):
                fluxcut = 1.9*10.**-17.
                mcut = 24.
                band = 'MegaCam-i-atmos'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'VVDS-WIDE'):
                fluxcut = 3.5*10.**-17.
                mcut = 22.5
                band = 'MegaCam-i-atmos'

                mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                sel0 = (mag < mcut)

            elif (survey == 'eBOSS'):
                fluxcut = 10.**-16. #erg/s/cm^2

                g = group['mag_DES-g_o_tot_ext'].value + tomag 
                r = group['mag_DES-r_o_tot_ext'].value + tomag 
                z = group['mag_DES-z_o_tot_ext'].value + tomag 
                rz = r-z ; gr = g-r

                sel0 = (g>21.825) & (g<22.825) & \
                    (gr>-0.068*rz + 0.457) & \
                    (gr< 0.112*rz + 0.773) & \
                    (rz> 0.218*gr + 0.571) & \
                    (rz<-0.555*gr + 1.901)

            elif (survey == 'DESI'):
                fluxcut = 8.*10.**-17. #erg/s/cm^2

                g = group['mag_DES-g_o_tot_ext'].value + tomag 
                r = group['mag_DES-r_o_tot_ext'].value + tomag 
                z = group['mag_DES-z_o_tot_ext'].value + tomag 
                rz = r-z ; gr = g-r

                sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                    (rz>0.9*gr+0.12) & (rz<1.345-0.85*gr)
#####################
        ind = np.where(sel)

        mass = np.log10(mass1[ind])
        sfr = np.log10(sfr1[ind])
 
        # GSMF
        H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
        gsmf[index,:] = gsmf[index,:] + H

        # SFR
        H, bins_edges = np.histogram(sfr,bins=np.append(sbins,smax))
        sfrf[index,:] = sfrf[index,:] + H

        # SFR-GSMF
        H, xedges, yedges = np.histogram2d(ssfr,mass,bins=[np.append(sbins,smax),np.append(mbins,mmax)])
        smf[index,:,:] = smf[index,:,:] + H
            
if (volume>0.):
    gsmf[index,:] = gsmf[index,:]/volume/dm  # In Mpc^3/h^3
    sfrf[index,:] = sfrf[index,:]/volume/ds  
    smf[index,:,:] = smf[index,:,:]/volume/dm/ds  
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

## Figure http://matplotlib.org/users/gridspec.html
#fig = plt.figure(figsize=(8.5,9.))
#gs = gridspec.GridSpec(3, 3)
#gs.update(wspace=0., hspace=0.)
#ax = plt.subplot(gs[1:,:-1])
#cols = get_distinct(ntypes + 1) 
#colors = cols #; g = ['grey'] ; colors.extend(g)
#colors[len(colors)-1] = 'grey'
#
## SFRF vs M
#xtit="$log_{10}(\\rm M_{*}/M_{\odot}h^{-1})$"
#ytit="$log_{10}(\\rm SFR/M_{\odot}h^{-1}Gyr^{-1})$"
#xmin=8.5 ; xmax=11.9 ; ymin=3. ; ymax=11.9
#ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
#ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
#
#
## Plot Franx+08 cut
#y = [np.log10(slim),np.log10(slim)] ; x= [xmin,xmax]
#ax.plot(x,y,'k',linewidth=1)
#y = [np.log10(0.3*slim),np.log10(0.3*slim)] ; x= [xmin,xmax]
#ax.plot(x,y,'k',linewidth=2)
#
#lsty = ['solid','dashed','dashed','solid']
#lwdt = [3.,1.5,1.5,1.5]
#for ii in range(ntypes):
#    matplotlib.rcParams['contour.negative_linestyle'] = lsty[ii]
#    zz = np.zeros(shape=(len(shist),len(mhist))) 
#    pz = smf[ii,:,:] ; ind = np.where(pz>0.)
#    zz = np.log10(pz) 
#
#    # Heat map
#    #x,y = np.meshgrid(np.append(mbins,mmax),np.append(sbins,smax))
#    #ax.pcolormesh(x, y, zz, cmap=plt.get_cmap(ptmap))
#    # Contours
#    xx,yy = np.meshgrid(mbins,sbins)
#    #zz = 100.*zz ; al = [5.,16.,68.,95.] 
#    al = [-3.5,-2.5,-1.5] 
#    # Smooth
#    #zz = ndimage.gaussian_filter(zz, sigma=0.5)
#    cs = ax.contour(xx, yy, zz, levels=al, \
#                        colors=cols[ii],linewidths=lwdt[ii])
#    cs.levels = [nf(val) for val in cs.levels]
#    ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
#                  fontsize=10,fmt='%r')#fmt='%r %%')
#
## GSMF ###################################
#axm = plt.subplot(gs[0, :-1],sharex=ax)
#ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
#axm.set_autoscale_on(False) ;  axm.minorticks_on()
#axm.set_ylim(-5.5,-1.) 
#plt.setp(axm.get_xticklabels(), visible=False)
#
## Plot observations from Baldry+2012
#dobs = '/cosma/home/violeta/Galform2/galform-2.6.0/Obs_Data2/'
#file = dobs+'mf/mf_baldry_2012.txt'
#oh = 0.7 
#lm,p3,dp3 = np.loadtxt(file,usecols=[0,1,2],unpack=True)            
#xobs = lm + np.log10(oh)
#yobs = xobs*0. - 999.
#indx = np.where( p3 > 0)
#yobs[indx] = np.log10(p3[indx]) -3. -3.*np.log10(oh)
#lerr = yobs*0. - 999. 
#indx = np.where( (p3-dp3) > 0)
#lerr[indx]  = np.log10(p3[indx] - dp3[indx])-3. -3.*np.log10(oh)
#herr = yobs*0. + 999.
#indx = np.where( (p3+dp3) > 0)
#herr[indx]  = np.log10(p3[indx] + dp3[indx])-3. -3.*np.log10(oh)
#plt.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], fmt='o', ecolor='grey',color='grey', mec='grey',label="Baldry+2012, z<0.06")
#
## Models
##lsty = ['solid','dashed','dotted','solid']
#lsty = ['-','--',':','-']
#lwdt = [3.,1.5,1.5,1.5]
#for ii in range(ntypes):
#    py = gsmf[ii,:] ; ind = np.where(py>0.)
#    x = mhist[ind] ; y = np.log10(py[ind])
#    ind = np.where(y < 0.)
#    axm.plot(x[ind],y[ind],color=cols[ii],\
#                 linestyle=lsty[ii],linewidth=lwdt[ii], \
#                 label=inleg[ii])
#
#
## SFRF
#axs = plt.subplot(gs[1:, 2],sharey=ax)
#xtit="$log_{10}(\Phi)$" ; axs.set_xlabel(xtit)
#axs.set_autoscale_on(False) ;  axs.minorticks_on()
#axs.set_xlim(-4.4,0.0)
#start, end = axs.get_xlim()
#axs.xaxis.set_ticks(np.arange(-4., end, 1.))
#plt.setp(axs.get_yticklabels(), visible=False)
#
#for ii in range(ntypes):
#    px = sfrf[ii,:] ; ind = np.where(px>0.)
#    y = shist[ind] ; x = np.log10(px[ind])
#    ind = np.where(x < 0.)
#    axs.plot(x[ind],y[ind],color=cols[ii],label=inleg[ii],\
#                 linestyle=lsty[ii],linewidth=lwdt[ii])
#
## Legend
#leg = axs.legend(bbox_to_anchor=(1., 1.4),fontsize='small')
#for color,text in zip(cols,leg.get_texts()):
#    text.set_color(color)
#    leg.draw_frame(False)
#
#
## Save figures
#plotfile = plotdir + 'ssfr_m.pdf'
#fig.savefig(plotfile)
#print 'Output: ',plotfile
