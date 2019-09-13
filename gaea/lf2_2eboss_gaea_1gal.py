#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from Cosmology import * 
from dust import *
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gaea/'
snap_g = [0.83,1.0] 
nvol_g = 2#5
gaea = '/cosma5/data/durham/violeta/gaea/NebCat-H17_FIRE-H16_z'

# Gaea cosmology
h0_g = 0.73 ; lambda0_g = 0.75
omega0_g = 0.25 ; omegab_g = 0.045
set_cosmology(omega0=omega0_g,omegab=omegab_g, \
                  lambda0=lambda0_g,h0=h0_g, \
                  universe="Flat",include_radiation=False)

vol1 = (500.**3)/5.
av = 1

llsun = np.log10(3.839) + 33.

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'lfs/lf2_ebossmod'
plotfile = outdir+line+'.pdf'
############################# Obs
obsh0 = 0.677
obs_dir = '/cosma5/data/durham/violeta/lines/desi_hod_o2/lf_obs_data/'
#############################

obsnom = ['DEEP2','VVDSDEEP','VVDSWIDE']
obands = ['R24.2','I24','I22.5']

##########

surveys = ['All','DEEP2','VVDS-DEEP','VVDS-Wide','DESI','eBOSS-SGC','eBOSSmod']
ntypes = len(surveys)
zleg = []
cols = get_distinct(ntypes-2)
cols.insert(0,'grey')

# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
fig = plt.figure(figsize=(6.5,14.))

xtit = "${\\rm log}_{10}(L\\rm{"+lline+"}/h^{-2}erg\, s^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dex}^{-1})$"

xmin = 40.2 ; xmax = 43.7
ymin = -5.9 ; ymax = -1.

# Loop over the redshifts of interest
jj = 410
for iz,zz in enumerate(snap_g):
    jj = jj + 1

    lf = np.zeros(shape=(ntypes,len(lhist)))
    lf_ext = np.zeros(shape=(ntypes,len(lhist)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol_g):
        gfile = gaea+str(zz)+'_'+str(ivol+1)+'.dat'
        if (os.path.isfile(gfile)):
            tomag = band_corrected_distance_modulus(zz) 

            ff = open(gfile, 'r')
            for iline, line in enumerate(ff):
                g_nodust = float(line.split()[12])
                r_nodust = float(line.split()[13])
                i_nodust = float(line.split()[14])
                z_nodust = float(line.split()[15])
                
                g = gaea_calzetti_mag(g_nodust,4770,av)
                r = gaea_calzetti_mag(r_nodust,6231,av)
                i = gaea_calzetti_mag(r_nodust,7625,av)
                z = gaea_calzetti_mag(z_nodust,9134,av)

                rz = r-z ; gr = g-r
                
                lum1 = float(line.split()[22]) # log10(L/Lsolar)
                lum = lum1 + llsun + 2*np.log(h0_g)  # h^2/erg/s
                lum_ext = gaea_calzetti_lum(lum,3727,av)

                # Testing-----------
                if (iline > 10000):
                    break
                #-------------------

                for index,survey in enumerate(surveys):
                    if (survey == 'All'):
                        if (lum_ext>0.):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if (lum>0.):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1

                    elif (survey == 'DEEP2'):
                        fluxcut = 2.7*10.**-17
                        mcut = 24.1
                        mag = r + tomag

                        lcut = emission_line_luminosity(fluxcut,zz)
                        if((mag<mcut) & (lum_ext>lcut)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((mag<mcut) & (lum>lcut)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1
                        
                    elif (survey == 'VVDS-DEEP'):
                        fluxcut = 1.9*10.**-17.
                        mcut = 24.
                        mag = i + tomag

                        lcut = emission_line_luminosity(fluxcut,zz)
                        if((mag<mcut) & (lum_ext>lcut)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((mag<mcut) & (lum>lcut)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1
                        
                    elif (survey == 'VVDS-WIDE'):
                        fluxcut = 3.5*10.**-17.
                        mcut = 22.5
                        mag = i + tomag

                        lcut = emission_line_luminosity(fluxcut,zz)
                        if((mag<mcut) & (lum_ext>lcut)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((mag<mcut) & (lum>lcut)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1
                        
                    elif (survey == 'DESI'):
                        fluxcut = 8.*10.**-17. #erg/s/cm^2

                        lcut = emission_line_luminosity(fluxcut,zz)
                        if((r<23.4) & \
                           (rz>0.3) & (gr>-0.3) & \
                           (rz>0.9*gr+0.12) & \
                           (rz<1.345-0.85*gr) & \
                           (lum_ext>lcut)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((r<23.4) & \
                           (rz>0.3) & (gr>-0.3) & \
                           (rz>0.9*gr+0.12) & \
                           (rz<1.345-0.85*gr) & \
                           (lum>lcut)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1

                    elif (survey == 'eBOSS-SGC'):
                        fluxcut = 10.**-16. #erg/s/cm^2

                        lcut = emission_line_luminosity(fluxcut,zz)                        
                        if((lum_ext>lcut) & \
                           (g>21.825) & (g<22.825) & \
                           (gr>-0.068*rz + 0.457) & \
                           (gr<0.112*rz + 0.773) & \
                           (rz>0.218*gr + 0.571) & \
                           (rz<-0.555*gr + 1.901)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((lum>lcut) & \
                           (g>21.825) & (g<22.825) & \
                           (gr>-0.068*rz + 0.457) & \
                           (gr<0.112*rz + 0.773) & \
                           (rz>0.218*gr + 0.571) & \
                           (rz<-0.555*gr + 1.901)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1

                    elif (survey == 'eBOSSmod'):
                        fluxcut = 10.**-16. #erg/s/cm^2

                        lcut = emission_line_luminosity(fluxcut,zz)                        
                        if((lum_ext>lcut) & \
                           (g>21.825) & (g<22.825) & \
                           (gr>-0.068*rz + 0.457) & \
                           (gr<0.112*rz + 0.773) & \
                           #(rz>0.218*gr + 0.571) & \
                           (rz>0.218*gr + 0.85) & \
                           (rz<-0.555*gr + 1.901)):
                            ibin = np.digitize(lum_ext,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf_ext[index,ibin] = lf_ext[index,ibin] + 1
                        if((lum>lcut) & \
                           (g>21.825) & (g<22.825) & \
                           (gr>-0.068*rz + 0.457) & \
                           (gr<0.112*rz + 0.773) & \
                           #(rz>0.218*gr + 0.571) & \
                           (rz>0.218*gr + 0.85) & \
                           (rz<-0.555*gr + 1.901)):
                            ibin = np.digitize(lum,lbins) - 1
                            if (ibin >=0. and ibin <= len(lhist)-1):
                                lf[index,ibin] = lf[index,ibin] + 1

    lf = lf/dl/volume
    lf_ext = lf_ext/dl/volume
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    if (iz == 0):
        ax1 = fig.add_subplot(jj) ; ax1.set_autoscale_on(False)
        ax1.set_xlim([xmin,xmax]) ; ax1.set_ylim([ymin,ymax]) 
        ax1.set_autoscale_on(False) ;  ax1.minorticks_on()
        ax1.set_xlabel(xtit) ; ax1.set_ylabel(ytit)
        #ax1.tick_params(labelsize=fs-2) 
        ax1.text(42.5, -1.7, zleg[iz])
    else:
        ax = fig.add_subplot(jj,sharex=ax1,sharey=ax1)
        ax.set_autoscale_on(False) ;  ax.minorticks_on()
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.text(42.5, -1.7, zleg[iz])

    # Plot all observations
    ox, oy, el, eh = jc.read_jc_lf(obs_dir+'lf_may16_comparat/',zz,h0=obsh0,\
                                       infile=\
                                       'O2_3728-data-summary-Planck15.txt')
    ind = np.where(oy>-5) 
    oxr = ox[ind] ; oyr = oy[ind]
    arrinds = oxr.argsort()
    oxr = oxr[arrinds]
    oyr = oyr[arrinds]

    if(isinstance(ox, (np.ndarray))):
        if (iz == 0):
            ax1.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                             ecolor='grey',color='grey',mec='grey')
        else:
            ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                            ecolor='grey',color='grey',mec='grey')
    # Plot the observations  O2_3728-*-z*.txt
    for i,isurvey in enumerate(obsnom):        
        ox, oy, el, eh = jc.read_jc_indlf(obs_dir+'lf_may16_comparat/individual_LF/',\
                                              zz,h0=obsh0,\
                                              line='O2_3728',\
                                              survey=isurvey,band=obands[i])
        if(isinstance(ox, (np.ndarray))):
            col = cols[i+1]
            if (iz == 0):
                ax1.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                                 ecolor=col,color=col,mec=col)
            else:
                ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                                ecolor=col,color=col,mec=col)

    # Plot Prabhakar Tiwari's LF
    if (zsnap == 41):
        oxh, oyh, oeh = np.loadtxt(obs_dir+'tiwari/OII_3728_LF_v11.dat',
                                 usecols=(0,1,2),unpack=True,skiprows=1)

        ox = oxh + 2*np.log10(h0)
        oy = oyh + 3*np.log10(obsh0) - 3*np.log10(h0)
        oe = oeh + 3*np.log10(obsh0) - 3*np.log10(h0)

        i = inleg.index('eBOSS-SGC') 

        if (iz == 0):
            ax1.errorbar(ox,oy,yerr=[oe,oe],fmt='o',\
                             ecolor=cols[i],color=cols[i],mec=cols[i])
        else:
            ax.errorbar(ox,oy,yerr=[oe,oe],fmt='o',\
                            ecolor=cols[i],color=cols[i],mec=cols[i])

    # Plot the model predictions
    for index in range(ntypes):
        # Attenuated
        py = 0. ; py = lf_ext[index,:]
        ind = np.where(py > 0)
        x = lhist[ind]
        y = np.log10(py[ind])
        ind = np.where(y < 0.)       
        if (iz == 0):
            if (index<ntypes-1):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                         label=inleg[index])
            else:
                ax1.plot(x[ind],y[ind],color=cols[index-1],linestyle='--')
        else:
            if (index<ntypes-1):
                ax.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                        label=inleg[index])
            else:
                ax.plot(x[ind],y[ind],color=cols[index-1],linestyle='--')

        # Ratios
        if (index == 0):
            my = np.interp(oxr,x,y) 
            diff = abs(my-oyr) ; ratio = 10.**(diff)
            print max(ratio),' diff(min,max)',min(diff),max(diff)
            print oxr,ratio

        # Intrinsic
        if (index == 0):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle=':')
            else:
                ax.plot(x[ind],y[ind],color=cols[index],linestyle=':')

        # Intrinsic eBOSS
        if (index == 5):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle=':')
            else:
                ax.plot(x[ind],y[ind],color=cols[index],linestyle=':')

        # Intrinsic eBOSSmod
        if (index == 6):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index-1],linestyle='-.')
            else:
                ax.plot(x[ind],y[ind],color=cols[index-1],linestyle='-.')

        # Legend
        if (iz == len(snap_list)-1):
            leg = plt.legend(loc=1, handlelength=0, handletextpad=0)
            renderer = fig.canvas.get_renderer()
            shift1 = max([t.get_window_extent(renderer).width for t in leg.get_texts()])
            shift2 = min([t.get_window_extent(renderer).width for t in leg.get_texts()])
            for item in leg.legendHandles:
                item.set_visible(False)
            for color,text in zip(cols,leg.get_texts()):
                text.set_color(color)
                text.set_ha('right') ; text.set_position((shift1-shift2,0))
                leg.draw_frame(False)           

    
# Save figures
fig.subplots_adjust(hspace=0)
fig.savefig(plotfile)
print 'Output: ',plotfile
###########################################
