#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from Cosmology import * 
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp18/'
#gaea_calzetti_lum(llum,leff,av)
#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/cosmicweb/plots/'+model+'lfs/compared2_'
plotfile = outdir+line+'.pdf'
############################# Obs
obsh0 = 0.677
obs_dir = '/gpfs/data/violeta/lines/desi_hod_o2/lf_obs_data/lf_may16_comparat/'
#############################

snap_list = [41,39] #MillGas
nvol = 64

gaea = '/gpfs/data/violeta/gaea/'
gaea_root = gaea+'NebCat-H17_FIRE-H16_z' 
zlist = [0.83,1.0] #Gaea
# Cosmology
h0_g=0.73 ; omega0_g=0.25 ; omegab_g=0.045 ; lambda0_g=0.75
#(log10L+lsun2lgal) erg/s h-2, Lsun taken from the IAU
lsun2lgal = np.log10(3.8275) + 33. + 2*np.log10(h0_g) 

#####
obsnom = ['DEEP2','VVDSDEEP','VVDSWIDE']
obands = ['R24.2','I24','I22.5']

bands = ['DEIMOS-R','MegaCam-i-atmos','MegaCam-i-atmos','eBOSS-SGC','DESI']
mcuts = [24.1, 24, 22.5]
fcuts = [2.7*10.**-17., 1.9*10.**-17., 3.5*10.**-17.,10.**-16.,8.*10.**-17.] #erg/s/cm2

inleg = ['All','DEEP2','VVDS-DEEP','VVDS-Wide','eBOSS-SGC','DESI']
##########

ntypes = len(inleg)
zleg = []
cols = get_distinct(ntypes-1)
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
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    lf = np.zeros(shape=(ntypes,len(lhist)))
    lf_ext = np.zeros(shape=(ntypes,len(lhist)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            zz   = f['Output001/redshift'].value
            group = f['Parameters']
            vol1 = group['volume'].value ; volume = volume + vol1
            h0 = group['h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 = f['Parameters/lambda0'].value
            f.close()

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
                firstpass = False
            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                              h0=h0, universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz) 

            efile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r') #; print (efile)
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                for index in range(ntypes):
                    if index==0:
                        ind  = np.where(lum_ext>0.)
                        indi = np.where(lum>0.)
                    elif ((inleg[index] == 'eBOSS-SGC') or (inleg[index] == 'DESI')):
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)

                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                        z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        if (inleg[index] == 'eBOSS-SGC'): 
                            ind = np.where((lum_ext>lcut) & \
                                               (g>21.825) & (g<22.825) & \
                                               (gr>-0.068*rz + 0.457) & \
                                               (gr<0.112*rz + 0.773) & \
                                               (rz>0.218*gr + 0.571) & \
                                               (rz<-0.555*gr + 1.901))
                            indi = np.where((lum>lcut) & \
                                               (g>21.825) & (g<22.825) & \
                                               (gr>-0.068*rz + 0.457) & \
                                               (gr<0.112*rz + 0.773) & \
                                               (rz>0.218*gr + 0.571) & \
                                               (rz<-0.555*gr + 1.901))

                        elif (inleg[index] == 'DESI'): 
                            ind  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum_ext>lcut))

                            indi  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum>lcut))

                    else:
                        ib = bands[index-1]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value\
                            + tomag
                        icut = mcuts[index-1]
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)

                        ind  = np.where((mag<icut) & (lum_ext>lcut))
                        indi = np.where((mag<icut) & (lum>lcut))

                        
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H

                    if (np.shape(indi)[1] > 0.):
                        ll = np.log10(lum[indi]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                f.close()

    if (volume>0.):
        lf = lf/dl/volume
        lf_ext = lf_ext/dl/volume
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
    else:
        print('STOP: Nothing read, last file= ',gfile) ; sys.exit()

    #GAEA
    lf_gaea = np.zeros(shape=(ntypes,len(lhist)))
    volume_g = 0.
    for ii in range(5):
        zz = zlist[iz]
        ivol = str(ii+1)
        volume_g = volume_g + (500.**3.)/5.
        gfile = gaea_root+str(zz)+'_'+ivol+'.dat'
        if(not os.path.isfile(gfile)):
            print ('STOP: GAEA file not found ',gfile) ; sys.exit()

        # Read file
        ff = open(gfile, 'r') 
        for iline, ffline in enumerate(ff):
            gmag = float(ffline.split()[12]) #Intrinsic
            rmag = float(ffline.split()[13])
            imag = float(ffline.split()[14])
            zmag = float(ffline.split()[15])
            rz = rmag-zmag ; gr = gmag-rmag

            l02 = float(ffline.split()[22]) + lsun2lgal
            iflux = logL2flux(l02,zz) #erg/s/cm^2

            for index in range(ntypes):
                if index==0:
                    if(l02>-9.):
                        jb = np.digitize(l02, np.append(lbins,lmax))
                        if (jb>0 and jb<len(lbins)+1):
                            jbin = jb-1
                            lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 
                else:
                    fluxcut = fcuts[index-1]

                    if (inleg[index] == 'eBOSS-SGC'): 
                        if((iflux>fluxcut) & \
                               (gmag>21.825) & (gmag<22.825) & \
                               (gr>-0.068*rz + 0.457) & \
                               (gr<0.112*rz + 0.773) & \
                               (rz>0.218*gr + 0.571) & \
                               (rz<-0.555*gr + 1.901)):
                            jb = np.digitize(l02, np.append(lbins,lmax))
                            if (jb>0 and jb<len(lbins)+1):
                                jbin = jb-1
                                lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 

                    elif (inleg[index] == 'DESI'): 
                        if((iflux>fluxcut) & \
                               (rz>0.3) & (gr>-0.3) & \
                               (gr<1.1*rz-0.13) & \
                               (gr<-1.18*rz+1.6)):
                            jb = np.digitize(l02, np.append(lbins,lmax))
                            if (jb>0 and jb<len(lbins)+1):
                                jbin = jb-1
                                lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 
                    elif (inleg[index] == 'DEEP2'):
                        if((iflux>fluxcut) & \
                               (rmag<24.1)):
                            jb = np.digitize(l02, np.append(lbins,lmax))
                            if (jb>0 and jb<len(lbins)+1):
                                jbin = jb-1
                                lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 
                    elif (inleg[index] == 'VVDS-DEEP'):
                        if((iflux>fluxcut) & \
                               (imag<24.)):
                            jb = np.digitize(l02, np.append(lbins,lmax))
                            if (jb>0 and jb<len(lbins)+1):
                                jbin = jb-1
                                lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 
                    elif (inleg[index] == 'VVDS-Wide'):
                        if((iflux>fluxcut) & \
                               (imag<22.5)):
                            jb = np.digitize(l02, np.append(lbins,lmax))
                            if (jb>0 and jb<len(lbins)+1):
                                jbin = jb-1
                                lf_gaea[index,jbin] = lf_gaea[index,jbin]+1 
                        
        
            ##Testing-----------
            #if (iline > 10000):
            #    break
            #------------------        
        ff.close()
    if (volume_g >0.):
        lf_gaea = lf_gaea/dl/volume_g
    else:
        print('STOP: Nothing read, last file= ',gfile) ; sys.exit()

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
    ox, oy, el, eh = jc.read_jc_lf(obs_dir,zz,h0=obsh0,\
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
        ox, oy, el, eh = jc.read_jc_indlf(obs_dir+'individual_LF/',\
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

    # Plot the model predictions
    for index in range(ntypes):
        # Attenuated
        py = 0. ; py = lf_ext[index,:]
        ind = np.where(py > 0)
        x = lhist[ind]
        y = np.log10(py[ind])
        ind = np.where(y < 0.)       
        if (iz == 0):
            ax1.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                         label=inleg[index])
        else:
            ax.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                        label=inleg[index])

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

        # Intrinsic GAEA
        py = 0. ; py = lf_gaea[index,:]
        ind = np.where(py > 0)
        x = lhist[ind]
        y = np.log10(py[ind])
        ind = np.where(y < 0.)
        if (iz == 0):
            ax1.plot(x[ind],y[ind],color=cols[index],linestyle='--')
        else:
            ax.plot(x[ind],y[ind],color=cols[index],linestyle='--')

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


#On splitting the legend
#line1, = plt.plot([1,2,3], label="Line 1", linestyle='--')
#line2, = plt.plot([3,2,1], label="Line 2", linewidth=4)
#
## Create a legend for the first line.
#first_legend = plt.legend(handles=[line1], loc=1)
#
## Add the legend manually to the current Axes.
#ax = plt.gca().add_artist(first_legend)
#
## Create another legend for the second line.
#plt.legend(handles=[line2], loc=4)
    
# Save figures
fig.subplots_adjust(hspace=0)
fig.savefig(plotfile)
print 'Output: ',plotfile
