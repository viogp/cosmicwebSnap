#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from Cosmology import * 
from distinct_colours import get_distinct 
import mpl_style
plt.style.use(mpl_style.style1)

#path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
#model = 'MillGas/gp15newmg/' #model = 'MillGas/gp14/'

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/'
#model = 'MillGas/gp17/' 
model = 'MillGas/gp17.spin/' 

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'
plotfile = outdir+'plots/cuts/'+model+'hod_'+line+'_bot_centrals.pdf'
#############################

snap_list = [44, 42, 40, 37, 34] #MillGas
nvol = 64

bands = ['DEIMOS-R','MegaCam-i-atmos','MegaCam-i-atmos']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17., 10.**-17., 8*10.**-17.]

surveys = ['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']    
ntypes = len(surveys)
cols = get_distinct(ntypes) 

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15
fig = plt.figure(figsize=(8.5,9.))

xtit = "${\\rm log}_{10}(M_{\\rm halo}/M_{\odot}h^{-1})$"
ytit = "$\\langle N_M\\rangle$"

xmin = 10. ; xmax = 15.
ymin = -3. ; ymax = 2.

# Loop over the redshifts of interest
jj = 320
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    nm = np.zeros(shape=(len(surveys),len(lhist)))
    cen= np.zeros(shape=(len(surveys),len(lhist)))
    sat= np.zeros(shape=(len(surveys),len(lhist)))
    cs = np.zeros(shape=(len(surveys),len(lhist)))
    cd = np.zeros(shape=(len(surveys),len(lhist)))
    ss = np.zeros(shape=(len(surveys),len(lhist)))
    sd = np.zeros(shape=(len(surveys),len(lhist)))
    nh = np.zeros(shape=(len(lhist)))
    mhcd = np.zeros(shape=(len(surveys)))
    mhcs = np.zeros(shape=(len(surveys)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        #print gfile
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
            sfr = f['Output001/mstardot'].value +\
                f['Output001/mstardot_burst'].value
            mass = f['Output001/mstars_disk'].value +\
                f['Output001/mstars_bulge'].value
            lssfr = np.zeros(shape=(len(mass))) ; lssfr.fill(-999.)
            ind = np.where((sfr>0.) & (mass>0.))                   
            lssfr[ind] = np.log10(sfr[ind]) - np.log10(mass[ind])
                                                                   
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0, h0=h0, \
                              universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)
            slim = 0.3/tHubble(zz) # Franx+08
            if (firstpass):
                szz = "{:.2f}".format(zz)     

            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                # Number of haloes per bin
                gtype  = f['Output001/type'].value
                mhhalo = f['Output001/mhhalo'].value
                BoT = f['Output001/BoT'].value
                ind = np.where(gtype == 0)
                ll = np.log10(mhhalo[ind])
                H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                nh = nh + H
            
                for index,ib in enumerate(bands):
                    mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                    lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                    icut = mcuts[index] ; fluxcut = fcuts[index]
                    lcut = emission_line_luminosity(fluxcut,zz)

                    # All
                    ind = np.where((mag<icut) & (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[index,:] = nm[index,:] + H

                    # Centrals 
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype==0))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cen[index,:] = cen[index,:] + H

                    # Centrals Spheroids
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype==0) & (BoT>0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cs[index,:] = cs[index,:] + H

                        av = np.mean(mhhalo[ind])
                        if firstpass:
                            mhcs[index] = av
                        else:
                            mhcs[index] = (mhcs[index]+av)/2.

                    # Satellites 
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype>0))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        sat[index,:] = sat[index,:] + H

                    # Satellites Spheroids
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype>0) & (BoT>0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        ss[index,:] = ss[index,:] + H

                    # Centrals Disk
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype==0) & (BoT<=0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cd[index,:] = cd[index,:] + H

                        av = np.mean(mhhalo[ind])
                        if firstpass:
                            mhcd[index] = av
                        else:
                            mhcd[index] = (mhcd[index]+av)/2.

                    # Satellites Disk
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype>0) & (BoT<=0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        sd[index,:] = sd[index,:] + H

                #eBOSS & DESI
                g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag 
                rz = r-z ; gr = g-r
                                  
                fluxcut = 10.**-17.
                lcut = emission_line_luminosity(fluxcut,zz)
                  
                #eBOSS decam180 selection all 
                index = 3
                # All
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nm[index,:] = nm[index,:] + H

                # Centrals 
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype==0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cen[index,:] = cen[index,:] + H

                # Centrals spheroids
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype==0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cs[index,:] = cs[index,:] + H

                    av = np.mean(mhhalo[ind])
                    if firstpass:
                        mhcs[index] = av
                    else:
                        mhcs[index] = (mhcs[index]+av)/2.


                # Satellites
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype>0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sat[index,:] = sat[index,:] + H

                # Satellites spheroids
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype>0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    ss[index,:] = ss[index,:] + H

                # Centrals disks
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype==0) & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cd[index,:] = cd[index,:] + H

                    av = np.mean(mhhalo[ind])
                    if firstpass:
                        mhcd[index] = av
                    else:
                        mhcd[index] = (mhcd[index]+av)/2.

                # Disk Satellites
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & (lum_ext>lcut) & \
                                   (gtype>0)  & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sd[index,:] = sd[index,:] + H

                #DESI
                index = 4
                # All
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nm[index,:] = nm[index,:] + H

                # Centrals 
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype==0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cen[index,:] = cen[index,:] + H

                # Centrals Spheroid
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype==0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cs[index,:] = cs[index,:] + H

                    av = np.mean(mhhalo[ind])
                    if firstpass:
                        mhcs[index] = av
                    else:
                        mhcs[index] = (mhcs[index]+av)/2.

                # Satellites
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype>0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sat[index,:] = sat[index,:] + H


                # Satellites spheroid
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype>0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    ss[index,:] = ss[index,:] + H

                # Centrals disk
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype==0) & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cd[index,:] = cd[index,:] + H

                    av = np.mean(mhhalo[ind])
                    if firstpass:
                        mhcd[index] = av
                    else:
                        mhcd[index] = (mhcd[index]+av)/2.

                # Satellites disk
                ind = np.where((r<23.4) &\
                                   (rz>0.3) & (gr>-0.3) &\
                                   (rz>0.9*gr+0.12) &\
                                   (rz<1.345-0.85*gr) &\
                                   (lum_ext>lcut) & \
                                   (gtype>0) & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sd[index,:] = sd[index,:] + H

                if firstpass:
                    firstpass = False
                f.close()

                
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    ax = fig.add_subplot(jj)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)
    ax.text(xmax-(xmax-xmin)*0.2, ymin+(ymax-ymin)*0.05, 'z='+szz)

    # Plot the model predictions
    for index,survey in enumerate(surveys):
        ig = survey
        indh = np.where(nh > 0) ; nhalos = sum(nh)
        x = lhist[indh]

        py = nm[index,:] ; nall = sum(py)
        yall = py[indh]/nh[indh] 
        y = np.log10(yall) 
        #ax.plot(x,y,color=cols[index],linestyle=':')

        py = cen[index,:] ; ncen = sum(py)
        ycen = py[indh]/nh[indh] 
        y = np.log10(ycen) 
        #ax.plot(x,y,color=cols[index],linestyle='-')

        py = cs[index,:] ; ncs = sum(py)
        ycs = py[indh]/nh[indh]
        y = np.log10(ycs) 
        ax.plot(x,y,color=cols[index],linestyle='-')

        py = cd[index,:]  ; ncd = sum(py)
        ycd = py[indh]/nh[indh]
        y = np.log10(ycd)
        extras = "<Md>=%.1f, <Ms>=%.1f" % (np.log10(mhcd[index]),np.log10(mhcs[index]))
        ax.plot(x,y,color=cols[index],\
                    linestyle='--',\
                    label=ig+extras)
        
        py = sat[index,:] ; nsat = sum(py)
        ysat = py[indh]/nh[indh]
        y = np.log10(ysat) 
        #ax.plot(x,y,color=cols[index],linestyle='--')

        py = ss[index,:] ; nss = sum(py)
        yss = py[indh]/nh[indh]
        y = np.log10(yss) 
        #ax.plot(x,y,color=cols[index],linestyle='--')
        
        py = sd[index,:]  ; nsd = sum(py)
        ysd = py[indh]/nh[indh]
        y = np.log10(ysd) 
        #ax.plot(x,y,color=cols[index],\
        #            marker='x',markeredgecolor=cols[index])

        # Legend
        leg = plt.legend(loc=2,prop={'size':10})
        for color,text in zip(cols,leg.get_texts()):
            text.set_color(color)
            leg.draw_frame(False)
                
        # Output files
        tofile = zip(x,yall,ycen, ysat, ycs,ycd,yss,ysd)
        output = outdir+'hod/'+model+'z'+str(zsnap)+'_'+line+'_'+survey+'_bot_central.txt'   
        with open(output, 'w') as outf:
            outf.write('# Mean HOD: '+model+', snap='+str(zsnap)+', '+line+', '+ig+' \n')
            outf.write('# \n')
            outf.write('# Total number of haloes='+str(nhalos)+' , Number of ELGs='+str(nall)+' \n')
            outf.write('# Percentages of: Central ELGs='+str(100.*ncen/nall)+' , Spheroids ones='+str(100.*ncs/ncen)+' \n')
            outf.write('# Percentages of: Satellite ELGs='+str(100.*nsat/nall)+' , SF ones='+str(100.*nss/nsat)+' \n')
            outf.write('# \n')
            outf.write('# log10(M*/Msunh-1), Total, cen, sat, Central spheroids, Central Disks, Satellite spheroids, Satellite disks \n')
            np.savetxt(outf,tofile,fmt=('%.5f'), delimiter=',')
        outf.closed


# Save figures
fig.tight_layout()
fig.savefig(plotfile)
print 'Output: ',plotfile
print '        ',output
