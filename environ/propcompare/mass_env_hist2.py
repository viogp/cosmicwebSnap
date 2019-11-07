import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

propname = 'lmass'
ytit = ''
ymin = 8.5 ; ymax = 15.

##########################################

model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
##########################################

# Separation of environments in the plot
sep = 0.85 

# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15

xlabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

cols = ['darkred','dodgerblue']
hatching = [' ','/','o',' ','//','O'] 

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
##########################################
# Output fraction summary
envsumfile = path+'env_files/'+model+'env_fractions.txt'
sumfile = open(envsumfile,'w')
sumfile.write('File : Fraction in Voids,Sheets,Filaments,Knots \n')
sumfile.close()
sumfile = open(envsumfile,'a')

# Loop over the different files
for cw in ['Vweb','Pweb']:
    epath = path+'env_files/'+model+cw+'/'

    for iis,survey in enumerate(surveys1):
        inleg = ['Mass cut, All','Mass cut, '+survey,'Mass cut, '+surveys2[iis],
                 'SFR cut, All','SFR cut, '+survey,'SFR cut, '+surveys2[iis]]
        numinleg = len(inleg)
        lbar = sep/numinleg

        for iz in ['39','41']:
            for nd in ['-2.0','-3.0','-4.2']:
                # Initialize the parameters for the figures
                fig = plt.figure(figsize=(8.5,9.))
                jj = 111 ; ax = fig.add_subplot(jj)
                plt.xticks(xlabels,elabels)
                ax.tick_params(axis='x',which='minor',bottom=False)
                ax.set_ylabel(ytit,fontsize=fs) ; ax.set_ylim(ymin,ymax)

                if(iz == '41'):
                    ztext = 'z = 0.83; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
                elif(iz == '39'):
                    ztext = 'z = 0.99; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
                ax.text(1.7, 0.95, ztext)

                ii = -1
                for ic, cut in enumerate(['m','sfr']):
                    end = cut+'cut_All_nd'+nd+'_sn'+iz+'.dat'
                    allfile = epath+end
                    allprop = proppath+end

                    end = cut+'cut_'+survey+'_nd'+nd+'_sn'+iz+'.dat'
                    elgfile1 = epath+end
                    elg1prop = proppath+end

                    end = cut+'cut_'+surveys2[iis]+'_nd'+nd+'_sn'+iz+'.dat'
                    elgfile2 = epath+end
                    elg2prop = proppath+end

                    files = [allfile,elgfile1,elgfile2]
                    fprop = [allprop,elg1prop,elg2prop]

                    for iif,efile in enumerate(files):
                        ii += 1
                        pfile = fprop[iif] 

                        # Check if files exist and has more than one line
                        if (not os.path.isfile(efile) or not os.path.isfile(pfile)):
                            continue
                        wcl_line = subprocess.check_output(["wc", "-l",efile])
                        wcl = int(wcl_line.split()[0])
                        pcl_line = subprocess.check_output(["wc", "-l",pfile])
                        pcl = int(pcl_line.split()[0])
                        if (wcl <= 1 or pcl <= 1):
                            continue

                        # Read the file
                        xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
                        env = fenv.astype(int)
                        
                        # Read property file
                        # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
                        # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
                        # lum 9,lum_ext 10 (10^40 h^-2 erg/s),
                        # type 11 (0= Centrals; 1,2= Satellites) 
                        px, py, pz, lmass = np.loadtxt(pfile, usecols=(0,1,2,7), unpack=True)

                        # Chech that the coordinates have the same size
                        if ((len(xx) != len(px)) or 
                            (len(yy) != len(py)) or
                            (len(zz) != len(pz))):
                            print('STOP! Different lengths coordinates: {}\n {}\n'.
                                  format(efile,pfile)) ; sys.exit()
                        # Chech that the coordinates are ordered in the same way
                        if ((not np.allclose(xx,px,atol=1e-08,equal_nan=True)) or 
                            (not np.allclose(yy,py,atol=1e-08,equal_nan=True)) or
                            (not np.allclose(zz,pz,atol=1e-08,equal_nan=True))):
                            print('STOP! Files with different coordinates: {}\n {}\n'.
                                  format(efile,pfile)) ; sys.exit()

                        # Loop over type of environment
                        for ienv in np.unique(env):
                            ind = np.where(env == ienv)
                            if (np.shape(ind)[1] <= 1):
                                continue
                            prop = lmass[ind]
                            # Plot
                            xenv = ienv + 0.5*(1.-sep) + 0.5*lbar + ii*lbar
                            ax.boxplot(prop,positions=[xenv])
###here
                            #ax.bar(xenv, frac, lbar, \
                            #       color=cols[ic], label=inleg[ii])
                            #ax.bar(xenv, frac, lbar, \
                            #       color=cols[ic], label=inleg[ii],\
                            #       hatch=hatching[ii],fill=False,edgecolor=cols[ic])

                newnd = nd.replace('.','p')
                plotfile = path+'plots/'+model+'environ/props/'+cw+\
                           '/'+propname+survey+'_nd'+newnd+'_sn'+iz+'_env2.pdf'
                print(plotfile)
                sys.exit()
                # Legend
                leg = plt.legend(loc=2)
                #for color,text in zip(zip(cols,cols),leg.get_texts()):
                #  text.set_color(cols)
                leg.draw_frame(False)

                # Save figure
                fig.savefig(plotfile)
                print 'Plot: ',plotfile
                plt.close()

sumfile.close()
print 'Enviroment fractions: ',envsumfile
