# Pretty plot with spatial distribution
# compare to ~/lines/desi_hod_o2/pretty_image/millgas_o2.py

import os.path, sys
from glob import glob
from numpy import *
import h5py
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distinct_colours import get_distinct

model = 'gp18'

snap = '41' ; survey = 'eBOSS'
#snap = '39' ; survey = 'DEEP2'

endings = ['_m8.5_sn'+snap+'.dat','_s3.0_sn'+snap+'.dat']
root = ['mcut_','scut_']

##############################
plotdir =  '/gpfs/data/violeta/lines/cosmicweb/selections/'
path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
line = 'OII3727' ; lline = '[OII]'

files = [] 
for i, ending in enumerate(endings):
    censf = survey+'_centrals'+ending
    satf  = survey+'_satellites'+ending
    allf  = root[i]+survey+ending
    files = files + [censf,satf,allf]
    #files = files + (glob(plotdir+model+'/ascii_files/*'+survey+'*'+ending))

ntypes = len(files)
cols = get_distinct(ntypes)
symb = ['o','^','x','o','^','+']

zup = ; zlow = 
volume = 500.*500.*(zup-zlow)

# Loop over files
for fil in files:
    infil = plotdir+model+'/ascii_files/'+fil
    if (not os.path.isfile(infil)):
        print('STOP: {} not found'.format(infil)) ; sys.exit()

    xgal,ygal,zgal,lum_ext = np.loadtxt(infil, usecols=(0,1,2,10),\
                                            unpack=True)





# DEEP2
print 'DEEP2 density=',len(xd2)/volume
for i in range(len(xd2)):
    x = xd2[i]-shift
    y = yd2[i]-shift
    plt.plot(x,y,"o",c=cols[0],markeredgecolor=cols[0],\
                 markersize=ld2[i]*4.,alpha=0.8)

# VVDSWIDE
print 'VVDSWIDE density=',len(xvw)/volume
for i in range(len(xvw)):
    x = xvw[i]-shift
    y = yvw[i]-shift
    #plt.plot(x,y,"^",c=cols[1],markeredgecolor=cols[1],\
    #             markersize=lvw[i]*5.,alpha=0.8)

# Haloes
print 'Haloes density=',len(xh)/volume, min(lh),max(lh)
for i in range(len(xh)):
    x = xh[i]-shift
    y = yh[i]-shift
    plt.plot(x,y,"o",c='none',markeredgecolor='r',\
                 markersize=lh[i]/1.5,linewidth=3)


# Save figures
plotfile = plotdir+model+'deep2z0.9.pdf'
plt.savefig(plotfile)
print 'Output: ',plotfile
