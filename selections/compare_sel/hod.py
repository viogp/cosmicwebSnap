import sys, os.path
import subprocess
import numpy as np

model = 'gp19/'

sn_list = ['41','39']
surveys = ['All','DEEP2','VVDS-DEEP','eBOSS-SGC','DESI']
nds = ['-2.0','-3.0','-4.2']
cuts = ['test'] #cuts = ['m','sfr']

verbose = False

# Paths
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
ndpath = inpath+'selections/'
hodpath = inpath+'hod/'

# Read the information from the different surveys
for sn in sn_list: 
    # Read the number of haloes
#    hmffile = hodpath+model+'hmf_sn'+sn+'.txt'
    hmffile = hodpath+model+'hmf_sn'+sn+'_2vols.txt'
    mhist,mlow,mhigh,nhs = np.loadtxt(hmffile,unpack='True')

    for nd in nds:
        for survey in surveys:
            for cut in cuts:
                infile = ndpath+model+'ascii_files/'+cut+\
                         'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'
                # Check the existance of the files
                if (not os.path.isfile(infile)):
                    if verbose:
                        print('WARNING: {} not found'.format(infile))
                    continue
                # Jump files with only a one line header
                wcl_line = subprocess.check_output(["wc", "-l",infile])
                wcl = int(wcl_line.split()[0])
                if (wcl <= 1):
                    if verbose:
                        print('WARNING: {} has too few lines'.format(infile))
                    continue

                massh,mass,sfr,lum_ext,gtype = np.loadtxt(infile,usecols=(6,7,8,10,11),unpack=True)

                # Initialize array for hod
                hod = np.zeros(shape=len(nhs))

                for ii,nh in enumerate(nhs):
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]))
                    num = np.shape(ind)[1] 
                    if (num > 1):
                        hod[ii] = float(num)/nh 

                # Output file
                hfile = hodpath+model+cut+\
                         'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'
                print(hfile)


