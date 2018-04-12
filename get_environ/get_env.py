import os, sys
sys.path.append('/data1/users/weiguang/Violeta/get_environments/')
import glob, time
import numpy as np
from getE import return_env

'''
def return_env(Pos, threshold=0.1, webf="../VPweb_data/VP_039_DM_MHD.000256.Vweb", Vweb=True):
        """
    Parameters:
    ----------
    Pos: must be 3D data array in shape of (N,3) in units of Mpc/h
    threshold: the \lambda_th for classifying environments, default 0.1
    webf: the file location for Vweb output, default: ../VPweb_data/VP_039_DM_MHD.000256.Vweb
    Vweb: If you want to use Vweb method, Set False for the Pweb method.
    
    Return:
    An 1D array indicates the position environments: 0 Void; 1 sheet; 2 filament; 3 Knots.
    """
'''
#----LSE parameters---------------
boxsize = 500.  # Side of the simulation box (Mpc/h)
threshold = 0.1  # For classifying environments
pathVweb = '/data1/users/weiguang/Violeta/VPweb_data/VP_039_DM_MHD.000256.Vweb'
Vweb = True 
#---------------------------------

# Read files from gal_files
# The name should be: *.dat
files = glob.glob('gal_files/'+'*.dat')
for infile in files:
    time0 = time.clock()
    
    pos = np.loadtxt(infile, usecols=(0,1,2), unpack= True)
    pos = np.transpose(pos)
    print 'Spent time w loadtxt = ',time.clock()-time0, np.shape(pos)
    time0 = time.clock()

    count = -1
    with open(infile) as ff:
        for line in ff:
            if (line[0] != '#'):
                count += 1 #; print count           
                x = float(line.split()[0])
                y = float(line.split()[1])
                z = float(line.split()[2])
                newpos = np.array([x,y,z])
                newpos = np.expand_dims(newpos, axis=0)

                if (count == 0):
                    pos = newpos
                else:
                    pos = np.vstack((pos, newpos))
    print 'Spent time line by line = ',time.clock()-time0, np.shape(pos)

    env = return_env(pos, boxsize, threshold = threshold, \
                     webf = pathVweb, Vweb = Vweb)
    
    # Write output file with positions and environment
    env2f = np.expand_dims(env, axis=1)
    tofile = np.concatenate((pos,env2f),axis=1)
    
    nom = infile.split('/')[-1:] # File name
    outfile = 'env_files/'+nom[0]
    if os.path.exists(outfile):
        os.remove(outfile)

    outf = open(outfile,'w')
    outf.write('# x,y,z (Mpc/h), environments: 0 Void; 1 sheet; 2 \
filament; 3 Knots \n')
    np.savetxt(outf, tofile, fmt='%.5f %.5f %.5f %1i')
    outf.close()
    
    print('Output: ',outfile)
