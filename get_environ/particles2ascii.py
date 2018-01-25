import os, sys
import h5py
import numpy as np

basedir = '/data1/users/weiguang/Violeta/' # taurus
basedir = '/gpfs/data/Millgas/data/dm/500/' # cosma
outdir = '/gpfs/data/violeta/Junk/dm'

isnap = 39 #z=1.
basename = basedir+"snapdir_%03d/500_dm_%03d" % (isnap, isnap)

for ifile in range(2):#512):

    # Read coordinates from one file
    fname = "%s.%d.hdf5" % (basename,ifile)
    if (not os.path.isfile(fname)):
        print "NOT found: %s" % fname
        continue

    f = h5py.File(fname, "r")
    # Coordinates = comoving position (Mpc/h)
    pos_x = f["PartType1/Coordinates"].value[:,0]
    pos_y = f["PartType1/Coordinates"].value[:,1]
    pos_z = f["PartType1/Coordinates"].value[:,2]

    # Velocities = sqrt(a)*dx/dt  ; a = 1/(1+z)
    # Peculiar velocity = sqrt(a)*Velocities (km/s)
    vel_x = f["PartType1/Velocities"].value[:,0]
    vel_y = f["PartType1/Velocities"].value[:,1]
    vel_z = f["PartType1/Velocities"].value[:,2]
    f.close()
    
    #Write the information into an ascii file
    output = "%s.%d.dat" % (outdir,ifile)

    tofile = zip(pos_x,pos_y,pos_z,vel_x,vel_y,vel_z)
    with open(output, 'w') as outf:                            
        outf.write('# x,y,z (Mpc/h), vel_x,vel_y,vel_z \n')                    
        np.savetxt(outf,tofile,fmt=('%.8f'))    
        outf.closed      

          

