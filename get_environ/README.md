# README #

At the moment, we are able to get the environment for objects in snapshot 39 (MS-W7 simulation), given an ascii table with x,y,z coordinates.

STEPS:

1. **xyz2ascii.py** Output the x,y,z coordinates of galaxies into an ascii file named: [cut]\_sn[#]\_[all/centrals/satellites].dat 

2. Copy to taurus the resulting ascii table (maybe use rsync?):
   scp [name] vgonzalez@taurus.ft.uam.es:/home2/vgonzalez/cosmicweb/gal_files/

3. Run Wiguang's code with the script:
   vgonzalez@taurus.ft.uam.es:/home2/vgonzalez/cosmicweb/get_env.py

## Weiguang's code  ##

**getE.py** reads a lits of x,y,z and gets the environment. This is a copy of the module placed at taurus:

/data1/users/weiguang/Violeta/get_environments/getE.py

This code needs numpy, to install it locally:
'''pip install --user numpy'''

This code runs over the Vweb data stored in the files:

/data1/users/weiguang/Violeta/VPweb_data/VP_039_DM_MHD.000256.Vweb

/data1/users/weiguang/Violeta/VPweb_data/VP_041_DM_MHD_512.000512.Vweb

## Particles ##

**particles2ascii.py** Example program to read the hdf5 files and produce ascii files.

The snapshot with particle information was copied from cosma (/gpfs/data/Millgas/data/dm/500/snapdir_039/ 178GB) to taurus:

> rsync -avzhe ssh snapdir_039/ vgonzalez@taurus.ft.uam.es:/data1/users/weiguang/Violeta/

