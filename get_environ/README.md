# README on 'get_environ' folder #

We get the environment for objects in snapshots 39, z=1, and 41, z=0.83 (MS-W7 simulation), given an ascii table with x,y,z coordinates.

STEPS:

1. **../selections/mass_selections.py sfr_selections.py elg_mselections.py elg_sselections.py** Output the x,y,z coordinates of galaxies into an ascii file named: [cut]\_sn[#]\_[all/centrals/satellites].dat 

**xyz2ascii.py** Used for producing a halo catalogue, before Weiguang had access to the simulation particles.

2. Copy to taurus the resulting ascii table (maybe use rsync?):
   scp [name] vgonzalez@taurus.ft.uam.es:/home2/vgonzalez/cosmicweb/gal_files/

3. Run Weiguang's code with the script for both sn39 and sn41 (note that the path to Vweb should be changed accordingly):
   vgonzalez@taurus.ft.uam.es:/home2/vgonzalez/cosmicweb/get_env.py

   > nohup python get_env.py >& out &

   > ps ax | grep get_env

4. Copy into cosma the environment files:
   scp -r vgonzalez@taurus.ft.uam.es:/home2/vgonzalez/cosmicweb/env_files/ /cosma6/data/dp004/dc-gonz3/lines/cosmicweb/env_files

## 'cosmicweb' folder at taurus ##

**get_env.py** This program reads the positions of galaxies from the files 'gal_files/*_sn*.dat' and find the environment. This information is stored in 'env_files'. Those files that have already been processed or are not being process in the run should be moved to the folder 'gal_files/processed/'. A copy of this program can also be found in this cosma folder.


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

> rsync -avzhe ssh snapdir_039/ vgonzalez@taurus.ft.uam.es:/data1/users/weiguang/Violeta/

From /cosma5/data/GRID/jch/millgas

    rsync -avzhe ssh snapdir_041/ vgonzalez@taurus.ft.uam.es:/data1/users/weiguang/Violeta/snapdir_041/

