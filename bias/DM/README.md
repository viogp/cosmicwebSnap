Snapshot: 39, 41 (MS-W7)
rmin,rmax (Mpc/h) = 0.02, 60.

## Information from Nuala

I am using the DM particle information for MS-W7, which are in HDF5 files in here:
/cosma5/data/GRID/jch/millgas/snapdir_###

I edited the CUTE code to be able to read in the HDF5 files for MS-W7 and do the subsampling. Peder asked me a few months ago to measure the DM correlation function for several snapshots, so that’s why I had it all set up already! I just checked and apparently it took 4 days for the code to finish then (!) so hopefully I will send you something next week.

 the first time I tried running them (downsampling the particles by a factor of 100 #sampling_rate#) it took too long and the job timed out. I just ran it again downsampling by a factor of 1000, and it finished much quicker. I can try running it again with more particles but these look pretty smooth, so maybe they are ok for what you need. 

My CUTE stuff is in:
/cosma/home/nmccull/code/CUTE-1.2/CUTE_box

currently I have NB_R=84, R_MAX = 65, N_LOGINT = 12. I have run one snapshot with this setting (snapshot 44), and it is in:
/gpfs/data/nmccull/mr7corr/violeta/

DEFINEOPTIONS += -D_LOGBIN uncommented (log binning for the correlation functions)
ats
-----------------------
~/lines/cosmic_web/bias/DM/hdf5_cute_box/paramfiles/millgas
params_sn39.txt
params_sn41.txt