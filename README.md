# README #

At the moment, we are able to get the environment for objects in snapshot 39 (MS-W7 simulation), given an ascii table with x,y,z coordinates.

## Models ##
Relevant plots can be found also in the Galform wiki: [for the Spin part](http://galform.pbworks.com/w/page/123049548/Spin-up%20on%20and%20off) and for the [exploration of the passive fraction plot](http://galform.pbworks.com/w/page/123760410/The%20passive%20fraction). 

### From v2.6.0 to v2.7.0 ###
Comparing:
* /gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/gp15newmg/
* /gpfs/data/violeta/Galform_Out/v2.7.0/stable/MillGas/gp17

v2.6.0: 
SMBH_spinup =  0 
SMBH_delta_m = 1 (default value for self-gravity, 0 for warp)

v2.7.0: 
SMBH_merger_spinup = 1 ; SMBH_accretion_spinup = 0
SMBH_delta_m = 2 (default value for self-gravity, 1 for warp)

Other new parameters in gp17:
* spin_from_main_subhalo_only = F
* SMBH_default_spin =   1.00000E-01
* spin_magnitude_only = T 

### From gp17 to gp17.spin ###

Differnces between these 2 models:

* SMBH_accretion_spinup = 0 (gp17), 2 (gp17.spin)
* SMBH_Mseed =   0 (gp17) ; 10 (gp17.spin)

### New model ###
As the ELG results appeared to be very sensitive to the number of star forming satellite galaxies, we have tried to better reproduced the passive fraction at z=0, as a first attempt to deal with this issue.

## Cosmic web information ##
Cosmic web information will be taken from taurus.

The snapshot with particle information was copied from cosma (/gpfs/data/Millgas/data/dm/500/snapdir_039/ 178GB) to taurus:

> rsync -avzhe ssh snapdir_039/ vgonzalez@taurus.ft.uam.es:/data1/users/weiguang/Violeta/

From /cosma5/data/GRID/jch/millgas
> rsync -avzhe ssh snapdir_041/ vgonzalez@taurus.ft.uam.es:/data1/users/weiguang/Violeta/snapdir_041/
