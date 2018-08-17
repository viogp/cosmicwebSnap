#! /bin/tcsh 

set testing = true

set logdir = /gpfs/data/$user/Junk/

set path = /cosma/home/violeta/lines/cosmic_web/bias/DM/
set p2p = ${path}hdf5_cute_box/paramfiles/millgas/
#set exec = ${path}hdf5_cute_box/CUTE_box
set exec = ${path}CUTE/CUTE_box/CUTE_box

foreach iparam (params_sn41.txt params_sn39.txt)
    if ( "$testing" == true )then
	./run_cute.csh $exec ${p2p}$iparam
    else
	set logfile = ${logdir}$iparam
	echo bsub -P durham -n 1 -q cordelia -J "$iparam" -o "$logfile.%I" run_cute.csh $exec ${p2p}$iparam 
    endif
end
