#! /bin/tcsh -f

set name    = selections
set logpath = /gpfs/data/$user/Junk/$name


time bsub -P durham -n 1 -q cordelia -J "$name" -o $logpath.%J.%I run.csh 


