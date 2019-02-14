#!/bin/tcsh

# ifort reformat125.f90 -L/gpfs/data/violeta/f90lib  -lvio
set exec = /cosma5/data/durham/violeta/CUTE/CUTE_box/CUTE_box

limit datasize unlimited
limit stacksize unlimited
limit coredumpsize 0k
limit vmemoryuse 5000m
###########################################################################
# Read in INPUT ARGUMENTS   ~/lines/desi_hod_o2/o2/cuts/xi
###########################################################################

echo "Submitting job to Cordelia, param_file =" $parfil

$exec $parfil


