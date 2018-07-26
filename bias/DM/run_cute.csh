#!/bin/tcsh
limit datasize unlimited
limit stacksize unlimited
limit coredumpsize 0k
limit vmemoryuse 5000m

############################################################################
## Read in INPUT ARGUMENTS   
############################################################################
set nargs = $#argv
echo "number of script arguments = $nargs"
if( $nargs == 2 )then
    set exec = $1
    echo 'CUTE =' $exec
    set paramfile = $2
    echo 'Parameter file =' $paramfile
else
    echo 'ERROR - USAGE: CUTE paramfile' 
    exit 1
endif

$exec $paramfile


