#!/bin/tcsh -f

#set exec = lf_K.py
#set exec = lf_bJ.py
#set exec = passive.py
set exec = ssfr_m.py
#set exec = lf2_cuts.py
#set exec = decam_eboss.py

python $exec

echo 'The end'

