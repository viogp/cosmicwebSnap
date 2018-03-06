#!/bin/tcsh -f

#set exec = ssfr_m.py
#set exec = lf_K.py
#set exec = lf_bJ.py
#set exec = passive.py
set exec = lf4_cuts.py

python $exec

echo 'The end'

