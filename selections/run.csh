#!/bin/tcsh -f

#set exec = mass_cum.py
#set exec = sfr_cum.py
#set exec = elg_mselections.py
set exec = elg_sselections.py
#set exec = mass_selections.py
#set exec = sfr_selections.py

python $exec

echo 'The end'

