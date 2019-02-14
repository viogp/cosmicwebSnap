#!/bin/tcsh -f

set exec =  ${nom}'_contributions.py'

echo $exec $model
python $exec $model


