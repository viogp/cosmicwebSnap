#! /bin/tcsh -f

set exec = 'shuffleMSI.py'

set Testing = False

set nom = 'shuffle'
set logpath = /cosma6/data/dp004/dc-gonz3/Junk/
set logname = ${logpath}/props_$nom.%A.%a.log
set job_file = ${logpath}/${nom}.job

if ($Testing == 'True') then
    echo 'Testing'
    echo $exec
    python3 $exec
else
    cat > $job_file <<EOF
#! /bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${nom}
#SBATCH -o ${logname}
#SBATCH -p cosma6
#SBATCH -A dp004
#SBATCH -t 12:00:00
#
echo $exec
python3 $exec
#
EOF

    sbatch $job_file
    rm $job_file
endif

echo 'All submited'
