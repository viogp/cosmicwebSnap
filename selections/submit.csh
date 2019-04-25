#! /bin/tcsh -f

set script = run.csh
set jobname = selections
set logpath = /gpfs/data/${user}/Junk
set logname = ${logpath}/${jobname}.%A.%a.log

# Construct a batch script and submit it to SLURM as a single job -
# the scritpt consists of the Slurm header below followed by the 
# contents of ${script}
cat << EOF - ${script} | sbatch 
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${jobname}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 72:00:00

# Run script follows
EOF



