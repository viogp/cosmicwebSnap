#! /bin/tcsh -f

#set noms = (deep2 vvdsdeep vvdswide eboss desi)
set noms = ( desi eboss f16)

set model = 'gp19'  
#set model = 'gp19.font'  
#set model = 'gp19.starvation'  

set script = run.csh

foreach nom ($noms)
    set logname = /cosma5/data/durham/$user/Junk/contributions_$nom.%A.%a.log  

    #time bsub -P durham -n 1 -q cordelia -J "$nom" -o $logpath.%J.%I run.csh $nom $model

    cat << EOF - ${script} | sbatch
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${nom}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 48:00:00
#
# Set parameters
set nom = ${nom}
set model = ${model}

# Run script follows
EOF

	    
    # Testing
    #./run.csh $nom $model
end

echo 'All submited'
