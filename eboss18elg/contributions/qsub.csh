#! /bin/tcsh -f

#set noms = (deep2 vvdsdeep vvdswide eboss desi)
set noms = ( deep2)

#set model = 'gp18'  
#set model = 'gp18.font'  
set model = 'gp18.starvation'  

foreach nom ($noms)
    set logpath = /gpfs/data/$user/Junk/contributions_$nom  

    time bsub -P durham -n 1 -q cordelia -J "$nom" -o $logpath.%J.%I run.csh $nom $model
	    
    # Testing
    #./run.csh $nom $model
end

echo 'All submited'
