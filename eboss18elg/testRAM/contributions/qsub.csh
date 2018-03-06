#! /bin/tcsh -f

#set noms = (deep2 vvdsdeep vvdswide eboss desi)
set noms = ( deep2)

#set model = 'MillGas/gp17/'  
#set model = 'MillGas/gp17.spin/'
#set model = 'MillGas/gp17.spin.ramt0.01/'
#set model = 'MillGas/gp17.spin.ramt0.01.griffinBH/'
set model = 'MillGas/gp17.spin.ramt0.01.griffinBH.stb075/'
#set model = 'MillGas/gp17.spin.ramt0.01.stabledisk0.75.ac085/'


foreach nom ($noms)
    set logpath = /gpfs/data/$user/Junk/contributions_$nom  

    time bsub -P durham -n 1 -q cordelia -J "$nom" -o $logpath.%J.%I run.csh $nom $model
	    
    # Testing
    #./run.csh $nom $model
end

echo 'All submited'
