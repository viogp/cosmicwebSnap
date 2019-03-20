Copied from ~/lines/desi_hod_o2/o2/cuts/contributions/


11th March from Johan
Can you plot LOII unextincted vs observed g magnitude and 
LOII extincted vs observed g magnitude. 
and the same vs. r-band at redshift 0.83.  at z=0.83, lambda_[OII] = 3727*1.83 =  6820.41 AA in the middle of the r-band and the NUV flux (2000-2500AA) is in the g-band. So if you use the r-band for selection band, you may have less difference due to extinction law ?  

Another question (or a different way to re phrase the previous statement) is how are g and r band computed: do they take line fluxes into account ?  It may be a bug as well :  - at z=0.83, the curves do not seem to match the total number:  --  maybe you can cut brighter in the g band 22.8 or 22.7 to see ?  -- one of the red dashed curves has the right shape. So may also be a unit problem in x or y ? i.e. do you predict between 200 and 300 / deg2 over the range 0.6-1.2. (the measured curve looks like 200/deg2 after a rough calculation by hand if dlogL=0.05) 

 Cheers, J