import sys
import numpy as np

model = 'gp19'

infile = '/cosma5/data/durham/violeta/lines/cosmicweb/env_files/gp19/env_fractions.txt'

sns = ['39','41']
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr']
nelgs = ['DEEP2','DESI','VVDS-DEEP','eBOSS-SGC']
lse = ['Voids','Sheets','Filaments','Knots']
###########################
dum = -999.

# Initialize the arrays
alls = np.zeros(shape=(len(sns),len(nds),len(cuts),len(lse)))
alls.fill(dum)
elgs = np.zeros(shape=(len(sns),len(nds),len(cuts),len(nelgs),len(lse)))
elgs.fill(dum)

with open(infile) as ff:
    for line in ff:
        if (not line.rstrip() or line[0]=='#'):
            continue

        # Data
        array1 = line.split(':')[1]
        array2 = array1.split('[')[1]
        array3 = array2.split(']')[0]
        array = array3.split()

        # Name
        name = line.split(':')[0]
        cut = name.split('cut')[0]
        ic = cuts.index(cut)
        nd = (name.split('_')[2]).split('nd')[1]
        ind = nds.index(nd)
        sn = (name.split('sn')[1]).split('.')[0]
        isn = sns.index(sn)
        survey = name.split('_')[1]
        if (survey == 'All'):
            alls[isn,ind,ic,:] = array
        else:
            ielg = nelgs.index(survey)
            elgs[isn,ind,ic,ielg,:] = array

for isn in range(len(sns)):
    for ind in range(len(nds)):
        for ilse in range(len(lse)):
            val = alls[isn,ind,1,ilse]
            print('### {}, nd{}, sn{}: All SFR frac={:.2f}'
                  .format(lse[ilse],nds[ind],sns[isn],val))

            for  ic in range(len(cuts)):
                per = np.zeros(shape=(len(nelgs)))
                per.fill(dum)
                diffs = np.zeros(shape=(len(nelgs)))
                diffs.fill(dum)

                for ielg in range(len(nelgs)):
                    if (elgs[isn,ind,ic,ielg,ilse]>0.):
                        per[ielg] = abs(val-elgs[isn,ind,ic,ielg,ilse])*100./val
                        diffs[ielg] = abs(val-elgs[isn,ind,ic,ielg,ilse])
                        
                print('   {}cut'.format(cuts[ic]))                
                newper = per[np.where(per>dum)]
                newdiffs = diffs[np.where(diffs>dum)]

                iis = np.where(per == np.amin(newper))
                ii = iis[0][0]
                jjs = np.where(diffs == np.amin(newdiffs))
                jj = jjs[0][0]
                print('   Min diff={:.2f} {} ({:.2f} % {})'.
                      format(np.amin(newdiffs),nelgs[jj],np.amin(newper),nelgs[ii])) 
                iis = np.where(per == np.amax(newper))
                ii = iis[0][0]
                jjs = np.where(diffs == np.amax(newdiffs))
                jj = jjs[0][0]
                print('   Max diff={:.2f} {} ({:.2f} % {})'.
                      format(np.amax(newdiffs),nelgs[jj],np.amax(newper),nelgs[ii]))

print('\n Frac(Sheets+filaments) ')
for isn in range(len(sns)):
    for ind in range(len(nds)):
        for ielg in range(len(nelgs)):
            for  ic in range(len(cuts)):
                addfs = 0.
                for ilse in range(len(lse)):
                    if (lse[ilse] in ['Sheets','Filaments']):
                        if (elgs[isn,ind,ic,ielg,ilse]>0.):
                            addfs = addfs + elgs[isn,ind,ic,ielg,ilse]
                print('* {}, {}, {}cut, nd{}, sn{}'
                  .format(addfs,nelgs[ielg],cuts[ic],nds[ind],sns[isn]))

                
