import numpy as np
import sys
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

rmin = 0.02

sn =['41','39']
path ='/gpfs/data/violeta/mr7corr/'
froot=path+'xi_real_sn0'

#Plot
fig = plt.figure(figsize=(8.,9.))
xmin = -2. ; xmax = 2.
xtit = "${\\rm log}_{10}(\\rm{r/Mpc}\, h^{-1})$"
ymin = -3. ; ymax = 4
ytit = "${\\rm log}_{10}\\xi (\\rm{r})$" 
ax = plt.subplot()
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

#Read the data and plot it
for s in sn:
    ff = froot+s+'.txt'
    r,xi,error,dd = np.loadtxt(ff,unpack=True)

    ind=np.where((r>rmin) & (xi>0.))
    lerr = error[ind]*np.log10(np.exp(1))/r[ind]

    #ax.errorbar(np.log10(r[ind]),np.log10(xi[ind]),yerr=lerr,label=s)
    ax.errorbar(np.log10(r[ind]),np.log10(xi[ind]),label=s)


# Plot Test
#ff='/cosma/home/violeta/lines/cosmic_web/bias/DM/hdf5_cute_box/test/corr128.dat'
#ff='/gpfs/data/violeta/CUTE/CUTE_box/test/corr128.dat'
ff='/cosma/home/violeta/lines/cosmic_web/bias/DM/CUTE/CUTE_box/test/corr128.dat'
r,xi,error,dd = np.loadtxt(ff,unpack=True)
ind=np.where((r>rmin) & (xi>0.))
lerr = error[ind]*np.log10(np.exp(1))/r[ind]
ax.errorbar(np.log10(r[ind]),np.log10(xi[ind]),label=s)

# Save figure
plotfile = path+'xiDM.pdf'
fig.savefig(plotfile)
print('Output: ',plotfile)
