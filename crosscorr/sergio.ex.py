# https://corrfunc.readthedocs.io/en/master/
# https://arxiv.org/pdf/1911.03545.pdf
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

def _RR(binfile, lbox, N1, N2=None):
    if N2 is None: N2=N1 #Autocorrelation Function
    n2 = N2/float(lbox**3)

    vshell = (binfile[1:]**3-binfile[:-1]**3)*4/3.*np.pi
    
    return N1*n2*vshell
# _var: _ hints that var is intended for internal use (defined in PEP 8).


# For DeLucia+2007 SAM
boxsize = 62.5
NBin = 30
log_rmin = -2
log_rmax = 1

# The number of OpenMP threads to use if OpenMP enabled during compilation.
nthreads = 4 # Change this on a cluster

# Binfile: A sequence of ``r`` values that provides the bin-edges
binfile = np.logspace(log_rmin,log_rmax,NBin+1,True)
bin_middle = (binfile[1:] + binfile[:-1])/2.

x,y,z = np.loadtxt('DL.txt').T

#--------Auto-correlations
# 3-D correlation function in a periodic cosmological box
# print(help(xi))
cf = xi(boxsize, nthreads, binfile, x,y,z)
# cf: A numpy structured array containing, for each binin binfile,
#     [rmin, rmax, ravg, xi, npairs, weightavg] 
_CF = cf['xi']

# 3-D pair-counts corresponding to the real-space correlation function
#print(help(DD))
autocorr = True
cf1 = DD(autocorr, nthreads, binfile, x,y,z, boxsize=boxsize)
# cf1: A numpy structured array containing, for each binin binfile,
#      [rmin, rmax, ravg, npairs, weightavg]
print('Expected array of 0s: {} \n'.format(cf1['npairs']-cf['npairs']))

#--------Cross-correlations
# Compute the 2 CCF and the full CF to compare. For the DD, the input parameters are:
mask1 = x > boxsize/2.
mask2 = x <= boxsize/2.

autocorr = False
DD1 = DD(autocorr, nthreads, binfile,
         x[mask1], y[mask1], z[mask1], 
         X2=x, Y2=y, Z2=z,boxsize=boxsize)
DD2 = DD(autocorr, nthreads, binfile,
         x[mask2], y[mask2], z[mask2], 
         X2=x, Y2=y, Z2=z,boxsize=boxsize)

# Get the cross-correlation functions
_CF1 = DD1['npairs']/_RR(binfile,boxsize,len(x),len(x[mask1])) - 1
_CF2 = DD2['npairs']/_RR(binfile,boxsize,len(x),len(x[mask2])) - 1

#------------Test
# Check that the numbers of pairs of the full CF are equal to
# the sum of the two CCFs
print('#r \t DD1 \t DD2 \t DD1+DD2 \t CC_total \t (DD1+DD2)-CC_total')
for i in range(NBin):
    print('{:.2f}\t{}\t{}\t{}\t{}\t{}'.format(bin_middle[i],
                                              DD1['npairs'][i],DD2['npairs'][i],
                                              DD1['npairs'][i]+DD2['npairs'][i],
                                              cf['npairs'][i],
                                              DD1['npairs'][i]+DD2['npairs'][i]-
                                              cf['npairs'][i]))

# Plot
# The auto correlation is expected to be between the 2 halves
plt.figsize=(6.5,6.5)
plt.loglog(bin_middle, _CF, 'k', label=r'Full correlation function')
plt.loglog(bin_middle, _CF1, 'r', label=r'CCF first half')
plt.loglog(bin_middle, _CF2, 'b', label=r'CCF second half')
# r+string: all escape codes to be ignored
 
plotfig = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/CCF.ex.pdf'
plt.legend(loc='best')
plt.ylabel(r'$\xi(r)$')
plt.xlabel(r'$r (h^{-1}{\rm Mpc})$')
plt.savefig(plotfig)
print('Output: {}'.format(plotfig))
