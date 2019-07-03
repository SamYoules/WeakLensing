## SY 7/6/19
## Plots a single kappa power spectrum.

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys


# Use input, eg midpoint xnoisy
esttype = sys.argv[1]
maptype = sys.argv[2]

# Open kappa fits files
kxi = fits.open('maps/true_corr/kappa_noisy_lensed_true_rt70_rp100.fits.gz') \
                                                              [1].data.kappa
kinput = fits.open('maps/input/kappa_input1.fits')[1].data.I
kest = fits.open('maps/{}/{}/kappa1.fits.gz'.format(esttype, maptype)) \
                                                             [1].data.kappa
wkest = fits.open('maps/{}/{}/kappa1.fits.gz'.format(esttype, maptype)) \
                                                             [1].data.wkappa

##- Reset resolution of input maps to match estimated maps
nside=int(hp.npix2nside(kest.size))
kinput = hp.ud_grade(kinput, nside)
kxi = hp.ud_grade(kxi, nside)

kest*=-1
##- Mask area outside footprint
w = (wkest == 0)
mask_size = 1 - wkest[w].size/wkest.size
kxi[w] = hp.UNSEEN
kest[w] = hp.UNSEEN
kinput = kinput*(w)+hp.UNSEEN*(~w) 
kinput[~w]=hp.UNSEEN

##- Get Cls
cli   = hp.sphtfunc.anafast(kinput, lmax=3*nside-1)
clxi  = hp.sphtfunc.anafast(kxi, lmax=3*nside-1)
clest = hp.sphtfunc.anafast(kest, kinput, lmax=3*nside-1)

# Get the theoretical values for ells and C_ells
#th = kappa_lya.Theory()
#ell, cell = th.get_cl_kappa(2.1)

# kappa plot
P.rcParams.update({'font.size':20})
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
#P.plot(ell, cell * mask_size, label='Theory')

P.plot(cli, label=r'Input')
P.plot(clxi, label=r'$\xi - \xi_L$')
P.plot(clest, label='Noisy Auto x Input')

P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
P.xlabel('$l$', fontsize=18)
ax.set_xlim([-2, 800])
ax.set_yscale('log')
P.legend(fontsize=22, loc = "center right")
#P.title(r'DESI $\kappa {\rm \ Correlations}$')
#P.savefig('plots/kappa_corr_noise.png')
#P.close()

