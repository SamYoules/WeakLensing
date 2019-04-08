## SY 19/7/18
## Plots auto-correlations of kappa power spectra for input and estimated maps
## (true-correlation, noisy and noiseless mocks).

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys


# Open kappa fits files
ktrue = fits.open('kappa_input.fits')[1].data.I
wkxi = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.wkappa
kxi = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.kappa
knoisy = fits.open('kappa-gaussian-v3.fits.gz')[1].data.kappa
kquiet = fits.open('kappa-noiseless.fits.gz')[1].data.kappa


# Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(kxi.size))
alm = hp.sphtfunc.map2alm(ktrue, lmax=3*NSIDE-1)
ktrue = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Mask off areas outside eBOSS footprint
w = (wkxi == 0)
mask_size = 1 - wkxi[w].size/wkxi.size
kxi[w] = hp.UNSEEN
ktrue[w] = hp.UNSEEN
knoisy[w] = hp.UNSEEN
kquiet[w] = hp.UNSEEN

# Get the C_ells for auto- and cross-correlations
Cls1=hp.sphtfunc.anafast(ktrue, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(kxi, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(knoisy, lmax=3*NSIDE-1)
Cls4=hp.sphtfunc.anafast(kquiet, lmax=3*NSIDE-1)

# Get the theoretical values for ells and C_ells
th = kappa_lya.Theory()
ell, cell = th.get_cl_kappa(2.1)

# kappa plot
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
P.plot(ell, cell * mask_size, label='Theory')

P.plot(Cls1, label='Input')
P.plot(Cls2, label='True Correlation')
P.plot(Cls3, label='Noisefull')
P.plot(Cls4, label='Noiseless')

P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
P.xlabel('$l$', fontsize=18)
ax.set_xlim([-2, 210])
ax.set_ylim([10e-10, 10e-2])
ax.set_yscale('log')
P.legend(fontsize=18, loc = "center right")
P.title(r'$\kappa {\rm \ Correlations}$')
#P.savefig('plots/kappa_corr_noise.png')
#P.close()


