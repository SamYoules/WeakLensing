## SY 17/8/18
## Plots auto-correlations of kappa power spectra for noisy maps with varying R_p max to see which is optimal.

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys


# Open kappa fits files
wkxi = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.wkappa
kxi = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.kappa
knoisy1 = fits.open('kappa-noiseless-40-10.fits.gz')[1].data.kappa
knoisy2 = fits.open('kappa-noiseless-40-20.fits.gz')[1].data.kappa
knoisy3 = fits.open('kappa-noiseless-40-30.fits.gz')[1].data.kappa
knoisy4 = fits.open('kappa-noiseless-40-40.fits.gz')[1].data.kappa
knoisy5 = fits.open('kappa-noiseless-40-50.fits.gz')[1].data.kappa
knoisy6 = fits.open('kappa-noiseless-40-60.fits.gz')[1].data.kappa
knoisy7 = fits.open('kappa-noiseless-40-70.fits.gz')[1].data.kappa


# Get nside (to use for all kappa maps)
NSIDE=int(hp.npix2nside(kxi.size))

# Mask off areas outside eBOSS footprint
w = (wkxi == 0)
mask_size = 1 - wkxi[w].size/wkxi.size
kxi[w] = hp.UNSEEN
knoisy1[w] = hp.UNSEEN
knoisy2[w] = hp.UNSEEN
knoisy3[w] = hp.UNSEEN
knoisy4[w] = hp.UNSEEN

# Get the C_ells for auto- and cross-correlations
Cls1=hp.sphtfunc.anafast(knoisy1, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(knoisy2, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(knoisy3, lmax=3*NSIDE-1)
Cls4=hp.sphtfunc.anafast(knoisy4, lmax=3*NSIDE-1)
Cls5=hp.sphtfunc.anafast(knoisy5, lmax=3*NSIDE-1)
Cls6=hp.sphtfunc.anafast(knoisy6, lmax=3*NSIDE-1)
Cls7=hp.sphtfunc.anafast(knoisy7, lmax=3*NSIDE-1)
#Cls5=hp.sphtfunc.anafast(kxi, lmax=3*NSIDE-1)

# Get the theoretical values for ells and C_ells
th = kappa_lya.Theory()
ell, cell = th.get_cl_kappa(2.1)

# kappa plot
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)

P.plot(Cls1, label='d1*d2 rtmax=40 rpmax=10')
P.plot(Cls2, label='d1*d2 rtmax=40 rpmax=20')
P.plot(Cls3, label='d1*d2 rtmax=40 rpmax=30')
P.plot(Cls4, label='d1*d2 rtmax=40 rpmax=40')
P.plot(Cls5, label='d1*d2 rtmax=40 rpmax=50')
P.plot(Cls6, label='d1*d2 rtmax=40 rpmax=60')
P.plot(Cls7, label='d1*d2 rtmax=40 rpmax=70')
#P.plot(Cls5, label=r'$\xi - \xi_L$')

P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
P.xlabel('$l$', fontsize=18)
ax.set_xlim([-2, 800])
ax.set_yscale('log')
leg = P.legend(fontsize=18, bbox_to_anchor=(0.3,0.3), loc = "center left")
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)

P.title(r'${\rm Noiseless} \ \kappa {\rm \ Auto-Correlations}$')
#P.savefig('plots/kappa_corr_noise.png')
#P.close()
