## SY: Produces plots of ells vs C_ells for cross correlation of power spectra
##	of kappa healpix maps. map_1 is produced using an estimator after moving ##	mock forests with map_2, which is a gaussian realisation of kappa.

import healpy as hp
import numpy as np
import pylab as P
import sys
import fitsio
import kappa_lya
from kappa_lya import *
from astropy.io import fits

def rebin_cl(Cls):
    '''Smooth curves by rebinning Cls'''
    binwidth = 2
    l_bin = []
    Cl_bin = []
    cl_temp = 0.
    cl_weight = 0.
    for i,Cl in enumerate(Cls):
        l = i + 1
        cl_temp += Cl*l*(l+1)
        cl_weight += l*(l+1)
        if i % binwidth == 0:
            l_bin.append(l)
            Cl_bin.append(cl_temp / cl_weight)
            cl_temp = 0.
            cl_weight = 0.
    l_bin = np.array(l_bin)
    Cl_bin = np.array(Cl_bin)
    return l_bin, Cl_bin


## Get the maps and put them into the right format

# input map names 
map_1 = 'kappa-gaussian-true-v3.fits.gz'
#map_1 = 'kappa-wsky.fits.gz'
map_2 = 'kappa_input_gaussian.fits'
#map_3 = 'kappa-gaussian-v3.fits.gz'

map1 = fits.open(map_1)[1].data.kappa
wmap1 = fits.open(map_1)[1].data.wkappa
map2 = fits.open(map_2)[1].data.I
#map3 = fits.open(map_3)[1].data.kappa

# degrade kappa_true map to match number of pixels of other
NSIDE=int(hp.npix2nside(map1.size))
alm = hp.sphtfunc.map2alm(map2, lmax=3*NSIDE-1)
map2 = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Trim kappa_true to match eBOSS footprint
w = (wmap1 == 0)
mask_size = 1 - wmap1[w].size/wmap1.size
map1[w] = hp.UNSEEN
map2[w] = hp.UNSEEN
#map3[w] = hp.UNSEEN


## Get Cls from maps for auto and cross correlations

# kappa_v3 (estimator) auto-correlation
Cls1=hp.sphtfunc.anafast(map1, lmax=3*NSIDE-1)
#l1, Cls1 = rebin_cl(Cls1)
l1 = np.arange(Cls1.size)+1

# kappa_input.fits (true) auto-correlation
Cls2=hp.sphtfunc.anafast(map2, lmax=3*NSIDE-1)
#l2, Cls2 = rebin_cl(Cls2)

# true/estimator cross-correlation
Cls12 = hp.sphtfunc.anafast(map1, map2, lmax=3*NSIDE-1)
#l12, Cls12 = rebin_cl(Cls12)

# true/noisy cross-correlation
#Cls23 = hp.sphtfunc.anafast(map2, map3, lmax=3*NSIDE-1)
#l23, Cls23 = rebin_cl(Cls23)

# Get the theoretical values for ells and C_ells
th = kappa_lya.Theory()
ell, cell = th.get_cl_kappa(2.1)

P.rcParams.update({'font.size': 22})

# signal as a function of scale
#rho = Cl12/np.sqrt(Cl1*Cl2)
P.ion()
#P.figure()
#P.plot(ell, cell, label='Theory')
P.plot(ell, cell * mask_size, label='Theory')
P.plot(l1, Cls1, label='$\kappa_{\kappa \kappa}$')
P.plot(l1, Cls2, label='$\kappa_{tt}$')
P.plot(l1, Cls12, label='$\kappa_{t \kappa}$')
#P.plot(l1, Cls23, label='$\kappa_{tn}$')
P.plot(l1, Cls1-Cls2, label='$\kappa_{\kappa \kappa} - \kappa_{tt}$')
#P.plot(l1, rho, label='$\\rho$')
P.xlim(10,200)
P.yticks(np.arange(0.0, 2.0e-8, step=0.5e-8))
P.ylabel('$C_{l}^{\kappa \kappa}$', fontsize=22)
P.xlabel('$l$', fontsize=22)
P.legend()
#P.title(r'${\rm Cross \ Correlations \ of \ C}_{l}s \rm \ Power \ Spectrum$')
#P.savefig('plots/kappa_whole_sky.png')
#P.close()

