## SY 14/10/18
## Plots auto-correlations of kappa power spectra for noisy maps with R_t max. Clips off outliers which were creating odd results.

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys
import os

def read_kappa_input(nside=256):
    kappa_input = fits.open('kappa_input_gaussian.fits')[1].data.I.ravel()
    #kappa_input = fits.open('kappa_input_noiseless.fits')[1].data.I.ravel()
    alm_input = hp.map2alm(kappa_input, lmax=3*nside-1)
    kappa_input = hp.alm2map(alm_input, nside=nside)
    return kappa_input

def rebin_cl(Cls):
    '''Smooth curves by rebinning Cls'''
    binwidth = 20
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

# Open kappa fits files
nside=256
lw=0.5
kappa_input = read_kappa_input(nside=nside)
ell = np.arange(3*nside)

P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)

data = fits.open('kappa-noiseless-70-40.fits.gz')[1].data
#data = fits.open('kappa-noiseless-70-100.fits.gz')[1].data
kappa = data.kappa
mask = data.wkappa!=0
mask &= (data.kappa>np.percentile(data.kappa[mask], 0.5)) & \
            (data.kappa<np.percentile(data.kappa[mask], 99.5))
kappa[~mask]=hp.UNSEEN
f_sky = sum(mask)/mask.size
kappa_input_masked = kappa_input*(mask)+hp.UNSEEN*(~mask) 
cl_data = hp.anafast(kappa)
cl_input = hp.anafast(kappa_input_masked)
cl_cross = hp.anafast(kappa, kappa_input_masked)
lcross, clcross = rebin_cl(cl_cross)
P.plot(ell, cl_data, lw = 2, label=r'$\kappa_{ee}$')
P.plot(ell, cl_input, lw = 2, label=r'$\kappa_{tt}$')
P.plot(lcross, clcross, lw = 2, label=r'$\kappa_{et}$')

#ax.set_yscale('log')
ax.set_yscale('symlog', linthreshy=1e-7)

P.ylabel('$C_l^{\kappa \kappa}$', fontsize=22)
P.xlabel('$l$', fontsize=22)
leg = P.legend(fontsize=22, bbox_to_anchor=(0.95,0.55), loc = "center right")
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
#P.title(r'${\rm Noiseless} \ \kappa {\rm \ Auto-Correlations}$')


