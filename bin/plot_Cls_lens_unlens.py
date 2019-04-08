## SY 8/12/18
## Plots cross-correlations of qso-lya kappa power spectrum.

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys

def rebin(Cls,ellmin=40, ellmax=768, nell=7):
    '''Smooth curves by rebinning Cls'''
    ell = np.arange(Cls.size)
    weights = 2.*ell+1.
    if not ellmin:
        ellmin = min(ell)
    if not ellmax:
        ellmax = max(ell)
    if nell==0:
        nell = ellmax-ellmin

    w = (ell>=ellmin)&(ell<=ellmax)
    index = np.floor( (ell[w]-ellmin)*1./(ellmax-ellmin)*nell ).astype(int)
    well = np.bincount( index, weights=weights[w])
    sell = np.bincount( index, weights=weights[w]*ell[w])
    scl  = np.bincount( index, weights=weights[w]*Cls[w])
    ell = sell/well
    Cls = scl/well
        
    return ell, Cls, well 


## Open kappa fits files
wkunlensed = fits.open('kappa-xcut-unlensed-true-rt70-rp100.fits.gz')[1].data.wkappa
kunlensed = fits.open('kappa-xcut-unlensed-true-rt70-rp100.fits.gz')[1].data.kappa
wklensed = fits.open('kappa-xcut-lensed-true-rt70-rp100.fits.gz')[1].data.wkappa
klensed = fits.open('kappa-xcut-lensed-true-rt70-rp100.fits.gz')[1].data.kappa
kinput = fits.open('est_maps_xcut/kappa_input1.fits')[1].data.I


# Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(klensed.size))
alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Mask off areas outside survey footprint

mask = wklensed==0
#mask_size = 1 - wklensed[mask].size/wklensed.size
mask &= (klensed>np.percentile(klensed[mask], 0.5)) & \
           (klensed<np.percentile(klensed[mask], 99.5))
klensed[mask]=hp.UNSEEN
kunlensed[mask]=hp.UNSEEN
kinput_masked = kinput*(~mask)+hp.UNSEEN*(mask) 

# Get the C_ells for auto- and cross-correlations
Cls1=hp.sphtfunc.anafast(klensed, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(kunlensed, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(kinput_masked, lmax=3*NSIDE-1)

Cls13=hp.sphtfunc.anafast(klensed, kinput_masked, lmax=3*NSIDE-1)
Cls23=hp.sphtfunc.anafast(kunlensed, kinput_masked, lmax=3*NSIDE-1)

P.rcParams.update({'font.size':20})
# kappa plot
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)

P.plot(Cls1, label=r'$\langle \kappa_{lensed} \kappa_{lensed} \rangle$')
P.plot(Cls2, label=r'$\langle \kappa_{unlensed} \kappa_{unlensed} \rangle$')
P.plot(Cls4, label=r'$\langle \kappa_{lya} \kappa_{lya} \rangle$')
P.plot(Cls3, label=r'$\langle \kappa_{input} \ \kappa_{input} \rangle$')
#P.plot(Cls23, label=r'$\langle \kappa_{input} \ \kappa_{unlensed} \rangle$')
#P.plot(Cls13, label=r'$\langle \kappa_{input} \ \kappa_{lensed} \rangle$')


P.ylabel('$C_l^{\kappa \kappa}$', fontsize=24)
P.xlabel('$l$', fontsize=24)
ax.set_xlim([-2, 800])
#ax.set_ylim([-4e-6, 1.4e-6])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
P.legend(fontsize=22, loc = "upper right")

#P.title(r'DESI Noiseless $\kappa_{{\rm \ QSO \ and \ Ly}\alpha {\rm \ True \ Correlations}$')

