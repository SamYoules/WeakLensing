## SY 8/12/18
## Plots cross-correlations of qso-lya kappa power spectrum.

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys

def rebin(Cls,ellmin=40, ellmax=768, nell=20):
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
wk1 = fits.open('kappa-xcut-lensed-true-rt20-rp100.fits.gz')[1].data.wkappa
k1 = fits.open('kappa-xcut-lensed-true-rt20-rp100.fits.gz')[1].data.kappa
k2 = fits.open('kappa-xcut-lensed-true-rt40-rp100.fits.gz')[1].data.kappa
k3 = fits.open('kappa-xcut-lensed-true-rt70-rp100.fits.gz')[1].data.kappa
k4 = fits.open('kappa-xcut-lensed-true-rt100-rp100.fits.gz')[1].data.kappa
k5 = fits.open('kappa-xcut-lensed-true-rt140-rp100.fits.gz')[1].data.kappa
kinput = fits.open('est_maps_xcut/kappa_input1.fits')[1].data.I


# Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(k1.size))
alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Mask off areas outside survey footprint

mask = wk1==0
#mask_size = 1 - wk1[mask].size/wk1.size
mask &= (k1>np.percentile(k1[mask], 0.5)) & \
           (k1<np.percentile(k1[mask], 99.5))
k1[mask]=hp.UNSEEN
k2[mask]=hp.UNSEEN
k3[mask]=hp.UNSEEN
k4[mask]=hp.UNSEEN
k5[mask]=hp.UNSEEN
kinput_masked = kinput*(~mask)+hp.UNSEEN*(mask) 

# Get the C_ells for auto- and cross-correlations
Cls1=hp.sphtfunc.anafast(k1, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(k2, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(k3, lmax=3*NSIDE-1)
Cls4=hp.sphtfunc.anafast(k4, lmax=3*NSIDE-1)
Cls5=hp.sphtfunc.anafast(k5, lmax=3*NSIDE-1)
Cls6=hp.sphtfunc.anafast(kinput_masked, lmax=3*NSIDE-1)

ell1, Cell1, well = rebin(Cls1)
ell2, Cell2, well = rebin(Cls2)
ell3, Cell3, well = rebin(Cls3)
ell4, Cell4, well = rebin(Cls4)
ell5, Cell5, well = rebin(Cls5)
ell6, Cell6, well = rebin(Cls6)

P.rcParams.update({'font.size':20})
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)

#P.plot(Cls1, label='rtmax=20')
#P.plot(Cls2, label='rtmax=40')
#P.plot(Cls3, label='rtmax=70')
#P.plot(Cls4, label='rtmax=100')
#P.plot(Cls5, label='rtmax=140')
#P.plot(Cls6, label='Input')

P.plot(ell1, Cell1, label='rtmax=20')
P.plot(ell2, Cell2, label='rtmax=40')
P.plot(ell3, Cell3, label='rtmax=70')
P.plot(ell4, Cell4, label='rtmax=100')
P.plot(ell5, Cell5, label='rtmax=140')
P.plot(ell6, Cell6, label='Input')


P.ylabel('$C_l^{\kappa \kappa}$', fontsize=24)
P.xlabel('$l$', fontsize=24)
ax.set_xlim([-2, 800])
#ax.set_ylim([-4e-6, 1.4e-6])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
P.legend(fontsize=22, loc = "upper right")

P.title(r'DESI Noiseless $\langle \kappa_{{\rm \ QSO}} \kappa_{{\rm \ Ly}\alpha} \rangle {\rm \ True \ Correlations}$')

