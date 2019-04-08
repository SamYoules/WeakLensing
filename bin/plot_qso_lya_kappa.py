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

#wkqsolya = fits.open('kappa-xnoiseless.fits.gz')[1].data.wkappa
#kqsolya = fits.open('kappa-xnoiseless.fits.gz')[1].data.kappa
#klyalya = fits.open('kappa-noiseless-70-100.fits.gz')[1].data.kappa
#kinput = fits.open('kappa_input_noiseless.fits')[1].data.I

#wkqsolya = fits.open('kappa-xnoiseless-true.fits.gz')[1].data.wkappa
#kqsolya = fits.open('kappa-xnoiseless-true.fits.gz')[1].data.kappa
#klyalya = fits.open('kappa-noiseless-70-100.fits.gz')[1].data.kappa
#kinput = fits.open('kappa_input_noiseless.fits')[1].data.I

#---------------

#wkqsolya = fits.open('kappa-xnoisy.fits.gz')[1].data.wkappa
#kqsolya = fits.open('kappa-xnoisy.fits.gz')[1].data.kappa
#klyalya = fits.open('est_maps_xnoisy/kappa1.fits.gz')[1].data.kappa
#kinput = fits.open('est_maps_xnoisy/kappa_input1.fits')[1].data.I

wkqsolya = fits.open('kappa-xnoisy-true.fits.gz')[1].data.wkappa
kqsolya = fits.open('kappa-xnoisy-true.fits.gz')[1].data.kappa
klyalya = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.kappa
#klyalya = fits.open('est_maps_noisy/kappa1-true.fits.gz')[1].data.kappa
kinput = fits.open('est_maps_xnoisy/kappa_input1.fits')[1].data.I

#---------------

#wkqsolya = fits.open('kappa-xcut.fits.gz')[1].data.wkappa
#kqsolya = fits.open('kappa-xcut.fits.gz')[1].data.kappa
#klyalya = fits.open('est_maps_cut_noiseless/kappa1.fits.gz')[1].data.kappa
#kinput = fits.open('est_maps_cut_noiseless/kappa_input1.fits')[1].data.I

#wkqsolya = fits.open('kappa-xcut-true.fits.gz')[1].data.wkappa
#kqsolya = fits.open('kappa-xcut-true.fits.gz')[1].data.kappa
#klyalya = fits.open('est_maps_cut_noiseless/kappa1.fits.gz')[1].data.kappa
#kinput = fits.open('est_maps_cut_noiseless/kappa_input1.fits')[1].data.I

#---------------

# Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(kqsolya.size))
alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Mask off areas outside survey footprint

mask = wkqsolya==0
#mask_size = 1 - wkqsolya[mask].size/wkqsolya.size
mask &= (kqsolya>np.percentile(kqsolya[mask], 0.5)) & \
           (kqsolya<np.percentile(kqsolya[mask], 99.5))
kqsolya[mask]=hp.UNSEEN
klyalya[mask]=hp.UNSEEN
kinput_masked = kinput*(~mask)+hp.UNSEEN*(mask) 

# Get the C_ells for auto- and cross-correlations
Cls1=hp.sphtfunc.anafast(kqsolya, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(klyalya, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(kinput_masked, lmax=3*NSIDE-1)

Cls13=hp.sphtfunc.anafast(kqsolya, kinput_masked, lmax=3*NSIDE-1)
Cls23=hp.sphtfunc.anafast(klyalya, kinput_masked, lmax=3*NSIDE-1)

P.rcParams.update({'font.size':20})
# kappa plot
P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)

#P.plot(Cls1, label=r'$\langle \kappa_{q} \kappa_{q} \rangle$')
#P.plot(Cls2, label=r'$\langle \kappa_{\alpha} \kappa_{\alpha} \rangle$')
P.plot(Cls3, label=r'$\langle \kappa_{input} \ \kappa_{input} \rangle$')
P.plot(Cls23, label=r'$\langle \kappa_{input} \ \kappa_{\alpha \alpha} \rangle$')
P.plot(Cls13, label=r'$\langle \kappa_{input} \ \kappa_{\alpha q} \rangle$')


P.ylabel('$C_l^{\kappa \kappa}$', fontsize=24)
P.xlabel('$l$', fontsize=24)
ax.set_xlim([-2, 800])
#ax.set_ylim([-4e-6, 1.4e-6])
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
P.legend(fontsize=22, loc = "upper right")

#P.title(r'DESI Noiseless $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ Correlations}$')
#P.title(r'DESI Noisy $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ Correlations}$')
#P.title(r'DESI Noiseless Cut $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ Correlations}$')

#P.title(r'DESI Noiseless $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ True \ Correlations}$')
P.title(r'DESI Noisy $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ True \ Correlations}$')
#P.title(r'DESI Noiseless Cut $\kappa {\rm \ QSO \ and \ Ly}\alpha {\rm \ True \ Correlations}$')

