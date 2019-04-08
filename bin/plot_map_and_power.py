import numpy as np
import pylab as plt
import healpy as hp
from astropy.io import fits
import sys

def rebin(Cls,ellmin=2, ellmax=400, nell=200):
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

mapin = sys.argv[1]
# kappalists/kappa_cut.npz
# kappalists/kappa_noisy.npz
# kappalists/map.npz
maptype = sys.argv[2]
# Noiseless
# Noisy
# Anze

plt.ion()
## Load estimated map
da=np.load(mapin); Nside=int(da['arr_0']); ids=da['arr_1']; mp=da['arr_2'][0];
pix=np.zeros(Nside**2*12)
pix[ids]=mp
hp.mollview(pix,rot=(180,0,0),min=-1e-4,max=1e-4)
#hp.mollview(pix,rot=(180,0,0),min=-0.05,max=0.05)
#hp.mollview(pix,rot=(180,0,0),min=-0.01,max=0.01)

#k_input = fits.open('est_maps_noiseless/kappa_input1.fits')[1].data.I
k_input = fits.open('kappalists/kappa_input.fits')[1].data.kappa
k_input = np.ravel(k_input)
k_input = hp.ud_grade(k_input, Nside)

## Set Nside to match estimated map
#alm = hp.sphtfunc.map2alm(k_input, lmax=3*Nside-1)
#k_input = hp.sphtfunc.alm2map(alm, nside=Nside)

## Get Cls
Cl_auto = hp.sphtfunc.anafast(pix, lmax=3*Nside-1)
Cl_cross = hp.sphtfunc.anafast(pix, k_input, lmax=3*Nside-1)

##- Rebin
#ell_a, Cl_a, w_a = rebin(Cl_auto)
#ell_c, Cl_c, w_c = rebin(Cl_cross)

ell = np.arange(Cl_auto.size)

fig=plt.figure()
#plt.plot(ell_c,Cl_c, lw=2, color='b', label='{} x input'.format(maptype))
#plt.plot(ell_a,Cl_a, lw=2, color='r', label='{} auto'.format(maptype))
plt.plot(ell,Cl_cross, color='b', label='{} x input'.format(maptype))
plt.plot(ell,Cl_auto, color='r', label='{} auto'.format(maptype))
plt.xlabel(r'$\ell$', fontsize=18)
plt.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
plt.ylim([-1e-9,1e-9])
#plt.ylim([-1e-7,1e-7])
plt.legend()
