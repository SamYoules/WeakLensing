import kappa_lya
from astropy.io import fits
import numpy as np
import pylab as plt
import healpy as hp
import os

def read_kappa_input(nside=256):

    kappa_input = fits.open('kappa_input.fits')[1].data.kappa.ravel()
    alm_input = hp.map2alm(kappa_input, lmax=3*nside-1)
    kappa_input = hp.alm2map(alm_input, nside=nside)
    return kappa_input

def get_cl_cross(kappa_input, kappa, mask=None):
    kappa_input_masked = kappa_input*1.
    if mask is not None:
        kappa_input_masked[~mask] = hp.UNSEEN
    cl_cross = hp.anafast(kappa_input_masked, kappa)
    return cl_cross

def plot_cl_vs_r(rr='rt', rp=10, rt=70, nside=256, lw=0.5):
   
    if rr=='rt':
        root = 'kappa-noiseless-cnstwe-fromdeltas-nside256-rt%d-'+'rp%d.fits.gz' %rp
    else:
        root = 'kappa-noiseless-cnstwe-fromdeltas-nside256-rt%d-'%rt+'rp%d.fits.gz'

    kappa_input = read_kappa_input(nside=nside)
    ell = np.arange(3*nside)
    f, ax = plt.subplots(1, 3, figsize=(14, 5))
    c = 0
    for i, r in enumerate(range(10, 170, 10)):
        if not os.path.exists(root%r): continue
        data = fits.open(root%r)[1].data
        kappa = data.kappa
        mask = data.wkappa!=0
        #mask &= (data.kappa>np.percentile(data.kappa[mask], 0.5)) & \
        #        (data.kappa<np.percentile(data.kappa[mask], 99.5))
        mask &= data.npairs>np.percentile(data.npairs[mask], 0.5) 
        kappa[~mask]=hp.UNSEEN
        f_sky = sum(mask)/mask.size
        kappa_input_masked = kappa_input*(mask)+hp.UNSEEN*(~mask) 
        cl_data = hp.anafast(kappa)
        cl_input = hp.anafast(kappa_input_masked)
        cl_cross = get_cl_cross(kappa, kappa_input_masked)

        ax[0].plot(ell, cl_data, lw=lw,  label=rr+'max=%d'%r)
        ax[1].plot(ell, cl_input, lw=lw, label=rr+'max=%d'%r)
        ax[2].plot(ell, cl_data/cl_input, lw=lw, label=rr+'max=%d'%r)
        #ax[3].plot(ell, cl_cross, lw=lw, label=rr+'=%d'%r)
        c+=1
    #ax[3].axhline(0., color='k', ls=':')
    for a in ax:
        a.set_yscale('log')
        a.legend(loc=0, fontsize=8)
    
    ax[0].set_ylabel(r'Estimated $C_{\ell}$')
    ax[1].set_ylabel(r'Input (masked) $C_{\ell}$')
    ax[2].set_ylabel('Ratio of estimated/input')
    plt.tight_layout()

