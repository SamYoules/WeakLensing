## SY 9/10/18
## Plots auto-correlations of kappa power spectra for noisy maps with varying R_p max to see which is optimal. Clips off outliers which were creating odd results.

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
    alm_input = hp.map2alm(kappa_input, lmax=3*nside-1)
    kappa_input = hp.alm2map(alm_input, nside=nside)
    return kappa_input

# Open kappa fits files
root = 'kappa-noiseless-70-%d.fits.gz'
#root = 'kappa-noiseless-40-%d.fits.gz'
nside=256
lw=0.5
kappa_input = read_kappa_input(nside=nside)
ell = np.arange(3*nside)

P.ion()
fig = P.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1)
ncolors=11
colors = plt.cm.jet(np.linspace(0,1,ncolors))
c = 0

for i, r in enumerate([10,20,30,40,50,60,70,80,90,100]):
    if not os.path.exists(root%r): continue
    data = fits.open(root%r)[1].data
    kappa = data.kappa
    mask = data.wkappa!=0
    mask &= (data.kappa>np.percentile(data.kappa[mask], 0.5)) & \
                (data.kappa<np.percentile(data.kappa[mask], 99.5))
    #mask &= data.npairs>np.percentile(data.npairs[mask], 0.5) 
    kappa[~mask]=hp.UNSEEN
    f_sky = sum(mask)/mask.size
    kappa_input_masked = kappa_input*(mask)+hp.UNSEEN*(~mask) 
    cl_data = hp.anafast(kappa)
    cl_input = hp.anafast(kappa_input_masked)
    P.plot(ell, cl_data, lw=lw,  color=colors[c], label='rp max=%d'%r)
    c+=1
    #if c==ncolors:
    #    break

ax.set_yscale('log')
P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
P.xlabel('$l$', fontsize=18)
P.plot(ell, cl_input, lw=lw, color=colors[c], label='input')
P.legend()
#leg = P.legend(fontsize=18, bbox_to_anchor=(0.95,0.45), loc = "center right")
#for legobj in leg.legendHandles:
#    legobj.set_linewidth(2.0)
P.title(r'${\rm Noiseless} \ \kappa {\rm \ Auto-Correlations, rt = 70}$')


