from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys

nside=64

wkappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside{}.fits.gz'.format(nside))[1].data.wkappa
kappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside{}.fits.gz'.format(nside))[1].data.kappa

##- Mask off areas outside DESI footprint
mask = wkappa_mid!=0
mask &= (kappa_mid>np.percentile(kappa_mid[mask], 0.5)) & \
                (kappa_mid<np.percentile(kappa_mid[mask], 99.5))
kappa_mid[~mask]=hp.UNSEEN

##- Open input map and reset nside
kinput = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
kappa_input = hp.ud_grade(kinput, nside)

##- Get Midpoint Cls
cla_mid = hp.sphtfunc.anafast(kappa_mid)
clx_mid = hp.sphtfunc.anafast(kappa_input, kappa_mid)

np.savetxt('kappa_opt_srad/cls/cl_auto_cut_midpoint_100_{}.txt'.format(nside), cla_mid)
np.savetxt('kappa_opt_srad/cls/cl_cross_cut_midpoint_100_{}.txt'.format(nside), clx_mid)

