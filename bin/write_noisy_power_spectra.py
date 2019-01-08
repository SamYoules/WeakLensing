## SY 15/11/18
## Write power spectra for auto-cf, cross-cf and input-auto-cfs for 100 realisations of estimated kappa,
## to be used for plotting errors.

from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys

def get_correlations(filenumber, mask, NSIDE):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    ##- Open maps
    knoisy = fits.open('est_maps/kappa{}.fits.gz'.format \
                                                   (filenumber))[1].data.kappa
    kinput = fits.open('est_maps/kappa_input{}.fits'.format \
                                                   (filenumber))[1].data.I

    knoisy *= (-1)
    ##- Reset resolution of input map to match estimated maps
    alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
    kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

    ##- Mask area outside footprint
    kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
    kinput[~mask]=hp.UNSEEN
    knoisy[~mask]=hp.UNSEEN

    ##- Get the power spectra from the maps
    Cl_auto = hp.sphtfunc.anafast(knoisy, lmax=3*NSIDE-1)
    Cl_cross = hp.sphtfunc.anafast(knoisy, kinput, lmax=3*NSIDE-1)
    Cl_input = hp.sphtfunc.anafast(kinput, lmax=3*NSIDE-1)
    return Cl_auto, Cl_cross, Cl_input

##- Open kappa true-correlation files
kxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.kappa
wkxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.wkappa

##- Get resolution of map
NSIDE=int(hp.npix2nside(kxi.size))

##- Mask off areas outside DESI footprint
mask = wkxi!=0
mask &= (kxi>np.percentile(kxi[mask], 0.5)) & \
                (kxi<np.percentile(kxi[mask], 99.5))

##- Get the C_ells for each estimated kappa realisation
##- (auto and crossed with input)
C_ells = []
N_files = 100

for i in range(N_files):
    j = str(i+1)
    C_ells.append(get_correlations(j, mask, NSIDE))

Cl_array = np.asarray(C_ells)
Cl_autos = Cl_array[:,0]
Cl_crosses = Cl_array[:,1]
Cl_inputs = Cl_array[:,2]
np.savetxt('Cl_autos_noisy.txt', Cl_autos)
np.savetxt('Cl_crosses_noisy.txt', Cl_crosses)
np.savetxt('Cl_inputs_noisy.txt', Cl_inputs)



