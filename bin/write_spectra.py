## SY 16/11/18
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

def get_correlations(filenumber, maptype, NSIDE):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    ##- Open maps
    kest = fits.open('est_maps_{}/kappa{}.fits.gz'.format \
                                                   (maptype, filenumber))[1].data.kappa
    wkest = fits.open('est_maps_{}/kappa{}.fits.gz'.format \
                                                   (maptype, filenumber))[1].data.wkappa
    kinput = fits.open('est_maps_{}/kappa_input{}.fits'.format \
                                                   (maptype, filenumber))[1].data.I
    kest *= (-1)
    ##- Reset resolution of input maps to match estimated maps
    alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
    kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

    ##- Mask area outside footprint
    mask = wkest!=0
    kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
    kinput[~mask]=hp.UNSEEN
    mask &= (kest>np.percentile(kest[mask], 0.5)) & \
                (kest<np.percentile(kest[mask], 99.5))
    kest[~mask]=hp.UNSEEN

    ##- Get the power spectra from the maps
    Cl_auto = hp.sphtfunc.anafast(kest, lmax=3*NSIDE-1)
    Cl_cross = hp.sphtfunc.anafast(kest, kinput, lmax=3*NSIDE-1)
    Cl_input = hp.sphtfunc.anafast(kinput, lmax=3*NSIDE-1)
    return Cl_auto, Cl_cross, Cl_input

##- Get map type name (e.g. xnoisy)
maptype = sys.argv[1]

##- Open a kappa file to get nside and create mask
kxi = fits.open('est_maps_{}/kappa1.fits.gz'.format(maptype))[1].data.kappa
wkxi = fits.open('est_maps_{}/kappa1.fits.gz'.format(maptype))[1].data.wkappa

##- Get resolution of map
NSIDE=int(hp.npix2nside(kxi.size))

##- Mask off areas outside DESI footprint
#mask = wkxi!=0
#mask &= (kxi>np.percentile(kxi[mask], 0.5)) & \
                #(kxi<np.percentile(kxi[mask], 99.5))

##- Get the C_ells for each estimated kappa realisation
##- (auto and crossed with input)
C_ells = []
N_files = 100

for i in range(N_files):
    j = str(i+1)
    C_ells.append(get_correlations(j, maptype, NSIDE))

Cl_array = np.asarray(C_ells)
Cl_autos = Cl_array[:,0]
Cl_crosses = Cl_array[:,1]
Cl_inputs = Cl_array[:,2]
np.savetxt('Cls/Cl_autos_{}.txt'.format(maptype), Cl_autos)
np.savetxt('Cls/Cl_crosses_{}.txt'.format(maptype), Cl_crosses)
np.savetxt('Cls/Cl_inputs_{}.txt'.format(maptype), Cl_inputs)



