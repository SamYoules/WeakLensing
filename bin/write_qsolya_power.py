from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys


##- Open a kappa file to get nside and create mask

##- Open maps
kest = fits.open('kappa-cross.fits.gz')[1].data.kappa
wkest = fits.open('kappa-cross.fits.gz')[1].data.wkappa
kinput = fits.open('est_maps_10_100/kappa_input1.fits')[1].data.I
kest *= (-1)

##- Get resolution of map
NSIDE=int(hp.npix2nside(kest.size))

##- Mask off areas outside DESI footprint
mask = wkest!=0
#mask &= (kest>np.percentile(kest[mask], 0.5)) & \
               # (kest<np.percentile(kest[mask], 99.5))

##- Reset resolution of input maps to match estimated maps
alm = hp.sphtfunc.map2alm(kinput, lmax=3*NSIDE-1)
kinput = hp.sphtfunc.alm2map(alm, nside=NSIDE)

##- Mask area outside footprint
kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
kinput[~mask]=hp.UNSEEN
kest[~mask]=hp.UNSEEN

##- Get the power spectra from the maps
Cl_auto = hp.sphtfunc.anafast(kest, lmax=3*NSIDE-1)
Cl_cross = hp.sphtfunc.anafast(kest, kinput, lmax=3*NSIDE-1)
Cl_input = hp.sphtfunc.anafast(kinput, lmax=3*NSIDE-1)

Cl_autos = Cl_auto
Cl_crosses = Cl_cross
Cl_inputs = Cl_input
np.savetxt('Cl_autos_qsolya.txt', Cl_autos)
np.savetxt('Cl_crosses_qsolya.txt', Cl_crosses)
np.savetxt('Cl_inputs_qsolya.txt', Cl_inputs)

