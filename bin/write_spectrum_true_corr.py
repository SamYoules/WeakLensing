## SY 16/11/18
## Write power spectra for auto-cf, cross-cf and input-auto-cfs for 100 realisations of estimated kappa,
## to be used for plotting errors.

from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
#import kappa_lya
#from kappa_lya import *
import sys

def get_correlations(esttype, maptype, suffix):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    ##- Open maps
    kest = fits.open('maps/{}/true_corr/kappa_{}_{}.fits.gz'.format \
                                  (esttype, maptype, suffix))[1].data.kappa
    wkest = fits.open('maps/{}/true_corr/kappa_{}_{}.fits.gz'.format \
                                 (esttype, maptype, suffix))[1].data.wkappa
    kinput = fits.open('maps/input/kappa_input1.fits')[1].data.I
    kest *= (-1)

    ##- Get resolution of map
    nside=int(hp.npix2nside(kest.size))

    ##- Reset resolution of input maps to match estimated maps
    kinput = hp.ud_grade(kinput, nside)

    ##- Mask area outside footprint
    mask = wkest!=0
    kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
    kinput[~mask]=hp.UNSEEN
    mask &= (kest>np.percentile(kest[mask], 0.5)) & \
                (kest<np.percentile(kest[mask], 99.5))
    kest[~mask]=hp.UNSEEN

    ##- Get the power spectra from the maps
    Cl_auto = hp.sphtfunc.anafast(kest, lmax=3*nside-1)
    Cl_cross = hp.sphtfunc.anafast(kest, kinput, lmax=3*nside-1)

    #Cl_input = hp.sphtfunc.anafast(kinput, lmax=3*nside-1)
    #return Cl_auto, Cl_cross, Cl_input
    return Cl_auto, Cl_cross

##- Input (e.g.  midpoint xnoisy rt70  )
esttype = sys.argv[1]
maptype = sys.argv[2]
suffix  = sys.argv[3]

Cl_autos, Cl_crosses =get_correlations(esttype, maptype, suffix)

np.savetxt('maps/{}/true_corr/Cls/Cl_autos_{}_{}.txt'.format(esttype, maptype, suffix), Cl_autos)
np.savetxt('maps/{}/true_corr/Cls/Cl_crosses_{}_{}.txt'.format(esttype, maptype, suffix), Cl_crosses)

