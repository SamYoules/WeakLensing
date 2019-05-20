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

def get_correlations(filenumber, maptype):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    ##- Open maps
    data = np.load('est_opt_{}/kappa{}.npz'.format(maptype, filenumber))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas*(-1)

    kinput = fits.open('est_maps_{}/kappa_input{}.fits'.format(maptype, \
                                                   filenumber))[1].data.I

    ##- Reset resolution of input map to match estimated map
    kappa_input  = hp.ud_grade(kinput, nside_opt)

    ##- Mask area outside footprint
    wkest = fits.open('est_maps_{}/kappa{}.fits.gz'.format \
                                     (maptype, filenumber))[1].data.wkappa
    mask = wkest!=0
    kappa_input = kappa_input*(mask)+hp.UNSEEN*(~mask) 
    kappa_input[~mask]=hp.UNSEEN

    ##- Get the power spectra from the maps
    cl_input     = hp.sphtfunc.anafast(kappa_input, lmax=3*nside_opt-1)
    cl_auto      = hp.sphtfunc.anafast(kappa_opt, lmax=3*nside_opt-1)
    cl_cross     = hp.sphtfunc.anafast(kappa_opt, kappa_input, lmax=3*nside_opt-1)
    return cl_auto, cl_cross, cl_input

##- Get map type name (e.g. xnoisy)
maptype = sys.argv[1]

c_ells = []
n_files = 100

for i in range(n_files):
    j = str(i+1)
    c_ells.append(get_correlations(j, maptype))

cl_array = np.asarray(c_ells)
cl_autos = cl_array[:,0]
cl_crosses = cl_array[:,1]
cl_inputs = cl_array[:,2]
np.savetxt('est_opt_{}/Cls/Cl_autos_{}.txt'.format(maptype,maptype), cl_autos)
np.savetxt('est_opt_{}/Cls/Cl_crosses_{}.txt'.format(maptype,maptype), cl_crosses)
np.savetxt('est_opt_{}/Cls/Cl_inputs.txt'.format(maptype,maptype), cl_inputs)



