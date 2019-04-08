from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys

def load_kappa_maps(nside_opt, sradius, damping, k0):
    '''Reads in optimal estimated kappa maps for different search radii, and gets the Cls
       normalised by pixel area.'''

    data = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside{}_srad{}_{}.npz'.format(nside_opt, sradius, damping))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    cl_auto      = hp.sphtfunc.anafast(kappa_opt)
    cl_cross     = hp.sphtfunc.anafast(kappa_opt, k0)
    return cl_auto, cl_cross

nside = 64
rtmax = 100
srad  = '01'
damping = '10'

##- Open input map and reset nside
kinput = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
kappa_input = hp.ud_grade(kinput, nside)

##- Open optimal map and get cls
cl_auto, cl_cross = load_kappa_maps(nside, srad, damping, kappa_input)
#cl_input = hp.sphtfunc.anafast(kappa_input)

np.savetxt('kappa_opt_srad/cls/cl_auto_cut_{}_{}_{}_{}.txt'.format(nside, rtmax, srad, damping), cl_auto)
np.savetxt('kappa_opt_srad/cls/cl_cross_cut_{}_{}_{}_{}.txt'.format(nside, rtmax, srad, damping), cl_cross)
#np.savetxt('kappa_opt_srad/cls/cl_input.txt', cl_input)

