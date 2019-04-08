from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np

nside = 64
kinput = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
kappa_input = hp.ud_grade(kinput, nside)
cl_input = hp.sphtfunc.anafast(kappa_input)
np.savetxt('kappa_opt_srad/cls/cl_input_64.txt', cl_input)

