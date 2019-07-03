# Author: Sam Youles
# 31/5/19

import numpy as np
import healpy as hp
import sys
import glob
import os
import fitsio
from kappa_lya import *
import argparse


#-- Input arguments

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute the convergence (kappa) between qsos and deltas.')

    parser.add_argument('--mapdir', required=True, type=str, \
           help='folder containing input maps')
    parser.add_argument('--mapnumber', required=False, type=int, default=1, \
           help='index number of input map')
    args, unknown = parser.parse_known_args()

    mapdir = args.mapdir
    mapnumber = args.mapnumber
    mapname = '{}/kappa_input{}.fits'.format(mapdir, mapnumber)

    #-- Create angular power spectrum of kappa
    theory = Theory()
    ell, cell = theory.get_cl_kappa(2.1, kmax=100., nz=100, lmax=10000)

    #-- Create kappa map
    nside=1024
    npix=nside**2*12
    seed=int(mapnumber)
    np.random.seed(seed)
    kappa = create_gaussian_kappa(ell, cell, nside=nside, seed=seed)
    hp.fitsfunc.write_map(mapname, kappa.A, column_names='I', fits_IDL=False, overwrite=True)

