# Author: Sam Youles
# Modified by Julian Bautista

import numpy as np
import healpy as hp
import sys
import glob
import os
import fitsio
from kappa_lya import *


# input directory name containing delta files
indir = sys.argv[1]
outdir = sys.argv[2]
mapdir = sys.argv[3]
mapnumber = sys.argv[4]
mapname = '{}/kappa_input{}.fits'.format(mapdir, mapnumber)

#-- Create angular power spectrum of kappa
theory = Theory()
ell, cell = theory.get_cl_kappa(2.1, kmax=100., nz=100, lmax=10000)

#nside = 256      # SY 27/11/18 Smooth maps to match nside of estimated maps (Was 1024, this reduces noise)
#npix=nside**2*12
#seed=int(mapnumber)
#np.random.seed(seed)
#kappa = create_gaussian_kappa(ell, cell, nside=nside, seed=seed)
#hp.fitsfunc.write_map(mapname, kappa.A, fits_IDL=False)

#-- Use pre-existing kappa input file
kappa = fitsio.FITS(mapname)

# Amend DEC and RA in each of the delta files by the bend angle from alpha map
alldeltas = glob.glob(indir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
for filename in alldeltas:
    #hdus = fits.open(filename)
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    out = fitsio.FITS(outdir+"/"+os.path.basename(filename),'rw',clobber=True)

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ra = header['RA']
        dec = header['DEC']

        # Add bend angles to ra and dec
        theta_lens, phi_lens = kappa.displace_objects(np.pi/2-dec, ra) 

        theta = np.pi/2-dec
        phi = ra
        ipix = hp.ang2pix(nside, theta, phi)
        dtheta = kappa[ipix]
        dphi = self.dphi_map.A[ipix]/np.sin(theta)
        dd = np.sqrt(dtheta**2+dphi**2)
        #alpha = np.arctan2(dphi,dtheta) + np.pi
        alpha = np.arctan2(dphi,dtheta)  ## SY 19/6/18
        #-- Equation A15 from 0502469
        thetap = np.arccos(np.cos(dd)*np.cos(theta) - 
                           np.sin(dd)*np.sin(theta)*np.cos(alpha))
        phip = phi+np.arcsin(np.sin(alpha)*np.sin(dd)/np.sin(thetap))
        return thetap, phip




        # Rewrite new delta file with new values
        header['RA'] = phi_lens
        header['DEC'] = np.pi/2-theta_lens
        header['RA0'] = ra
        header['DEC0'] = dec
      
        #-- Re-create columns (maybe there's a better way to do this?) 
        ll = hdu['LOGLAM'][:]
        de = hdu['DELTA'][:]
        we = hdu['WEIGHT'][:]
        co = hdu['CONT'][:] 
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols, names=names, header=header, \
                  extname=str(header['THING_ID']))

    out.close()
