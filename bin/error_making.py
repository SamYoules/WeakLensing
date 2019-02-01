## Author: Ben Mawdsley
## Jan 2019
## Create residual errors fits file for forward fitting kappa maps, by comparing 60 input maps with corresponding estimated maps

import numpy as np
from astropy.io import fits
import healpy as hp


NSIDE=256
nmaps=60
all_residuals=np.zeros((hp.nside2npix(NSIDE), nmaps))

for o in range(1, nmaps):
    maptrue=hp.fitsfunc.read_map('kappa_input'+str(o)+'.fits')
    maptrue_deg=hp.pixelfunc.ud_grade(maptrue, nside_out=NSIDE)
    maprecon=hp.fitsfunc.read_map('kappa-noiseless'+str(o)+'.fits')
    print(hp.npix2nside(len(maprecon)))
    all_residuals[:, o]=np.absolute(maprecon-maptrue_deg)
    print(o)

errormap=np.mean(all_residuals, axis=1)
errormap[np.isnan(errormap)]=0.0
errormap[np.isinf(errormap)]=0.0
print(errormap.shape)
print(hp.npix2nside(len(errormap)))
hp.fitsfunc.write_map('ResidualEstimatedErrors.fits', errormap)

