## Author: Ben Mawdsley
## Modified by Sam Youles
## Jan 2019
## Calls fitting_kappa.py
## Forward fits a kappa map. Reads in estimated kappa map to get nside and footprint. Reads in error fits file (created in error_making.py). Generates random kappa map of correct resolution and footprint. Calls fitting_kappa.py, which uses perturbation theory to find map with best chi2. Outputs fitted map.
## Command line: python testing_fitting.py map-in error-map map-out-suffix

import numpy as np
import healpy as hp
import fitting_kappa
import sys

#map_data=hp.fitsfunc.read_map('est_maps_noiseless/kappa-noiseless1.fits.gz'
#map_errors=hp.fitsfunc.read_map('ResidualEstimatedErrors.fits')
map_data=hp.fitsfunc.read_map(sys.argv[1])
map_errors=hp.fitsfunc.read_map(sys.argv[2])

mask=np.zeros(map_data.shape)
mask[np.where(map_data!=0)]=1.0

map_errors=mask*map_errors


map_errors[np.where(map_errors<0)]=0.0
map_errors[np.isinf(map_errors)]=0.0
map_errors[np.isnan(map_errors)]=0.0
print(map_errors)
print(map_errors[np.where(map_errors==0)])

#NSIDE=256
NSIDE=64 ## SY: lower resolution for blob kappa

randommap=np.random.normal(0.0, 1.0, hp.nside2npix(NSIDE))

print('random map go', randommap)
random_alm=hp.sphtfunc.map2alm(randommap, lmax=2*NSIDE-1)


fitting=fitting_kappa.fit_map(random_alm, map_data, map_errors)

np.savetxt('OutputElm'+str(sys.argv[3])+'.txt', fitting)

outmap1 = fitting_kappa.ealm2kappamap(fitting, NSIDE)
outmap = np.reshape(outmap1, 12*NSIDE**2)
hp.fitsfunc.write_map('FittedKappaMap'+str(sys.argv[3])+'.fits', outmap)
