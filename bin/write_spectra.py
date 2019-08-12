## SY 16/11/18
## Write power spectra for auto-cf, cross-cf and input-auto-cfs for 100 realisations of estimated kappa,
## to be used for plotting errors.

from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import sys

def get_correlations(filenumber, esttype, maptype):
    '''Open the estimated and input kappa files and get the auto and cross
       power spectrum'''

    print (filenumber, "of 100")
    ##- Open maps
    kest = fits.open('maps/{}/{}/kappa{}.fits.gz'.format \
                                  (esttype, maptype, filenumber))[1].data.kappa
    wkest = fits.open('maps/{}/{}/kappa{}.fits.gz'.format \
                                 (esttype, maptype, filenumber))[1].data.wkappa
    kinput = fits.open('maps/input/kappa_input{}.fits'.format \
                                                        (filenumber))[1].data.I
    #kest *= (-1)

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

##- Input (e.g.  midpoint xnoisy)
esttype = sys.argv[1]
maptype = sys.argv[2]

C_ells = []
N_files = 100

for i in range(N_files):
    j = str(i+1)
    C_ells.append(get_correlations(j, esttype, maptype))

Cl_array = np.asarray(C_ells)
Cl_autos = Cl_array[:,0]
Cl_crosses = Cl_array[:,1]
#Cl_inputs = Cl_array[:,2]
np.savetxt('maps/{}/{}/Cls/Cl_autos.txt'.format(esttype, maptype), Cl_autos)
np.savetxt('maps/{}/{}/Cls/Cl_crosses.txt'.format(esttype, maptype), Cl_crosses)
#np.savetxt('maps/input/Cl_inputs.txt', Cl_inputs)



