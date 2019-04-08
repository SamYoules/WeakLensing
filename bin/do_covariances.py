## SY 31/10/18
## Reads in power spectra of 100 realisations of estimated kappa for LYAxLYA and QSOxLYA and input kappa.
## Plots cross-cf error bars (for noisy and noiseless datasets) for input and estimated kappa.

from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys
from collections import OrderedDict
from astropy.table import Table

def rebin(Cls,ellmin=40, ellmax=768, nell=7):
    '''Smooth curves by rebinning Cls'''
    ell = np.arange(Cls.size)
    weights = 2.*ell+1.
    if not ellmin:
        ellmin = min(ell)
    if not ellmax:
        ellmax = max(ell)
    if nell==0:
        nell = ellmax-ellmin

    w = (ell>=ellmin)&(ell<=ellmax)
    index = np.floor( (ell[w]-ellmin)*1./(ellmax-ellmin)*nell ).astype(int)
    well = np.bincount( index, weights=weights[w])
    sell = np.bincount( index, weights=weights[w]*ell[w])
    scl  = np.bincount( index, weights=weights[w]*Cls[w])
    ell = sell/well
    Cls = scl/well
        
    return ell, Cls, well 


def get_vals(Cl_crosses_in):

    ##- Trim off smallest and largest scales (gives 7 bins of ell=104)
    Cl_crosses = Cl_crosses_in[:,40:768]

    ##- Rebin
    cross_rebin = []
    cross_weights = []

    for i,j in enumerate(Cl_crosses):
        ell, crosses, cross_wei = rebin(Cl_crosses[i])
        cross_rebin.append(crosses)
        cross_weights.append(cross_wei)

    cross_cl = np.asarray(cross_rebin)
    cross_we = np.asarray(cross_weights)
    cross_mean = (np.sum(cross_cl*cross_we,axis=0) / np.sum(cross_we,axis=0))

    ##- Calculate covariance matrices, variance and errors
    QW = (cross_cl-cross_mean)*cross_we

    return QW, cross_we


##- Open kappa true-correlation files
kxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.kappa
wkxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.wkappa
#kxi = fits.open('kappa-noiseless-cnstwe-lensed-truecorr-rtmax70-rpmax100.fits.gz')[1].data.SKAPPA
#wkxi = fits.open('kappa-noiseless-cnstwe-lensed-truecorr-rtmax70-rpmax100.fits.gz')[1].data.WKAPPA
kxi_input = fits.open('est_maps_noiseless/kappa_input1.fits')[1].data.I

##- Get resolution of map
NSIDE=int(hp.npix2nside(kxi.size))

##- Mask off areas outside DESI footprint
mask = wkxi!=0
mask &= (kxi>np.percentile(kxi[mask], 0.5)) & \
                (kxi<np.percentile(kxi[mask], 99.5))
kxi[~mask]=hp.UNSEEN

##- Reset resolution of input maps to match estimated maps
alm = hp.sphtfunc.map2alm(kxi_input, lmax=3*NSIDE-1)
kxi_input = hp.sphtfunc.alm2map(alm, nside=NSIDE)

##- Mask area outside footprint
kxi_input = kxi_input*(mask)+hp.UNSEEN*(~mask) 
kxi_input[~mask]=hp.UNSEEN

##- Get input Cls
Cl_xi_input = hp.sphtfunc.anafast(kxi_input, lmax=3*NSIDE-1)

##- Get File names for auto and cross power spectra
Clx = ['Cl_crosses_noisy.txt', 'Cl_crosses_xnoisy.txt', 'Cl_crosses_noiseless.txt', 'Cl_crosses_xnoiseless.txt', 'Cl_crosses_cut.txt', 'Cl_crosses_xcut.txt']

QW1 = []
cwe = []

##- Read in cross- and input- power spectra
for i,j in enumerate(Clx):
    Cl_crosses_in = np.loadtxt(j)
    QW, cr_weights = get_vals(Cl_crosses_in)
    QW1.append(QW)
    cwe.append(cr_weights)

##- Write 14x14 covariance and correlation matrices for auto and cross
i=[0,2,4]
k=['noisy','noiseless','cut']
fout1=['cov_noisy.txt', 'cov_noiseless.txt', 'cov_cut.txt']
fout2=['corr_noisy.txt', 'corr_noiseless.txt', 'corr_cut.txt']
for l,j in enumerate(i):
    QW2 = np.concatenate([QW1[j], QW1[j+1]], axis=1)
    cwe2 = np.concatenate([cwe[j], cwe[j+1]], axis=1)

    cross_cov = QW2.T.dot(QW2)/cwe2.T.dot(cwe2)
    np.savetxt(fout1[l], cross_cov)

    rho = cross_cov/np.sqrt(np.outer(np.diag(cross_cov),np.diag(cross_cov)))
    np.savetxt(fout2[l], rho)

    #cross_var = sp.diagonal(cross_cov)
    #cross_stdev = np.sqrt(cross_var)

