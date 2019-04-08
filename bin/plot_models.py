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
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText

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


def cov(cl, we):
    '''Computes weighted covariance matrix'''
    mcl = (cl*we).sum(axis=0)
    swe = we.sum(axis=0)
    w = swe>0.
    mcl[w] /= swe[w]
    wcl = we*(cl-mcl)
    print("Computing cov...")
    co = wcl.T.dot(wcl)
    sswe = swe*swe[:,None]
    w = sswe>0.
    co[w] /= sswe[w]
    return co

def get_vals(Cl_crosses_in, Cl_inputs):

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

    input_mean = (np.sum(Cl_inputs, axis=0) / 100)

    ##- Calculate covariance matrices, variance and errors
    QW = (cross_cl-cross_mean)*cross_we
    cross_cov = QW.T.dot(QW)/cross_we.T.dot(cross_we)
    cross_var = sp.diagonal(cross_cov)
    cross_stdev = np.sqrt(cross_var)
    #cross_stdev = np.sqrt(cross_var/100)

    ##- Calc chi2
    chi2 = np.sum((cross_mean - 0.)**2/cross_stdev**2 )

    return ell, cross_mean, cross_stdev, input_mean, chi2

##- Open xkappa true-correlation files
kxicut = fits.open('kappa_xcut_true.fits.gz')[1].data.kappa
wkxicut = fits.open('kappa_xcut_true.fits.gz')[1].data.wkappa
kxinoisy = fits.open('kappa_xnoisy_true.fits.gz')[1].data.kappa
wkxinoisy = fits.open('kappa_xnoisy_true.fits.gz')[1].data.wkappa
kxinoiseless = fits.open('kappa_xnoiseless_true.fits.gz')[1].data.kappa
wkxinoiseless = fits.open('kappa_xnoiseless_true.fits.gz')[1].data.wkappa
kxi_input = fits.open('est_maps_noiseless/kappa_input1.fits')[1].data.I

##- Get resolution of map
NSIDE=int(hp.npix2nside(kxicut.size))

##- Mask off areas outside DESI footprint
mask = wkxicut!=0
mask &= (kxicut>np.percentile(kxicut[mask], 0.5)) & \
                (kxicut<np.percentile(kxicut[mask], 99.5))
kxicut[~mask]=hp.UNSEEN
kxinoisy[~mask]=hp.UNSEEN
kxinoiseless[~mask]=hp.UNSEEN

##- Reset resolution of input maps to match estimated maps
alm = hp.sphtfunc.map2alm(kxi_input, lmax=3*NSIDE-1)
kxi_input = hp.sphtfunc.alm2map(alm, nside=NSIDE)

##- Mask area outside footprint
kxi_input = kxi_input*(mask)+hp.UNSEEN*(~mask) 
kxi_input[~mask]=hp.UNSEEN

##- Get true_corr Cls
Cl_xi_input = hp.sphtfunc.anafast(kxi_input, lmax=3*NSIDE-1)
Cl_xicut = hp.sphtfunc.anafast(kxicut, lmax=3*NSIDE-1)
Cl_xinoisy = hp.sphtfunc.anafast(kxinoisy, lmax=3*NSIDE-1)
Cl_xinoiseless = hp.sphtfunc.anafast(kxinoiseless, lmax=3*NSIDE-1)

##- Get File names

Clx = ['Cl_crosses_noisy.txt', 'Cl_crosses_xnoisy.txt', 'Cl_crosses_cnoisy.txt', 'Cl_crosses_noiseless.txt', 'Cl_crosses_xnoiseless.txt', 'Cl_crosses_cnoiseless.txt', 'Cl_crosses_cut.txt', 'Cl_crosses_xcut.txt', 'Cl_crosses_ccut.txt']

Cli = ['Cl_inputs_noisy.txt', 'Cl_inputs_xnoisy.txt', 'Cl_inputs_cnoisy.txt', 'Cl_inputs_noiseless.txt', 'Cl_inputs_xnoiseless.txt', 'Cl_inputs_cnoiseless.txt', 'Cl_inputs_cut.txt', 'Cl_inputs_xcut.txt', 'Cl_inputs_ccut.txt']

cmean = []
cstdev = []
imean = []
chi_2 = []

##- Read in cross- and input- power spectra
for i,j in enumerate(Clx):
    Cl_crosses_in = np.loadtxt(j)
    Cl_inputs = np.loadtxt(Cli[i])
    ell, cr_mean, cr_stdev, in_mean, chi = get_vals(Cl_crosses_in, Cl_inputs)
    cmean.append(cr_mean)
    cstdev.append(cr_stdev)
    imean.append(in_mean)
    chi_2.append(chi)

chi2 = np.asarray(chi_2)
cross_stdev = np.asarray(cstdev)
cross_mean = np.asarray(cmean)
input_mean = np.asarray(imean)


##- fit for large-angles
ycut = (Cl_crosses_in[7] / Cl_inputs[7])
xcut = np.arange(ycut.size)
z = np.polyfit(xcut, ycut, 5)
modelcut = np.polyval(z, xcut)

ynoisy = (Cl_crosses_in[1] / Cl_inputs[1])
xnoisy = np.arange(ynoisy.size)
z = np.polyfit(xnoisy, ynoisy, 5)
modelnoisy = np.polyval(z, xnoisy)

ynoiseless = (Cl_crosses_in[4] / Cl_inputs[4])
xnoiseless = np.arange(ynoiseless.size)
z = np.polyfit(xnoiseless, ynoiseless, 5)
modelnoiseless = np.polyval(z, xnoiseless)

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

P.figure()
P.plot(input_mean[1]*modelnoisy*xnoisy, color=colors[0], lw=2, linestyle="-",label='Noisy')
P.plot(input_mean[7]*modelcut*xcut, color=colors[1], lw=2, linestyle="-",label='Noiseless')
P.plot(input_mean[4]*modelnoiseless*xnoiseless, color=colors[2], lw=2, linestyle="-",label='High Density')
P.legend()
