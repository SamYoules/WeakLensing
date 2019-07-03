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


def rebin(cls,ellmin=40, ellmax=768, nell=20):
    '''Smooth curves by rebinning cls'''
    ell = np.arange(cls.size)
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
    scl  = np.bincount( index, weights=weights[w]*cls[w])
    ell = sell/well
    cls = scl/well
        
    return cls


##- Open kappa true-correlation files and corresponding input map
kinput = fits.open('maps/input/kappa_input1.fits')[1].data.I
k_maps = ['noisy', 'xnoisy', 'cut', 'xcut', 'noiseless', 'xnoiseless']
model = []
ell = []
cl_in_masked = []

for i in k_maps:
    kxi = fits.open('maps/midpoint/true_corr/kappa_{}_rt70.fits.gz'.format(i)) \
                                                                 [1].data.kappa
    wkxi = fits.open('maps/midpoint/true_corr/kappa_{}_rt70.fits.gz'.format(i)) \
                                                                [1].data.wkappa
    nside=int(hp.npix2nside(kxi.size))
    kin = hp.ud_grade(kinput, nside)
    mask = wkxi!=0
    mask &= (kxi>np.percentile(kxi[mask], 0.5)) & \
                                           (kxi<np.percentile(kxi[mask], 99.5))
    kxi[~mask]=hp.UNSEEN
    kin = kin*(mask)+hp.UNSEEN*(~mask) 
    kin[~mask]=hp.UNSEEN
    kin = hp.ud_grade(kin, nside)
    cl_in = hp.sphtfunc.anafast(kin, lmax=3*nside-1)
    cl_xi = hp.sphtfunc.anafast(kxi, lmax=3*nside-1)
    y = (cl_xi / cl_in)
    x = np.arange(y.size)
    z = np.polyfit(x, y, 5)
    mod = np.polyval(z, x)
    model.append(mod)
    ell.append(x)
    cl_in_masked.append(cl_in)

##- Get input cls for creating model
cl_inputs = np.loadtxt('maps/input/Cl_inputs.txt')
input_mean = (np.sum(cl_inputs, axis=0) / 100)

##- Setup figures
P.rcParams.update({'font.size':18})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
#colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

##- Plot bias figure
P.figure(figsize=(8.2,6))
#P.plot(input_mean[0]*x, color=colors[1], linewidth=2.0, linestyle="-",label='Masked Input')

#P.plot(cl_in_masked[0]*model[0]*ell[0], color=colors[2], lw=2, linestyle="-",label='Noisy Auto')
#P.plot(cl_in_masked[1]*model[1]*ell[1], color=colors[2], lw=2, linestyle="--",label='Noisy Cross')
#P.plot(cl_in_masked[2]*model[2]*ell[2], color=colors[3], lw=2, linestyle="-",label='Noiseless Auto')
#P.plot(cl_in_masked[3]*model[3]*ell[3], color=colors[3], lw=2, linestyle="--",label='Noiseless Cross')
#P.plot(cl_in_masked[4]*model[4]*ell[4], color=colors[4], lw=2, linestyle="-",label='High Density Auto')
#P.plot(cl_in_masked[5]*model[5]*ell[5], color=colors[4], lw=2, linestyle="--",label='High Density Cross')
P.plot(input_mean*model[0]*ell[0], color=colors[2], lw=2, linestyle="-",label='Noisy Auto')
P.plot(input_mean*model[1]*ell[1], color=colors[2], lw=2, linestyle="--",label='Noisy Cross')
P.plot(input_mean*model[2]*ell[2], color=colors[3], lw=2, linestyle="-",label='Noiseless Auto')
P.plot(input_mean*model[3]*ell[3], color=colors[3], lw=2, linestyle="--",label='Noiseless Cross')
P.plot(input_mean*model[4]*ell[4], color=colors[4], lw=2, linestyle="-",label='High Density Auto')
P.plot(input_mean*model[5]*ell[5], color=colors[4], lw=2, linestyle="--",label='High Density Cross')

P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
#P.ylim([-0.8e-6, 1.1e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'upper right', fontsize=16)
P.title('Damped Masked Input x Large Angle Damping')

