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


##- Open kappa true-correlation files
kxi = fits.open('maps/true_corr/kappa_noisy_lensed_true_rt70_rp100.fits.gz')[1].data.kappa
wkxi = fits.open('maps/true_corr/kappa_noisy_lensed_true_rt70_rp100.fits.gz')[1].data.wkappa
kxi_input = fits.open('maps/input/kappa_input1.fits')[1].data.I

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

##- Get File names

Clx = ['maps/midpoint/xnoisy/Cls/Cl_crosses.txt', '../desi_old/Cls/Cl_crosses_noisy.txt']
Cli = ['maps/input/Cl_inputs.txt', '../desi_old/Cls/Cl_inputs_noisy.txt']

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

chi21a='%.3f'%(chi2[0])
chi21b='%.3f'%(chi2[1])

chi2new=r'$\chi^2_{cross} = $' + chi21a + r', $\chi^2_{auto} = $' + chi21b

##- fit for large-angles
Cl_true_corr = hp.sphtfunc.anafast(kxi, lmax=3*NSIDE-1)
y = (Cl_true_corr / Cl_xi_input)
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model = np.polyval(z, x)

##- Setup figures
P.rcParams.update({'font.size':18})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
#colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

##- Plot bias figure
P.figure(figsize=(8.2,6))
fig, ax = plt.subplots()
P.plot(input_mean[0]*x, color=colors[1], linewidth=2.0, linestyle="-",label='Masked Input')
P.plot(input_mean[0]*model*x, color=colors[3], lw=2, linestyle="-",label='Masked Input x large-angle damping')
P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[0], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[0], label='QSO x Lya x Input')
P.errorbar(ell + 8, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[2], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[2], label='Auto x Input')

props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.55, 0.95, chi2new, transform=ax.transAxes, fontsize=18,
        verticalalignment='top', bbox=props)

P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=22)
P.xlabel(r'$\ell$', fontsize=22)
P.xlim([0, 800])
P.ylim([-0.8e-6, 1.1e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=18, fancybox=True, framealpha=0.5)
#P.tight_layout()
#P.savefig('plots/bias.pdf')

