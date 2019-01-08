## SY 31/10/18
## Reads in power spectra of 100 realisations of estimated kappa.
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
kxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.kappa
wkxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.wkappa
kxi_input = fits.open('kappa_input.fits')[1].data.I

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
#Cl1 = ['Cl_crosses.txt','Cl_crosses_noisy.txt', 'Cl_crosses_cut.txt', 'Cl_crosses_noisy_smoothed.txt','Cl_crosses_noisy_10_100.txt']
#Cl2 = ['Cl_inputs.txt', 'Cl_inputs_noisy.txt', 'Cl_inputs_cut.txt', 'Cl_inputs_noisy_smoothed.txt', 'Cl_inputs_noisy_10_100.txt']
Cl1 = ['Cl_crosses.txt','Cl_crosses_noisy.txt', 'Cl_crosses_cut.txt', 'Cl_crosses_noisy_smoothed.txt']
Cl2 = ['Cl_inputs.txt', 'Cl_inputs_noisy.txt', 'Cl_inputs_cut.txt', 'Cl_inputs_noisy_smoothed.txt']

cmean = []
cstdev = []
imean = []
chi_2 = []

##- Read in cross power spectra
for i,j in enumerate(Cl1):
    Cl_crosses_in = np.loadtxt(j)
    Cl_inputs = np.loadtxt(Cl2[i])
    ell, cr_mean, cr_stdev, in_mean, chi = get_vals(Cl_crosses_in, Cl_inputs)
    cmean.append(cr_mean)
    cstdev.append(cr_stdev)
    imean.append(in_mean)
    chi_2.append(chi)

chi2 = np.asarray(chi_2)
cross_stdev = np.asarray(cstdev)
cross_mean = np.asarray(cmean)
input_mean = np.asarray(imean)

print('chi2 (noiseless, 725,816)', chi2[0])
print('chi2 (noisy, 290,298)', chi2[1])
print('chi2 (noiseless, 290,823)', chi2[2])
print('chi2 (smoothed noisy, 290,298)', chi2[3])
#print('chi2 (noisy rtmax = 10, 290,298)', chi2[4])

##- fit for large-angles
Cl_true_corr = hp.sphtfunc.anafast(kxi, lmax=3*NSIDE-1)
y = (Cl_true_corr / Cl_xi_input)
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model = np.polyval(z, x)

##- Setup figures
P.rcParams.update({'font.size':20})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Plot bias figure
P.figure()

P.plot(input_mean[0]*x, color='k', linewidth=2.0, linestyle="-",label='Masked Input')
P.plot(input_mean[0]*model*x, color=colors[0], lw=2, linestyle="-",label='Masked Input x Large-angle damping')
#P.plot(input_mean[0]*model*x, color='k', lw=2, linestyle="-",label='Masked Input x Large-angle damping')
P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color='#006edb', marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor='#006edb', label='Noiseless Estimated x Input')

P.axhline(0., color='k', ls=':')
P.ylabel('$l \ C_l^{true, est}$', fontsize=28)
P.xlabel('$l$', fontsize=28)
P.xlim([0, 800])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=18)
P.title(r'True and Estimated ${\kappa}$ Correlations')

##- Plot errors figure
P.figure()

P.plot(input_mean[0]*model*x, color='k', lw=2, linestyle="-",label='Masked Input x Large-angle damping')

#P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color='#006edb', marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor='#006edb', label='Noiseless Estimated x Input')

#P.errorbar(ell + 2, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color='#6e00db', marker = 'o', markeredgecolor='#6e00db', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label='Noisy Estimated x Input')

#P.errorbar(ell + 4, cross_mean[2]*ell, xerr=None, yerr=cross_stdev[2]*ell, color='r', marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor='r', label='Cut Noiseless Estimated x Input')

#P.errorbar(ell + 6, cross_mean[3]*ell, xerr=None, yerr=cross_stdev[3]*ell, color='g', marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor='g', label='Smoothed Noisy Estimated x Input')

#P.errorbar(ell + 8, cross_mean[4]*ell, xerr=None, yerr=cross_stdev[4]*ell, color='magenta', marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor='magenta', label='Noisy rtmax=10')

P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[0], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[0], label='Noiseless Estimated x Input')

P.errorbar(ell + 2, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[1], marker = 'o', markeredgecolor=colors[1], fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label='Noisy Estimated x Input')

P.errorbar(ell + 4, cross_mean[2]*ell, xerr=None, yerr=cross_stdev[2]*ell, color=colors[2], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[2], label='Cut Noiseless Estimated x Input')

P.errorbar(ell + 6, cross_mean[3]*ell, xerr=None, yerr=cross_stdev[3]*ell, color=colors[3], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[3], label='Smoothed Noisy Estimated x Input')

#P.errorbar(ell + 8, cross_mean[4]*ell, xerr=None, yerr=cross_stdev[4]*ell, color=colors[4], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[4], label='Noisy rtmax=10')

P.axhline(0., color='k', ls=':')
P.ylabel('$l \ C_l^{true, est}$', fontsize=24)
P.xlabel('$l$', fontsize=24)
P.xlim([0, 800])
P.ylim([-1.0e-5,1.2e-5])
#P.yscale('symlog', linthreshy=1e-7)
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'upper right', fontsize=18)
P.title(r'True and Estimated ${\kappa}$ Correlations')


