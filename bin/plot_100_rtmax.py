## SY 13/6/19
## Reads in power spectra of 100 realisations of estimated kappa for LYAxLYA and QSOxLYA and input kappa.
## Plots cross-cf error bars for input and estimated kappa, for different values of rt_max.

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

    w = (ell>=ellmin)&(ell<ellmax)
    index = np.floor( (ell[w]-ellmin)*1./(ellmax-ellmin)*nell ).astype(int)
    well = np.bincount( index, weights=weights[w])
    sell = np.bincount( index, weights=weights[w]*ell[w])
    scl  = np.bincount( index, weights=weights[w]*Cls[w])
    ell = sell/well
    Cls = scl/well
        
    return ell, Cls, well 


def get_vals(Cl_crosses_in, Cl_inputs):

    ##- Trim off smallest and largest scales (gives 7 bins of ell=104)
    Cl_crosses = Cl_crosses_in

    ##- Rebin
    cross_rebin = []
    cross_weights = []

    for i,j in enumerate(Cl_crosses):
        ell, crosses, cross_wei = rebin(Cl_crosses[i], ellmin=50, ellmax=750, nell=7)
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
    chi2_cov = cross_mean.dot( np.linalg.inv(cross_cov).dot(cross_mean))
    return ell, cross_mean, cross_stdev, input_mean, chi2_cov


##- Get true-correlation Cls
Clxi = ['maps/midpoint/true_corr/Cls/Cl_crosses_xnoisy_rt30.txt', 'maps/midpoint/true_corr/Cls/Cl_crosses_xnoisy_rt70.txt', 'maps/midpoint/true_corr/Cls/Cl_crosses_xnoisy_rt100.txt', 'maps/midpoint/true_corr/Cls/Cl_crosses_xnoisy_rt150.txt']
Clxi30 = np.loadtxt(Clxi[0])
Clxi70 = np.loadtxt(Clxi[1])
Clxi100 = np.loadtxt(Clxi[2])
Clxi150 = np.loadtxt(Clxi[3])


##- Get input Cls
Cl_xii = 'maps/input/Cl_xi_input.txt'
Cl_xi_input = np.loadtxt(Cl_xii)
Cli = 'maps/input/Cl_inputs.txt'

##- Get Cls for 100 realisations
Clx = ['maps/midpoint/xnoisy30/Cls/Cl_crosses.txt', 'maps/midpoint/xnoisy/Cls/Cl_crosses.txt', 'maps/midpoint/xnoisy100/Cls/Cl_crosses.txt', 'maps/midpoint/xnoisy150/Cls/Cl_crosses.txt']

cmean = []
cstdev = []
imean = []
chi_2 = []

##- Read in cross- and input- power spectra
for i,j in enumerate(Clx):
    Cl_crosses_in = np.loadtxt(j)
    Cl_inputs = np.loadtxt(Cli)
    ell, cr_mean, cr_stdev, in_mean, chi = get_vals(Cl_crosses_in, Cl_inputs)
    cmean.append(cr_mean)
    cstdev.append(cr_stdev)
    imean.append(in_mean)
    chi_2.append(chi)

chi2 = np.asarray(chi_2)
cross_stdev = np.asarray(cstdev)
cross_mean = np.asarray(cmean)
input_mean = np.asarray(imean)

chi2a='%.3f'%(chi2[0])
chi2b='%.3f'%(chi2[1])
chi2c='%.3f'%(chi2[2])
chi2d='%.3f'%(chi2[3])

chi2=r'$\chi^2_{30} = $' + chi2a + r'$,\ \ \chi^2_{70} = $' + chi2b + '\n' + r'$\chi^2_{100} = $' + chi2c + r', $\chi^2_{150} = $' + chi2d

##- fit for large-angles
y = (Clxi30 / Cl_xi_input)*-1
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model30 = np.polyval(z, x)
y = (Clxi70 / Cl_xi_input)*-1
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model70 = np.polyval(z, x)
y = (Clxi100 / Cl_xi_input)*-1
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model100 = np.polyval(z, x)
y = (Clxi150 / Cl_xi_input)*-1
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model150 = np.polyval(z, x)

##- Setup figures
P.rcParams.update({'font.size':18})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
#colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

##- Plot bias figure
#P.figure(figsize=(8.2,6))
fig, ax = P.subplots(figsize=(10,6))
P.plot(input_mean[0]*model30*x, color=colors[0], lw=2, linestyle="-")
P.plot(input_mean[1]*model70*x, color=colors[1], lw=2, linestyle="-")
P.plot(input_mean[2]*model100*x, color=colors[2], lw=2, linestyle="-")
P.plot(input_mean[3]*model150*x, color=colors[3], lw=2, linestyle="-")
P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[0], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[0], label='rt_max = 30')
P.errorbar(ell + 5, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[1], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[1], label='rt_max = 70')
P.errorbar(ell + 10, cross_mean[2]*ell, xerr=None, yerr=cross_stdev[2]*ell, color=colors[2], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[2], label='rt_max = 100')
P.errorbar(ell + 15, cross_mean[3]*ell, xerr=None, yerr=cross_stdev[3]*ell, color=colors[3], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[3], label='rt_max = 150')


props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.5, 0.95, chi2, transform=ax.transAxes, fontsize=16, verticalalignment='top', bbox=props)

P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=22)
P.xlabel(r'$\ell$', fontsize=22)
P.xlim([0, 800])
P.ylim([-0.8e-6, 1.1e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=16, fancybox=True, framealpha=0.5, ncol=2)
P.tight_layout()
P.savefig('plots/compare_rtmax_xnoisy_cov.pdf')

