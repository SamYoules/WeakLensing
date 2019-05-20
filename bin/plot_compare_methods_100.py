## SY 20/5/19
## Reads in power spectra of 100 realisations of estimated kappa for midpoint & optimal estimators.
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
kxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.kappa
wkxi = fits.open('kappa-gaussian-true-70-100.fits.gz')[1].data.wkappa
kxi_input = fits.open('est_maps_noiseless/kappa_input1.fits')[1].data.I

##- Get resolution of map
nside=int(hp.npix2nside(kxi.size))

##- Mask off areas outside DESI footprint
mask = wkxi!=0
mask &= (kxi>np.percentile(kxi[mask], 0.5)) & \
                (kxi<np.percentile(kxi[mask], 99.5))
kxi[~mask]=hp.UNSEEN

##- Reset resolution of input maps to match estimated maps
kxi_input = hp.ud_grade(kxi_input, nside)

##- Mask area outside footprint
kxi_input = kxi_input*(mask)+hp.UNSEEN*(~mask) 
kxi_input[~mask]=hp.UNSEEN

##- Get input Cls
Cl_xi_input = hp.sphtfunc.anafast(kxi_input, lmax=3*nside-1)

##- Get File names

Clx = ['Cls/Cl_crosses_noisy.txt', 'est_opt_noisy/Cls/Cl_crosses_noisy.txt', 'Cls/Cl_crosses_cut.txt', 'est_opt_cut/Cls/Cl_crosses_cut.txt']
Cli = 'Cls/Cl_inputs_noisy.txt'
cmean = []
cstdev = []
imean = []
chi_2 = []

Cl_inputs = np.loadtxt(Cli)

##- Get normalisation factor for optimal by matching it to midpoint at ell=400
Cl_mid_noisy = np.loadtxt(Clx[0])
Cl_opt_noisy = np.loadtxt(Clx[1])
#norm_noisy = Cl_mid_noisy[0][300]/Cl_opt_noisy[0][300]
norm_noisy = np.sum(Cl_mid_noisy)/np.sum(Cl_opt_noisy)
Cl_mid_cut = np.loadtxt(Clx[2])
Cl_opt_cut = np.loadtxt(Clx[3])
norm_cut = np.sum(Cl_mid_cut)/np.sum(Cl_opt_cut)


##- Read in cross- and input- power spectra
for i,j in enumerate(Clx):
    Cl_crosses_in = np.loadtxt(j)
    if i==1:
        Cl_crosses_in *= norm_noisy
    elif i==3:
        Cl_crosses_in *= norm_cut
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
chi22a='%.3f'%(chi2[2])
chi22b='%.3f'%(chi2[3])

chi2noisy=r'$\chi^2_{midpoint} = $' + chi21a + ', ' + r'$\chi^2_{optimal} = $' + chi21b

chi2cut=r'$\chi^2_{midpoint} = $' + chi22a + ', ' + r'$\chi^2_{optimal} = $' + chi22b

##- fit for large-angles
Cl_true_corr = hp.sphtfunc.anafast(kxi, lmax=3*nside-1)
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
P.plot(input_mean[0]*x, color=colors[2], linewidth=2.0, linestyle="-",label='Masked Input')
P.plot(input_mean[0]*model*x, color=colors[3], lw=2, linestyle="-",label='Masked Input x large-angle damping')
P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[0], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[0], label='Midpoint Noisy Estimated x Input')
P.errorbar(ell + 4, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[1], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[1], label='Optimal Noisy Estimated x Input')

P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.ylim([-1.0e-6, 1.5e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=16, fancybox=True, framealpha=0.5)
#P.tight_layout()
P.savefig('plots/compare_bias.pdf')


##- Plot errors figure
# Create 2x1 sub plots and set right-hand axis labels
P.rcParams.update({'font.size':14})
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

gs = gridspec.GridSpec(2, 1, hspace = 0.2, wspace=0.3)
fig = P.figure(figsize=(7,8))
ax1 = P.subplot(gs[0, 0]) # row 0, col 0
ax1.text(1.02, 0.5, "Noisy",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

ax2 = P.subplot(gs[1, 0]) # row 1, col 0
ax2.text(1.02, 0.5, "Noiseless",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

##- Subplot 1
##- Masked Input x Large-angle damping
ax1.plot(input_mean[0]*model*x, color=colors[3], lw=2, linestyle="-",label='Damped Masked Input')

ax1.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[0], marker = 'o', markeredgecolor=colors[0], fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label=r'$\langle \kappa_{midpoint} \kappa_{\rm{input}} \rangle$')

ax1.errorbar(ell + 4, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[1], marker = 'o', markeredgecolor=colors[1], fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label=r'$\langle \kappa_{optimal} \kappa_{\rm{input}} \rangle$')

ax1.axhline(0., color='k', ls=':')
ax1.set_ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
#ax1.set_xlabel(r'$\ell$', fontsize=18)
ax1.set_xlim([0, 800])
ax1.set_ylim([-0.08e-5,0.2e-5])
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
anchored_text = AnchoredText(chi2noisy, loc=1)
ax1.add_artist(anchored_text)

ax1.legend(numpoints = 1, bbox_to_anchor=(0.5, 1.25), loc="upper center", borderaxespad=0., fontsize=14, ncol=3, handletextpad=0.5, handlelength=1, columnspacing=1)

##- Subplot 2
ax2.plot(input_mean[0]*model*x, color=colors[3], lw=2, linestyle="-",label='Damped Masked Input')

ax2.errorbar(ell, cross_mean[2]*ell, xerr=None, yerr=cross_stdev[2]*ell, color=colors[0], marker = 'o', markeredgecolor=colors[0], fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label=r'$\langle \kappa_{midpoint} \kappa_{\rm{input}} \rangle$ Noiseless')

ax2.errorbar(ell + 4, cross_mean[3]*ell, xerr=None, yerr=cross_stdev[3]*ell, color=colors[1], marker = 'o', markeredgecolor=colors[1], fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, label=r'$\langle \kappa_{optimal} \kappa_{\rm{input}} \rangle$ Noiseless')

ax2.axhline(0., color='k', ls=':')
ax2.set_xlabel(r'$\ell$', fontsize=18)
ax2.set_ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
ax2.set_xlim([0, 800])
ax2.set_ylim([-0.08e-5,0.2e-5])
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
anchored_text = AnchoredText(chi2cut, loc=1)
ax2.add_artist(anchored_text)

P.savefig('plots/compare_100.pdf', bbox_inches='tight')

