## SY 28/6/19

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

def rebin(cls,ellmin=40, ellmax=768, nell=7):
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
        
    return ell, cls, well 

def get_vals(cl_crosses_in):

    ##- Trim off smallest and largest scales (gives 7 bins of ell=104)
    cl_crosses = cl_crosses_in[:,40:768]

    ##- Rebin
    cross_rebin = []
    cross_weights = []

    for i,j in enumerate(cl_crosses):
        ell, crosses, cross_wei = rebin(cl_crosses[i])
        cross_rebin.append(crosses)
        cross_weights.append(cross_wei)

    cross_cl = np.asarray(cross_rebin)
    cross_we = np.asarray(cross_weights)
    cross_mean = (np.sum(cross_cl*cross_we,axis=0) / np.sum(cross_we,axis=0))

    ##- Calculate covariance matrices, variance and errors
    QW = (cross_cl-cross_mean)*cross_we
    cross_cov = QW.T.dot(QW)/cross_we.T.dot(cross_we)
    cross_var = sp.diagonal(cross_cov)
    cross_stdev = np.sqrt(cross_var)
    #cross_stdev = np.sqrt(cross_var/100)

    ##- Calc chi2
    chi2 = np.sum((cross_mean - 0.)**2/cross_stdev**2 )

    return ell, cross_mean, cross_stdev, chi2

def get_model(kappa_type):
    cl_true_corr = hp.sphtfunc.anafast(kappa_type, lmax=3*nside-1)
    y = (cl_true_corr / cl_xi_input)
    x = np.arange(y.size)
    z = np.polyfit(x, y, 5)
    model = np.polyval(z, x)
    return model, x


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

##- Get File names
esttype = 'midpoint'
maptype = ['noisy', 'xnoisy', 'cnoisy', 'cut', 'xcut', 'ccut',
                    'noiseless', 'xnoiseless', 'cnoiseless', ]
cl_inputs = np.loadtxt('maps/input/Cl_inputs.txt')
input_mean = (np.sum(cl_inputs, axis=0) / 100)

cmean = []
cstdev = []
imean = []
chi_2 = []

##- Read in cross- and input- power spectra
for i,j in enumerate(maptype):
    cl_crosses_in = np.loadtxt('maps/{}/{}/Cls/Cl_crosses.txt'.format(esttype,j))
    ell, cr_mean, cr_stdev, chi = get_vals(cl_crosses_in)
    cmean.append(cr_mean)
    cstdev.append(cr_stdev)
    chi_2.append(chi)

chi2 = np.asarray(chi_2)
cross_stdev = np.asarray(cstdev)
cross_mean = np.asarray(cmean)

chi21a='%.3f'%(chi2[0])
chi21b='%.3f'%(chi2[1])
chi21c='%.3f'%(chi2[2])
chi22a='%.3f'%(chi2[3])
chi22b='%.3f'%(chi2[4])
chi22c='%.3f'%(chi2[5])
chi23a='%.3f'%(chi2[6])
chi23b='%.3f'%(chi2[7])
chi23c='%.3f'%(chi2[8])

##- Setup figures
P.rcParams.update({'font.size':18})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
gs = gridspec.GridSpec(2, 1, hspace = 0.2, wspace=0.3)
fig = P.figure(figsize=(8,8))
ax1 = P.subplot(gs[0, 0]) # row 0, col 0
ax1.text(1.02, 0.5, "LYAxLYA",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)
ax2 = P.subplot(gs[1, 0]) # row 1, col 0
ax2.text(1.02, 0.5, "LYAxQSO",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)
ax1.plot(ell, cstdev[0]*ell, color=colors[1], linewidth=2.0, linestyle="-",label='Noisy')
ax1.plot(ell, cstdev[2]*ell, color=colors[2], linewidth=2.0, linestyle="-",label='Noiseless')
ax1.plot(ell,cstdev[4]*ell, color=colors[3], linewidth=2.0, linestyle="-",label='High Density')
ax1.axhline(0., color='k', ls=':')
ax1.set_ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=18)
ax1.set_xlabel(r'Errors', fontsize=18)
#ax1.set_xlim([0, 800])
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=16)
ax2.plot(ell, cstdev[1]*ell, color=colors[1], linewidth=2.0, linestyle="-",label='Noisy')
ax2.plot(ell, cstdev[3]*ell, color=colors[2], linewidth=2.0, linestyle="-",label='Noiseless')
ax2.plot(ell, cstdev[5]*ell, color=colors[3], linewidth=2.0, linestyle="-",label='High Density')
ax2.axhline(0., color='k', ls=':')
ax2.set_ylabel(r'Errors', fontsize=18)
ax2.set_xlabel(r'$\ell$', fontsize=18)
#ax2.set_xlim([0, 800])
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax2.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'lower left', fontsize=16)



