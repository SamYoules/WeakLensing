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

def rebin(cls,ellmin=50, ellmax=750, nell=7):
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

def get_vals(cl_crosses_in):

    ##- Trim off smallest and largest scales (gives 7 bins of ell=104)
    #cl_crosses = cl_crosses_in[:,40:768]

    ##- Rebin
    cross_rebin = []
    cross_weights = []

    for i,j in enumerate(cl_crosses_in):
        ell, crosses, cross_wei = rebin(cl_crosses_in[i].ravel(), ellmin=48, ellmax=768, nell=7)
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
    chi2_cov = cross_mean.dot( np.linalg.inv(cross_cov).dot(cross_mean))
    return ell, cross_mean, cross_stdev, chi2_cov

def get_model(kappa_type, cl_xi_input):
    cl_true_corr = hp.sphtfunc.anafast(kappa_type, lmax=3*nside-1)
    y = (cl_true_corr / cl_xi_input)
    x = np.arange(y.size)
    z = np.polyfit(x, y, 5)
    model = np.polyval(z, x)
    return model, x


esttype = 'midpoint'

##- Open kappa true-correlation files and corresponding input map

kappa_type = ['noisy', 'xnoisy', 'cut', 'xcut', 'noiseless', 'xnoiseless']
wkxi_noisy = fits.open('maps/midpoint/true_corr/kappa_noisy_rt70.fits.gz')[1].data.wkappa

kxi=[]
for i, k in enumerate(kappa_type):
    kxi.append(fits.open('maps/midpoint/true_corr/kappa_{}_rt70.fits.gz'.format(k))[1].data.kappa)

kinput = fits.open('maps/input/kappa_input1.fits')[1].data.I

##- Get resolution of map
nside=int(hp.npix2nside(kxi[0].size))
kinput = hp.ud_grade(kinput, nside)

##- Mask off areas outside DESI footprint
mask = wkxi_noisy!=0
mask &= (kxi[0]>np.percentile(kxi[0][mask], 0.5)) & \
                (kxi[0]<np.percentile(kxi[0][mask], 99.5))

for i, k in enumerate(kxi):
    kxi[i][~mask]=hp.UNSEEN

kinput = kinput*(mask)+hp.UNSEEN*(~mask) 
kinput[~mask]=hp.UNSEEN

##- Create models for damped input
cl_xi_input = hp.sphtfunc.anafast(kinput, lmax=3*nside-1)
x = []
model = []
for i in kxi:
    mod, ell_true = get_model(i, cl_xi_input)
    model.append(mod)
    x.append(ell_true)

cl_inputs = np.loadtxt('maps/input/Cl_inputs.txt')
input_mean = (np.sum(cl_inputs, axis=0) / 100)

cmean = []
cstdev = []
chi_2 = []
maptype = ['noisy', 'xnoisy', 'cnoisy', 'cut', 'xcut', 'ccut',
                    'noiseless', 'xnoiseless', 'cnoiseless', ]

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

chi2text=[]
j=0
for i in range(3):
    chi2text.append(r'$\chi^2_{\alpha \alpha} = $' + '%.3f'%(chi2[j]) + ', '
                      + r'$\chi^2_{q\alpha} = $' + '%.3f'%(chi2[j+1]) + ', '
                      + r'$\chi^2_{comb} = $' + '%.3f'%(chi2[j+2]))
    j+=3

x = []
model = []
for i in kxi:
    mod, ell_true = get_model(i, cl_xi_input)
    model.append(mod)
    x.append(ell_true)
    

##- Setup figures
P.rcParams.update({'font.size':14})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
#colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

##- Plot bias figure
fig, ax = P.subplots(figsize=(10,8))

P.plot(input_mean*model[0]*x[0], color="k", lw=2, linestyle="-",label='Masked Input x large-angle damping')
P.errorbar(ell, cross_mean[0]*ell, xerr=None, yerr=cross_stdev[0]*ell, color=colors[2], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[2], label='Noisy Auto x Input')
P.errorbar(ell + 6, cross_mean[1]*ell, xerr=None, yerr=cross_stdev[1]*ell, color=colors[3], marker = 'o', fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, markeredgecolor=colors[3], label='Noisy QSO-Cross x Input')
anchored_text = AnchoredText(chi2text[1], loc=3)
ax.add_artist(anchored_text)

P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
#P.ylim([-0.5e-6, 1.3e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'upper right', fontsize=16)
P.savefig('plots/bias_190808.pdf')


##- Plot errors figure
mapnames = ['Noisy', 'Noiseless', 'High Density']

# Create 2x2 sub plots and set right-hand axis labels
#fig = P.subplots(figsize=(13,8))
#gs = gridspec.GridSpec(2, 2, hspace = 0.2, wspace=0.3)

# Create 3x1 sub plots and set right-hand axis labels
fig = P.figure(figsize=(7,12))
gs = gridspec.GridSpec(3, 1, hspace = 0.3, wspace=0.)

axes=[]
axes.append(P.subplot(gs[0, 0])) # row 0, col 0
axes[0].text(1.02, 0.5, mapnames[0],
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)
axes.append(P.subplot(gs[1, 0])) # row 1, col 0
axes[1].text(1.02, 0.5, mapnames[1],
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)
#axes.append(P.subplot(gs[1, 1])) # row 1, col 1
axes.append(P.subplot(gs[2, 0])) # row 2, col 0
axes[2].text(1.02, 0.5, mapnames[2],
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

j=0
for i,k in enumerate([0,3,6]):
    axes[i].errorbar(ell, cross_mean[k]*ell,
         xerr=None, yerr=cross_stdev[k]*ell, 
	 color=colors[0], marker = 'o', markeredgecolor=colors[0], 
	 fmt='.', capsize=5, elinewidth=2, markeredgewidth=2, 
	 label=r'$\langle \kappa_{\alpha \alpha} \kappa_{\rm{input}} \rangle$')
    axes[i].errorbar(ell + 4, cross_mean[k+1]*ell,
         xerr=None, yerr=cross_stdev[k+1]*ell,
         color=colors[1], marker = 'o', markeredgecolor=colors[1],
         fmt='.', capsize=5, elinewidth=2, markeredgewidth=2,
         label=r'$\langle \kappa_{q \alpha} \kappa_{\rm{input}} \rangle$')
    axes[i].errorbar(ell + 8, cross_mean[k+2]*ell,
         xerr=None, yerr=cross_stdev[k+2]*ell,
         color=colors[2], marker = 'o', markeredgecolor=colors[2],
         fmt='.', capsize=5, elinewidth=2, markeredgewidth=2,
         label=r'$\langle \kappa_{comb} \kappa_{\rm{input}} \rangle$')
    axes[i].plot(input_mean*model[j]*x[j],
         color=colors[4], lw=2, linestyle="-",
         label='Damped Masked Input (Auto)')
    axes[i].plot(input_mean*model[j+1]*x[j+1],
         color=colors[3], lw=2, linestyle="-",
         label='Damped Masked Input (Cross)')
    axes[i].plot(input_mean*model[j+1]*x[j+1],
         color=colors[3], lw=2, linestyle="-",
         alpha=0, label=' ') # this is a dirty hack to better format the legend
    j+=2

    axes[i].axhline(0., color='k', ls=':')
    axes[i].set_ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
    axes[i].set_xlim([0, 800])
    axes[i].set_ylim([-0.08e-5,0.2e-5])
    axes[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
    anchored_text = AnchoredText(chi2text[i], loc=1)
    axes[i].add_artist(anchored_text)
    #if i!=0:
    #    axes[i].set_xlabel(r'$\ell$', fontsize=18)
    axes[2].set_xlabel(r'$\ell$', fontsize=18)

    #axes[0].legend(numpoints=1, bbox_to_anchor=(1.7, 1), loc="upper center",
    #     borderaxespad=0., fontsize=14, ncol=1, handletextpad=0.5,
    #     handlelength=1, columnspacing=1)

    axes[0].legend(numpoints=1, bbox_to_anchor=(0.5, 1.5), loc="upper center",
         borderaxespad=0., fontsize=14, ncol=2, handletextpad=0.5,
         handlelength=1, columnspacing=1)

P.savefig('plots/100_realisations.pdf', bbox_inches="tight")

