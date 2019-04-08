# SY 3/12/18
# Plots C_ells for estimated and input maps for a range of maximum r_transverse
# separations. r_parallel is set at 100 Mpc/h for all. The third sublot shows
# the ratio of est/input, which is used to model the large angle effect in
# get_errors6.py. Output is rtmax.png.

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

def get_Cls(rtmax, k_input, NSIDE):
    '''Open estimated kappa files, mask them and calculate their C_ells'''
    ##- Open kappa true-correlation files
    kest = fits.open('kappa-gaussian-true-{}-100.fits.gz'.format(rtmax))[1].data.kappa
    wkest = fits.open('kappa-gaussian-true-{}-100.fits.gz'.format(rtmax))[1].data.wkappa

    ##- Mask off areas outside DESI footprint
    mask = wkest!=0
    mask &= (kest>np.percentile(kest[mask], 0.5)) & \
                (kest<np.percentile(kest[mask], 99.5))
    kest[~mask]=hp.UNSEEN
    k_input = k_input*(mask)+hp.UNSEEN*(~mask) 
    k_input[~mask]=hp.UNSEEN
    Cl_input = hp.sphtfunc.anafast(k_input, lmax=3*NSIDE-1)
    Cl_estim = hp.sphtfunc.anafast(kest, lmax=3*NSIDE-1)
    return Cl_input, Cl_estim


k_input = fits.open('est_maps_10_100/kappa_input1.fits')[1].data.I
NSIDE=int(hp.npix2nside(k_input.size))
Cl_in = []
Cl_est = []
a = [10, 20, 30, 40, 50, 60, 70, 80, 100, 120]
for i in a:
    inp, est = get_Cls(i, k_input, NSIDE)
    Cl_in.append(inp)
    Cl_est.append(est)

##- Setup figures
P.rcParams.update({'font.size':14})
P.ion()
ncolors=10
colors = P.cm.jet(np.linspace(0,1,ncolors))

##- Setup 3 subplots with right-hand axis labels
gs = gridspec.GridSpec(3, 1, hspace = 0.2, wspace=0.3)
fig = P.figure(figsize=(7,11))

ax1 = P.subplot(gs[0, 0]) # row 0, col 0
ax1.text(1.02, 0.5, "Estimated",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

ax2 = P.subplot(gs[1, 0]) # row 1, col 0
ax2.text(1.02, 0.5, "Input (masked)",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

ax3 = P.subplot(gs[2, 0]) # row 2, col 0
ax3.text(1.02, 0.5, "Estimated / Input",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)
P.setp(ax3, xticks=[0, 200, 400, 600, 800], xlim=([0, 800]))


##- Plot figure
for i in np.arange(10):
    ax1.plot(Cl_est[i], color=colors[i], lw=2, linestyle="-", label="rtmax={}".format(a[i]))
    ax2.plot(Cl_in[i], color=colors[i], lw=2, linestyle="-", label="rtmax={}".format(a[i]))
    ax3.plot(Cl_est[i]/Cl_in[i], color=colors[i], lw=2, linestyle="-", label="rtmax={}".format(a[i]))

ax1.set_ylabel(r'$\ell C_{\ell}^{\kappa \kappa}$', fontsize=18)
ax2.set_ylabel(r'$\ell C_{\ell}^{\kappa \kappa}$', fontsize=18)
ax3.set_ylabel(r'$\ell C_{\ell}^{\kappa \kappa}$', fontsize=18)
ax1.set_ylim([0,2.5e-8])
ax2.set_ylim([0,2.5e-8])
ax3.set_ylim([0,3.5])
ax1.set_xticks([0, 200, 400, 600, 800])
ax2.set_xticks([0, 200, 400, 600, 800])
ax3.set_xticks([0, 200, 400, 600, 800])
ax1.set_xlim([0, 800])
ax2.set_xlim([0, 800])
ax3.set_xlim([0, 800])
ax3.set_xlabel(r'$\ell$', fontsize=18)

#axs[0].set_yscale('symlog', linthreshy=1e-8)
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'upper right', fontsize=14, ncol=2)

#P.savefig('plots/rtmax_vertical.pdf', bbox_inches='tight')

