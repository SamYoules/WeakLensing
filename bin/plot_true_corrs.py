## SY 27/6/19
## Plots true_corr power spectra for different datasets.

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


def get_model(c_ell):
    y = (c_ell / cl_input)
    x = np.arange(y.size)
    z = np.polyfit(x, y, 5)
    model = np.polyval(z, x)
    return model, x


esttype = 'midpoint'

##- Open kappa true-correlation and input c_ells
cl_noisy = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_autos_noisy_rt70.txt')
clx_noisy = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_crosses_noisy_rt70.txt')

cl_cut = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_autos_cut_rt70.txt')
clx_cut = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_crosses_cut_rt70.txt')

cl_noiseless = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_autos_noiseless_rt70.txt')
clx_noiseless = np.loadtxt('maps/midpoint/true_corr/Cls/Cl_crosses_noiseless_rt70.txt')

cl_input = np.loadtxt('maps/input/Cl_xi_input.txt')
input_mean = np.loadtxt('maps/input/Cl_input_mean.txt')

x = []
model = []
kappa_true = [cl_noisy, clx_noisy, cl_cut, clx_cut, cl_noiseless, clx_noiseless]
for i in kappa_true:
    mod, ell_true = get_model(i)
    model.append(mod)
    x.append(ell_true)
    

##- Setup figures
P.rcParams.update({'font.size':18})
P.ion()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))
#colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

##- Plot figure
P.figure(figsize=(8.2,6))
#P.plot(input_mean[0], color=colors[1], linewidth=2.0, linestyle="-",label='Masked Input')

P.plot(cl_noisy, color=colors[2], lw=2, linestyle="-",label='Noisy Auto')
P.plot(clx_noisy, color=colors[2], lw=2, linestyle="--",label='Noisy Cross')
P.plot(cl_cut, color=colors[3], lw=2, linestyle="-",label='Noiseless Auto')
P.plot(clx_cut, color=colors[3], lw=2, linestyle="--",label='Noiseless cross')
P.plot(cl_noiseless, color=colors[4], lw=2, linestyle="-",label='High Density Auto')
P.plot(clx_noiseless, color=colors[4], lw=2, linestyle="--",label='High Density Cross')

P.title('True_corr')
P.axhline(0., color='k', ls=':')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{true, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
#P.ylim([0.0, 1.1e-6])
P.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
handles, labels = P.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
P.legend(by_label.values(), by_label.keys(), numpoints = 1, loc = 'upper right', fontsize=16)


