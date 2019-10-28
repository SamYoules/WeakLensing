#!/usr/bin/env python

#-- run mc to calculate delta r

import numpy as np
from numpy import linalg
import scipy as sp
from scipy import interpolate
import fitsio
import h5py
import argparse
import pylab as plt
from astropy.table import Table

def load_model(file_xi, file_fit, nbins=50) :

    h = fitsio.FITS(file_xi)
    ff = h5py.File(file_fit, 'r')
    base = file_xi 
    fit = ff[base+'/fit'][...] 
    data_rp = h[1]['RP'][:]
    data_rt = h[1]['RT'][:]
    hh = h[1].read_header()
    rpmin = hh['RPMIN']
    rpmax = hh['RPMAX']
    rtmin = 0 
    rtmax = hh['RTMAX']
    h.close()
    ff.close()

    rpmin = data_rp.reshape(50, 50)[0].max()
    rpmax = data_rp.reshape(50, 50)[-1].min()
    rtmin = data_rt.reshape(50, 50)[:, 0].max()
    rtmax = data_rt.reshape(50, 50)[:, -1].min()

    #-- create the regular grid for griddata
    rp = np.linspace(rpmin, rpmax, nbins)
    rt = np.linspace(rtmin, rtmax, nbins)

    xim = sp.interpolate.griddata((data_rt, data_rp), fit,
                (np.outer(np.ones(rp.size), rt).ravel(),
                np.outer(rp, np.ones(rt.size)).ravel()), method='cubic')

    #-- create interpolator object
    xi2d = sp.interpolate.RectBivariateSpline(rt, rp, \
               xim.reshape((50, 50)).T )

    return xi2d

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Improve fit of xcf for lya-lya.')

    parser.add_argument('--rp', required=True, type=int, \
               help='r parallel')
    parser.add_argument('--rt', required=True, type=int, \
               help='r transverse')
    parser.add_argument('--nreal', required=True, type=int, default=1000000, \
               help='number of delta pair realisations')
    parser.add_argument('--mcreal', required=True, type=int, default=100, \
               help='number of mc realisations')
    parser.add_argument('--var', required=True, type=float, default=0.005, \
               help='variance of deltas (approx 0.04 at z=2 to 0.15 at z=3)')
    parser.add_argument('--xi', required=True, \
               help='text file containing model')
    parser.add_argument('--fit', required=True, \
               help='text file containing correlation function')

    args, unknown = parser.parse_known_args()

    rp_in = args.rp
    rt_in = args.rt
    nrealisations = args.nreal
    mc_realisations = args.mcreal
    var1 = args.var
    var2 = args.var
    xi_file = args.xi
    fit_file = args.fit

    # Make range of drt to plot
    drt_in = [-0.02, -0.014, -0.008, -0.002, 0.004, 0.01, 0.016, 0.022]

    # Get xi_model interpolator object
    print('Reading model from ', fit_file)
    xi2d = load_model(xi_file, fit_file)

    # Calculate delta_rt using Monte Carlo for lensed & unlensed xi
    drt_l_mean = []
    drt_u_mean = []
    drt_l_std  = []
    drt_u_std  = []

    for drt in drt_in:
        drt_u_out = []
        drt_l_out = []
        for n in range(mc_realisations):
            # Construct two columns of Gaussian random variables
            g = np.random.randn(2, nrealisations)

            # Construct Covariance matrix
            xi_un = xi2d(rt_in, rp_in, grid=False)
            C = np.array([ [var1, xi_un  ],
                           [xi_un,   var2] ])

            # Cholesky decomposition, see 2.1 in 1108.5606
            L = np.linalg.cholesky(C)

            # Apply transformation to correlate variables
            delta = L.dot(g)

            # Apply lensing (approximation ignoring shear)
            rt_lens = rt_in + drt

            # Estimate drt_out
            xi_lensed = xi2d(rt_lens, rp_in)
            xi_unlensed = xi2d(rt_in, rp_in)
            R_l = 1 / rt_lens / xi2d(rt_lens, rp_in, dx=1)
            R_u = 1 / rt_lens / xi2d(rt_in, rp_in, dx=1)
           
            drt_u_out.append(np.mean((delta[0]*delta[1] - xi_unlensed) * R_u))
            drt_l_out.append(np.mean((delta[0]*delta[1] - xi_lensed) * R_l))

        drt_u_out = np.array(drt_u_out)
        drt_l_out = np.array(drt_l_out)

        drt_u_mean.append(np.mean(drt_u_out))
        drt_u_std.append(np.std(drt_u_out))
        drt_l_mean.append(np.mean(drt_l_out))
        drt_l_std.append(np.std(drt_l_out))

    drt_u_mean = np.array(drt_u_mean)
    drt_l_mean = np.array(drt_l_mean)
    drt_u_std = np.array(drt_u_std)
    drt_l_std = np.array(drt_l_std)
    residuals_u = (drt_u_mean - drt_in)/drt_u_std
    residuals_l = (drt_l_mean - drt_in)/drt_l_std

    ##- Setup figures
    plt.ion()
    plt.rcParams.update({'font.size':14})
    colors=['#396AB1','#DA7C30','#3E9651','#CC2529','#535154','#6B4C9A','#922428','#948B3D']

    ##- Plot bias figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8,8))
    fig.subplots_adjust(hspace=0.5)
    ax[0].set_title(fr'$\Delta r$, 1000000 realisations, rt={rt_in}, rp={rp_in}')
    ax[0].plot(drt_in, drt_in, color=colors[4], ls=':')
    ax[0].errorbar(drt_in, drt_u_mean, yerr=drt_u_std,
         color=colors[0], marker = 'o', markeredgecolor=colors[0],
         fmt='.', capsize=5, elinewidth=2, markeredgewidth=2,
         label=r'Unlensed $\xi$')
    ax[0].errorbar(drt_in, drt_l_mean, yerr=drt_l_std,
         color=colors[1], marker = 'o', markeredgecolor=colors[1],
         fmt='.', capsize=5, elinewidth=2, markeredgewidth=2,
         label=r'Lensed $\xi$')
    ax[0].legend()
    ax[0].set_xlabel(r'drt in $[h^{-1}Mpc]$')
    ax[0].set_ylabel(r'drt out $[h^{-1}Mpc]$')
    ax[1].set_title(r'Residuals')
    ax[1].plot(drt_in, residuals_u, '-o', color=colors[2], label=r'Unlensed $\xi$') 
    ax[1].plot(drt_in, residuals_l, '-o', color=colors[3], label=r'Lensed $\xi$') 
    ax[1].axhline(0, color=colors[4], ls=':')
    ax[1].set_xlabel(r'drt in $[h^{-1}Mpc]$')
    ax[1].set_ylabel(r'drt out - drt in $[h^{-1}Mpc]$')
    ax[1].legend() 














