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
    parser.add_argument('--xi', required=True, \
               help='text file containing model')
    parser.add_argument('--fit', required=True, \
               help='text file containing corrrelation function')

    args, unknown = parser.parse_known_args()

    rp_in = args.rp
    rt_in = args.rt
    xi_file = args.xi
    fit_file = args.fit

    # Make range of drt to plot
    drt_in = [-0.2, -0.14, -0.08, -0.02, 0.04, 0.1, 0.16, 0.22]

    # Get xi_model interpolator object
    print('Reading model from ', fit_file)
    xi2d = load_model(xi_file, fit_file)

    drt_out = []
    for drt in drt_in:
        # Apply lensing
        rt_lens = rt_in + drt

        # Estimate drt_out
        xi_unlensed = xi2d(rt_in, rp_in)
        xi_lensed = xi2d(rt_lens, rp_in)
        drt_out.append((xi_lensed - xi_unlensed) / xi2d(rt_in, rp_in, dx=1))

    # Plot 2D figure of drt
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(8,10))
    fig.subplots_adjust(hspace=0.5)
    ax[0].plot(drt_in, ravel(drt_out), 'o')
    ax[0].plot(drt_in, drt_in, color='g', ls=':')
    ax[0].set_xlabel(r'drt in $[h^{-1}Mpc]$')
    ax[0].set_ylabel(r'drt out $[h^{-1}Mpc]$')
    ax[0].set_title(fr'$\Delta r$, 1 realisation, rt={rt_in}, rp={rp_in}')

    # Plot residual drt
    ax[1].plot(drt_in, ravel(drt_out) - drt_in, 'o')
    ax[1].set_title(r'Residuals')
    ax[1].axhline(0, color='g', ls=':')
    ax[1].set_xlabel(r'drt in $[h^{-1}Mpc]$')
    ax[1].set_ylabel(r'drt out - drt in $[h^{-1}Mpc]$')

    # Plot delta r over r
    ax[2].plot(drt_in, ravel(drt_out) / rt_in, '-bo', label=r'$drt / rt_{Unlensed}$')
    ax[2].plot(drt_in, ravel(drt_out) / rt_lens, '--ro', label=r'$drt / rt_{Lensed}$')
    ax[2].set_title(r'$\Delta r / rt$')
    ax[2].axhline(0, color='g', ls=':')
    ax[2].set_xlabel(r'drt in $[h^{-1}Mpc]$')
    ax[2].set_ylabel(r'drt / rt $[h^{-1}Mpc]$')
    ax[2].legend()
    plt.savefig(f'plots/mc/delta_r_1.png')

