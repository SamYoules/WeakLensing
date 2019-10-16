#!/usr/bin/env python

#-- plot 2d residuals (kappa_in - kappa_out)/std as function of rt, rp

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

    return xi2d, rt, rp

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Improve fit of xcf for lya-lya.')

    parser.add_argument('--kappa_in', required=True, type=float, \
               help='kappa value to model')
    parser.add_argument('--nreal', required=True, type=int, default=1000000, \
               help='number of delta pair realisations')
    parser.add_argument('--mcreal', required=True, type=int, default=100, \
               help='number of mc realisations')
    parser.add_argument('--var', required=True, type=float, default=0.005, \
               help='variance of deltas (approx 0.04 at z=2 to 0.15 at z=3)')
    parser.add_argument('--xi', required=True, \
               help='text file containing model')
    parser.add_argument('--fit', required=True, \
               help='text file containing corrrelation function')

    args, unknown = parser.parse_known_args()

    kappa_in = args.kappa_in
    nrealisations = args.nreal
    mc_realisations = args.mcreal
    var1 = args.var
    var2 = args.var
    xi_file = args.xi
    fit_file = args.fit

    # Get xi_model interpolator object
    xi2d, rt, rp = load_model(xi_file, fit_file)

    residuals = np.zeros(shape=(50,50))
    for i, r_t in enumerate(rt):
        for j, r_p in enumerate(rp):
            kappas_out = []
            for r in range(mc_realisations):

                # Construct two columns of gaussian random variables
                g = np.random.randn(2, nrealisations)

                # Construct Covariance matrix
                xi_un = xi2d(r_t, r_p, grid=False)
                C = np.array([ [var1, xi_un  ],
                               [xi_un,   var2] ])

                # Cholesky decomposition, see 2.1 in 1108.5606
                L = np.linalg.cholesky(C)

                # Apply transformation to correlate variables
                delta = L.dot(g)

                # Apply lensing (approximation ignoring shear)
                rt_lens = r_t * (1 - kappa_in)

                # Estimate kappa
                xi_model = xi2d(rt_lens, r_p)
                R = 1 / rt_lens / xi2d(rt_lens, r_p, dx=1)
                kappa_out = np.mean(  (delta[0]*delta[1] - xi_model)*R  )
                kappas_out.append(kappa_out)

            kappas_out = np.array(kappas_out)
            kappa_mean = np.mean(kappas_out)
            kappa_std = np.std(kappas_out)
            residuals[i,j] = (kappa_mean - kappa_in)/kappa_std

    # Save mean kappa and error to file
    t = Table([kappa_mean, kappa_std], names=('MEAN', 'STD'))
    t.write(f'monte_carlo/mc_kappa_{kappa_in}.fits', format='fits', overwrite=True)

    # Plot 2D figure of residuals
    plt.figure()
    plt.pcolormesh(residuals, vmin=-2, vmax=2)
    plt.xlabel('rt')
    plt.ylabel('rp')
    plt.title(f'Auto kappa residuals, kappa = {kappa_in}')
    plt.colorbar()
    plt.savefig(f'kappa_{kappa_in}.png')

















