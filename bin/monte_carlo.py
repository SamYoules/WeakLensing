#!/usr/bin/env python

import numpy as np
from numpy import linalg
import scipy as sp
from scipy import interpolate
from picca.data import forest, delta
import fitsio
import h5py
import argparse
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
        description='Improve fit of xcf for qso-lya.')

    parser.add_argument('--rt', required=True, type=float, \
               help='r transverse separation between qso and delta')
    parser.add_argument('--nreal', required=True, type=int, default=1000000, \
               help='number of realisations')
    parser.add_argument('--var', required=True, type=float, default=0.005, \
               help='variance of deltas (approx 0.04 at z=2 to 0.15 at z=3)')
    parser.add_argument('--xi', required=True, \
               help='text file containing model')
    parser.add_argument('--fit', required=True, \
               help='text file containing corrrelation function')

    args, unknown = parser.parse_known_args()

    rt = args.rt
    nrealisations = args.nreal
    var1 = args.var
    var2 = args.var
    xi_file = args.xi
    fit_file = args.fit

    # Construct two columns of gaussian random variables
    g = np.random.randn(2, nrealisations)

    # Get xi_model interpolator object
    xi2d = load_model(xi_file, fit_file)

    # Set range of kappa and rp to test
    kappa_in = [-0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0.,
                 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]

    # kappa range for the blob
    #kappa_in = [0.0001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001]
    r_p = [10,20,30,40,50,60,70]

    # Open template file
    t = Table.read('monte_carlo/monte_carlo_template.fits')
    for rp in r_p:
        # Construct Covariance matrix
        xi_un = xi2d(rp,rt, grid=False)
        C = np.array([ [var1, xi_un  ],
                       [xi_un,   var2] ])

        # Cholesky decomposition, see 2.1 in 1108.5606
        L = np.linalg.cholesky(C)

        # Apply transformation to correlate variables
        delta = L.dot(g)

        # Check covariance matrix of delta is similar to C + noise
        print(np.cov(delta)) 

        for k in kappa_in:

            # Apply lensing (approximation ignoring shear)
            rt_lens = rt * (1 - k)

            # Estimate kappa
            xi_model = xi2d(rt_lens, rp)
            R = - 1 / rt_lens / xi2d(rt_lens, rp, dx=1)
            kappa_out = np.mean(  (delta[0]*delta[1] - xi_model)*R  )

            # Write output
            t.add_row([k, rp, rt, var1, var2, nrealisations, kappa_out])
            
    t.remove_row(0)
    t.write('monte_carlo/monte_carlo_results{}.fits'.format(int(rt)), format='fits', overwrite=True)


















