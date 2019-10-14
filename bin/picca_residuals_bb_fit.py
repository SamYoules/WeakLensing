#!/usr/bin/env python

#-- Plot residuals from cf fit with added broadband terms and chi2

import picca.wedgize
import os
import h5py
import fitsio
import numpy as np
import pylab as plt
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--xi', required=True, help='Input correlation function')
parser.add_argument('--fit', required=False, default=None, help='Input HDF5 file with best fit')
parser.add_argument('--scale_r', default=2, type=int, help='Scale xi by r**scale_r')
parser.add_argument('--mus', nargs='+', help='Array with mus defining wedges (for negative values add quotes and a space before the number e.g. " -1"')
parser.add_argument('--savefig', help='Name of figure to be saved in PDF format')

args = parser.parse_args()
if args.mus is None:
    args.mus = [1., 0.95, 0.8, 0.5, 0]
else:
    args.mus = np.array(args.mus).astype(float)

def plot_bestfit(file_xi, file_fit=None, power=2, mus=[1., 0.95, 0.8, 0.5, 0], 
                 figsize=(6, 8), label=None):

    #- Read correlation function and covariance matrix
    h = fitsio.FITS(file_xi)
    da = h[1]['DA'][:]
    co = h[1]['CO'][:]
    hh = h[1].read_header()
    rpmin = hh['RPMIN']
    rpmax = hh['RPMAX']
    rtmin = 0 
    rtmax = hh['RTMAX']
    nrp = hh['NP']
    nrt = hh['NT']
    h.close()

    #-- Read fit h5 file if there's one
    if file_fit:
        ff = h5py.File(file_fit, 'r')
        keys = list(ff.keys())
        for k in keys:
            if k != 'best fit':
                base = k
                break
        fit = ff[base+'/fit'][...] 
        attr = dict(ff['best fit'].attrs)
        chi2 = attr['fval']
        ndata = attr['ndata']
        npars = attr['npar']
        rchi2 = chi2/(ndata-npars)
        print(f'chi2/(ndata-npars) = {chi2}/({ndata}-{npars}) = {rchi2}')
        ff.close()

    for i, (mumax,mumin) in enumerate(zip(mus[:-1],mus[1:])):
        b = picca.wedgize.wedge(mumin=mumin, mumax=mumax, 
                                rpmin=rpmin, rpmax=rpmax, rtmin=rtmin, rtmax=rtmax, nrt=nrt, 
                                nrp=nrp,absoluteMu=False, 
                                rmin=0., rmax=min(rpmax, rtmax), nr=min(nrt, nrp))
        r,d,c = b.wedge(da,co)

        w = r<=100  # sy
        r = r[w]    # sy
        d = d[w]    # sy
        c = c[w]    # sy

        nrows = 1 
        fig, ax = plt.subplots(nrows=nrows, ncols=1, figsize=figsize)
        fig.subplots_adjust(hspace=0)

        #-- Best model and Residuals
        y = d*r**power
        dy = np.sqrt(c.diagonal())*r**power
        if file_fit:
            r, model, _ = b.wedge(fit, co)
            w = r<=100  # sy
            r = r[w]    # sy
            model = model[w]    # sy

            y = (d-model)*r**power
            ax[0].errorbar(r, y, dy, fmt='o')
            ax[0].axhline(0, color='k', ls=':')
            ax[0].set_ylabel(r"$r^{power} [\xi-\xi_{{\rm model}}](r)$".format(power=power))
        ax[0].set_xlabel(r"$r \, [h^{-1}\, \mathrm{Mpc}]$")
        plt.tight_layout()
        if args.savefig:
            plt.savefig(args.savefig.format(i=i))
    plt.show()

plot_bestfit(args.xi, file_fit=args.fit, power=args.scale_r, mus=args.mus)

