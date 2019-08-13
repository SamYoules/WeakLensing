#!/usr/bin/env python

## 17/12/18
## Calculate kappa for midpoint between quasars using estimator for lya x qso

import numpy as np
import pylab as P
import scipy as sp
import fitsio
from astropy.io import fits
import argparse
import glob
import healpy
import sys
from scipy import random
from scipy import interpolate ## SY 7/6/19
import copy
from picca import constants, xcf, io_lens
from picca.data_lens import delta
from multiprocessing import Pool,Process,Lock,Manager,cpu_count,Value
import h5py ## SY 5/6/19


class kappa:

    #nside = 256
    nside_data = 32
    rot = healpy.Rotator(coord=['C', 'G'])
    lambda_abs = 1215.67

    fid_Om=0.31
    cosmo = constants.cosmo(fid_Om)

    rt_min=0.
    rt_max=40.
    rp_min=0.
    rp_max=10.
    nt = 1
    np = 1

    xi2d = None

    angmax = None
    z_min_pix = None

    counter = Value('i',0)
    lock = Lock()

    data={}
    ndata=0


    ## - SY 4/6/19 Changed to use new picca model and fit files
    @staticmethod
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

        rpmin = data_rp.reshape(100, 50)[0].max()
        rpmax = data_rp.reshape(100, 50)[-1].min()
        rtmin = data_rt.reshape(100, 50)[:, 0].max()
        rtmax = data_rt.reshape(100, 50)[:, -1].min()

        #-- create the regular grid for griddata
        rp = np.linspace(rpmin, rpmax, nbins*2)
        rt = np.linspace(rtmin, rtmax, nbins)

        #xim = sp.interpolate.griddata((data_rt, data_rp), fit, \
        #            (rt[:, None], rp[None, :]), method='cubic')

        xim = sp.interpolate.griddata((data_rt, data_rp), fit,
                    (np.outer(np.ones(rp.size), rt).ravel(),
                    np.outer(rp, np.ones(rt.size)).ravel()), method='cubic')

        #-- create interpolator object
        xi2d = sp.interpolate.RectBivariateSpline(rt, rp, \
                   xim.reshape((100, 50)).T )

        kappa.xi2d = xi2d
        return xi2d

    @staticmethod
    def read_deltas(in_dir, nspec=None):
        data = {}
        ndata = 0
        dels = []

        fi = glob.glob(in_dir+"/*.fits.gz")
        for i,f in enumerate(fi):
            sys.stderr.write("\rread {} of {} {}".format(i,len(fi),ndata))
            hdus = fitsio.FITS(f)
            dels += [delta.from_fitsio(h) for h in hdus[1:]]
            ndata+=len(hdus[1:])
            hdus.close()
            if nspec and ndata > nspec:
                break

        phi = [d.ra for d in dels]
        th = [sp.pi/2.-d.dec for d in dels]
        pix = healpy.ang2pix(kappa.nside_data, th, phi)

        z_min_pix = 10**dels[0].ll[0]/kappa.lambda_abs-1

        for d,p in zip(dels,pix):
            if p not in data:
                data[p]=[]
            data[p].append(d)

            z = 10**d.ll/kappa.lambda_abs-1.
            z_min_pix = sp.amin( sp.append([z_min_pix], z) ) 
            d.z = z
            d.r_comov = kappa.cosmo.r_comoving(z)
            d.we *= ((1.+z)/(1.+2.25))**(2.9-1.)
            d.project()

        kappa.z_min_pix = z_min_pix
        kappa.angmax = 2.*\
                sp.arcsin(kappa.rt_max/(2.*kappa.cosmo.r_comoving(z_min_pix)))
        kappa.ndata = ndata
        kappa.data = data

        return data

    @staticmethod
    def fill_neighs():
        data = kappa.data
        objs = kappa.objs
        print('\n Filling neighbors')
        for ipix in data.keys():
            for d in data[ipix]:
                npix = healpy.query_disc(kappa.nside_data, \
                        [d.xcart, d.ycart, d.zcart], \
                        kappa.angmax, inclusive = True)
                npix = [p for p in npix if p in objs]
                neighs = [q for p in npix for q in objs[p] if q.thid != d.thid]
                ang = d^neighs
                w = ang<kappa.angmax
                neighs = sp.array(neighs)[w]
                d.qneighs = [q for q in neighs]

        
    @staticmethod
    def get_kappa(pixels):

        id1 = []
        id2 = []
        skappa = [] 
        wkappa = []

        for ipix in pixels:
            for i,d in enumerate(kappa.data[ipix]):
                sys.stderr.write("\rcomputing kappa: {}%".format(\
                                 round(kappa.counter.value*100./kappa.ndata,2)))
                with kappa.lock:
                    kappa.counter.value += 1
                for q in d.qneighs:

                    #-- angle between skewers
                    ang = d^q
                    if kappa.true_corr:
                        ang_delensed = d.delensed_angle(q)
 
                    if kappa.true_corr:
                        sk, wk = kappa.fast_kappa_true(\
                                d.z, d.r_comov, \
                                q.zqso, q.r_comov, \
                                ang_delensed, ang)
                    else:
                        sk, wk = kappa.fast_kappa(\
                                d.z, d.r_comov, d.we, d.de, \
                                q.zqso, q.r_comov, q.we, ang)

                    if wk != 0:
                        id1.append(q.thid)
                        id2.append(d.thid)
                        skappa.append(sk)
                        wkappa.append(wk)

                setattr(d, "neighs", None)
                
        return id1, id2, skappa, wkappa

    @staticmethod
    def fast_kappa(z1,r1,w1,d1,zq,rq,wq,ang):
        rp = (r1[:,None]-rq)*sp.cos(ang/2)
        rt = (r1[:,None]+rq)*sp.sin(ang/2)
        
        we = w1[:,None]*wq
        de = d1[:,None]

        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & \
            (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        we = we[w]
        de = de[w]
 
        #-- getting model and first derivative
        xi_model = kappa.xi2d(rt, rp, grid=False)
        xip_model = kappa.xi2d(rt, rp, dx=1, grid=False)

        #-- weight of estimator
        R = -1/(xip_model*rt)

        ska = sp.sum( (de - xi_model)/R*we )
        wka = sp.sum( we/R**2 ) 

        return ska, wka

    @staticmethod
    def fast_kappa_true(z1, r1, zq, rq, ang, ang_lens):
        
        rp      = (r1[:,None]-rq)*sp.cos(ang/2)
        rt      = (r1[:,None]+rq)*sp.sin(ang/2)
        rp_lens = (r1[:,None]-rq)*sp.cos(ang_lens/2)
        rt_lens = (r1[:,None]+rq)*sp.sin(ang_lens/2)
        
        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & \
            (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        rp_lens = rp_lens[w]
        rt_lens = rt_lens[w]

        #-- getting model and first derivative
        xi_model  = kappa.xi2d(rt,      rp,       grid=False)
        xi_lens   = kappa.xi2d(rt_lens, rp_lens,  grid=False)
        xip_model = kappa.xi2d(rt,      rp, dx=1, grid=False)
        R = -1/(xip_model*rt)

        ska = sp.sum( (xi_lens - xi_model)/R )
        wka = sp.sum( 1/R**2  )

        return ska, wka


def compute_kappa(p):
    id1, id2, skappa, wkappa = kappa.get_kappa(p)
    return id1, id2, skappa, wkappa
            

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute the convergence (kappa) between qsos and deltas.')

    parser.add_argument('--deltas', required=True, type=str, \
               help='folder containing deltas in pix format')
    parser.add_argument('--xi', required=True, \
               help='text file containing model') # SY 5/6/19
    parser.add_argument('--fit', required=True, \
               help='text file containing fit')  # SY 5/6/19
    parser.add_argument('--out', required=True, \
               help='output fits file with kappa map')
    parser.add_argument('--nproc', required=False, type=int, default=1, \
               help='number of procs used in calculation')
    parser.add_argument('--nspec', required=False, type=int, default=None, \
               help='number of spectra to process')
    parser.add_argument('--drq', type=str, default=None, \
               required=True, help='Catalog of objects in DRQ format')
    parser.add_argument('--z_min_obj', type=float, default=None, \
               required=False, help='Min redshift for object field')
    parser.add_argument('--z_max_obj', type=float, default=None, \
               required=False, help='Max redshift for object field')
    parser.add_argument('--z_cut_min', type=float, default=0., required=False,
        help='Use only pairs of forest x object with the mean of the last absorber \
        redshift and the object redshift larger than z_cut_min') #SY 13/2/19
    parser.add_argument('--z_cut_max', type=float, default=10., required=False,
        help='Use only pairs of forest x object with the mean of the last absorber \
        redshift and the object redshift smaller than z_cut_max') #SY 13/2/19
    parser.add_argument('--z_evol_obj', type=float, default=1., \
               required=False, help='Exponent of the redshift evolution of \
               the object field')
    parser.add_argument('--z_ref', type=float, default=2.25, required=False, \
               help='Reference redshift')
    parser.add_argument('--fid_Om', type=float, default=0.315, \
               required=False, help='Omega_matter(z=0) of fiducial LambdaCDM \
               cosmology')
    parser.add_argument('--rt_min', required=False, type=float, default=3., \
               help='minimum transverse separation')
    parser.add_argument('--rp_min', required=False, type=float, default=3., \
               help='minimum radial separation')
    parser.add_argument('--rt_max', required=False, type=float, default=40., \
               help='maximum transverse separation')
    parser.add_argument('--rp_max', required=False, type=float, default=10., \
               help='maximum radial separation')
    parser.add_argument('--true_corr', required=False, default=False,\
               action='store_true', help='use actual lensed correlation')
    args, unknown = parser.parse_known_args()

    ### Read objects
    cosmo = constants.cosmo(args.fid_Om)
    objs,zmin_obj = io_lens.read_objects(args.drq, kappa.nside_data, \
                                   args.z_min_obj, args.z_max_obj,\
                                   args.z_evol_obj, args.z_ref,cosmo)
    sys.stderr.write("\n")

    kappa.objs = objs
    kappa.true_corr = args.true_corr
    kappa.rt_min = args.rt_min
    kappa.rp_min = args.rp_min
    kappa.rt_max = args.rt_max
    kappa.rp_max = args.rp_max
    kappa.z_cut_max = args.z_cut_max #SY 13/2/19
    kappa.z_cut_min = args.z_cut_min #SY 13/2/19
    kappa.load_model(args.xi, args.fit) # SY 5/6/19
    kappa.read_deltas(args.deltas, nspec=args.nspec)
    kappa.fill_neighs()

    cpu_data = {}
    for p in kappa.data.keys():
        cpu_data[p] = [p]

    print(' ', len(kappa.data.keys()), 'pixels with data')
    pool = Pool(processes=args.nproc)
    results = pool.map(compute_kappa, cpu_data.values())
    pool.close()
    print('')

    id1 = np.empty(0, dtype=int)
    id2 = np.empty(0, dtype=int)
    skappa = np.empty(0)
    wkappa = np.empty(0)
    
    for r in results:
        id1 = np.append(id1, np.array(r[0]).astype(int))
        id2 = np.append(id2, np.array(r[1]).astype(int))
        skappa = np.append(skappa, np.array(r[2]))
        wkappa = np.append(wkappa, np.array(r[3]))

    
    out = fitsio.FITS(args.out, 'rw', clobber=True)
    head = {}
    head['RPMIN']=kappa.rp_min
    head['RPMAX']=kappa.rp_max
    head['RTMIN']=kappa.rt_min
    head['RTMAX']=kappa.rt_max
    head['NT']=kappa.nt
    head['NP']=kappa.np
    out.write([id1, id2, skappa, wkappa],
              names=['THIDQ', 'THID', 'SKAPPA', 'WKAPPA'],
              header=head)
    out.close()

