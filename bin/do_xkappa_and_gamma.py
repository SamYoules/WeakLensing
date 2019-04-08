## SY 17/12/18
## Calculate kappa for midpoint between quasars using estimator for lya x qso

#!/usr/bin/env python

import numpy as N
import pylab as P
import scipy as sp
import fitsio
from astropy.io import fits
import argparse
import glob
import healpy
import sys
from scipy import random
import copy

from picca import constants, xcf, io
from picca.data import delta

from multiprocessing import Pool,Process,Lock,Manager,cpu_count,Value

import configargparse



class kappa:

    nside = 256
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

    @staticmethod
    def load_model(modelfile, nbins=50) :

        data_rp, data_rt, xi_dist = N.loadtxt(modelfile, unpack=1)

        #-- get the larger value of the first separation bin to make a grid
        rp_min = data_rp.reshape(100, 50)[0].max()
        rp_max = data_rp.reshape(100, 50)[-1].min()
        rt_min = data_rt.reshape(100, 50)[:, 0].max()
        rt_max = data_rt.reshape(100, 50)[:, -1].min()
        #-- create the regular grid for griddata
        rp = N.linspace(rp_min, rp_max, nbins*2)
        rt = N.linspace(rt_min, rt_max, nbins)
        xim = sp.interpolate.griddata((data_rt, data_rp), xi_dist, \
                    (rt[:, None], rp[None, :]), method='cubic')

        #-- create interpolator object
        xi2d = sp.interpolate.RectBivariateSpline(rt, rp, xim)

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
                d.neighs = [q for q in neighs if d.ra > q.ra]
                #d.neighs = sp.array([q for q in neighs if (10**(d.ll[-1]- \
                #         sp.log10(lambda_abs))-1 + q.zqso)/2. >= z_cut_min \
                #         and (10**(d.ll[-1]- sp.log10(lambda_abs))-1 + \
                #         q.zqso)/2. < z_cut_max])

        
    @staticmethod
    def get_kappa(pixels):

        ikappa = []
        skappa = {} 
        wkappa = {}
        gamma1 = {}
        gamma2 = {}

        for ipix in pixels:
            for i,d in enumerate(kappa.data[ipix]):
                sys.stderr.write("\rcomputing kappa: {}%".format(\
                                 round(kappa.counter.value*100./kappa.ndata,2)))
                with kappa.lock:
                    kappa.counter.value += 1
                for q in d.neighs:
                    #--  compute the cartesian mid points and convert back 
                    #--  to ra, dec
                    mid_xcart = 0.5*(d.xcart+q.xcart)
                    mid_ycart = 0.5*(d.ycart+q.ycart)
                    mid_zcart = 0.5*(d.zcart+q.zcart)
                    mid_ra, mid_dec = get_radec(\
                        N.array([mid_xcart, mid_ycart, mid_zcart]))

                    #-- apply rotation into Galactic coordinates
                    #th, phi = kappa.rot(sp.pi/2-mid_dec, mid_ra)
                    #-- keeping without rotation
                    th, phi = sp.pi/2-mid_dec, mid_ra
                    th1, phi1 = sp.pi/2-d.dec, d.ra
                    thq, phiq = sp.pi/2-q.dec, q.ra

                    #-- check if pair of skewers belong to same spectro
                    same_half_plate = (d.plate == q.plate) and\
                            ( (d.fid<=500 and q.fid<=500) or \
                            (d.fid>500 and q.fid>500) )

                    #-- angle between skewers
                    ang = d^q
                    if kappa.true_corr:
                        ang_delensed = d.delensed_angle(q)

                    #-- getting pixel in between 
                    mid_pix = healpy.ang2pix(kappa.nside, \
                                  th, phi) 

                    if kappa.true_corr:
                        sk, wk, g1, g2 = kappa.fast_kappa_true(\
                                d.z, d.r_comov, \
                                q.zqso, q.r_comov, \
                                ang_delensed, ang, \
                                th1, phi1, thq, phiq)
                    else:
                        sk, wk, g1, g2 = kappa.fast_kappa(\
                                d.z, d.r_comov, d.we, d.de, \
                                q.zqso, q.r_comov, q.we, ang, \
                                th1, phi1, thq, phiq)

                    if mid_pix in ikappa:
                        skappa[mid_pix]+=sk
                        wkappa[mid_pix]+=wk
                        gamma1[mid_pix]+=g1
                        gamma2[mid_pix]+=g2
                    else:
                        ikappa.append(mid_pix)
                        skappa[mid_pix]=sk
                        wkappa[mid_pix]=wk
                        gamma1[mid_pix]=g1
                        gamma2[mid_pix]=g2

                setattr(d, "neighs", None)
                
        return ikappa, skappa, wkappa, gamma1, gamma2

    @staticmethod
    def fast_kappa(z1,r1,w1,d1,zq,rq,wq,ang, \
                   th1, phi1, thq, phiq):
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

        x = (phiq - phi1) * sp.sin(th1)
        y = thq - th1
        gamma_ang = sp.arctan2(y,x)

        #-- getting model and first derivative
        xi_model = kappa.xi2d(rt, rp, grid=False)
        xip_model = kappa.xi2d(rt, rp, dx=1, grid=False)

        #-- weight of estimator
        R = 1/(xip_model*rt)
       
        ska = sp.sum( (de - xi_model)/R*we )
        wka = sp.sum( we/R**2 ) 
        gam1 = 4*sp.cos(2*gamma_ang) * sp.sum( (de - xi_model)/R*we )
        gam2 = 4*sp.sin(2*gamma_ang) * sp.sum( (de - xi_model)/R*we )

        return ska, wka, gam1, gam2

    @staticmethod
    def fast_kappa_true(z1, r1, zq, rq, ang, ang_lens, \
                        th1, phi1, thq, phiq):
        
        rp      = abs(r1[:,None]-rq)*sp.cos(ang/2)
        rt      = (r1[:,None]+rq)*sp.sin(ang/2)
        rp_lens = abs(r1[:,None]-rq)*sp.cos(ang_lens/2)
        rt_lens = (r1[:,None]+rq)*sp.sin(ang_lens/2)
        
        w = (rp>=kappa.rp_min) & (rp<=kappa.rp_max) & \
            (rt<=kappa.rt_max) & (rt>=kappa.rt_min)

        rp = rp[w]
        rt = rt[w]
        rp_lens = rp_lens[w]
        rt_lens = rt_lens[w]

        x = (phiq - phi1) * sp.sin(th1)
        y = thq - th1
        gamma_ang = sp.arctan2(y,x)

        #-- getting model and first derivative
        xi_model  = kappa.xi2d(rt,      rp,       grid=False)
        xi_lens   = kappa.xi2d(rt_lens, rp_lens,  grid=False)
        xip_model = kappa.xi2d(rt,      rp, dx=1, grid=False)
        R = 1/(xip_model*rt)

        ska = sp.sum( (xi_lens - xi_model)/R )
        wka = sp.sum( 1/R**2  )
        gam1 = 4*sp.cos(2*gamma_ang) * sp.sum( (xi_lens - xi_model)/R )
        gam2 = 4*sp.sin(2*gamma_ang) * sp.sum( (xi_lens - xi_model)/R )

        return ska, wka, gam1, gam2


def get_radec(pos):
    ra = N.arctan(pos[1]/pos[0]) + N.pi + N.pi*(pos[0]>0)
    ra -= 2*N.pi*(ra>2*N.pi)
    dec = N.arcsin(pos[2]/N.sqrt(pos[0]**2+pos[1]**2+pos[2]**2))
    return ra, dec



def compute_kappa(p):
    tmp = kappa.get_kappa(p)
    return tmp
            

if __name__=='__main__':

    parser = configargparse.ArgParser()

    parser.add('--deltas', required=True, type=str, \
               help='folder containing deltas in pix format')
    parser.add('--model', required=True, \
               help='text file containing model')
    parser.add('--out', required=True, \
               help='output fits file with kappa map')
    parser.add('--gout', required=True, \
               help='output fits file with gamma maps')
    parser.add('--nproc', required=False, type=int, default=1, \
               help='number of procs used in calculation')
    parser.add('--nspec', required=False, type=int, default=None, \
               help='number of spectra to process')
    parser.add_argument('--drq', type=str, default=None, \
               required=True, help='Catalog of objects in DRQ format')
    parser.add_argument('--nside', type=int, default=256, required=False, \
               help='Healpix nside')
    parser.add_argument('--z-min-obj', type=float, default=None, \
               required=False, help='Min redshift for object field')
    parser.add_argument('--z-max-obj', type=float, default=None, \
               required=False, help='Max redshift for object field')
    parser.add_argument('--z-evol-obj', type=float, default=1., \
               required=False, help='Exponent of the redshift evolution of \
               the object field')
    parser.add_argument('--z-ref', type=float, default=2.25, required=False, \
               help='Reference redshift')
    parser.add_argument('--fid-Om', type=float, default=0.315, \
               required=False, help='Omega_matter(z=0) of fiducial LambdaCDM \
               cosmology')
    parser.add('--rt_min', required=False, type=float, default=3., \
               help='minimum transverse separation')
    parser.add('--rp_min', required=False, type=float, default=3., \
               help='minimum radial separation')
    parser.add('--rt_max', required=False, type=float, default=40., \
               help='maximum transverse separation')
    parser.add('--rp_max', required=False, type=float, default=10., \
               help='maximum radial separation')
    parser.add('--true_corr', required=False, default=False,\
               action='store_true', help='use actual lensed correlation')
    args, unknown = parser.parse_known_args()

    ### Read objects
    cosmo = constants.cosmo(args.fid_Om)
    objs,zmin_obj = io.read_objects(args.drq, kappa.nside_data, args.z_min_obj, args.z_max_obj,\
                                args.z_evol_obj, args.z_ref,cosmo)
    sys.stderr.write("\n")

    kappa.objs = objs
    kappa.true_corr = args.true_corr
    kappa.rt_min = args.rt_min
    kappa.rp_min = args.rp_min
    kappa.rt_max = args.rt_max
    kappa.rp_max = args.rp_max
    kappa.load_model(args.model)
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

    #-- compiling results from pool
    kap = sp.zeros(12*kappa.nside**2)
    wkap = sp.zeros(12*kappa.nside**2)
    ga1 = sp.zeros(12*kappa.nside**2)
    ga2 = sp.zeros(12*kappa.nside**2)
    wgam = sp.zeros(12*kappa.nside**2)
    for i, r in enumerate(results):
        print(i, len(results))
        index = N.array(r[0])
        for j in index:
            kap[j]  += r[1][j]
            wkap[j] += r[2][j]
            ga1[j]  += r[3][j]
            ga2[j]  += r[4][j]
            wgam[j] += r[2][j]*2. # SY: factor of 2 comes from squaring the
			          #  prefactors of r in the gamma estimator
    
    w = wkap>0
    kap[w]/= wkap[w]
    ga1[w]/= wgam[w]
    ga2[w]/= wgam[w]
    
    out = fitsio.FITS(args.out, 'rw', clobber=True)
    head = {}
    head['RPMIN']=kappa.rp_min
    head['RPMAX']=kappa.rp_max
    head['RTMIN']=kappa.rt_min
    head['RTMAX']=kappa.rt_max
    head['NT']=kappa.nt
    head['NP']=kappa.np
    head['NSIDE']=kappa.nside
    out.write([kap, wkap], names=['kappa', 'wkappa'], header=head)
    out.close()

    # SY: Write gamma estimator fits file
    out = fitsio.FITS(args.gout, 'rw', clobber=True)
    head = {}
    head['RPMIN']=kappa.rp_min
    head['RPMAX']=kappa.rp_max
    head['RTMIN']=kappa.rt_min
    head['RTMAX']=kappa.rt_max
    head['NT']=kappa.nt
    head['NP']=kappa.np
    head['NSIDE']=kappa.nside
    out.write([ga1, ga2, wgam], names=['gamma1', 'gamma2', 'wgamma'], header=head)
    out.close()

