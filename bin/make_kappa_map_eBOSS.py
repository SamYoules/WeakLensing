#!/usr/bin/env python

## Read in kappalist file (qso pairwise kappa) and produce midpoint-estimator kappa file

from astropy.table import Table, join
import numpy as np
import healpy as hp
import argparse
import configargparse
import fitsio

def get_radec(thid_in, ra_in, dec_in, thid_out):
    ''' Function to obtain RA, DEC from THING_ID values'''

    win = np.argsort(thid_in)
    wout = np.searchsorted(thid_in[win], thid_out)
    return ra_in[win][wout], dec_in[win][wout]

def radec_to_cart(ra, dec):

    ra_r = np.radians(ra)
    dec_r = np.radians(dec)
    x = np.cos(dec_r)*np.cos(ra_r)
    y = np.cos(dec_r)*np.sin(ra_r)
    z = np.sin(dec_r)
    return x, y, z

def cart_to_radec(x, y, z):

    dist = np.sqrt(x**2+y**2+z**2)
    dec = 90 - np.degrees(np.arccos(z / dist))
    ra = np.degrees(np.arctan2(y, x))
    ra[ra < 0] += 360
    return ra, dec

def make_map(nside, ra1, dec1, ra2, dec2, skappa, wkappa):

    x1, y1, z1 = radec_to_cart(ra1, dec1)
    x2, y2, z2 = radec_to_cart(ra2, dec2)

    x_mid = 0.5*(x2+x1)
    y_mid = 0.5*(y2+y1)
    z_mid = 0.5*(z2+z1)

    ra_mid, dec_mid = cart_to_radec(x_mid, y_mid, z_mid)

    th_mid, phi_mid = np.pi/2-np.radians(dec_mid), np.radians(ra_mid)
    pix_mid = hp.ang2pix(nside, th_mid, phi_mid)


    smap = np.bincount(pix_mid, weights=skappa, 
                    minlength=hp.nside2npix(nside))
    wmap = np.bincount(pix_mid, weights=wkappa, 
                    minlength=hp.nside2npix(nside))

    kappa_map = smap*0
    w = wmap !=0
    kappa_map[w] = smap[w]/wmap[w]

    return kappa_map, wmap


if __name__=='__main__':

    parser = configargparse.ArgParser()

    parser.add('--mapin', required=True, type=str, \
               help='pathname to kappalist fits file')
    parser.add('--mapout', required=True, type=str, \
               help='output fits file with kappa values') 
    parser.add('--rt_min', required=False, type=float, default=3., \
               help='minimum transverse separation')
    parser.add('--rp_min', required=False, type=float, default=3., \
               help='minimum radial separation')
    parser.add('--rt_max', required=True, type=float, default=40., \
               help='maximum transverse separation')
    parser.add('--rp_max', required=True, type=float, default=10., \
               help='maximum radial separation')
    parser.add('--nside', required=False, type=int, default=256, \
               help='resolution of map')
    parser.add('--zcat', required=False, type=str, \
               help='catalogue file')
    args, unknown = parser.parse_known_args()

    mapin  = args.mapin
    mapout = args.mapout
    rt_min = args.rt_min
    rp_min = args.rp_min
    rt_max = args.rt_max
    rp_max = args.rp_max
    nside  = args.nside
    zcat   = args.zcat
 
    #-- read list of kappas
    kap = Table.read(mapin)

    #-- read QSO catalog with THING_ID, RA and DEC information 
    drq = Table.read(zcat)
    drq.remove_columns(['Z', 'PLATE', 'MJD', 'FIBERID'])

    ra1, dec1 = get_radec(drq['THING_ID'].data, drq['RA'].data, drq['DEC'].data,
                      kap['THID1'].data)
    ra2, dec2 = get_radec(drq['THING_ID'].data, drq['RA'].data, drq['DEC'].data,
                      kap['THID2'].data)
    skappa = kap['SKAPPA'].data
    wkappa = kap['WKAPPA'].data

    kappa_map, wmap = make_map(nside, ra1, dec1, ra2, dec2, skappa, wkappa)

    ##- Write map
    out = fitsio.FITS(mapout, 'rw', clobber=True)
    head = {}
    head['RPMIN']=rp_min
    head['RPMAX']=rp_max
    head['RTMIN']=rt_min
    head['RTMAX']=rt_max
    head['NT']=1
    head['NP']=1
    head['NSIDE']=nside
    out.write([kappa_map, wmap], names=['kappa', 'wkappa'], header=head)
    out.close()


