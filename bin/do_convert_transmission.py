#!/usr/bin/env python

import argparse
import picca
from picca import utils
import fitsio
import scipy as sp
import glob
import healpy
import sys
from picca.data import delta

def desi_convert_transmission_to_delta_files(zcat,outdir,indir=None,infiles=None,lObs_min=3600.,lObs_max=5500.,lRF_min=1040.,lRF_max=1200.,dll=3.e-4,nspec=None):
    """Convert desi transmission files to picca delta files

    Args:
        zcat (str): path to the catalog of object to extract the transmission from
        indir (str): path to transmission files directory
        outdir (str): path to write delta files directory
        lObs_min (float) = 3600.: min observed wavelength in Angstrom
        lObs_max (float) = 5500.: max observed wavelength in Angstrom
        lRF_min (float) = 1040.: min Rest Frame wavelength in Angstrom
        lRF_max (float) = 1200.: max Rest Frame wavelength in Angstrom
        dll (float) = 3.e-4: size of the bins in log lambda
        nspec (int) = None: number of spectra, if 'None' use all

    Returns:
        None

    """

    ### Catalog of objects
    h = fitsio.FITS(zcat)
    key_val = sp.char.strip(sp.array([ h[1].read_header()[k] for k in h[1].read_header().keys()]).astype(str))
    if 'TARGETID' in key_val:
        zcat_thid = h[1]['TARGETID'][:]
    elif 'THING_ID' in key_val:
        zcat_thid = h[1]['THING_ID'][:]
    w = h[1]['Z'][:]>max(0.,lObs_min/lRF_max -1.)
    w &= h[1]['Z'][:]<max(0.,lObs_max/lRF_min -1.)
    zcat_ra = h[1]['RA'][:][w].astype('float64')*sp.pi/180.
    zcat_dec = h[1]['DEC'][:][w].astype('float64')*sp.pi/180.
    zcat_thid = zcat_thid[w]
    h.close()

    ### List of transmission files
    if (indir is None and infiles is None) or (indir is not None and infiles is not None):
        print("ERROR: No transmisson input files or both 'indir' and 'infiles' given")
        sys.exit()
    elif indir is not None:
        fi = glob.glob(indir+'/*/*/transmission*.fits') + glob.glob(indir+'/*/*/transmission*.fits.gz')
        h = fitsio.FITS(sp.sort(sp.array(fi))[0])
        #in_nside = h[1].read_header()['NSIDE'] # SY
        in_nside = 16 #SY
        nest = True
        h.close()
        in_pixs = healpy.ang2pix(in_nside, sp.pi/2.-zcat_dec, zcat_ra, nest=nest)
        fi = sp.sort(sp.array([ indir+'/'+str(int(f/100))+'/'+str(f)+'/transmission-'+str(in_nside)+'-'+str(f)+'.fits' for f in sp.unique(in_pixs)]))
    else:
        fi = sp.sort(sp.array(infiles))

    ### Stack the transmission
    lmin = sp.log10(lObs_min)
    lmax = sp.log10(lObs_max)
    nstack = int((lmax-lmin)/dll)+1
    T_stack = sp.zeros(nstack)
    n_stack = sp.zeros(nstack)

    deltas = {}

    ### Read
    for nf, f in enumerate(fi):
        sys.stderr.write("\rread {} of {} {}".format(nf,fi.size,sp.sum([ len(deltas[p]) for p in list(deltas.keys())])))
        h = fitsio.FITS(f)
        thid = h[1]['MOCKID'][:]
        if sp.in1d(thid,zcat_thid).sum()==0:
            h.close()
            continue
        ra = h[1]['RA'][:].astype('float64')*sp.pi/180.
        dec = h[1]['DEC'][:].astype('float64')*sp.pi/180.
        z = h[1]['Z'][:]
        ll = sp.log10(h[2].read())
        trans = h[3].read()
        nObj = z.size
        pixnum = f.split('-')[-1].split('.')[0]

        if trans.shape[0]!=nObj:
            trans = trans.transpose()

        bins = sp.floor((ll-lmin)/dll+0.5).astype(int)
        tll = lmin + bins*dll
        lObs = (10**tll)*sp.ones(nObj)[:,None]
        lRF = (10**tll)/(1.+z[:,None])
        w = sp.zeros_like(trans).astype(int)
        w[ (lObs>=lObs_min) & (lObs<lObs_max) & (lRF>lRF_min) & (lRF<lRF_max) ] = 1
        nbPixel = sp.sum(w,axis=1)
        cut = nbPixel>=50
        cut &= sp.in1d(thid,zcat_thid)
        if cut.sum()==0:
            h.close()
            continue

        ra = ra[cut]
        dec = dec[cut]
        z = z[cut]
        thid = thid[cut]
        trans = trans[cut,:]
        w = w[cut,:]
        nObj = z.size
        h.close()

        deltas[pixnum] = []
        for i in range(nObj):
            tll = ll[w[i,:]>0]
            ttrans = trans[i,:][w[i,:]>0]

            bins = sp.floor((tll-lmin)/dll+0.5).astype(int)
            cll = lmin + sp.arange(nstack)*dll
            cfl = sp.bincount(bins,weights=ttrans,minlength=nstack)
            civ = sp.bincount(bins,minlength=nstack).astype(float)

            ww = civ>0.
            if ww.sum()<50: continue
            T_stack += cfl
            n_stack += civ
            cll = cll[ww]
            cfl = cfl[ww]/civ[ww]
            civ = civ[ww]
            deltas[pixnum].append(delta(thid[i],ra[i],dec[i],z[i],thid[i],thid[i],thid[i],cll,civ,None,cfl,1,None,None,None,None,None,None))
        if not nspec is None and sp.sum([ len(deltas[p]) for p in list(deltas.keys())])>=nspec: break

    print('\n')

    ### Get stacked transmission
    w = n_stack>0.
    T_stack[w] /= n_stack[w]

    ### Transform transmission to delta and store it
    for nf, p in enumerate(sorted(list(deltas.keys()))):
        sys.stderr.write("\rwrite {} of {} ".format(nf,len(list(deltas.keys()))))
        out = fitsio.FITS(outdir+'/delta-{}'.format(p)+'.fits.gz','rw',clobber=True)
        for d in deltas[p]:
            bins = sp.floor((d.ll-lmin)/dll+0.5).astype(int)
            d.de = d.de/T_stack[bins] - 1.
            d.we *= T_stack[bins]**2

            hd = {}
            hd['RA'] = d.ra
            hd['DEC'] = d.dec
            hd['Z'] = d.zqso
            hd['PMF'] = '{}-{}-{}'.format(d.plate,d.mjd,d.fid)
            hd['THING_ID'] = d.thid
            hd['PLATE'] = d.plate
            hd['MJD'] = d.mjd
            hd['FIBERID'] = d.fid
            hd['ORDER'] = d.order

            cols = [d.ll,d.de,d.we,sp.ones(d.ll.size)]
            names = ['LOGLAM','DELTA','WEIGHT','CONT']
            out.write(cols,names=names,header=hd,extname=str(d.thid))
        out.close()

    print('\n')

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Script to convert noiseless transmission files to delta picca files')

    parser.add_argument('--object-cat', type = str, default = None, required=True,
            help = 'Path to a catalog of objects to get the transmission from')

    parser.add_argument('--in-dir',type = str,default=None,required=False,
            help='Desi formated data directory to transmission files')

    parser.add_argument('--in-files',type = str,default=None,required=False,
            help='List of transmission files.', nargs='*')

    parser.add_argument('--out-dir',type = str,default=None,required=True,
            help='Output directory')

    parser.add_argument('--lambda-min',type = float,default=3600.,required=False,
            help='Lower limit on observed wavelength [Angstrom]')

    parser.add_argument('--lambda-max',type = float,default=5500.,required=False,
            help='Upper limit on observed wavelength [Angstrom]')

    parser.add_argument('--lambda-rest-min',type = float,default=1040.,required=False,
            help='Lower limit on rest frame wavelength [Angstrom]')

    parser.add_argument('--lambda-rest-max',type = float,default=1200.,required=False,
            help='Upper limit on rest frame wavelength [Angstrom]')

    parser.add_argument('--dll',type = float,default=3.e-4,required=False,
            help='Size of the rebined pixels in log lambda')

    parser.add_argument('--nspec',type = int,default=None,required=False,
            help="Number of spectra to fit, if None then run on all files")

    args = parser.parse_args()

desi_convert_transmission_to_delta_files(args.object_cat, args.out_dir, indir=args.in_dir, infiles=args.in_files, lObs_min=args.lambda_min, lObs_max=args.lambda_max, lRF_min=args.lambda_rest_min, lRF_max=args.lambda_rest_max, dll=args.dll, nspec=args.nspec)
