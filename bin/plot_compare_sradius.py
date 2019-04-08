## SY 8/3/19
## Plots power spectra for optimal map-making method with different search radii.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

def rebin(Cls,ellmin=1, ellmax=200, nell=8):
    '''Smooth curves by rebinning Cls'''
    ell = np.arange(Cls.size)
    weights = 2.*ell+1.
    if not ellmin:
        ellmin = min(ell)
    if not ellmax:
        ellmax = max(ell)
    if nell==0:
        nell = ellmax-ellmin

    w = (ell>=ellmin)&(ell<=ellmax)
    index = np.floor( (ell[w]-ellmin)*1./(ellmax-ellmin)*nell ).astype(int)
    well = np.bincount( index, weights=weights[w])
    sell = np.bincount( index, weights=weights[w]*ell[w])
    scl  = np.bincount( index, weights=weights[w]*Cls[w])
    ell = sell/well
    Cls = scl/well
        
    return ell, Cls 

def load_kappa_maps(nside_opt, sradius, k0):
    '''Reads in optimal estimated kappa maps for different search radii, and gets the Cls
       normalised by pixel area.'''

    data = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside{}_srad{}_nogaussian.npz'.format(nside_opt, sradius))
    #data = np.load('kappa_opt_dd/kappa_cut_rt100_rp100_nside{}_srad{}.npz'.format(nside_opt, sradius))
    #data = np.load('kappa_opt_srad/kappa_eBOSS_rt70_rp40_nside{}_srad{}.npz'.format(nside_opt, sradius))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    cls          = hp.sphtfunc.anafast(kappa_opt, k0)/pixel_area
    ell          = np.arange(cls.size)
    #ell, cls    = rebin(cls)
    return ell, cls


nside_opt = 128

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Input Cls
kappa_input = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
#kappa_input = fits.open('kappalists/kappa_input.fits')[1].data.kappa.ravel()
Cl0 = hp.sphtfunc.anafast(kappa_input, kappa_input)
ell0 = np.arange(Cl0.size)
P.plot(ell0, Cl0*ell0, color='k', lw=2, linestyle="-", label='Input')
k00 = hp.ud_grade(kappa_input, nside_opt)

##- Get Optimal Cls
suffix = ['10', '05', '01']
sradii = [1.0, 0.5, 0.1]
#suffix = ['1.5', '1', '08', '06', '05', '04', '02', '01', '005']
#sradii = [1.5, 1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.1, 0.05]
for i,j in enumerate(suffix):
    ell, Cl = load_kappa_maps(nside_opt, j, k00)
    P.plot(ell, Cl*ell, color=colors[i], lw=2, linestyle="-", \
                       label='sradius={}'.format(sradii[i]))


##- Open Midpoint kappa
#wkappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside64.fits.gz')[1].data.wkappa
#kappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside64.fits.gz')[1].data.kappa
wkappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside128.fits.gz')[1].data.wkappa
kappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside128.fits.gz')[1].data.kappa
nside_mid=int(hp.npix2nside(kappa_mid.size))

##- Mask off areas outside DESI footprint
mask = wkappa_mid!=0
mask &= (kappa_mid>np.percentile(kappa_mid[mask], 0.5)) & \
                (kappa_mid<np.percentile(kappa_mid[mask], 99.5))
kappa_mid[~mask]=hp.UNSEEN

##- Reset midpoint map resolution if necessary
if (nside_mid != nside_opt):
    kappa_mid = hp.ud_grade(kappa_mid, nside_opt)

##- Get Midpoint Cls
clx_mid = hp.sphtfunc.anafast(k00, kappa_mid)
ell_mid = np.arange(clx_mid.size)
#ell_mid, clx_mid = rebin(clx_mid)

P.plot(ell_mid, clx_mid*ell_mid, color=colors[i+1], lw=2, linestyle="-", label='Midpoint')

#P.title('Estimated DESI kappa, nside={}, rtmax=100'.format(nside_opt))
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
#P.annotate(r'DESI noiseless $\delta \delta$ Cross, nside={}, rtmax=rpmax=100', xy=(0.15, 0.95), xycoords='axes fraction'.format(nside_opt))
P.annotate(r'DESI noiseless, nside={}, rtmax=rpmax=100, no Gaussian'.format(nside_opt), xy=(0.05, 0.95), xycoords='axes fraction')
#P.annotate(r'eBOSS noiseless, nside={}, rtmax=70, rpmax=40'.format(nside_opt), xy=(0.15, 0.95), xycoords='axes fraction')

#P.savefig('plots/compare_eBOSS_srad_{}.pdf'.format(nside_opt), bbox_inches='tight')
#P.savefig('plots/compare_desi_{}.pdf'.format(nside_opt), bbox_inches='tight')

