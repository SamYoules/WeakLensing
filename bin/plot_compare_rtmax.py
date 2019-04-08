## SY 8/3/19
## Plots power spectra for optimal map-making method with different rtmax.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

def load_kappa_maps(nside_opt, rtmax, sradius, k0):
    '''Reads in optimal estimated kappa maps for different search radii, and gets the Cls
       normalised by pixel area.'''

    data = np.load('kappa_opt_srad/kappa_cut_rt{}_rp100_nside{}_srad{}.npz'.format(rtmax, nside_opt, sradius))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    Cls          = hp.sphtfunc.anafast(kappa_opt, k0)/pixel_area
    ell          = np.arange(Cls.size)
    return ell, Cls

nside_opt = 128
sradius   = '05'

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

##- Open Midpoint kappa
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
    alm = hp.sphtfunc.map2alm(kappa_mid, lmax=3*nside_opt-1)
    kappa_mid = hp.sphtfunc.alm2map(alm, nside=nside_opt)

##- Get Midpoint Cls
clx_mid = hp.sphtfunc.anafast(k00, kappa_mid)
ell_mid = np.arange(clx_mid.size)

rtmax = [40, 70, 100]

for i,j in enumerate(rtmax):
    ell, Cl = load_kappa_maps(nside_opt, j, sradius, k00)
    P.plot(ell, Cl*ell, color=colors[i], lw=2, linestyle="-", \
                       label='rtmax={}'.format(rtmax[i]))

P.plot(ell_mid, clx_mid*ell_mid, color=colors[i+1], lw=2, linestyle="-", label='Midpoint')

#P.title('Estimated DESI kappa, nside={}, rtmax=100'.format(nside_opt))
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
P.annotate(r'DESI noiseless, nside={}, sradius={}'.format(nside_opt,sradius), xy=(0.15, 0.95), xycoords='axes fraction')

