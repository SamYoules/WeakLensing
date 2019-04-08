## SY 8/3/19
## Plots power spectra for optimal map-making method with and without Gaussian damping.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

nside = 128
srad  = 0.1

##- Get Cls
#k0           = fits.open('kappalists/kappa_input.fits')[1].data.kappa.ravel()
k0           = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
Cls0         = hp.sphtfunc.anafast(k0, k0)
ell0         = np.arange(Cls0.size)
k00          = hp.ud_grade(k0, nside)

#data         = np.load('kappa_opt_eBOSS/kappa_eBOSS_rt70_rp40_nside_128.npz')
data         = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside128_srad01_nogaussian.npz')
nside_k1     = int(data['arr_0'])
ids          = data['arr_1']
pixel_kappas = data['arr_2']
pixel_area   = hp.nside2pixarea(nside_k1)
k1           = np.zeros(nside_k1**2*12)
k1[ids]      = pixel_kappas
Cls1         = hp.sphtfunc.anafast(k1, k00)/pixel_area
ell1         = np.arange(Cls1.size)

data         = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside128_srad01.npz')
nside_k2     = int(data['arr_0'])
ids          = data['arr_1']
pixel_kappas = data['arr_2']
pixel_area   = hp.nside2pixarea(nside_k2)
k2           = np.zeros(nside_k2**2*12)
k2[ids]      = pixel_kappas
Cls2         = hp.sphtfunc.anafast(k2, k00)/pixel_area
ell2         = np.arange(Cls2.size)

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
if (nside_mid != nside):
    kappa_mid = hp.ud_grade(kappa_mid, nside)

##- Get Midpoint Cls
clx_mid = hp.sphtfunc.anafast(k00, kappa_mid)
ell_mid = np.arange(clx_mid.size)

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

P.plot(ell0, Cls0*ell0, color='k', lw=2, linestyle="-", label='Input')
P.plot(ell1, Cls1*ell1, color=colors[0], lw=2, linestyle="-", label='No damping')
P.plot(ell2, Cls2*ell2, color=colors[1], lw=2, linestyle="-", label='Gaussian damping')
P.plot(ell_mid, clx_mid*ell_mid, color=colors[2], lw=2, linestyle="-", label='Midpoint')

P.title('Estimated x Input DESI kappa, nside={}, srad={}'.format(nside, srad))
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
P.savefig('plots/compare_damping.pdf', bbox_inches='tight')

