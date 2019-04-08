## SY 7/3/19
## Plots power spectra for optimal map-making method with different nsides.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

def load_kappa_maps(nside, srad, k0):
    '''Reads in optimal estimated kappa maps for different nsides, and gets the Cls
       normalised by pixel area.'''

    #data = np.load('kappa_opt_eBOSS/kappa_eBOSS_rt70_rp40_nside_{}.npz'.format(nside))
    data = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside{}_srad{}_nogaussian.npz'.format(nside, srad))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    Cl1          = hp.sphtfunc.anafast(kappa_opt, k0)/pixel_area ## SY 12/3/19
    ell          = np.arange(Cl1.size)
    return ell, Cl1

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Cls
#k0 = fits.open('kappalists/kappa_input.fits')[1].data.kappa.ravel()
k0 = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
cl0 = hp.sphtfunc.anafast(k0, k0)
ell_input = np.arange(cl0.size)
P.plot(ell_input, cl0*ell_input, color='k', lw=2, linestyle="-", label='Input')

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
k00 = hp.ud_grade(k0, 128)
if (nside_mid != 128):
    kappa_mid = hp.ud_grade(kappa_mid, 128)

##- Get Midpoint Cls
clx_mid = hp.sphtfunc.anafast(k00, kappa_mid)
ell_mid = np.arange(clx_mid.size)

#nsides = [256, 128, 64, 32, 16]
nsides   = [128,64]
srad     = '01'
sradname = '0.1'
for i,j in enumerate(nsides):
    k00 = hp.ud_grade(k0, j)
    ell, Cl1 = load_kappa_maps(j, srad, k00)
    P.plot(ell, Cl1*ell, color=colors[i], lw=2, linestyle="-", label='nside={}'.format(j))

P.plot(ell_mid, clx_mid*ell_mid, color=colors[i+1], lw=2, linestyle="-", label='Midpoint')

P.title('Estimated DESI kappa, rtmax=100, srad={}'.format(sradname))
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
#P.savefig('plots/compare_nsides_norm.pdf', bbox_inches='tight')

