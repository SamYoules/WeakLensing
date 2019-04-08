## SY 8/3/19
## Plots power spectra for optimal map-making method with different search radii.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

def load_kappa_maps(nside, k00):
    '''Reads in optimal estimated kappa maps for different nsides, and gets the Cls
       normalised by pixel area.'''

    data = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside{}_srad01.npz'.format(nside))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    Cls          = hp.sphtfunc.anafast(kappa_opt, k00)/pixel_area**0.5
    ell          = np.arange(Cls.size)
    return ell, Cls

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Cls
k0 = fits.open('est_maps_cut/kappa_input1.fits')[1].data.I
#k0 = fits.open('kappalists/kappa_input.fits')[1].data.kappa.ravel()
Cl0 = hp.sphtfunc.anafast(k0, k0)
ell0 = np.arange(Cl0.size)
P.plot(ell0, Cl0*ell0/400, color='k', lw=2, linestyle="-", label='Input/400')

#nsides = [16, 32, 64, 128, 256]
nsides = [256]
for i,j in enumerate(nsides):
    k00 = hp.ud_grade(k0, j)
    ell, Cl = load_kappa_maps(j,k00)
    ell, Cl = load_kappa_maps(j,k00)
    #P.plot(ell, Cl*ell, color=colors[i], lw=2, linestyle="-", label='nside={}'.format(j))
    P.plot(ell, Cl*ell, color=colors[i], lw=2, linestyle="-", label='Optimal')

##- Get Midpoint Cls
kappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside256.fits.gz')[1].data.kappa
wkappa_mid = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100-nside256.fits.gz')[1].data.wkappa
nside_mid=int(hp.npix2nside(kappa_mid.size))

##- Reset resolution of input map to match midpoint maps
k_in = hp.ud_grade(k0, nside_mid)

##- Mask off areas outside DESI footprint
mask = wkappa_mid!=0
mask &= (kappa_mid>np.percentile(kappa_mid[mask], 0.5)) & \
                (kappa_mid<np.percentile(kappa_mid[mask], 99.5))
kappa_mid[~mask]=hp.UNSEEN

##- Mask area outside footprint
k_in = k_in*(mask)+hp.UNSEEN*(~mask) 
k_in[~mask]=hp.UNSEEN

cl_mid=hp.sphtfunc.anafast(kappa_mid, k_in, lmax=3*nside_mid-1)
ell_mid = np.arange(cl_mid.size)
#P.plot(ell_mid, cl_mid*ell_mid, color=colors[6], lw=2, linestyle="-", label='Midpoint, 256')
P.plot(ell_mid, cl_mid*ell_mid/130, color=colors[2], lw=2, linestyle="-", label='Midpoint/130')

##- fit for large-angles
cl_in=hp.sphtfunc.anafast(k_in, lmax=3*nside_mid-1)
y = (cl_mid / cl_in)
x = np.arange(y.size)
z = np.polyfit(x, y, 5)
model = np.polyval(z, x)
#P.plot(x, cl_in*x*model, color=colors[3], lw=2, linestyle="-", label='Masked Damped Input')

#P.title('DESI: kappa estimators, nside=256, rtmax=rpmax=100, sradius=0.1')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower left', fancybox=True, framealpha=0.5)
#P.annotate(r'DESI noiseless, sradius=0.1, rtmax=rpmax=100', xy=(0.12, 0.95), xycoords='axes fraction')

#P.savefig('plots/compare_eBOSS_srad_128.pdf', bbox_inches='tight')
#P.savefig('plots/compare_desi_128.pdf', bbox_inches='tight')

