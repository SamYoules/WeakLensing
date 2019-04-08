## SY 8/3/19
## Plots power spectra for optimal map-making method with different search radii.
## Maps made with get_map.py

import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits

def load_kappa_maps(cutoff, k0):
    '''Reads in optimal estimated kappa maps for different cutoffs, and gets the Cls
       normalised by pixel area.'''

    data = np.load('kappa_opt_cutoff/kappa_eBOSS_rt70_rp40_nside128_cutoff{}.npz'.format(cutoff))
    nside_opt    = int(data['arr_0'])
    pixel_ids    = data['arr_1']
    pixel_kappas = data['arr_2']
    pixel_area   = hp.nside2pixarea(nside_opt)
    kappa_opt    = np.zeros(nside_opt**2*12)
    kappa_opt[pixel_ids] = pixel_kappas
    Cls          = hp.sphtfunc.anafast(kappa_opt, k0)/pixel_area
    ell          = np.arange(Cls.size)
    return ell, Cls

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Cls
k0 = fits.open('kappalists/kappa_input.fits')[1].data.kappa.ravel()
Cl0 = hp.sphtfunc.anafast(k0, k0)
ell0 = np.arange(Cl0.size)
P.plot(ell0, Cl0*ell0, color='k', lw=2, linestyle="-", label='Input')
k00 = hp.ud_grade(k0, 128)

cutoff = [2, 4, 10, 100]
for i,j in enumerate(cutoff):
    ell, Cl = load_kappa_maps(j,k00)
    P.plot(ell, Cl*ell, color=colors[i], lw=2, linestyle="-", \
                       label='cutoff={}'.format(j))

P.title('Estimated eBOSS kappa')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
P.savefig('plots/compare_cutoffs.pdf', bbox_inches='tight')

