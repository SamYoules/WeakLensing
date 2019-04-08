import healpy as hp
import numpy as np
import pylab as P
from astropy.io import fits


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

k00 = hp.ud_grade(k0, 256)

data = np.load('kappa_opt_srad/kappa_cut_rt100_rp100_nside256_srad01.npz')
nside_opt    = int(data['arr_0'])
pixel_ids    = data['arr_1']
pixel_kappas = data['arr_2']
pixel_area   = hp.nside2pixarea(nside_opt)
kappa_opt    = np.zeros(nside_opt**2*12)
kappa_opt[pixel_ids] = pixel_kappas
Cls          = hp.sphtfunc.anafast(kappa_opt, k00)
ell          = np.arange(Cls.size)
cl1 = Cls/pixel_area
cl2 = Cls/pixel_area**0.5

P.plot(ell0, Cl0*ell0/3, color='k', lw=2, linestyle="-", label='Input/3')
P.plot(ell, cl1*ell, color=colors[2], lw=2, linestyle="-", label='pixarea')
P.plot(ell0, Cl0*ell0/731, color=colors[4], lw=2, linestyle="-", label='Input/731')
P.plot(ell, cl2*ell, color=colors[3], lw=2, linestyle="-", label='sqrt(pixarea)')


#P.title('Estimated DESI kappa, nside=256, rtmax=100, sradius=0.1')
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower left', fancybox=True, framealpha=0.5)
P.annotate(r'DESI noiseless, nside=256, rtmax=rpmax=100, sradius=0.1', xy=(0.05, 0.95), xycoords='axes fraction')

