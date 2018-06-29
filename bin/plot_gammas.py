## SY

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys
import configargparse


def open_gaussian():

    g1 = fits.open('gamma-gaussian-true-v3.fits.gz')[1].data.gamma1
    g2 = fits.open('gamma-gaussian-true-v3.fits.gz')[1].data.gamma2
    kap = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.kappa
    wkap = fits.open('kappa-gaussian-true-v3.fits.gz')[1].data.wkappa
    ktrue = fits.open('kappa_input.fits')[1].data.I
    return g1, g2, kap, wkap, ktrue


def open_blob():

    g1 = fits.open('gamma-blob.fits.gz')[1].data.gamma1
    g2 = fits.open('gamma-blob.fits.gz')[1].data.gamma2
    kap = fits.open('kappa-blob.fits.gz')[1].data.kappa
    wkap = fits.open('kappa-blob.fits.gz')[1].data.wkappa
    ktrue = fits.open('kappa_blob_input.fits')[1].data.I
    return g1, g2, kap, wkap, ktrue


def plot_mollviews(ga1, ga2, k_g1, k_g2, kt, kg):

    hp.mollview(ga1, min=-0.5, max=0.5, title = 'gamma1 (blob) estimator')
    P.savefig('plots/gamma1_estim_blob.png')
    P.close()

    hp.mollview(k_g1, min=-0.5, max=0.5, title = 'gamma1 (blob) input')
    P.savefig('plots/gamma1_input_blob.png')
    P.close()

    hp.mollview(ga2, min=-0.5, max=0.5, title = 'gamma2 (blob) estimator')
    P.savefig('plots/gamma2_estim_blob.png')
    P.close()

    hp.mollview(k_g2, min=-0.5, max=0.5, title = 'gamma2 (blob) input')
    P.savefig('plots/gamma2_input_blob.png')
    P.close()

    hp.mollview(kt, min=-0.5, max=0.5, title = 'kappa (blob) input')
    P.savefig('plots/kappa_input_blob.png')
    P.close()

    hp.mollview(kg, min=-0.5, max=0.5, title = 'kappa (blob) estimator')
    P.savefig('plots/kappa_estim_blob.png')
    P.close()


def plot_gamma(ga1, ga2, k_g1, k_g2):

    # Get the theoretical values for ells and C_ells
    th = kappa_lya.Theory()
    ell, cell = th.get_cl_kappa(2.1)

    # Adjust kappa C_ells to get gamma C_ells
    # from equations 48 - 51 in Castro et al
    for i, c_ell in enumerate(cell):
        l = i + 1
        cell[i] = (c_ell*(l+2)*(l-1))/(l*(l+1))

    # Get the angular power spectrum (Cls) for auto- and cross-correlations
    Cls1=hp.sphtfunc.anafast(ga1, lmax=3*NSIDE-1)
    Cls2=hp.sphtfunc.anafast(ga2, lmax=3*NSIDE-1)
    Cls3=hp.sphtfunc.anafast(k_g1, lmax=3*NSIDE-1)
    Cls4=hp.sphtfunc.anafast(k_g2, lmax=3*NSIDE-1)
    Cls12=hp.sphtfunc.anafast(ga1, ga2, lmax=3*NSIDE-1)
    Cls13=hp.sphtfunc.anafast(ga1, k_g1, lmax=3*NSIDE-1)
    Cls24=hp.sphtfunc.anafast(ga2, k_g2, lmax=3*NSIDE-1)

    # gamma1 plot (divide cell by 4 to match eBOSS footprint
    #              and by 2 because gamma?)
    P.plot(ell, cell/8, label='Theory')
    P.plot(Cls1, label='$\gamma_{1 est}$')
    P.plot(Cls3, label='$\gamma_{1 true}$')
    P.plot(Cls13, label='$\gamma_{1 est}\gamma_{1 true}$')
    P.xlim(10,200)
    P.ylabel('$C_l$', fontsize=18)
    P.xlabel('$l$', fontsize=18)
    P.legend()
    P.title(r'$\gamma {\rm \ 1 \ Correlations}$')
    P.savefig('plots/gamma1_corr.png')
    P.close()

    # gamma2 plot
    P.plot(ell, cell/8, label='Theory')
    P.plot(Cls2, label='$\gamma_{2 est}$')
    P.plot(Cls4, label='$\gamma_{2 true}$')
    P.plot(Cls24, label='$\gamma_{2 est}\gamma_{2 true}$')
    P.xlim(10,200)
    P.ylabel('$C_l$', fontsize=18)
    P.xlabel('$l$', fontsize=18)
    P.legend()
    P.title(r'$\gamma {\rm \ 2 \ Correlations}$')
    P.savefig('plots/gamma2_corr.png')
    P.close()


def plot_kappa(ga1, ga2, k_k, k_t):
    '''Derive kappa_g from the estimated gamma1 and gamma2, and then plot
       auto- and cross-correlations'''

    # Make a 3D array using the shear components
    a3darraymap = np.array([ga1*0.0, ga1, ga2]) 

    # Find nside from input map
    nside = hp.npix2nside(len(a3darraymap[1])) 
    lmax=3*nside-1
    nell = hp.sphtfunc.Alm.getsize(lmax=lmax)
    
    # Put ell modes in an array
    lmode = hp.sphtfunc.Alm.getlm(lmax=lmax, i=np.arange(nell))[0] 

    # ell mode coefficients to go from E_alm to kappa
    lfactor_kappa = np.sqrt((lmode*(lmode+1)) / ((lmode+2)*(lmode-1)))
    lfactor_kappa[lmode<2]=0

    # Get E-mode alms from input map
    alms = hp.sphtfunc.map2alm(a3darraymap, lmax=lmax, pol=True)

    # Compute kappa alms
    Kalm = lfactor_kappa*alms[1]

    # Get kappa_g map from kappa alms
    k_g = hp.sphtfunc.alm2map(Kalm, nside=nside, lmax=lmax)
    # Fix amplitude problem coming from the gammas
    k_g=k_g*2

    # Get the C_ells for auto- and cross-correlations
    Cls1=hp.sphtfunc.anafast(k_g, lmax=3*NSIDE-1)
    Cls2=hp.sphtfunc.anafast(k_k, lmax=3*NSIDE-1)
    Cls3=hp.sphtfunc.anafast(k_t, lmax=3*NSIDE-1)
    Cls12=hp.sphtfunc.anafast(k_g, k_k, lmax=3*NSIDE-1)
    Cls13=hp.sphtfunc.anafast(k_g, k_t, lmax=3*NSIDE-1)
    Cls23=hp.sphtfunc.anafast(k_k, k_t, lmax=3*NSIDE-1)

    # Get the theoretical values for ells and C_ells
    th = kappa_lya.Theory()
    ell, cell = th.get_cl_kappa(2.1)

    # kappa plots
    P.figure(figsize=(10,8))
    # The factor of 4 below accounts for eBOSS footprint ~ 1/4 full sky
    #P.plot(ell, cell/4, label='Theory')

    #P.plot(Cls1, label='$\kappa_{\gamma \gamma} $')
    #P.plot(Cls2, label='$\kappa_{\kappa \kappa} $')
    #P.plot(Cls3, label='$\kappa_{tt} $')

    #P.plot(Cls12, label='$\kappa_{\kappa \gamma} $')
    #P.plot(Cls13, label='$\kappa_{\gamma t} $')
    #P.plot(Cls23, label='$\kappa_{\kappa t} $')

    #P.plot(Cls23/Cls3, label='$\kappa_{\kappa t} / \kappa_{tt} $')
    #P.plot(Cls13/Cls3, label='$\kappa_{\gamma t} / \kappa_{tt} $')

    P.plot(Cls2 - Cls3, label='$\kappa_{\kappa \kappa} - \kappa_{tt} $')
    P.plot(Cls1 - Cls3, label='$\kappa_{\gamma \gamma} - \kappa_{tt} $')
    P.plot(Cls12 - Cls3, label='$\kappa_{\kappa \gamma} - \kappa_{tt} $')

    P.xlim(10,200)
    #P.ylim(-0.5e-8,2.1e-8)
    P.ylabel('$C_l$', fontsize=18)
    P.xlabel('$l$', fontsize=18)
    P.legend(fontsize=22)
    P.title(r'$\kappa {\rm \ Correlations}$')
    #P.savefig('plots/kappa_corr_b.png')
    #P.close()



# input: enter 'g' or 'b' for gaussian / blob
filetype = sys.argv[1]
# input: enter 'g' or 'k' for gamma / kappa or 'n' for none
outtype = sys.argv[2]

# Open files
if (filetype == 'b'):
    g1, g2, kap, wkap, ktrue = open_blob()
else:
    g1, g2, kap, wkap, ktrue = open_gaussian()
    #g1 = g1*2
    #g2 = g2*2

# Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(kap.size))
alm = hp.sphtfunc.map2alm(ktrue, lmax=3*NSIDE-1)
ktrue = hp.sphtfunc.alm2map(alm, nside=NSIDE)

# Mask off areas outside eBOSS footprint
w = (wkap == 0)
kap[w] = hp.UNSEEN
ktrue[w] = hp.UNSEEN
g1[w] = hp.UNSEEN
g2[w] = hp.UNSEEN

# Get gamma1 and gamma2 from input kappa file
kappa_alm = hp.map2alm(ktrue, lmax=3*NSIDE-1)

# Get coefficients in ell and em
LMAX=3*NSIDE-1
lmode , em =hp.sphtfunc.Alm.getlm(lmax=LMAX)

# Inverse of coefficients that relate the kappa field and the "E mode type"
# measurement from lensing shear.
LFACSEB=(lmode*(lmode+1)/((lmode+2)*(lmode-1)))**(-0.5)  

LFACSEB[np.isinf(LFACSEB)]=0
LFACSEB[np.isnan(LFACSEB)]=0

# Get an E mode in spherical space from this kappa
Etruealm=LFACSEB*kappa_alm 

# Get gamma1 and gamma2 from E modes
[kdummy, g1t, g2t] = hp.alm2map([Etruealm*0.0, Etruealm, Etruealm*0.0], \
                       nside=NSIDE, pol=True)

# Mask off areas outside eBOSS footprint
g1t[w] = hp.UNSEEN
g2t[w] = hp.UNSEEN

# Make plots
if (outtype == 'g'):
    plot_gamma(g1, g2, g1t, g2t)
elif (outtype == 'k'):
    plot_kappa(g1, g2, kap, ktrue)
else:
    plot_mollviews(g1, g2, g1t, g2t, ktrue, kap)




