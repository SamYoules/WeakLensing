## SY
## Plot gamma and/or kappa auto- and cross-correlations

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P
import kappa_lya
from kappa_lya import *
import sys


def open_gaussian():

    g1 = fits.open(gfile)[1].data.gamma1
    g2 = fits.open(gfile)[1].data.gamma2
    kap = fits.open(kfile)[1].data.kappa
    wkap = fits.open(kfile)[1].data.wkappa
    ktrue = fits.open(infile)[1].data.I
    return g1*(-1), g2, kap, wkap, ktrue  # SY 22/8/18


def plot_gamma(ga1, ga2, ktrue, ell, cell, NSIDE):
    '''Plot auto-correlations of estimated-gamma and gamma-made-from-input-
       kappa, cross-correlation of both, and a computed theoretical value.'''

    #-- Get gamma1 and gamma2 from input kappa file

    #-- Get kappa alms from map
    kappa_alm = hp.map2alm(ktrue, lmax=3*NSIDE-1)

    #-- Get coefficients in ell and em
    LMAX=3*NSIDE-1
    lmode , em =hp.sphtfunc.Alm.getlm(lmax=LMAX)

    #-- Inverse of coefficients that relate the kappa field and the
    #-- "E mode type" measurement from lensing shear.
    LFACSEB=(lmode*(lmode+1)/((lmode+2)*(lmode-1)))**(-0.5)  

    LFACSEB[np.isinf(LFACSEB)]=0
    LFACSEB[np.isnan(LFACSEB)]=0

    #-- Get an E mode in spherical space from this kappa
    Etruealm=LFACSEB*kappa_alm 

    #-- Get gamma1 and gamma2 from E modes
    [kdummy, g1t, g2t] = hp.alm2map([Etruealm*0.0, Etruealm, Etruealm*0.0], \
                       nside=NSIDE, pol=True)

    #-- Mask off areas outside eBOSS footprint
    g1t[w] = hp.UNSEEN
    g2t[w] = hp.UNSEEN

    #-- Adjust kappa C_ells to get gamma C_ells
    #-- from equations 48 - 51 in Castro et al (astro-ph/0503479v2)
    for i, c_ell in enumerate(cell):
        l = i + 1
        cell[i] = (c_ell*(l+2)*(l-1))/(l*(l+1))

    #-- Get the angular power spectrum (Cls) for auto- and cross-correlations
    Cls1=hp.sphtfunc.anafast(ga1, lmax=3*NSIDE-1)
    Cls2=hp.sphtfunc.anafast(ga2, lmax=3*NSIDE-1)
    Cls3=hp.sphtfunc.anafast(g1t, lmax=3*NSIDE-1)
    Cls4=hp.sphtfunc.anafast(g2t, lmax=3*NSIDE-1)
    Cls13=hp.sphtfunc.anafast(ga1, g1t, lmax=3*NSIDE-1)
    Cls24=hp.sphtfunc.anafast(ga2, g2t, lmax=3*NSIDE-1)

    #-- gamma1 plot
    P.ion()
    P.figure()
    P.plot(ell, cell, label='Theory')
    P.plot(Cls1, label='$\gamma_{1 est, est}$')
    P.plot(Cls3, label='$\gamma_{1 true, true}$')
    P.plot(Cls13, label='$\gamma_{1 est, true}$')
    P.xlim(1,200)
    P.ylabel('$C_l^{\gamma \gamma}$', fontsize=18)
    P.xlabel('$l$', fontsize=18)
    P.legend()
    P.title(r'$\gamma {\rm \ 1 \ DESI \ Double \ Density}$')
    #P.title(r'$\gamma {\rm \ 1 \ DESI \ Half \ Density}$')
    #P.title(r'$\gamma {\rm \ 1 \ DESI}$')
    #P.title(r'$\gamma {\rm \ 1 \ Whole \ Sky, \ Double \ Density}$')
    #P.title(r'$\gamma {\rm \ 1 \ Whole \ Sky, \ eBOSS \ Density}$')
    #P.title(r'$\gamma {\rm \ 1 \ eBOSS \ Footprint}$')
    #P.title(r'$\gamma {\rm \ 1 \ Whole \ Sky, \ 1.5 Density}$')
    #P.title(r'$\gamma {\rm \ 1 \ eBOSS \ Footprint \ Half \ Density}$')
    #P.savefig('plots/gamma1_corr_ws.png')
    #P.close()

    #-- gamma2 plot
    P.figure()
    P.plot(ell, cell, label='Theory')
    P.plot(Cls2, label='$\gamma_{2 est, est}$')
    P.plot(Cls4, label='$\gamma_{2 true, true}$')
    P.plot(Cls24, label='$\gamma_{2 est, true}$')
    P.xlim(1,200)
    P.ylabel('$C_l^{\gamma \gamma}$', fontsize=18)
    P.xlabel('$l$', fontsize=18)
    P.legend()
    P.title(r'$\gamma {\rm \ 2 \ DESI \ Double \ Density}$')
    #P.title(r'$\gamma {\rm \ 2 \ DESI \ Half \ Density}$')
    #P.title(r'$\gamma {\rm \ 2 \ DESI}$')
    #P.title(r'$\gamma {\rm \ 2 \ Whole \ Sky, \ Double \ Density}$')
    #P.title(r'$\gamma {\rm \ 2 \ Whole \ Sky, \ eBOSS \ Density}$')
    #P.title(r'$\gamma {\rm \ 2 \ Whole \ Sky, \ 1.5 Density}$')
    #P.title(r'$\gamma {\rm \ 2 \ eBOSS \ Footprint}$')
    #P.title(r'$\gamma {\rm \ 2 \ eBOSS \ Footprint \ Half \ Density}$')
    #P.savefig('plots/gamma2_corr_ws.png')
    #P.close()


def plot_kappa(ga1, ga2, k_k, k_t, ell, cell):
    '''Derive kappa_g from the estimated gamma1 and gamma2, and then plot
       auto- and cross-correlations'''

    #-- Make a 3D array using the shear components
    a3darraymap = np.array([ga1*0.0, ga1, ga2]) 

    #-- Find nside from input map
    nside = hp.npix2nside(len(a3darraymap[1])) 
    lmax=3*nside-1
    nell = hp.sphtfunc.Alm.getsize(lmax=lmax)
    
    #-- Put ell modes in an array
    lmode = hp.sphtfunc.Alm.getlm(lmax=lmax, i=np.arange(nell))[0] 

    #-- ell mode coefficients to go from E_alm to kappa
    lfactor_kappa = np.sqrt((lmode*(lmode+1)) / ((lmode+2)*(lmode-1)))
    lfactor_kappa[lmode<2]=0

    #-- Get E-mode alms from input map
    alms = hp.sphtfunc.map2alm(a3darraymap, lmax=lmax, pol=True)

    #-- Compute kappa alms
    Kalm = lfactor_kappa*alms[1]

    #-- Get kappa_g map from kappa alms
    k_g = hp.sphtfunc.alm2map(Kalm, nside=nside, lmax=lmax)

    #-- Get the C_ells for auto- and cross-correlations
    Cls1=hp.sphtfunc.anafast(k_g, lmax=3*NSIDE-1)
    Cls2=hp.sphtfunc.anafast(k_k, lmax=3*NSIDE-1)
    Cls3=hp.sphtfunc.anafast(k_t, lmax=3*NSIDE-1)
    Cls12=hp.sphtfunc.anafast(k_g, k_k, lmax=3*NSIDE-1)
    Cls13=hp.sphtfunc.anafast(k_g, k_t, lmax=3*NSIDE-1)
    Cls23=hp.sphtfunc.anafast(k_k, k_t, lmax=3*NSIDE-1)

    #-- kappa plots
    P.ion()
    P.figure(figsize=(10,8))
    P.plot(ell, cell, label='Theory')
    P.plot(Cls1, label='$\kappa_{\gamma \gamma} $')
    P.plot(Cls2, label='$\kappa_{\kappa \kappa} $')
    P.plot(Cls3, label='$\kappa_{tt} $')
    P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
    P.xlim(1,200)
    P.xlabel('$l$', fontsize=18)
    P.legend(fontsize=22, loc = 'upper right')
    #P.title(r'$\kappa {\rm DESI \ Quarter \ Density \ Auto-Correlations}$')
    P.title(r'$\kappa {\rm DESI \ Double \ Density \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm DESI \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ Double \ Density, \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ eBOSS \ Density, \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ 1.5 Density, \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Half \ Density \ Auto-Correlations}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Auto-Correlations}$')

    P.figure(figsize=(10,8))
    P.plot(ell, cell, label='Theory')
    P.plot(Cls12, label='$\kappa_{\kappa \gamma} $')
    P.plot(Cls13, label='$\kappa_{\gamma t} $')
    P.plot(Cls23, label='$\kappa_{\kappa t} $')
    P.ylabel('$C_l^{\kappa \kappa}$', fontsize=18)
    P.xlim(1,200)
    P.xlabel('$l$', fontsize=18)
    P.legend(fontsize=22, loc = 'upper right')
    P.title(r'$\kappa {\rm DESI \ Double \ Density \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm DESI \ Half \ Density \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ DESI \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ Double \ Density, \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ eBOSS \ Density, \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ 1.5 Density, \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Half \ Density \ Cross-Correlations}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Cross-Correlations}$')

    P.figure(figsize=(10,8))
    P.plot(Cls23/Cls3, label='$\kappa_{\kappa t} / \kappa_{tt} $')
    P.plot(Cls13/Cls3, label='$\kappa_{\gamma t} / \kappa_{tt} $')
    P.ylabel('$C_l^{\kappa \kappa}: \ Estimator \ / \ True$', \
            fontsize=18)
    P.xlim(1,200)
    P.xlabel('$l$', fontsize=18)
    P.legend(fontsize=22, loc = 'upper right')
    P.title(r'$\kappa {\rm \ DESI, \ Double \ Density \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ DESI, \ Half \ Density \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ DESI, \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ Double \ Density, \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ eBOSS \ Density, \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ 1.5 Density, \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Ratio \ of \ C_l^{\kappa \kappa}}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Half \ Density, \ Ratio \ of \ C_l^{\kappa \kappa}}$')

    P.figure(figsize=(10,8))
    P.plot(Cls2 - Cls3, label='$\kappa_{\kappa \kappa} - \kappa_{tt} $')
    P.plot(Cls1 - Cls3, label='$\kappa_{\gamma \gamma} - \kappa_{tt} $')
    P.plot(Cls12 - Cls3, label='$\kappa_{\kappa \gamma} - \kappa_{tt} $')
    P.ylabel('$C_l^{\kappa \kappa} \ Estimator \ - \ True$', fontsize=18)
    P.xlim(1,200)
    P.xlabel('$l$', fontsize=18)
    P.legend(fontsize=22, loc = 'upper right')
    P.title(r'$\kappa {\rm \ DESI, \ Double \ Density \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ DESI, \ Half \ Density \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ DESI, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ Double \ Density, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ eBOSS \ Density, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ Whole \ Sky, \ 1.5 Density, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')
    #P.title(r'$\kappa {\rm \ eBOSS \ Footprint, \ Half \ Density, \ C_l^{\kappa \kappa} \ Estimator \ - \ True}$')

    ###win = hp.pixwin(nside=256, pol = False)
    ###for i,j in enumerate(Cls23):
    ###    j = j/win[i]**4
    ###    Cls13[i] = Cls13[i]/win[i]**4

    ###P.plot(Cls23/Cls3, label='$\kappa_{\kappa t} / \kappa_{tt} $')
    ###P.plot(Cls13/Cls3, label='$\kappa_{\gamma t} / \kappa_{tt} $')
    ###P.plot(win, label='$window-function$')
    ###P.plot(win**4, label='$window-function^4$')
    ###P.ylabel('$Ratio \ of \ C_l^{\kappa \kappa} \ Estimator \ and \ True$', \
    ###        fontsize=18)

    #P.savefig('plots/kappa_test.png')
    #P.close()


#-- input: enter 'g' or 'k' for gamma / kappa or 'n' for none
outtype = sys.argv[1]
infile = sys.argv[2]
kfile = sys.argv[3]
gfile = sys.argv[4]

#-- open kappa & gamma fits files
g1, g2, kap, wkap, ktrue = open_gaussian()

#-- Reset nside of input kappa file to match lower resolution maps
NSIDE=int(hp.npix2nside(kap.size))
alm = hp.sphtfunc.map2alm(ktrue, lmax=3*NSIDE-1)
ktrue = hp.sphtfunc.alm2map(alm, nside=NSIDE)

#-- Mask off areas outside eBOSS footprint
w = (wkap == 0)
mask_size = 1 - wkap[w].size/wkap.size
kap[w] = hp.UNSEEN
ktrue[w] = hp.UNSEEN
g1[w] = hp.UNSEEN
g2[w] = hp.UNSEEN

#-- Get the theoretical values for ells and C_ells
th = kappa_lya.Theory()
ell, cell = th.get_cl_kappa(2.1)
cell *= mask_size

print(mask_size)

#-- Make plots
if (outtype == 'g'):
    cell /= 2.
    plot_gamma(g1, g2, ktrue, ell, cell, NSIDE)
elif (outtype == 'k'):
    plot_kappa(g1, g2, kap, ktrue, ell, cell)




