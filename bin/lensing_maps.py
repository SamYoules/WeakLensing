import healpy as hp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl


def make_maps(a3darraymap): # a function that makes maps from gamma 1 and gamma 2 shear measurements - needs a 3d array of an empty map, a g1 map and a g2map

        NSIDE=hp.npix2nside(len(a3darraymap[1])) #finds nside of a map from just input map

        LMAX=2*NSIDE-1

        ells=hp.sphtfunc.Alm.getsize(lmax=LMAX)

        lmode=np.zeros(ells)
        for jh in range(0, ells):
                lmode[jh]=hp.sphtfunc.Alm.getlm(lmax=LMAX, i=jh)[0] #puts ell modes in an array, doesnt need to be a loop
        lfacsphi=-2/np.sqrt((lmode+2)*(lmode+1)*(lmode)*(lmode-1)) #ell mode coefficients to go from E_alm to gravitational potential
        lfacsphi[np.where(lmode<2)]=0
        lfacsphi[np.isnan(lfacsphi)]=0
        lfacsphi[np.isinf(lfacsphi)]=0

        lfacskappa=(lmode*(lmode+1)*((lmode+2)**-1)*((lmode-1)**-1))**0.5 ##ell mode coefficients to go from E_alm to kappa

        lfacskappa[np.where(lmode<2)]=0
        lfacskappa[np.isinf(lfacskappa)]=0
        lfacskappa[np.isnan(lfacskappa)]=0

        ## SY 30/5/18 Added minus sign below
        alphafacs=-2/np.sqrt((lmode+2)*(lmode-1)) #ell mode coefficients to go from E_alm to alpha, bend angle
        alms=hp.sphtfunc.map2alm(a3darraymap, lmax=LMAX, pol=True)
        Ealm=alms[1]
        Kalm=lfacskappa*Ealm
        Palm=lfacsphi*Ealm

        alphaplus=alphafacs*alms[2]
        alphaminus=alphafacs*alms[1]

        alphaplus[np.where(lmode<2)]=0
        alphaminus[np.where(lmode<2)]=0 #2 different bend angles. 
        Kmap=hp.sphtfunc.alm2map(Kalm, nside=NSIDE, lmax=LMAX)
        Pmap=hp.sphtfunc.alm2map(Palm, nside=NSIDE, lmax=LMAX)
        Amap=hp.alm2map_spin([alphaplus, alphaminus], nside=NSIDE, spin=1, lmax=LMAX) # spin1 transform to get the 2 bend angles
        return Kmap, Pmap, Amap



NSIDE=256 #defines pixel scale, can use hp.nside2resol to get a pixel scale
LMAX=2*NSIDE-1 #defines the biggest l mode used in some spherical approaches (used in some functions)

Cls=np.loadtxt('c_ells.txt') #loading in a CAMB C_ls power spectrum file
intercl=np.zeros(int(np.amax(Cls[:,0])))

cli=np.interp(x=np.arange(1, np.amax(Cls[:,0])), xp=Cls[:,0],fp=Cls[:,1]) #doing an interpolation so that the input Cl for synfast is in the right format
#intercl[0:]=cli
intercl[1:]=cli

obsK=hp.sphtfunc.synfast(intercl, nside=NSIDE, lmax=LMAX, pol=True) #Using synfast to generate a Gaussian random field using the power spectrum, chose to define this as observed kappa as opposed to observed E mode. Pol=False means a scalar field is output only. Can play with these to get more types of fields out
lmode , em =hp.sphtfunc.Alm.getlm(lmax=LMAX) #gives coefficients in ell and em


obsKalm=hp.sphtfunc.map2alm(obsK, lmax=LMAX, pol=False) # Kappa map into harmonic space

LFACSEB=(lmode*(lmode+1)/((lmode+2)*(lmode-1)))**(-0.5) #inverse of coefficients that relate the kappa field and the "E mode type" measurement from lensing shear. 

LFACSEB[np.isinf(LFACSEB)]=0
LFACSEB[np.isnan(LFACSEB)]=0

Etruealm=LFACSEB*obsKalm #Getting a E mode in spherical space from this kappa 

[kk, g1in, g2in]=hp.sphtfunc.alm2map([Etruealm*0.0, Etruealm, Etruealm*0.0], nside=NSIDE, lmax=LMAX, pol=True) #Getting a shear (spin-2 ) field from the E modes

# SY 13/4/18
a3darraymap = np.array([g1in*0.0, g1in, g2in]) # first array needs to be empty & same dimension as others
Kmap, Pmap, Amap = make_maps(a3darraymap)
hp.fitsfunc.write_map('outputs/kmap.fits', Kmap)
hp.fitsfunc.write_map('outputs/amap.fits', Amap)
hp.fitsfunc.write_map('outputs/pmap.fits', Pmap)
