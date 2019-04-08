## SY

from astropy.io import fits
import healpy as hp
import numpy as np
import pylab as P

def plot_mollviews(ga1, ga2, g1true, g2true, ktrue, kap):

    hp.mollview(ga1, min=-0.5, max=0.5, title = 'gamma1 estimator')
    P.savefig('plots/gamma1_estim_blob.png')
    P.close()

    hp.mollview(g1true, min=-0.5, max=0.5, title = 'gamma1 input')
    P.savefig('plots/gamma1_input_blob.png')
    P.close()

    hp.mollview(ga2, min=-0.5, max=0.5, title = 'gamma2 estimator')
    P.savefig('plots/gamma2_estim_blob.png')
    P.close()

    hp.mollview(g2true, min=-0.5, max=0.5, title = 'gamma2 input')
    P.savefig('plots/gamma2_input_blob.png')
    P.close()

    hp.mollview(ktrue, min=-0.5, max=0.5, title = 'kappa input')
    P.savefig('plots/kappa_input_blob.png')
    P.close()

    hp.mollview(kap, min=-0.5, max=0.5, title = 'kappa estimator')
    P.savefig('plots/kappa_estim_blob.png')
    P.close()


g1 = fits.open('gamma-blob.fits.gz')[1].data.gamma1
g2 = fits.open('gamma-blob.fits.gz')[1].data.gamma2
kap = fits.open('kappa-blob.fits.gz')[1].data.kappa
wkap = fits.open('kappa-blob.fits.gz')[1].data.wkappa
ktrue = fits.open('kappa_blob_input.fits')[1].data.I

plot_mollviews(g1, g2, g1t, g2t, ktrue, kap)




