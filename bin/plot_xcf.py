import numpy as N
import pylab as plt
from astropy.io import fits
import picca.wedgize
import sys

def plot_xcf_file(xcf_file):

    a = fits.open(xcf_file)

    #-- these are 100x50=2500 bins vectors containing the separations of each bin
    rp = a[1].data.RP  #-- this is radial separation r_parallel
    rt = a[1].data.RT   #-- this is transverse separation r_perp
    z = a[1].data.Z      #-- this is the mean redshift (almost constant)
    we = a[2].data.WE  #-- this is the weight associated to each bin, and each pixel
    cf = a[2].data.DA    #-- this is the correlation function in each bin and each pixel

    np = a[1].header['NP']
    nt = a[1].header['NT']
    rpmax = a[1].header['RPMAX']
    rtmax = a[1].header['RTMAX']

    rper = (N.arange(nt)+0.5)*rtmax/nt
    rpar = (N.arange(np)+0.5)*rpmax/np

    #-- this shows the number of pixels and number of bins
    print (cf.shape)
    npix = cf.shape[0]
    nbins = cf.shape[1]

    #-- doing a weighted average to get the final correlation function
    mcf = N.sum(cf*we, axis=0)/N.sum(we, axis=0)
    mwe = N.sum(we, axis=0)

    r = N.sqrt(rp**2+rt**2)

    #-- making a 2D plot of the correlation function
    plt.figure()
    plt.ion()
    plt.pcolormesh( rper, rpar, (r*mcf).reshape(100, 50))
    plt.xlabel(r'$r_\perp$', fontsize=20)
    plt.ylabel(r'$r_\parallel$', fontsize=20)
    plt.colorbar()
    plt.savefig('corr_func.png')

    #-- computing covariance matrix (important!)
    coss=N.zeros([nbins, nbins])
    for i in range(npix):
        if i%50==0:
            print (i, 'of',  npix)
        coss+= N.outer(we[i]*(cf[i]-mcf), we[i]*(cf[i]-mcf))
    coss /= N.outer(mwe, mwe)

    #-- plotting wedges like in Fig. 8

    #for mus in [[0., 0.5], [0.5, 0.8], [0.8, 0.95], [0.95, 1.]]:
    for mus in [[-1., -0.5], [-0.5, 0.0], [0.0, 0.05], [0.05, 1.]]:
        w = picca.wedgize.wedge(rpmin=-200.0, rpmax=200.0, nrp=100, \
                rtmin=0.0, rtmax=200.0, nrt=50, \
                rmin=0.0, rmax=200.0, nr=50, \
                mumin=mus[0], mumax=mus[1], ss=10)
        r, wed, wedcov = w.wedge(mcf, coss)

        #-- errors are the square root of the diagonal elements of the covariance matrix
        dwed = N.sqrt(N.diag(wedcov))

        #-- we multiply the wedges by r**2 so we can see the BAO peak
        plt.errorbar(r, wed*r**2, dwed*r**2, fmt='o')
        plt.title(r'$%.1f < \mu < %.1f$'%(mus[0], mus[1]))
        plt.xlabel('r [Mpc/h]')
        #plt.savefig('wedges_{}-{}.png'.format(mus[0], mus[1]))

xcf_file = sys.argv[1]
plot_xcf_file(xcf_file)
