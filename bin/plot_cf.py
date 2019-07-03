import numpy as N
import pylab as P
from astropy.io import fits
import picca.wedgize

def plot_cf_file():

    a = fits.open('cf.noisy.fits.gz')

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
    P.figure()
    P.ion()
    P.pcolormesh( rper, rpar, (r*mcf).reshape(50, 50))
    P.xlabel(r'$r_\perp$', fontsize=20)
    P.ylabel(r'$r_\parallel$', fontsize=20)
    P.colorbar()
    P.savefig('corr_func.png')

    #-- computing covariance matrix (important!)
    coss=N.zeros([nbins, nbins])
    for i in range(npix):
        print (i, 'of',  npix)
        coss+= N.outer(we[i]*(cf[i]-mcf), we[i]*(cf[i]-mcf))
    coss /= N.outer(mwe, mwe)

    #-- plotting wedges like in Fig. 8
    for mus in [[0., 0.5], [0.5, 0.8], [0.8, 0.95], [0.95, 1.]]:
        w = picca.wedgize.wedge(rpmin=0.0, rpmax=200.0, nrp=50, \
                rtmin=0.0, rtmax=200.0, nrt=50, \
                rmin=0.0, rmax=200.0, nr=50, \
                mumin=mus[0], mumax=mus[1], ss=10)
        r, wed, wedcov = w.wedge(mcf, coss)

        #-- errors are the square root of the diagonal elements of the covariance matrix
        dwed = N.sqrt(N.diag(wedcov))

        P.figure()
        #-- we multiply the wedges by r**2 so we can see the BAO peak
        P.errorbar(r, wed*r**2, dwed*r**2, fmt='o')
        P.title(r'$%.1f < \mu < %.1f$'%(mus[0], mus[1]))
        P.xlabel('r [Mpc/h]')
        #P.savefig('wedges_{}-{}.png'.format(mus[0], mus[1]))

plot_cf_file()
