#!/usr/bin/env python

##-- Compute and plot weighted mean of kappa rings

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


def compute_mean(j, band, mapnum):
    '''compute weighted mean and standard deviation for ring kappas'''

    skxa=0.
    skaa=0.

    #kx=fits.open('maps/midpoint/blob/rp_bands/rings_broad/kappa{}_ring{}.fits.gz'
    #                           .format(band,mapnum))[1].data['kappa']
    #kxw=fits.open('maps/midpoint/blob/rp_bands/rings_broad/kappa{}_ring{}.fits.gz'
    #                          .format(band,mapnum))[1].data['wkappa']
    kx=fits.open('maps/midpoint/blob/rp_bands/rings/kappa{}_ring{}.fits.gz'
                               .format(band,mapnum))[1].data['kappa']
    kxw=fits.open('maps/midpoint/blob/rp_bands/rings/kappa{}_ring{}.fits.gz'
                              .format(band,mapnum))[1].data['wkappa']
    mask  = kx!=0
    kx    = kx[mask]
    kxw   = kxw[mask]

    #- weighted average, variance
    m   = sp.sum(kxw*kx) / sp.sum(kxw)
    var = sp.sum(kxw*kx*kx) / sp.sum(kxw) - m**2 

    #- error
    xsigma = np.sqrt(var)

    #- error in the mean
    msigma = np.sqrt(var)/np.sqrt(kxw.size)

    if j == 0:
        ki = fits.open('maps/midpoint/blob/input/kappa_input_ring{}.fits.gz'
                                   .format(mapnum))[1].data['T']
        ki = ki[mask]
        kimean = np.mean(ki)
    else:
        kimean = 0
    return(m, msigma, kimean)


band = ['_20_16', '_16_12', '_12_8', '_8_4', '_4_0', '_0_4', '_4_8', '_8_12', '_12_16', '_16_20']
narrowband = ['-20-16', '-16-12', '-12-8', '-8-4', '-4-0', '0-4', '4-8', '8-12', '12-16', '16-20']
#band = [1,2,3,4,5,6,7,8,9,10]
#fatband = ['-100 to -80','-80 to -60','-60 to -40','-40 to -20','-20 to 0','0 to 20','20 to 40','40 to 60','60 to 80','80 to 100']

input_mean    = []

x=[1,3,5,7,9,11,13,15]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
 '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plt.figure()

for j, b in enumerate(band):
    weighted_xmean = []
    yerr           = []
    for i in np.arange(len(x)):
        wxm, msig, kimean = compute_mean(j, b, i+1)
        weighted_xmean.append(wxm)
        yerr.append(msig)
        if j==0:
            input_mean.append(kimean)

    plt.errorbar(x, weighted_xmean, yerr=yerr, label=narrowband[j], marker = '.', fmt='-o',
           color=colors[j], capsize=5, elinewidth=2, markeredgewidth=2,
           markeredgecolor=colors[j])

plt.plot(x, input_mean, 'o', color='k', label='Input')
plt.plot((0,16),(0,0), 'g--', alpha=0.5)

plt.xticks(x)
plt.xlabel('Ring Radius (degrees)')
plt.ylabel('Mean kappa')
plt.title('Annuli around blob for different bands of rp')
plt.legend()
plt.show()
#plt.savefig('plots/rp_bands/rings_narrow_rp.png')
