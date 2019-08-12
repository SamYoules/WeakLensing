#!/usr/bin/env python
  
import numpy as np
import pylab as P
import scipy as sp
import fitsio
from scipy import random
from scipy import interpolate
import h5py

file_xi='xcf-exp-noisy.out.gz'
file_fit='xcf_noisy.h5'
nbins=50

h = fitsio.FITS(file_xi)
ff = h5py.File(file_fit, 'r')
base = file_xi 
fit = ff[base+'/fit'][...] 
data_rp = h[1]['RP'][:]
data_rt = h[1]['RT'][:]
hh = h[1].read_header()
rpmin = hh['RPMIN']
rpmax = hh['RPMAX']
rtmin = 0 
rtmax = hh['RTMAX']
h.close()
ff.close()

rpmin = data_rp.reshape(100, 50)[0].max()
rpmax = data_rp.reshape(100, 50)[-1].min()
rtmin = data_rt.reshape(100, 50)[:, 0].max()
rtmax = data_rt.reshape(100, 50)[:, -1].min()

#-- create the regular grid for griddata
rp = np.linspace(rpmin, rpmax, nbins*2)
rt = np.linspace(rtmin, rtmax, nbins)
xim = sp.interpolate.griddata((data_rt, data_rp), fit,
           (np.outer(np.ones(rp.size), rt).ravel(),
           np.outer(rp, np.ones(rt.size)).ravel()), method='cubic')

#-- create interpolator object
xi2d = sp.interpolate.RectBivariateSpline(rt, rp, \
           xim.reshape((50, 100)) )

P.figure()
P.pcolormesh([data_rp, data_rt])
P.savefig('xi_fit.png')

xim1 = xim.reshape(50,100)
P.figure()
P.pcolormesh([xim1[0], xim1[1]])
P.savefig('xi_xim.png')

