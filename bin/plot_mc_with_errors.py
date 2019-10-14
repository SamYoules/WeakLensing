#!/usr/bin/env python

#-- Plot mean of 100 realisations of residuals (kappa_out - kappa_in) with RMS errors

import numpy as np
import pylab as plt
from astropy.table import Table
import glob

rt = 30
rp = 10
mcreal = 1000000
nreal = 100
var = 0.05
kin = np.array([-0.14, -0.12, -0.1, -0.08, -0.06, -0.04, -0.02, 0.,
             0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14])
residual = np.zeros(kin.size)
indir = 'monte_carlo_1'

allfiles = glob.glob(indir+'/*.fits') 
for i, filename in enumerate(allfiles):
    t = Table.read(filename)
    residual+=(t['kappa_out']-kin)

#-- Get RMSE of residual
rmse = np.sqrt(residual**2 / nreal)
mean_resid = residual / nreal

#-- plot the mean residuals with the error bars
plt.figure()
plt.title('rt = {}, rp = {}, nreal = {}, variance = {}'.format(rt, rp, mcreal, var))
plt.errorbar(kin, mean_resid, yerr=rmse, fmt='o', label='Mean residual')
plt.plot((-0.14,0.14),(0,0), 'g--', alpha=0.5)
plt.legend()
plt.xlabel('Kappa input')
plt.ylabel('Mean residuals')
plt.show()
plt.savefig('plots/mc/mc.png')