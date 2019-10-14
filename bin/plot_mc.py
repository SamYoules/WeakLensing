#!/usr/bin/env python

import numpy as np
import pylab as plt
import sys
from astropy.table import Table

rt = sys.argv[1]
t = Table.read('monte_carlo/monte_carlo_results{}.fits'.format(rt))
n = t['nrealisations'][0]
var = t['var1'][0]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plt.figure()
plt.title('rt = {}, nreal = {}, variance = {}'.format(rt, n, var))
#rp = np.unique(t['rp'])
rp = [10,20,30,40,50,60,70]
for i, j in enumerate(rp):
    w = j == t['rp']
    kin = t['kappa_in'][w]
    kout = t['kappa_out'][w]
    plt.plot(kin, kout, color=colors[i], label='rp = {}'.format(j))
plt.legend()
plt.xlabel('kappa input')
plt.ylabel('kappa output')
plt.show()
plt.savefig('plots/mc/mc_rt{}.png'.format(rt))
