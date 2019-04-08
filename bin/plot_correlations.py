## SY 31/10/18
## Reads in power spectra of 100 realisations of estimated kappa for LYAxLYA and QSOxLYA and input kappa.
## Plots cross-cf error bars (for noisy and noiseless datasets) for input and estimated kappa.

import numpy as np
import pylab as P

tnoisy = np.loadtxt('corr_noisy.txt', skiprows=1)
tnoiseless = np.loadtxt('corr_noiseless.txt', skiprows=1)
tcut = np.loadtxt('corr_cut.txt', skiprows=1)


cax1=P.matshow(tnoisy, cmap='viridis_r')
#P.title('Noisy')
groups = ['1a','2a','3a','4a','5a', '6a', '7a', '1x','2x','3x','4x','5x', '6x', '7x', ]
 
x_pos = np.arange(len(groups))
P.xticks(x_pos,groups)
 
y_pos = np.arange(len(groups))
P.yticks(y_pos,groups)
P.colorbar(cax1)
P.show()

cax2=P.matshow(tnoiseless, cmap='viridis_r')
#P.title('High Density')
groups = ['1a','2a','3a','4a','5a', '6a', '7a', '1x','2x','3x','4x','5x', '6x', '7x', ]
 
x_pos = np.arange(len(groups))
P.xticks(x_pos,groups)
 
y_pos = np.arange(len(groups))
P.yticks(y_pos,groups)
P.colorbar(cax2)
P.show()

cax3=P.matshow(tcut, cmap='viridis_r')
#P.title('Noiseless')
groups = ['1a','2a','3a','4a','5a', '6a', '7a', '1x','2x','3x','4x','5x', '6x', '7x', ]
 
x_pos = np.arange(len(groups))
P.xticks(x_pos,groups)
 
y_pos = np.arange(len(groups))
P.yticks(y_pos,groups)
P.colorbar(cax3)
P.show()

