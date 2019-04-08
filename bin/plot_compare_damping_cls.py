## SY 25/3/19
## Plots power spectra for optimal map-making method with different search radii.
## Maps made with get_map.py, cls made with write_spectrum.py

import numpy as np
import pylab as P

nside   = 64
rtmax   = 100
sradius = '05'
srad    = 0.5
damping = ['_01', '_05', '']
damp    = [0.1, 0.5, 1.0]

##- Setup figure
P.rcParams.update({'font.size':14})
P.ion()
P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Input Cls
cl_input  = np.loadtxt('kappa_opt_srad/cls/cl_input_{}.txt'.format(nside))
ell_input = np.arange(cl_input.size)
P.plot(ell_input, cl_input*ell_input, color='k', lw=2, linestyle="-", label='Input')

##- Get Midpoint Cls
cl_mid  = np.loadtxt('kappa_opt_srad/cls/cl_cross_cut_midpoint_100_{}.txt'.format(nside))
ell_mid = np.arange(cl_mid.size)
P.plot(ell_mid, cl_mid*ell_mid, color=colors[0], lw=2, linestyle="-", label='Midpoint')

##- Get Optimal Cls
for i,j in enumerate(damping):
    cl  = np.loadtxt('kappa_opt_srad/cls/cl_cross_cut_{}_{}_{}{}.txt'.format(nside,rtmax,sradius,j))
    ell = np.arange(cl.size)
    P.plot(ell, cl*ell, color=colors[i+1], lw=2, linestyle="-", label='damping={}'.format(damp[i]))

#P.title('Estimated DESI kappa, nside={}, rtmax=100'.format(nside))
P.ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
P.xlabel(r'$\ell$', fontsize=18)
P.xlim([0, 800])
P.yscale('log')
P.legend(loc='lower right', fancybox=True, framealpha=0.5)
P.annotate(r'DESI noiseless, nside={}, rtmax=rpmax={}, sradius={}'.format(nside,rtmax,srad), xy=(0.05, 0.95), xycoords='axes fraction')

