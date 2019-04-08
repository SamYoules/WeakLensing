## SY 25/3/19
## Plots power spectra for optimal map-making method with different search radii.
## Maps made with get_map.py, cls made with write_spectrum.py

import numpy as np
import pylab as P

nside   = 64
rtmax   = 100
sradius = ['01', '05', '10']
srad    = [0.1, 0.5, 1.0]
damping = ['_01', '_05', '_10']
damp    = [0.1, 0.5, 1.0]

##- Setup figure
P.ion()
fig, [ax1, ax2, ax3] = P.subplots(1, 3, sharex=True, figsize=(16,7))
fig.suptitle('DESI noiseless, nside={}, rtmax=rpmax={}, xcf coefficient'.format(nside,rtmax))
#fig.rcParams.update({'font.size':14})
#fig.ion()
#P.figure()
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

##- Get Input Cls
cl_input  = np.loadtxt('kappa_opt_srad/cls/cl_input_{}.txt'.format(nside))

##- Get Optimal Cls
for i,j in enumerate(sradius):
    cla  = np.loadtxt('kappa_opt_srad/cls/cl_auto_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[0]))
    clx  = np.loadtxt('kappa_opt_srad/cls/cl_cross_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[0]))
    cl   = clx / np.sqrt(cla*cl_input)
    ell = np.arange(cl.size)
    ax1.plot(ell, cl, color=colors[i+1], lw=2, linestyle="-", label='sradius={}'.format(srad[i]))
    ax1.set_xlabel(r'$\ell$', fontsize=18)
    ax1.set_xlim([0, 200])
    #ax1.set_ylim([10e-9, 10e-6])
    #ax1.set_yscale('log')
    ax1.annotate(r'damping={}'.format(damp[0]), xy= (0.05, 0.95), xycoords='axes fraction')
    ax1.set_ylabel(r'$C_{\ell}^{\rm{t,e}}/\sqrt{C_{\ell}^{\rm{t,t}}C_{\ell}^{\rm{e,e}}}}$', fontsize=18)
    ax1.legend(loc='lower right', fancybox=True, framealpha=0.5)

for i,j in enumerate(sradius):
    cla  = np.loadtxt('kappa_opt_srad/cls/cl_auto_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[1]))
    clx  = np.loadtxt('kappa_opt_srad/cls/cl_cross_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[1]))
    cl   = clx / np.sqrt(cla*cl_input)
    ell = np.arange(cl.size)
    ax2.plot(ell, cl, color=colors[i+1], lw=2, linestyle="-", label='sradius={}'.format(srad[i]))
    ax2.set_xlabel(r'$\ell$', fontsize=18)
    ax2.set_xlim([0, 200])
    #ax2.set_ylim([10e-9, 10e-6])
    #ax2.set_yscale('log')
    ax2.annotate(r'damping={}'.format(damp[1]), xy= (0.05, 0.95), xycoords='axes fraction')

for i,j in enumerate(sradius):
    cla  = np.loadtxt('kappa_opt_srad/cls/cl_auto_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[2]))
    clx  = np.loadtxt('kappa_opt_srad/cls/cl_cross_cut_{}_{}_{}{}.txt'.format(nside,rtmax,j,damping[2]))
    cl   = clx / np.sqrt(cla*cl_input)
    ell = np.arange(cl.size)
    ax3.plot(ell, cl, color=colors[i+1], lw=2, linestyle="-", label='sradius={}'.format(srad[i]))
    ax3.set_xlabel(r'$\ell$', fontsize=18)
    ax3.set_xlim([0, 200])
    #ax3.set_ylim([10e-9, 10e-6])
    #ax3.set_yscale('log')
    ax3.annotate(r'damping={}'.format(damp[2]), xy= (0.05, 0.95), xycoords='axes fraction')

