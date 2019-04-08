##- by Anze, 6/4/19

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from  scipy.stats import chi2 as schi2

def loadps(fn):
    '''Reads in Cl files for auto, cross and input, then reshapes the data'''
    data=np.loadtxt(fn)[Noffs:]
    data=data.reshape((-1,Navg)).mean(axis=1)
    return data  

##- setup fig
plt.ion()
plt.figure()

##- fraction of sky covered by survey footprint
fsky=0.36

##- variables for selecting and shaping input data
Navg=7
Noffs=34

##- Get number of ells, modes and degrees of freedom
full_ell=np.arange(128*3)[Noffs:]
nmodes=fsky*(2*full_ell+1).reshape((-1,Navg)).sum(axis=1)
ell=full_ell.reshape((-1,Navg)).mean(axis=1)
dof=len(ell)

plotcolours=['#396AB1','#CC2529','#3E9651','#DA7C30','#535154','#6B4C9A','#922428','#948B3D']
barcolours=['#7293CB','#D35E60','#84BA5B','#E1974C','#808585','#9067A7','#AB6857','#CCC210']

nside = 128
rtmax = 100
#srad = '01'
#sradius = 0.1
srad = ['01', '05', '10']
sradius = [0.1, 0.5, 1.0]
#damp = ['01', '05', '10']
#damping = ['0.1', '0.5', '1.0']
damp = '01'
damping = 0.1
dir = 'kappa_opt_srad/cls/'

##- Get midpoint noise
mid = [dir+'cl_auto_cut_midpoint_{}_{}.txt'.format(rtmax,nside), dir+'cl_cross_cut_midpoint_{}_{}.txt'.format(rtmax,nside), dir+'cl_input_{}.txt'.format(nside)]
auto,cross,autoin=[loadps(fn) for fn in mid]
itrans=cross/autoin
fitfunc_trans=lambda x:x[0]+x[1]*ell#+x[2]*ell*ell
fit=opt.least_squares(lambda x:fitfunc_trans(x)-itrans,[0,0])
transfit=1/fitfunc_trans(fit.x)
print (fit.x)
noise=auto*transfit**2-autoin
plt.plot(ell,noise, color=plotcolours[3], lw=2, linestyle="-", label='midpoint')
fitfunc_noise=lambda x:1e-7*(x[0]+x[1]*np.exp(-ell/x[2]))
fit=opt.least_squares(lambda x:1e8*(fitfunc_noise(x)-noise),[0.3,1.2,100])
noisefit=fitfunc_noise(fit.x)
print (fit.x)


#for i,j in enumerate(damp):
#    fsuffix ='_cut_{}_{}_{}_{}.txt'.format(nside, rtmax, srad, j)
for i,j in enumerate(srad):
    fsuffix ='_cut_{}_{}_{}_{}.txt'.format(nside, rtmax, j, damp)
    fin = [dir+'cl_auto{}'.format(fsuffix), dir+'cl_cross{}'.format(fsuffix), dir+'cl_input_{}.txt'.format(nside)]

    ##- Get the Cls
    auto,cross,autoin=[loadps(fn) for fn in fin]

    ##- Get the transfer function (input Cls -> estimated cross-correlation Cls)
    ##- and fit with a smooth, linear function (to minimize coupling between bins 
    ##- due to finite sky area)
    itrans=cross/autoin
    fitfunc_trans=lambda x:x[0]+x[1]*ell#+x[2]*ell*ell
    fit=opt.least_squares(lambda x:fitfunc_trans(x)-itrans,[0,0])
    transfit=1/fitfunc_trans(fit.x)
    print (fit.x)

    ##- Calculate and plot the noise power spectrum
    noise=auto*transfit**2-autoin

    ##- Fit the noise with another smooth function (const + decaying exponential)
    #plt.plot(ell,noise, color=plotcolours[i], lw=2, linestyle="-", label='damp={}'.format(damping[i]))
    plt.plot(ell,noise, color=plotcolours[i], lw=2, linestyle="-", label='srad={}'.format(sradius[i]))
    fitfunc_noise=lambda x:1e-7*(x[0]+x[1]*np.exp(-ell/x[2]))
    fit=opt.least_squares(lambda x:1e8*(fitfunc_noise(x)-noise),[0.3,1.2,100])
    noisefit=fitfunc_noise(fit.x)
    plt.plot(ell,noisefit, color=barcolours[i], lw=2, linestyle="-")
    plt.plot(ell,noisefit, 'y-')
    print (fit.x)

#plt.legend(loc='lower left')
#plt.annotate(r'Noise power spectrum, nside={}, sradius={}'.format(nside, sradius), xy= (0.15, 0.95), xycoords='axes fraction')
plt.legend(loc='center right')
plt.annotate(r'Noise power spectrum, nside={}, damping={}'.format(nside, damping), xy= (0.15, 0.95), xycoords='axes fraction')


