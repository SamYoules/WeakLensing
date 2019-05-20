##- by Anze, 6/4/19
##- modified by Sam, 16/5/19

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
plt.rcParams.update({'font.size':14})

##- fraction of sky covered by survey footprint
fsky=0.36

##- variables for selecting and shaping input data
Navg=7    ##- how many ells to average over
Noffs=34  ##- offset for ells

##- Get number of ells, modes and degrees of freedom
full_ell=np.arange(128*3)[Noffs:]            ##- remove first 34 values in ells
nmodes=fsky*(2*full_ell+1).reshape((-1,Navg)).sum(axis=1)
ell=full_ell.reshape((-1,Navg)).mean(axis=1)   ##- 50 ell values from 37 to 380
dof=len(ell)

##- Get the Cls for auto, cross and input
dir = 'kappa_opt_srad/cls/'
fin = [dir+'cl_auto_cut_128_100_01_01.txt', dir+'cl_cross_cut_128_100_01_01.txt',
                                                         dir+'cl_input_128.txt']
auto,cross,autoin=[loadps(fn) for fn in fin]

##- Get the transfer function (input Cls -> estimated cross-correlation Cls)
##- and fit with a smooth, linear function (to minimize coupling between bins 
##- due to finite sky area)
itrans=cross/autoin
fitfunc_trans=lambda x:x[0]+x[1]*ell#+x[2]*ell*ell
fit=opt.least_squares(lambda x:fitfunc_trans(x)-itrans,[0,0])
transfit=1/fitfunc_trans(fit.x)
plt.figure()
plt.plot(ell,itrans, label='cross/input')
plt.plot(ell,1/transfit, label='transfer function')
plt.xlabel('$\ell$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
plt.legend(loc = "upper left")
print (fit.x)

##- Calculate and plot the noise power spectrum
noise=auto*transfit**2-autoin

##- Fit the noise with another smooth function (const + decaying exponential)
plt.figure()
plt.plot(ell,noise, label='Auto Noise')
fitfunc_noise=lambda x:1e-7*(x[0]+x[1]*np.exp(-ell/x[2]))
fit=opt.least_squares(lambda x:1e8*(fitfunc_noise(x)-noise),[0.3,1.2,100])
noisefit=fitfunc_noise(fit.x)
plt.plot(ell,noisefit,'y-')
plt.xlabel('$\ell$')
plt.ylabel('$C_\ell^{\kappa \kappa}$')
plt.legend(loc = "upper right")
print (fit.x)

##- plot the errors on individual measurements (auto)
auto_err=np.sqrt(2*(autoin+noisefit)**2/nmodes)
auto_corr=auto*transfit**2-noisefit ## corr here means corrected for transfer function and noise bias
plt.figure()
plt.errorbar(ell,auto_corr,yerr=auto_err,fmt='bo', label='Auto')
plt.plot(ell,autoin,'r-')
plt.xlabel('$\ell$')
plt.ylabel('$C_\ell^{\kappa \kappa}$')
chi2=(((auto_corr-autoin)/auto_err)**2).sum()
print(dof)
print ("chi2 =",chi2, "dof=",dof,"PTE=",1-schi2(df=dof).cdf(chi2))
print ("chi2 of null model =",((auto_corr/auto_err)**2).sum())
chi2text1 = "chi2="+str(chi2)
chi2text2 = "dof="+str(dof)
chi2text3 = "PTE="+str(1-schi2(df=dof).cdf(chi2))
chi2text4 = "chi2 of null model ="+str(((auto_corr/auto_err)**2).sum())
plt.text(150,5.2e-8, chi2text1)
plt.text(150,4.7e-8, chi2text2)
plt.text(150,4.2e-8, chi2text3)
plt.text(150,3.7e-8, chi2text4)
plt.legend(loc = "upper right")

##- plot the errors on individual measurements (cross)
cross_err=np.sqrt(((autoin+noisefit)*autoin+autoin**2)/nmodes)/2 ## i don't understand the factor of 2 here
cross_corr=cross*transfit
plt.figure()
plt.errorbar(ell,cross_corr,yerr=cross_err,fmt='bo', label='Cross')
plt.plot(ell,autoin,'r-')
plt.xlabel('$\ell$')
plt.ylabel('$C_\ell^{\kappa \kappa}$')
chi2=(((cross_corr-autoin)/cross_err)**2).sum()
print ("chi2 =",chi2, "dof=",dof,"PTE=",1-schi2(df=dof).cdf(chi2))
print ("chi2 of null model =",((cross_corr/cross_err)**2).sum())
chi2text1 = "chi2="+str(chi2)
chi2text2 = "dof="+str(dof)
chi2text3 = "PTE="+str(1-schi2(df=dof).cdf(chi2))
chi2text4 = "chi2 of null model ="+str(((cross_corr/cross_err)**2).sum())
plt.text(150,5.2e-8, chi2text1)
plt.text(150,4.7e-8, chi2text2)
plt.text(150,4.2e-8, chi2text3)
plt.text(150,3.7e-8, chi2text4)
plt.legend(loc = "upper right")


