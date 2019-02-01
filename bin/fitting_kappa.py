## Author: Ben Mawdsley
## Jan 2019
## Called by testing_fitting.py
## Forward fits a kappa map by perturbing random l-modes, creating new kappa map and testing chi2 of new map. Rinse & repeat until achieve 500 iterations with no improvement of chi2.

import numpy as np
import healpy as hp
import os

def chi2(f1, f2, e):
#calculates chi2
    e[np.where(e<0)]=0.0
    aa=np.where(e!=0)
    c2=np.sum((f1[aa]-f2[aa])**2/e[aa]**2)
    return c2

def ealm2kappamap(ealm, n_side):
#makes a kappa map from e mode coefficients
    
    ells=hp.sphtfunc.Alm.getsize(lmax=2*n_side-1)
    lmode=hp.sphtfunc.Alm.getlm(lmax=2*n_side-1)[0]
    lfacskappa=(lmode*(lmode+1.)*((lmode+2.)**-1)*((lmode-1.)**-1))**0.5
    lfacskappa[np.where(lmode<3)]=0.0
#    print(lmode)
#    print(lfacskappa)
#    quit()
    lfacskappa[np.isinf(lfacskappa)]=0
    lfacskappa[np.isnan(lfacskappa)]=0
    kappa_alm=lfacskappa*ealm
    kmap=hp.sphtfunc.alm2map(kappa_alm, nside=n_side, lmax=2*n_side-1)
    return kmap 


def make_guess(harm_coeff,kick_size,  L_MAX):
#perturbs coefficients of a randomly chosen l mode
    ellchange=np.random.randint(2, L_MAX+1)
    ell_coords=hp.sphtfunc.Alm.getlm(lmax=L_MAX)[0]
    kicks=np.random.normal(0.0, kick_size, ell_coords[np.where(ell_coords==ellchange)].shape)
    hyp=harm_coeff.copy()
    locs=np.where(ell_coords==ellchange)
    hyp[locs]= harm_coeff[locs]+(kicks*harm_coeff[locs])
#    print('hyp changed',ellchange , hyp-harm_coeff, hyp[locs]-harm_coeff[locs])
    return hyp

def guess2obs(harmonics, n_side, field='kappa'):
#a function to take e mode coefficients and produce a map to compare against observables. Will edit this to have a shear option too
    obs_field=0.0
    if field =='kappa' :
#        print('making kappa map')
        obs_field=ealm2kappamap(harmonics, n_side)

    return obs_field


def fit_map(start_guess, data_fitting, error_field):
#Over-arching fitting loop. Should be the only routine you need to call to produce a fitted map
    NSIDE=hp.pixelfunc.get_nside(data_fitting)
    start_map=ealm2kappamap(start_guess, NSIDE)
    start_chi=chi2(start_map, data_fitting, error_field)
    print('gogogo', start_map)
    best_guess=start_guess.copy()
    chi_ref=start_chi*1.0 #first chi2 that we need to improve upon
    dof=np.sum(error_field!=0.0) #number of pixels fitting to
    chi_target=dof # want a reduced chi2 =1
    kk=1.0 #amplitude of kicks to give new hypotheses
    n_loops=0 # iterator measuring number of attempts to improve hypothesis without improvement
    n_quit=500 #maximum number of times to try a new guess without improvement before giving up and outputting
    print(start_chi, dof, hp.nside2npix(NSIDE))
    allchi=[]
    n_all=0
    while chi_target < chi_ref:
        new_guess=make_guess(best_guess, kk, 2*NSIDE-1)
        ab=new_guess-start_guess
        new_map=guess2obs(new_guess, NSIDE)
        new_chi2=chi2(new_map, data_fitting, error_field)
        print('new chi2', new_chi2)
        if new_chi2 < chi_ref:
            best_guess=new_guess.copy()
            chi_ref=new_chi2
            n_loops=0
        if n_loops>n_quit :
            chi_ref=chi_ref*-1.
        if n_loops > 20:
            kk=kk*0.99
        if n_loops >400:
            kk=kk*1.1
        n_loops+=1
        n_all+=1
        print(n_loops, chi_ref, new_chi2)
        allchi.append(chi_ref/dof)
        np.savetxt('RunningChi_step'+str(n_all)+'.txt', allchi)
        if n_all>1:        
            os.remove('RunningChi_step'+str(n_all-1)+'.txt')

    return best_guess







