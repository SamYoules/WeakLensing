## SY 27/2/19
## Plots power spectra of true correlation for both map-making methods for rtmax 40, 70 and 100.

from astropy.io import fits
import healpy as hp
import scipy as sp
import numpy as np
import pylab as P
import matplotlib.gridspec as gridspec

def Get_Cls(fin):
    da=np.load(fin); Nside=int(da['arr_0']); ids=da['arr_1']; mp=da['arr_2'][0];
    k=np.zeros(Nside**2*12)
    k[ids]=mp
    #mask = k!=0
    #mask &= (k>np.percentile(k[mask], 0.5)) & (k<np.percentile(k[mask], 99.5))
    #k[~mask]=hp.UNSEEN
    Cls=hp.sphtfunc.anafast(k, lmax=3*Nside-1)
    return Cls

##- Open midpoint kappa files
k_input = fits.open('est_maps_noiseless/kappa_input1.fits')[1].data.I
wk1 = fits.open('kappa_midpoints/kappa-noisy-lensed-truecorr-rt40-rp100.fits.gz')[1].data.wkappa
k1 = fits.open('kappa_midpoints/kappa-noisy-lensed-truecorr-rt40-rp100.fits.gz')[1].data.kappa
k2 = fits.open('kappa_midpoints/kappa-noisy-lensed-truecorr-rt70-rp100.fits.gz')[1].data.kappa
k3 = fits.open('kappa_midpoints/kappa-noisy-lensed-truecorr-rt100-rp100.fits.gz')[1].data.kappa
k4 = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt40-rp100.fits.gz')[1].data.kappa
k5 = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt70-rp100.fits.gz')[1].data.kappa
k6 = fits.open('kappa_midpoints/kappa-cut-lensed-truecorr-rt100-rp100.fits.gz')[1].data.kappa

##- Get resolution of estimated map
NSIDE=int(hp.npix2nside(k1.size))

##- Mask off areas outside DESI footprint
mask = wk1!=0
mask &= (k1>np.percentile(k1[mask], 0.5)) & \
                (k1<np.percentile(k1[mask], 99.5))
k1[~mask]=hp.UNSEEN
k2[~mask]=hp.UNSEEN
k3[~mask]=hp.UNSEEN
k4[~mask]=hp.UNSEEN
k5[~mask]=hp.UNSEEN
k6[~mask]=hp.UNSEEN

##- Reset resolution of input map to match estimated maps
alm = hp.sphtfunc.map2alm(k_input, lmax=3*NSIDE-1)
k_input = hp.sphtfunc.alm2map(alm, nside=NSIDE)

##- Mask area outside footprint
k_input = k_input*(mask)+hp.UNSEEN*(~mask) 
k_input[~mask]=hp.UNSEEN

##- Get input Cls
Cl_input = hp.sphtfunc.anafast(k_input, lmax=3*NSIDE-1)

# Get the C_ells for auto-correlations
Cls0=hp.sphtfunc.anafast(k1, lmax=3*NSIDE-1)
Cls1=hp.sphtfunc.anafast(k2, lmax=3*NSIDE-1)
Cls2=hp.sphtfunc.anafast(k3, lmax=3*NSIDE-1)
Cls3=hp.sphtfunc.anafast(k4, lmax=3*NSIDE-1)
Cls4=hp.sphtfunc.anafast(k5, lmax=3*NSIDE-1)
Cls5=hp.sphtfunc.anafast(k6, lmax=3*NSIDE-1)


## Load optimal kappa files
opt_kappas=['kappa_opt/kappa_cut_rt40_rp100.npz','kappa_opt/kappa_cut_rt70_rp100.npz','kappa_opt/kappa_cut_rt100_rp100.npz','kappa_opt/kappa_noisy_rt40_rp100.npz','kappa_opt/kappa_noisy_rt70_rp100.npz','kappa_opt/kappa_noisy_rt100_rp100.npz']

Cls=[]
for i in opt_kappas:
    Cls.append(Get_Cls(i))

ell = np.arange(Cls0.size)

##- Setup figures
P.rcParams.update({'font.size':14})
P.ion()
# Create 2x1 sub plots and set right-hand axis labels
ncolors=9
colors = P.cm.Set1(np.linspace(0,1,ncolors))

gs = gridspec.GridSpec(3, 1, hspace = 0.2, wspace=0.3)
fig = P.figure(figsize=(7,10))
ax1 = P.subplot(gs[1, 0]) # row 0, col 0
ax1.text(1.02, 0.5, "Noisy",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)

ax2 = P.subplot(gs[2, 0]) # row 1, col 0
ax2.text(1.02, 0.5, "Noiseless",
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=P.gca().transAxes)


##- Subplot 1
ax1.plot(Cl_input*ell, color=colors[3], lw=2, linestyle="-",label='Input')

ax1.plot(Cls0*ell, color=colors[0], lw=2, linestyle="-", label='Midpoint, rt_max = 40')
ax1.plot(Cls[0]*ell, color=colors[0], lw=2, linestyle="--", label='Optimal, rt_max = 40')
ax1.plot(Cls1*ell, color=colors[1], lw=2, linestyle="-", label='Midpoint, rt_max = 70')
ax1.plot(Cls[1]*ell, color=colors[1], lw=2, linestyle="--", label='Optimal, rt_max = 70')
ax1.plot(Cls2*ell, color=colors[2], lw=2, linestyle="-", label='Midpoint, rt_max = 100')
ax1.plot(Cls[2]*ell, color=colors[2], lw=2, linestyle="--", label='Optimal, rt_max = 100')

ax1.axhline(0., color='k', ls=':')
ax1.set_ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
#ax1.set_xlabel(r'$\ell$', fontsize=18)
ax1.set_xlim([0, 800])
#ax1.set_ylim([0.,1.2e-4])
ax1.set_yscale('log')
#ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)

#handles, labels = ax1.get_legend_handles_labels()
#order = [1,2,3,0]
#ax1.legend([handles[idx] for idx in order],[labels[idx] for idx in order], bbox_to_anchor=(0.5, 1.5), loc="upper center", borderaxespad=0., fontsize=14, ncol=3, handletextpad=0.5, handlelength=1, columnspacing=1)

#ax1.legend()
ax1.legend(bbox_to_anchor=(0.5, 2.6), loc="upper center", borderaxespad=0., fontsize=14, ncol=1, handletextpad=0.0, handlelength=2, columnspacing=1)

##- Subplot 2
ax2.plot(Cl_input*ell, color=colors[3], lw=2, linestyle="-",label='Input')

ax2.plot(Cls3*ell, color=colors[0], lw=2, linestyle="-", label='Midpoint, rt_max = 40')
ax2.plot(Cls[3]*ell, color=colors[0], lw=2, linestyle="--", label='Optimal, rt_max = 40')
ax2.plot(Cls4*ell, color=colors[1], lw=2, linestyle="-", label='Midpoint, rt_max = 70')
ax2.plot(Cls[4]*ell, color=colors[1], lw=2, linestyle="--", label='Optimal, rt_max = 70')
ax2.plot(Cls5*ell, color=colors[2], lw=2, linestyle="-", label='Midpoint, rt_max = 100')
ax2.plot(Cls[5]*ell, color=colors[2], lw=2, linestyle="--", label='Optimal, rt_max = 100')

ax2.axhline(0., color='k', ls=':')
ax2.set_ylabel(r'$\ell \ C_{\ell}^{\rm{input, est}}$', fontsize=18)
ax2.set_xlabel(r'$\ell$', fontsize=18)
ax2.set_xlim([0, 800])
#ax2.set_ylim([0.,5e-4])
ax2.set_yscale('log')
#ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)

#handles, labels = ax1.get_legend_handles_labels()
#order = [1,2,3,0]
#ax2.legend([handles[idx] for idx in order],[labels[idx] for idx in order], bbox_to_anchor=(0.5, 1.5), loc="upper center", borderaxespad=0., fontsize=14, ncol=3, handletextpad=0.5, handlelength=1, columnspacing=1)

#ax2.legend()

P.savefig('plots/compare_maps.pdf', bbox_inches='tight')

