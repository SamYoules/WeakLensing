# SY 15/5/18
# 
# Code lifted from:
#    https://github.com/slosar/wldisplace/blob/master/transforms_check.ipynb
#    by Anze Slosar

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
#%matplotlib inline
from SphericalDiff import *

Nside=256
Npix=Nside**2*12

#get kappa map
kappa=SphericalMap(hp.fitsfunc.read_map('outputs/kmap-blob.fits'))
#kappa=SphericalMap(hp.fitsfunc.read_map('kmap.fits'))

#generate some points to displace
Ndots=1000
theta=np.arccos(np.random.uniform(-np.sin(1),np.sin(1),Ndots))
phi=np.random.uniform(-1,1,Ndots)
plt.figure(figsize=(8,8))
plt.plot(phi,np.cos(theta),'b.')
plt.savefig('plots/mock-qsos-test.png')
plt.close()

thetap,phip=kappa.DisplaceObjects(theta,phi)

## plot displacements, there are issues with projection hence weird at the edges
plt.figure(figsize=(8,8))
for t,p,t2,p2 in zip(np.cos(theta),phi,np.cos(thetap),phip):
    plt.plot([p,p2],[t,t2],'r-')
    plt.title('Displacements using kappa map')
plt.savefig('plots/displacements-kappa-blob.png')
plt.close()

## SY 14/5/18: alternative displacement using bend angles instead of kappa
amap_plus = hp.fitsfunc.read_map('outputs/amap-blob.fits', field=1)   # change in dec
amap_minus = hp.fitsfunc.read_map('outputs/amap-blob.fits', field=0)  # change in RA
ipix=hp.ang2pix(Nside,theta,phi)
theta_sy = []
phi_sy = []
for i, j in enumerate(ipix):
    theta_sy.append(amap_plus[j])
    phi_sy.append(-amap_minus[j])             #24/5/18 converted ra to phi
theta_a = np.array(theta_sy) + theta
phi_a = np.array(phi_sy) + phi

plt.figure(figsize=(8,8))
for t,p,t2,p2 in zip(np.cos(theta),phi,np.cos(theta_a),phi_a):
    plt.plot([p,p2],[t,t2],'r-')
    plt.title('Displacements using alpha map')
plt.savefig('plots/displacements-alpha-blob.png')
plt.close()

# SY: Create output file to compare displacements with kappa and displacements with alpha
fout = open('outputs/displacements-blob.txt', 'w')
print('Differences in theta & phi from using kappa map, and using alpha map', file = fout)
print('theta kappa', 'theta_alpha', 'phi kappa', 'phi alpha', file = fout)

fout_2 = open('outputs/displacements_2-blob.txt', 'w')
print('Differences in displacements from using kappa map, and using alpha map', file = fout_2)
print('d_kappa           d_alpha             diff(rad)        diff(arcsec)    (dk-da)/dk', file = fout_2)

d_kappa = np.sqrt((theta - thetap)**2 + (phi - phip)**2)
d_alpha = np.sqrt((theta - theta_a)**2 + (phi - phi_a)**2)
diff = d_kappa - d_alpha
diff_arcsec = diff / np.pi * 180. * 3600.
for j in range(len(ipix)):
    print(thetap[j], theta_a[j], phip[j], phi_a[j], file = fout)
    print(d_kappa[j], d_alpha[j], diff[j], diff_arcsec[j], diff[j]/d_kappa[j], file = fout_2)
fout.close()
fout_2.close()

