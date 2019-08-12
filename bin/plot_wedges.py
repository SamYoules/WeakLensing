import picca.wedgize
from astropy.io import fits
import numpy as np
import pylab as plt

d0 = fits.open('xcf-exp-noisy-unlensed.out.gz')[1].data
d1 = fits.open('xcf-exp-noisy.out.gz')[1].data

f, (axs) = plt.subplots(nrows=2, ncols=2, figsize=(12,7))
i=0

for mus in [[0., 0.5], [0.5, 0.8], [0.8, 0.95], [0.95, 1.]]:
        w = picca.wedgize.wedge(rpmin=-100.0, rpmax=100.0, nrp=100, \
                rtmin=3.0, rtmax=100.0, nrt=50, \
                rmin=3.0, rmax=100.0, nr=50, \
                mumin=mus[0], mumax=mus[1], ss=10)
        r, wed0, wedcov0 = w.wedge(d0.DA, d0.CO)
        r, wed1, wedcov1 = w.wedge(d1.DA, d1.CO)

        dwed0 = np.sqrt(np.diag(wedcov0))
        dwed1 = np.sqrt(np.diag(wedcov1))

        axs[i//2][i%2].errorbar(r, wed0*r**2, dwed0*r**2, fmt='o', color='#e74c3c', label='Unlensed')
        axs[i//2][i%2].errorbar(r, wed1*r**2, dwed1*r**2, fmt='o', color='#3498db', label='Lensed')
        #axs[i//2][i%2].plot(r,f*r**power,color='red',linewidth=2, label=os.path.basename(file_fit)[:-3])
        axs[i//2][i%2].set_ylabel(r"$r^{2}\xi(r)$")
        if i//2==1:
            axs[i//2][i%2].set_xlabel(r"$r \, [h^{-1}\, \mathrm{Mpc}]$")
        axs[i//2][i%2].set_title(r'$%.1f < \mu < %.1f$'%(mus[0], mus[1]))
        axs[i//2][i%2].legend(fontsize=10)
        i+=1
plt.tight_layout()
plt.show()




