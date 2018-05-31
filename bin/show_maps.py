import healpy as hp
import matplotlib.pyplot as P

a = hp.fitsfunc.read_map('outputs/kmap.fits')
hp.mollview(a, title = 'Kmap Mollweide View')
P.savefig('plots/kmap.png')
P.close()

b = hp.fitsfunc.read_map('outputs/amap.fits')
hp.mollview(b, title = 'Amap Mollweide View')
P.savefig('plots/amap.png')
P.close()

c = hp.fitsfunc.read_map('outputs/pmap.fits')
hp.mollview(c, title = 'Pmap Mollweide View')
P.savefig('plots/pmap.png')
P.close()
