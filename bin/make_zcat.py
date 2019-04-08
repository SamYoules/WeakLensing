# SY 13/2/19
# Make zcat with RA0 = RA and DEC) = DEC.

import numpy as np
from astropy.table import Table, join

t = Table.read('zcat_desi_drq_lensed.fits')

t['DEC'] = t['DEC0']
t['RA'] = t['RA0']

t.write('zcat_desi_drq_fake.fits', format='fits')


