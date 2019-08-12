from astropy.table import Table, join

zcat='zcats/zcat_desi_drq_unlensed.fits'
t = Table.read(zcat)
t['DEC'] = t['DEC0']
t['RA'] = t['RA0']
t.write(zcat, format='fits', overwrite=True)

