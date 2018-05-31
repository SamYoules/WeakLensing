# move_forests.py
# Author: Sam Youles
# 19/4/18
# Uses amap.fits (created in lensing_maps.py) to add bend angles from a lens to
# the RA and DEC of forests in the delta-{}.fits.gz files (from do_deltas.py)
# Writes textfile of old and new RA and DEC
# user input: run move_forests.py delta-files-directory-name

from astropy.io import fits
import numpy as N
import healpy as hp
import sys
import glob


# User should input directory name after run command
dir = sys.argv[1]

# Create output file for recording old & new RA & DEC
fout = open('outputs/test_lensingRAandDEC.txt', 'w')
print('Old DEC, Old RA, New DEC, New RA, alpha_plus, alpha_minus', file = fout)

# Load bend angles from alpha map into an array
amap_plus = hp.fitsfunc.read_map('outputs/amap.fits', field=1)  # dec
amap_minus = hp.fitsfunc.read_map('outputs/amap.fits', field=0) # RA
NSIDE = 256

# Amend DEC and RA in each of the delta files by the bend angle from alpha map
for filename in glob.glob(dir+'/*.fits.gz'):
    try:
        a = fits.open(filename)
    except IOError:
        print("error reading {}\n".format(filename))
        continue

    ra = a[1].header[22]
    dec = a[1].header[18]

    # Convert RA and DEC to index to find corresponding pixel in amap
    # Note: theta is defined in the range [0,pi] and declination is defined in
    # the range [-pi/2,pi/2].
    index = hp.pixelfunc.ang2pix(NSIDE,-dec + N.pi/2., 2.*N.pi - ra)

    # Add bend angles to dec and ra
    dec += amap_plus[index]
    ra += amap_minus[index]

    # Write old & new RA & DEC to textfile
    print(a[1].header[18], a[1].header[22], dec, ra, amap_plus[index], amap_minus[index], file=fout) 
    
    # Rewrite new delta file with new values
    a[1].header[22] = ra
    a[1].header[18] = dec

    a.writeto('lensed_' + filename)

fout.close()
