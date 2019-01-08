# Date: 20 Dec 2018
# Author: Sam Youles
# make_qsoinfo.py
# Create qso.txt and qsoradec.txt files containing: thingid, RA, Dec, z

import glob
from astropy.io import fits
import numpy as N
import sys

DR = sys.argv[1]
f1 = open('qso.txt', 'w')
f2 = open('qsoradec.txt', 'w')

allfiles = glob.glob('{}/*.fits.gz'.format(DR))
print(allfiles)
m = 0
anze_list=[]
qso_list=[]

for f in allfiles:
    try:
        tab = fits.open(f)
    except:
        print ("Couldn't open ", f)
        continue

    for i in range(1, len(tab)):
        anze_list.append([tab[i].header['THING_ID'],tab[i].header['RA'], /
              tab[i].header['DEC'],tab[i].header['Z']])
        qso_list.append([tab[i].header['THING_ID'])

        #thid.append(tab[i].header['THING_ID'])
        #ra.append(tab[i].header['RA'])
        #dec.append(tab[i].header['DEC'])
        #z.append(tab[i].header['Z'])
        #f1.write(thid)
        #f2.write([thid, ra, dec, z])

    # Counter to see progress
    m +=1
    if (m % 100 == 0):
        print (m)

N.savetxt('qsoradec.txt', anze_list, fmt='%d %.16e %.16e %d')
N.savetxt('qso.txt', qso_list, fmt='%d')

