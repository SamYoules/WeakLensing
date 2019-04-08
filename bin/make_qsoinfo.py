# Date: 20 Dec 2018
# Author: Sam Youles
# make_qsoinfo.py
# Create qso.txt and qsoradec.txt files containing: thingid, RA, Dec, z

import glob
from astropy.io import fits
import numpy as N
import sys

DR = sys.argv[1]
#f1 = open('qso.txt', 'w')
#f2 = open('qsoradec.txt', 'w')

allfiles = glob.glob('{}/*.fits.gz'.format(DR))
#print(allfiles)
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
        anze_list.append([tab[i].header['THING_ID'],tab[i].header['RA'], \
              tab[i].header['DEC'],tab[i].header['Z']])
        qso_list.append(tab[i].header['THING_ID'])

    # Counter to see progress
    m +=1
    if (m % 100 == 0):
        print (m)

N.savetxt('qsoradec_cut.txt', anze_list, fmt='%d %.14f %.14f %.14f')
N.savetxt('qso_cut.txt', qso_list, fmt='%d')

