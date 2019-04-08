import glob
import os
from shutil import copyfile

outdir = '../deltas_half/'

fi = glob.glob("*.fits.gz")
for i,f in enumerate(fi):
    if i%2 != 0:
        copyfile(f, outdir+f)

