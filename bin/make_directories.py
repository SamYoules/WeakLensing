import os

## to make directories
for i in range(1,101):
    os.makedirs('deltas_lensed{}'.format(i))

## to remove files
for i in range(1,34):
    os.remove('do_xkappa{}_errfile'.format(i))

## to remove empty directories
for i in range(1,101):
    os.rmdir('deltas_lensed{}'.format(i))

## to remove directories + contents	
import shutil
for i in range(1,101):
    shutil.rmtree('deltas_lensed{}'.format(i))
