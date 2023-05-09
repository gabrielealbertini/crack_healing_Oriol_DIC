import numpy as np
import os
import glob

files = glob.glob('*.jpg')

for f in files:
    n = int(f.split('frame')[1].split('.jpg')[0])
    new_f = '{}_frame_{:0>5}.jpg'.format(f.split('frame')[0],n)
    os.rename(f,new_f)
