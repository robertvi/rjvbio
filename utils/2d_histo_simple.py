#!/usr/bin/python

'''
plot a 2D histogram using matplotlib
from a file with two columns of paired data
http://matplotlib.org/examples/pylab_examples/hist2d_log_demo.html
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='data file, two values per line')
ap.add_argument('--bins',default=40,type=int,help='bins')
conf = ap.parse_args()

from matplotlib.colors import LogNorm
from pylab import *

x = []
y = []

f = open(conf.inp)
for line in f:
    if line.startswith('#'): continue
    
    tok = line.strip().split()
    if len(tok) != 2: continue
    
    try:
        xval = float(tok[0])
        yval = float(tok[1])
    except:
        continue
        
    x.append(xval)
    y.append(yval)

f.close()

hist2d(x, y, bins=conf.bins, norm=LogNorm())
colorbar()
show()
