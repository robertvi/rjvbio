#!/usr/bin/python

'''
plot a 2D histogram using matplotlib from a text format matrix of count data
http://matplotlib.org/examples/pylab_examples/colorbar_tick_labelling_demo.html
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='matrix of integer count data')
conf = ap.parse_args()

import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

#cts = np.loadtxt(conf.inp,dtype=np.int64)
cts = np.loadtxt(conf.inp,dtype=np.float32)

minct = np.amin(np.amin(cts))
maxct = np.amax(np.amax(cts))
print minct,maxct

cts += 1
cts = np.log10(cts)
minct = np.amin(np.amin(cts))
maxct = np.amax(np.amax(cts))
cts -= minct
cts /= maxct - minct

minct = np.amin(np.amin(cts))
maxct = np.amax(np.amax(cts))
print minct,maxct
#print fcts

#fcts = np.clip(np.random.randn(1000,1000),-1,1)

fig, ax = plt.subplots()

#cax = ax.imshow(cts, interpolation='nearest', cmap=cm.coolwarm)
cax = ax.imshow(cts, interpolation='nearest')

#cbar = fig.colorbar(cax, ticks=[minct,maxct])

plt.show()
