#!/usr/bin/python

'''
given corresponding per-base coverage and percent gc files
make a 2D histogram data structure and save to file
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpgc',required=True,type=str,help='input percent GC values, assumed to be floats in the range 0...100, and where < 0 means undefined')
ap.add_argument('--inpcov',required=True,type=str,help='input coverage values, integers >= 0')
ap.add_argument('--gcbins',default=100,type=int,help='how many bins to divide GC into')
ap.add_argument('--maxcov',default=10000,type=int,help='coverage greater than this will be clipped')
ap.add_argument('--covbins',default=1000,type=int,help='input coverage values')
ap.add_argument('--out',required=True,type=str,help='save 2D histogram matrix to this file')
conf = ap.parse_args()

import sys
import numpy as np

cts = np.zeros([conf.gcbins,conf.covbins],dtype=np.int64)

fgc = open(conf.inpgc)
fcov = open(conf.inpcov)
linect = 0

while True:
    gc = fgc.readline()
    cov = fcov.readline()
    
    if gc == '' or cov == '': break #end of file
    
    tokgc = gc.strip().split()
    tokcov = cov.strip().split()

    assert tokgc[0] == tokcov[0], "tokgc[0] %s, tokcov[0] %s"%(tokgc[0],tokcov[0])
    assert tokgc[1] == tokcov[1], "tokgc[1] %s, tokcov[1] %s"%(tokgc[1],tokcov[1])

    if float(tokgc[2]) < 0.0: continue #skip where gc is undefined

    #map 0.0 -> 0, maxcov -> covbins-1
    x = int(round(float(tokcov[2]) / conf.maxcov * (conf.covbins - 1.0)))
    if x > conf.covbins-1: x = conf.covbins - 1

    #map 0.0->0, maxgc->gcbins-1
    y = int(round(float(tokgc[2]) / 100.0 * (conf.gcbins - 1.0)))
    if y > conf.gcbins-1: y = conf.gcbins - 1

    cts[y][x] += 1
    linect += 1
    
    if linect % 10000000 == 0:
        print linect
        np.savetxt(conf.out,cts,fmt='%d')
    #print x,y,cts[x][y]

np.savetxt(conf.out,cts,fmt='%d')
