#!/usr/bin/python

'''
plot a 2D histogram using matplotlib
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input percent GC values, assumed to be in the range 0...100, and where < 0 means undefined')
conf = ap.parse_args()

import sys
from collections import defaultdict

cts = defaultdict(defaultdict(int))

fgc = open(conf.inpgc)
fcov = open(conf.inpcov)

while True:
    gc = fgc.readline()
    cov = fcov.readline()
    
    if gc == '' or cov == '': break #end of file
    
    tokgc = gc.strip().split()
    tokcov = cov.strip().split()

    assert tokgc[0] == tokcov[0]
    assert tokgc[1] == tokcov[1]

    if float(tokgc[2]) < 0.0: continue #skp where gc is undefined

    x = int(tokcov[2])
    y = int(float(tokgc[2])*10.0)
