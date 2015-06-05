#!/usr/bin/python

'''
sample the quality values used to help deduce the encoding system used
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input FASTQ file')
ap.add_argument('--samples',default=100000,type=int,help='how many reads to sample')
ap.add_argument('--verbose',action='store_true',help='print frequency break down instead of min,max')
conf = ap.parse_args()

from rjvbio.seq import generate_fastq

#iterate fastq records loading the raw quality characters
cts = {}
minval=1000
maxval=-1000
for i,rec in enumerate(generate_fastq(conf.inp)):
    if i >= conf.samples: break
    #count frequency of the characters
    for x in rec['qual']:
        if not x in cts: cts[x] = 1
        else:            cts[x] += 1
        c = ord(x)
        minval = min(minval,c)
        maxval = max(maxval,c)
    
#print the results
if conf.verbose:
    results = [[ord(x),cts[x]] for x in cts]
    results.sort(key=lambda x:x[0])
        
    for row in results:
        print row[0],row[1]
else:
    print minval,maxval
