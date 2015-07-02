#!/usr/bin/python

'''
convert raw integer values into histogram type counts

awk based alternative one liner:
awk '{ct[$3]=ct[$3]+1} END{for(x in ct)print x,ct[x]}'
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',default='STDIN',type=str,help='raw integer values, one per line')
ap.add_argument('--out',default='STDOUT',type=str,help='output file: value count')
conf = ap.parse_args()

import sys,collections

if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

cts = collections.defaultdict(int)

if conf.inp == 'STDIN':
    f = sys.stdin
else:
    f = open(conf.inp)
    
for line in f:
    cts[int(line.strip())] += 1
    
if conf.inp != 'STDIN':
    f.close()

for key in cts:
    print key,cts[key]

if conf.out != 'STDOUT':
    fout.close()
