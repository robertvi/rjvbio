#!/usr/bin/python

'''
calcule percent GC content using a sliding window
output as a simple per-base flat file comparable to
the -d output file of bedtools genomecov
'''

import argparse,sys

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',default='STDIN',type=str,help='input FASTA')
ap.add_argument('--window',default=200,type=int,help='sliding window size')
ap.add_argument('--mincalled',default=100,type=int,help='minimum called (non-N) bases to calculate a percent GC value')
ap.add_argument('--out',default='STDOUT',type=str,help='output per-base percent GC file (NOT bed or bedgraph format)')
conf = ap.parse_args()

import sys
from Bio.SeqIO import parse

#open files
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')
if conf.inp == 'STDIN': f = sys.stdin
else:                   f = open(conf.inp)

offset = int(conf.window / 2)

for rec in parse(f,'fasta'):
    seqid = rec.id
    length = len(rec.seq)
    if length < conf.window: continue
    
    for i in xrange(length):
        if i-offset < 0 or i-offset+conf.window > length:
            #window is off either end of the chromosome/scaffold
            percent = -1.0
            
        else:
            seq = str(rec.seq[i-offset:i-offset+conf.window]).upper()
            at = seq.count('A') + seq.count('T')
            gc = seq.count('G') + seq.count('C')
            total = gc + at
            
            if total < conf.mincalled:
                #the sequence contains too few known bases
                percent = -1.0
            else:
                percent = float(gc) / total * 100.0
        
        fout.write("%s %d %.1f\n"%(seqid,i+1,percent))

#close files
if conf.inp != 'STDIN': f.close()
if conf.out != 'STDOUT': fout.close()
