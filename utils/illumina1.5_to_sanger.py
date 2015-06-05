#!/usr/bin/python

'''
FASTQ quality values: convert Illumina 1.5 to Sanger 

incomplete script
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input FASTQ file')
ap.add_argument('--out',required=True,type=str,help='output FASTQ file')
conf = ap.parse_args()

from rjvbio.seq import generate_fastq

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

#iterate fastq records loading the raw quality characters
for rec in generate_fastq(conf.inp):
    
    
if conf.out != 'STDOUT': fout.close()
