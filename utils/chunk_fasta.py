#!/usr/bin/python

'''
output one chunk (subset) of a fasta which must be samtools-indexed already
eg:

chunk_fasta.py --inp INPUT.fa --chunks 10 --thischunk 1 | other_script.sh

'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input fasta file')
ap.add_argument('--chunks',required=True,type=int,help='how many chunks to split the file into')
ap.add_argument('--thischunk',required=True,type=int,help='which chunk to output, 1-based')
conf = ap.parse_args()

#from pyfaidx import Fasta
#import pysam,
from rjvbio.seq import fasta
import math,os

#index must already exist (create with samtool faidx / pysam / pyfaidx etc)
#do not make the index on-the-fly as we expect to be accessing the file in parallel
#with other processes
assert os.path.isfile(conf.inp+'.fai')
assert conf.thischunk <= conf.chunks

#get list of records from index file
seqs = fasta(conf.inp)

records = len(seqs)

per_chunk = float(records) / float(conf.chunks)

if conf.thischunk == 1:
    first = 0
else:
    first = int(round(float(conf.thischunk-1)*per_chunk))
    
if conf.thischunk == conf.chunks:
    last = records
else:
    last = int(round(float(conf.thischunk)*per_chunk))

for i in xrange(first,last):
    print '>' + seqs[i].uid
    for line in seqs.iterseq(i):
        print line
