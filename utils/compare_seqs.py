#!/usr/bin/python

'''
compare two fasta files
report any differences
based on seqids not the ordering in the file
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp1',required=True,type=str,help='input FASTA file 1')
ap.add_argument('--inp2',required=True,type=str,help='input FASTA file 2')
conf = ap.parse_args()

import pyfaidx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq,translate

fasta1 = pyfaidx.Fasta(conf.inp1)
fasta2 = pyfaidx.Fasta(conf.inp2)

for rec in fasta1.records.iterkeys():
    if not rec in fasta2:
        print conf.inp2,'missing',rec

for rec in fasta2.records.iterkeys():
    if not rec in fasta1:
        print conf.inp1,'missing',rec

for rec in fasta1.records.iterkeys():
    if not rec in fasta2: continue
    if fasta1[rec][:].seq != fasta2[rec][:].seq:
        print 'different sequence',rec
        print fasta1[rec][:].seq
        print
        print fasta2[rec][:].seq
