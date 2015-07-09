#!/usr/bin/python

'''
output a subset of the sequences in a fasta file
sequence ids specified on the command line or in file(s)
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input FASTA file')
ap.add_argument('--seqid',nargs='*',type=str,help='seqid(s)')
ap.add_argument('--seqidfile',nargs='*',type=str,help='filename(s) containing seqid(s)')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args()

import rjvbio.seq
import Bio.SeqIO,Bio.SeqRecord
from Bio.Seq import translate, Seq

if conf.seqid == None and conf.seqidfile == None:
    print 'no seqids specified'
    exit()

#compile the set of seqids from command line
id_dict = {}
if conf.seqid:
    for seqid in conf.seqid:
        id_dict[seqid] = True
    
#compile the set of seqids from seqid files
if conf.seqidfile:
    for fname in conf.seqidfile:
        f = open(fname)
        for line in f:
            seqid = line.strip()
            id_dict[seqid] = True
        f.close()

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')
    
for rec in Bio.SeqIO.parse(conf.inp,'fasta'):
    if rec.id in id_dict:
        Bio.SeqIO.write(rec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
