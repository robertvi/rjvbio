#!/usr/bin/python

'''
merge all sequences in a fasta file into a single pseudomolecule
separate original sequences by a run of Ns of the desired length
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input fasta file')
ap.add_argument('--spacer',default=5000,type=int,help='how many Ns to use when joining adjacent sequences together')
ap.add_argument('--linewidth',default=80,type=int,help='bases per line before wrapping')
ap.add_argument('--newseqid',default=None,type=str,help='sequence id to give to the fused sequence (default: input filename)')
ap.add_argument('--out',default='STDOUT',type=str,help='output fasta file')
conf = ap.parse_args()

import sys,os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq,translate

#open output file
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

if conf.newseqid == None: conf.newseqid = os.path.basename(conf.inp)

fout.write('>' + conf.newseqid + '\n')

class seq_adder:
    def __init__(self,fout,wrap,spacing):
        self.fout = fout
        self.wrap = wrap
        self.spacing = spacing
        self.seq = ''
        
    def add_seq(self,line):
        self.seq += line
        
        while len(self.seq) >= self.wrap:
            fout.write(self.seq[:self.wrap] + '\n')
            self.seq = self.seq[self.wrap:]
        
    def add_spacer(self):
        self.add_seq('N' * self.spacing)
        
    def finish(self):
        fout.write(self.seq + '\n')

sa = seq_adder(fout,conf.linewidth,conf.spacer)

f = open(conf.inp)
for i,line in enumerate(f):
    if line.startswith('>'):
        if i == 0: continue #first sequence does not need a space prepending to it
        sa.add_spacer()
        continue
        
    sa.add_seq(line.strip())
    
sa.finish()
    
if conf.out != 'STDOUT': fout.close()
