#!/usr/bin/python

'''
split the input sequences whereever a run of NNNs is found
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input FASTA file')
ap.add_argument('--minNs',default=1,type=int,help='minimum number of Ns to trigger splitting')
ap.add_argument('--minlen',default=1,type=int,help='minimum length of sequence to output')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args()

import rjvbio.seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO,Bio.SeqRecord
from Bio.Seq import translate, Seq
import re

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')
    
for rec in Bio.SeqIO.parse(conf.inp,'fasta'):
    seq = str(rec.seq)
    m = [[x.start(),x.end()] for x in re.finditer('N{%d,}'%conf.minNs,seq.upper())]
    for i in xrange(len(m)+1):
        if i == 0 and len(m) == 0:
            start = 0
            end = len(seq)
        elif i == 0:
            start = 0
            end = m[i][0]
        elif i == len(m):
            start = m[i-1][1]
            end = len(seq)
        else:
            start = m[i-1][1]
            end = m[i][0]
        
        newseq = Seq(seq[start:end])
        if len(newseq) < conf.minlen: continue
        
        newrec = SeqRecord(newseq, id=rec.id+'_'+str(i+1), description='')

        #write out to file
        Bio.SeqIO.write(newrec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
