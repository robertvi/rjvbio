#!/usr/bin/python

'''
output a subset of the sequences in a fasta file
sequence ids specified on the command line or in file(s)
'''

import sys,argparse
import rjvbio.seq
import Bio.SeqIO,Bio.SeqRecord
from Bio.Seq import translate, Seq

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input FASTA file')
ap.add_argument('--seqid',required=True,type=str,help='seqid(s)')
ap.add_argument('--start',required=True,type=int,help='starting base, 1-based inclusive')
ap.add_argument('--end',required=True,type=int,help='final base, 1-based inclusive')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args()

if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')
    
for rec in Bio.SeqIO.parse(conf.inp,'fasta'):
    if rec.id == conf.seqid:
        uid = conf.seqid+'::'+str(conf.start)+'-'+str(conf.end)
        seq = rec.seq[conf.start-1:conf.end]
        newrec = Bio.SeqRecord.SeqRecord(seq, id=uid,description='')

        Bio.SeqIO.write(newrec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
