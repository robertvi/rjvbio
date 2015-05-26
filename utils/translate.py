#!/usr/bin/python

'''
translate nucleotide sequences into protein sequences 
'''

import sys,argparse
import rjvbio.seq
import Bio.SeqIO,Bio.SeqRecord
from Bio.Seq import translate, Seq

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',nargs='+',required=True,type=str,help='input sequence file(s)')
ap.add_argument('--inpformat',default='fasta',type=str,help='format of input file(s), eg fasta,fastq,genbank see http://biopython.org/wiki/SeqIO#File_Formats for details')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args()

if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')
    
for fname in conf.inp:
    for rec in Bio.SeqIO.parse(fname,conf.inpformat):
        seq = translate(rec.seq)
        newrec = Bio.SeqRecord.SeqRecord(seq, id=rec.id,description='')
        Bio.SeqIO.write(newrec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
