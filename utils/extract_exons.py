#!/usr/bin/python

'''
extract exons from an annotated genome
extract just genes listed in an id file
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpfasta',required=True,type=str,help='input FASTA file of genome assembly (masked or unmasked)')
ap.add_argument('--inpdb',required=True,type=str,help='input gffutils DB of genome annotation')
ap.add_argument('--inpids',required=True,type=str,help='file containing ids of gene to extract exons from')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file containing exonic sequences')
conf = ap.parse_args() #sys.argv

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq,translate
import gffutils
from pyfaidx import Fasta

#open output file
if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

db = gffutils.FeatureDB(conf.inpdb)
fasta = Fasta(conf.inpfasta)

f = open(conf.inpids)
for line in f:
    uid = line.strip()
    
    for child in db.children(uid,order_by='start'):
        if child.featuretype != 'exon': continue
        seq = fasta[child.seqid][child.start:child.end].seq
        #print child.id,child.start,child.end
        #print seq
        #print

        seq = Seq(seq)
        newrec = SeqRecord(seq, id=child.id,description=uid)

        #write out to file
        SeqIO.write(newrec,fout,"fasta")
        
f.close()

if conf.out != 'STDOUT':
    fout.close()
