#!/usr/bin/python

'''
Read a FASTA and a GFF file, output the sequence of all mRNAs to a FASTA file
output as protein or transcript
run check_gff.py on the GFF file first to check for errors
'''

import sys,argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq,translate
from BCBio import GFF
from rjvbio.seq import parse_gff

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpfasta',required=True,type=str,help='input FASTA file')
ap.add_argument('--inpgff',required=True,type=str,help='input GFF file')
ap.add_argument('--protein',action='store_true',help='if true output protein else output transcript')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args() #sys.argv

#read in vesca pseudomolecules as a dictionary name:seqrecord
seq_dict = SeqIO.to_dict(SeqIO.parse(conf.inpfasta,"fasta"))

data = parse_gff(conf.inpgff)

#open output file
if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

if conf.protein:
    required_type = 'CDS'
else:
    required_type = 'exon'

for rec in data.items:
    #only use mRNA items
    if rec.type != 'mRNA': continue
    part_list = []
    
    #get a list of its exons or CDSs
    for kid in data.kids[rec.id]:
        sub = data.ids[kid]
        assert sub.strand == rec.strand
        if sub.type == required_type: part_list.append(sub)

    #sort into order by start bp
    part_list.sort(key=lambda sub:sub.start)
    
    #assemble exons
    seq = ''
    for sub in part_list:
        seq += str(seq_dict[rec.seqid].seq[sub.start:sub.end])
        
    seq = Seq(seq)
    if rec.strand == '-': seq = seq.reverse_complement()
    if conf.protein: seq = translate(seq)
    newrec = SeqRecord(seq, id=rec.id,description='')

    #write out to file
    SeqIO.write(newrec,fout,"fasta")

if conf.out != 'STDOUT':
    fout.close()
