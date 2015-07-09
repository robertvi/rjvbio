#!/usr/bin/python

'''
Read a FASTA and a gffutils database file, output the sequence of all mRNAs to a FASTA file
output as protein or transcript
this version uses gffutils so that it can handle gtf as well as gff3
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpfasta',required=True,type=str,help='input FASTA file')
ap.add_argument('--gffdb',required=True,type=str,help='input gffutils database')
ap.add_argument('--featuretype',default='mRNA',type=str,help='featuretype denoting transcripts')
ap.add_argument('--subfeaturetype',default=None,type=str,help='featuretype denoting exon or CDS')
ap.add_argument('--protein',action='store_true',help='if true output protein else output transcript')
ap.add_argument('--excludestop',action='store_true',help='if true remove final * from proteins')
ap.add_argument('--minlen',default=1,type=int,help='discard transcripts below this length')
ap.add_argument('--minunmasked',default=1,type=int,help='discard transcripts below this number of unmasked (non N/X) bases/residues')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args() #sys.argv

import gffutils,pyfaidx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq,translate

#genome fasta file
fasta = pyfaidx.Fasta(conf.inpfasta)

#gffutils database
db = gffutils.FeatureDB(conf.gffdb)

#open output file
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

if conf.protein:
    if conf.subfeaturetype == None:
        conf.subfeaturetype = 'CDS'
else:
    if conf.subfeaturetype == None:
        conf.subfeaturetype = 'exon'

#get all transcript items
for rec in db.features_of_type(conf.featuretype,order_by=['seqid','start']):
    #get list of relevant sub features (exon or CDS)
    part_list = [sub for sub in db.children(rec,featuretype=conf.subfeaturetype)]

    #sort into order by start bp
    part_list.sort(key=lambda sub:sub.start)
    
    #assemble exons / CDSs
    seq = ''.join([str(fasta[rec.seqid][sub.start-1:sub.end].seq) for sub in part_list])
        
    seq = Seq(seq,generic_dna)
    if rec.strand == '-' and len(seq) > 0: seq = seq.reverse_complement()
    
    if conf.protein:
        if len(seq) == 0: seq = Seq('',generic_protein)
        else:             seq = translate(seq)
        if len(seq) < conf.minlen: continue
        if len(seq) - seq.count('X') < conf.minunmasked: continue
        if len(seq) > 0 and conf.excludestop and seq[-1] == '*': seq = seq[:-1]
    else:
        if len(seq) < conf.minlen: continue
        if len(seq) - seq.count('N') < conf.minunmasked: continue
            
    newrec = SeqRecord(seq, id=rec.id,description='')

    #write out to file
    SeqIO.write(newrec,fout,"fasta")

if conf.out != 'STDOUT': fout.close()
