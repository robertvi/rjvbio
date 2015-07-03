#!/usr/bin/python

'''
mask a fasta file at the regions indicted by a GFF file
see also bedtools maskfasta
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpfasta',required=True,type=str,help='input FASTA file')
ap.add_argument('--gffdb',required=True,type=str,help='input GFF database file')
ap.add_argument('--featuretypes',default='ALL',type=str,help='which feature type(s) to mask')
ap.add_argument('--out',default='STDOUT',type=str,help='output masked FASTA file')
conf = ap.parse_args()

import gffutils
import rjvbio.seq
import Bio.SeqIO,Bio.SeqRecord
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')
    
db = gffutils.FeatureDB(conf.gffdb)

if conf.featuretypes == 'ALL': conf.featuretypes = None

for rec in Bio.SeqIO.parse(conf.inpfasta,'fasta'):
    seqid = rec.id.strip()
    seq = MutableSeq(str(rec.seq))
    length = len(seq)

    for feature in db.all_features(limit=[seqid,0,length],completely_within=False,featuretype=conf.featuretypes):
        start = min(feature.start-1,feature.end)
        end = max(feature.start-1,feature.end)
        flength = end - start
        seq[start:end] = 'N' * flength
        assert len(seq) == length

    newrec = Bio.SeqRecord.SeqRecord(seq, id=seqid, description='')
    Bio.SeqIO.write(newrec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
