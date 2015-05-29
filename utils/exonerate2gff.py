#!/usr/bin/python

'''
convert exonerate hits from a protein versus genome search
into GFF annotations defining the location of the CDS / exon in the genome

exonerate should be run including sugar and targetgff in the output something like this:

exonerate --showsugar --showtargetgff --model protein2genome PROTEINS.fa GENOME.fa > OUTPUT_FILE

'''

import sys,argparse
from rjvbio.seq import generate_exonerate

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input exonerate file')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
conf = ap.parse_args() #sys.argv

#open output file
if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

for rec in generate_exonerate(conf.inp):
    #print rec.qid,rec.sid
    for tok in rec.lines:
        fout.write('\t'.join(tok) + '\n')

if conf.out != 'STDOUT':
    fout.close()
