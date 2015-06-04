#!/usr/bin/python

'''
convert a multifasta alignment to nexus format
ready for mybayes
'''

'''
example NEXUS file

#NEXUS
Begin data;
Dimensions ntax=4 nchar=15;
Format datatype=dna missing=? gap=-;
Matrix
Species1   atgctagctagctcg
Species2   atgcta??tag-tag
Species3   atgttagctag-tgg
Species4   atgttagctag-tag           
;
End;

'''

import argparse
ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input multi FASTA alignment file')
ap.add_argument('--out',default='STDOUT',type=str,help='output NEX file')
ap.add_argument('--datatype',default='dna',type=str,help='datatype')
ap.add_argument('--missing',default='N',type=str,help='symbol denoting missing data')
ap.add_argument('--gap',default='-',type=str,help='symbol denoting a gap in the alignment')
conf = ap.parse_args()

import sys
from Bio import SeqIO
ntax = 0
nchar = None

#count taxa and sequence length
for rec in SeqIO.parse(conf.inp, "fasta"):
    ntax += 1
    if nchar == None:
        nchar = len(rec.seq)
    else:
        assert nchar == len(rec.seq)

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

#output nexus to stdout
fout.write('#NEXUS\n')
fout.write('Begin data;\n')
fout.write('Dimensions ntax=%d nchar=%d;\n'%(ntax,nchar))
fout.write('Format datatype=dna missing=N gap=-;\n')
fout.write('Matrix\n')

for rec in SeqIO.parse(conf.inp, "fasta"):
    fout.write(rec.id + '   ' + str(rec.seq) + '\n')
    
fout.write(';\n')
fout.write('End;\n')

if conf.out != 'STDOUT': fout.close()
