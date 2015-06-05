#!/usr/bin/python

'''
convert cuffmerge GTF into augustus exon hints file
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='cuffmerge GTF file')
ap.add_argument('--hint_type',default='exon',type=str,help='hint type, eg exon, exonpart')
ap.add_argument('--pri',default=4,type=int,help='hint priority')
ap.add_argument('--src',default='W',type=str,help='type of hint (W=RNASeq)')
ap.add_argument('--out',default='STDOUT',type=str,help='output Augustus hints file')
conf = ap.parse_args() #sys.argv

import sys

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

f = open(conf.inp)
for line in f:
    if not '\texon\t' in line: continue
    tok = line.strip().split('\t')
    tok2 = tok[8].split(';')[1].strip().split()
    assert tok2[0] == 'transcript_id'
    tran_id = tok2[1].replace('"','')
    tok[8] = 'grp=%s;pri=%d;src=%s'%(tran_id,conf.pri,conf.src)
    tok[2] = conf.hint_type
    fout.write('\t'.join(tok) + '\n')
    
f.close()

if conf.out != 'STDOUT': fout.close()
