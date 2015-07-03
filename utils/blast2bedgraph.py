#!/usr/bin/python

'''
Convert blast hits in tabular format into BED format
suitable for viewing in a genome browser
'''

import argparse,sys

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',default='STDIN',type=str,help='blast tabular file created with --outfmt 6 from blast+')
ap.add_argument('--query',action='store_true',help='if set, output annotation against the query sequences, otherwise against the subject')
ap.add_argument('--scorecol',default=10,type=int,help='which column offset to use as the score column')
ap.add_argument('--out',default='STDOUT',type=str,help='output bed file')
conf = ap.parse_args()

#open files
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')
if conf.inp == 'STDIN': f = sys.stdin
else:                   f = open(conf.inp)

#fout.write('track type=bedGraph\n')
for line in f:
    line = line.strip()
    if line == '': continue # ignore blank lines
    if line.startswith('#'): continue #ignore header lines or comments
    tok = line.strip().split()

    query = tok[0]
    subject = tok[1]
    qstart = min(int(tok[6]),int(tok[7]))
    qend = max(int(tok[6]),int(tok[7]))
    sstart = min(int(tok[8]),int(tok[9]))
    send = max(int(tok[8]),int(tok[9]))

    if conf.query:
        chrom = query
        start = qstart - 1 #convert to 0-based
        end  = qend
    else:
        chrom = subject
        start = sstart - 1 #convert to 0-based
        end  = send

    score = tok[conf.scorecol].strip()
    fout.write('\t'.join([str(x) for x in [chrom,start,end,score]]) + '\n')
    
f.close()


if conf.out != 'STDOUT': fout.close()
if conf.inp != 'STDIN': f.close()
