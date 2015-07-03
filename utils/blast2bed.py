#!/usr/bin/python

'''
Convert blast hits in tabular format into BED format
suitable for viewing in a genome browser
'''

import argparse,sys

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='blast tabular file created with --outfmt 6 from blast+')
ap.add_argument('--query',action='store_true',help='if set, output annotation against the query sequences, otherwise against the subject')
ap.add_argument('--out',default='STDOUT',type=str,help='output bed file')
conf = ap.parse_args()

#open output file
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

fout.write('track name=%s\n'%conf.inp)

f = open(conf.inp)
for line in f:
    line = line.strip()
    if line == '': continue # ignore blank lines
    if line.startswith('#'): continue #ignore header lines or comments
    tok = line.strip().split('\t')

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
        name = '%s=>%s:%d-%d'%(query,subject,sstart,send)
    else:
        chrom = subject
        start = sstart - 1 #convert to 0-based
        end  = send
        name = '%s=>%s:%d-%d'%(subject,query,qstart,qend)
        
    fout.write('\t'.join([str(x) for x in [chrom,start,end,name]]) + '\n')
    
f.close()


if conf.out != 'STDOUT': fout.close()
