#!/usr/bin/python

'''
convert augustus gff file to GFF3
'''

import argparse,sys,copy
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input augustus GFF file')
ap.add_argument('--prefix',default='',type=str,help='prefix to apply to feature IDs')
ap.add_argument('--suffix',default='',type=str,help='suffix to apply to feature IDs')
ap.add_argument('--out',default='STDOUT',type=str,help='output GFF3 file')
conf = ap.parse_args()

'''
LG1	AUGUSTUS	gene	2176	3500	0.14	+	.	Gene g1
LG1	AUGUSTUS	mRNA	2176	3500	0.14	+	.	mRNA g1.t1 ; Gene g1
LG1	AUGUSTUS	tss	2176	2176	.	+	.	mRNA g1.t1
LG1	AUGUSTUS	exon	2176	2762	.	+	.	mRNA g1.t1
LG1	AUGUSTUS	start_codon	2245	2247	.	+	0	mRNA g1.t1
LG1	AUGUSTUS	intron	2763	3119	0.86	+	.	mRNA g1.t1
LG1	AUGUSTUS	CDS	2245	2762	0.88	+	0	mRNA g1.t1
'''

class augrec:
    def __init__(self,tok):
        self.seqid = tok[0]
        self.source = tok[1]
        self.ftype = tok[2]
        self.start = int(tok[3])
        self.end = int(tok[4])
        self.score = tok[5]
        self.strand = tok[6]
        self.phase = tok[7]
        self.id = tok[8].split()[1]
        
#open output file
if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')

#load all
geneid = None
f = open(conf.inp)
for line in f:
    tok = line.strip().split('\t')
    rec = augrec(tok)
    
    if rec.ftype == 'gene':
        geneid = conf.prefix+rec.id+conf.suffix
        rec.attributes = 'ID=%s'%geneid
    elif rec.ftype == 'mRNA':
        rec.attributes = 'ID=mRNA_%s;Parent=%s'%(geneid,geneid)
    elif rec.ftype in ['CDS','exon']:
        rec.attributes = 'Parent=mRNA_%s'%(geneid)
    else:
        #ignore other feature types
        continue
        
    items = [rec.seqid,rec.source,rec.ftype,rec.start,rec.end,rec.score,rec.strand,rec.phase,rec.attributes]
    fout.write('\t'.join( (str(x) for x in items) ) + '\n' )
    
f.close()

    
if conf.out != 'STDOUT':
    fout.close()
