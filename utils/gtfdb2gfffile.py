#!/usr/bin/python

'''
convert GTF database (made by gffutils) into a gff3 file

rationale: gffutils seems good at loading gtf files, but
I've not so far found a way to easily convert the output
into a gff3 file, because gffutils works by storing the attributes as
raw json in the database, therefore when converting back to a flat file
it becomes gtf again, even if you change the dialect to be gff3-like

rather than write a direct gtf -> gff3 convert I'll use the gffutils loader
to deal with variations in gtf format, then convert from the gffutils
database into gff3 flat file format
'''

import argparse,sys,copy

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpdb',required=True,type=str,help='input gffutils GTF database')
ap.add_argument('--conv',default="{'exon':[('Parent','transcript_id')],'transcript':[(2,'mRNA'),('ID','transcript_id'),('Parent','gene_id')],'gene':[('ID','gene_id')]}",type=str,help='method of converting (see source code for details)')
conf = ap.parse_args()

import sqlite3,simplejson

db = sqlite3.connect(conf.inpdb)

c = db.cursor()

#how to convert attributes
#convstr=

conv=eval(conf.conv)

for row in c.execute('SELECT * FROM features'):
    tok = [str(x) for x in row[1:9]]
    attr = simplejson.loads(row[9])
    
    ftype = tok[2]
    
    if not ftype in conv: continue
    
    newattr = {}
    for x in conv[ftype]:
        if type(x[0]) == int:
            tok[x[0]] = x[1]
        else:
            newattr[x[0]] = attr[x[1]][0]
        
    newattrstr = ';'.join([ '%s=%s'%(k,v) for k,v in newattr.iteritems() ])
    print '\t'.join(tok) + '\t' + newattrstr
