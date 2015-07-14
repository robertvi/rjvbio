#!/usr/bin/python

'''
Read a gffutils database file
output only the longest isoform for each gene
output written to stdout
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--gffdb',required=True,type=str,help='input gffutils database')
ap.add_argument('--type1',default='gene',type=str,help='featuretype denoting a gene')
ap.add_argument('--type2',default='mRNA',type=str,help='featuretype denoting a transcript')
ap.add_argument('--type3',default='exon',type=str,help='featuretype denoting an exon/CDS')
conf = ap.parse_args() #sys.argv

import gffutils

#gffutils database
db = gffutils.FeatureDB(conf.gffdb)

#get all transcript items
for gene in db.features_of_type(conf.type1,order_by=['seqid','start']):
    #get list of relevant sub features (exon or CDS)
    transcript_list = []
    for transcript in db.children(gene,featuretype=conf.type2):
        transcript_length=0
        for exon in db.children(transcript,featuretype=conf.type3,order_by=['seqid','start']):
            transcript_length += exon.end - exon.start + 1
        transcript_list.append([transcript,transcript_length])
        
    if len(transcript_list) == 0: continue
    transcript_list.sort(key=lambda x:x[1],reverse=True)
    transcript = transcript_list[0][0]

    print db[gene]
    print db[transcript]
    for exon in db.children(transcript,featuretype=conf.type3,order_by=['seqid','start']):
        print db[exon]
