#!/usr/bin/python

'''
extract only the requested geneids from the gff
create db file using make_gff_database.py first
'''

import sys,argparse,re
from Bio import SeqIO
import gffutils

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--db',required=True,type=str,help='gffutils db file')
ap.add_argument('--ids',default=None,type=str,help='ids of features of interest')
ap.add_argument('--idsfile',default=None,type=str,help='file listing ids of features of interest, one per line')
ap.add_argument('--out',default='STDOUT',type=str,help='output gff')
conf = ap.parse_args() #sys.argv

db = gffutils.FeatureDB(conf.db)

id_list = []

if conf.ids != None:
    id_list += conf.ids
    
if conf.idsfile != None:
    f = open(conf.idsfile)
    for line in f:
        uid = line.strip()
        id_list.append(uid)
    f.close()

if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

def print_features(db,uid,fout,printed):
    if uid in printed: return
    printed[uid] = True
    item = db[uid]
    fout.write(str(item) + '\n')
    
    
    for child in db.children(item,level=1,order_by='start'):
        print_features(db,child,fout,printed)

printed = {}

for uid in id_list: print_features(db,uid,fout,printed)

if conf.out != 'STDOUT': fout.close()
