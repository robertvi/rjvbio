#!/usr/bin/python

'''
list parental ids of a list of feature ids
eg from a set of mRNA ids, output the parental gene ids
remove duplicated parental ids
if no parent return the original id itself
'''

import sys,argparse,re
from Bio import SeqIO
import gffutils

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--db',required=True,type=str,help='gffutils db file')
ap.add_argument('--ids',default=None,type=str,help='ids of features of interest')
ap.add_argument('--idsfile',default=None,type=str,help='file listing ids of features of interest, one per line')
ap.add_argument('--featuretypes',default=None,nargs='*',help='list of feature types to output')
ap.add_argument('--out',default='STDOUT',type=str,help='output file containing ids')
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

parental_dict = {}

for uid in id_list:
    for parent in db.parents(uid,featuretype=conf.featuretypes):
        if not parent in parental_dict:
            fout.write(str(parent.id) + '\n')
        parental_dict[parent] = True

if conf.out != 'STDOUT': fout.close()
