#!/usr/bin/python

'''
make a gffutils db from a gff file 

merge strategies: merge,create_unique,error,warning,replace
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input gff file')
ap.add_argument('--merge',default="create_unique",type=str,help='merge strategy for features with the same ID (merge,create_unique,error,warning)')
ap.add_argument('--db',required=True,type=str,help='output db file')
conf = ap.parse_args() #sys.argv

import gffutils
db = gffutils.create_db(conf.inp, force=True, dbfn=conf.db, merge_strategy=conf.merge)
