#!/usr/bin/python

'''
make a gffutils db from a gff file
'''

import sys,argparse
import gffutils

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input gff file')
ap.add_argument('--db',required=True,type=str,help='output db file')
conf = ap.parse_args() #sys.argv

db = gffutils.create_db(conf.inp, force=True, dbfn=conf.db)
