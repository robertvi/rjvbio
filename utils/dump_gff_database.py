#!/usr/bin/python

'''
dump a gffutils database to gff3 format

not working yet!

note also gffutils-cli provides an command line interface to gffutils
in my case its installed at ~/.local/bin/gffutils-cli
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inpdb',required=True,type=str,help='gffutils database file')
ap.add_argument('--out',default='STDOUT',type=str,help='output gff3 file')
conf = ap.parse_args() #sys.argv

import gffutils

#gffutils database
db = gffutils.FeatureDB(conf.inpdb)

#open output file
if conf.out == 'STDOUT': fout = sys.stdout
else:                    fout = open(conf.out,'wb')

for x in db: print db[x]

