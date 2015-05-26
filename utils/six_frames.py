#!/usr/bin/python

'''
translate nucleotide sequences into six frames protein sequences 
--help for help info
'''

import sys,argparse
import rjvbio.seq
import Bio.SeqIO,Bio.SeqRecord

if len(sys.argv) < 2: sys.argv.append('--help')
ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',nargs='+',required=True,type=str,help='input sequence file(s)')
ap.add_argument('--inpformat',default='fasta',type=str,help='format of input file(s), eg fasta,fastq,genbank see http://biopython.org/wiki/SeqIO#File_Formats for details')
ap.add_argument('--out',default='STDOUT',type=str,help='output FASTA file')
#ap.add_argument('--suffix',default=r'_frame%%+d',type=str,help='suffix to append to original gene id, %%d will be replaced by the frame number')
conf = ap.parse_args()

suffix1 = '_frame%+d'
suffix2 = '_frame-%d'

if conf.out == 'STDOUT':
    fout = sys.stdout
else:
    fout = open(conf.out,'wb')
    
for fname in conf.inp:
    for rec in Bio.SeqIO.parse(fname,conf.inpformat):
        frames = rjvbio.seq.six_frames(rec.seq)
        
        for i in [0,1,2]:
            newrec = Bio.SeqRecord.SeqRecord(frames['%+d'%i], id=rec.id + suffix1%i,description='')
            Bio.SeqIO.write(newrec,fout,"fasta")
            
        for i in [0,1,2]:
            newrec = Bio.SeqRecord.SeqRecord(frames['-%d'%i], id=rec.id + suffix2%i,description='')
            Bio.SeqIO.write(newrec,fout,"fasta")
            
if conf.out != 'STDOUT': fout.close()
