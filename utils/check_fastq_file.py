#!/usr/bin/python

'''
check fastq file(s) for obvious signs of corruption,
written to deal with files where end of line
characters may be missing
'''

import sys,argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inps',required=True,nargs='+',type=str,help='fastq filename(s)')
ap.add_argument('--newline',default='unix',type=str,help='newline: unix, dos or mac')
conf = ap.parse_args() #sys.argv

buffsize = 10000
newline = {'unix':'\n','dos':'\r\n','mac':'\r'}[conf.newline]
fout = sys.stdout

for fname in conf.inps:
    readct = 0
    f = open(fname)

    while True:
        readct += 1
        pos = f.tell()
        
        #read in a chunk of data ignoring line endings
        buff = f.read(buffsize)
        if buff == '':
            fout.write('%s %d reads ok\n'%(fname,readct))
            fout.flush()
            break #end of file
            
        #split into lines
        read = buff.split(newline)
        
        if len(read) < 4:
            fout.write('%s read %d less than four lines available %s\n'%(fname,readct,read[0]))
            fout.flush()
            break
            
        read = read[:4]
        readsize = len(newline.join(read) + newline)
        f.seek(pos+readsize)
        
        if not read[0].startswith('@'):
            fout.write('%s read %d header invalid: %s\n'%(fname,readct,read[0]))
            fout.flush()
            break
            
        if not read[2].startswith('+'):
            fout.write('%s read %d second header invalid: %s\n'%(fname,readct,read[2]))
            fout.flush()
            break
            
        if len(read[1]) != len(read[3]):
            fout.write('%s read %d seq and qual lengths differ: %d %d\n'%(fname,readct,len(read[1]),len(read[3])))
            fout.flush()
            break
            
        errflag = False
        for x in read[1]:
            if x not in 'ATCGNatcgn':
                fout.write('%s read %d invalid seq character: %s\n'%(fname,readct,x))
                fout.flush()
                errflag = True
                break
                
        if errflag: break
                
        for x in read[3]:
            if not (33 <= ord(x) <= 74):
                fout.write('%s read %d invalid qual character: %s\n'%(fname,readct,x))
                fout.flush()
                errflag = True
                break
        
        if errflag: break
        
    f.close()
