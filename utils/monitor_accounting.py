#!/usr/bin/python

'''
monitor gridengine accounting file
add new records to sqlite database as they appear
see: man 5 accounting

example usage:
nohup ~/git_repos/rjvbio/utils/monitor_accounting.py --db ./gridengine.db --sleep 10 2> /dev/null > /dev/null &
'''

import argparse

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--db',required=True,type=str,help='sqlite db file to output data into')
ap.add_argument('--init',action='store_true',help='if true init new database (delete old database if present) and exit')
ap.add_argument('--cols',default='/home/vicker/git_repos/rjvbio/colnames.csv',type=str,help='csv file listing column names and datatypes')
ap.add_argument('--file',default='/var/lib/gridengine/blacklace/common/accounting',type=str,help='accounting file to read input from')
ap.add_argument('--sleep',default=0,type=int,help='poll file every x second, 0 => run once then exit')
conf = ap.parse_args() #sys.argv

import sqlite3,os,time

columns = 45

#read in the gridengine column names and types
f = open(conf.cols)
col_list = [line.strip().split(',') for line in f]
f.close()

assert len(col_list) == columns

#initialise a new database then quit
if conf.init == True:
    #remove any existing database file
    try:
        os.remove(conf.db)
    except OSError:
        pass
        
    #create and populate database
    db = sqlite3.connect(conf.db)
    cur = db.cursor()
    cmd = 'create table data (fposn integer primary key,' + ','.join(['[%s] %s'%(col[0],col[1]) for col in col_list])  + ')'
    cur.execute(cmd)
    db.commit()
    db.close()
        
    assert os.path.isfile(conf.file)
    exit()


prev_mtime = None

while True:
    mtime = os.stat(conf.file).st_mtime
    
    if mtime == prev_mtime:
        #file didn't change
        print 'sleeping'
        time.sleep(conf.sleep)
        continue
        
    prev_mtime = mtime
    db = sqlite3.connect(conf.db)
    cur = db.cursor()

    #find where the scan got up to last time it ran
    cur.execute('select max(fposn) from data')
    maxfposn = cur.fetchone()[0]

    f = open(conf.file)
    #skip to last known record + 1
    if maxfposn != None:
        f.seek(maxfposn)
        line = f.readline()

    ct = 0

    while True:
        fposn = f.tell()
        line = f.readline()
        if line == '': break #end of file
        if line.startswith('#'): continue
        tok = line.strip().split(':')
        if len(tok) != columns: continue
        
        items = [str(fposn)]
        for j in xrange(columns):
            if tok[j] == 'NONE':
                items.append('null')
            elif col_list[1] == 'integer':
                items.append('%d'%tok[j])
            elif col_list[1] == 'real':
                items.append('%f'%tok[j])
            else:
                items.append("'%s'"%tok[j])
                
        cmd = 'insert into data values ('\
            + ','.join(items)\
            + ')'
        cur.execute(cmd)
        ct += 1
        
        if ct % 1000000 == 0:
            db.commit()
            print ct
        
    f.close()
    db.commit()
    db.close()

    print 'added',ct,'records'
    
    #do not loop if we're only doing one scan
    if conf.sleep == 0: break
