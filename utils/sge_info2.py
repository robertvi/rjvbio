#!/usr/bin/python

#
# overview of cluster resource utilisation
# relies on batch_top script to update info
# should be more efficient than running ssh 'top' for each update!

import subprocess,collections,os,time,calendar

batch_top_dir='/home/vicker/gridengine'

cols = collections.defaultdict(float)
hosts= []
cpu_list=[]
mem_list=[]
users= {}

#users to ignore
ignore = ['Debian-e','root','sgeadmin','statd','messageb','sshd',
          'apt-cach','ntp','mysql','dnsmasq','www-data','avahi']
#ignore=[]

#create fake entry for head node, giving name, ncpus, mem (not used yet)
qhost = ['blacklace00.blacklace - 16 - 15.7G']

#capture output of qhost
out = subprocess.check_output('qhost',shell=True)
qhost += [tok for tok in out.strip().split('\n')[3:]]

#for each host
for line in qhost:
    tok = line.split()
    host      = tok[0].split('.')[0]
    cpu       = int(tok[2])      #number of cores
    
    hosts.append(host)

    #read output of top created by the batch_top script
    #only works if batch_top is running on all the nodes!
    fname = batch_top_dir + '/' + host + '.top'
    if not os.path.isfile(fname): continue

    mtime = os.path.getmtime(fname)
    now = calendar.timegm(time.gmtime())
    if now - mtime > 100: continue
    
    f = open(fname)
    top = [x.strip() for x in f]
    f.close()
    
    total_mem = 0.0
    total_cpu = 0.0
    
    for item in top:
        row  = item.split()
        user = row[0]
        icpu  = float(row[1])/cpu  #%CPU -> %total CPU
        imem  = float(row[2])      #%MEM
        
        total_mem += imem
        total_cpu += icpu
        
        if user in ignore: continue
        
        users[user] = True

        cols[(host,user,'cpu')] += icpu
        cols[(host,user,'mem')] += imem

    mem_list.append(total_mem)
    cpu_list.append(total_cpu)

def fcpu(x):
    if x < 0.5: return ' '*8+'.'
    return '%9.0f'%x
    
def fmem(x):
    if x < 0.5: return ' '*5+'.'
    return '%6.0f'%x

#heading (host names, cpu%, mem%)
print '%10s'%' ' + ''.join('%16s'%x for x in hosts)
print ' '*10 + ''.join('%10s'%'%cpu'+'%6s'%'%mem' for x in hosts)
print

print '%10s'%'TOTAL',
for i in xrange(len(hosts)):
    print fcpu(cpu_list[i])+fmem(mem_list[i]),
print
print

user_list = users.keys()
user_list.sort()

for user in user_list:
    print '%10s'%user,
    for host in hosts:
        print fcpu(cols[(host,user,'cpu')])+fmem(cols[(host,user,'mem')]),
    print
    print
