#!/usr/bin/python

#
# overview of cluster resource utilisation
#

import subprocess,collections

cols = collections.defaultdict(float)
hosts= []
cpu_list=[]
mem_list=[]
users= {}

#users to ignore
ignore = ['Debian-e','root','sgeadmin','statd','messageb','sshd',
          'apt-cach','ntp','mysql','dnsmasq','www-data','avahi']
#ignore=[]

#capture output of qhost
out = subprocess.check_output('qhost',shell=True)
qhost = [tok for tok in out.strip().split('\n')[3:]]

#capture head node top 
out = subprocess.check_output("top -b -n 1",shell=True)
top = [x for x in out.strip().split('\n')[7:]]

total_mem = 0.0
total_cpu = 0.0
cpu = 16
host = 'head'
hosts.append(host)

for item in top:
    row  = item.split()
    user = row[1]
    icpu  = float(row[8])/cpu  #%CPU -> %total CPU
    imem  = float(row[9])      #%MEM
    
    total_mem += imem
    total_cpu += icpu
    
    if user in ignore: continue
    
    users[user] = True

    cols[(host,user,'cpu')] += icpu
    cols[(host,user,'mem')] += imem

mem_list.append(total_mem)
cpu_list.append(total_cpu)

#for each host
for line in qhost:
    tok = line.split()
    host      = tok[0].split('.')[0]
    cpu       = int(tok[2])      #number of cores
    
    hosts.append(host)

    #capture output of top run on this host
    out = subprocess.check_output("ssh %s 'top -b -n 1'"%host,shell=True)
    top = [x for x in out.strip().split('\n')[7:]]
    
    total_mem = 0.0
    total_cpu = 0.0
    
    for item in top:
        row  = item.split()
        user = row[1]
        icpu  = float(row[8])/cpu  #%CPU -> %total CPU
        imem  = float(row[9])      #%MEM
        
        total_mem += imem
        total_cpu += icpu
        
        if user in ignore: continue
        
        users[user] = True

        cols[(host,user,'cpu')] += icpu
        cols[(host,user,'mem')] += imem

    mem_list.append(total_mem)
    cpu_list.append(total_cpu)

#heading (host names, cpu%, mem%)
print '%10s'%' ' + ''.join('%16s'%x for x in hosts)
print ' '*10 + ''.join('%10s'%'%cpu'+'%6s'%'%mem' for x in hosts)
print

print '%10s'%'TOTAL',
for i in xrange(len(hosts)):
    print '%9.1f'%cpu_list[i]+'%6.1f'%mem_list[i],
print
print

user_list = users.keys()
user_list.sort()

for user in user_list:
    print '%10s'%user,
    for host in hosts:
        print '%9.1f'%cols[(host,user,'cpu')]+'%6.1f'%cols[(host,user,'mem')],
    print
    print
