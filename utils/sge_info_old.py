#!/usr/bin/python

#
# overview of cluster resource utilisation
#

import subprocess,collections

#convert 2.34G -> 2.34e9, 3.45M -> 3.45e6
def conv(x): return float(x.upper().replace('G','e9').replace('M','e6').replace('K','e3'))

#capture output of qstat -u "*"
#out = subprocess.check_output('qstat -u "*"',shell=True)
#qstat = [tok for tok in out.strip().split('\n')[2:]]

cols = collections.defaultdict(float)
hosts= []
cpu_list=[]
mem_list=[]
users= {}

#users to ignore
ignore = ['Debian-e','root','sgeadmin','statd','messageb','apt-cach','ntp','mysql','dnsmasq','www-data','avahi']
#ignore = []

#capture output of qhost
out = subprocess.check_output('qhost',shell=True)
qhost = [tok for tok in out.strip().split('\n')[3:]]

#head node
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

    cols[host+':'+user+':cpu'] += icpu
    cols[host+':'+user+':mem'] += imem

mem_list.append(total_mem)
cpu_list.append(total_cpu)

#for each host
for line in qhost:
    tok = line.split()
    host      = tok[0].split('.')[0]
    cpu       = int(tok[2])      #number of cores
    #cpu_used  = float(tok[3])    #load average (cores)
    #mem       = conv(tok[4])     #bytes
    #mem_used  = conv(tok[5])
    #swap      = conv(tok[6])
    #swap_used = conv(tok[7])
    #print mem,mem_used,mem_used/mem
    #exit()
    
    hosts.append(host)

    #cpu_list.append(cpu_used/cpu*100)
    #cols[host+':swap'] = swap_used/swap*100
    
    #capture output of top run on the host
    
    out = subprocess.check_output("ssh %s 'top -b -n 1'"%host,shell=True)
    top = [x for x in out.strip().split('\n')[7:]]
    
    total_mem = 0.0
    total_cpu = 0.0
    
    for item in top:
        row  = item.split()
        user = row[1]
        icpu  = float(row[8])/cpu  #%CPU -> %total CPU
        #if not 'g' in row[4] and not 'm' in row[4]: row[4] += 'k'
        #imem  = conv(row[4])   #virtual mem
        imem  = float(row[9])   #%MEM
        
        total_mem += imem
        total_cpu += icpu
        
        if user in ignore: continue
        
        #print row[4],str(imem)
        
        users[user] = True

        cols[host+':'+user+':cpu'] += icpu
        cols[host+':'+user+':mem'] += imem

    mem_list.append(total_mem)
    cpu_list.append(total_cpu)

#heading (host names, cpu%, mem%)
print '%10s'%' ' + ''.join('%16s'%x for x in hosts)
print ' '*10 + ''.join('%10s'%'%cpu'+'%6s'%'%mem' for x in hosts)
print
print '%10s'%'total',
for i in xrange(len(hosts)): print '%9.1f'%cpu_list[i]+'%6.1f'%mem_list[i],
print
print

for user in users:
    print '%10s'%user,
    for host in hosts:
        key = '%s:%s'%(host,user)
        print '%9.1f'%cols[key+':cpu']+'%6.1f'%cols[key+':mem'],
    print
    print
