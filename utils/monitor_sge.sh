#!/bin/bash

#
# run using eg:
# watch --interval 5 monitor_sge.sh
#

echo =====QHOST 
qhost
echo

echo =====QSTAT -f
qstat -f -u "*"
echo

echo =====QSTAT
qstat -u "*"
echo

echo =====QSTAT -j
for uid in $(qstat | grep '^[0-9]' | cut -f 1 -d ' ' | uniq)
do
    qstat -j ${uid}\
        | grep -e '^usage' -e '^hard' -e '^job_number' -e '^job_name'
    echo
done
