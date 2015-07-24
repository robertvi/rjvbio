#!/bin/bash

#
# kill batch_top on all nodes (inc head)
#

set -eu

kill $(ps -ef | grep batch_top | grep -v grep | awk '{print $2}')


#for each worker node
for x in $(seq 1 5)
do
    hname=$(printf "blacklace%02d" ${x})
    ssh ${hname} 'kill $(ps -ef | grep batch_top | grep -v grep | awk "{print \$2}")'
done

