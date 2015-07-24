#!/bin/bash

#
# launch batch_top on all nodes (inc head)
#

set -eu

delay=10
cmd='/home/vicker/git_repos/rjvbio/utils/batch_top.sh'
out='/home/vicker/gridengine/blacklace00.top'

#for head node
nohup ${cmd} ${out} ${delay} </dev/null >/dev/null 2>/dev/null &

#for each worker node
for x in $(seq 1 5)
do
    hname=$(printf "blacklace%02d" ${x})
    out="/home/vicker/gridengine/${hname}.top"
    ssh ${hname} "nohup ${cmd} ${out} ${delay} </dev/null >/dev/null 2>/dev/null &"
done

