#!/bin/bash

#
# dump top output to file $1 every $2 seconds
# save only username, %CPU, %MEM
#

#usage example: nohup batch_top.sh hostname.top_output 10&

set -eu
set -o pipefail

while [ 1 ]
do
    top -b -n 1 | tail -n +8 | awk '{print $2,$9,$10}' > $1
    sleep $2
done
