#!/bin/bash

#
# show info on the latest jobs to have finished
# should be used in conjunction with monitor_accounting.py
# which turns the accounting flat file into an sqlite database
#

set -eu
set -o pipefail

source /home/vicker/git_repos/rjvbio/assign_args.sh
ARGS='|--db|--owner|--limit|'
limit=30
db=/home/vicker/gridengine/gridengine.db
assign_args ${ARGS} "$@"

tmpfile=~/tmp/${RANDOM}_${RANDOM}.tmp
mkdir -p ~/tmp

cat > ${tmpfile} <<XXX
.headers on
.mode column
.width 11 40 8 4 16 4 3
XXX

echo -n 'select hostname,job_name,job_number,task_number,maxvmem,exit_status,failed' >> ${tmpfile}
echo " from data where owner == '${owner}' order by fposn desc limit ${limit};" >> ${tmpfile}

cat ${tmpfile} | sqlite3 ${db}

rm ${tmpfile}
