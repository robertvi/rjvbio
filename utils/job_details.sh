#!/bin/bash

#
# show detailed info on a job selected by job and task number
# should be used in conjunction with monitor_accounting.py
# which turns the accounting flat file into an sqlite database
#

set -eu
set -o pipefail

source ~/git_repos/rjvbio/assign_args.sh
ARGS='|--db|--jobid|--taskid|'
taskid=0
db=~/gridengine/gridengine.db
assign_args ${ARGS} "$@"

tmpfile=~/tmp/${RANDOM}_${RANDOM}.tmp

cat > ${tmpfile} <<XXX
.headers on
.mode column
XXX

if [ ${taskid} == 0 ]; then
    echo "select hostname,job_name,job_number,task_number,maxvmem,exit_status,failed,category from data where job_number == ${jobid};" >> ${tmpfile}
else
    echo "select hostname,job_name,job_number,task_number,maxvmem,exit_status,failed,category from data where job_number == ${jobid} and task_number == ${taskid};" >> ${tmpfile}
fi

cat ${tmpfile} | sqlite3 ${db}

rm ${tmpfile}
