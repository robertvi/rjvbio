#!/bin/bash

#
# automatically set the job name to the name of the script file
# so that qacct information contains a record of the script name
# to facilitate tracing errors in failed jobs
#

# usage example: myqsub.sh ./scripts/augustus_vesca.sh
# becomes: qsub -N augustus_vesca.sh ./scripts/augustus_vesca.sh
# this only work if there is no -N option inside the script file
# I cannot find any SGE variable which contains the script filename

set -eu

qsub -terse -cwd -V -N "$(basename $1)[$(readlink -f $1 | tr '/' '_')]" $1

