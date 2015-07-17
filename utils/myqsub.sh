#!/bin/bash

#
# automatically set the job name to the name of the script file
# plus the directory name (with / converted to -)
# so that qacct information contains a record of the full script name
# to facilitate tracing errors in failed jobs
# this only works if there is no -N option inside the script file
# I implemented this because I cannot find any SGE variable which
# contains the script filename
#

# usage example: myqsub.sh ./scripts/augustus_vesca.sh [script_args...]

set -eu
#set -o pipefail

#make sure a logs directory exists
mkdir -p logs

qsub -V -terse -cwd -e ./logs -o ./logs -N "$(basename $1)_$(dirname $(readlink -f $1) | tr '/' '-')" $@

