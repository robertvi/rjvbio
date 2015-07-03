#!/bin/bash

#
# automatically set the job name to the name of the script file
# plus the directory name
# so that qacct information contains a record of the full script name
# to facilitate tracing errors in failed jobs
# only work if there is no -N option inside the script file
# I cannot find any SGE variable which contains the script filename
#

# usage example: myqsub.sh ./scripts/augustus_vesca.sh

set -eu
set -o pipefail

qsub -V -terse -cwd -e ./logs -o ./logs -N "$(basename $1)_$(dirname $(readlink -f $1) | tr '/' '_')" $1

