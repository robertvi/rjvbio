#!/bin/bash

#
# remove softmasking: ie convert all sequence data to upper case
#

source ~/git_repos/rjvbio/assign_args.sh
ARGS='|--inp|--out|'
assign_args ${ARGS} "$@"

set -eu
set -o pipefail

#convert all sequence to uppercase
cat ${inp} | awk '/^>/ {print $0} /^[^>]/ {print toupper($0)}' > ${out}
