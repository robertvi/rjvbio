#!/bin/bash

#
# find scaffold positions in pseudomolecules using blastn
# blastn must be in the path

source ~/git_repos/rjvbio/assign_args.sh
ARGS='|--scaffs|--pseudo|--out|'
assign_args ${ARGS} "$@"

#allow only 100% ident, no mismatches or gaps, qstart==1
#penalise gaps and mismatches
#retain best hit per scaffold
blastn -query ${scaffs} -subject ${pseudo} -outfmt 6\
       -ungapped -dust no -max_target_seqs 1\
    | awk '$3==100&&$5==0&&$6==0&&$7==1{print} '\
    > ${out}

