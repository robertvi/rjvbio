#!/bin/bash

#
# search for contaminant sequences in a fastq file
# output a csv giving number of blastn hits per adapter file
# per fastq
#

source ~/git_repos/rjvbio/assign_args.sh
ARGS='|--adapters|--fastqs|--tmpdir|--outcsv|--evalue|--reads|'
tmpdir=./tmp
evalue=1e-2
assign_args ${ARGS} "$@"

#create header: names of adapter fasta files
echo -n 'fastq' > ${outcsv}
for x in $(ls -1 ${adapters}) ; do echo -n ,$(basename ${x}) >> ${outcsv} ; done
echo >> ${outcsv}

for x in $(ls -1 ${fastqs})
do
    #print name of fastq
    echo -n $(basename ${x}) >> ${outcsv}
    
    #convert a sample of reads into fasta format in a temp file
    head -n $((reads*4)) ${x}\
        | awk 'NR%4==1{print ">" $0} NR%4==2{print $0}'\
        > ${tmpdir}/tmp_fastq_sample.fa
    
    #use blastn to search for adapter sequences from each file
    for y in $(ls -1 ${adapters})
    do
        blastn -query ${tmpdir}/tmp_fastq_sample.fa\
               -subject ${y}\
               -outfmt 6 \
               -evalue ${evalue}\
        | wc --lines \
        | cut -d ' ' -f 1\
        | xargs echo -n ,\
        >> ${outcsv}
        #print the number of raw hits
    done
    
    #end the line of output
    echo >> ${outcsv}
done

#clean up
rm ${tmpdir}/tmp_fastq_sample.fa
