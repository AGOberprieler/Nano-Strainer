#!/bin/bash

# Author: Ulrich Lautenschlager
# Usage: Execute in the folder that contains folders "reads" and "results" from the previous mapping step
# names of all individuals/sequenced samples have to be given - these need to be identical with the filenames in the "reads" folder, without the .fastq suffix

individuals=(individual_1 individual_2 individual3)

####

mkdir -p reads_unmapped

for ind in ${individuals[*]}; do

    # count SAM files
    n_ref=$(ls -1q results/"${ind}"_* | wc -l)
    
    if [ $n_ref -eq 0 ]; then
        echo "ERROR: no SAM file found for individual ${ind}! Individual skipped." > /dev/stderr
        continue
    fi

    # get IDs of reads without any map
    cat results/"${ind}"_* | awk -v "n=$n_ref" '
    
        !/^@/ && and($2, 0x4) {
            unmapped[$1]++
        }
        
        END {
            for (read_id in unmapped) {
                if (unmapped[read_id] == n) {
                    print read_id
                }
            }
        }
        
    ' > ${ind}_unmapped.txt

    # prepare filtering
    sed 's/\(.*\)/@\1 /' ${ind}_unmapped.txt > ${ind}_unmapped_regex.txt
    # more general, but much slower to grep:
    # sed 's/\(.*\)/^@\1[[:space:]]/' ${ind}_unmapped.txt > ${ind}_unmapped_regex.txt

    # extract fastq entries
    grep --no-group-separator -A3 -f ${ind}_unmapped_regex.txt reads/${ind}.fastq > reads_unmapped/${ind}_unmapped.fastq
    
done
