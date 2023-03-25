#!/bin/bash

# Author: Ulrich Lautenschlager
# Usage: Replace "/path/to/flagstat_results" with the directory containing the results files from samtools flagstat and run: bash 9_extract_mapping_statistics.sh

for f in /path/to/flagstat_results/*; do
        printf "${f/*\//}\t" # filename without directory
        awk '/primary mapped/ {p=$1} /secondary/ {s=$1} /supplementary/ {u=$1} END {print p"\t"s"\t"u}' $f
done > mapping_statistics.csv
