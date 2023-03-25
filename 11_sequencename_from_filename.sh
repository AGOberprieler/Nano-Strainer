#!/bin/bash

# Author: Ulrich Lautenschlager
# Usage: Replace [input_folder] and [output_folder] with your own input and output folder names/paths, preserving the quotation marks.
# In line 15, "substr(id, 1, 12)" will extract the first 12 characters from the filename; it can be modified as needed.

indir="[input_folder]"
outdir="[output_folder]"

for f in "$indir"/*; do

    awk '
    BEGIN {
        id = ARGV[1]  # input file ($f)
        sub(/.*\//, "", id)  # remove directories from path
        id = substr(id, 1, 12)  # extract first 12 characters
    }
    /^>/ {
        print ">"id"_"i+1; i++  # replace sequence name; add sequential no.
    }
    !/^>/ {
        print
    }' $f > "$outdir/${f/*\//}"

done
