#!/bin/bash

# Author: Ulrich Lautenschlager
# Usage: This script is executed once WITHIN a directory containing one subfolder per amplified marker, each containing one subfolder per individual, each containing the output of the respective CANU assembly run.
# Please ensure that the subdirectories are given meaningful names because the output files will also be named accordingly.

mkdir -p summary
rm -f summary/*

for dir in */*/; do

    prefix="$(basename $dir)"
    f_out="summary/$(dirname $dir).txt"
    
    # input file names
    f_report="$dir/${prefix}.report"
    f_contigs="$dir/${prefix}.contigs.fasta"
    f_unassembled="$dir/${prefix}.unassembled.fasta"
    f_sserr="$dir/${prefix}.seqStore.err"

    # get assembly info
    n_reads=""
    n_assembled=""
    if test -f "$f_report"; then
        n_reads=$(awk '/^\[CORRECTION\/READS\]/ || p==1 {print; p=1} /^$/ {p=0}' "$f_report" | awk '/Found [0-9]* reads./ {print $3}')
        n_assembled=$(awk '/^\[UNITIGGING\/READS\]/ || p==1 {print; p=1} /^$/ {p=0}' "$f_report" | awk '/Found [0-9]* reads./ {print $3}')
    fi
    
    n_contigs=""
    n_reads_per_contig=""
    if test -f "$f_contigs"; then
        n_contigs=$(grep -c "^>" "$f_contigs")
        n_reads_per_contig=$(grep "^>" "$f_contigs" | sed 's/.*reads=\([0-9]*\).*/\1/' | tr "\n" "\t" | sed 's/\t$//')
    fi
    
    n_unassembled=""
    if test -f "$f_unassembled"; then
        n_unassembled=$(grep -c "^>" "$f_unassembled")
    fi

    n_short=""
    p_short=""
    if test -f "$f_sserr"; then
        # Note: Only the first "Short"-line of the file is processed!
        n_short=$(awk '$1=="Short" {print $2; exit}' "$f_sserr")
        p_short=$(awk '$1=="Short" {print $3; exit}' "$f_sserr" | tr -d "%")
    fi
    
    # add placeholders
    if test -z "$n_reads"; then
        n_reads="-"
    fi
    if test -z "$n_assembled"; then
        n_assembled="-"
    fi
    if test -z "$n_contigs"; then
        n_contigs="-"
    fi
    if test -z "$n_reads_per_contig"; then
        n_reads_per_contig="-"
    fi
    if test -z "$n_unassembled"; then
        n_unassembled="-"
    fi
    if test -z "$n_short"; then
        n_short="-"
    fi
    if test -z "$p_short"; then
        p_short="-"
    fi

    echo -e "$prefix\t$n_short\t$p_short\t$n_reads\t$n_assembled\t$n_unassembled\t$n_contigs\t$n_reads_per_contig" >> "$f_out"

done
