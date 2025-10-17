#!/bin/bash

bc_pattern=$1
fastq1=$2
fastq2=$3
fastq1_output=$4
fastq2_output=$5
white_list=$6

umi_tools extract --bc-pattern=$bc_pattern \
    --stdin $fastq1 \
    --stdout $fastq1_output \
    --read2-in=$fastq2 \
    --read2-out=$fastq2_output \
    --whitelist=$white_list 
