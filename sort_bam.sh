#!/bin/bash
results=$1;
bam_file=$2;

hash samtools 2>/dev/null || { echo "I require samtools in path but it's not installed.  Aborting."; exit 1; }
echo "LINKYX_PATH:" $LINKYX_PATH
echo "SAMTOOLS:" $SAMTOOLS
echo "PATH:" $PATH
samtools sort $bam_file $results 2>log.txt

