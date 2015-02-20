#!/bin/bash
set -e;

#variables
mkdir -p scripts/statistics; #create directory if it doesn't exist

echo "Calculation of statistics started."

#get statistics about mapped and unmapped reads for every contig

	idxstats_bam_file () 
	{
	  for ARG in "$@";
	    do
	      samtools idxstats bam/aln_cleaned_sorted_deduplicated_${ARG}_male.bam >scripts/statistics/stat_${ARG} &
	    done;
	}

	idxstats_bam_file mother father daughter son;

wait;
echo "Idxstats successfully created."




