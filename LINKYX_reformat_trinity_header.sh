#!/bin/bash

results=$1 #reformatted fasta sequence written here
trinity_fasta=${2}

perl -pe 's/(>[\S]+).+/$1/' $trinity_fasta >$results
samtools faidx $results


