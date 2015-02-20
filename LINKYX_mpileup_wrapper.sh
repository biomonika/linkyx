#!/bin/bash
results=$1;
bam_file=$2;
reference=$3;

sorted_bam=${bam_file}_sorted

samtools mpileup -uf $reference $bam_file 2>log.txt | bcftools view -p 0.85 -cgv - 2>log.txt >$results
