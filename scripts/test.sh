#!/bin/bash
#example usage: ./test.sh comp100628_c0_seq1 5 A/contig_pos_ref_alt_father mother.bam

contig_to_check=$1
percentage=$2
forbidden=$3
bam=$4
dir=$5;

hits=0;

grep $contig_to_check $forbidden | (while read line; do 

	#tresholds described in table 1C
		result=`bash ${LINKYX_PATH}/scripts/contains_these_variants.sh $percentage $bam ${line}`
	
	echo "result: " $result  >>kontrola;
	hits=$(($hits + $result))
	echo "hits: " $hits  >>kontrola;
done;

echo $hits;
)

