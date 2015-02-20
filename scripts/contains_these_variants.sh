#!/bin/bash
#input example: perc% bam comp100134_c0_seq1	28	T	C
#returns true if the alternative allele occurs in >= than perc (%) cases

perc=$1; bam=$2; contig=$3; position=$4; ref=$5; alt=$6; 

echo contig $contig position $position ref $ref alt $alt >>kontrola;

depth=0; 

perl ${LINKYX_PATH}/scripts/bam_analysis.pl $bam reference.fasta $contig $position; #writes to mpileup;
pro_hits=`grep -w "$alt" mpileup | cut -f2 -d":"`;

if [ -z "$pro_hits" ]; then
	hits=0;
else
	hits=$pro_hits;
fi 

depth=`grep -w "DEPTH" mpileup | cut -f2 -d":"`;

echo "hits: " $hits " depth: " $depth >>kontrola;

if [[ $depth -eq $zero ]]; then
  echo "0";
  exit;
fi

if [ $depth -lt 5 ]; then
	#depth <=4 
	if [ $hits -gt 0 ]; then
		echo "1"; #in this depth, none of variants can be tolerated
	else
		echo "0"; #ok
	fi
else
	if [ $depth -lt 25 ]; then
		#depth 5..24
		if [ $hits -gt 1 ]; then
			echo "1"; #in this depth, one variant can be tolerated
		else
			echo "0"; #ok
		fi
	else
		#depth >=25
		hits=$[hits*100];
		percentage=`echo "scale=3; $hits/$depth" | bc`;

		if (( $(echo "$perc < $percentage"|bc -l) )); then
			echo "1"; #too much of WRONG variants
		else
			echo "0"; #in tolerance
		fi
	fi
fi


