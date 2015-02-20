#!/bin/bash
#get_average_contig_coverage $perc $contig $bam
#returns "1" -> #low coverage, "0" -> #in tolerance
perc=$1; contig=$2; bam=$3; reference=$4;

pro_coverage=`samtools mpileup -f $reference -r $contig $bam | cut -f 4 | awk '{s+=$1}END{print s}' 2>/dev/null`;

if [ -z "$pro_coverage" ]; then
	coverage=0;
else
	coverage=$pro_coverage;
fi 

pro_positions=`samtools mpileup -f $reference -r $contig $bam | wc -l 2>/dev/null`

if [ -z "$pro_positions" ]; then
	positions=0;
else
	positions=$pro_positions;
fi 

if [[ $positions -eq $zero ]]; then
  echo "1";
  exit;
fi


average_coverage=`echo "scale=3; $coverage/$positions" | bc`;

		if (( $(echo "$average_coverage < $perc"|bc -l) )); then
			echo "1"; #low coverage
		else
			echo "0"; #in tolerance
		fi




