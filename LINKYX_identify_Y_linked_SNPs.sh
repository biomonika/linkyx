#!/bin/bash

#LINKYX_PATH set during toolshed installation

results=$1 #X-linked candidates written here

reference=${2}
mother_male_bam=${3}
father_male_bam=${4}
daughter_male_bam=${5}
son_male_bam=${6}

results_fasta=${7}

mkdir -p Y references; dir="Y/";
grep ">" $reference | sed 's/>//g' >references/male_contig_names

file_exists () 
	{
	for ARG in "$@";
	do
		if [ ! -f ${ARG} ]; then
		echo "File $ARG not found"; 
		fi
	done;
	}

#index bam files
	
index_bam_file () 
	{
	for ARG in "$@";
		do
		bam_file_index=${ARG}".bai"
		if [ -e $bam_file_index ];
		then	
			echo "Index already exists. Not being created again.";
		else
			samtools index ${ARG};
		fi
	    done;
	}

#get statistics about mapped and unmapped reads for every contig

idxstats_bam_file () 
	{
		#get statistics about mapped and unmapped reads for every contig
		echo "Idxstats" $1 $2;
		stat_male_name=statistics/stat_male_${2};
		samtools idxstats $1 >$stat_male_name #statistics for bam files mapping on male reference

		#with MAPQ>0
		samtools view -bq 1 $1 >statistics/high_score_male_${2}.bam
		samtools index statistics/high_score_male_${2}.bam

		stat_male_high_name=statistics/stat_male_${2}_high
		samtools idxstats statistics/high_score_male_${2}.bam >$stat_male_high_name #statistics for bam files mapping on male reference with MAPQ>0
	}

echo "Checking existence of input bam files."
file_exists $mother_male_bam $father_male_bam $son_male_bam $daughter_male_bam
index_bam_file $mother_male_bam $father_male_bam $son_male_bam $daughter_male_bam

#retrieve statistics
mkdir -p statistics;
idxstats_bam_file $mother_male_bam mother
idxstats_bam_file $father_male_bam father
idxstats_bam_file $daughter_male_bam daughter
idxstats_bam_file $son_male_bam son

#read statistics
reads_number_mother=`awk '{read_number+=($3+$4)} END {print read_number}' statistics/stat_male_mother`;
reads_number_father=`awk '{read_number+=($3+$4)} END {print read_number}' statistics/stat_male_father`;
reads_number_daughter=`awk '{read_number+=($3+$4)} END {print read_number}' statistics/stat_male_daughter`;
reads_number_son=`awk '{read_number+=($3+$4)} END {print read_number}' statistics/stat_male_son`;

echo "Reads numbers: " $reads_number_mother $reads_number_father $reads_number_daughter $reads_number_son

#filter contigs, which are candidate Y genes - very divergent from their X counterparts
perl ${LINKYX_PATH}/scripts/filter_diverged_Y_contigs.pl $reads_number_mother $reads_number_father $reads_number_daughter $reads_number_son $reference ${LINKYX_PATH} >Y_results.txt;

echo "Y_results: " >>$results

final_contigs=`cat Y_list | wc -l`
#have we found any genes?

if [ $final_contigs -gt 0 ]; then
	#we found SOME genes
	#retrieve sequences in fasta
	perl ${LINKYX_PATH}/scripts/get_sequences_based_on_ids.pl $reference Y_list >Y_results_sequences.fasta;

	get_results_in_bam_file () 
	{
	#candidates full_bam_file candidates_bam_file
	cat Y_list | xargs samtools view -F 4 -b $1 > $2;
	samtools index $2
	}

	get_results_in_bam_file $mother_male_bam ${8} #argument for output bam file mother
	get_results_in_bam_file $father_male_bam ${9} #argument for output bam file father
	get_results_in_bam_file $daughter_male_bam ${10} #argument for output bam file daughter
	get_results_in_bam_file $son_male_bam ${11} #argument for output bam file son
	wait;
else
		#we found NO genes
		echo "No sequences found." >Y_results.txt;
		echo "No sequences found." >Y_results_sequences.fasta;
	fi

#output results to Galaxy
	cat Y_results.txt >>$results
	cat Y_results_sequences.fasta >>$results_fasta
	ls -lah */* >>$results
