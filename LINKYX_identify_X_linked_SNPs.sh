#!/bin/bash

#LINKYX_PATH set during toolshed installation

results=$1 #X-linked candidates written here

reference=${2}
mother_bam=${3}
mother_vcf=${4}
father_bam=${5}
father_vcf=${6}
daughter_bam=${7}
daughter_vcf=${8}
son_bam=${9}
son_vcf=${10}

results_fasta=${11}

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
			samtools index ${ARG} &
		fi
	    done;
	}

idxstats_bam_file () 
	{
		#get statistics about mapped and unmapped reads for every contig
		echo "Idxstats" $1 $2;
		stat_name=statistics/stat_${2};
		samtools idxstats $1 >$stat_name #statistics for bam files mapping on male reference

		#with MAPQ>0
		samtools view -bq 1 $1 >statistics/high_score_${2}.bam
		samtools index statistics/high_score_${2}.bam

		stat_high_name=statistics/stat_${2}_high
		samtools idxstats statistics/high_score_${2}.bam >$stat_high_name #statistics for bam files mapping on male reference with MAPQ>0
	}


file_exists $testing $results $reference $mother_bam $mother_vcf $father_bam $father_vcf $daughter_bam $daughter_vcf $son_bam $son_vcf #check if all input files are present
index_bam_file $mother_bam $father_bam $daughter_bam $son_bam #index bam files

#retrieve statistics
mkdir -p statistics;
idxstats_bam_file $mother_bam mother
idxstats_bam_file $father_bam father
idxstats_bam_file $daughter_bam daughter
idxstats_bam_file $son_bam son

newline () {
	echo -e "\n---";
}

#variables

dir="A/";
echo -e "A variant"; newline; mkdir -p A references bam;

grep ">" $reference | sed 's/>//g' >references/reference_contig_names
mother_contigs=`cat references/reference_contig_names | wc -l`;

#DAUGHTER PREPARE FILES

	grep -E '#|AF1=0.' $daughter_vcf | cut -f1 | uniq >${dir}daughter_on_mother_one_allele_heterozygous_ids.vcf
	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/daughter_on_mother_one_allele_heterozygous_ids.vcf"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if exists $a{$res}' $daughter_vcf | sort | uniq | sed '/#/d' >${dir}daughter_on_mother_heterozygous.vcf

	#be aware that if there are more snps, one with allele frequency '0.' is enough for the whole contig to be in contig_pos_ref_heterozygous
	sed '/#/d' ${dir}daughter_on_mother_heterozygous.vcf | perl -pe 's/(comp\S+)\t(\d+)\t\.\t(\S+)\t(\S+).+/$1\t$2\t$3\t$4/' | sort >${dir}contig_pos_ref_alt_daughter_heterozygous
	cut -f 1 ${dir}contig_pos_ref_alt_daughter_heterozygous | sort | uniq >${dir}contigs_where_daughter_heterozygous

	sed '/#/d' $daughter_vcf | perl -pe 's/(comp\S+)\t(\d+)\t\.\t(\S+)\t(\S+).+/$1\t$2\t$3\t$4/' | sort >${dir}contig_pos_ref_alt_daughter

#FATHER PREPARE FILES

	sed '/#/d' $father_vcf | perl -pe 's/(comp\S+)\t(\d+)\t\.\t(\S+)\t(\S+).+/$1\t$2\t$3\t$4/' | sort >${dir}contig_pos_ref_alt_father

#SON PREPARE FILES

	sed '/#/d' $son_vcf | perl -pe 's/(comp\S+)\t(\d+)\t\.\t(\S+)\t(\S+).+/$1\t$2\t$3\t$4/' | sort >${dir}contig_pos_ref_alt_son


echo -e "Step 0 of 3.\nMother must be homozygote.";
	sed '/#/d' $mother_vcf | cut -f1 | sort | uniq >${dir}contigs_where_mother_differs_from_mother;
	differing=`cat ${dir}contigs_where_mother_differs_from_mother | wc -l`;

	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/contigs_where_mother_differs_from_mother"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if !exists $a{$res}' references/reference_contig_names | sort | uniq >${dir}contigs_where_mother_same_as_mother_unchecked

	#check whether mother really does not contain any father variants
	wrong_hits=0;
	cat ${dir}contigs_where_mother_same_as_mother_unchecked | (while read line; do 
		contig_check=`bash ${LINKYX_PATH}/scripts/test.sh ${line} 4 ${dir}contig_pos_ref_alt_father $mother_bam A`
		if [[ $contig_check -gt 0 ]]; then
			echo "contig: " ${line} "contig_check:" $contig_check >>kontrola;
			echo "DID NOT PASS "${line}; 
			wrong_hits=`expr $wrong_hits + 1`
		else
			#echo "0K"; 
			echo ${line} >>${dir}contigs_where_mother_same_as_mother
		fi
	done;
	echo "There is this number of wrong variants: "$wrong_hits;
	)

	same=`cat ${dir}contigs_where_mother_same_as_mother | wc -l`;
	echo -e "Number of contigs, where mother differs: $differing.";
	echo -e "NUMBER of contigs AFTER just this filter: $same/$mother_contigs";

newline;

echo -e "Step 1 of 3.\nSons must be same as mother.";
	#First we have to get rid of the intersection with father, possible Y variants do not matter in the calculation
	comm -12 ${dir}contig_pos_ref_alt_son ${dir}contig_pos_ref_alt_father | sort | uniq >${dir}possible_Yvariants
	cut -f 1 ${dir}possible_Yvariants | sort | uniq >${dir}contigs_with_possible_Yvariants
	
	#if these possible_Yvariants are in the daughter, they are not Y variants
	comm -12 ${dir}contig_pos_ref_alt_daughter ${dir}possible_Yvariants | sort | uniq >${dir}NOT_Yvariants
	cut -f 1 ${dir}NOT_Yvariants  | sort | uniq >${dir}contigs_with_NOT_Yvariants

	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/contigs_with_NOT_Yvariants"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if !exists $a{$res}' ${dir}contigs_with_possible_Yvariants | sort | uniq >${dir}contigs_with_Yvariants

	#contig_pos_ref_alt_Yvariants_confirmed; we deplet just those contigs that contain Y variants and then decide if sons are same as mother

	touch ${dir}contig_pos_ref_alt_Yvariants_confirmed_pro #todo - needs to be created?
cat ${dir}contigs_with_Yvariants | (while read line; do 
	grep ${line} ${dir}contig_pos_ref_alt_father | sort >>${dir}contig_pos_ref_alt_Yvariants_confirmed_pro
done;
)
	cat ${dir}contig_pos_ref_alt_Yvariants_confirmed_pro | sort >${dir}contig_pos_ref_alt_Yvariants_confirmed;
	
	comm -23 ${dir}contig_pos_ref_alt_son ${dir}contig_pos_ref_alt_Yvariants_confirmed | cut -f 1 | sort | uniq >${dir}contigs_where_sons_variation_not_caused_by_father
	#rm ${dir}contig_pos_ref_alt_Yvariants_confirmed_pro;

	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/contigs_where_sons_variation_not_caused_by_father"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if !exists $a{$res}' references/reference_contig_names | sort | uniq >${dir}contigs_where_sons_same_as_mother

	same=`cat ${dir}contigs_where_sons_same_as_mother | wc -l`;
	echo -e "Number of contigs, where sons differ: $differing.";
	echo -e "NUMBER of contigs AFTER just this filter: $same/$mother_contigs";


	newline;
echo -e "Step 2 of 3.\nFather must differ from mother.";
	sed '/#/d' $father_vcf | cut -f1 | sort | uniq >${dir}/contigs_where_father_differs_from_mother;
	differing=`cat ${dir}contigs_where_father_differs_from_mother | wc -l`;
	echo -e "NUMBER of contigs AFTER just this filter: $differing/$mother_contigs";

	newline;

echo -e "Step 3 of 3.\nDaughter must have mother and father alleles.";
	#we need daughter to have alleles from both parents (AF stands for 'allele frequency')

	#unique lines in first file mean that those are not subsets of father variants - otherwise they would be in common
	comm -23 ${dir}contig_pos_ref_alt_daughter_heterozygous ${dir}contig_pos_ref_alt_father | cut -f 1 | sort | uniq >${dir}daughter_NOT_subset_of_father

	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/daughter_NOT_subset_of_father"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if !exists $a{$res}' ${dir}contigs_where_daughter_heterozygous | sort | uniq >${dir}/contigs_where_daughter_from_parents
	

same=`cat ${dir}contigs_where_daughter_from_parents | wc -l`;
	echo -e "NUMBER of contigs AFTER just this filter: $same/$mother_contigs";
	
	newline;

echo -e "Step 4 of 3.\nDaughters and sons can not have any intersection.";
	comm -12 ${dir}contig_pos_ref_alt_daughter ${dir}contig_pos_ref_alt_son >${dir}sons_and_daughters_similiar
	cut -f 1 ${dir}sons_and_daughters_similiar | sort | uniq >${dir}contigs_where_sons_and_daughters_similiar

	perl -MFile::Slurp -ne 'BEGIN{@a{split(/\s/,read_file("A/contigs_where_sons_and_daughters_similiar"))}=()}; $b=$_; @data = split(/\s/, $b); $res=$data[0]; print if !exists $a{$res}' ${dir}contigs_where_daughter_heterozygous | sort | uniq >${dir}contigs_where_daughters_and_sons_different_pro

	#check whether they are really different (small amount of variants is also wrong)

	echo -e "daughter contains son variants?"
	#check whether daughter really does not contain any of son variants
	wrong_hits=0;
	cat ${dir}contigs_where_daughters_and_sons_different_pro | (while read line; do 
		contig_check=`bash ${LINKYX_PATH}/scripts/test.sh ${line} 4 ${dir}contig_pos_ref_alt_daughter $son_bam A`
		if [[ $contig_check -gt 0 ]]; then
			echo "DID NOT PASS "${line}; 
			echo "daughter contains son variants DID NOT PASS "${line} >>kontrola
			wrong_hits=`expr $wrong_hits + 1`
		else
			#echo "0K"; 
			echo ${line} >>${dir}contigs_where_daughters_and_sons_different_un1
		fi
	done;
	echo "There is this number of wrong variants: "$wrong_hits;
	)

	echo -e "son contains daughter variants?"
	#check whether son really does not contain any of daughter variants
	wrong_hits=0;
	cat ${dir}contigs_where_daughters_and_sons_different_pro | (while read line; do 
		contig_check=`bash ${LINKYX_PATH}/scripts/test.sh ${line} 4 ${dir}contig_pos_ref_alt_son $daughter_bam A`
		if [[ $contig_check -gt 0 ]]; then
			echo "DID NOT PASS "${line}; 
			echo "son contains daughter variants DID NOT PASS "${line} >>kontrola
			wrong_hits=`expr $wrong_hits + 1`
		else
			#echo "0K"; 
			echo ${line} >>${dir}contigs_where_daughters_and_sons_different_un2
		fi
	done;
	echo "There is this number of wrong variants: "$wrong_hits;
	)

	comm -12 ${dir}contigs_where_daughters_and_sons_different_un1 ${dir}contigs_where_daughters_and_sons_different_un2 | sort >${dir}contigs_where_daughters_and_sons_different_un

	cat ${dir}contigs_where_daughters_and_sons_different_un | sort | uniq >${dir}contigs_where_daughters_and_sons_different;
	
#rm ${dir}contigs_where_daughters_and_sons_different_un1 ${dir}contigs_where_daughters_and_sons_different_un2 ${dir}contigs_where_daughters_and_sons_different_un;

	same=`cat ${dir}contigs_where_daughters_and_sons_different | wc -l`;
	echo -e "NUMBER of contigs AFTER just this filter: $same/$mother_contigs";

	echo -e "Calculating intersection of candidate contigs from steps 0-3.";

	comm -12 ${dir}contigs_where_mother_same_as_mother ${dir}contigs_where_sons_same_as_mother | sort >${dir}first_couple;
        comm -12 ${dir}contigs_where_father_differs_from_mother ${dir}contigs_where_daughter_from_parents | sort >${dir}second_couple;
	comm -12 ${dir}first_couple ${dir}second_couple | sort >${dir}set1
	comm -12 ${dir}contigs_where_daughters_and_sons_different ${dir}set1 | sort >${dir}set2; 

cat ${dir}set2 >>results; #todo remove
echo -e "Number of contigs after intersection: " `cat ${dir}set2 | wc -l`;

echo -e "Step -1 of 3.\nCoverage must be sufficient.";
	touch ${dir}contigs_with_sufficient_coverage;

cat ${dir}set2 | (while read line; do 

	coverage_check1=`bash ${LINKYX_PATH}/scripts/get_average_contig_coverage.sh 0.75 ${line} $mother_bam $reference 2>/dev/null &`;
	coverage_check2=`bash ${LINKYX_PATH}/scripts/get_average_contig_coverage.sh 0.75 ${line} $father_bam $reference 2>/dev/null &`;
	coverage_check3=`bash ${LINKYX_PATH}/scripts/get_average_contig_coverage.sh 0.75 ${line} $son_bam $reference 2>/dev/null &`;
	coverage_check4=`bash ${LINKYX_PATH}/scripts/get_average_contig_coverage.sh 0.75 ${line} $daughter_bam $reference 2>/dev/null &`;

	#echo "coverage check: " $coverage_check1 $coverage_check2 $coverage_check3 $coverage_check4

	coverage_check=$(($coverage_check1 + $coverage_check2 + $coverage_check3 + $coverage_check4))
		if [[ $coverage_check -gt 0 ]]; then
			echo "LOW COVERAGE "${line}; 
		else
			echo "COVERAGE 0K"; 
			echo ${line} >>${dir}contigs_with_sufficient_coverage
		fi
done;
)

	cat ${dir}contigs_with_sufficient_coverage | sort >X_results.txt; 

#have we found any genes?
	final_contigs=`cat X_results.txt | wc -l`;
	echo -e "$final_contigs contigs found fulfilling criteria.";

	if [ $final_contigs -gt 0 ]; then
		#we found SOME genes

		#retrieve sequences in fasta
		perl ${LINKYX_PATH}/scripts/get_sequences_based_on_ids.pl $reference X_results.txt >X_results_sequences.fasta;
		#sort sequences
		perl ${LINKYX_PATH}/scripts/sorting.pl A/ $reference ${LINKYX_PATH} >>$results

	get_results_in_bam_file () 
	{
	#candidates full_bam_file candidates_bam_file
	cat X_results.txt | xargs samtools view -F 4 -b $1 > $2;
	samtools index $2
	}

	get_results_in_bam_file $mother_bam ${12} #argument for output bam file mother
	get_results_in_bam_file $father_bam ${13} #argument for output bam file father
	get_results_in_bam_file $daughter_bam ${14} #argument for output bam file daughter
	get_results_in_bam_file $son_bam ${15} #argument for output bam file son
	wait;

	#send BAM files with corresponding contigs
	#mutt -s "LINKYX: Results of computation for X-linked genes" $email_address_of_user -a bam/mother.bam -a bam/mother.bam.bai -a bam/father.bam -a bam/father.bam.bai -a bam/daughter.bam -a bam/daughter.bam.bai -a bam/son.bam -a bam/son.bam.bai < messages_to_user/message_results_sentBamX;

	else
		#we found NO genes
		echo "No sequences found." >X_results.txt;
		echo "No sequences found." >X_results_sequences.fasta;
	fi

	newline;
echo -e "Writing output to X_results.txt";
echo -e "Total number of found X-linked genes: $final_contigs";

#TEST OUTPUT
	cat X_results.txt >>$results
	cat X_results_sequences.fasta >>$results_fasta
	ls -lah */* >>$results


