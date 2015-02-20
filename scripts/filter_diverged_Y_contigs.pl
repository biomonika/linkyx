#!/usr/bin/perl -w
use strict;

my @inputs = @ARGV;
chomp(@inputs);

my $reads_number_mother=$inputs[0];
my $reads_number_father=$inputs[1]; 
my $reads_number_daughter=$inputs[2]; 
my $reads_number_son=$inputs[3]; 
my $reference=$inputs[4]; 
my $LINKYX_PATH=$inputs[5]; 

my $mother="statistics/stat_male_mother";
my $father="statistics/stat_male_father";
my $daughter="statistics/stat_male_daughter";
my $son="statistics/stat_male_son";

open(MOTHER_FILE, $mother) or die "Impossible to open source file: $!\n";
open(FATHER_FILE, $father) or die "Impossible to open source file: $!\n";
open(DAUGHTER_FILE, $daughter) or die "Impossible to open source file: $!\n";
open(SON_FILE, $son) or die "Impossible to open source file: $!\n";

my %hash = (); #contig <-> overall coverage hash

my $mother_line; my $mother_cell; my @mother_cells;
my $father_line;  my $father_cell; my @father_cells;
my $daughter_line;  my $daughter_cell; my @daughter_cells;
my $son_line; my $son_cell; my @son_cells;

my $check;
my $comp; my $contig;
my $comp1; my $comp2; my $comp3; my $comp4; my $comp_overall; my $count=0; my $coverage_check; my $girls_check; my $girls;

print "read numbers: mother ".$reads_number_mother." father ".$reads_number_father." daughter ".$reads_number_daughter." son ".$reads_number_son."\n";

#variables to set
my $read_limit=15; #at least 15 reads
my $girl_limit=3; #each of the females can have 3 reads at most


	while (<MOTHER_FILE>) { 


		$mother_line = $_; chomp($mother_line);
		$father_line = <FATHER_FILE>; chomp($father_line);
		$daughter_line = <DAUGHTER_FILE>; chomp($daughter_line);
		$son_line = <SON_FILE>; chomp($son_line);

		#contig length - mapped - unmapped
		#example: comp1260_c0_seq1	273	17	10
		@mother_cells = split("\t", $mother_line); $mother_cell=$mother_cells[2];
		@father_cells = split("\t", $father_line); $father_cell=$father_cells[2];
		@daughter_cells = split("\t", $daughter_line); $daughter_cell=$daughter_cells[2];
		@son_cells = split("\t", $son_line); $son_cell=$son_cells[2];

		#contig being processed
		$contig=$mother_cells[0];

		#check if files are sorted the same way
		$check=(($mother_cells[0] eq $father_cells[0]) eq ($daughter_cells[0] eq $son_cells[0]));
		if ($check==1) {
			#print "ok";
		} else {
			die("Check input files, wrong order of contigs.");
		}
		
		#normalization
		$father_cell = $father_cell/($reads_number_father/$reads_number_mother);
		$daughter_cell = $daughter_cell/($reads_number_daughter/$reads_number_mother);
		$son_cell = $son_cell/($reads_number_son/$reads_number_mother);
		
		#minimal number of reads for males
		$coverage_check=(($father_cell>=$read_limit) && ($son_cell>=$read_limit));

		#for diverged Y-LINKED gene, males should have much greater coverage then females
	
		#mother must have less coverage than father and son
		$comp1=$mother_cell<$father_cell;
		$comp2=$mother_cell<$son_cell;

		#daughter must have less coverage than father and son
		$comp3=$daughter_cell<$father_cell;
		$comp4=$daughter_cell<$son_cell;

		#females have to have low overall number of reads
		$girls_check=(($mother_cell<=$girl_limit) && ($daughter_cell<=$girl_limit));

		#together, daughter and mother must be at most 1/3 of overall coverage
		$girls=$mother_cell+$daughter_cell; my $guys=$father_cell+$son_cell;

		my $fraction;
		if (($girls+$guys)>0) {
			$fraction=$girls/($girls+$guys);		
		} else {
			$fraction=0;
		}

		$comp_overall=$fraction<0.03;
		
		#all conditions satisfied
		$comp=$comp1 && $comp2 && $comp3 && $comp4 && $comp_overall && $coverage_check && $girls_check;	

		#print contig name
		if ($comp==1) {
			#contig fulfilled all the criteria, added to hash
				$count=$count+1;
				$hash{ $contig } = $father_cell + $son_cell; 
				system("echo $contig >>Y_list");
		} 
	}

	my $iso;
	foreach (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
	  print "$_ $hash{$_}\n\n";

			#get sequence for each contig;
			system("perl ${LINKYX_PATH}/scripts/get_sequences_based_on_ids.pl $reference $_");

			$iso = `bash $LINKYX_PATH/scripts/isoform_exists.sh $_ references/male_contig_names`;

			if ($iso eq '') {
 				print "No isoforms.\n";
			} else {
				print "Isoforms exist:\n";
				print $iso;
			}
			print "----------\n";
	}

	print "Number of candidate diverged Y contigs found: ".$count."\n";









