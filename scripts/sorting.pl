#!/usr/bin/perl -w
use strict;
#usage: perl sorting.pl A/

my @inputs = @ARGV;
chomp(@inputs);

my $dir = $inputs[0];
my $reference = $inputs[1];
my $LINKYX_PATH = $inputs[2];
my $results = "X_results.txt"; 
my %hash = (); #contig <-> family member with the maximal number of not uniquely mapped reads
my $contig; my $count=0; my $iso;

open(RESULTS_FILE, $results) or die "(sorting.pl) Impossible to open source file: $!\n";

	while (<RESULTS_FILE>) { 
		$contig = $_; chomp($contig);
		$count = `bash $LINKYX_PATH/scripts/number_of_not_uniq_mapped_reads.sh $dir $contig`;
		$hash{$contig} = $count; 
	}

foreach (sort {$hash{$a} <=> $hash{$b}} keys %hash) {
	  print "$_ $hash{$_}";

			#get sequence
			print "% of not uniquely mapped reads for contig ".$_.": ".$hash{$_}."\n";
			system("perl $LINKYX_PATH/scripts/get_sequences_based_on_ids.pl $reference $_");

			$iso = `bash $LINKYX_PATH/scripts/isoform_exists.sh $_ references/reference_contig_names`;
			if ($iso eq '') {
 				print "No isoforms.\n";
			} else {
				print "Isoforms exists:\n";
				print $iso;
			}
			print "----------\n";
	}