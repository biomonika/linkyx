#!/usr/bin/perl 
#NOTE:  The following code is adapted from the post http://www.biostars.org/p/44867/#44906

use strict;
use warnings;

use Bio::DB::Sam;
open (MP, '>mpileup');

# all necessary inputs (I use Getopt::Long to obtain input parameters)
my $opt_i     = $ARGV[0]; #"mother.bam"; # should have an index file named input.bam.bai
my $opt_f     = $ARGV[1]; #"reference.fasta"; 
my $chr_id    = $ARGV[2]; #'comp1_c0_seq1'; # your chromosome-id
my $snp_pos = $ARGV[3]; # position you want to call all variants in this position 
# create the object
my $sam = Bio::DB::Sam->new( -fasta => $opt_f, -bam => $opt_i );

# get all reads that cover this SNP -- therefore, start and end set to $snp_pos
my @alignments = $sam->get_features_by_location(-seq_id => $chr_id, -start => $snp_pos, -end => $snp_pos );

	my $depth=@alignments;
	print MP "DEPTH:".$depth."\n";

my %res; # output hash that'll contain the count of each available nucleotide (or blank if the read covering the SNP is spliced in this position).
# loop over all alignments
for my $cAlign (@alignments) { 
    # parameters we'll need from our bam file

    my $start     = $cAlign->start;            # get start coordinate
    my $end       = $cAlign->end;               # get end coordinate
    my $ref_seq   = $cAlign->dna;               # get reference DNA sequence
    my $read_seq  = $cAlign->query->dna;        # get query DNA sequence

	#print $read_seq."\n";
	#print "start: ".$start." end: ".$end."\n";

    my $cigar     = $cAlign->cigar_str;     # get CIGAR sequence
    my $cigar_ref = $cAlign->cigar_array;       # probably the important useful of all. splits cigar to array of array reference

	#print $cigar."\n";                            # Ex: $cigar = 20M100N50M => $cigar_ref = [ [ 'M' 20 ] [ 'N' 100] ['M' 50] ]

    my $ref_cntr     = $start;       # $ref_cntr = assign start to counter variable for reference. This will ultimately
                                       # keep track of the current position on the chromosome
    my $read_cntr    = 0;           # $read_cntr = computes offset on the read
    my $read_snp_pos = "";              # variable to hold base at $snp_pos from current read

    foreach my $deref ( @$cigar_ref ) {
        my $cigar_chr = $deref->[0];        # cigar character (ex: M, N, I, D etc...)
        my $len_chr   = $deref->[1];        # number corresponding to `this` cigar character ( ex: 20, 100 )

        # NOTE: I => insertion in to the reference, meaning the read has the base(s) but the REFERENCE fasta does NOT
        #       D => deletion from the reference, meaning the READ does NOT have the base(s), but the reference does        

        # modify reference counter only on M, N or D occurrence
        if( $cigar_chr eq "M"  || $cigar_chr eq "N" || $cigar_chr eq "D") {
            $ref_cntr += $len_chr;
        }

	if( $cigar_chr eq "S" ) {
		#print "soft_clipped.\n";
	}

        # modify offset only on M or I occurrence
        if( $cigar_chr eq "M" || $cigar_chr eq "I" || $cigar_chr eq "S") {
            $read_cntr += $len_chr;
        }


        # now, check if the current operation takes ref_cntr past the SNP position. If it does,
        # 1) If the current operation is NOT "M", then its either "N" or "D". Both of them mean
        # that the read doesn't have a base at the SNP. So, this read is either spliced or has 
        # a deletion at that point and is not useful for your SNP location.
        # 2) If the current position IS "M", then the current operation has gotten is past SNP
        # location. So, we FIRST SUBTRACT what we added in this operation for ref_cntr and read_cntr, 
        # and then just add the difference ($snp_pos - $ref_cntr + 1)

	#print "read_cntr: ".$read_cntr."\nref_cntr: ".$ref_cntr." cigar: ".$cigar_chr."\n";

        if( $ref_cntr > $snp_pos ) { 
            if( $cigar_chr eq "M" || $cigar_chr eq "S") {
                $ref_cntr  -= $len_chr;
                $read_cntr -= $len_chr;
		
		#print "backing\n";
		#print "read_cntr: ".$read_cntr."\nref_cntr: ".$ref_cntr."\n";

		#IMPORTANT
                my $pos = $snp_pos - $ref_cntr + $read_cntr;
		
		#print "SNP position: ".$snp_pos;
		#print "position: ".$pos;
		
                $read_snp_pos = substr( $read_seq, $pos, 1 );

			#print "This is what we add: ".$read_snp_pos."\n";
            }
            # if $cigar_chr is "N" or "D", do nothing. $read_snp_pos is set to ""
	
	    if( $cigar_chr eq "N" || $cigar_chr eq "D") { 
		print MP "NAN\n";
	    }

            # add value of $read_snp_pos to hash and get out of loop - to next read
            $res{$read_snp_pos}++;
            last;
        }
    }
}

#Here, I am just printing to output.
#print "Count\n";
foreach my $key ( keys %res ) {
    print MP "$key:$res{$key}\n";
}
#print "--\n";
close MP;
