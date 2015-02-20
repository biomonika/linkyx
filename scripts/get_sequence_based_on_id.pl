use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
    $seq = $queryFile;
    chomp($queryFile);

    my $sequence = $db->seq($seq);
    if  (!defined( $sequence )) {
            die "Sequence $seq not found. \n" 
    }   
    print ">$seq\n", "$sequence\n";