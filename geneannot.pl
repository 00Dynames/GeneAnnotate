#!/usr/bin/perl -w

use strict;

my $query_seq = $ARGV[0];
my $query_name = (split /\./, $query_seq)[0];

# open output file
open(OUTPUT, ">annotateout.gff") or die "Cannot open file\n";
print OUTPUT "##gff-version 3\n";


# Do ORF finding
system("getorf -sequence $query_seq -outseq $query_name.gff -osformat2 gff3
");

open(ORF, "<$query_name.gff") or die "Can't open getorf output";
	
# parse ORF file    
while(my $line = <ORF>){
    if($line =~ /^>/){
        # Set index offset for (REVERSE STRING)
        my $index = ($line =~ /reverse/i) ? 2 : 0;
        
        my @gff_row = split /\s/, $line;
        (my $result = $gff_row[7 + $index]) =~ s/,//g; # col1 => seqid
        
        (my $start = $gff_row[1]) =~ s/[^0-9]//g;
        (my $end = $gff_row[3]) =~ s/[^0-9]//g;

        # No type 
        $result = $result." ORF . $start $end .";
        $result = ($line =~ /reverse/i) ? $result." -" : $result." +";
        $result = $result." . $query_seq"; 
    
        print OUTPUT "$result\n" or die "File not open";

    }
}

close(ORF);

# do web blast




close(OUTPUT);
