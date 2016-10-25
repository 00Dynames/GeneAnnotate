#!/usr/bin/perl 

#use strict;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use List::Util qw( min max );

# Option flags
my $ps = 0;
my $f = 0;

if (($#ARGV + 1) < 1){ # Print usage 
    print "usage: geneannot.pl query\n";
    exit;
} else { # Check options
    foreach my $option (@ARGV){
        if ($option =~ /^-ps/gi){ # Run gff output through gff2ps to get ps output
            $ps = 1;
        }
        
        if ($option =~ /^-f/gi){ # Collect protein sequences of BLASTx hits in the fasta format
            $f = 1;
        }
        
        if ($option =~ /^-h/gi){ # Print help information
            print "=== Geneannotate ===\n Take a sequence in the fasta format and reutrn a gff file with associated ORF's and BLASTx matches\n\n -h => Print this help page\n -ps => Use gff2ps on gff output produced byu this script\n -f => Collect protein sequences of BLASTx hits in the fasta format\n\n";
            exit 0;                   
        }
        shift @ARGV;
    }
}

my $user_agent = LWP::UserAgent->new;
my $query_seq = $ARGV[0];
my $query_name = (split /\./, $query_seq)[0];
my $query_id = "";

# open output file
open(OUTPUT, ">annotate$query_name.gff") or die "Cannot open file\n";
print OUTPUT "##gff-version 3\n";

# Do ORF finding
system("getorf -sequence $query_seq -outseq $query_name.gff -osformat2 gff3
");

open(ORF, "<$query_name.gff") or die "Can't open getorf output";
	
# Parse ORF file    
while(my $line = <ORF>){
    if($line =~ /^>/){
        # Set index offset for (REVERSE STRING)
        my $index = ($line =~ /reverse/i) ? 2 : 0;
        
        my @gff_row = split /\s/, $line;
        (my $result = $gff_row[7 + $index]) =~ s/,//g; # col1 => seqid
        ($query_id = $gff_row[7 + $index]) =~ s/,//g;
        
        (my $start = $gff_row[1]) =~ s/[^0-9]//g;
        (my $end = $gff_row[3]) =~ s/[^0-9]//g;

        if ($start > $end){
            my $tmp = $start;
            $start = $end;
            $end = $tmp;
        }

        $result = $result."\tORF\t.\t$start\t$end\t.";
        $result = ($line =~ /reverse/i) ? $result."\t-" : $result."\t+";
        $result = $result."\t.\t$query_seq"; 
    
        print OUTPUT "$result\n" or die "File not open";
    }
}

close(ORF);

# Run web BLAST with BLASTx on Swissprot
my $encoded_query = "";    
open(QUERY,$query_seq);
while(<QUERY>){
    $encoded_query = $encoded_query . uri_escape($_);
} 
close(QUERY);

# build the request
my $args = "CMD=Put&PROGRAM=blastx&DATABASE=swissprot&QUERY=" . $encoded_query;
 
# make the request
my $req = new HTTP::Request POST => 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi';
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);

# get the response
my $response = $user_agent->request($req);

# parse out the request id
$response->content =~ /^    RID = (.*$)/m;
my $rid=$1;
#$rid = "" . $rid;

# parse out the estimated time to completion
$response->content =~ /^    RTOE = (.*$)/m;
my $rtoe=$1;
 
# wait for search to complete
sleep $rtoe;
 
# poll for results => Based on code from https://blast.ncbi.nlm.nih.gov/docs/web_blast.pl
while (0){ # while true:
    #sleep 5;
    
    $req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&F    ORMAT_OBJECT=SearchInfo&RID=$rid";
    $response = $user_agent->request($req);
    if ($response->content =~ /Status=WAITING/m){
        print "Searching...\n";
        next;
    }
 
    if ($response->content =~ /\s+Status=FAILED/m){
        print "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
    }
  
    if ($response->content =~ /\s+Status=UNKNOWN/m){
       print STDERR "Search $rid expired.\n";
       exit 3;
    }
 
    if ($response->content =~ /\s+Status=READY/m){
        if ($response->content =~ /\s+ThereAreHits=yes/m){
            print "Search complete, retrieving results...\n";
            last;
        } else {
             print "No hits found.\n";
             exit 2;
        }
    }

    # if we get here, something unexpected happened.
    exit 5;
} 

# retrieve and display results
$req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
$response = $user_agent->request($req);
#print $response->content;

# parse out BLAST results to GFF output
my $result = "";
my @query = (); # Query sequence indexes
my @subject = (); # Subject sequence indexes
my $subject_seq = "";
my $first_loop = 1; #HACK


foreach my $line (split(/\n/, $response->content)){

    chomp $line;
    my @line_split = split(/ +|, +/, $line);

    if ($line =~ /^>/gi){
        my $e_val = (split(/ /, $result))[2];
        
        if ($e_val < 1e-3){ # Print result to GFF output
            $result = "$result ". min(@query) . " " . max(@query) . " " . min(@subject) . " " . max(@subject);
            my @result_split = split(/ /, $result);      
            my @out = ();
            push @out, ($query_id, "BLASTx", "EST_match", min(@query), max(@query), $result_split[1], $result_split[3], ".", "target=$result_split[0] " . min(@subject) . " " . max(@subject));
            
            $result = join("\t", @out);
            print OUTPUT "$result\n" if ($first_loop == 0); # HACK
            $subject_seq =~ s/\n//g;
            print OUTPUT "###FASTA\n>$result_split[0]\n$subject_seq\n" if ($f == 1);
            $first_loop = 0; # HACK
        }
        
        # Reset results everytime there is a new match
        ($result = $line_split[0]) =~ s/>//gi;
        @query = ();
        @subject = ();
        $subject_seq = "";
    } elsif ($line =~ /score/gi){
        $result = "$result $line_split[3] $line_split[8]";
    } elsif ($line =~ /^query /gi){ # Collect query indexes
        push(@query, $line_split[1]);     
        push(@query, $line_split[3]);
    } elsif ($line =~ /sbjct /gi){ # Collect subject indexes
        push(@subject, $line_split[1]);
        $subject_seq = $subject_seq . $line_split[2];
        push(@subject, $line_split[3]);       
    } elsif ($line =~ /frame/gi){ # Parse strand
        (my $add = $line_split[3]) =~ s/[0-9]//gi;
        $result = "$result $add";
    }
}

if ($ps == 1){ # Execute ps optional action
    `gff2ps -v annotate$query_name.gff 1> annotate$query_name.ps 2> STDIN`
}

# Remove intermediary files
`rm $query_name.gff STDIN`;

close(OUTPUT);
