#!/usr/bin/perl 

#use strict;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use List::Util qw( min max );

my $user_agent = LWP::UserAgent->new;

my $query_seq = $ARGV[0];
my $query_name = (split /\./, $query_seq)[0];
my $query_id = "";

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
        ($query_id = $gff_row[7 + $index]) =~ s/,//g;
        
        (my $start = $gff_row[1]) =~ s/[^0-9]//g;
        (my $end = $gff_row[3]) =~ s/[^0-9]//g;

        # No type 
        $result = $result."\tORF\t.\t$start\t$end\t.";
        $result = ($line =~ /reverse/i) ? $result."\t-" : $result."\t+";
        $result = $result."\t.\t$query_seq"; 
    
        print OUTPUT "$result\n" or die "File not open";

    }
}

close(ORF);

# do web blast => ./web_blast.pl blastx swissprot geobacter.fasta
# maybe make this a function?

my $encoded_query = "";

#foreach $query (@ARGV){ # possible extra feature => do stuff on multiple sequences
    open(QUERY,$query_seq);
    while(<QUERY>){
        $encoded_query = $encoded_query . uri_escape($_);
    } 
    close(QUERY);
#}

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
 
# poll for results
while (0){ # while true:
    sleep 5;
    
    $req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&F    ORMAT_OBJECT=SearchInfo&RID=$rid";
    $response = $user_agent->request($req);
    if ($response->content =~ /\s+Status=WAITING/m){
        print "Searching...\n";
        next;
    }
 
    if ($response->content =~ /\s+Status.+searching/im){
        print "Searching...\n";
        next;
    }
 
 
    if ($response->content =~ /\s+Status=FAILED/m){
        print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
    }
  
    if ($response->content =~ /\s+Status=UNKNOWN/m){
       print STDERR "Search $rid expired.\n";
       exit 3;
    }
 
    if ($response->content =~ /\s+Status=READY/m){
        if ($response->content =~ /\s+ThereAreHits=yes/m){
            #  print STDERR "Search complete, retrieving results...\n";
            last;
        } else {
             print STDERR "No hits found.\n";
             exit 2;
        }
    }

    # if we get here, something unexpected happened.
    exit 5;
} 

# retrieve and display results
$req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
$response = $user_agent->request($req);


# parse out BLAST results to GFF output
my $result = "";
my @query = ();
my @subject = ();
my $first_loop = 1; #HACK


foreach my $line (split(/\n/, $response->content)){

    chomp $line;
    my @line_split = split(/ +|, +/, $line);

    if ($line =~ /^>/gi){
    
        my $e_val = (split(/ /, $result))[2];
        #print "$e_val\n";   #if not ($result eq ""){#or (split(/ /, $result))[2] gt 1e-3){
        
        if ($e_val < 1e-3){
            # print result to GFF output
            $result = "$result ". min(@query) . " " . max(@query) . " " . min(@subject) . " " . max(@subject);
            #my $out = $query_id . "\t" . "BLASTX\tEST_match\t" . min(@query) . "\t" . max(@query) . "\t" .   
            my @result_split = split(/ /, $result);      
            my @out = ();
            push @out, ($query_id, "BLASTx", "EST_match", min(@query), max(@query), $result_split[1], $result_split[3], ".", "target=$result_split[0] " . min(@subject) . " " . max(@subject));
            #print "@out\n";
            
            $result = join("\t", @out);
            print OUTPUT "$result\n" if ($first_loop == 0); # HACK
            $first_loop = 0; # HACK
        }
        
        # reset results everytime
        ($result = $line_split[0]) =~ s/>//gi;
        @query = ();
        @subject = ();
    } elsif ($line =~ /score/gi){
        $result = "$result $line_split[3] $line_split[8]";
    } elsif ($line =~ /^query /gi){
        #print "===> $line_split[1] $line_split[3]";
        push(@query, $line_split[1]);     
        push(@query, $line_split[3]);
    } elsif ($line =~ /sbjct /gi){
        #print "<=== $line_split[1] $line_split[3]";
        push(@subject, $line_split[1]);
        push(@subject, $line_split[3]);       
    } elsif ($line =~ /frame/gi){
        (my $add = $line_split[3]) =~ s/[0-9]//gi;
        $result = "$result $add";
    }
}

close(OUTPUT);
