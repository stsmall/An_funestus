#!/usr/bin/perl

my $usage = "
# # # # # #
# GapDistrFromFasta.pl
# written by LinnÃ©a Smeds 19 oct 2010
# ======================================================
# Makes a list of all gap sequences (N:s) in a fasta 
# file, and also prints the length distribution of
# gaps. Can be run in \"fast\" mode, and then only return
# the total number of N:s in the input file.
# ======================================================
# Usage: perl NDistrFromFasta.pl <fastafile> <slow|fast>
#			<outpref>
#
# Example 1: perl NDistrFromFasta.pl mySeqs.fa fast
# 	(Returns the total number of Ns in file)
# Example 2: perl NDistrFromFasta.pl mySeqs.fa slow prefix
# 	(Returns prefix.gaps with a list of all gaps (start,
# 	stop and length), and a prefix.gaphist with the gap 
# 	size distribution: column1 = gap size,
#	column2 = #gap of this size)
";

use strict;
use warnings;

# Input parameters
my $fasta = $ARGV[0]; 
my $flag = $ARGV[1];
my $outpref = $ARGV[2];

my ($listOut, $histOut);

if($flag eq "slow"){
	if(defined $outpref) {
		$listOut = $outpref.".gaps";
		$histOut = $outpref.".gaphist";
	}
	else {
		$listOut = $fasta.".gaps";
		$histOut = $fasta.".gaphist";
	}
}
elsif($flag eq "fast") {
}
else {
	die "Flag must be either \"slow\" or \"fast\"\n\n$usage";
	
}

open(IN, $fasta);
if($flag eq "slow") {
	open(OUT, ">$listOut");
}
my %hist = ();
my ($head, $seq) = ("","");
my $sum = 0;
while(<IN>) {
	if($_ =~ m/^>/){
		my @tab = split(/\s+/, $_);
		$head = $tab[0];
		chomp($head);
		$head =~ s/>//;
		$seq = "";

		my $next = <IN>;
		while ($next !~ m/^>/) {
			chomp($next),
			$seq.= $next;
			if(eof(IN)) {
				last;
			}	
			$next = <IN>;
		}
		seek(IN, -length($next), 1);


		#Only sum up the number of N:s
		if($flag eq "fast") {		
			my @hits = $seq =~ m/(N+)/g;
			for(@hits) {
				$sum+=length($_);
			}	
		}
		
		#Print each gap with start, end and size to a .gaps file
		elsif($flag eq "slow") {
			my @seq = split(//,$seq);
			my $Ntemp = "";	
			my ($start, $end) = ("",""); 		
			for(my $i=0; $i<scalar(@seq); $i++) {
				if($seq[$i] eq 'N' || $seq[$i] eq 'n') {
					if($Ntemp eq "") {
						$start = $i+1;
					}
					$Ntemp .= 'N';
				}
				else {
					if($Ntemp ne "") {
						$end = $i;
						my $size = $end-$start+1;
						print OUT $head ."\t".$start."\t".$end."\t".$size."\n";
						$sum += $size;
						($Ntemp,$start,$end) = ("","","");
						if(defined $hist{$size}) {
							$hist{$size}++;
						}
						else {
							$hist{$size}=1;
						}
					}
				}
			}
		}
	}
}
close(IN);

# Print the distribution of gap sizes to a second file (.gaphist).
if($flag eq "slow") {
	close(OUT);
	open(OUT, ">$histOut");
	foreach my $key (sort {$a<=>$b} keys %hist) {
		print OUT $key."\t".$hist{$key}."\n";
	}
	close(OUT);
}

print "Total number of N: $sum bp\n";