#!/usr/bin/perl
#Filter annotations by coverage of hmm

use strict; use warnings;

my $hmmerfile = $ARGV[0];
my $hmmlength = $ARGV[1];
open(HMM, '<', $hmmlength);
open(NHMMER, '<', $hmmerfile);

my %hmms;
while(<HMM>){
	my @sRNA = split(' ',$_);
	$hmms{$sRNA[0]} = $sRNA[1];
}

close HMM;

while(<NHMMER>){
	if ($_ =~/^\S+\s+nhmmer\s+(\S+)\s+(.*?)hmmfrom=(\d+);hmmto=(\d+)/){
		if(($4-$3) >= (0.7 * $hmms{$1})){
			print $_;
		}
	}
}

close NHMMER;


	
