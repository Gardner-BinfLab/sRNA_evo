#!/usr/bin/perl
#a script to stop me screwing up esl-sfetch 
#takes gff from STDIN, and pre-built esl-sfetch index
#-g names sequences "seqname_feature" for multi-genome, multi-gene gff files 
#dies if seqfetch fails (numbers out of range)

use strict; use warnings;
use Getopt::Long;

my ($help,$name,$index);
&GetOptions(
	"h|help"	=>	\$help,
	"g|genome"	=>	\$name,
	"i|index=s"	=>	\$index
	);

if($help){
	&help();
	exit(1);
}


if (-t STDIN){
	print "Error: No input from STDIN\n";
	&help();
        exit(1);
}



### Check gff
my $errCount = 0;
my $fetches =0;

while(<STDIN>){
	if ($_ =~/^\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t[+-]\t\S+\t\S+/){
		my $line = $_;
		my $err;
		my @columns = split('\t',$line);
		if($name){
			$columns[2] = $columns[0]."\_".$columns[2];
			}
		if($columns[6] eq "+"){
			my $cmd = "esl-sfetch -n $columns[2] -c $columns[3]\.\.$columns[4] $index $columns[0]";
			if (system($cmd) > 0){
				$err =1;
				}
			}
		 elsif($columns[6] eq "-"){
			my $cmd = "esl-sfetch -r -n $columns[2] -c $columns[3]\.\.$columns[4] $index $columns[0]";
			if (system($cmd) > 0){
				$err=1;
				}
			}
		if($err){
			$errCount ++;
		}
		$fetches ++;
	}
}

if ($errCount){
	die "$errCount out of $fetches attempted fetches failed. Something wrong with your gff?\n";
}

sub help{
	print "fetchGff.pl\nUsage [Gff from STDIN} | fetchGff.pl index_file\n\nRequires pre-built esl-sfetch index.\n\nOptions:\n-------\n-h --help\tDisplay this help\n\n";
}
