#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my(
    @inFiles,
    $range,
    $name,
    $help
);

&GetOptions(
    "i|infiles=s@" => \@inFiles,
    "r|range=s" => \$range,
    "n|name=s" => \$name,
    "h|help" => \$help
);

if( $help ) {
    &help();
    exit(1);
}

if(defined $range && defined $name){
	print "-r and -n options not compatible\n";
	&help();
	exit(1);
}

if(defined $range){
	my($start,$stop)=split(",",$range);
	my($newStart,$newStop);	
	my $seqs;
	foreach my $i (@inFiles){
	    print $i,"\n";
	    my $awkCmd = "awk 'BEGIN{RS=\">\"}NR>";
	    	if((substr($start,0,1) eq "-")){
		   $seqs = "grep \"^>\" ".$i." | wc -l";
		   $seqs=`$seqs`;
		   $newStart = $seqs + $start;
		}
		else{
		    $newStart=$start;
		}
	        $awkCmd = $awkCmd.$newStart."{sub(\"\\n\",\"\\t\"); gsub(\"\\n\",\"\"); print RS\$0}' ".$i." | awk 'NR<";
		if((substr($stop,0,1) eq "-")){
		   $seqs = "grep \"^>\" ".$i." | wc -l";
		   $seqs=`$seqs`;
		   $newStop = $seqs + $stop - $newStart + 1;    
		}
		else{    
		    $newStop=$stop - $newStart + 1;
		}
		$awkCmd = $awkCmd.$newStop." {print \$1 \"\\n\" \$2}' ";
	        if($newStop < 1){
		    next;
		}
	    
	    else{
		print $awkCmd,"\n";
		system($awkCmd);
	    }
	}
}
