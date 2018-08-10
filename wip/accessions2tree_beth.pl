#!/usr/bin/perl -w


use warnings;
use strict; 
use Getopt::Long;


my(
	$alignment,
	$species_summary,
	$help,
	$verbose
	);

&GetOptions(
	"a|alignment=s"	=>	\$alignment,
	"s|summary=s"	=>	\$species_summary,
	"h|help"	=>	\$help,
	"v|verbose"	=>	\$verbose
	);

if($help){
	&help();
	exit(1);
}


if (-t STDIN){
	print "Error: No accessions provided\n";
        &help();
        exit(1);
}

my @accessions;

open(NAMES, ">tmpNames.txt");
while(<STDIN>){
        my $input = $_;
        chomp $input;
        push @accessions, $input;
        print NAMES $_;
        }
close(NAMES);


unless ($accessions[0]){
        print "Error: No accessions provided\n";
        system("rm tmpNames.txt");
	&help();
        exit(1);
        }




if(not defined($alignment)){
	print "Error: Cannot find alignment file -a";
	&help();
	exit(1);
	}
elsif(not defined($species_summary)){
	print "Error: Cannot find annotation file -s";
	&help();
	exit(1);
	}
else{
	print "[@accessions]\n" if (defined($verbose)); 

        
	#grab seqs of accessions from alignment file
	my $egrepStr = join('|', @accessions);
	print "[$egrepStr]\n" if (defined($verbose));
		
	open(IN, "egrep \'$egrepStr\' $alignment | grep -v \'^\#' | "); 
	open(OUT, "> tmpOut.txt");
	print OUT "# STOCKHOLM 1.0\n\n";
	while(my $in=<IN>){
    				if($in=~/^([A-Z]+[0-9]+)\/\d+\-\d+\s+(\S+)/){
					my ($acc,$seq) = ($1,$2);
					print OUT "$acc	$seq\n";
				}
   
  				 print $in;
    
	}
	print OUT "\/\/\n";
	close(OUT);
	
	## generate tree with accession names
	
	system("esl-reformat --mingap phylip  tmpOut.txt > infile && rm -f outfile outtree && echo \"Y\" | dnadist ");
	system("mv outfile infile && echo \"Y\" | neighbor");
	
	#system("mv outtree intree && printf \"Y\nY\n\Y\nY\nY\nY\n\" | drawgram");
	

	## replace accessions in newick tree with ACCESSION_LONG_NAME (destroys old tree)
	## generates new file "formatted_tree" in easyfig format

        system("while read name; do grep \$name $species_summary | cut -f 8 | cut -d \",\" -f 1 | sed 's/complete genome//;s/complete genome.//;s/DNA//;s/genomic//;s/genome assembly//;s/ genome.//;s/main chromosome//;s/draft assembly//;s/ genome//;s/ sequence//;s/ sequence.//;s/main chrosome//;s/draft genome//;s/draft//' | sed 's/[^[:alpha:][:digit:] \t]/ /g' >> tmpSp.txt; done < tmpNames.txt && paste tmpNames.txt tmpSp.txt > tmpLookup.txt");
	open(OUTFILE, "> formatted_tree");
	print "Renaming tip labels...\n\n";
	my $ntax = @accessions;
        print OUTFILE "#NEXUS\nbegin taxa;\n\tdimensions\n\tntax=$ntax;\n\ttaxlabels\n";

	
	foreach my $acc (@accessions){
		my $sum = `grep ^$acc tmpLookup.txt | sed 's/[[:space:]]/\_/g'`;
		chomp($sum);
		print OUTFILE "\t$sum\n";
		system("sed -i 's/$acc/$sum/g' outtree");
		}
	
	print OUTFILE ";\nend;\n\nbegin trees;\n\ttree tree_1 = [&R] ";
	open(INFILE, "< outtree");
	
	while(<INFILE>){
		print OUTFILE $_;
		}
	
	print OUTFILE "end;";
	
	close(OUTFILE);
	close(INFILE);
	
	print "Renamed tree successfully written to formatted_tree\n\n"; 
	
	##cleanup of temp files
	#system("rm infile outfile outtree tmpLookup.txt tmpNames.txt tmpOut.txt tmpSp.txt");
	
	exit(0);
}

######################################################################

sub help{
	print "\n\n$0: Generates labelled phlyogenetic trees from a set of aligned sequences. Currently based on Phylip using DNAdist, NJ and produces Newick tree with Easyfig-readable tip labels.\n\nLabels produced in long format beginning with accession number in format ACCESSION_LONG_NAME \n\nUsage: $0 [accessions from <STDIN>] [options]\n\nOptions:\n	-h|--help	help\n	-a|--alignment	path to sequence alignment (unwrapped Stockholm format)\n        -s|--summary	path to summary text file\n\nsummary text file format example:\n\nEMBL_ACC        LEN_(MB)        LEN_(BPs)       NUM_GENES       GC_CONTENT      NCBI_TAXONOMY   SPECIES_NAME    DESCRIPTION
AB097150.1      0.0MB   5045    610     0.24    Bacteria; Tenericutes; Mollicutes; Acholeplasmatales; Acholeplasmataceae;Candidatus Phytoplasma; Candidatus Phytoplasma asteris.        Onion yellows phytoplasma       Onion yellows phytoplasma extrachromosomal DNA, complete sequence.\n\nTo-do:\n		- Add flag to read accession file?\n";
}
 
	
