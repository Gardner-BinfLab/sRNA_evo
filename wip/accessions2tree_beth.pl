#!/usr/bin/perl -w
#accessions2tree.pl - takes a list of bacterial genome accessions and generates a phylogenetic tree with Phylip from a pre-computed alignment of 16S sequences
#output format is a newick tree formatted for use with Figtree

use warnings;
use strict; 
use Getopt::Long;
use Cwd qw(cwd);

my(
	$alignment,
	$species_summary,
	$help,
    $verbose,
    $outdir,
    @accessions
	);

my $dir=cwd;

&GetOptions(
    "a|alignment=s"	=>	\$alignment,
    "s|summary=s"	=>	\$species_summary,
    "h|help"	=>	\$help,
    "v|verbose"	=>	\$verbose
    "o|output" => \$outdir
	);

if($help){
	&help();
	exit(1);
}

sub setup{
    if (-t STDIN){
	print "Error: No accessions provided\n";
        &help();
        exit(1);
    }

    system("mkdir -p $dir/$output");
    while(<STDIN>){
        my $input = $_;
        chomp $input;
        push @accessions, $input;
    }

    unless ($accessions[0]){
        print "Error: No accessions provided\n";
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
}

sub grabSeqs{
	print "[@accessions]\n" if (defined($verbose)); 
       
	#grab seqs of accessions from alignment file
	my $egrepStr = join('|', @accessions);
	print "[$egrepStr]\n" if (defined($verbose));		
	open(IN, "egrep \'$egrepStr\' $alignment | grep -v \'^\#' | "); 
	open(OUT, ">  $dir/$output/$$_tmpOut.txt");
	#parse CMalign stockholm
	print OUT "# STOCKHOLM 1.0\n\n";
	while(my $in=<IN>){
	    #Expected file format: ACCESSION/from-to Aligned_sequence
    				if($in=~/^(\S+)\/\d+\-\d+\s+(\S+)/){
					my ($acc,$seq) = ($1,$2);
					print OUT "$acc	$seq\n";
				}
   
  				 print $in;
    
	}
	print OUT "\/\/\n";
	close(OUT);
}	


sub makeTree{
    ## generate tree from new alignment
    system("esl-reformat --mingap phylip $dir/$output/$$_tmpOut.txt > infile && echo \"Y\" | dnadist ");
    system("mv outfile infile && echo \"Y\" | neighbor");

    ## replace accessions in newick tree with ACCESSION_LONG_NAME (destroys old tree)
    ## generates new file "formatted_tree" in easyfig format
    open(OUTFILE, "> formatted_tree");
    print "Renaming tip labels...\n\n";
    my $ntax = @accessions;
    my %newNames;
    print OUTFILE "#NEXUS\nbegin taxa;\n\tdimensions\n\tntax=$ntax;\n\ttaxlabels\n";
    for(my $i=0,$i<$ntax,$i++){
	my $acc = $accessions[$i];
	my $spName = `grep $acc $species_summary | rev | cut -f2 | rev`
	chomp($spName);
	$spName = $acc,"\t",$spName
	$newNames{$acc} = $spName 
    }
    
    print OUTFILE ";\nend;\n\nbegin trees;\n\ttree tree_1 = [&R] ";
    open(INFILE, "< outtree");
    my $infile = "";
    while(<INFILE>){
	my $line = $_;
	$infile = $infile.$line;
    }
    close(INFILE);
    for my $acc in keys($newNames){
	my $name = $newNames{$acc};
	$infile =~ 's/$acc/$name/';
    }
    print OUTFILE $infile;
    print OUTFILE "end;";
    close(OUTFILE);
    
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
 
	
