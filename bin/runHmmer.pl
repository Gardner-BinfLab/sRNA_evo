#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use POSIX;
########
#TO-DO
#Check for install of easel/squid
#Implement restart
########

my( 
	@inFiles,
	$method,
	$jobs,
	$queryFile,
	$restart,
	$eValue,
	$nice,
	$dirty,
	$verbose,
	$nosplit,
	$help,
        $outDir
    );

$jobs = 1;
$eValue = 0.01;

&GetOptions( 
	"db|infiles=s@"	=> \@inFiles,
	"q|qfile=s"	=> \$queryFile,
	"m|method=s"	=> \$method,
	"j|jobs=i"	=> \$jobs,
	"e|evalue=f"	=> \$eValue,
        "o|out=s"       => \$outDir,
	"d|dirty"	=> \$dirty,
	"v|verbose"	=> \$verbose,
	"x|nosplit"	=> \$nosplit,
	"h|help"	=> \$help,
	"n|nice"	=> \$nice
    );

if( $help ) {
    &help();
    exit(1);
}

print $jobs,"\t",$eValue,"\n";
### Remove or alter this if statement if you want more CPUS. Set to 36 to be nice to the Otago compute server. :)
if($jobs > 36){
	print "Maximum 36 jobs allowed. Mess with the script if you want more CPUs\n";
	&help();
	exit(1);
} 

if(not defined($restart)){
	if(not -e $inFiles[0]){
		print "FATAL: missing fasta files, -db is required!\n";
		&help();
		exit(1);
	}  
	system("mkdir -p $outDir") and die("ERROR: FAILED TO RUN [mkdir -p $outDir]\n[$!]"); 
	if(not defined($nosplit)){
		my $newInFile  = catFiles(\@inFiles);
		my ($fileString, $dbsize, $numSeqs)=splitSeq($newInFile);
		runMethod($method,$fileString,$queryFile,$dbsize,$numSeqs);
	}
	else{
		my $fileString = join(" ",@inFiles);
		runMethod($method,$fileString, $queryFile);
	}
	parse($method);
}

##Cleanup all the files:
if (not defined($dirty)){
    system("rm $outDir/*.fa $outDir/*.tblout $outDir/*.out") and die "FATAL: failed to clean up [$outDir] files\n[$!]";
}

exit(0); 

############################################################################################
sub runMethod{
	my ($method, $fileString, $queryFile, $dbsize,$numSeqs)=@_; 
	my $opt;
	print $method,".\n";
	if (defined $verbose){
		print "Running: [$method] on [splitSeq($jobs,$outDir/$$\.*\.fa";
		print ",$queryFile" if (defined $queryFile);
		print ")]\n";
		}
	if ($method eq "nhmmer"){ 
		my $hmmfile = $queryFile;
		###Input QC - is the query file a valid HMM?
		if($queryFile !~ /fa$/i and $queryFile!~/fasta$/i and $queryFile !~ /hmm$/i){
			$hmmfile = "$outDir/tmp\.stk";
			system("esl-reformat pfam $queryFile > $hmmfile") and die "FATAL: failed to execute [esl-reformat pfam $queryFile > $hmmfile] (query appears to be a fasta file)\n[$!]"; 
			}
		if($queryFile !~ /hmm$/i){
			system("hmmbuild $queryFile\.hmm $queryFile > /dev/null") and die "FATAL: failed to execute [hmmbuild $queryFile\.hmm $queryFile]\n[$!]"; 
			$hmmfile = "$queryFile\.hmm";
			system("hmmstat $hmmfile > /dev/null") and die "FATAL: failed to execute [hmmstat $hmmfile]\n[$!]";
			}
		###Options for different methods
		if($method eq "nhmmer"){
		    $opt = "-E $eValue -A {}.msa --tblout {}.tblout -o {}.out ";
			$dbsize = $dbsize/10**6;
		    	if($dbsize>0){
				$opt = "-Z $dbsize ".$opt;
	    		}
		}		
	}
	my $paraCmd = "$method $opt $queryFile {}";
	submitJobs($paraCmd,$fileString);
}

######################################################################
#parse: parse method results to gff
###########################################################
sub parse{
    my @outFiles = glob("$outDir/*.tblout"); 
    my $outFile="$outDir/nhmmer_$$.gff";
    my $fileString = join(" ",@inFiles);
    open(OUT, "> $outFile"); 
    print OUT "##gff-version 3\n##$method results on [$queryFile] vs [$fileString]\n";
    foreach my $File (@outFiles){
	open(IN, "< $File"); 
	my ($name,$comment,$score,$f,$t,$strand)=("","",0,0,0,"+");
	while(my $in=<IN>){
	    if($method eq "nhmmer"){
			chomp($in);
			next if ($in=~/^#/);
           		my @preds = split(/\s+/, $in); 
            		if(defined($preds[6]) and defined($preds[7]) and $preds[6]=~/^\d+$/ and $preds[7]=~/^\d+$/){
                		$name=$preds[0];
                		$score=$preds[13];
                		$strand=$preds[11];
                		($f,$t)=($preds[6],$preds[7]);
                		if($f>$t){
                    			($f,$t)=($t,$f);
                		}		
                		$comment.="E-value=$preds[12];query=$preds[2];hmmfrom=$preds[4];hmmto=$preds[5];";
                		print OUT "$name\t$method\t$preds[2]\t$f\t$t\t$score\t$strand\t\.\t$comment\n";
                		($comment,$score,$f,$t,$strand)=("",0,0,0,"+");
            		}
		}
		elsif($method eq "hmmsearch"){
			chomp($in);
			next if ($in=~/^#/);
       			my @preds = split(/\s+/, $in); 
			if(defined($preds[17]) and defined($preds[18]) and $preds[17]=~/^\d+$/ and $preds[18]=~/^\d+$/){
 				$name=$preds[0];
               		  	$score=$preds[7];
               		  	$strand='+';
                          	($f,$t)=($preds[17],$preds[18]);
               	             	if($f>$t){
                                	($f,$t)=($t,$f);
                                	}
                                $comment.="E-value=$preds[6];query=$preds[3];hmmfrom=$preds[15];hmmto=$preds[16];";
                                print OUT "$name\t$method\t$preds[3]\t$f\t$t\t$score\t$strand\t\.\t$comment\n";
                                ($comment,$score,$f,$t,$strand)=("",0,0,0,"+");
                           }
		}
	    }
	close(IN);
	}
   close(OUT);
}

######################################################################
#submitJobs: Run method with GNU parallel
sub submitJobs {
	my ($paraCmd,$input) = @_;
	my $pOpt =" --jobs $jobs ";
	if (defined ($nice)){
		$pOpt = $pOpt." --nice 19 ";
	}
	if (defined ($verbose)){
		$pOpt = $pOpt." --eta ";
	}
	my $paraExe = "parallel $pOpt '$paraCmd' ::: $input";
#	print $paraExe,"\n";
	system($paraExe);
}

#####################################################################
#Split up large sequence files into chunks

sub splitSeq {  
    my ($inFile) = @_; 
    print "Splitting $inFile into $jobs chunks\n";
    #initialise 
    my %newInFiles;
    for (my $i=1; $i<$jobs+1; $i++){
	$newInFiles{"$outDir/$$\.$i\.fa"}=0;
    }
    
    my $seqStat = seqStat($inFile); 
    my ($fileCount,$seqCount)=(1,0);
    open(IN,  "< $inFile");
    open(OUT, "> $outDir/$$\.$fileCount\.fa" );
    while(my $fa = <IN>){	
	if($fa =~ /^>/){
	    $seqCount++;
	    close(OUT);
	    open(OUT, ">> $outDir/$$\.$fileCount\.fa" );
	    $fileCount++;
	    $fileCount=1 if ($fileCount>$jobs);
	    $newInFiles{"$outDir/$$\.$fileCount\.fa"}++;
	}	
	print OUT $fa;     
   }    
    close(OUT);
    close(IN);
    
    my $numSeqs = $seqStat->{'numSeqs'}; 
    my $dbsize = $seqStat->{'numRes'};
    printf "WARNING: the number of printed sequences [$seqCount] is not equal to the number of counted sequences [%d]\n", $seqStat->{'numSeqs'} if ($seqStat->{'numSeqs'} ne $seqCount);
    printf "WARNING: created too many files [$fileCount] should be [$jobs]\n" if ($fileCount>$jobs);
    my $fileString = "";
    foreach my $newFile (keys %newInFiles){
	$fileString = $fileString." ".$newFile;
	}
    return ($fileString, $dbsize, $numSeqs); 
}
#######################################################################
sub seqStat {
	my $inFile = shift; 
	my %seqStat;  
	#Find number of sequences and total sequence length:
	open(SEQSTAT, "esl-seqstat $inFile |")  or die "FATAL: failed to open pipe: [seqstat -a $inFile]\n[$!]";
	while (my $stat = <SEQSTAT>){
		if($stat =~ /^Number.*\s+(\d+)$/){
			$seqStat{'numSeqs'}=$1;
			printf "[%d] seqs\n", $seqStat{'numSeqs'};
			}
		elsif ($stat =~ /^Total.*\s+(\d+)$/){
		    $seqStat{'numRes'}=$1;
		}
	    }
	close(SEQSTAT);
	return \%seqStat; 
}

######################################################################
#catFiles: concatenate a lot of seq Files
sub catFiles {
    my $inFiles = shift; 
    my $inFile = "$outDir/" . $$ . ".fasta"; 
    
    #TODO: CHECK FILE NAMES ARE UNIQUE
    #concatenate files
    foreach my $in (@{$inFiles}){
	system("sreformat -r -u fasta $in >> $inFile") and die "FATAL: failed to run [sreformat -r -u fasta $in >> $inFile]!\n[$!]"; 
    }
    
    return $inFile;    
}

###################################################################
sub help {
    print STDERR <<EOF;
    
$0: runs nhmmer over large sequence databases using GNU parallel 
Usage: $0 fasta_file(s)


Script options:
	       -h |--help	show this help
               -d |--dirty	Don't cleanup temporary files
               -v |--verbose    Print lots of stuff
	       -x |--nosplit	Don't split seqs (just run my method in parallel and parse it for me)
	       -o |--out        Output directory
Method options:
               -db|--infiles	<str>  Sequence database (Fasta files)
               -m |--method	<str>  Run the given method. (Currently supported: nhmmer)
	       -q |--qfile	<str>  Query file used to search sequence database
	       -e |--evalue	<str>  E-value setting [DEFAULT 0.01]

Job options:
		-j |--jobs	<int>  Total number of jobs to run in parallel - will split sequence database into this many chunks. [DEFAULT:1, Max: 36]
       	        -n |--nice	Set job priority to very nice (19)
		-v |--verbose	Monitor job status

Examples:
runHmmer.pl -d -n -e 0.001 -j 30 -m nhmmer -db [*.fasta] -q [*.hmm] -o data/nhmmer

EOF
}

