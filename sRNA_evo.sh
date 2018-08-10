#!/bin/sh

usage(){
    echo "
    	 Wrapper script for running analysis pipeline for Salmonella sRNA conservation project.
	 Usage: sRNA_evo.sh -d [data directory] -g [starting gff] -f [starting genome.fasta] -1 [training database multi-fasta] -2 [target database multi-fasta]
	 "
}

GFF=""
GENOME=""
DB1=""
DB2=""
DIR=""


#get options 

while getopts :g:f:d:1:2:h opt
    do
        case "${opt}" in
	    h  ) usage;exit;;
	    g  ) GFF=${OPTARG};;
	    f  ) GENOME=${OPTARG};;
	    d  ) DIR=${OPTARG};;
	    1  ) DB1=${OPTARG};;
	    2  ) DB2=${OPTARG};;
	    \? ) echo "Unknown option: -${OPTARG}" >&2; exit 1;;
	    :  ) echo "Missing option argument for -${OPTARG}" >&2; exit 1;;
	    *  ) echo "Unimplemented option: -${OPTARG}" >&2; exit;;
	esac
done
shift $((${OPTIND}-1))

setup(){
    if [ -z ${GENOME} ] || [ -z ${GFF} ] || [ -z ${DB1} ] || [ -z ${DB2} ] || [ -z ${DIR} ] ; then
	echo "Error: Missing input(s)." >&2
	usage
	exit 1
    fi

    #check dependencies are installed
    declare -a dependencies=(
	"nhmmer"
	"esl-sfetch"
	"esl-seqstat"
	"esl-alimanip"
	"hmmbuild"
	"makehmmerdb"
	"mafft-qinsi"
	"sreformat"
	"parallel"
    )

    for i in "${dependencies[@]}"
        do
            if ! [ -x "$(command -v $i)" ]
                then
                    echo "Error: Couldn't find $i in PATH. Is $i installed? Check README for required dependencies." >&2
                    exit 1
            fi
    done

    #check FASTA files are ok and copy to working directory
    mkdir -v -p ${DIR}/hmms/training/flank ${DIR}/nhmmer/final ${DIR}/hmms/starting/flank ${DIR}/QC ${DIR}/data 

    for i in ${GENOME} ${DB1} ${DB2}; do
	if ! [[ $(esl-seqstat $i 2> /dev/null) ]]; then
	    echo "Error: $i is not a valid FASTA file." >&2
	    usage
	    exit 1
	else
	    cp --verbose $i ${DIR}/data/	
	fi
    done

    cut -f3 $GFF | sort | uniq > ${DIR}/data/sRNA_list.txt
}

buildHmms(){
    FASTA=$1
    if [[ $(esl-seqstat ${FASTA}.fasta | grep Number | awk '{print $4}') > 1 ]]; then
        mafft-qinsi --retree 1 --thread 20 --quiet ${FASTA}.fasta > ${FASTA}.maf
	esl-reformat clustal ${FASTA}.maf > ${FASTA}.clustal
	hmmbuild ${FASTA}.hmm ${FASTA}.clustal >/dev/null
    else
        sreformat -u stockholm ${FASTA}.fasta > ${FASTA}.stk
	hmmbuild ${FASTA}.hmm ${FASTA}.stk >/dev/null
    fi
}

getSeqs(){
    echo "Indexing database (this may take a while)..."
    esl-sfetch --index $2 >/dev/null

    echo -n "Fetching sequences and building hmms..."
    while read i; do
	if [[ $(grep -E "\s$i\s" $1 | fetchGff.pl -i $2 > ${DIR}/hmms/$3/$i.fasta) -ne 0 ]] ; then
	    echo "Error: Couldn't fetch all sequences. Does your genome accession match FASTA header? (i.e CP002487 vs CP002487.1) \n"
	    usage
	    exit 1
	else
	    ###Grab 150nt of sequence either side of genes, build 1-sequence HMMs for filtering nhmmer results
            grep -E "\s+$i\s+" $1 | awk -v OFS='\t' '{ if ($7 =="+") {$3 = $3"_1.fasta"; $5=$4; $4 = $4-150; print} else {$3 = $3"_1.fasta"; $4=$5; $5 = $5+150; print}}' | fetchGff.pl -i $2 > ${DIR}/hmms/$3/flank/$i\_1.fasta 
            grep -E "\s+$i\s+" $1 | awk -v OFS='\t' '{ if ($7 =="+") {$3 = $3"_2.fasta"; $4=$5; $5 = $5+150; print} else {$3 = $3"_2.fasta"; $5=$4; $4 = $4-150; print}}' | fetchGff.pl -i $2 > ${DIR}/hmms/$3/flank/$i\_2.fasta
	    buildHmms  ${DIR}/hmms/$3/$i
	    buildHmms  ${DIR}/hmms/$3/flank/$i\_1 
	    buildHmms  ${DIR}/hmms/$3/flank/$i\_2 
	fi
    done < ${DIR}/data/sRNA_list.txt
    
    cat ${DIR}/hmms/$3/*.hmm > ${DIR}/hmms/$3Seqs.hmm
    cat ${DIR}/hmms/$3/flank/*.hmm > ${DIR}/hmms/$3Flank.hmm
    echo "Done!"
}

bigSearch(){
    echo -n "Setting up hmmer database (this may take a while)..."
    makehmmerdb $1 $1.db >/dev/null
    echo "Done!"
    
    JOBS=$5
    echo "Running nhmmer with ${JOBS} jobs."

    runHmmer.pl -e $4 -j ${JOBS} -m nhmmer -db $1 -q ${DIR}/hmms/$2Seqs.hmm -o ${DIR}/nhmmer/$3
    mv ${DIR}/nhmmer/$3/nhmmer*.gff ${DIR}/nhmmer/$2Seqs.gff
    runHmmer.pl -e 0.001 -j ${JOBS} -m nhmmer -db $1 -q ${DIR}/hmms/$2Flank.hmm -o ${DIR}/nhmmer/$3
    mv ${DIR}/nhmmer/$3/nhmmer*.gff ${DIR}/nhmmer/$2Flank.gff
    echo "Done!"
}

filterGff(){
    echo -n "Filtering sequences.."
    ## Remove results with low hmm coverage
    ls ${DIR}/hmms/$1/*.hmm | parallel 'hmmstat {}' | grep "^1" | awk '{print $2"\t"$6}' > ${DIR}/data/hmm_length.txt
    ls ${DIR}/hmms/$1/flank/*.hmm | parallel 'hmmstat {}' | grep "^1" | awk '{print $2"\t"$6}' > ${DIR}/data/flank_hmm_length.txt
    hmm_coverage.pl ${DIR}/nhmmer/$1Seqs.gff ${DIR}/data/hmm_length.txt > ${DIR}/nhmmer/tmp && mv ${DIR}/nhmmer/tmp ${DIR}/nhmmer/$1Seqs.gff
    hmm_coverage.pl ${DIR}/nhmmer/$1Flank.gff ${DIR}/data/flank_hmm_length.txt > ${DIR}/nhmmer/tmp && mv ${DIR}/nhmmer/tmp ${DIR}/nhmmer/$1Flank.gff
    
    ##Filter results with flanking regions and bitscore
    filter.pl -1 ${DIR}/nhmmer/$1Seqs.gff -2 ${DIR}/nhmmer/$1Flank.gff > ${DIR}/data/FILT_$1Seqs.gff

    echo "Done!"
}

#setup directories and check dependencies
setup


GENOME=$(basename ${GENOME})
DB1=$(basename ${DB1})
DB2=$(basename ${DB2})
echo -e "Starting with files: \n${GFF}\n${GENOME}\n${DB1}\n${DB2}\n"

#Iteration 1:
#    Grab starting sequences using GFF and GENOME inputs and build 1-sequence HMMs for model training
getSeqs ${GFF} ${GENOME} starting

#    Run starting sequence and flanking sequence hmms over training database with nhmmer and parse results to gff format. Split into multiple jobs if DB1 is very large using GNU parallel
bigSearch ${DIR}/data/${DB1} starting training 0.001 10

#    Filter results, grab sequences and build new models
filterGff starting

#Iteration 2:
getSeqs ${DIR}/data/FILT_startingSeqs.gff ${DIR}/data/${DB1} training
bigSearch ${DIR}/data/${DB2} training final 0.01 30
filterGff training
