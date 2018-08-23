#!/bin/sh
#Downloads fastq files from SRA, trims, maps and generates plotfiles for visualisation in artemis

#Dependencies:
#curl
#sratoolkit
#samtools 1.6 (older versions may not work for generating plotfiles)
#bowtie2
#trimmomatic 0.36

usage(){
    echo "sra2plot.sh is a wrapper script for downloading, mapping and visualising RNA-seq data from the NCBI Sequence Read Archive (SRA). Currently assumes paired end reads with TruSeq3 adaptors. Path for trimmomatic needs to be set to run. 
Usage sra2plot [opts] [input]
    
    Options:
		-h	Display this help

		Input	       
	       	-r	Reference genome accession (required)
		-s	SRA run accession or name of split fastq files (Required. Format: FILE_1.fastq FILE_2.fastq)
    		-n	Number of cores

		Turn off defaults
		-d	Turn off download. Default: download genome and SRA from NCBI if not found in working directory. 
			(Genome accession must be in Genbank nucleotide format: https://www.ncbi.nlm.nih.gov/Sequin/acc.html)
		-t	Turn off trimming
		-m	Turn off mapping
		-p	Don't make plotfiles
		-x	Don't cleanup files"
}
TPATH= #path for trimmomatic (no quotes)
OUTDIR=""
SRA=""
GENOME=""
TRIM=true
MAP=true
PLOT=true
CLEAN=true
DOWNLOAD=true
THREADS=1

while getopts :s:r:n:thdmpx opt; do
    case "${opt}" in
	h) usage;exit;;
	t) TRIM=false;;
	s) SRA=${OPTARG};;
	r) GENOME=${OPTARG};;
	d) DOWNLOAD=false;;
	m) MAP=false;;
	p) PLOT=false;;
	x) CLEAN=false;;
	n) THREADS=${OPTARG};;
	\?) echo "Unknown option: -${OPTARG}" >&2; exit 1;;
	:) echo "Missing option argument for -${OPTARG}" >&2; exit 1;;
	*) echo "Unimplemented option: -${OPTARG}" >&2; exit;;
    esac
done
shift $((${OPTIND}-1))

if [ -z ${GENOME} ] || [ -z ${SRA} ]; then
    echo "Error: No input specified." >&2
    usage
    exit 1
fi

if [ -z ${TPATH} ]; then
    echo "Error: Path to trimmomatic install folder is not set.\n" >&2
    exit 1
fi

if $DOWNLOAD;then
    if [ ! -f ${GENOME}.fasta ];then
	curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=${GENOME}&amp;rettype=fasta" > ${GENOME}.fasta
    fi
    if [ ! -f ${SRA}_*.fastq ];then
	fastq-dump --split-3 ${SRA}
    fi
fi

if $TRIM;then 
    if [ ! -d trimmed ];then
	mkdir trimmed
    fi
    java -jar ${TPATH}/trimmomatic-0.36.jar PE -threads `echo $((2*${THREADS}))` ${SRA}_1.fastq ${SRA}_2.fastq trimmed/${SRA}_1_paired.fastq trimmed/${SRA}_1_unpaired.fastq trimmed/${SRA}_2_paired.fastq trimmed/${SRA}_2_unpaired.fastq ILLUMINACLIP:${TPATH}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

if $MAP; then
    #Build index of genome if necessary
    if [ ! -d index ]; then
	bowtie2-build ${GENOME}.fasta ${GENOME} &&
	mkdir index && mv *.bt2* index/
    fi
    bowtie2 -p `echo "$((2*${THREADS}))"` -x index/${GENOME} -1 trimmed/${SRA}_1_paired.fastq -2 trimmed/${SRA}_2_paired.fastq -S ${SRA}.sam
fi

if $PLOT;then
    samtools view -bS -@ ${THREADS} ${SRA}.sam > ${SRA}.bam
    samtools sort -@ ${THREADS} ${SRA}.bam > ${SRA}.sorted.bam
    # Forward strand.
    #alignments of the second in pair if they map to the forward strand
    samtools view -b -f 128 -F 16 -@ ${THREADS} ${SRA}.sorted.bam > ${SRA}.fwd1.bam
    samtools index ${SRA}.fwd1.bam
    #alignments of the first in pair if they map to the reverse strand
    samtools view -b -f 80 -@ ${THREADS} ${SRA}.sorted.bam > ${SRA}.fwd2.bam
    samtools index ${SRA}.fwd2.bam
    #combine alignments that originate on the forward strand
    samtools merge -f $DATA.fwd.bam $DATA.fwd1.bam $DATA.fwd2.bam
    samtools index $DATA.fwd.bam

    # Reverse strand
    #alignments of the second in pair if they map to the reverse strand
    samtools view -b -f 144 -@ ${THREADS} ${SRA}.sorted.bam > ${SRA}.rev1.bam
    samtools index ${SRA}.rev1.bam
    #alignments of the first in pair if they map to the forward strand
    samtools view -b -f 64 -F 16 -@ ${THREADS} ${SRA}.sorted.bam > ${SRA}.rev2.bam
    samtools index ${SRA}.rev2.bam
    #combine alignments that originate on the reverse strand.
    samtools merge -f ${SRA}.rev.bam ${SRA}.rev1.bam ${SRA}.rev2.bam
    samtools index ${SRA}.rev.bam

    #Generate plotfiles
    samtools mpileup -aa ${SRA}.fwd.bam > ${SRA}.fwd.mpileup
    samtools mpileup -aa ${SRA}.rev.bam > ${SRA}.rev.mpileup
    cat ${SRA}.fwd.mpileup | cut -f4 > ${SRA}.fwd.plot
    cat ${SRA}.rev.mpileup | cut -f4 > ${SRA}.rev.plot
    paste ${SRA}.rev.plot ${SRA}.fwd.plot > ${SRA}.plot   
fi

if $CLEAN; then
    rm *.bam *.mpileup *.bai
	if [ ! -d fastq ]; then
	    mkdir fastq
	fi
	mv ${SRA}_*.fastq fastq/
fi



#To-do
#add install checks
#add opts for directory outputs
#write readme
#make logs/verbose?
