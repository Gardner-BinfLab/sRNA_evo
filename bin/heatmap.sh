#!/bin/bash

#produces summary stats for heatmap figure, to show sequence conservation and annotation success for each genera
#
usage(){
    echo "
heatmap.sh: generate summary stats for heatmap figure, to show sequence conservation and annotation success for each genus

usage: heatmap.sh [working directory containing files from sRNA_evo.sh and sRNA_QC.sh]
    
    Options:
		-h	Display this help
"
}
while getopts :h opt; do
    case "${opt}" in
	h) usage;exit;;
	\?) echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done
shift $((${OPTIND}-1))

if [ -z "$1" ]
  then
      echo "ERROR! No argument supplied for directory"
      usage
      exit 1
fi

DIR=$1

if [ -z ${DIR}/data/FILT_trainingSeqs.gff ] || [ -z ${DIR}/data/FILT_startingSeqs.gff ]
    then
      echo "ERROR! Can't find annotation files in ${DIR}/data/*.gff"
      usage
      exit 1
fi

##make opt?
if [ -z ~/data/summary.txt ]
    then
      echo "ERROR! Can't find summary.txt in ~/data/summary.txt (a bad hardcoded path)"
      usage
      exit 1
fi


    
##From annotations, grab genome accesions that are from Enterobacteriaceae
cut -f1 ${DIR}/data/FILT_trainingSeqs.gff ${DIR}/data/FILT_startingSeqs.gff | sort | uniq | xargs -ifoo grep foo ~/data/summary.txt | grep Enterobact* | awk -F"\t|;." '{print $1" "$(NF-1)}' | cut -d ' ' -f1-2 > ${DIR}/data/Ents.txt

cut -f3 ${DIR}/data/*Seqs.gff | sort | uniq > ${DIR}/data/sRNA_list.txt

getSeqs(){
    echo -n "Indexing ${DIR}/QC files..."
    ##index fasta files generated by sRNA_QC
    ls ${DIR}/QC/*.fasta | parallel --jobs 20 'esl-sfetch --index {} >/dev/null'
    mkdir ${DIR}/data/heatmap
    #if 3 or more results for a genus, make a file of genome accessions for that genus
    for i in `cut -d ' ' -f2 ${DIR}/data/Ents.txt | sort | uniq -d -c | awk '{if($1>2){print $2}}'`; do
	grep -E "$i$" ${DIR}/data/Ents.txt | cut -d " " -f1 | sed 's/\..*$//g' > ${DIR}/data/heatmap/$i.genus
    done
    echo "Done!"
}


make_heatmap(){    
    ##For each genus:
    ##for each genome/gene, grab sequences, align each sequence to the original seed sequence 
    ##create csv file of % identities for each genus
    
    echo -n "Getting sequence identities..."
    for i in ${DIR}/data/heatmap/*.genus; do
	#make column headers (genome,gene names)
	echo -n "," >> $i.identity
	cat ${DIR}/data/sRNA_list.txt | xargs echo | sed "s/ /,/g" >> $i.identity
	##for each genome/gene, grab sequences, align each sequence to the original seed sequence 
	while read genome; do
	    #rowname (genome)
	    echo -n $genome >> $i.identity
	    while read gene; do
                esl-sfetch ${DIR}/QC/$gene.fasta $genome\_$gene 1> $i.fasta 2>/dev/null
		#if annotation exists for this gene, align & add % identity as value. else record "NA".
                if [[ -s $i.fasta ]]; then
		    cat ${DIR}/hmms/starting/$gene.fasta >> $i.fasta
                    hmmalign ${DIR}/hmms/starting/$gene.hmm $i.fasta > $i.stk
                    echo -n "," >> $i.identity && echo -n `esl-alistat $i.stk | grep identity | cut -d " " -f6 | tr -d %` >> $i.identity
		else
                    echo -n ",NA" >> $i.identity
                fi
	    done < ${DIR}/data/sRNA_list.txt 
            echo "" >> $i.identity
	done < $i && rm $i.stk $i.fasta
    done
    echo "Done!"
    ##For each genus:
    ##for each genome/gene, count annotations 
    ##create csv file of annotation success for all genera
    echo -n "Getting annotation counts..."
    #set up colnames (sRNAs, Total_Annotations, Missing_in)
    echo -n "sRNA," >> ${DIR}/data/heatmap/all.csv
    while read line; do
	echo -n "$line," >> ${DIR}/data/heatmap/all2.csv
    done < ${DIR}/data/sRNA_list.txt
    echo "Total_annotations,Missing_in" >> ${DIR}/data/heatmap/all.csv

    #get rownames (genera)
    ls ${DIR}/data/heatmap/*.genus | sed 's/.genus$//g' > ${DIR}/data/heatmap/list.txt

    #for each genus and gene, count number of annotations. Also track total annotations and total genomes with no annotations (missing_in)
    while read gene; do
	annotationTotal=0
	genusSum=0
	echo -n "$gene," 
	while read genus; do
	    genusTotal=0
            while read genome; do
		genusSum=$(( $genusSum + 1 ))
		if grep -q $genome ${DIR}/QC/$gene.fasta ; then
		    genusTotal=$(( $genusTotal + 1 ))
		fi
	    done < $genus.genus
	    annotationTotal=$(( $annotationTotal + $genusTotal ))
	    echo -n "$genusTotal,"
	done < ${DIR}/data/heatmap/list.txt
	echo "$annotationTotal,$(( $genusSum - $annotationTotal ))"
    done < ${DIR}/data/sRNA_list.txt >> ${DIR}/data/heatmap/all.csv
    echo "Done!"
}

cleanup(){
    echo -n "Cleaning up..."
    rm ${DIR}/data/sRNA_list.txt ${DIR}/data/heatmap/Ents.txt ${DIR}/heatmap/*.genus ${DIR}/data/heatmap/list.txt
    echo "Done!"
}

getSeqs
make_heatmap
cleanup
