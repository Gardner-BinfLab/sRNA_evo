#!/bin/bash

#Takes maf alignments produced by sRNA_evo.sh, runs RNA QC programs and produces summary stats in QC.csv
for i in *.maf; do
    esl-reformat clustal $i | RNAcode > $i.RNAcode
done
for i in `ls *.maf | sed 's/.maf//'`; do
    weight -f 0.95 $i.maf > $i.maf.strip
    esl-reformat clustal $i.maf.strip | RNAalifold --id-prefix=$i --aln-stk=$i > $i.alifold
    esl-reformat clustal $i.maf.strip | RNAcode > $i.RNAcode
    R-scape $i.stk > $i.rscape
    ./runr2r.sh $i.stk
done


#Parse output files from various RNA QC programs. 
echo "name,ali_energy,ali_MFE,ali_cov,Rsc_cov_min,Rsc_conv_max,code_length,code_score,code_P"
for i in `ls *.maf | sed 's/.maf//'`; do
    echo -n "$i,"
    #RNAalifold
    #hey jude
    #if file exists
    if [[ -s $i.alifold ]]; then
	    #take that output\
		#and make it smaller
        echo -n `sed 's| (|\t|;s|)||g' $i.alifold | cut -f2 | tail -n1 | awk 'BEGIN{OFS=","} {print $1,$3,$5}'`
    else
	echo -n ",NA,NA,NA"
    fi
    #R-scape
    if [[ -s $i.out $(( `grep GTp $i.out | wc -l` )) -gt 0 ]]; then
	echo -n"," `grep GTp $i.out | tr -d '[]|#GTp' |  sed -E 's/^\s+//g;s/\s+$//g;s/\s+/,/g' | cut -d "," -f2-3`
    else
	echo -n ",NA,NA"
    fi
    #RNAcode
    if [[ -s "$i.maf.RNAcode" && $(( `grep HSS $i.RNAcode | wc -l` )) -gt 0 ]]; then
	    grep -P '^\s+\d+' $i.maf.RNAcode | head -n 1 | perl -lane 'print ",",$F[2],",",$F[6],",",$F[8]; '
    else
	echo ",NA,NA,NA"
    fi
done > QC.csv
