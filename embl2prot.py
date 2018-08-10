#!/usr/bin/python3
#embl2prot.py - Takes embl files and a .gff file of annotation
#gff seqname must be in name of embl file
#For each annotation, finds and prints the 6 nearest proteins (3 upstream, 3 downstream)
#Usage: embl2prot.py [path to gff]
#run in same directory as embl files
#todo: make proper opts, take arg for embl directory

import Bio.SeqIO
import os
import sys

path=os.getcwd()
files = os.listdir(path)

gff = sys.argv[1]

with open(gff, "r") as annotations:
    #split gff annotations by field
    for line in annotations:
        line = line.split()
        genome = line[0]
        gene = line[2]
        start = line[3]
        end = line[4]
        emblFile = ''
        after =0
        #find relevant embl file for annotation
        for embl in files:
            if genome in embl:
                if embl.endswith(".embl"):
                    emblFile = embl
                    break
        embl = path + "/" + emblFile
        seqList = []
        for record in Bio.SeqIO.parse(embl,"embl"):
            #open embl and store all proteins in list 'seqList'
            for feat in record.features:
                if feat.type=='CDS':
                    startCDS=feat.location.start
                    seq = feat.qualifiers.get("translation")
                    seqList.append(seq)
                    #Find first protein downstream of annotation & record list index in 'after'
                    if after== 0 and int(startCDS) >= int(end):
                        after = (len(seqList)-1)
            seqNum = 1 #name prots 1-6 (1-3 upstream, 4-6 downstream)
            #Stop out-of-range calls
            rangeMax = after+3
            rangeMin = after-3
            if rangeMax > (len(seqList)-1):
                rangeMax = len(seqList)-1
            if rangeMin < 0:
                rangeMin = 0
            #print protein sequence in fasta file
            #assumes 1 annotation of each gene per genome
            #naming format: >GENOME_GENE_SEQ[1-6]
            for i in range(rangeMin, rangeMax,1):
                header=">"+genome+"_"+gene+"_seq"+str(seqNum)
                print(header)
                num += 1
                print(str(seqList[i]).strip("[]'"))
        
    
