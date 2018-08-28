rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggthemes)
library(viridis)
library(plotrix)
library(grid)

setwd("/mnt/9386241f-599d-4e2c-b4fc-60498ed920a1/projects/salmonella_srna/data/heatmap")

#all_2.csv - all.csv but columns ordered manually to reflect phylogeny
all <- read.csv("all_2.csv", header = TRUE)

percIdentity <- data.frame(stringsAsFactors = FALSE)
numAnnotations <- data.frame(data.frame(stringsAsFactors = FALSE))

#Make empty dataframe to store sum of average conservation to rank plot
avgsum <- data.frame(matrix(0,nrow=280, ncol=1))

#Get genera names
#genera ordering is already in correct order 
genera <- colnames(all)
for(x in c(1, 16, 16)){
  genera <- genera[-x]
}

for(x in genera){
  #For each genus in "all", open relevant file
  myFile <- paste(x,".genus.identity", sep='')
  myFile <- read.csv(myFile, header=TRUE, row.names=1)
  
  #percIdentity: each genus is row, each column is average seq conservation for that genus
  percIdentity <- rbind(colMeans(myFile, na.rm = TRUE), percIdentity)
  #numAnnotations: each genus is row, each column is average conservation (ie % of genomes in genus with annotation)
  numAnnotations <-rbind(all[,x]/max(all[,x], na.rm=TRUE), numAnnotations)
  #Avgsum = sum of annotations for each gene
  avgsum <- avgsum + (all[,x]/max(all[,x], na.rm=TRUE))
}


#Sort out column and row names
colnames(percIdentity) <- all$sRNA
colnames(numAnnotations) <- all$sRNA
rownames(percIdentity) <- rev(genera)
rownames(numAnnotations) <- rev(genera)
colnames(avgsum) <- "average"

#Subset of genes that have passed QC
srna_names <- read.table(file = "intergenic_only/intergenic2.txt")
tmpnames <- c()
for(x in srna_names$V1){
  tmpnames <- c(tmpnames,x)
}
srna_names <-tmpnames

#rankedGenes = give sRNAs values based on num of annotations * sum of conservation for ranking
rankedGenes <-all$Total_annotations*avgsum$average
rankedGenes <-cbind.data.frame(rankedGenes/max(rankedGenes))
rownames(rankedGenes) <- all$sRNA
rankedGenes <-rankedGenes[srna_names,]

#Transpose numAnnotations & percIdentity and rank sRNAs by rankedGenes values
percIdentity <-t(percIdentity)
percIdentity <- percIdentity[srna_names,]
percIdentity <-percIdentity[order(rankedGenes),]

percIdentity.m <- melt(percIdentity)
numAnnotations <- t(numAnnotations)
numAnnotations <- numAnnotations[srna_names,]
numAnnotations <- numAnnotations[order(rankedGenes),]
numAnnotations.m <- melt(numAnnotations)

#Plot as heatmap:
#1st layer - presence/absence of annotations from sRNA, genus (all2.m$Var2, all2.m$Var1)
#2nd layer - fill with colour as % sequence identity (identity.m), set alpha to % of genomes annotated for that genus (all2.m$value)

#to-do: fix legends
plot <- (ggplot(numAnnotations.m, aes(x=numAnnotations.m$Var2, y=numAnnotations.m$Var1))
          + geom_tile(aes(fill=percIdentity.m$value, alpha = numAnnotations.m$value), size= 0.3)
            + scale_alpha_continuous(range = c(0,1.0), name= "sRNA presence in clade (%)", guide=guide_legend(keywidth = 0.4, keyheight=0.4,legend.position = "right", title.position = "left", label.theme = element_text(size=1, angle =90), title.theme = element_text(size=1.2,angle=90)))
            + scale_fill_gradient(high= "Blue", low = "Red",name = "Average % sequence identity\n(relative to starting sequence)",guide=guide_colourbar(barwidth = 0.3, barheight=1.7, legend.position="right", title.position="left",label.theme = element_text(size=0.7, angle =90), title.theme = element_text(size=1, angle=90), title.hjust = 0.01))
            + coord_equal(ratio =0.5) + labs(y="sRNA", x="")
            + theme_tufte(base_family = "Helvetica")
            + theme(axis.ticks.x = element_blank(), axis.ticks.y= element_line(size = 0.1))
            + theme(axis.text.y= element_text(size=1, angle=20))
            + theme(axis.text.x= element_text(size = 3,angle = 90))
            + theme(legend.position = c(1,0.3))
            + theme(legend.margin = unit(0.25, "cm"))
            )
print(plot, vp=viewport(angle=-90))