library(RColorBrewer)
library(ggplot2)
library(stringr)
library(reshape2)
library(sqldf)
library(viridis)
library(gtools)
library(DESeq2)


###########################################################################
##################### Load the data #######################################
###########################################################################


## Read the main data
ob_cond   <- read.csv("data/cond.csv")
rownames(ob_cond) <- ob_cond[,1]
ob_cond <- ob_cond[,-1]
ob_counts <- read.csv("data/counts.csv")
rownames(ob_counts) <- ob_counts[,1]
ob_counts <- ob_counts[,-1]

ensid <- read.csv("genesym_human.csv", stringsAsFactors = FALSE,sep="\t")
colnames(ensid) <- c("geneid","symbol")

ensid_covid <- read.csv("genesym_covid.csv", stringsAsFactors = FALSE,sep="\t")
colnames(ensid_covid) <- c("geneid","symbol")

ensid <- rbind(ensid,ensid_covid)


## reduced ensid list for our count table
ensid_red <- ensid[ensid$geneid %in% rownames(ob_counts),]
ensid_red <- ensid_red[!duplicated(ensid_red$geneid),]
rownames(ensid_red) <- ensid_red$geneid



###################################################################################################################
##################### DE testing ##################################################################################
###################################################################################################################



####################################################   Impact of infection, on untreated cells --------- fig 4c

keep <- ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]

ggplot(r[str_starts(r$geneid,"ENSG"),], aes(log2FoldChange, -log10(pvalue))) + geom_point(size=0.1) 


####################################################   impact of betamethasone, on infected cells ------- fig 4d

keep <- ob_cond$infected=="TRUE"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~Treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("Treatment","beta03","none")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]

ggplot(r[str_starts(r$geneid,"ENSG"),], aes(log2FoldChange, -log10(pvalue))) + geom_point(size=0.1)


