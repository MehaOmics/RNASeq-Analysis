library(DESeq2) 
library(ggplot2)
library(tidyverse)
library(reshape)
library(RColorBrewer)
library(pheatmap)

##Creating countData
read.counts <- read.table("Feature_Counts_Sorted.txt",header = TRUE)
## gene IDs should be stored as row.names
row.names(read.counts)<-read.counts$Geneid

##Remove first column
read.counts<- read.counts[,-c(1)]
view(read.counts)

names (read.counts) <- c ("Sample01","Sample02","Sample03","Sample04" ,"Sample05", "Sample06","Sample07","Sample08","Sample09" ,"Sample10")
##give meaningful sample names
names(read.counts) <- c( paste("Control", c(1:5), sep = "_"),
                        paste("Inoculated", c(1:5), sep = "_") )
###data frame str converts values to string form. To check data
str(read.counts)
head(read.counts)
##Creating colData:this should be a data.frame, where the rows directly match the columns of the count data
sample.info <- data.frame (condition = c(rep("Control",5),rep("Inoculated",5)),row.names = names(read.counts))
str(sample.info)
view(sample.info)

## Generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData=read.counts,colData = sample.info,design = ~ condition )
DESeq.ds
head(counts(DESeq.ds))
##To know number of reads generated in each library
colSums(counts(DESeq.ds))
colSums(counts(DESeq.ds)) |> barplot()
dim(DESeq.ds)
##Filteration:Remove genes with no reads. Here we perform a minimal pre-filtering to keep only rows that have at least 0 reads total?
keep_genes <- rowSums(counts(DESeq.ds)) > 0
dim(DESeq.ds)
DESeq.ds <- DESeq.ds[ keep_genes, ]
dim(DESeq.ds)
counts(DESeq.ds) |> str()
assay(DESeq.ds) |> str()

##Diffrential Gene Expression
##set control condition as reference
dds <- relevel (DESeq.ds$condition, ref="Control")
dds <- DESeq(DESeq.ds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
res
resultsNames(dds)
##Plots to check data
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup = "condition")
plotDispEsts(dds)
##Two ways to proceed from here either with LFC Shrinkage or without LFC shrinkage:
## 1). Without shrinkage, Adjusted P value
#resOrdered <- res[order(res$pvalue),]
#padj.cutoff <- 0.05
#lfc.cutoff <- 0.58
#threshold <- resOrdered$padj < padj.cutoff & abs(resOrdered$log2FoldChange) > lfc.cutoff
#length(which(threshold))
#summary(resOrdered, threshold)

##To add this vector to our results table we can use the $ notation to create the column on the left hand side of the assignment operator, and the assign the vector to it instead of using cbind()

#resOrdered$threshold <- threshold
##Now we can easily subset the results table to only include those that are significant using the subset() function
#sigOE <- data.frame(subset(resOrdered, threshold==TRUE))
#sigOE
#write.csv(as.data.frame(sigOE), file="DESeq2_Stranded_Results.csv")

## With Shrinkage(Used in current analysis)
## Two arguments to the results function allow for threshold-based Wald tests: lfcThreshold, which takes a numeric of a non-negative threshold value, and altHypothesis, which specifies the kind of test
## Log fold change shrinkage for visualization and ranking 
library(apeglm)
resLFC <- lfcShrink(dds, coef="condition_Inoculated_vs_Control", type="apeglm")
## In case you want to filter results. But in final analysis i did filtration on final output while making plots
resLFC
resOrdered <- resLFC[order(resLFC$pvalue),]
resOrdered
sum(resOrdered$padj < 0.05, na.rm=TRUE)
## To know summary of DEGs at FRD 0.05, summary() function doesn't have an argument for fold change threshold.
summary(results(dds, alpha=0.05))
## For all the genes make table from resLFC. This will make table for all the genes without any filteration
write.csv(as.data.frame(resLFC), file="DESeq2_Results_with_Shrinkage_2.csv")

##Heatmap: resLFC and res_tableOE are same
resLFC <- lfcShrink(dds, coef="condition_Inoculated_vs_Control", type="apeglm")
##Let's first create variables that contain our threshold criteria:
padj.cutoff <- 0.05
lfc.cutoff <- 1
##Let's create vector that helps us identify the genes that meet our criteria:
threshold <- resLFC$padj < padj.cutoff & abs(resLFC$log2FoldChange) > lfc.cutoff
length(which(threshold))

##To add this vector to our results table we can use the $ notation to create the column on the left hand side of the assignment operator, and the assign the vector to it instead of using cbind()

resLFC$threshold <- threshold
##Now we can easily subset the results table to only include those that are significant using the subset() function
sigOE <- data.frame(subset(resLFC, threshold==TRUE))
sigOE
write.csv(as.data.frame(sigOE), file="DESeq2_Filtered_Output.csv")


##Order significant results by padj values
sigOE_ordered <- sigOE[order(sigOE$padj), ]
##For top 40 genes
sigOE_orderedTOP <- as.data.frame(sigOE_ordered)[1:40, ]
write.csv(as.data.frame(sigOE_orderedTOP), file="DESeq2_LFCShrink_Top40.csv")

##we can extract the normalized count values:
normalized_counts <- counts(dds, normalized=T)


