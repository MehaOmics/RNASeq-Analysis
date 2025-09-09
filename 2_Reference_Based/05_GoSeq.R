library(goseq)
## Load Tables as data frame
DEG <- read.table("DE_Info.txt",header = TRUE)
LN <- read.table("Length_Info.txt",header = TRUE)
GO <- read.table("GO_Terms.txt",header = TRUE)

##Convert to named Vector
DEG.vector <- as.vector(DEG$DE_Status)
names(DEG.vector)<- DEG$Gene_id
LN.vector <- as.vector(LN$Length)
names(LN.vector)<- LN$Length


##Make empty list
GO.list <- list()

##To populate list from GO Dataframe
for (i in 1:nrow(GO)){
  GO.list[[i]]<-unlist(strsplit(GO[i,2], split = ","))
}

##Same names to GO ID
names(GO.list)<-GO$Gene_ID
View(GO.list)

## 
pwf=nullp(DEG.vector,bias.data=LN.vector)
GO.wall = goseq(pwf, gene2cat = GO.list, use_genes_without_cat=TRUE)
View(GO.wall)

#this gave table with p-values...now correct for multiple testing using FDR
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

# add new column with over represented GO terms padj
GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
head(enriched.GO)

# add new column with under represented GO terms padj
GO.wall$under_rep_padj=p.adjust(GO.wall$under_represented_pvalue, method="BH")
dep.GO=GO.wall$category[p.adjust(GO.wall$under_represented_pvalue, method="BH")<.05]
write.csv(GO.wall, file = "Upregulated_GOSEq.csv")




