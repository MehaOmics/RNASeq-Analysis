library(DESeq2)
library(tximport)
library(readr)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
# ========================
#Input Files
# ========================
# Define Salmon quantification files
files <- c(
  "LB1" = "LB1_quant.sf", "LB2" = "LB2_quant.sf", "LB3" = "LB3_quant.sf",
  "EL1" = "EL1_quant.sf", "EL2" = "EL2_quant.sf", "EL3" = "EL3_quant.sf",
  "ML1" = "ML1_quant.sf", "ML2" = "ML2_quant.sf", "ML3" = "ML3_quant.sf",
  "S1" = "S1_quant.sf",   "S2" = "S2_quant.sf",   "S3" = "S3_quant.sf",
  "YR1" = "YR1_quant.sf", "YR2" = "YR2_quant.sf", "YR3" = "YR3_quant.sf"
)


# ========================
# Sample Metadata
# ========================
sample.info <- data.frame(
  SampleID = c("LB1", "LB2", "LB3", "EL1", "EL2", "EL3", "ML1", "ML2", "ML3",
               "S1", "S2", "S3","YR1","YR2","YR3"), 
  condition = factor(c("Condition1", "Condition1", "Condition1",
                       "Condition2", "Condition2", "Condition2", 
                       "Condition3", "Condition3", "Condition3",
                       "Condition4", "Condition4", "Condition4","Condition5", "Condition5", "Condition5")),
  Group = factor(c("Treatment", "Treatment", "Treatment",
                   "Treatment", "Treatment", "Treatment",
                   "Treatment", "Treatment", "Treatment",
                   "Treatment", "Treatment", "Treatment",
                   "Control", "Control", "Control" 
                     ))
)
rownames(sample.info) <- sample.info$SampleID

# ========================
# Import Quantification Data
# ========================
txi <- tximport(files, type="salmon", txOut=TRUE)

# ========================
#  DESeq2 Analysis
# ========================
dds <- DESeqDataSetFromTximport(txi, colData = sample.info, design = ~ Group)
dds <- DESeq(dds)

# Shrink Log2 Fold Changes
resLFC <- lfcShrink(dds, coef = "Group_Control_vs_Treatment", type = "apeglm")
resLFC_df <- as.data.frame(resLFC) %>%
  rownames_to_column("transcript_id")

# ========================
# Filter DEGs
# ========================
padj.cutoff <- 0.05
lfc.cutoff <- 1

resLFC_df <- resLFC_df %>%
  mutate(Significance = ifelse(!is.na(padj) & padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff, 1, 0))

# Save DEGs and all results
write.csv(resLFC_df, "./DEG-Output/Transcript_level/ControlvsTreatment/results.csv", row.names = FALSE)

# Significant DEGs only
sig_degs <- resLFC_df %>% filter(Significance == 1)
write.csv(sig_degs, "./DEG-Output/Transcript_level/ControlvsTreatment/Significant_DESeq2.csv", row.names = FALSE)

# ========================
# 📈 DEG Summary
# ========================
num_upregulated <- sum(sig_degs$log2FoldChange > lfc.cutoff, na.rm = TRUE)
num_downregulated <- sum(sig_degs$log2FoldChange < -lfc.cutoff, na.rm = TRUE)

cat("Upregulated DEGs:", num_upregulated, "\n")
cat("Downregulated DEGs:", num_downregulated, "\n")

# ========================
# 🔎 Normalized Counts
# ========================
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "./DEG-Output/Transcript_level/ControlvsTreatment/normalized_counts.csv")
