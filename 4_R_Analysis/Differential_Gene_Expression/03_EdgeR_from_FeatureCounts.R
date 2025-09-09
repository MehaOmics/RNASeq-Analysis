library(edgeR)

# Load count data
x <- read.delim("Feature_count.txt", header = TRUE, sep = "\t", row.names = 1)

# Define groups
group <- factor(c(rep("control",5), rep("Inoculated",5)))
y <- DGEList(counts = x, group = group, genes = row.names(x))

# Filter lowly expressed genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
y <- calcNormFactors(y)

# Create design matrix
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# Fit GLM and perform QL test
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit)

# Top genes
topGenes <- topTags(qlf, n = Inf)
write.table(topGenes$table, file = "EdgeR_AllGenes.tsv", sep = "\t", quote = FALSE)

# Identify significant DEGs (FDR < 0.05 & |log2FC| ≥ 1.5)
de.HL <- decideTestsDGE(qlf, lfc = log2(1.5))
summary(de.HL)
write.csv(as.data.frame(de.HL), file = "EdgeR_DEGs_FDR0.05_LFC.csv")

# Optional: significant downregulated genes ordered by logFC
sigDownReg <- topGenes$table[topGenes$table$FDR < 0.05 & topGenes$table$logFC < -log2(1.5), ]
sigDownReg <- sigDownReg[order(sigDownReg$logFC), ]
write.table(sigDownReg, file = "EdgeR_SigDownregulated.tsv", sep = "\t", quote = FALSE)
