library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# ========================
#Volcano Plot (from DeSeq2)
# ========================
resLFC_df <- resLFC_df %>%
  mutate(significance = case_when(
    padj < padj.cutoff & log2FoldChange > lfc.cutoff ~ "Upregulated",
    padj < padj.cutoff & log2FoldChange < -lfc.cutoff ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))
png("Volcano_plot.png", width = 4000, height = 3000, res = 600)

volcano <- ggplot(resLFC_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "orange", "Downregulated" = "#AEC6CF", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed") +
  labs(
    title = "Volcano Plot (Below vs Above Ground Tissues)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  theme_classic()
# ========================
#  Heatmap of Significant Genes (from DeSeq2)
# ========================
sig_genes <- sig_degs$transcript_id
heat_data <- normalized_counts[sig_genes, ]

# Optional: pre-scale for speed
heat_data_scaled <- t(scale(t(heat_data)))

# Save PNG at high resolution
png("./Heatmapviridis.png",
    width = 5000, height = 3500, res = 600)

heatmap <- pheatmap(heat_data_scaled,
                    col = viridis(100),        # use viridis color palette
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,      # no clustering for tissue types
                    show_rownames = FALSE,
                    border_color = NA,
                    fontsize = 10,
                    height = 20)

print(heatmap)
dev.off()
