setwd("/home/ibab/Bioinformatics/NGS/Term_Paper_RA/Phase2/final_counts_matrices/")
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

early_ra_F = read.csv("ra_early_f_coding_final.csv", row.names = 1)
# length(early_ra_F$X)
# length(unique(early_ra_F$X))
# rownames(early_ra_F) <- early_ra_F$X
# early_ra_F<- early_ra_F[, -1] 
colnames(early_ra_F) <- paste("early_", colnames(early_ra_F), sep = "")

early_ra_M = read.csv("ra_early_m_coding_final.csv", row.names = 1)
# rownames(early_ra_F) <- early_ra_F$X
# early_ra_F<- early_ra_F[, -1] 
colnames(early_ra_M) <- paste("early_", colnames(early_ra_M), sep = "")
RA_F = read.csv("ra_estd_f_coding_final.csv", row.names = 1)
colnames(RA_F) <- paste("estd_", colnames(RA_F), sep = "")
RA_M = read.csv("ra_estd_m_coding_final.csv", row.names = 1)
colnames(RA_M) <- paste("estd_", colnames(RA_M), sep = "")
control_F = read.csv("synovium_controls_f_coding_final.csv", row.names = 1)
colnames(control_F) <- paste("control_", colnames(control_F), sep = "")
control_M = read.csv("synovium_controls_m_coding_final.csv", row.names = 1)
colnames(control_M) <- paste("control_", colnames(control_M), sep = "")

control_v_early_metadata_m = data.frame(sample = c(colnames(control_M)[2:ncol(control_M)], colnames(early_ra_M)[2:ncol(early_ra_M)]) , 
                                        condition = c(rep("control", length(colnames(control_M)) -1 ),
                                                      rep("early_M", length(colnames(early_ra_M)) -1 )))
control_v_estd_metadata_m = data.frame(sample = c(colnames(control_M)[2:ncol(control_M)], colnames(RA_M)[2:ncol(RA_M)]), 
                                       condition =c(rep("control", length(colnames(control_M))-1),
                                                    rep("RA_M", length(colnames(RA_M))-1)))
control_v_early_metadata_f = data.frame(sample = c(colnames(control_F)[2:ncol(control_F)], colnames(early_ra_F)[2:ncol(early_ra_F)]), 
                                        condition = c(rep("control", ncol(control_F) -1 ),
                                                      rep("early_F", ncol(early_ra_F)-1 )))
control_v_estd_metadata_f = data.frame(sample = c(colnames(control_F)[2:ncol(control_F)], colnames(RA_F)[2:ncol(RA_F)]), 
                                       condition =c(rep("control", length(colnames(control_F))-1),
                                                    rep("RA_F", length(colnames(RA_F))-1)))

control_v_early_f = merge(control_F, early_ra_F, by.x='control_X', by.y ='early_X')
control_v_estd_m = merge(control_M, RA_M, by.x = 'control_X', by.y = 'estd_X')
control_v_early_m = merge(control_M, early_ra_M, by.x='control_X', by.y ='early_X')
control_v_estd_f = merge(control_F, RA_F, by.x = 'control_X', by.y = 'estd_X')


# control_v_early_f[, -1] <- lapply(control_v_early_f[,-1], as.integer)
control_v_early_f <- control_v_early_f[!duplicated(control_v_early_f$control_X),]
rownames(control_v_early_f) <- control_v_early_f$control_X
control_v_early_f <- control_v_early_f[, -1]

# control_v_early_m[, -1] <- lapply(control_v_early_m[,-1], as.integer)
control_v_early_m <- control_v_early_m[!duplicated(control_v_early_m$control_X),]
rownames(control_v_early_m) <- control_v_early_m$control_X
control_v_early_m <- control_v_early_m[, -1]

# control_v_estd_f[, -1] <- lapply(control_v_estd_f[,-1], as.integer)
control_v_estd_f <- control_v_estd_f[!duplicated(control_v_estd_f$control_X),]
rownames(control_v_estd_f) <- control_v_estd_f$control_X
control_v_estd_f <- control_v_estd_f[, -1]

# control_v_estd_m[, -1] <- lapply(control_v_estd_m[,-1], as.integer)
control_v_estd_m <- control_v_estd_m[!duplicated(control_v_estd_m$control_X),]
rownames(control_v_estd_m) <- control_v_estd_m$control_X
control_v_estd_m <- control_v_estd_m[, -1]


# at_ce <- at_ce[, !colnames(at_ce) %in% c("ce2", "ce4", "at6", "ce6", "ce5", "ce7", "at5", "at1")] # df of counts matrix
# at_vs_ce_metadata = read.csv("at_vs_ce_rem2.csv", row.names = 1) # metadata file (attached). Make sure the colnames of df is same as that in Metadat
#print(ncol(control_v_early_f))
#nrow(control_v_early_metadata_f)
control_v_early_f[] <- lapply(control_v_early_f, as.integer)
control_v_early_m[] <- lapply(control_v_early_m, as.integer)
control_v_estd_f[] <- lapply(control_v_estd_f, as.integer)
control_v_estd_m[] <- lapply(control_v_estd_m, as.integer)

control_v_early_metadata_f$condition <- factor(control_v_early_metadata_f$condition)
control_v_early_metadata_m$condition <- factor(control_v_early_metadata_m$condition)
control_v_estd_metadata_f$condition <- factor(control_v_estd_metadata_f$condition)
control_v_estd_metadata_m$condition <- factor(control_v_estd_metadata_m$condition)

dds_cef = DESeqDataSetFromMatrix(countData = control_v_early_f,
                             colData = control_v_early_metadata_f,
                             design = ~ condition)

dds_cem = DESeqDataSetFromMatrix(countData = control_v_early_m,
                                 colData = control_v_early_metadata_m,
                                 design = ~ condition)

dds_craf = DESeqDataSetFromMatrix(countData = control_v_estd_f,
                                 colData = control_v_estd_metadata_f,
                                 design = ~ condition)

dds_cram = DESeqDataSetFromMatrix(countData = control_v_estd_m,
                                 colData = control_v_estd_metadata_m,
                                 design = ~ condition)

gdds_cef = dds_cef[rowSums(counts(dds_cef)) > 10, ]
gdds_cem = dds_cem[rowSums(counts(dds_cem)) > 10, ]
gdds_craf = dds_craf[rowSums(counts(dds_craf)) > 10, ]
gdds_cram = dds_cram[rowSums(counts(dds_cram)) > 10, ]

# nrow(dds)
vsd_cef <- vst(dds_cef, blind = TRUE)# Variance stabilizing transformation
vsd_cem <- vst(dds_cem, blind = TRUE)
vsd_craf <- vst(dds_craf, blind = TRUE)
vsd_cram <- vst(dds_cram, blind = TRUE)

pca_data_cef = plotPCA(vsd_cef, intgroup = "condition", returnData = TRUE)
pca_data_cem = plotPCA(vsd_cem, intgroup = "condition", returnData = TRUE)
pca_data_craf = plotPCA(vsd_craf, intgroup = "condition", returnData = TRUE)
pca_data_cram = plotPCA(vsd_cram, intgroup = "condition", returnData = TRUE)

# Extract percentage variance explained
percentVar_cef = round(100 * attr(pca_data_cef, "percentVar"))
percentVar_cem = round(100 * attr(pca_data_cem, "percentVar"))
percentVar_craf = round(100 * attr(pca_data_craf, "percentVar"))
percentVar_cram = round(100 * attr(pca_data_cram, "percentVar"))


# Create PCA plot with sample names and variance info
ggplot(pca_data_cef, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_data_cef))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar_cef[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar_cef[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers Control vs Early in Females") +
  theme_minimal()

ggplot(pca_data_cem, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_data_cem))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar_cem[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar_cem[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers Control vs Early in Males") +
  theme_minimal()

ggplot(pca_data_craf, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_data_craf))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar_craf[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar_craf[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers Control vs Established in Females") +
  theme_minimal()

ggplot(pca_data_cram, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_data_cram))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar_cram[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar_cram[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers Control vs Established in Males") +
  theme_minimal()


dds_cef = DESeq(dds_cef)  
dds_cem = DESeq(dds_cem)  
dds_craf = DESeq(dds_craf)  
dds_cram = DESeq(dds_cram)  

resultsNames(dds_cef)
dds_vec = list(dds_cef, dds_cem, dds_craf, dds_cram)
normalized_counts = list(control_v_early_f, control_v_early_m, control_v_estd_f, control_v_estd_m)
i = 0
for (dds in dds_vec){
  i = i + 1
#dds = dds_cef
  split_string <- strsplit(resultsNames(dds), "_")[[2]]
  part1 <- split_string[1]    # "condition"
  part2 <- paste(split_string[2],split_string[3],sep="_")    
  part3 <- split_string[5] 
  res = results(dds, contrast = c(part1, part3, part2))
  res = res[order(res$padj, na.last = NA), ]
  res_df = as.data.frame(res)
  res_p_val_cut = subset(res_df, res_df$padj <= 0.05)
  at_vs_ce_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
  at_vs_ce_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)
  
  at_vs_ce_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$padj <= 0.01)
  at_vs_ce_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$padj <= 0.01)
  
  at_vs_ce_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$padj <= 0.01)
  at_vs_ce_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$padj <= 0.01)
  print(at_vs_ce_downgenes_0_5)
  
  write.csv(at_vs_ce_upgenes_1_5, file = paste(resultsNames(dds)[2],"1.5_up.csv"))
  write.csv(at_vs_ce_downgenes_1_5, file = paste(resultsNames(dds)[2],"1.5_down.csv"), row.names = T)
  
  write.csv(at_vs_ce_upgenes_1, file = paste(resultsNames(dds)[2],"1_up.csv"), row.names = T)
  write.csv(at_vs_ce_downgenes_1, file = paste(resultsNames(dds)[2],"1_down.csv"), row.names = T)
  
  write.csv(at_vs_ce_upgenes_0_5, file = paste(resultsNames(dds)[2],"0.5_up.csv"), row.names = T)
  write.csv(at_vs_ce_downgenes_0_5, file = paste(resultsNames(dds)[2],"0.5_down.csv"), row.names = T)
  print("__________")
  # Get top 25 upregulated genes
  top_25_up_ctrl_mer <- at_vs_ce_upgenes_1_5 %>%
    filter(log2FoldChange >= 1.5 & pvalue <= 0.05) %>%
    arrange(desc(log2FoldChange)) %>%
    head(25)
  print("__________")
  # Get top 25 downregulated genes
  top_25_down_ctrl_mer <- at_vs_ce_downgenes_1_5 %>%
    filter(log2FoldChange <= -1.5 & pvalue <= 0.05) %>%
    arrange(log2FoldChange) %>%
    head(25)
  print("__________")
  normalized_counts_df_up = subset(normalized_counts[[i]], rownames(normalized_counts[[i]]) %in% 
                                     rownames(top_25_up_ctrl_mer))
  print(normalized_counts_df_up)
  print("__________")
  normalized_counts_df_down = subset(normalized_counts[[i]], rownames(normalized_counts[[i]]) %in% 
                                       rownames(top_25_down_ctrl_mer))
  print("__________")
  heatmap_data_up <- as.matrix(normalized_counts_df_up)
  heatmap_data_down <- as.matrix(normalized_counts_df_down)
  print("__________")
  print(heatmap_data_up)
  pheatmap(heatmap_data_up, 
           color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
           cluster_rows = T,  # Cluster rows
           cluster_cols = T,  # Cluster columns
           scale = "row",        # Normalize by row (optional)
           show_rownames = TRUE,  # Show row labels
           show_colnames = TRUE,  # Show column labels
           main = paste("Heatmap of Top 25 Upregulated Genes",resultsNames(dds)[2]))
  pheatmap(heatmap_data_down, 
           color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
           cluster_rows = T,  # Cluster rows
           cluster_cols = T,  # Cluster columns
           scale = "row",        # Normalize by row (optional)
           show_rownames = TRUE,  # Show row labels
           show_colnames = TRUE,  # Show column labels
           main = paste("Heatmap of Top 25 Downregulated Genes",resultsNames(dds)[2]))
  
  res_p_val_cut$neg_logp_val = -log10(res_p_val_cut$pvalue)
  res_p_val_cut$Expression = ifelse(res_p_val_cut$log2FoldChange >= 1.0 & res_p_val_cut$pvalue <= 0.05, "Upregulated",
                                   ifelse(res_p_val_cut$log2FoldChange <= -1.0 & 
                                          res_p_val_cut$pvalue <= 0.05, "Downregulated", "Not Significant"))
  
  volcano_plot = ggplot(res_p_val_cut, aes(x=log2FoldChange, y=neg_logp_val, colour = Expression)) +
    geom_point() +
    scale_color_manual(values=c("red", "black", "blue")) +
    geom_vline(xintercept = c(-1, 1), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
    labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)") +
    ggtitle(paste("Differential Expression Analysis ",resultsNames(dds)[2]))
  
  print(volcano_plot) 
}

# Plot heatmaps with top 25 upreg genes and top 25 downreg genes
# Plot Volcano Plot

# Get top 25 upregulated genes
top_25_up_ctrl_mer <- at_vs_ce_upgenes_1_5 %>%
  filter(logFC >= 1.5 & PValue <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ctrl_mer <- ctrl_mer_res %>%
  filter(logFC <= -1.5 & PValue <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts, rownames(normalized_counts) %in% 
                                   rownames(top_25_up_ctrl_mer))

normalized_counts_df_down = subset(normalized_counts, rownames(normalized_counts) %in% 
                                     rownames(top_25_down_ctrl_mer))

# Convert data frame to matrix (if not already in matrix format)
heatmap_data_up <- as.matrix(normalized_counts_df_up)

# Create the heatmap
pheatmap(heatmap_data_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                        "C10", "C11", "C12", "C13", "C14", "Mer1", "Mer2", "Mer3", 
                        "Mer4", "Mer5", "Mer6", "Mer7", "Mer8", "Mer9", "Mer10", 
                        "Mer11", "Mer12", "Mer13", "Mer14", "Mer15", "Mer16", "Mer17", 
                        "Mer18", "Mer19", "Mer20", "Mer21", "Mer22", "Mer23", "Mer24"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Male Early RA")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                        "C10", "C11", "C12", "C13", "C14", "Mer1", "Mer2", "Mer3", 
                        "Mer4", "Mer5", "Mer6", "Mer7", "Mer8", "Mer9", "Mer10", 
                        "Mer11", "Mer12", "Mer13", "Mer14", "Mer15", "Mer16", "Mer17", 
                        "Mer18", "Mer19", "Mer20", "Mer21", "Mer22", "Mer23", "Mer24"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Male Early RA")


ctrl_mer_res$neg_logp_val = -log10(ctrl_mer_res$PValue)
ctrl_mer_res$Expression = ifelse(ctrl_mer_res$logFC >= 1.0 & ctrl_mer_res$PValue <= 0.05, "Upregulated",
                                 ifelse(ctrl_mer_res$logFC <= -0.5 & 
                                          ctrl_mer_res$PValue <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(ctrl_mer_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)") +
  ggtitle("Differential Expression Analysis of Cardioembolic vs Embolic Stroke of Undetermined Source")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Male Early RA")

write.csv(at_vs_ce_downgenes_0_5, file = "at_vs_rest_downgenes_0_5.csv", row.names = T)