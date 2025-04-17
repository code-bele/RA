# Limma-Voom
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Male Early RA vs Male Controls

male_early_ra = read.csv(file = "ra_early_m_coding_final.csv", row.names = 1)
length(unique(male_early_ra$X))

# Remove Duplicates
male_early_ra = male_early_ra[!duplicated(male_early_ra$X), ]
rownames(male_early_ra) = male_early_ra$X
male_early_ra$X = NULL

ctrl_m = read.csv(file = "synovium_controls_m_coding_final.csv", row.names = 1)

# Remove Duplicates
ctrl_m= ctrl_m[!duplicated(ctrl_m$X), ]
rownames(ctrl_m) = ctrl_m$X
ctrl_m$X = NULL

# combine the 2 dfs
ra_early_m_vs_ctrl_m = cbind(ctrl_m, male_early_ra)
colnames(ra_early_m_vs_ctrl_m)[1:14] = paste0("C", 1:14)
colnames(ra_early_m_vs_ctrl_m)[15:38] = paste0("Mer", 1:24)

# Run Limma:Voom
d0 =  DGEList(ra_early_m_vs_ctrl_m)
cutoff = 1
drop = which(apply(cpm(d0), 1, max) < cutoff)
d = d0[-drop,] 
dim(d) # number of genes left
snames = colnames(d) # Sample names
snames
group = as.factor(c(rep("C", 14), rep("Mer", 24)))
head(group)
plotMDS(d, col = as.numeric(group))
mm = model.matrix(~0 + group)
y = voom(d, mm, plot = T)
normalized_counts = y$E

fit = lmFit(y, mm)
head(coef(fit))
contr = makeContrasts(groupC - groupMer, levels = colnames(coef(fit)))
tmp = contrasts.fit(fit, contr)
head(tmp)
######Empirical Bayes smoothing of standard errors 
tmp = eBayes(tmp)
head(tmp)
ctrl_mer_res = topTable(tmp, sort.by = "P", n = Inf)
ctrl_mer_significant = subset(ctrl_mer_res, ctrl_mer_res$adj.P.Val <= 0.05)
write.csv(ctrl_mer_significant, file = "male_control_vs_male_early_ra_limma_allgenes.csv")
ctrl_mer_up_1_5 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC >= 1.5)
ctrl_mer_down_1_5 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC <= -1.5)
ctrl_mer_up_1 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC >= 1)
ctrl_mer_down_1 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC <= -1)
ctrl_mer_up_0_5 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC >= 0.5)
ctrl_mer_down_0_5 = subset(ctrl_mer_significant, ctrl_mer_significant$logFC <= -0.5)

# Get top 25 upregulated genes
top_25_up_ctrl_mer <- ctrl_mer_res %>%
  filter(logFC >= 1.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ctrl_mer <- ctrl_mer_res %>%
  filter(logFC <= -1.5 & adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts, rownames(normalized_counts) %in% 
                                   rownames(top_25_up_ctrl_mer))

normalized_counts_df_down = subset(normalized_counts, rownames(normalized_counts) %in% 
                                     rownames(top_25_down_ctrl_mer))
library(pheatmap)

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


ctrl_mer_res$neg_logp_val = -log10(ctrl_mer_res$adj.P.Val)
ctrl_mer_res$Expression = ifelse(ctrl_mer_res$logFC >= 0.5 & ctrl_mer_res$adj.P.Val <= 0.05, "Upregulated",
                                 ifelse(ctrl_mer_res$logFC <= -0.5 & 
                                          ctrl_mer_res$adj.P.Val <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(ctrl_mer_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Male Early RA")












# Female Early RA vs Female Controls

female_early_ra = read.csv(file = "ra_early_f_coding_final.csv", row.names = 1)
length(unique(female_early_ra$X))

# Remove Duplicates
female_early_ra = female_early_ra[!duplicated(female_early_ra$X), ]
rownames(female_early_ra) = female_early_ra$X
female_early_ra$X = NULL

ctrl_f = read.csv(file = "synovium_controls_f_coding_final.csv", row.names = 1)

# Remove Duplicates
ctrl_f= ctrl_f[!duplicated(ctrl_f$X), ]
rownames(ctrl_f) = ctrl_f$X
ctrl_f$X = NULL

# combine the 2 dfs
ra_early_f_vs_ctrl_f = cbind(ctrl_f, female_early_ra)
colnames(ra_early_f_vs_ctrl_f)[1:14] = paste0("C", 1:14)
colnames(ra_early_f_vs_ctrl_f)[15:47] = paste0("Fer", 1:33)

d0 =  DGEList(ra_early_f_vs_ctrl_f)
cutoff = 1
drop = which(apply(cpm(d0), 1, max) < cutoff)
d = d0[-drop,] 
dim(d) # number of genes left
snames = colnames(d) # Sample names
snames
group = as.factor(c(rep("C", 14), rep("Fer", 33)))
head(group)
plotMDS(d, col = as.numeric(group))
mm = model.matrix(~0 + group)
y = voom(d, mm, plot = T)
normalized_counts = y$E

fit = lmFit(y, mm)
head(coef(fit))
contr = makeContrasts(groupC - groupFer, levels = colnames(coef(fit)))
tmp = contrasts.fit(fit, contr)
head(tmp)
######Empirical Bayes smoothing of standard errors 
tmp = eBayes(tmp)
head(tmp)
ctrl_fer_res = topTable(tmp, sort.by = "P", n = Inf)
ctrl_fer_significant = subset(ctrl_fer_res, ctrl_fer_res$adj.P.Val <= 0.05)
write.csv(ctrl_fer_significant, file = "female_control_vs_female_early_ra_limma_allgenes.csv")
ctrl_fer_up_1_5 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC >= 1.5)
ctrl_fer_down_1_5 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC <= -1.5)
ctrl_fer_up_1 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC >= 1)
ctrl_fer_down_1 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC <= -1)
ctrl_fer_up_0_5 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC >= 0.5)
ctrl_fer_down_0_5 = subset(ctrl_fer_significant, ctrl_fer_significant$logFC <= -0.5)

# Get top 25 upregulated genes
top_25_up_ctrl_fer <- ctrl_fer_res %>%
  filter(logFC >= 1.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ctrl_fer <- ctrl_fer_res %>%
  filter(logFC <= -1.5 & adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts, rownames(normalized_counts) %in% 
                                   rownames(top_25_up_ctrl_fer))

normalized_counts_df_down = subset(normalized_counts, rownames(normalized_counts) %in% 
                                     rownames(top_25_down_ctrl_fer))
library(pheatmap)

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
                        "C10", "C11", "C12", "C13", "C14", "Fer1", "Fer2", "Fer3",
                        "Fer4", "Fer5", "Fer6", "Fer7", "Fer8", "Fer9", "Fer10", 
                        "Fer11", "Fer12", "Fer13", "Fer14", "Fer15", "Fer16", "Fer17",
                        "Fer18", "Fer19", "Fer20", "Fer21", "Fer22", "Fer23", "Fer24", 
                        "Fer25", "Fer26", "Fer27", "Fer28", "Fer29", "Fer30", "Fer31",
                        "Fer32", "Fer33"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Female Early RA")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                        "C10", "C11", "C12", "C13", "C14", "Fer1", "Fer2", "Fer3",
                        "Fer4", "Fer5", "Fer6", "Fer7", "Fer8", "Fer9", "Fer10", 
                        "Fer11", "Fer12", "Fer13", "Fer14", "Fer15", "Fer16", "Fer17",
                        "Fer18", "Fer19", "Fer20", "Fer21", "Fer22", "Fer23", "Fer24", 
                        "Fer25", "Fer26", "Fer27", "Fer28", "Fer29", "Fer30", "Fer31",
                        "Fer32", "Fer33"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Female Early RA")


ctrl_fer_res$neg_logp_val = -log10(ctrl_fer_res$adj.P.Val)
ctrl_fer_res$Expression = ifelse(ctrl_fer_res$logFC >= 1.0 & ctrl_fer_res$adj.P.Val <= 0.05, "Upregulated",
                                 ifelse(ctrl_fer_res$logFC <= -1 & 
                                          ctrl_fer_res$adj.P.Val <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(ctrl_fer_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Female Early RA")












# Male Estd RA vs Male Controls

male_estd_ra = read.csv(file = "ra_estd_m_coding_final.csv", row.names = 1)
length(unique(male_estd_ra$X))

# Remove Duplicates
male_estd_ra = male_estd_ra[!duplicated(male_estd_ra$X), ]
rownames(male_estd_ra) = male_estd_ra$X
male_estd_ra$X = NULL

ctrl_m = read.csv(file = "synovium_controls_m_coding_final.csv", row.names = 1)

# Remove Duplicates
ctrl_m= ctrl_m[!duplicated(ctrl_m$X), ]
rownames(ctrl_m) = ctrl_m$X
ctrl_m$X = NULL

# combine the 2 dfs
ra_estd_m_vs_ctrl_m = cbind(ctrl_m, male_estd_ra)
colnames(ra_estd_m_vs_ctrl_m)[1:14] = paste0("C", 1:14)
colnames(ra_estd_m_vs_ctrl_m)[15:36] = paste0("Mes", 1:22)

# Run Limma:Voom
d0 =  DGEList(ra_estd_m_vs_ctrl_m)
cutoff = 1
drop = which(apply(cpm(d0), 1, max) < cutoff)
d = d0[-drop,] 
dim(d) # number of genes left
snames = colnames(d) # Sample names
snames
group = as.factor(c(rep("C", 14), rep("Mes", 22)))
head(group)
plotMDS(d, col = as.numeric(group))
mm = model.matrix(~0 + group)
y = voom(d, mm, plot = T)
normalized_counts = y$E

fit = lmFit(y, mm)
head(coef(fit))
contr = makeContrasts(groupC - groupMes, levels = colnames(coef(fit)))
tmp = contrasts.fit(fit, contr)
head(tmp)
######Empirical Bayes smoothing of standard errors 
tmp = eBayes(tmp)
head(tmp)
ctrl_mes_res = topTable(tmp, sort.by = "P", n = Inf)
ctrl_mes_significant = subset(ctrl_mes_res, ctrl_mes_res$adj.P.Val <= 0.05)
write.csv(ctrl_mes_significant, file = "male_control_vs_male_estd_ra_limma_allgenes.csv")
ctrl_mes_up_1_5 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC >= 1.5)
ctrl_mes_down_1_5 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC <= -1.5)
ctrl_mes_up_1 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC >= 1)
ctrl_mes_down_1 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC <= -1)
ctrl_mes_up_0_5 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC >= 0.5)
ctrl_mes_down_0_5 = subset(ctrl_mes_significant, ctrl_mes_significant$logFC <= -0.5)

# Get top 25 upregulated genes
top_25_up_ctrl_mes <- ctrl_mes_res %>%
  filter(logFC >= 1.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ctrl_mes <- ctrl_mes_res %>%
  filter(logFC <= -1.5 & adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts, rownames(normalized_counts) %in% 
                                   rownames(top_25_up_ctrl_mes))

normalized_counts_df_down = subset(normalized_counts, rownames(normalized_counts) %in% 
                                     rownames(top_25_down_ctrl_mes))
library(pheatmap)

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
                        "C10", "C11", "C12", "C13", "C14", "Mes1", "Mes2", "Mes3", 
                        "Mes4", "Mes5", "Mes6", "Mes7", "Mes8", "Mes9", "Mes10", 
                        "Mes11", "Mes12", "Mes13", "Mes14", "Mes15", "Mes16", 
                        "Mes17", "Mes18", "Mes19", "Mes20", "Mes21", "Mes22"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Male Estd RA")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                        "C10", "C11", "C12", "C13", "C14", "Mes1", "Mes2", "Mes3", 
                        "Mes4", "Mes5", "Mes6", "Mes7", "Mes8", "Mes9", "Mes10", 
                        "Mes11", "Mes12", "Mes13", "Mes14", "Mes15", "Mes16", 
                        "Mes17", "Mes18", "Mes19", "Mes20", "Mes21", "Mes22"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Male Estd RA")


ctrl_mes_res$neg_logp_val = -log10(ctrl_mes_res$adj.P.Val)
ctrl_mes_res$Expression = ifelse(ctrl_mes_res$logFC >= 1.0 & ctrl_mes_res$adj.P.Val <= 0.05, "Upregulated",
                                 ifelse(ctrl_mes_res$logFC <= -1.0 & 
                                          ctrl_mes_res$adj.P.Val <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(ctrl_mes_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Male Estd RA")

















# Female Estd RA vs Male Controls

female_estd_ra = read.csv(file = "ra_estd_f_coding_final.csv", row.names = 1)
length(unique(female_estd_ra$X))

# Remove Duplicates
female_estd_ra = female_estd_ra[!duplicated(female_estd_ra$X), ]
rownames(female_estd_ra) = female_estd_ra$X
female_estd_ra$X = NULL

ctrl_f = read.csv(file = "synovium_controls_f_coding_final.csv", row.names = 1)

# Remove Duplicates
ctrl_f = ctrl_f[!duplicated(ctrl_f$X), ]
rownames(ctrl_f) = ctrl_f$X
ctrl_f$X = NULL

# combine the 2 dfs
ra_estd_f_vs_ctrl_f = cbind(ctrl_f, female_estd_ra)
colnames(ra_estd_f_vs_ctrl_f)[1:14] = paste0("C", 1:14)
colnames(ra_estd_f_vs_ctrl_f)[15:87] = paste0("Fes", 1:73)

# Run Limma:Voom
d0 =  DGEList(ra_estd_f_vs_ctrl_f)
cutoff = 1
drop = which(apply(cpm(d0), 1, max) < cutoff)
d = d0[-drop,] 
dim(d) # number of genes left
snames = colnames(d) # Sample names
snames
group = as.factor(c(rep("C", 14), rep("Fes", 73)))
head(group)
plotMDS(d, col = as.numeric(group))
mm = model.matrix(~0 + group)
y = voom(d, mm, plot = T)
normalized_counts = y$E

fit = lmFit(y, mm)
head(coef(fit))
contr = makeContrasts(groupC - groupFes, levels = colnames(coef(fit)))
tmp = contrasts.fit(fit, contr)
head(tmp)
######Empirical Bayes smoothing of standard errors 
tmp = eBayes(tmp)
head(tmp)
ctrl_fes_res = topTable(tmp, sort.by = "P", n = Inf)
ctrl_fes_significant = subset(ctrl_fes_res, ctrl_fes_res$adj.P.Val <= 0.05)
write.csv(ctrl_fes_significant, file = "female_control_vs_female_estd_ra_limma_allgenes.csv")
ctrl_fes_up_1_5 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC >= 1.5)
ctrl_fes_down_1_5 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC <= -1.5)
ctrl_fes_up_1 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC >= 1)
ctrl_fes_down_1 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC <= -1)
ctrl_fes_up_0_5 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC >= 0.5)
ctrl_fes_down_0_5 = subset(ctrl_fes_significant, ctrl_fes_significant$logFC <= -0.5)

# Get top 25 upregulated genes
top_25_up_ctrl_fes <- ctrl_fes_res %>%
  filter(logFC >= 1.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ctrl_fes <- ctrl_fes_res %>%
  filter(logFC <= -1.5 & adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts, rownames(normalized_counts) %in% 
                                   rownames(top_25_up_ctrl_fes))

normalized_counts_df_down = subset(normalized_counts, rownames(normalized_counts) %in% 
                                     rownames(top_25_down_ctrl_fes))
library(pheatmap)

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
                        "C10", "C11", "C12", "C13", "C14", "Fes1", "Fes2", "Fes3", 
                        "Fes4", "Fes5", "Fes6", "Fes7", "Fes8", "Fes9", "Fes10", 
                        "Fes11", "Fes12", "Fes13", "Fes14", "Fes15", "Fes16", "Fes17", 
                        "Fes18", "Fes19", "Fes20", "Fes21", "Fes22", "Fes23", "Fes24", 
                        "Fes25", "Fes26", "Fes27", "Fes28", "Fes29", "Fes30", "Fes31", 
                        "Fes32", "Fes33", "Fes34", "Fes35", "Fes36", "Fes37", "Fes38", 
                        "Fes39", "Fes40", "Fes41", "Fes42", "Fes43", "Fes44", "Fes45", 
                        "Fes46", "Fes47", "Fes48", "Fes49", "Fes50", "Fes51", "Fes52", 
                        "Fes53", "Fes54", "Fes55", "Fes56", "Fes57", "Fes58", "Fes59", 
                        "Fes60", "Fes61", "Fes62", "Fes63", "Fes64", "Fes65", "Fes66", 
                        "Fes67", "Fes68", "Fes69", "Fes70", "Fes71", "Fes72", "Fes73"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Female Estd RA")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                        "C10", "C11", "C12", "C13", "C14", "Fes1", "Fes2", "Fes3", 
                        "Fes4", "Fes5", "Fes6", "Fes7", "Fes8", "Fes9", "Fes10", 
                        "Fes11", "Fes12", "Fes13", "Fes14", "Fes15", "Fes16", "Fes17", 
                        "Fes18", "Fes19", "Fes20", "Fes21", "Fes22", "Fes23", "Fes24", 
                        "Fes25", "Fes26", "Fes27", "Fes28", "Fes29", "Fes30", "Fes31", 
                        "Fes32", "Fes33", "Fes34", "Fes35", "Fes36", "Fes37", "Fes38", 
                        "Fes39", "Fes40", "Fes41", "Fes42", "Fes43", "Fes44", "Fes45", 
                        "Fes46", "Fes47", "Fes48", "Fes49", "Fes50", "Fes51", "Fes52", 
                        "Fes53", "Fes54", "Fes55", "Fes56", "Fes57", "Fes58", "Fes59", 
                        "Fes60", "Fes61", "Fes62", "Fes63", "Fes64", "Fes65", "Fes66", 
                        "Fes67", "Fes68", "Fes69", "Fes70", "Fes71", "Fes72", "Fes73"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Female Estd RA")


ctrl_fes_res$neg_logp_val = -log10(ctrl_fes_res$adj.P.Val)
ctrl_fes_res$Expression = ifelse(ctrl_fes_res$logFC >= 1.0 & ctrl_fes_res$adj.P.Val <= 0.05, "Upregulated",
                                 ifelse(ctrl_fes_res$logFC <= -1.0 & 
                                          ctrl_fes_res$adj.P.Val <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(ctrl_fes_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Female Estd RA")
