#################################################
# Part 2: UMAP, Heatmap and Spearman-Kendall matrix
#################################################

# =========================================
# 1. TRPM4 UMAP plot
# =========================================
# Download and then load the RData available from HuggingFace: https://huggingface.co/datasets/kahkengwong/GAM_REML_PRSS_Project/tree/main
load("GSE185344_Seurat_processed.RData") # Filesize 9.516 GB

# Plot UMAP with log2 gene expression to visualize gene levels
plot_umap_gene_expression <- function(seurat_obj, gene, dataset_label, assay = "RNA") {
    DefaultAssay(seurat_obj) <- assay  # Set assay (RNA default)
    gene_colors_alpha <- c(scales::alpha("lightgray", 0.85), scales::alpha("lightpink", 0.85), scales::alpha("#FF6666", 0.85), 
                           scales::alpha("#BC2727", 0.85), scales::alpha("#660000", 0.85))  # Color palette (semi-transparent)
    umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))  # UMAP embeddings (coords)
    gene_expression <- FetchData(seurat_obj, vars = gene)[, 1]; umap_data$expression <- log2(gene_expression + 1)  # Log2 expression (pseudo-count 1)
    umap_data <- umap_data[order(umap_data$expression), ]  # Sort by expression (high expressers front)
    feature_plot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = expression)) +  # UMAP plot (gene-colored)
        geom_point(size = 0.2, alpha = 0.85) + scale_color_gradientn(colors = gene_colors_alpha, limits = range(umap_data$expression, na.rm = TRUE)) +
        theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
        labs(title = paste("UMAP log2", gene, "expression", dataset_label), x = "umap_1", y = "umap_2", color = "log2")
    print(feature_plot)
}

# Marker genes, TRPM4 and KLK4 that serves as an independent marker 
markers <- c("TRPM4", "KLK4")

# Plot UMAPs for prostate ca
for (marker in markers) {
    plot_umap_gene_expression(prostate_results$seurat_obj, marker, "_PCa", assay = "RNA")
}

# Plot UMAPs for non-cancerous
for (marker in markers) {
    plot_umap_gene_expression(non_cancerous_results$seurat_obj, marker, "_Non-cancerous", assay = "RNA")
}


# ======================================================
# 2. Spearman's r and Kendall's tau Matrix between Different Clusters
# ======================================================
# Get average expression for clusters
get_cluster_averages <- function(seurat_obj, clusters) {
    Idents(seurat_obj) <- "seurat_clusters"
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)  # Log-normalize (scale 10000)
    avg_expr <- AverageExpression(seurat_obj, idents = clusters, slot = "data")$RNA  # Average expression (RNA assay)
    return(avg_expr)
}

# Compute averages (non-cancerous cluster 3, prostate ca 0-22)
non_ca_avg <- get_cluster_averages(non_cancerous_results$seurat_obj, "3")
prostate_ca_avg <- get_cluster_averages(prostate_results$seurat_obj, 0:22)

# Combine data (matrix of averages)
all_clusters_avg <- cbind(non_ca_avg[,"g3"], prostate_ca_avg)
all_clusters_avg <- as.matrix(all_clusters_avg)
colnames(all_clusters_avg) <- c("Non-Ca_3", paste0("PCa_", 0:22))

# Filter low-expression genes (keep finite, >0)
genes_to_keep <- rowSums(is.finite(all_clusters_avg) & all_clusters_avg > 0) == ncol(all_clusters_avg)
all_clusters_avg <- all_clusters_avg[genes_to_keep, ]

# Spearman's r (non-parametric)
cor_matrix_spearman <- cor(all_clusters_avg, method = "spearman")

# Kendall's tau (rank-based) with parallel processing by using all CPU cores (except 1) to speed up the computations
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
n <- ncol(all_clusters_avg)
cluster_names <- colnames(all_clusters_avg)
kendall_cor <- function(i, j, data) { cor(data[,i], data[,j], method = "kendall") }  # Kendall function (pairwise)
cor_matrix_kendall_upper <- foreach(i = 1:(n-1), .combine = 'c') %dopar% {
    sapply((i+1):n, function(j) kendall_cor(i, j, all_clusters_avg))
}
stopCluster(cl)
cor_matrix_kendall <- matrix(0, nrow = n, ncol = n)
index <- 1
for (i in 1:(n-1)) {
    for (j in (i+1):n) {
        cor_matrix_kendall[i,j] <- cor_matrix_kendall_upper[index]
        index <- index + 1
    }
}
cor_matrix_kendall <- cor_matrix_kendall + t(cor_matrix_kendall)
diag(cor_matrix_kendall) <- 1
rownames(cor_matrix_kendall) <- cluster_names; colnames(cor_matrix_kendall) <- cluster_names

# Plot correlation matrix to visualize relationships
plot_correlation_matrix <- function(cor_matrix, title, palette, is_reversed = FALSE) {
    colors <- brewer.pal(11, palette); if (is_reversed) colors <- rev(colors)
    corrplot(cor_matrix, method = "color", col = colors, addCoef.col = "black", tl.col = "black", tl.srt = 45, 
             title = title, mar = c(0,0,1,0))  # Plot (color + coefficients)
}

# Save matrices as CSVs
write.csv(cor_matrix_spearman, "cluster_correlation_matrix_spearman.csv")
write.csv(cor_matrix_kendall, "cluster_correlation_matrix_kendall.csv")


# ==============================================================================
# 3. Spearman's r and Kendall's tau for TRPM4 in PCa Combined Clusters and Adjacent BP/Non-ca
# ==============================================================================
# Calculate correlations of TRPM4 vs all genes
calculate_correlations <- function(integrated_obj, original_obj, cluster_ids, gene_of_interest, output_file) {
    cluster_cells <- WhichCells(integrated_obj, idents = cluster_ids)  # Select cluster cells (integrated)
    cluster_obj <- subset(original_obj, cells = cluster_cells)  # Subset original object (pre-integration data)
    DefaultAssay(cluster_obj) <- "RNA"  # Set assay (RNA)
    gene_expression <- FetchData(cluster_obj, vars = gene_of_interest)[, 1]  # Gene expression (TRPM4)
    normalized_data <- GetAssayData(cluster_obj, slot = "data", assay = "RNA")  # Normalized data (all genes)
    num_cores <- max(1, parallel::detectCores() - 1)  # Parallel setup (max cores - 1)
    correlations <- pbapply(normalized_data, 1, function(gene_expr) {  # Compute correlations (with progress)
        c(Kendall_Correlation = cor(gene_expression, gene_expr, method = "kendall"),
          Spearman_Correlation = cor(gene_expression, gene_expr, method = "spearman"))
    }, cl = num_cores)
    correlations_df <- data.frame(Gene = rownames(normalized_data), Kendall_Correlation = correlations["Kendall_Correlation", ],
                                  Spearman_Correlation = correlations["Spearman_Correlation", ])  # Results dataframe
    write.xlsx(correlations_df, file = output_file, rowNames = FALSE)  # Export (Excel)
}

# Compute correlations for clusters focusing on TRPM4
combined_clusters <- c(6, 9, 11, 14, 19)
output_file_combined <- "TRPM4_Correlations_PCa_Clusters_Combined.xlsx"
calculate_correlations(prostate_results$seurat_obj, prostate_ca_seurat, combined_clusters, "TRPM4", output_file_combined)
output_file_cluster16 <- "TRPM4_Correlations_PCa_Cluster_16.xlsx"
calculate_correlations(prostate_results$seurat_obj, prostate_ca_seurat, 16, "TRPM4", output_file_cluster16)
output_file_non_cancerous <- "TRPM4_Correlations_BPNonCa.xlsx"
calculate_correlations(non_cancerous_results$seurat_obj, non_cancerous_seurat, 3, "TRPM4", output_file_non_cancerous)


# =========================================
# 4. Heatmap of TRPM4 vs Ribosomal and AR Genes
# =========================================
# Subset data (specific clusters for heatmap)
heatmap_prostate_subset <- subset(prostate_results$seurat_obj, idents = c(6, 9, 11, 14, 19, 16))
heatmap_noncancer_subset <- subset(non_cancerous_results$seurat_obj, idents = 3)
heatmap_obj <- merge(heatmap_prostate_subset, y = heatmap_noncancer_subset, add.cell.ids = c("PCa", "NonCa"))

# Define gene sets (ribosomal and AR-related)
genes_of_interest <- c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26", "KLK4", "KLK2", "KLK3", "PDLIM5", "ABHD2", "ALDH1A3", "SORD")
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(heatmap_obj)]
ribo_genes <- c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26"); ribo_genes <- ribo_genes[ribo_genes %in% rownames(heatmap_obj)]
ar_genes <- c("KLK4", "KLK2", "KLK3", "PDLIM5", "ABHD2", "ALDH1A3", "SORD"); ar_genes <- ar_genes[ar_genes %in% rownames(heatmap_obj)]

# Compute mean expression (aggregate ribosomal and AR)
ribo_expression <- colMeans(GetAssayData(heatmap_obj, slot = "data")[ribo_genes, ])
ar_expression <- colMeans(GetAssayData(heatmap_obj, slot = "data")[ar_genes, ])
trpm4_expression <- GetAssayData(heatmap_obj, slot = "data")["TRPM4", ]

# Prepare expression matrix (log2-transformed)
expression_matrix <- GetAssayData(heatmap_obj, slot = "data")[c("TRPM4", genes_of_interest), ]
expression_matrix_log <- log2(expression_matrix + 1)
expression_matrix_log <- rbind(expression_matrix_log["TRPM4", , drop = FALSE], Ribo = log2(ribo_expression + 1),
                               AR = log2(ar_expression + 1), expression_matrix_log[setdiff(rownames(expression_matrix_log), "TRPM4"), ])

# Scale to z-scores (cap at -1 to 2.3)
scale_to_zscore <- function(x) {
    z <- scale(x)
    return(pmin(pmax(z, -1), 2.3))
}
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Order cells by cluster and TRPM4 expression (ascending)
heatmap_grouping <- Idents(heatmap_obj)
heatmap_grouping <- factor(heatmap_grouping, levels = c("6", "9", "11", "14", "19", "16", "3"))
cell_order <- order(heatmap_grouping, trpm4_expression)
expression_matrix_log <- expression_matrix_log[, cell_order]
heatmap_grouping <- heatmap_grouping[cell_order]
trpm4_expression <- trpm4_expression[cell_order]
ribo_expression <- ribo_expression[cell_order]
ar_expression <- ar_expression[cell_order]
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Prepare heatmap annotations (cluster labels with counts)
n_groups <- length(levels(heatmap_grouping))
group_colors <- viridis(n_groups, option = "mako"); names(group_colors) <- levels(heatmap_grouping)
cell_counts <- table(heatmap_grouping); cell_percentages <- round(prop.table(cell_counts) * 100, 1)
group_labels <- paste0(names(cell_counts), "\n", cell_counts, "\n(", cell_percentages, "%)")
column_annotation <- HeatmapAnnotation(Cluster = anno_block(gp = gpar(fill = group_colors), labels = group_labels,
                                                            labels_gp = gpar(col = "white", fontsize = 8)),
                                       show_legend = FALSE, show_annotation_name = FALSE)

# Create heatmap of TRPM4 vs gene sets
heatmap <- Heatmap(expression_matrix_scaled, name = "zscore", column_title = NULL, row_title = "Genes",
                   show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                   top_annotation = column_annotation, col = viridis(100, option = "mako"), row_names_gp = gpar(fontsize = 8),
                   column_split = heatmap_grouping, row_gap = unit(1, "mm"), column_gap = unit(0.5, "mm"), border = TRUE,
                   use_raster = TRUE, raster_quality = 2)

# Save and display heatmap
pdf("TRPM4_heatmap.pdf", width = 12, height = 8)
draw(heatmap)
dev.off()
draw(heatmap)



