###############################################
# >>> Dataset GSE185344 scRNA-seq analysis <<<
###############################################
#############################################################
# Part 1: scRNA-seq Dataset Pre-processing and UMAP Clusters
#############################################################
setwd("C:/...") # Set the local directory
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(dorothea)
library(dplyr)
library(earth)
library(future)
library(ggplot2)
library(ggplot2) 
library(ggrepel)
library(glmGamPoi) 
library(Hmisc)
library(Matrix)
library(matrixStats)
library(mgcv)
library(monocle3)
library(openxlsx)
library(parallel)
library(pbapply)
library(pheatmap)
library(presto)
library(purrr)
library(RColorBrewer)
library(scales)
library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)
library(SparseArray)
library(tibble)
library(tidyr)
library(viper)
library(viridis)
library(writexl)
library(FSA)
library(rstatix)
library(readxl)
library(magick)
library(pdftools)
library(rmarkdown)


# =========================================
# 1. scRNA-seq Dataset Pre-processing
# =========================================
# Load data, GSE185344 prostate cancer scRNA-seq dataset
loaded_df <- readRDS("C:/...directory.../GSE185344_PH_scRNA.final.rds")

# Extract Seurat object, the core data structure
seurat_obj <- loaded_df$obj

# Define sample names (tumor vs benign)
prostate_ca_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                         "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                         "HYW_5755_Tumor")
non_cancerous_samples <- c("HYW_4701_Benign", "HYW_4847_Benign", "HYW_4880_Benign", 
                           "HYW_4881_Benign", "HYW_5386_Benign", "HYW_5742_Benign", 
                           "HYW_5755_Benign")

# Subset and process Seurat object for normalization and feature selection
process_seurat <- function(seurat_obj, sample_names, project_name) {
    subset_obj <- subset(seurat_obj, subset = orig.ident %in% sample_names)  # Subset by sample ID
    subset_obj <- NormalizeData(subset_obj, normalization.method = "LogNormalize", scale.factor = 10000)  # Log-normalize (expression scaling)
    subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000)  # Top 2000 variable features (VST method)
    return(subset_obj)
}

# Process samples to split into tumor and benign
prostate_ca_seurat <- process_seurat(seurat_obj, prostate_ca_samples, "prostate-ca")
non_cancerous_seurat <- process_seurat(seurat_obj, non_cancerous_samples, "NonCancerous")

# Filter by feature and count thresholds to remove low-quality cells
prostate_ca_seurat <- subset(prostate_ca_seurat, subset = nFeature_RNA > 500 & nCount_RNA > 0)
non_cancerous_seurat <- subset(non_cancerous_seurat, subset = nFeature_RNA > 500 & nCount_RNA > 0)

# Filter high ribosomal content to mitigate bias from cells with highest expression of ribosomal genes (top 10th percentile)
filter_ribosomal <- function(seurat_obj, method = "fixed", cutoff = 10) {
    rp_genes <- grep("^RP[SL]|^MRP[SL]", rownames(seurat_obj), value = TRUE)  # Ribosomal genes (RP/MRP prefixes)
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = rp_genes)  # % ribosomal expression (per cell)
    threshold <- if (method == "percentile") quantile(seurat_obj$percent.ribo, probs = cutoff) else cutoff  # Threshold (percentile or fixed)
    plot <- ggplot(seurat_obj@meta.data, aes(x = percent.ribo)) +  # Plot distribution (with threshold line)
        geom_histogram(bins = 100) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        ggtitle("Distribution of Ribosomal Gene Percentage")
    print(plot)
    genes_before <- nrow(seurat_obj); cells_before <- ncol(seurat_obj)  # Pre-filter counts (genes, cells)
    seurat_obj_filtered <- subset(seurat_obj, subset = percent.ribo < threshold)  # Apply filter (below threshold)
    genes_after <- nrow(seurat_obj_filtered); cells_after <- ncol(seurat_obj_filtered)  # Post-filter counts
    cat("Threshold:", threshold, "\nCells before:", cells_before, "\nCells after:", cells_after, 
        "\nRemoved:", round((cells_before - cells_after) / cells_before * 100, 2), "%\n",
        "Genes before:", genes_before, "\nGenes after:", genes_after, "\n")
    return(seurat_obj_filtered)
}

# Apply ribosomal filter (90th percentile cutoff)
cat("Filtering prostate cancer samples (ribosomal)\n")
prostate_ca_seurat <- filter_ribosomal(prostate_ca_seurat, method = "percentile", cutoff = 0.90)
cat("Filtering non-cancerous samples (ribosomal)\n")
non_cancerous_seurat <- filter_ribosomal(non_cancerous_seurat, method = "percentile", cutoff = 0.90)

# Filter high mitochondrial content to mitigate dying cells
filter_mitochondrial <- function(seurat_obj, method = "fixed", cutoff = 10) {
    mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)  # Mitochondrial genes (MT- prefix)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes)  # % mitochondrial expression (per cell)
    threshold <- if (method == "percentile") quantile(seurat_obj$percent.mt, probs = cutoff) else cutoff  # Threshold (percentile or fixed)
    plot <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +  # Plot distribution (with threshold line)
        geom_histogram(bins = 100) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        ggtitle("Distribution of Mitochondrial Gene Percentage")
    print(plot)
    genes_before <- nrow(seurat_obj); cells_before <- ncol(seurat_obj)  # Pre-filter counts (genes, cells)
    seurat_obj_filtered <- subset(seurat_obj, subset = percent.mt < threshold)  # Apply filter (below threshold)
    genes_after <- nrow(seurat_obj_filtered); cells_after <- ncol(seurat_obj_filtered)  # Post-filter counts
    cat("Threshold:", threshold, "\nCells before:", cells_before, "\nCells after:", cells_after, 
        "\nRemoved:", round((cells_before - cells_after) / cells_before * 100, 2), "%\n",
        "Genes before:", genes_before, "\nGenes after:", genes_after, "\n")
    return(seurat_obj_filtered)
}

# Apply mitochondrial filter (90th percentile cutoff)
cat("Filtering prostate cancer samples (mitochondrial)\n")
prostate_ca_seurat <- filter_mitochondrial(prostate_ca_seurat, method = "percentile", cutoff = 0.90)
cat("Filtering non-cancerous samples (mitochondrial)\n")
non_cancerous_seurat <- filter_mitochondrial(non_cancerous_seurat, method = "percentile", cutoff = 0.90)


# =========================================
# 2. Cell Cycle Regression
# =========================================
# Check pre-regression cell count (after the previous steps)
cat("Cells before cell cycle regression (prostate cancer):", ncol(prostate_ca_seurat), "\n")

# Default cell cycle genes (common S and G2M phase markers)
s_genes_default <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m_genes_default <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

# Score and regress out cell cycle to remove phase effects (prostate cancer)
prostate_ca_seurat <- CellCycleScoring(prostate_ca_seurat, s.features = s_genes_default, g2m.features = g2m_genes_default, set.ident = TRUE)
prostate_ca_seurat <- ScaleData(prostate_ca_seurat, vars.to.regress = c("S.Score", "G2M.Score"))

# Post-regression cell count to verify no cell loss
cat("Cells after cell cycle regression (prostate cancer):", ncol(prostate_ca_seurat), "\n") # 22796 (no changes)

# Score and regress out cell cycle (non-cancerous)
cat("Cells before cell cycle regression (non-cancerous):", ncol(non_cancerous_seurat), "\n")
non_cancerous_seurat <- CellCycleScoring(non_cancerous_seurat, s.features = s_genes_default, g2m.features = g2m_genes_default, set.ident = TRUE)
non_cancerous_seurat <- ScaleData(non_cancerous_seurat, vars.to.regress = c("S.Score", "G2M.Score"))


# =========================================
# 3. Doublets removal
# =========================================
# Set up parallel processing to speed up doublet detection
plan(multisession, workers = availableCores())

# Remove doublets using scDblFinder
remove_doublets <- function(seurat_obj) {
    sce_obj <- as.SingleCellExperiment(seurat_obj)  # Convert to SCE for scDblFinder
    samples <- seurat_obj@meta.data$orig.ident  # Batch info by sample IDs
    doublet_scores <- scDblFinder(sce_obj, samples = samples, k = 30, nfeatures = 2000)  # Run scDblFinder 
    batch_thresholds <- tapply(doublet_scores$scDblFinder.score, samples, function(x) quantile(x, probs = 0.95))  # Batch-specific thresholds (95th percentile)
    cat("Batch thresholds:\n"); print(batch_thresholds)
    doublet_cells <- colnames(sce_obj)[mapply(function(x, y) x > batch_thresholds[y], doublet_scores$scDblFinder.score, samples)]  # Identify doublets (above threshold)
    cat("Doublets:", length(doublet_cells), "\n")
    seurat_obj$doublet <- colnames(seurat_obj) %in% doublet_cells  # Mark doublets as TRUE/FALSE
    seurat_obj <- subset(seurat_obj, subset = doublet == FALSE)  # Remove doublets
    cat("Cells after removal:", ncol(seurat_obj), "\n")
    seurat_obj$filtered <- "filtered"  # Update metadata and flag filtered cells
    return(seurat_obj)
}

# Apply doublet removal to prostate ca and benign cases
cat("Processing prostate cancer samples\n")
prostate_ca_seurat <- remove_doublets(prostate_ca_seurat)
cat("Processing non-cancerous samples\n")
non_cancerous_seurat <- remove_doublets(non_cancerous_seurat)


# =========================================
# 4. Batch effects correction
# =========================================
# Disable parallel processing to avoid integration issues
plan(sequential)

# Correct batch effects by integrating across samples
correct_batch_effects <- function(seurat_obj) {
    cat("Metadata columns:\n"); print(colnames(seurat_obj@meta.data))
    cat("Unique orig.ident:\n"); print(unique(seurat_obj@meta.data$orig.ident))
    seurat_obj_before_integration <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 500)  # Pre-integration features (500 VST)
    seurat_obj_before_integration <- RunPCA(seurat_obj_before_integration, verbose = FALSE)  # PCA (dimensionality reduction)
    seurat_obj_before_integration <- RunUMAP(seurat_obj_before_integration, dims = 1:8, verbose = FALSE)  # UMAP (pre-integration visualization)
    plasma_colors <- viridis(n = length(unique(seurat_obj_before_integration$orig.ident)), option = "plasma")
    p1 <- DimPlot(seurat_obj_before_integration, group.by = "orig.ident", pt.size = 0.5, label = FALSE, repel = TRUE, cols = plasma_colors) + 
        ggtitle("UMAP Before Integration") + theme(legend.position = "right")  # Pre-integration UMAP (batch-colored)
    print(p1)
    sample_list <- SplitObject(seurat_obj_before_integration, split.by = "orig.ident")  # Split by batch (sample IDs)
    sample_list <- lapply(sample_list, function(x) {  # SCTransform and clean NAs (per sample)
        x <- SCTransform(x, verbose = FALSE, variable.features.n = 500, vst.flavor = "v2")
        x@meta.data <- x@meta.data[complete.cases(x@meta.data), ]
        x
    })
    anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:5, verbose = FALSE)  # Find anchors (dims 1-5)
    seurat_obj_integrated <- IntegrateData(anchorset = anchors, dims = 1:5, verbose = FALSE)  # Integrate (batch-corrected)
    seurat_obj_integrated <- ScaleData(seurat_obj_integrated, verbose = FALSE)  # Scale (post-integration)
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)  # PCA (integrated)
    seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, dims = 1:8)  # Neighbors (for clustering)
    seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.5)  # Clusters (res 0.5)
    seurat_obj_integrated <- RunUMAP(seurat_obj_integrated, dims = 1:8, verbose = FALSE, umap.method = "uwot", metric = "cosine")  # UMAP (post-integration)
    plasma_colors <- viridis(n = length(unique(seurat_obj_integrated$orig.ident)), option = "plasma")
    p3 <- DimPlot(seurat_obj_integrated, group.by = "orig.ident", pt.size = 0.5, label = FALSE, repel = TRUE, cols = plasma_colors) + 
        ggtitle("UMAP After Integration") + theme(legend.position = "right")  # Post-integration UMAP (batch-colored)
    print(p3)
    p4 <- DimPlot(seurat_obj_integrated, group.by = "orig.ident", pt.size = 0.5, label = TRUE, repel = TRUE) + 
        ggtitle("UMAP After Integration") + theme(legend.position = "right")  # Labeled UMAP by batch IDs
    print(p4)
    return(seurat_obj_integrated)
}

# Remove parallelization limits to ensure stability
options(future.globals.maxSize = Inf)

# Apply batch correction for prostate ca and benign cases
cat("Processing prostate cancer samples (batch effects)\n")
prostate_ca_seurat_integrated <- correct_batch_effects(prostate_ca_seurat)
cat("Processing non-cancerous samples (batch effects)\n")
non_cancerous_seurat_integrated <- correct_batch_effects(non_cancerous_seurat)


# =========================================
# 5. UMAP Clusters
# =========================================
# Generate elbow plots to assess PCA dimensionality reduction
generate_elbow_plot <- function(seurat_obj_integrated, output_prefix) {
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)  # PCA (dimensionality reduction)
    elbow_plot <- ElbowPlot(seurat_obj_integrated, ndims = 50) +  # Elbow plot
        labs(title = paste("Elbow Plot for", output_prefix), x = "Principal Components", y = "Standard Deviation") +
        theme(plot.title = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 12), 
              axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
    ggsave(paste0(output_prefix, "_ElbowPlot.pdf"), elbow_plot, width = 6.83, height = 6.41)
    print(elbow_plot)
}

# Generate elbow plots (prostate ca and benign)
generate_elbow_plot(prostate_ca_seurat_integrated, "prostate_ca")
generate_elbow_plot(non_cancerous_seurat_integrated, "non_cancerous")

# Downstream analyses with UMAP (clustering and markers)
downstream_analyses <- function(seurat_obj_integrated, gene_of_interest, output_prefix, dims = 15) {
    set.seed(10)
    DefaultAssay(seurat_obj_integrated) <- "RNA"
    seurat_obj_integrated <- FindVariableFeatures(seurat_obj_integrated, selection.method = "vst", nfeatures = 2000)  # Variable features (2000 VST)
    seurat_obj_integrated <- ScaleData(seurat_obj_integrated, verbose = FALSE)  # Scale (center and normalize)
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)  # PCA (dims reduction)
    seurat_obj_integrated <- RunUMAP(seurat_obj_integrated, dims = 1:dims, verbose = FALSE)  # UMAP (dims 1-15)
    set.seed(11); seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, dims = 1:dims)  # Neighbors (kNN graph)
    set.seed(12); seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.5)  # Clusters (Louvain, res 0.5)
    set.seed(13); cluster_markers <- FindAllMarkers(seurat_obj_integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)  # Marker genes (positive, logFC > 0.25)
    print(paste("Cluster markers:", nrow(cluster_markers)))
    if (nrow(cluster_markers) == 0) {
        print("Cluster levels:"); print(levels(Idents(seurat_obj_integrated)))
        print("Cells per cluster:"); print(table(Idents(seurat_obj_integrated)))
    }
    umap_data <- as.data.frame(Embeddings(seurat_obj_integrated, "umap")); umap_data$cluster_id <- Idents(seurat_obj_integrated)  # UMAP data (coords + clusters)
    umap_data_mean <- aggregate(. ~ cluster_id, data = umap_data, FUN = mean)  # Mean coords (per cluster)
    plasma_func <- colorRampPalette(viridis::viridis(100, direction = -1, option = "plasma")); portion <- 0.8  # Colors (plasma palette)
    n_colors <- round(length(unique(umap_data$cluster_id)) / portion); plasma_colors <- plasma_func(n_colors)
    set.seed(14); umap_plot_with_labels <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +  # Labeled UMAP (cluster IDs)
        geom_point(size = 0.3, alpha = 0.5) + scale_color_manual(values = plasma_colors) +
        geom_text(data = umap_data_mean, aes(label = cluster_id, x = umap_1, y = umap_2), color = "black", size = 3, fontface = "bold", check_overlap = TRUE) +
        theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
        labs(title = "UMAP plot colored by cluster (with labels)", x = "umap_1", y = "umap_2", color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3)))
    print(umap_plot_with_labels)
    set.seed(15); umap_plot_no_labels <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +  # Unlabeled UMAP (clusters only)
        geom_point(size = 0.3, alpha = 0.5) + scale_color_manual(values = plasma_colors) +
        theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
        labs(title = "UMAP plot colored by cluster (without labels)", x = "umap_1", y = "umap_2", color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3)))
    print(umap_plot_no_labels)
    if (gene_of_interest %in% rownames(seurat_obj_integrated)) {  # Gene expression UMAP (if gene exists)
        gene_colors_alpha <- c(scales::alpha("lightgray", 0.85), scales::alpha("lightpink", 0.85), scales::alpha("#FF6666", 0.85), 
                               scales::alpha("#BC2727", 0.85), scales::alpha("#660000", 0.85))
        set.seed(16); feature_plot <- FeaturePlot(seurat_obj_integrated, features = gene_of_interest, min.cutoff = 'q10', max.cutoff = 'q90',
                                                  pt.size = 0.2, cols = gene_colors_alpha) +
            theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
                  axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
            labs(title = paste("UMAP plot colored by", gene_of_interest, "expression"), x = "umap_1", y = "umap_2")
        print(feature_plot)
    } else {
        cat(paste("Warning: Gene", gene_of_interest, "not found.\nAvailable genes:\n"))
        print(head(rownames(seurat_obj_integrated), 20))
    }
    if (nrow(cluster_markers) > 0) {  # Top markers (50 per cluster)
        top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
    } else {
        top_markers <- data.frame()
        warning("No cluster markers found.")
    }
    write.table(top_markers, file = paste0(output_prefix, "_top_markers_for_each_cluster_vRibo.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)  # Save markers (TSV)
    return(list(seurat_obj = seurat_obj_integrated, cluster_markers = cluster_markers, top_markers = top_markers))
}

# Apply downstream analyses (TRPM4 focus)
set.seed(42)
prostate_results <- downstream_analyses(prostate_ca_seurat_integrated, "TRPM4", "prostate_ca", dims = 15)
non_cancerous_results <- downstream_analyses(non_cancerous_seurat_integrated, "TRPM4", "non_cancerous", dims = 15)

# Save workspace (~10GB, full analysis state)
save.image(file = "Dt2_scRNAseq_workspace_vRibo_v2.RData")

# For subsequent analysis, load the saved file with: load("Dt2_scRNAseq_workspace_vRibo_v2.RData")
