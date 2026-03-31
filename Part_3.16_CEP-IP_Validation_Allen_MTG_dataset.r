####################################################################
# Validation of CEP-IP Framework: Allen Human MTG Dataset
# Dataset: Allen Human MTG doi: 10.1038/s41586-019-1506-7
# 15,928 nuclei, 8 donors (age 24-66y), MTG, 75 cell types
####################################################################

#############################################################
# Part 1: scRNA-seq Dataset Pre-processing and UMAP Clusters
#############################################################
# =========================================
# 1. scRNA-seq Dataset Pre-processing
# =========================================
# Load required libraries
library(dplyr)
library(future)
library(ggplot2)
library(parallel)
library(scales)
library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
library(viridis)
library(rstatix)
library(writexl)
library(Matrix)  # ADDED: for sparse matrix construction

# ---------------------------------------------------------------------------
# Load Allen MTG SMART-seq dataset from raw CSVs
# NOTE: exon + intron matrices are combined (total RNA per gene per cell)
# Rows = genes, Columns = cells
# ---------------------------------------------------------------------------
cat("Loading Allen MTG SMART-seq data...\n")

# Load gene metadata (rows)
genes_df <- read.csv("human_MTG_2018-06-14_genes-rows.csv", row.names = 1)
cat("Genes loaded:", nrow(genes_df), "\n")

# Load sample/cell metadata (columns)
samples_df <- read.csv("human_MTG_2018-06-14_samples-columns.csv", row.names = 1)
cat("Cells loaded:", nrow(samples_df), "\n")

# Load exon and intron count matrices
# MODIFIED: Read raw CSVs instead of RDS; combine exon + intron for total counts
cat("Loading exon matrix (this may take a few minutes)...\n")
exon_mat <- read.csv("human_MTG_2018-06-14_exon-matrix.csv", row.names = 1, check.names = FALSE)

cat("Loading intron matrix...\n")
intron_mat <- read.csv("human_MTG_2018-06-14_intron-matrix.csv", row.names = 1, check.names = FALSE)

# Combine exon + intron counts (total gene-level expression)
cat("Combining exon + intron matrices...\n")
count_mat <- exon_mat + intron_mat

# Convert to sparse matrix for memory efficiency
count_mat <- as(as.matrix(count_mat), "sparseMatrix")

# Use gene_symbol as row names if available, else keep entrez IDs
if ("gene_symbol" %in% colnames(genes_df)) {
    gene_names <- make.unique(genes_df$gene_symbol)
} else {
    gene_names <- rownames(genes_df)
}
rownames(count_mat) <- gene_names
cat("Matrix dimensions (genes x cells):", nrow(count_mat), "x", ncol(count_mat), "\n")

# ---------------------------------------------------------------------------
# Build Seurat object from count matrix + sample metadata
# MODIFIED: Single cohort (healthy donors) - no tumor/benign split
# ---------------------------------------------------------------------------
seurat_obj <- CreateSeuratObject(
    counts = count_mat,
    meta.data = samples_df,
    project = "AllenMTG"
)
cat("Seurat object created:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")

# ---------------------------------------------------------------------------
# MODIFIED: Single cohort - process entire dataset (no sample split)
# orig.ident maps to donor (column name in samples-columns.csv) for batch correction downstream
# ---------------------------------------------------------------------------
# Ensure orig.ident reflects donor ID for batch correction later
if ("donor" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$orig.ident <- seurat_obj$donor
}

# Normalize and find variable features
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Filter by feature and count thresholds to remove low-quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 0)
cat("Cells after feature/count filter:", ncol(seurat_obj), "\n")

# Filter high ribosomal content (top 10th percentile)
filter_ribosomal <- function(seurat_obj, method = "fixed", cutoff = 10) {
    rp_genes <- grep("^RP[SL]|^MRP[SL]", rownames(seurat_obj), value = TRUE)
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = rp_genes)
    threshold <- if (method == "percentile") quantile(seurat_obj$percent.ribo, probs = cutoff) else cutoff
    plot <- ggplot(seurat_obj@meta.data, aes(x = percent.ribo)) +
        geom_histogram(bins = 100) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        ggtitle("Distribution of Ribosomal Gene Percentage")
    print(plot)
    genes_before <- nrow(seurat_obj); cells_before <- ncol(seurat_obj)
    seurat_obj_filtered <- subset(seurat_obj, subset = percent.ribo < threshold)
    genes_after <- nrow(seurat_obj_filtered); cells_after <- ncol(seurat_obj_filtered)
    cat("Threshold:", threshold, "\nCells before:", cells_before, "\nCells after:", cells_after,
        "\nRemoved:", round((cells_before - cells_after) / cells_before * 100, 2), "%\n",
        "Genes before:", genes_before, "\nGenes after:", genes_after, "\n")
    return(seurat_obj_filtered)
}

cat("Filtering MTG nuclei (ribosomal)\n")
seurat_obj <- filter_ribosomal(seurat_obj, method = "percentile", cutoff = 0.90)

# Mitochondrial content check
# NOTE: Mitochondria have no nucleus; snRNA-seq/nuclear RNA datasets yield near-zero
# MT- gene counts by definition. Applying a percentile filter would remove all cells
# (threshold collapses to 0). The filter is replaced with a diagnostic report only.
mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes)
cat("MT- genes detected:", length(mt_genes), "\n")
cat("percent.mt summary:\n"); print(summary(seurat_obj$percent.mt))
plot_mt <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
    geom_histogram(bins = 100) +
    ggtitle("Distribution of Mitochondrial Gene Percentage (diagnostic only - no filter applied)") +
    xlab("% Mitochondrial") + ylab("Count")
print(plot_mt)
cat("Mitochondrial filter skipped: nuclear RNA dataset (MT- counts are near-zero by definition)\n")
cat("Cells retained after ribosomal filter:", ncol(seurat_obj), "\n")


# =========================================
# 2. Cell Cycle Regression
# =========================================
cat("Cells before cell cycle regression:", ncol(seurat_obj), "\n")

s_genes_default <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m_genes_default <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

seurat_obj <- CellCycleScoring(seurat_obj, s.features = s_genes_default, g2m.features = g2m_genes_default, set.ident = TRUE)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"))

cat("Cells after cell cycle regression:", ncol(seurat_obj), "\n")


# =========================================
# 3. Doublets removal
# NOTE: SMART-seq is plate-based (one cell per well), so doublet rates are
# very low. scDblFinder is retained here as a conservative QC sanity check.
# =========================================
plan(multisession, workers = availableCores())

remove_doublets <- function(seurat_obj) {
    sce_obj <- as.SingleCellExperiment(seurat_obj)
    samples <- seurat_obj@meta.data$orig.ident
    doublet_scores <- scDblFinder(sce_obj, samples = samples, k = 30, nfeatures = 2000)
    batch_thresholds <- tapply(doublet_scores$scDblFinder.score, samples, function(x) quantile(x, probs = 0.95))
    cat("Batch thresholds:\n"); print(batch_thresholds)
    doublet_cells <- colnames(sce_obj)[mapply(function(x, y) x > batch_thresholds[y], doublet_scores$scDblFinder.score, samples)]
    cat("Doublets:", length(doublet_cells), "\n")
    seurat_obj$doublet <- colnames(seurat_obj) %in% doublet_cells
    seurat_obj <- subset(seurat_obj, subset = doublet == FALSE)
    cat("Cells after removal:", ncol(seurat_obj), "\n")
    seurat_obj$filtered <- "filtered"
    return(seurat_obj)
}

cat("Processing MTG nuclei (doublet removal)\n")
seurat_obj <- remove_doublets(seurat_obj)


# =========================================
# 4. Batch effects correction
# MODIFIED: Batch = donor (8 donors), single integrated object
# =========================================
plan(sequential)

correct_batch_effects <- function(seurat_obj) {
    cat("Metadata columns:\n"); print(colnames(seurat_obj@meta.data))
    cat("Unique orig.ident (donors):\n"); print(unique(seurat_obj@meta.data$orig.ident))
    seurat_obj_before_integration <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 500)
    seurat_obj_before_integration <- RunPCA(seurat_obj_before_integration, verbose = FALSE)
    seurat_obj_before_integration <- RunUMAP(seurat_obj_before_integration, dims = 1:8, verbose = FALSE)
    plasma_colors <- viridis(n = length(unique(seurat_obj_before_integration$orig.ident)), option = "plasma")
    p1 <- DimPlot(seurat_obj_before_integration, group.by = "orig.ident", pt.size = 0.5, label = FALSE, repel = TRUE, cols = plasma_colors) +
        ggtitle("UMAP Before Integration (by Donor)") + theme(legend.position = "right")
    print(p1)
    sample_list <- SplitObject(seurat_obj_before_integration, split.by = "orig.ident")
    sample_list <- lapply(sample_list, function(x) {
        x <- SCTransform(x, verbose = FALSE, variable.features.n = 500, vst.flavor = "v2")
        x@meta.data <- x@meta.data[complete.cases(x@meta.data), ]
        x
    })
    anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:5, verbose = FALSE)
    # Dynamically set k.weight to avoid error when smallest donor batch < default k.weight (100)
    # Uses 90% of smallest donor cell count, capped at 100, floored at 10
    donor_counts <- table(seurat_obj$orig.ident)
    cat("Cells per donor after QC:\n"); print(donor_counts)
    k_weight <- max(10, min(100, floor(min(donor_counts) * 0.9)))
    cat("Dynamic k.weight set to:", k_weight, "\n")
    seurat_obj_integrated <- IntegrateData(anchorset = anchors, dims = 1:5, k.weight = k_weight, verbose = FALSE)
    seurat_obj_integrated <- ScaleData(seurat_obj_integrated, verbose = FALSE)
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)
    seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, dims = 1:8)
    seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.5)
    seurat_obj_integrated <- RunUMAP(seurat_obj_integrated, dims = 1:8, verbose = FALSE, umap.method = "uwot", metric = "cosine")
    plasma_colors <- viridis(n = length(unique(seurat_obj_integrated$orig.ident)), option = "plasma")
    p3 <- DimPlot(seurat_obj_integrated, group.by = "orig.ident", pt.size = 0.5, label = FALSE, repel = TRUE, cols = plasma_colors) +
        ggtitle("UMAP After Integration (by Donor)") + theme(legend.position = "right")
    print(p3)
    # ADDED: UMAP colored by pre-annotated cell type (Allen MTG provides these labels)
    if ("cell_type_designation" %in% colnames(seurat_obj_integrated@meta.data)) {
        p5 <- DimPlot(seurat_obj_integrated, group.by = "cell_type_designation", pt.size = 0.5, label = FALSE, repel = TRUE) +
            ggtitle("UMAP After Integration (by Cell Type)") + theme(legend.position = "right")
        print(p5)
    }
    p4 <- DimPlot(seurat_obj_integrated, group.by = "orig.ident", pt.size = 0.5, label = TRUE, repel = TRUE) +
        ggtitle("UMAP After Integration (by Donor, labeled)") + theme(legend.position = "right")
    print(p4)
    return(seurat_obj_integrated)
}

options(future.globals.maxSize = Inf)

cat("Processing MTG nuclei (batch effects correction by donor)\n")
seurat_obj_integrated <- correct_batch_effects(seurat_obj)


# =========================================
# 5. UMAP Clusters
# MODIFIED: Single integrated object; gene_of_interest = TRPM4
# =========================================
generate_elbow_plot <- function(seurat_obj_integrated, output_prefix) {
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)
    elbow_plot <- ElbowPlot(seurat_obj_integrated, ndims = 50) +
        labs(title = paste("Elbow Plot for", output_prefix), x = "Principal Components", y = "Standard Deviation") +
        theme(plot.title = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
    ggsave(paste0(output_prefix, "_ElbowPlot.pdf"), elbow_plot, width = 6.83, height = 6.41)
    print(elbow_plot)
}

generate_elbow_plot(seurat_obj_integrated, "AllenMTG")

downstream_analyses <- function(seurat_obj_integrated, gene_of_interest, output_prefix, dims = 15) {
    set.seed(10)
    DefaultAssay(seurat_obj_integrated) <- "RNA"
    seurat_obj_integrated <- FindVariableFeatures(seurat_obj_integrated, selection.method = "vst", nfeatures = 2000)
    seurat_obj_integrated <- ScaleData(seurat_obj_integrated, verbose = FALSE)
    seurat_obj_integrated <- RunPCA(seurat_obj_integrated, verbose = FALSE)
    seurat_obj_integrated <- RunUMAP(seurat_obj_integrated, dims = 1:dims, verbose = FALSE)
    set.seed(11); seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated, dims = 1:dims)
    set.seed(12); seurat_obj_integrated <- FindClusters(seurat_obj_integrated, resolution = 0.5)
    # JoinLayers required in Seurat v5: after integration, RNA assay exists as split donor layers
    # (counts.1 ... counts.8); FindAllMarkers cannot run across split layers without joining first
    seurat_obj_integrated <- JoinLayers(seurat_obj_integrated)
    set.seed(13); cluster_markers <- FindAllMarkers(seurat_obj_integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    print(paste("Cluster markers:", nrow(cluster_markers)))
    if (nrow(cluster_markers) == 0) {
        print("Cluster levels:"); print(levels(Idents(seurat_obj_integrated)))
        print("Cells per cluster:"); print(table(Idents(seurat_obj_integrated)))
    }
    umap_data <- as.data.frame(Embeddings(seurat_obj_integrated, "umap")); umap_data$cluster_id <- Idents(seurat_obj_integrated)
    umap_data_mean <- aggregate(. ~ cluster_id, data = umap_data, FUN = mean)
    plasma_func <- colorRampPalette(viridis::viridis(100, direction = -1, option = "plasma")); portion <- 0.8
    n_colors <- round(length(unique(umap_data$cluster_id)) / portion); plasma_colors <- plasma_func(n_colors)
    set.seed(14); umap_plot_with_labels <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +
        geom_point(size = 0.3, alpha = 0.5) + scale_color_manual(values = plasma_colors) +
        geom_text(data = umap_data_mean, aes(label = cluster_id, x = umap_1, y = umap_2), color = "black", size = 3, fontface = "bold", check_overlap = TRUE) +
        theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
        labs(title = "UMAP plot colored by cluster (with labels)", x = "umap_1", y = "umap_2", color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3)))
    print(umap_plot_with_labels)
    set.seed(15); umap_plot_no_labels <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +
        geom_point(size = 0.3, alpha = 0.5) + scale_color_manual(values = plasma_colors) +
        theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
              axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
        labs(title = "UMAP plot colored by cluster (without labels)", x = "umap_1", y = "umap_2", color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3)))
    print(umap_plot_no_labels)
    # ADDED: UMAP colored by pre-annotated Allen cell type designations
    if ("cell_type_designation" %in% colnames(seurat_obj_integrated@meta.data)) {
        umap_data$cell_type <- seurat_obj_integrated@meta.data[rownames(umap_data), "cell_type_designation"]
        set.seed(17); umap_plot_celltype <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(cell_type))) +
            geom_point(size = 0.3, alpha = 0.5) +
            theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
                  axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white"),
                  legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm")) +
            labs(title = "UMAP colored by Allen MTG cell type designation", x = "umap_1", y = "umap_2", color = "Cell Type") +
            guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
        print(umap_plot_celltype)
    }
    # ADDED: UMAP colored by cortical layer
    if ("cortical_layer_label" %in% colnames(seurat_obj_integrated@meta.data)) {
        umap_data$layer <- seurat_obj_integrated@meta.data[rownames(umap_data), "cortical_layer_label"]
        set.seed(18); umap_plot_layer <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = as.factor(layer))) +
            geom_point(size = 0.3, alpha = 0.5) +
            theme(panel.border = element_rect(fill = NA, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"), axis.text.y = element_text(color = "black"),
                  axis.ticks.y = element_line(color = "black"), panel.background = element_rect(fill = "white")) +
            labs(title = "UMAP colored by cortical layer", x = "umap_1", y = "umap_2", color = "Layer") +
            guides(color = guide_legend(override.aes = list(size = 3)))
        print(umap_plot_layer)
    }
    if (gene_of_interest %in% rownames(seurat_obj_integrated)) {
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
    if (nrow(cluster_markers) > 0) {
        top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
    } else {
        top_markers <- data.frame()
        warning("No cluster markers found.")
    }
    write.table(top_markers, file = paste0(output_prefix, "_top_markers_for_each_cluster_vRibo.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    return(list(seurat_obj = seurat_obj_integrated, cluster_markers = cluster_markers, top_markers = top_markers))
}

# Apply downstream analyses (TRPM4 focus, single MTG cohort)
set.seed(42)
mtg_results <- downstream_analyses(seurat_obj_integrated, "TRPM4", "AllenMTG", dims = 20)


# =======================================================================
# 6. Boxplot plus jitter plot of TRPM4 expression across MTG clusters
# MODIFIED: Single cohort - one boxplot for all Louvain clusters
# =======================================================================
expression_values <- data.frame(
    clusters = Idents(mtg_results$seurat_obj),
    TRPM4 = GetAssayData(seurat_obj, slot = "data")["TRPM4", Cells(mtg_results$seurat_obj)]
)

cluster_order <- sort(unique(as.numeric(as.character(expression_values$clusters))))
expression_values$clusters <- factor(expression_values$clusters, levels = as.character(cluster_order))

plot_colors <- viridis(length(cluster_order), end = 1, option = "plasma")
plot_colors <- rev(plot_colors)
named_colors <- setNames(plot_colors, as.character(cluster_order))

p <- ggplot(expression_values, aes(x = clusters, y = TRPM4, fill = clusters)) +
    geom_jitter(aes(color = clusters), alpha = 0.5, width = 0.3, height = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    scale_fill_manual(values = named_colors) +
    scale_color_manual(values = named_colors) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
    ) +
    labs(
        x = "Allen MTG Cluster",
        y = "Expression of TRPM4",
        title = "Expression of TRPM4 by Allen MTG Cluster"
    )
print(p)

# ADDED: Boxplot by pre-annotated Allen cell type designation
if ("cell_type_designation" %in% colnames(mtg_results$seurat_obj@meta.data)) {
    expression_values_ct <- data.frame(
        cell_type = mtg_results$seurat_obj@meta.data[Cells(mtg_results$seurat_obj), "cell_type_designation"],
        TRPM4 = GetAssayData(seurat_obj, slot = "data")["TRPM4", Cells(mtg_results$seurat_obj)]
    )
    n_ct <- length(unique(expression_values_ct$cell_type))
    ct_colors <- viridis(n_ct, option = "plasma")
    p_ct <- ggplot(expression_values_ct, aes(x = cell_type, y = TRPM4, fill = cell_type)) +
        geom_jitter(aes(color = cell_type), alpha = 0.5, width = 0.3, height = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.2) +
        scale_fill_manual(values = ct_colors) +
        scale_color_manual(values = ct_colors) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
            legend.position = "none"
        ) +
        labs(
            x = "Allen MTG Cell Type",
            y = "Expression of TRPM4",
            title = "TRPM4 by Allen MTG Cell Type Designation"
        )
    print(p_ct)
}

# Kruskal-Wallis and Dunn's Tests across MTG clusters
get_expression_data <- function(seurat_obj, original_seurat, gene) {
    expression_values <- data.frame(
        clusters = Idents(seurat_obj),
        expression = GetAssayData(original_seurat, slot = "data")[gene, Cells(seurat_obj)]
    )
    expression_values %>%
        mutate(clusters = factor(clusters, levels = sort(unique(clusters))))
}

mtg_data <- get_expression_data(mtg_results$seurat_obj, seurat_obj, "TRPM4") %>%
    mutate(type = "AllenMTG")

stats_summary <- mtg_data %>%
    group_by(clusters) %>%
    summarise(
        median = median(expression),
        Q1 = quantile(expression, 0.25),
        Q3 = quantile(expression, 0.75),
        IQR = IQR(expression)
    )

kruskal_result <- kruskal_test(expression ~ clusters, data = mtg_data)
dunn_result <- dunn_test(expression ~ clusters, data = mtg_data, p.adjust.method = "BH")

kruskal_df <- as.data.frame(kruskal_result)
dunn_df <- as.data.frame(dunn_result)

write_xlsx(list(
    KruskalWallis = kruskal_df,
    DunnsTest = dunn_df,
    Statistics = stats_summary
), "AllenMTG_TRPM4_expression_analysis.xlsx")

# Save workspace
save.image(file = "AllenMTG_Seurat_processed.RData")

load("AllenMTG_Seurat_processed.RData")

###############################################################################
# Part 2: Spearman-Kendall Dual-Filter Correlation Analysis
# Genes: CARM1P1 (clusters 9 and 3) | PCP4 (clusters 13 and 4)
# Dataset: Allen Human MTG SMART-seq (2018)
#
# NOTE: Per-cluster median/IQR, KW, and Dunn's stats are already computed
#       in Part 1 Section 6 (boxplot block) - not duplicated here.
###############################################################################
library(Seurat)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(writexl)

# ---------------------------------------------------------------------------
# NOTE: Assumes objects in memory from Part 1. If starting fresh:
# load("AllenMTG_Seurat_processed.RData")
#   mtg_results$seurat_obj  - integrated Seurat object with cluster IDs
#   seurat_obj              - original (pre-integration) Seurat object
# ---------------------------------------------------------------------------

# ===========================================================================
# SECTION 1: Filter to expressed genes within a cluster
# ===========================================================================
get_expressed_matrix <- function(integrated_obj, original_obj, cluster_ids,
                                  gene_of_interest, min_cells = 10) {
    cluster_cells <- WhichCells(integrated_obj, idents = cluster_ids)
    cluster_obj   <- subset(original_obj, cells = cluster_cells)
    DefaultAssay(cluster_obj) <- "RNA"

    norm_mat      <- as.matrix(GetAssayData(cluster_obj, slot = "data", assay = "RNA"))
    gene_expr_vec <- norm_mat[gene_of_interest, ]

    # Keep genes expressed (>0) in at least min_cells cells - removes zero-variance genes
    genes_to_keep <- rowSums(is.finite(norm_mat) & norm_mat > 0) >= min_cells
    norm_mat_filt <- norm_mat[genes_to_keep, ]

    cat("  Genes before filter:", nrow(norm_mat),
        "| After (expressed in >=", min_cells, "cells):", nrow(norm_mat_filt),
        "| Cells:", ncol(norm_mat_filt), "\n")

    return(list(mat       = norm_mat_filt,
                gene_vec  = gene_expr_vec,
                cells     = ncol(norm_mat_filt),
                genes     = nrow(norm_mat_filt)))
}

# ===========================================================================
# SECTION 2: Spearman correlation - vectorised, no parallelisation needed
# cor(gene_rank, t(mat)) computes all gene correlations in one call
# ===========================================================================
run_spearman <- function(data_list, gene_of_interest, cluster_id, output_file) {
    mat      <- data_list$mat
    gene_vec <- data_list$gene_vec

    cat("  Running Spearman:", data_list$genes, "genes x", data_list$cells, "cells ... ")

    # Rank gene of interest once; cor() on t(mat) gives all gene correlations
    gene_rank     <- rank(gene_vec)
    spearman_vals <- cor(gene_rank, t(mat), method = "spearman")[1, ]

    spearman_df <- data.frame(
        Gene                 = names(spearman_vals),
        Spearman_Correlation = as.numeric(spearman_vals),
        stringsAsFactors     = FALSE
    ) %>% arrange(desc(Spearman_Correlation))

    cat("done\n")
    cat("  Tier 1 (Spearman >= 0.6):",
        sum(spearman_df$Spearman_Correlation >= 0.6, na.rm = TRUE), "genes\n")
    cat("  Tier 2 (Spearman >= 0.5, < 0.6):",
        sum(spearman_df$Spearman_Correlation >= 0.5 &
            spearman_df$Spearman_Correlation <  0.6, na.rm = TRUE), "genes\n")

    write_xlsx(spearman_df, path = output_file)
    cat("  Saved:", output_file, "\n\n")
    return(invisible(spearman_df))
}

# ===========================================================================
# SECTION 3: Kendall correlation - foreach/doParallel (same pattern as TRPM4 PCa)
# Each foreach iteration = one gene vs gene_of_interest
# ===========================================================================
run_kendall <- function(data_list, gene_of_interest, cluster_id, output_file) {
    mat      <- data_list$mat
    gene_vec <- data_list$gene_vec
    n_genes  <- nrow(mat)

    num_cores <- max(1, detectCores() - 1)
    cat("  Running Kendall:", n_genes, "genes x", data_list$cells,
        "cells,", num_cores, "cores\n")

    cl <- makeCluster(num_cores)
    registerDoParallel(cl)

    # foreach: one gene per iteration, distributed across workers
    # .export ensures gene_vec and mat are available in each worker environment
    kendall_vals <- foreach(i = seq_len(n_genes),
                            .combine   = 'c',
                            .export    = c("gene_vec", "mat")) %dopar% {
        cor(gene_vec, mat[i, ], method = "kendall")
    }
    stopCluster(cl)

    kendall_df <- data.frame(
        Gene                = rownames(mat),
        Kendall_Correlation = kendall_vals,
        stringsAsFactors    = FALSE
    ) %>% arrange(desc(Kendall_Correlation))

    cat("  Tier 1 (Kendall >= 0.5):",
        sum(kendall_df$Kendall_Correlation >= 0.5, na.rm = TRUE), "genes\n")
    cat("  Tier 2 (Kendall >= 0.4, < 0.5):",
        sum(kendall_df$Kendall_Correlation >= 0.4 &
            kendall_df$Kendall_Correlation <  0.5, na.rm = TRUE), "genes\n")

    write_xlsx(kendall_df, path = output_file)
    cat("  Saved:", output_file, "\n\n")
    return(invisible(kendall_df))
}

# ===========================================================================
# SECTION 4: CARM1P1 - cluster 9 (primary) and cluster 3 (secondary)
# ===========================================================================
cat("=================================================================\n")
cat("CARM1P1 DUAL-FILTER CORRELATIONS\n")
cat("=================================================================\n")

cat("\n-- CARM1P1 | Cluster 9 (Exc L3-5 RORB eC1 -- CARM1P1/COL22A1 subtype) --\n")
d_carm_9  <- get_expressed_matrix(mtg_results$seurat_obj, seurat_obj, 9, "CARM1P1")
sp_carm_9 <- run_spearman(d_carm_9, "CARM1P1", 9, "CARM1P1_Cluster9_spearman.xlsx")
kd_carm_9 <- run_kendall (d_carm_9, "CARM1P1", 9, "CARM1P1_Cluster9_kendall.xlsx")

cat("\n-- CARM1P1 | Cluster 3 (Exc L3 RORB CARTPT -- secondary non-zero mean) --\n")
d_carm_3  <- get_expressed_matrix(mtg_results$seurat_obj, seurat_obj, 3, "CARM1P1")
sp_carm_3 <- run_spearman(d_carm_3, "CARM1P1", 3, "CARM1P1_Cluster3_spearman.xlsx")
kd_carm_3 <- run_kendall (d_carm_3, "CARM1P1", 3, "CARM1P1_Cluster3_kendall.xlsx")


# ===========================================================================
# SECTION 5: PCP4 - cluster 13 (primary) and cluster 4 (secondary)
# ===========================================================================
cat("=================================================================\n")
cat("PCP4 DUAL-FILTER CORRELATIONS\n")
cat("=================================================================\n")

cat("\n-- PCP4 | Cluster 13 (Exc L4-6 FEZF2 eD Near-Projecting -- highest PCP4) --\n")
d_pcp4_13  <- get_expressed_matrix(mtg_results$seurat_obj, seurat_obj, 13, "PCP4")
sp_pcp4_13 <- run_spearman(d_pcp4_13, "PCP4", 13, "PCP4_Cluster13_spearman.xlsx")
kd_pcp4_13 <- run_kendall (d_pcp4_13, "PCP4", 13, "PCP4_Cluster13_kendall.xlsx")

cat("\n-- PCP4 | Cluster 4 (Exc L4-5 RORB FOLH1B -- secondary non-zero mean) --\n")
d_pcp4_4  <- get_expressed_matrix(mtg_results$seurat_obj, seurat_obj, 4, "PCP4")
sp_pcp4_4 <- run_spearman(d_pcp4_4, "PCP4", 4, "PCP4_Cluster4_spearman.xlsx")
kd_pcp4_4 <- run_kendall (d_pcp4_4, "PCP4", 4, "PCP4_Cluster4_kendall.xlsx")

cat("=================================================================\n")
cat("Part 2 complete. Output files:\n")
cat("  CARM1P1_Cluster9_spearman.xlsx\n")
cat("  CARM1P1_Cluster9_kendall.xlsx\n")
cat("  CARM1P1_Cluster3_spearman.xlsx\n")
cat("  CARM1P1_Cluster3_kendall.xlsx\n")
cat("  PCP4_Cluster13_spearman.xlsx\n")
cat("  PCP4_Cluster13_kendall.xlsx\n")
cat("  PCP4_Cluster4_spearman.xlsx\n")
cat("  PCP4_Cluster4_kendall.xlsx\n")
cat("=================================================================\n")



# =========================================
# Heatmap of CARM1P1 vs Axon Guidance and Housekeeping Genes
# Allen Human MTG SMART-seq (2018)
# =========================================
library(ComplexHeatmap)

# Subset data (specific clusters for heatmap)
# MODIFIED: single object, 8 MTG clusters; no merge needed
heatmap_obj <- subset(mtg_results$seurat_obj, idents = c(0, 3, 9, 15))

# Define gene sets
# MODIFIED: CARM1P1 correlated axon guidance genes + housekeeping genes
genes_of_interest <- c("CARM1P1", "CLMN", "EPHA3", "EPHA6", "LOC101928964", "ROBO2",
                       "ACTB", "GAPDH", "RPL13A")
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(heatmap_obj)]

# MODIFIED: extract CARM1P1 expression for sorting (replaces trpm4_expression)
carm1p1_expression <- GetAssayData(seurat_obj, slot = "data")["CARM1P1",
                          Cells(heatmap_obj)]

# Prepare expression matrix (log2-transformed)
# MODIFIED: no averaged composite rows (Ribo/AR) - all genes shown individually
expression_matrix     <- GetAssayData(seurat_obj, slot = "data")[genes_of_interest,
                             Cells(heatmap_obj)]
expression_matrix_log <- log2(expression_matrix + 1)

# Scale to z-scores (cap at -1 to 2.3)
scale_to_zscore <- function(x) {
    z <- scale(x)
    return(pmin(pmax(z, -1), 2.3))
}
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Order cells by cluster and CARM1P1 expression (ascending)
# MODIFIED: cluster levels = 8 MTG clusters; sort key = carm1p1_expression
heatmap_grouping <- Idents(heatmap_obj)
heatmap_grouping <- factor(heatmap_grouping,
                            levels = c("0", "3", "9", "15"))
cell_order <- order(heatmap_grouping, carm1p1_expression)
expression_matrix_log  <- expression_matrix_log[, cell_order]
heatmap_grouping       <- heatmap_grouping[cell_order]
carm1p1_expression     <- carm1p1_expression[cell_order]
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Prepare heatmap annotations (cluster labels with counts)
n_groups         <- length(levels(heatmap_grouping))
group_colors     <- viridis(n_groups, option = "mako"); names(group_colors) <- levels(heatmap_grouping)
cell_counts      <- table(heatmap_grouping); cell_percentages <- round(prop.table(cell_counts) * 100, 1)
group_labels     <- paste0(names(cell_counts), "\n", cell_counts, "\n(", cell_percentages, "%)")
column_annotation <- HeatmapAnnotation(
    Cluster = anno_block(gp = gpar(fill = group_colors), labels = group_labels,
                         labels_gp = gpar(col = "white", fontsize = 8)),
    show_legend = FALSE, show_annotation_name = FALSE)

# Create heatmap
# MODIFIED: output filename; width adjusted for 8 clusters
heatmap <- Heatmap(expression_matrix_scaled, name = "zscore", column_title = NULL, row_title = "Genes",
                   show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                   top_annotation = column_annotation, col = viridis(100, option = "mako"), row_names_gp = gpar(fontsize = 8),
                   column_split = heatmap_grouping, row_gap = unit(1, "mm"), column_gap = unit(0.5, "mm"), border = TRUE,
                   use_raster = TRUE, raster_quality = 2)

# Save and display heatmap
# MODIFIED: filename
pdf("AllenMTG_CARM1P1_heatmap.pdf", width = 14, height = 6)
draw(heatmap)
dev.off()
draw(heatmap)


##################################################################
# Allen MTG CARM1P1 GAM-REML-PRSS Analysis - Part 1
# x-axis: CARM1P1 (log2-transformed)
# y-axis: DFG = avg(CLMN, EPHA3, EPHA6, LOC101928964, ROBO2)
#         HKG = avg(ACTB, GAPDH, PPIA, TBP, RPL13A)
# Cluster 9 cells (Exc L3-5 RORB eC1), run separately per donor
# 8 donors x 2 gene sets = 16 GAM fits
# gamma = 1.5, method = REML, select = TRUE, bs = "tp"
#
# Output: AllenMTG_CARM1P1_GAM_Results_1.xlsx  (8 sheets)
##################################################################
setwd("C:/Users/Wong/Desktop/GO_preIP_&_postIP/Brain_datasets")
load("AllenMTG_Seurat_processed.RData")
library(mgcv)
library(Seurat)
library(dplyr)
library(writexl)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Per-donor data preparation for cluster 9
# ------------------------------------------------------------------
donors <- c("H200.1030", "H16.06.002", "H200.1025", "H200.1023",
            "H16.24.010", "H16.06.009", "H16.03.004", "H16.06.008")

dfg_genes <- c("CLMN", "EPHA3", "EPHA6", "LOC101928964", "ROBO2")
hkg_genes <- c("ACTB", "GAPDH", "TBP", "PPIA", "RPL13A")

cluster9_cells <- WhichCells(mtg_results$seurat_obj, idents = 9)

cat("Building per-donor data for cluster 9...\n")
donor_data <- lapply(donors, function(donor) {
    donor_cells <- colnames(seurat_obj)[seurat_obj$orig.ident == donor]
    cells       <- intersect(cluster9_cells, donor_cells)
    if (length(cells) < 10) {
        cat("  Donor", donor, ": only", length(cells),
            "cluster-9 cells -- skipping (min 10 required)\n")
        return(NULL)
    }
    obj <- subset(seurat_obj, cells = cells)
    DefaultAssay(obj) <- "RNA"
    mat <- as.matrix(GetAssayData(obj, slot = "data", assay = "RNA"))
    dfg_present <- dfg_genes[dfg_genes %in% rownames(mat)]
    hkg_present <- hkg_genes[hkg_genes %in% rownames(mat)]
    carm1p1 <- log2(mat["CARM1P1", ] + 1)
    dfg_avg  <- log2(colMeans(mat[dfg_present, , drop = FALSE]) + 1)
    hkg_avg  <- log2(colMeans(mat[hkg_present, , drop = FALSE]) + 1)
    cat("  Donor", donor, "| cells:", length(cells),
        "| unique CARM1P1:", length(unique(carm1p1)),
        "| DFG genes:", length(dfg_present),
        "| HKG genes:", length(hkg_present), "\n")
    list(
        dfg     = data.frame(CARM1P1 = carm1p1, Expression = dfg_avg),
        hkg     = data.frame(CARM1P1 = carm1p1, Expression = hkg_avg),
        carm1p1 = carm1p1,
        mat     = mat,
        dfg_present = dfg_present,
        hkg_present = hkg_present
    )
})
names(donor_data) <- donors
cat("Donors with sufficient cells:",
    sum(!sapply(donor_data, is.null)), "/", length(donors), "\n")

# ------------------------------------------------------------------
# SECTION 2: Helper functions
# ------------------------------------------------------------------

calculate_prss <- function(model, data) {
    RSS    <- sum(residuals(model)^2)
    lambda <- model$sp[1]
    S      <- model$smooth[[1]]$S[[1]]
    beta   <- coef(model)[model$smooth[[1]]$first.para:model$smooth[[1]]$last.para]
    f_dpi  <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

extract_equation <- function(model) {
    pt <- coef(model)[!grepl("s\\(CARM1P1\\)", names(coef(model)))]
    st <- coef(model)[ grepl("s\\(CARM1P1\\)", names(coef(model)))]
    lhs <- sprintf("%.9f", pt[1])
    if (length(pt) > 1)
        lhs <- paste0(lhs, " + ", sprintf("%.9f*x", pt[2]))
    rhs <- paste(mapply(function(v, i) sprintf("%.9f*\u03c6%d(x)", v, i),
                        st, seq_along(st)), collapse = " + ")
    paste0("f(x) = ", lhs, if (length(st) > 0) " + " else "", rhs)
}

extract_reml_iterations_enhanced <- function(model) {
    if (is.null(model) || is.null(model$outer.info))
        return(data.frame(iteration = numeric(0), lambda = character(0),
                          score = numeric(0), stringsAsFactors = FALSE))
    tryCatch({
        n  <- model$outer.info$iter
        df <- data.frame(iteration = seq_len(n), lambda = NA_character_,
                         score = NA_real_, stringsAsFactors = FALSE)
        if (!is.null(model$outer.info$sp)) {
            sp <- model$outer.info$sp
            if (is.matrix(sp)) {
                for (i in seq_len(min(nrow(sp), n)))
                    df$lambda[i] <- paste(format(sp[i,], digits=6), collapse=", ")
            } else if (is.vector(sp) && length(sp) >= n) {
                df$lambda <- format(sp[seq_len(n)], digits=6)
            } else if (length(sp) > 0) {
                df$lambda <- format(sp[1], digits=6)
            }
        }
        if (all(is.na(df$lambda)) && !is.null(model$sp))
            df$lambda <- paste(format(model$sp, digits=6), collapse=", ")
        if (!is.null(model$outer.info$score))
            df$score[seq_len(min(length(model$outer.info$score), n))] <-
                model$outer.info$score[seq_len(min(length(model$outer.info$score), n))]
        df
    }, error = function(e) {
        cat("  Warning - extract_reml_iterations:", conditionMessage(e), "\n")
        data.frame(iteration = numeric(0), lambda = character(0),
                   score = numeric(0), stringsAsFactors = FALSE)
    })
}

extract_convergence_details <- function(model) {
    blank <- list(convergence_criterion="Unknown", gradient_value=NA_character_,
                  gradient_norm=NA_real_, relative_score_change=NA_real_,
                  iterations=NA_integer_, max_iterations_reached=NA,
                  final_reml_score=NA_real_, optimizer=NA_character_,
                  monotonic_decrease="NA")
    if (is.null(model) || is.null(model$outer.info)) return(blank)
    tryCatch({
        iters <- model$outer.info$iter
        gv <- NA_character_; gn <- NA_real_
        if (!is.null(model$outer.info$grad)) {
            g <- model$outer.info$grad
            if (is.vector(g) && length(g) > 0) {
                fg <- g[length(g)]; gv <- format(fg, digits=8); gn <- abs(fg)
            } else if (is.matrix(g) && nrow(g) > 0) {
                fg <- g[nrow(g),]
                gv <- paste(format(fg, digits=8), collapse=", ")
                gn <- sqrt(sum(fg^2))
            }
        }
        rsc <- NA_real_
        if (!is.null(model$outer.info$score) && length(model$outer.info$score) > 1) {
            sc  <- model$outer.info$score
            rsc <- tail(abs(diff(sc)) / abs(sc[-length(sc)]), 1)
        }
        mono <- "NO"
        if (!is.null(model$outer.info$score) && length(model$outer.info$score) > 1)
            mono <- if (all(diff(model$outer.info$score) <= 1e-10)) "YES" else "NO"
        crit <- "Unknown"
        if (!is.na(iters) && !is.null(model$control$maxit) &&
            iters >= model$control$maxit) {
            crit <- "Maximum number of iterations reached"
        } else if (!is.na(gn) && gn < 1e-5) {
            crit <- "Gradient-based convergence"
        } else if (!is.na(rsc) && rsc < 1e-6) {
            crit <- "Score-based convergence"
        } else if (isTRUE(model$converged)) {
            crit <- "Fallback convergence"
        }
        list(convergence_criterion  = crit,
             gradient_value         = gv,
             gradient_norm          = gn,
             relative_score_change  = rsc,
             iterations             = iters,
             max_iterations_reached = if (!is.na(iters) && !is.null(model$control$maxit))
                                          iters >= model$control$maxit else NA,
             final_reml_score       = model$gcv.ubre,
             optimizer              = paste(model$optimizer, collapse=", "),
             monotonic_decrease     = mono)
    }, error = function(e) {
        cat("  Warning - extract_convergence:", conditionMessage(e), "\n")
        blank
    })
}

# ------------------------------------------------------------------
# SECTION 3: GAM fitting with two-phase PRSS optimisation
# ------------------------------------------------------------------
fit_gam_prss <- function(data, donor_id, gene_set_label, num_iterations = 100) {
    set.seed(123)
    unique_x <- length(unique(data$CARM1P1))
    max_k    <- min(10, unique_x - 1)
    cat("  [", donor_id, "|", gene_set_label, "]",
        "cells:", nrow(data), "| unique CARM1P1:", unique_x,
        "| max k:", max_k, "\n")

    best_prss <- Inf; best_model <- NULL; best_k <- NULL
    best_iter <- 0;   prss_values <- numeric(num_iterations)
    model_params <- list(); reml_iters_list <- list()
    early_stop_counter <- 0; early_stop_triggered <- FALSE

    run_one <- function(i, k) {
        set.seed(123 + i)
        gam_model <- gam(Expression ~ CARM1P1 + s(CARM1P1, bs = "tp", k = k),
                         data = data, method = "REML", select = TRUE, gamma = 1.5)
        pc <- calculate_prss(gam_model, data)
        ri <- extract_reml_iterations_enhanced(gam_model)
        if (nrow(ri) > 0) {
            ri$prss_iteration <- i; ri$k_value <- k
            ri$prss_score <- pc$PRSS; ri$is_best_model <- FALSE
            reml_iters_list[[i]] <<- ri
        }
        prss_values[i] <<- pc$PRSS
        reml_scores_i <- if (!is.null(gam_model$outer.info$score))
                             gam_model$outer.info$score else gam_model$gcv.ubre
        reml_n_i   <- if (!is.null(gam_model$outer.info$iter)) gam_model$outer.info$iter else 1L
        reml_min_i <- min(reml_scores_i, na.rm = TRUE)
        reml_max_i <- max(reml_scores_i, na.rm = TRUE)
        model_params[[i]] <<- list(
            iteration = i, k = k, lambda = pc$lambda,
            PRSS = pc$PRSS, RSS = pc$RSS,
            f_dpi = pc$f_double_prime_integral,
            dev_explained = summary(gam_model)$dev.expl,
            reml_n    = reml_n_i,
            reml_min  = reml_min_i,
            reml_max  = reml_max_i,
            smooth_edf = summary(gam_model)$edf,
            p_value = summary(gam_model)$s.table[4])
        list(model = gam_model, prss = pc$PRSS, k = k)
    }

    # ---- Phase 1: exploration (iterations 1-8) ----
    # k cycles through 3, 4, ..., max_k  (wraps via modulo)
    for (i in seq_len(min(8, num_iterations))) {
        k <- 3 + ((i - 1) %% (max_k - 2))
        tryCatch({
            res <- run_one(i, k)
            cat("    init iter", i, "| k =", k,
                "| PRSS:", round(res$prss, 4), "\n")
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_stop_counter <- 0
            } else { early_stop_counter <- early_stop_counter + 1 }
        }, error = function(e) {
            cat("    error iter", i, ":", conditionMessage(e), "\n")
            prss_values[i] <<- NA
        })
    }

    # ---- Set initial k_range from phase 1 best ----
    k_range <- if (best_k == max_k) c(max_k - 1, max_k) else
               if (best_k == 3)    c(3, 4) else
               c(best_k - 1, best_k, best_k + 1)
    cat("    best k after phase 1:", best_k,
        "| PRSS:", round(best_prss, 4), "\n")

    # ---- Phase 2: adaptive refinement (iterations 9-100) ----
    # Early stop after 20 consecutive iterations without improvement
    for (i in 9:num_iterations) {
        if (early_stop_counter >= 20) {
            cat("    early stopping at iter", i, "\n")
            early_stop_triggered <- TRUE; break
        }
        # Every 3rd iteration: sample k from best_k ± 1 range
        k <- if (i %% 3 == 0) sample(k_range, 1) else best_k
        tryCatch({
            res <- run_one(i, k)
            cat("    refine iter", i, "| k =", k,
                "| PRSS:", round(res$prss, 4), "\n")
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_stop_counter <- 0
                k_range <- if (best_k == max_k) c(max_k - 1, max_k) else
                           if (best_k == 3)     c(3, 4) else
                           c(best_k - 1, best_k, best_k + 1)
            } else { early_stop_counter <- early_stop_counter + 1 }
        }, error = function(e) {
            cat("    error iter", i, ":", conditionMessage(e), "\n")
            prss_values[i] <<- NA
            early_stop_counter <<- early_stop_counter + 1
        })
    }

    actual_iters <- max(which(!is.na(prss_values)))

    valid_reml <- Filter(Negate(is.null), reml_iters_list)
    all_reml   <- if (length(valid_reml) > 0) do.call(rbind, valid_reml) else NULL
    if (!is.null(all_reml) && nrow(all_reml) > 0)
        all_reml$is_best_model[all_reml$prss_iteration == best_iter] <- TRUE

    list(best_model    = best_model,
         best_iter     = best_iter,
         best_k        = best_k,
         best_prss     = best_prss,
         prss_values   = prss_values[seq_len(actual_iters)],
         model_params  = model_params[seq_len(actual_iters)],
         all_reml_iters = all_reml,
         early_stopped  = early_stop_triggered,
         actual_iters   = actual_iters,
         donor          = donor_id,
         gene_set       = gene_set_label,
         data           = data)
}

# ------------------------------------------------------------------
# SECTION 4: Run 16 fits (8 donors x 2 gene sets)
# ------------------------------------------------------------------
cat("\n========== Running GAM fits: 8 donors x 2 gene sets ==========\n")

all_results <- list()

for (donor in donors) {
    dd <- donor_data[[donor]]
    if (is.null(dd)) {
        cat("Skipping donor", donor, "(insufficient cells)\n")
        next
    }
    cat("\n--- Donor:", donor, "---\n")

    cat(" Fitting DFG...\n")
    key_dfg <- paste0(donor, "|DFG")
    all_results[[key_dfg]] <- tryCatch(
        fit_gam_prss(dd$dfg, donor, "DFG"),
        error = function(e) {
            cat("  FAILED:", conditionMessage(e), "\n"); NULL })

    cat(" Fitting HKG...\n")
    key_hkg <- paste0(donor, "|HKG")
    all_results[[key_hkg]] <- tryCatch(
        fit_gam_prss(dd$hkg, donor, "HKG"),
        error = function(e) {
            cat("  FAILED:", conditionMessage(e), "\n"); NULL })
}

n_success <- sum(!sapply(all_results, is.null))
cat("\nCompleted:", n_success, "/ 16 fits\n")

# ------------------------------------------------------------------
# SECTION 5: Sheet-building helpers
# ------------------------------------------------------------------

# ---- 5a. Convergence row (one per fit) ----
build_convergence_row <- function(result) {
    cv <- extract_convergence_details(result$best_model)
    data.frame(
        Donor                   = result$donor,
        Group                   = "Allen_MTG",
        Gene_Set                = result$gene_set,
        REML_monotonic_decrease = cv$monotonic_decrease,
        Convergence_criterion   = cv$convergence_criterion,
        Gradient_value          = cv$gradient_value,
        Gradient_norm           = cv$gradient_norm,
        Relative_score_change   = cv$relative_score_change,
        Iterations              = cv$iterations,
        Max_iterations_reached  = cv$max_iterations_reached,
        Final_REML_score        = cv$final_reml_score,
        Optimizer               = cv$optimizer,
        stringsAsFactors = FALSE)
}

# ---- 5b. PRSS iterations (all iterations per fit) ----
build_prss_rows <- function(result) {
    n     <- result$actual_iters
    iters <- seq_len(n)
    get_p <- function(i, key) {
        x <- result$model_params[[i]]
        if (is.null(x)) NA else x[[key]]
    }
    data.frame(
        Donor               = result$donor,
        Group               = "Allen_MTG",
        Gene_Set            = result$gene_set,
        Best_Model_Tag      = iters == result$best_iter,
        k_Value             = sapply(iters, function(i) get_p(i, "k")),
        Lambda_Value        = sapply(iters, function(i) get_p(i, "lambda")),
        PRSS_Iteration      = iters,
        PRSS_Value          = result$prss_values,
        REML_num_iterations = sapply(iters, function(i) get_p(i, "reml_n")),
        REML_min_score      = sapply(iters, function(i) get_p(i, "reml_min")),
        REML_max_score      = sapply(iters, function(i) get_p(i, "reml_max")),
        stringsAsFactors = FALSE)
}

# ---- 5c. REML iterations of best PRSS model ----
build_reml_best_rows <- function(result) {
    all_ri <- result$all_reml_iters
    if (is.null(all_ri) || nrow(all_ri) == 0) return(NULL)
    best_ri <- all_ri[all_ri$is_best_model == TRUE, ]
    if (nrow(best_ri) == 0) return(NULL)
    m <- result$best_model
    hess_str <- NA_character_
    if (!is.null(m$outer.info$hess)) {
        h <- m$outer.info$hess
        if (is.matrix(h))
            hess_str <- paste(apply(h, 1, function(r)
                paste(format(r, digits=4, scientific=TRUE), collapse=" ")),
                collapse="; ")
    }
    nr <- nrow(best_ri)
    grad_vec <- rep(NA_real_, nr)
    if (!is.null(m$outer.info$grad)) {
        g <- m$outer.info$grad
        if (is.matrix(g) && nrow(g) >= nr) {
            grad_vec <- g[seq_len(nr), 1]
        } else if (is.vector(g) && length(g) >= nr) {
            grad_vec <- g[seq_len(nr)]
        } else if (length(g) > 0) {
            grad_vec[seq_len(min(length(g), nr))] <-
                as.numeric(g)[seq_len(min(length(g), nr))]
        }
    }
    data.frame(
        Donor          = result$donor,
        Group          = "Allen_MTG",
        Gene_Set       = result$gene_set,
        Best_Model     = TRUE,
        REML_Iteration = best_ri$iteration,
        REML_Score     = best_ri$score,
        Lambda_Value   = best_ri$lambda,
        Gradient       = grad_vec,
        Hessian        = c(rep(NA_character_, nr - 1), hess_str),
        Hessian_Notes  = c(rep("Not calculated for intermediate steps", nr - 1),
                           "Final Hessian at convergence"),
        stringsAsFactors = FALSE)
}

# ---- 5d. GAM performance row - original columns + new EDF/basis columns ----
build_gam_row <- function(result) {
    m  <- result$best_model
    pc <- calculate_prss(m, result$data)
    sm <- summary(m)

    # ---- Original columns ----
    p_val <- sm$s.table[4]

    # ---- New: EDF components ----
    smooth_names  <- grepl("s\\(CARM1P1\\)", names(coef(m)))
    parametric_df <- sum(!smooth_names)
    total_edf     <- sum(sm$edf)    # intercept + parametric + smooth

    # ---- New: basis function counts ----
    smooth_coefs         <- coef(m)[smooth_names]
    nonzero_basis        <- sum(abs(smooth_coefs) > 1e-8)
    total_basis          <- length(smooth_coefs)
    smooth_coef_count    <- length(smooth_coefs)

    data.frame(
        Donor                   = result$donor,
        Group                   = "Allen_MTG",
        Gene_Set                = result$gene_set,
        PRSS                    = pc$PRSS,
        RSS                     = pc$RSS,
        Lambda_REML             = pc$lambda,
        f_double_prime_int      = pc$f_double_prime_integral,
        k                       = result$best_k,
        Model_deviance          = deviance(m),
        Null_deviance           = deviance(gam(Expression ~ 1, data = result$data)),
        DE                      = sm$dev.expl,
        Equation                = extract_equation(m),
        Smooth_EDF              = sm$edf[length(sm$edf)],
        Smooth_p_value          = p_val,
        # BH FDR added after pooling - placeholder NA, filled in Section 6
        Smooth_p_BH_FDR         = NA_real_,
        # New EDF columns
        Total_EDF               = total_edf,
        Parametric_DF           = parametric_df,
        # Nonlinearity flag filled in Section 6 after FDR is known
        Is_Significantly_Nonlinear = NA,
        # New basis function columns
        Nonzero_Basis_Functions = nonzero_basis,
        Total_Basis_Functions   = total_basis,
        Smooth_Coef_Count       = smooth_coef_count,
        stringsAsFactors = FALSE)
}

# ---- 5e. Observed / Predicted ----
build_obs_pred_rows <- function(result) {
    data.frame(
        Donor           = result$donor,
        Gene_Set        = result$gene_set,
        CARM1P1         = result$data$CARM1P1,
        Observed_value  = result$data$Expression,
        Predicted_value = as.numeric(predict(result$best_model,
                                             newdata = result$data)),
        stringsAsFactors = FALSE)
}

# ---- 5f. Basis functions ----
build_basis_rows <- function(result) {
    m           <- result$best_model
    X           <- predict(m, type = "lpmatrix")
    smooth_cols <- grep("s\\(CARM1P1\\)", colnames(X))
    bm          <- X[, smooth_cols, drop = FALSE]
    colnames(bm) <- paste0("phi_", seq_len(ncol(bm)))
    data.frame(Donor    = result$donor,
               Gene_Set = result$gene_set,
               CARM1P1  = result$data$CARM1P1,
               as.data.frame(bm), stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------
# SECTION 6: Collect all results into sheets
# ------------------------------------------------------------------
valid_results <- Filter(Negate(is.null), all_results)

sheet_conv      <- do.call(rbind, lapply(valid_results, build_convergence_row))
sheet_prss      <- do.call(rbind, lapply(valid_results, build_prss_rows))
sheet_reml_best <- do.call(rbind, lapply(valid_results, build_reml_best_rows))

# Split by gene set for separate performance sheets
dfg_results <- Filter(function(r) r$gene_set == "DFG", valid_results)
hkg_results <- Filter(function(r) r$gene_set == "HKG", valid_results)

sheet_perf_dfg <- do.call(rbind, lapply(dfg_results, build_gam_row))
sheet_perf_hkg <- do.call(rbind, lapply(hkg_results, build_gam_row))

# Apply cross-donor BH FDR per gene set, then compute nonlinearity flag
if (!is.null(sheet_perf_dfg) && nrow(sheet_perf_dfg) > 0) {
    sheet_perf_dfg$Smooth_p_BH_FDR <- p.adjust(sheet_perf_dfg$Smooth_p_value,
                                                 method = "BH")
    sheet_perf_dfg$Is_Significantly_Nonlinear <-
        sheet_perf_dfg$Smooth_EDF > 1 & sheet_perf_dfg$Smooth_p_BH_FDR < 0.05
}
if (!is.null(sheet_perf_hkg) && nrow(sheet_perf_hkg) > 0) {
    sheet_perf_hkg$Smooth_p_BH_FDR <- p.adjust(sheet_perf_hkg$Smooth_p_value,
                                                 method = "BH")
    sheet_perf_hkg$Is_Significantly_Nonlinear <-
        sheet_perf_hkg$Smooth_EDF > 1 & sheet_perf_hkg$Smooth_p_BH_FDR < 0.05
}

sheet_obs_pred  <- do.call(rbind, lapply(valid_results, build_obs_pred_rows))

# Basis functions: pad phi columns to uniform width before rbind
basis_list <- lapply(valid_results, build_basis_rows)
max_phi    <- max(sapply(basis_list, function(df) sum(grepl("^phi_", names(df)))))
pad_phi <- function(df) {
    ex <- sum(grepl("^phi_", names(df)))
    if (ex < max_phi)
        for (j in (ex + 1):max_phi)
            df[[paste0("phi_", j)]] <- NA_real_
    df
}
sheet_basis <- do.call(rbind, lapply(basis_list, pad_phi))

# ------------------------------------------------------------------
# SECTION 6b: Leave-one-out gene contribution analysis
# ------------------------------------------------------------------
fit_gam_loo <- function(data, best_k) {
    set.seed(123)
    k_use <- max(3L, min(best_k, length(unique(data$CARM1P1)) - 1L))
    tryCatch({
        m <- gam(Expression ~ CARM1P1 + s(CARM1P1, bs = "tp", k = k_use),
                 data = data, method = "REML", select = TRUE, gamma = 1.5)
        summary(m)$dev.expl
    }, error = function(e) NA_real_)
}

build_gene_contribution_sheet <- function(donor, dd, full_results) {
    rows <- list()
    for (gs in c("DFG", "HKG")) {
        key      <- paste0(donor, "|", gs)
        full_res <- full_results[[key]]
        if (is.null(full_res) || is.null(full_res$best_model)) next
        de_full <- summary(full_res$best_model)$dev.expl
        best_k  <- full_res$best_k
        genes   <- if (gs == "DFG") dd$dfg_present else dd$hkg_present
        mat     <- dd$mat
        carm1p1 <- dd$carm1p1
        if (length(genes) < 2) {
            rows[[length(rows) + 1]] <- data.frame(
                Donor            = donor,
                Group            = "Allen_MTG",
                Gene_Set         = gs,
                Removed_Gene     = "(only 1 gene present - LOO not applicable)",
                DE_Full          = de_full,
                DE_LOO           = NA_real_,
                Delta_DE         = NA_real_,
                Pct_Contribution = NA_real_,
                stringsAsFactors = FALSE)
            next
        }
        for (gene_out in genes) {
            genes_in <- setdiff(genes, gene_out)
            loo_avg  <- log2(colMeans(mat[genes_in, , drop = FALSE]) + 1)
            loo_data <- data.frame(CARM1P1 = carm1p1, Expression = loo_avg)
            de_loo   <- fit_gam_loo(loo_data, best_k)
            delta    <- de_full - de_loo
            pct_contrib <- if (!is.na(de_full) && de_full > 0)
                               100 * delta / de_full else NA_real_
            cat("    LOO [", donor, "|", gs, "] remove", gene_out,
                "| DE_full:", round(de_full, 4),
                "| DE_loo:", round(de_loo, 4),
                "| delta:", round(delta, 4), "\n")
            rows[[length(rows) + 1]] <- data.frame(
                Donor            = donor,
                Group            = "Allen_MTG",
                Gene_Set         = gs,
                Removed_Gene     = gene_out,
                DE_Full          = de_full,
                DE_LOO           = de_loo,
                Delta_DE         = delta,
                Pct_Contribution = pct_contrib,
                stringsAsFactors = FALSE)
        }
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

cat("\n========== Leave-one-out gene contribution analysis ==========\n")
cat("Fitting", length(donors), "donors x 2 gene sets x 5 genes = up to 80 LOO GAMs\n\n")

loo_list <- lapply(donors, function(donor) {
    dd <- donor_data[[donor]]
    if (is.null(dd)) return(NULL)
    cat("  LOO donor:", donor, "\n")
    tryCatch(
        build_gene_contribution_sheet(donor, dd, all_results),
        error = function(e) {
            cat("  LOO FAILED for", donor, ":", conditionMessage(e), "\n")
            NULL
        })
})
sheet_loo <- do.call(rbind, Filter(Negate(is.null), loo_list))
cat("\nLOO sheet rows:", nrow(sheet_loo), "\n")

# ------------------------------------------------------------------
# SECTION 7: Export
# ------------------------------------------------------------------
export_list <- list(
    REML_Convergence          = sheet_conv,
    All_Models_PRSS           = sheet_prss,
    REML_Iterations_BestModel = sheet_reml_best,
    GAM_Performance_DFG       = sheet_perf_dfg,
    GAM_Performance_HKG       = sheet_perf_hkg,
    Observed_Predicted        = sheet_obs_pred,
    Basis_Functions           = sheet_basis,
    Gene_Contribution_LOO     = sheet_loo
)

write_xlsx(export_list, path = "AllenMTG_CARM1P1_GAM_Results_1.xlsx")
cat("\nSaved: AllenMTG_CARM1P1_GAM_Results_1.xlsx\n")
cat("Sheets:\n")
cat("  REML_Convergence          ", nrow(sheet_conv),      "rows\n")
cat("  All_Models_PRSS           ", nrow(sheet_prss),      "rows\n")
cat("  REML_Iterations_BestModel ", nrow(sheet_reml_best), "rows\n")
cat("  GAM_Performance_DFG       ", nrow(sheet_perf_dfg),
    "rows | columns: PRSS, RSS, Lambda, f'', k, deviance, DE, Eq, EDF,",
    "p, BH-FDR, Total_EDF, Param_DF, Nonlinear, NZ-basis, Tot-basis, Coef-count\n")
cat("  GAM_Performance_HKG       ", nrow(sheet_perf_hkg),  "rows (same columns as DFG)\n")
cat("  Observed_Predicted        ", nrow(sheet_obs_pred),  "rows\n")
cat("  Basis_Functions           ", nrow(sheet_basis),     "rows\n")
cat("  Gene_Contribution_LOO     ", nrow(sheet_loo),
    "rows (8 donors x 2 gene sets x 5 genes)\n")

# Save model objects for Part 2 (plots)
saveRDS(list(all_results = all_results,
             donor_data  = donor_data,
             donors      = donors,
             sheet_loo   = sheet_loo),
        "AllenMTG_CARM1P1_GAM_models.rds")
cat("Models saved to: AllenMTG_CARM1P1_GAM_models.rds\n")


##################################################################
# Allen MTG CARM1P1 GAM-REML-PRSS Analysis - Part 2
# Plots and console summary statistics
##################################################################
library(mgcv)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(viridis)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data    <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results <- rds_data$all_results
donor_data  <- rds_data$donor_data
donors      <- rds_data$donors
sheet_loo   <- rds_data$sheet_loo

valid_results <- Filter(Negate(is.null), all_results)
cat("Valid fits loaded:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# SECTION 2: Shared theme
# ------------------------------------------------------------------
theme_gam <- function(base_size = 9) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1),
            plot.subtitle    = element_text(size = base_size - 1),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            legend.title     = element_text(size = base_size - 1),
            legend.text      = element_text(size = base_size - 2),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.35),
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
}

# ------------------------------------------------------------------
# SECTION 3: PRSS + REML plots - best model only
# ------------------------------------------------------------------
cat("========== Building PRSS + REML plots ==========\n")

PRSS_PLOT_MAX_ITER <- 28L   # show only first 28 PRSS iterations

# ---- 3a. PRSS plot builder ----
make_prss_plot <- function(result, compact = FALSE) {
    bs <- if (compact) 7 else 9
    
    # Step 1: hard-cap at PRSS_PLOT_MAX_ITER (absolute ceiling = 28)
    n_cap  <- min(PRSS_PLOT_MAX_ITER, result$actual_iters)
    iters  <- seq_len(n_cap)
    prss_v <- result$prss_values[iters]
    
    df <- data.frame(
        Iteration = iters,
        PRSS      = prss_v,
        Is_Best   = iters == result$best_iter
    )
    df <- df[!is.na(df$PRSS), ]
    
    # Step 2: find the first iteration where PRSS drops to near-zero and
    prss_max   <- max(df$PRSS, na.rm = TRUE)
    zero_thresh <- max(prss_max * 0.01, 1e-6)
    first_zero  <- which(df$PRSS < zero_thresh)[1]   # first near-zero index
    
    if (!is.na(first_zero) && first_zero > 1) {
        # Keep rows up to - but not including - the first near-zero iteration
        df <- df[seq_len(first_zero - 1), ]
    }
    # (If first_zero == 1 all values are near-zero; keep df as-is so plot is not blank)
    
    best_row  <- df[df$Is_Best, ]
    beyond_28 <- result$best_iter > PRSS_PLOT_MAX_ITER
    
    # If best_iter was trimmed away (rare), mark the minimum of what remains
    if (nrow(best_row) == 0 && nrow(df) > 0) {
        best_row  <- df[which.min(df$PRSS), ]
        beyond_28 <- TRUE   # treat as "beyond shown range" → grey marker
    }
    
    x_max <- max(df$Iteration, na.rm = TRUE)
    
    ggplot(df, aes(x = Iteration, y = PRSS)) +
        geom_line(linewidth = 0.6, color = "gray30") +
        geom_point(aes(color = Is_Best, size = Is_Best)) +
        { if (nrow(best_row) > 0)
            geom_vline(xintercept = best_row$Iteration[1],
                       linetype  = "dashed",
                       color     = if (beyond_28) "gray50" else "firebrick",
                       linewidth = 0.5) } +
        { if (nrow(best_row) > 0)
            geom_point(data = best_row,
                       aes(x = Iteration, y = PRSS),
                       color = if (beyond_28) "gray40" else "firebrick",
                       size  = 3.5, shape = 18,
                       inherit.aes = FALSE) } +
        scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "firebrick"),
                           guide = "none") +
        scale_size_manual(values  = c("FALSE" = 1.2,     "TRUE" = 3.5),
                          guide = "none") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
        coord_cartesian(xlim = c(1, x_max)) +
        labs(title    = paste0(result$donor, " | ", result$gene_set, " - PRSS"),
             subtitle = paste0("Best iter: ", result$best_iter,
                               "  k = ", result$best_k,
                               "  PRSS = ", round(result$best_prss, 3),
                               if (beyond_28 && result$best_iter <= PRSS_PLOT_MAX_ITER) "" else
                                   if (result$best_iter > PRSS_PLOT_MAX_ITER)
                                       "  [best beyond iter 28]" else ""),
             x = "PRSS Iteration",
             y = "PRSS") +
        theme_gam(base_size = bs)
}

# ---- 3b. REML plot builder (best PRSS iteration only) ----
make_reml_plot <- function(result, compact = FALSE) {
    bs     <- if (compact) 7 else 9
    all_ri <- result$all_reml_iters
    
    if (is.null(all_ri) || nrow(all_ri) == 0)
        return(ggplot() +
                   annotate("text", x = 0.5, y = 0.5, size = 3,
                            label = "REML data unavailable") +
                   labs(title = paste0(result$donor, " | ", result$gene_set,
                                       " - REML")) +
                   theme_void())
    
    # Filter to rows belonging to best PRSS iteration; fallback progressively
    best_ri <- all_ri[all_ri$is_best_model == TRUE, ]
    if (nrow(best_ri) == 0)
        best_ri <- all_ri[all_ri$prss_iteration == result$best_iter, ]
    if (nrow(best_ri) == 0) best_ri <- all_ri
    
    best_ri <- best_ri[order(best_ri$iteration), ]
    
    final_row <- best_ri[nrow(best_ri), ]
    
    ggplot(best_ri, aes(x = iteration, y = score)) +
        geom_line(linewidth = 0.6, color = "gray30") +
        geom_point(size = 2, color = "gray20") +
        geom_point(data = final_row,
                   aes(x = iteration, y = score),
                   color = "firebrick", size = 3.5, shape = 18,
                   inherit.aes = FALSE) +
        labs(title    = paste0(result$donor, " | ", result$gene_set, " - REML"),
             subtitle = paste0("PRSS iter ", result$best_iter,
                               "  |  Final REML = ",
                               round(final_row$score, 3)),
             x = "REML Iteration", y = "REML Score") +
        theme_gam(base_size = bs)
}

# ---- 3c. Build all panels then display as multi-page grid ----
all_panel_pairs <- list()

for (key in names(valid_results)) {
    res <- valid_results[[key]]
    cat(" Building panels for:", key, "\n")
    p_prss <- tryCatch(make_prss_plot(res, compact = TRUE),
                       error = function(e) {
                           cat("  PRSS error:", conditionMessage(e), "\n"); NULL })
    p_reml <- tryCatch(make_reml_plot(res, compact = TRUE),
                       error = function(e) {
                           cat("  REML error:", conditionMessage(e), "\n"); NULL })
    all_panel_pairs[[key]] <- list(prss = p_prss, reml = p_reml, key = key)
}

n_fits  <- length(all_panel_pairs)
fits_pp <- 4L
n_pages <- ceiling(n_fits / fits_pp)
keys_vec <- names(all_panel_pairs)

cat("\nDisplaying PRSS + REML grid (", n_fits, "fits,",
    n_pages, "page(s),", fits_pp, "fits/page)...\n")

for (pg in seq_len(n_pages)) {
    idx_start <- (pg - 1L) * fits_pp + 1L
    idx_end   <- min(pg * fits_pp, n_fits)
    page_keys <- keys_vec[idx_start:idx_end]
    
    # Flatten: col1_prss, col1_reml, col2_prss, col2_reml, …
    panel_list <- list()
    for (k in page_keys) {
        pp <- all_panel_pairs[[k]]
        panel_list <- c(
            panel_list,
            list(if (!is.null(pp$prss)) pp$prss else
                ggplot() + theme_void() +
                    annotate("text", 0.5, 0.5, label = paste(k, "PRSS N/A"))),
            list(if (!is.null(pp$reml)) pp$reml else
                ggplot() + theme_void() +
                    annotate("text", 0.5, 0.5, label = paste(k, "REML N/A"))))
    }
    
    n_cols <- length(page_keys)
    page_title <- paste0("PRSS (first 28 iters) & REML - Best Model per Fit  (page ",
                         pg, " of ", n_pages, ")")
    grid_obj <- tryCatch(
        gridExtra::grid.arrange(
            grobs = panel_list,
            ncol  = n_cols,
            nrow  = 2L,
            top   = grid::textGrob(
                page_title,
                gp = grid::gpar(fontface = "bold", fontsize = 11))),
        error = function(e) {
            cat("  grid.arrange error (page", pg, "):", conditionMessage(e), "\n")
            NULL
        }
    )
    if (!is.null(grid_obj)) print(grid_obj)
}

# ------------------------------------------------------------------
# SECTION 4: GAM curve plots - one per donor, DFG + HKG overlaid
# ------------------------------------------------------------------
cat("\n========== GAM curve plots ==========\n")

gene_set_colors <- c("DFG" = "#4B0082", "HKG" = "#4F75DE")
gene_set_fills  <- c("DFG" = "#4B0082", "HKG" = "#4F75DE")

make_gam_curve_plot <- function(donor, vres) {
    plot_data <- data.frame()
    for (gs in c("DFG", "HKG")) {
        res <- vres[[paste0(donor, "|", gs)]]
        if (is.null(res) || is.null(res$best_model)) next
        m    <- res$best_model
        data <- res$data
        pts  <- data.frame(
            CARM1P1    = data$CARM1P1,
            Expression = data$Expression,
            Gene_Set   = gs,
            se.fit     = NA_real_,
            type       = "data",
            stringsAsFactors = FALSE)
        set.seed(123)
        x_seq <- seq(min(data$CARM1P1), max(data$CARM1P1), length.out = 1000)
        pred  <- predict(m, newdata = data.frame(CARM1P1 = x_seq), se.fit = TRUE)
        fit_df <- data.frame(
            CARM1P1    = x_seq,
            Expression = pred$fit,
            Gene_Set   = gs,
            se.fit     = pred$se.fit,
            type       = "fit",
            stringsAsFactors = FALSE)
        plot_data <- rbind(plot_data, pts, fit_df)
    }
    if (nrow(plot_data) == 0)
        return(ggplot() +
                   annotate("text", x = 0.5, y = 0.5,
                            label = paste("No data for donor:", donor)) +
                   theme_void())
    
    de_labels <- sapply(c("DFG", "HKG"), function(gs) {
        res <- vres[[paste0(donor, "|", gs)]]
        if (is.null(res) || is.null(res$best_model)) return(NA_character_)
        sprintf("%s DE=%.1f%%", gs, summary(res$best_model)$dev.expl * 100)
    })
    de_labels <- de_labels[!is.na(de_labels)]
    
    set.seed(123)
    ggplot(plot_data, aes(x = CARM1P1, y = Expression, color = Gene_Set)) +
        geom_point(data = subset(plot_data, type == "data"),
                   position = position_jitter(width = 0.01, height = 0, seed = 123),
                   alpha = 0.25, size = 0.8) +
        geom_line(data  = subset(plot_data, type == "fit"), linewidth = 1.2) +
        geom_ribbon(data = subset(plot_data, type == "fit"),
                    aes(ymin = Expression - 1.96 * se.fit,
                        ymax = Expression + 1.96 * se.fit,
                        fill = Gene_Set),
                    alpha = 0.20, color = NA) +
        scale_color_manual(values = gene_set_colors, name = "Gene Set") +
        scale_fill_manual( values = gene_set_fills,  name = "Gene Set") +
        labs(title    = paste("GAM: Gene Set Expression vs CARM1P1 -", donor),
             subtitle = paste(de_labels, collapse = "  |  "),
             x = "CARM1P1 Expression (log2)",
             y = "Gene Set Mean Expression (log2)") +
        scale_x_continuous(labels = number_format(accuracy = 0.1)) +
        scale_y_continuous(labels = number_format(accuracy = 0.1)) +
        theme_gam() +
        theme(legend.position = "right")
}

for (donor in donors) {
    has_result <- any(sapply(c("DFG", "HKG"), function(gs)
        !is.null(valid_results[[paste0(donor, "|", gs)]])))
    if (!has_result) { cat("  Skipping", donor, "(no valid fits)\n"); next }
    cat(" ", donor, "\n")
    p <- tryCatch(make_gam_curve_plot(donor, valid_results),
                  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(p)) print(p)
}

# ------------------------------------------------------------------
# SECTION 5: LOO gene contribution plots
# ------------------------------------------------------------------
cat("\n========== LOO gene contribution plots ==========\n")

if (!is.null(sheet_loo) && nrow(sheet_loo) > 0) {
    
    # ---- 5a. Heatmap: all donors × genes, faceted by gene set ----
    loo_clean <- sheet_loo[!is.na(sheet_loo$Delta_DE), ]
    if (nrow(loo_clean) > 0) {
        loo_clean$Donor_short <- sub("^H", "", loo_clean$Donor)
        p_heat <- ggplot(loo_clean,
                         aes(x = Removed_Gene, y = Donor_short,
                             fill = Delta_DE * 100)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = sprintf("%.1f", Delta_DE * 100)),
                      size = 2.6, color = "black") +
            facet_wrap(~ Gene_Set, scales = "free_x") +
            scale_fill_distiller(palette = "Spectral", name = expression(Delta*"DE (%)"),
                     labels = number_format(accuracy = 0.1)) +
            labs(title    = "Leave-One-Out Gene Contribution to Deviance Explained",
                 subtitle = paste0("\u0394DE = DE(full 5-gene) \u2212 DE(4-gene LOO)  |  ",
                                   "red = gene helps fit,  blue = gene suppresses fit"),
                 x = "Removed Gene", y = "Donor") +
            theme_gam() +
            theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                  strip.text  = element_text(face = "bold", size = 10))
        print(p_heat)
    }
    
    # ---- 5b. Per-donor bar plots combined into one grid (2 columns) ----
    bar_plots <- list()
    for (donor in unique(sheet_loo$Donor)) {
        dd_loo <- sheet_loo[sheet_loo$Donor == donor & !is.na(sheet_loo$Delta_DE), ]
        if (nrow(dd_loo) == 0) next
        gene_order <- dd_loo %>%
            group_by(Removed_Gene) %>%
            summarise(mean_delta = mean(Delta_DE, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(mean_delta)) %>%
            pull(Removed_Gene)
        dd_loo$Removed_Gene <- factor(dd_loo$Removed_Gene, levels = gene_order)
        p_bar <- ggplot(dd_loo,
                        aes(x = Removed_Gene, y = Delta_DE * 100, fill = Gene_Set)) +
            geom_col(position = position_dodge(width = 0.7), width = 0.6) +
            geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
            geom_text(aes(label = sprintf("%.1f%%", Delta_DE * 100),
                          vjust = ifelse(Delta_DE >= 0, -0.4, 1.2)),
                      position = position_dodge(width = 0.7),
                      size = 2.2, color = "gray20") +
            scale_fill_manual(values = c("DFG" = "#4B0082", "HKG" = "#4F75DE"),
                              name = "Gene Set") +
            labs(title = sub("^H", "", donor), x = NULL,
                 y = expression(Delta*"DE (%)")) +
            theme_gam(base_size = 7) +
            theme(legend.position = "none",
                  axis.text.x     = element_text(angle = 30, hjust = 1))
        bar_plots[[donor]] <- p_bar
    }
    
    if (length(bar_plots) > 0) {
        cat("  Displaying LOO bar plots grid (", length(bar_plots), "donors)...\n")
        grid_bars <- tryCatch(
            gridExtra::grid.arrange(
                grobs = bar_plots,
                ncol  = 2,
                top   = grid::textGrob(
                    paste0("LOO Gene Contribution - \u0394DE per Gene per Donor\n",
                           "\u0394DE = DE(full) \u2212 DE(LOO);  positive = gene contributes"),
                    gp = grid::gpar(fontface = "bold", fontsize = 10))),
            error = function(e) {
                cat("  grid.arrange error:", conditionMessage(e), "\n"); NULL })
        if (!is.null(grid_bars)) print(grid_bars)
    }
    
} else {
    cat("  No LOO data found - skipping\n")
}

# ------------------------------------------------------------------
# SECTION 6: Console summary statistics
# ------------------------------------------------------------------
cat("\n========== Summary Statistics ==========\n")

for (gs in c("DFG", "HKG")) {
    gs_res <- Filter(function(r) r$gene_set == gs, valid_results)
    if (length(gs_res) == 0) next
    
    k_vals    <- sapply(gs_res, function(r) r$best_k)
    prss_vals <- sapply(gs_res, function(r) r$best_prss)
    de_vals   <- sapply(gs_res, function(r)
        if (!is.null(r$best_model)) summary(r$best_model)$dev.expl else NA_real_)
    edf_vals  <- sapply(gs_res, function(r)
        if (!is.null(r$best_model)) {
            sm <- summary(r$best_model)
            if (nrow(sm$s.table) > 0) sm$s.table[1, "edf"] else NA_real_
        } else NA_real_)
    iter_vals <- sapply(gs_res, function(r) r$best_iter)
    
    cat("\nGene Set:", gs, "(", length(gs_res), "donors)\n")
    cat("  Mean best k:       ", round(mean(k_vals,    na.rm=TRUE), 2), "\n")
    cat("  Mean best PRSS:    ", round(mean(prss_vals, na.rm=TRUE), 4), "\n")
    cat("  Mean DE:           ", round(mean(de_vals,   na.rm=TRUE)*100, 2), "%\n")
    cat("  Mean smooth EDF:   ", round(mean(edf_vals,  na.rm=TRUE), 4), "\n")
    cat("  Mean best iter:    ", round(mean(iter_vals, na.rm=TRUE), 1), "\n")
}

cat("\n========== Part 2 complete ==========\n")
cat("All plots displayed in RStudio Plots pane (no files saved).\n")
cat("All tabular results are in: AllenMTG_CARM1P1_GAM_Results_1.xlsx\n")

# Heatmap of CARM1P1 vs Axon Guidance and Housekeeping Genes
library(ComplexHeatmap)

# Subset data (specific clusters for heatmap)
# MODIFIED: single object, 8 MTG clusters; no merge needed
heatmap_obj <- subset(mtg_results$seurat_obj, idents = c(0, 3, 9, 15))

# Define gene sets
# MODIFIED: CARM1P1 correlated axon guidance genes + housekeeping genes
genes_of_interest <- c("CARM1P1", "CLMN", "EPHA3", "EPHA6", "LOC101928964", "ROBO2",
                       "ACTB", "GAPDH", "PPIA", "RPL13A", "TBP")
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(heatmap_obj)]

# MODIFIED: extract CARM1P1 expression for sorting (replaces trpm4_expression)
carm1p1_expression <- GetAssayData(seurat_obj, slot = "data")["CARM1P1",
                          Cells(heatmap_obj)]

# Prepare expression matrix (log2-transformed)
# MODIFIED: no averaged composite rows (Ribo/AR) - all genes shown individually
expression_matrix     <- GetAssayData(seurat_obj, slot = "data")[genes_of_interest,
                             Cells(heatmap_obj)]
expression_matrix_log <- log2(expression_matrix + 1)

# Scale to z-scores (cap at -1 to 2.3)
scale_to_zscore <- function(x) {
    z <- scale(x)
    return(pmin(pmax(z, -1), 2.3))
}
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Order cells by cluster and CARM1P1 expression (ascending)
# MODIFIED: cluster levels = 8 MTG clusters; sort key = carm1p1_expression
heatmap_grouping <- Idents(heatmap_obj)
heatmap_grouping <- factor(heatmap_grouping,
                            levels = c("0", "3", "9", "15"))
cell_order <- order(heatmap_grouping, carm1p1_expression)
expression_matrix_log  <- expression_matrix_log[, cell_order]
heatmap_grouping       <- heatmap_grouping[cell_order]
carm1p1_expression     <- carm1p1_expression[cell_order]
expression_matrix_scaled <- t(apply(expression_matrix_log, 1, scale_to_zscore))

# Prepare heatmap annotations (cluster labels with counts)
n_groups         <- length(levels(heatmap_grouping))
group_colors     <- viridis(n_groups, option = "mako"); names(group_colors) <- levels(heatmap_grouping)
cell_counts      <- table(heatmap_grouping); cell_percentages <- round(prop.table(cell_counts) * 100, 1)
group_labels     <- paste0(names(cell_counts), "\n", cell_counts, "\n(", cell_percentages, "%)")
column_annotation <- HeatmapAnnotation(
    Cluster = anno_block(gp = gpar(fill = group_colors), labels = group_labels,
                         labels_gp = gpar(col = "white", fontsize = 8)),
    show_legend = FALSE, show_annotation_name = FALSE)

# Create heatmap
# MODIFIED: output filename; width adjusted for 8 clusters
heatmap <- Heatmap(expression_matrix_scaled, name = "zscore", column_title = NULL, row_title = "Genes",
                   show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                   top_annotation = column_annotation, col = viridis(100, option = "mako"), row_names_gp = gpar(fontsize = 8),
                   column_split = heatmap_grouping, row_gap = unit(1, "mm"), column_gap = unit(0.5, "mm"), border = TRUE,
                   use_raster = TRUE, raster_quality = 2)

# Display heatmap
draw(heatmap)


##############################################################################
# Allen MTG CARM1P1 - Part 3: GAM Deviance Components & Explanatory Power
##############################################################################
library(mgcv)
library(writexl)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data    <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results <- rds_data$all_results
donors      <- rds_data$donors

valid_results <- Filter(Negate(is.null), all_results)
gene_sets     <- c("DFG", "HKG")

cat("Valid fits loaded:", length(valid_results), "\n")
cat("Donors:", paste(donors, collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# SECTION 2: Deviance component functions 
# ------------------------------------------------------------------

# ---- 2a. Extract all deviance components from one model ----
extract_deviance_components <- function(model, data) {
    set.seed(123)
    n             <- nrow(data)
    response_var  <- "Expression"           # fixed for CARM1P1 analysis
    y             <- data[[response_var]]

    fitted_values      <- fitted(model)
    residuals_response <- residuals(model, type = "response")

    edf_total   <- sum(model$edf)
    df_residual <- model$df.residual
    df_model    <- n - df_residual

    scale_parameter <- model$scale

    # Null model: intercept only
    null_model <- gam(Expression ~ 1, data = data, family = gaussian())
    y_null     <- fitted(null_model)
    y_mean     <- mean(y)

    cat("  n =", n, "| total EDF =", round(edf_total, 3),
        "| residual DF =", round(df_residual, 3), "\n")
    cat("  Null fitted =", round(unique(y_null)[1], 6),
        "| Mean(y) =", round(y_mean, 6), "\n")

    # Sum of squares
    ss_total    <- sum((y - y_mean)^2)
    ss_residual <- sum((y - fitted_values)^2)
    ss_model    <- ss_total - ss_residual

    # Null deviance
    null_deviance_raw    <- sum((y - y_mean)^2)
    null_deviance_scaled <- null_deviance_raw / scale_parameter
    null_deviance_mgcv   <- if (!is.null(null_model$null.deviance))
                                null_model$null.deviance else deviance(null_model)

    # Model deviance
    model_deviance_raw    <- sum((y - fitted_values)^2)
    model_deviance_scaled <- model_deviance_raw / scale_parameter
    model_deviance_alt    <- df_residual * scale_parameter
    model_deviance_mgcv   <- deviance(model)

    # Deviance explained - four methods
    dev_explained_raw     <- 1 - (model_deviance_raw    / null_deviance_raw)
    dev_explained_scaled  <- 1 - (model_deviance_scaled / null_deviance_scaled)
    dev_explained_mgcv    <- 1 - (model_deviance_mgcv   / null_deviance_mgcv)
    dev_explained_summary <- summary(model)$dev.expl

    list(
        n_observations        = n,
        response_variable     = response_var,
        edf_total             = edf_total,
        df_residual           = df_residual,
        df_model              = df_model,
        scale_parameter       = scale_parameter,
        y_mean                = y_mean,
        y_range               = range(y),
        fitted_range          = range(fitted_values),
        ss_total              = ss_total,
        ss_residual           = ss_residual,
        ss_model              = ss_model,
        null_deviance_raw     = null_deviance_raw,
        null_deviance_scaled  = null_deviance_scaled,
        null_deviance_mgcv    = null_deviance_mgcv,
        model_deviance_raw    = model_deviance_raw,
        model_deviance_scaled = model_deviance_scaled,
        model_deviance_alt    = model_deviance_alt,
        model_deviance_mgcv   = model_deviance_mgcv,
        dev_explained_raw     = dev_explained_raw,
        dev_explained_scaled  = dev_explained_scaled,
        dev_explained_mgcv    = dev_explained_mgcv,
        dev_explained_summary = dev_explained_summary,
        scale_check           = all.equal(model_deviance_raw / df_residual,
                                          scale_parameter),
        deviance_alt_check    = all.equal(model_deviance_alt, model_deviance_mgcv)
    )
}

# ---- 2b. Run extract_deviance_components across all 16 fits ----
analyze_all_deviances <- function(valid_res) {
    deviance_data <- list()
    for (key in names(valid_res)) {
        res <- valid_res[[key]]
        parts          <- strsplit(key, "\\|")[[1]]
        donor_id       <- parts[1]
        gene_set_label <- parts[2]
        cat("Processing:", key, "\n")
        tryCatch({
            dc <- extract_deviance_components(res$best_model, res$data)
            dc$donor      <- donor_id
            dc$gene_set   <- gene_set_label
            dc$group      <- "Allen_MTG"
            deviance_data[[key]] <- dc
            cat("  Dev explained (summary):",
                round(dc$dev_explained_summary, 4), "\n")
        }, error = function(e) {
            cat("  Error:", conditionMessage(e), "\n")
        })
    }
    deviance_data
}

# ---- 2c. Flatten deviance_data into a wide summary data frame ----
create_deviance_summary <- function(deviance_data) {
    if (length(deviance_data) == 0) return(data.frame())
    rows <- lapply(deviance_data, function(dc) {
        data.frame(
            Donor                 = dc$donor,
            Group                 = dc$group,
            Gene_Set              = dc$gene_set,
            N_Observations        = dc$n_observations,
            Response_Variable     = dc$response_variable,
            EDF_Total             = dc$edf_total,
            DF_Residual           = dc$df_residual,
            DF_Model              = dc$df_model,
            Scale_Parameter       = dc$scale_parameter,
            Y_Mean                = dc$y_mean,
            Y_Min                 = dc$y_range[1],
            Y_Max                 = dc$y_range[2],
            SS_Total              = dc$ss_total,
            SS_Residual           = dc$ss_residual,
            SS_Model              = dc$ss_model,
            Null_Deviance_Raw     = dc$null_deviance_raw,
            Null_Deviance_Scaled  = dc$null_deviance_scaled,
            Null_Deviance_mgcv    = dc$null_deviance_mgcv,
            Model_Deviance_Raw    = dc$model_deviance_raw,
            Model_Deviance_Scaled = dc$model_deviance_scaled,
            Model_Deviance_Alt    = dc$model_deviance_alt,
            Model_Deviance_mgcv   = dc$model_deviance_mgcv,
            Dev_Explained_Raw     = dc$dev_explained_raw,
            Dev_Explained_Scaled  = dc$dev_explained_scaled,
            Dev_Explained_mgcv    = dc$dev_explained_mgcv,
            Dev_Explained_Summary = dc$dev_explained_summary,
            Scale_Check           = as.character(dc$scale_check),
            Deviance_Alt_Check    = as.character(dc$deviance_alt_check),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

# ---- 2d. Formula reference sheet ----
create_deviance_explanations <- function() {
    data.frame(
        Component = c(
            "EDF_Total", "DF_Residual", "DF_Model", "Scale_Parameter",
            "Y_Mean", "Y_Min", "Y_Max",
            "SS_Total", "SS_Residual", "SS_Model",
            "Null_Deviance_Raw", "Null_Deviance_Scaled", "Null_Deviance_mgcv",
            "Model_Deviance_Raw", "Model_Deviance_Scaled",
            "Model_Deviance_Alt", "Model_Deviance_mgcv",
            "Dev_Explained_Raw", "Dev_Explained_Scaled",
            "Dev_Explained_mgcv", "Dev_Explained_Summary",
            "Scale_Check", "Deviance_Alt_Check"
        ),
        Formula = c(
            "Sum of effective degrees of freedom",
            "n - df_model",
            "n - df_residual",
            "phi = sum((yi - f(xi))^2) / (n - df)",
            "Mean of response variable",
            "Minimum value of response variable",
            "Maximum value of response variable",
            "sum((yi - y_mean)^2)",
            "sum((yi - f(xi))^2)",
            "SS_Total - SS_Residual",
            "sum((yi - y_mean)^2)",
            "sum((yi - y_mean)^2) / phi",
            "mgcv null deviance calculation",
            "sum((yi - f(xi))^2)",
            "sum((yi - f(xi))^2) / phi",
            "(n - df) x phi",
            "mgcv model deviance calculation",
            "1 - (Model_Deviance_Raw / Null_Deviance_Raw)",
            "1 - (Model_Deviance_Scaled / Null_Deviance_Scaled)",
            "1 - (Model_Deviance_mgcv / Null_Deviance_mgcv)",
            "mgcv summary deviance explained",
            "Verification: phi = Model_Deviance_Raw / DF_Residual",
            "Verification: Model_Deviance_Alt = Model_Deviance_mgcv"
        ),
        Description = c(
            "Total effective degrees of freedom used by smooth terms",
            "Residual degrees of freedom for error estimation",
            "Model degrees of freedom",
            "Scale parameter estimate (residual variance)",
            "Mean of the response variable values",
            "Minimum value in the response variable",
            "Maximum value in the response variable",
            "Total sum of squares (same as raw null deviance)",
            "Residual sum of squares (same as raw model deviance)",
            "Explained sum of squares by the model",
            "Raw null deviance: sum of squared deviations from mean",
            "Scaled null deviance: raw null deviance / scale parameter",
            "Null deviance as calculated by mgcv package",
            "Raw model deviance: sum of squared residuals from fitted model",
            "Scaled model deviance: raw model deviance / scale parameter",
            "Alternative model deviance: (n-df) times scale parameter",
            "Model deviance as calculated by mgcv package",
            "Deviance explained using raw (unscaled) deviances",
            "Deviance explained using scaled deviances (standard Gaussian)",
            "Deviance explained using mgcv deviance calculations",
            "Deviance explained from mgcv model summary",
            "Logical check: scale parameter calculation correct",
            "Logical check: alternative deviance matches mgcv"
        ),
        stringsAsFactors = FALSE
    )
}

# ---- 2e. Long-format detailed breakdown (one row per metric per fit) ----
create_deviance_detailed <- function(deviance_data) {
    if (length(deviance_data) == 0) return(data.frame())
    rows <- lapply(names(deviance_data), function(key) {
        dc <- deviance_data[[key]]
        # Scalar fields only (exclude list/vector fields like y_range, fitted_range)
        scalar_fields <- c(
            "n_observations", "response_variable",
            "edf_total", "df_residual", "df_model", "scale_parameter",
            "y_mean", "ss_total", "ss_residual", "ss_model",
            "null_deviance_raw", "null_deviance_scaled", "null_deviance_mgcv",
            "model_deviance_raw", "model_deviance_scaled",
            "model_deviance_alt", "model_deviance_mgcv",
            "dev_explained_raw", "dev_explained_scaled",
            "dev_explained_mgcv", "dev_explained_summary",
            "scale_check", "deviance_alt_check"
        )
        # Add range fields as separate named entries
        extra <- list(
            y_min         = dc$y_range[1],
            y_max         = dc$y_range[2],
            fitted_min    = dc$fitted_range[1],
            fitted_max    = dc$fitted_range[2]
        )
        vals <- c(dc[scalar_fields], extra)
        data.frame(
            Donor    = dc$donor,
            Group    = dc$group,
            Gene_Set = dc$gene_set,
            Metric   = names(vals),
            Value    = sapply(vals, function(x) as.character(x[1])),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

# ------------------------------------------------------------------
# SECTION 3: Explanatory power functions 
# ------------------------------------------------------------------

# ---- 3a. Cell-level EP extraction for one fit ----
extract_ep_for_fit <- function(result, donor_id, gene_set_label) {
    cat("  Extracting EP for:", donor_id, "|", gene_set_label, "\n")

    model      <- result$best_model
    data       <- result$data
    n_cells    <- nrow(data)

    # Overall model deviance explained
    model_dev_explained <- summary(model)$dev.expl

    # Null model fitted values (= mean of y)
    null_model    <- gam(Expression ~ 1, data = data)
    null_fitted   <- fitted(null_model)
    model_fitted  <- fitted(model)

    null_residuals  <- data$Expression - null_fitted
    model_residuals <- data$Expression - model_fitted

    null_sq_diff  <- null_residuals^2
    model_sq_diff <- model_residuals^2

    # Individual cell explanatory power
    ep_value <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )

    # Sort cells by EP (highest first)
    sorted_idx       <- order(ep_value, decreasing = TRUE)
    sorted_cell_names <- rownames(data)[sorted_idx]

    # Deviance cells: top (model_dev_explained × n) cells by EP
    target_cells      <- round(n_cells * model_dev_explained)
    target_cells      <- max(1L, min(target_cells, n_cells))  # guard bounds
    deviance_cells    <- sorted_cell_names[seq_len(target_cells)]
    non_deviance_cells <- sorted_cell_names[(target_cells + 1):n_cells]

    output_df <- data.frame(
        Donor                = donor_id,
        Group                = "Allen_MTG",
        Gene_Set             = gene_set_label,
        Cell_Index           = rownames(data),
        EP_Value             = ep_value,
        Sorted_Index         = match(rownames(data), sorted_cell_names),
        Deviance_Cell        = ifelse(rownames(data) %in% deviance_cells,
                                      "YES", "NO"),
        Non_Deviance_Cell    = ifelse(rownames(data) %in% non_deviance_cells,
                                      "YES", "NO"),
        CARM1P1_Expression   = data$CARM1P1,
        Expression           = data$Expression,
        Model_Dev_Explained  = model_dev_explained,
        stringsAsFactors = FALSE
    )

    # Sort output by EP (highest first)
    output_df <- output_df[order(output_df$EP_Value, decreasing = TRUE), ]
    rownames(output_df) <- NULL

    cat("    Cells:", n_cells,
        "| Dev explained:", round(model_dev_explained, 4),
        "| Target deviance cells:", target_cells, "\n")
    output_df
}

# ------------------------------------------------------------------
# SECTION 4: Run all analyses
# ------------------------------------------------------------------
cat("\n========== SECTION 4: Running deviance component analysis ==========\n")

deviance_data    <- analyze_all_deviances(valid_results)
deviance_summary <- create_deviance_summary(deviance_data)
deviance_expl    <- create_deviance_explanations()
deviance_detail  <- create_deviance_detailed(deviance_data)

cat("\nDeviance summary rows:", nrow(deviance_summary), "\n")

cat("\n========== SECTION 5: Running cell-level EP analysis ==========\n")

all_ep_list     <- list()   # all fits combined
per_fit_ep      <- list()   # named list for individual sheets

for (key in names(valid_results)) {
    parts          <- strsplit(key, "\\|")[[1]]
    donor_id       <- parts[1]
    gene_set_label <- parts[2]
    cat("Processing EP:", key, "\n")
    tryCatch({
        ep_df <- extract_ep_for_fit(valid_results[[key]], donor_id, gene_set_label)
        all_ep_list[[key]]  <- ep_df
        # Sheet name: "EP_<donor>_<gene_set>"  (max 31 chars for Excel)
        sheet_nm <- paste0("EP_",
                           gsub("\\.", "", donor_id), "_", gene_set_label)
        sheet_nm <- substr(sheet_nm, 1, 31)
        per_fit_ep[[sheet_nm]] <- ep_df
    }, error = function(e) {
        cat("  Error:", conditionMessage(e), "\n")
    })
}

ep_combined <- if (length(all_ep_list) > 0)
                   do.call(rbind, all_ep_list) else data.frame()
rownames(ep_combined) <- NULL

# EP summary: one row per fit
ep_summary <- if (length(all_ep_list) > 0) {
    do.call(rbind, lapply(names(all_ep_list), function(key) {
        df <- all_ep_list[[key]]
        data.frame(
            Donor               = df$Donor[1],
            Group               = df$Group[1],
            Gene_Set            = df$Gene_Set[1],
            Total_Cells         = nrow(df),
            Mean_EP_Value       = round(mean(df$EP_Value,  na.rm = TRUE), 4),
            Max_EP_Value        = round(max(df$EP_Value,   na.rm = TRUE), 4),
            Min_EP_Value        = round(min(df$EP_Value,   na.rm = TRUE), 4),
            Deviance_Cells_Count = sum(df$Deviance_Cell == "YES"),
            Model_Dev_Explained = round(unique(df$Model_Dev_Explained), 4),
            stringsAsFactors = FALSE
        )
    }))
} else data.frame()

cat("\nEP combined rows:", nrow(ep_combined), "\n")
cat("EP per-fit sheets:", length(per_fit_ep), "\n")

# ------------------------------------------------------------------
# SECTION 6: Assemble and export single Excel file
# ------------------------------------------------------------------
cat("\n========== SECTION 6: Exporting to Excel ==========\n")

export_list <- c(
    list(
        Deviance_Summary      = deviance_summary,
        Deviance_Explanations = deviance_expl,
        Deviance_Detailed     = deviance_detail,
        EP_All_Combined       = ep_combined,
        EP_Summary            = ep_summary
    ),
    per_fit_ep   # individual EP sheets: EP_<donor>_<gene_set>
)

out_file <- "AllenMTG_CARM1P1_GAM_Results_3.xlsx"

tryCatch({
    write_xlsx(export_list, path = out_file)
    cat("\nSaved:", out_file, "\n")
    cat("Sheets:\n")
    cat("  Deviance_Summary      :", nrow(deviance_summary),
        "rows (one per donor x gene set)\n")
    cat("  Deviance_Explanations :", nrow(deviance_expl),
        "rows (formula reference)\n")
    cat("  Deviance_Detailed     :", nrow(deviance_detail),
        "rows (long-format per-metric breakdown)\n")
    cat("  EP_All_Combined       :", nrow(ep_combined),
        "rows (cell-level EP, all fits)\n")
    cat("  EP_Summary            :", nrow(ep_summary),
        "rows (one per donor x gene set)\n")
    for (nm in names(per_fit_ep))
        cat(" ", nm, ":", nrow(per_fit_ep[[nm]]), "rows\n")
}, error = function(e) {
    cat("Error writing Excel:", conditionMessage(e), "\n")
})

# ------------------------------------------------------------------
# SECTION 7: Console summary statistics  (mirrors Part 3.09 output)
# ------------------------------------------------------------------
cat("\n========== Summary Statistics ==========\n")

if (nrow(deviance_summary) > 0) {
    cat("\nScale parameter (phi):\n")
    cat("  Min:", round(min(deviance_summary$Scale_Parameter,  na.rm=TRUE), 6), "\n")
    cat("  Max:", round(max(deviance_summary$Scale_Parameter,  na.rm=TRUE), 6), "\n")
    cat("  Mean:", round(mean(deviance_summary$Scale_Parameter, na.rm=TRUE), 6), "\n")

    cat("\nDeviance explained (Dev_Explained_Summary):\n")
    cat("  Min:", round(min(deviance_summary$Dev_Explained_Summary,  na.rm=TRUE), 4), "\n")
    cat("  Max:", round(max(deviance_summary$Dev_Explained_Summary,  na.rm=TRUE), 4), "\n")
    cat("  Mean:", round(mean(deviance_summary$Dev_Explained_Summary, na.rm=TRUE), 4), "\n")

    cat("\nTotal EDF:\n")
    cat("  Min:", round(min(deviance_summary$EDF_Total,  na.rm=TRUE), 3), "\n")
    cat("  Max:", round(max(deviance_summary$EDF_Total,  na.rm=TRUE), 3), "\n")
    cat("  Mean:", round(mean(deviance_summary$EDF_Total, na.rm=TRUE), 3), "\n")

    for (gs in gene_sets) {
        sub_df <- deviance_summary[deviance_summary$Gene_Set == gs, ]
        cat("\nGene Set:", gs, "(", nrow(sub_df), "donors)\n")
        cat("  Mean Dev_Explained_Summary:",
            round(mean(sub_df$Dev_Explained_Summary, na.rm=TRUE), 4), "\n")
        cat("  Mean Scale_Parameter:      ",
            round(mean(sub_df$Scale_Parameter, na.rm=TRUE), 6), "\n")
        cat("  Mean EDF_Total:            ",
            round(mean(sub_df$EDF_Total, na.rm=TRUE), 3), "\n")
    }
}

cat("\n========== Part 3 complete ==========\n")
cat("Output: AllenMTG_CARM1P1_GAM_Results_3.xlsx\n")



##############################################################################
# Allen MTG CARM1P1 - Part 4: Inflection Point (IP) Detection
#                              DFG and HKG gene sets
##############################################################################
library(ggplot2)
library(dplyr)
library(mgcv)
library(openxlsx)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data    <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results <- rds_data$all_results
donors      <- rds_data$donors

valid_results <- Filter(Negate(is.null), all_results)
gene_sets     <- c("DFG", "HKG")

cat("Valid fits loaded:", length(valid_results), "\n")
cat("Donors:", paste(donors, collapse = ", "), "\n\n")

# Gene set labels and members (adapt member lists to your actual gene sets)
GS_LABELS <- list(
    DFG = "DFG composite",
    HKG = "HKG composite"
)
GS_MEMBERS <- list(
    DFG = "CLMN, EPHA3, EPHA6, LOC101928964, ROBO2",   # update as needed
    HKG = "ACTB, GAPDH, PPIA, RPL13A, TBP"       # update as needed
)

# ------------------------------------------------------------------
# SECTION 2: Tuning parameters
# ------------------------------------------------------------------
STAT_ALPHA        <- 0.05
SHAPIRO_ALPHA     <- 0.05
POST_IP_MIN_CELLS <- 5L
MIN_PURPLE_CELLS  <- 10L

IPRS_VERY_STRONG <- 0.80
IPRS_STRONG      <- 0.60
IPRS_MODERATE    <- 0.40
IPRS_WEAK        <- 0.20

# ------------------------------------------------------------------
# SECTION 3: Helper functions
# ------------------------------------------------------------------

# ---- 3a. SW-gated two-group stat test ----
stat_test_zc <- function(pre_r, post_r) {
    out <- list(stat_test   = NA_character_,
                stat_p_raw  = NA_real_,
                stat_p_adj  = NA_real_,
                both_normal = NA)
    
    if (length(pre_r) < 3L || length(post_r) < POST_IP_MIN_CELLS)
        return(out)
    
    sw_pre  <- if (length(pre_r)  >= 3L)
        suppressWarnings(shapiro.test(pre_r))$p.value  else NA_real_
    sw_post <- if (length(post_r) >= 3L)
        suppressWarnings(shapiro.test(post_r))$p.value else NA_real_
    both_normal <- !is.na(sw_pre) && !is.na(sw_post) &&
        sw_pre > SHAPIRO_ALPHA && sw_post > SHAPIRO_ALPHA
    out$both_normal <- both_normal
    
    if (isTRUE(both_normal)) {
        tt <- suppressWarnings(
            t.test(post_r, pre_r, alternative = "greater", var.equal = FALSE))
        out$stat_test  <- "Welch_t"
        out$stat_p_raw <- tt$p.value
    } else {
        mw <- suppressWarnings(
            wilcox.test(post_r, pre_r, alternative = "greater", exact = FALSE))
        out$stat_test  <- "MannWhitneyU"
        out$stat_p_raw <- mw$p.value
    }
    out$stat_p_adj <- p.adjust(out$stat_p_raw, method = "BH")
    out
}

# ---- 3b. Brent CI on R(x) GAM confidence bands ----
find_ci_brent <- function(residual_gam, x_seq, r_hat, r_se,
                          multiplier, near_x) {
    band     <- r_hat + multiplier * r_se
    sign_chg <- which(diff(sign(band)) != 0L)
    if (length(sign_chg) == 0L) return(NA_real_)
    
    if (length(sign_chg) > 1L)
        sign_chg <- sign_chg[which.min(abs(x_seq[sign_chg] - near_x))]
    
    bracket <- c(x_seq[sign_chg], x_seq[sign_chg + 1L])
    
    tryCatch({
        uniroot(
            f        = function(x) {
                p <- predict(residual_gam,
                             newdata = data.frame(CARM1P1 = x), se.fit = TRUE)
                as.numeric(p$fit) + multiplier * as.numeric(p$se.fit)
            },
            interval = bracket,
            tol      = .Machine$double.eps^0.5
        )$root
    }, error = function(e) {
        x1 <- bracket[1L]; x2 <- bracket[2L]
        y1 <- band[sign_chg]; y2 <- band[sign_chg + 1L]
        x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
    })
}

# ---- 3b2. Proximity fallback CI bound ----
find_ci_proximity_fallback <- function(x_seq, r_hat, r_se, multiplier,
                                       ip_value, side = c("lower", "upper")) {
    side  <- match.arg(side)
    band  <- r_hat + multiplier * r_se
    
    if (side == "lower") {
        idx <- which(x_seq <= ip_value)
    } else {
        idx <- which(x_seq >= ip_value)
    }
    if (length(idx) < 3L) return(NA_real_)
    
    band_region <- band[idx]
    if (any(diff(sign(band_region)) != 0L)) return(NA_real_)
    
    x_seq[idx[which.min(abs(band_region))]]
}

# ---- 3c. IPRS component calculation ----
compute_iprs <- function(ip_value, ip_ci_lower, ip_ci_upper,
                         stat_p_adj, x_seq, r_hat, r_upper, r_lower,
                         residual_data) {
    
    # --- C1: Precision ---
    carm1p1_range <- max(x_seq) - min(x_seq)
    
    both_ci  <- !is.na(ip_ci_lower) && !is.na(ip_ci_upper)
    lower_ok <- !is.na(ip_ci_lower)
    upper_ok <- !is.na(ip_ci_upper)
    
    ci_width <- if (both_ci) abs(ip_ci_upper - ip_ci_lower) else NA_real_
    
    C1 <- if (both_ci && carm1p1_range > 0) {
        max(0, 1 - ci_width / carm1p1_range)
    } else if ((lower_ok || upper_ok) && !is.na(ip_value) && carm1p1_range > 0) {
        half_width <- if (lower_ok) abs(ip_value - ip_ci_lower)
        else          abs(ip_ci_upper - ip_value)
        0.5 * max(0, 1 - half_width / carm1p1_range)
    } else {
        0
    }
    
    C1_note <- if (both_ci)
        sprintf("Both CI bounds: CI_width=%.4f / CARM1P1_range=%.4f → C1=%.4f",
                ci_width, carm1p1_range, C1)
    else if (lower_ok || upper_ok) {
        known_bound <- if (lower_ok) ip_ci_lower else ip_ci_upper
        bound_label <- if (lower_ok) "lower" else "upper"
        half_w      <- abs(ip_value - known_bound)
        sprintf("Partial CI (%s only): half_width=%.4f / CARM1P1_range=%.4f → 0.5×max(0,1-ratio)=%.4f",
                bound_label, half_w, carm1p1_range, C1)
    } else
        "Both CI bounds NA → C1=0"
    
    # --- C2: Statistical ---
    GW_ALPHA  <- 0.05
    scale_c2  <- atanh(0.5) / (-log10(GW_ALPHA))
    C2 <- if (!is.na(stat_p_adj))
        tanh(scale_c2 * (-log10(max(stat_p_adj, 1e-300))))
    else
        0
    
    C2_note <- if (is.na(stat_p_adj))
        "stat_p_adj not available → C2=0"
    else
        sprintf("tanh(%.4f * -log10(%.4e)) → C2=%.4f",
                scale_c2, stat_p_adj, C2)
    
    # --- C3: Abruptness ---
    deriv_at_ip <- if (!is.na(ip_value) && length(x_seq) >= 3L) {
        i_ip <- which.min(abs(x_seq - ip_value))
        n    <- length(x_seq)
        if (i_ip > 1L && i_ip < n) {
            abs((r_hat[i_ip + 1L] - r_hat[i_ip - 1L]) /
                    (x_seq[i_ip + 1L] - x_seq[i_ip - 1L]))
        } else if (i_ip == 1L) {
            abs((r_hat[2L] - r_hat[1L]) / (x_seq[2L] - x_seq[1L]))
        } else {
            abs((r_hat[n] - r_hat[n - 1L]) / (x_seq[n] - x_seq[n - 1L]))
        }
    } else NA_real_
    
    C3 <- if (!is.na(deriv_at_ip)) tanh(deriv_at_ip) else 0
    
    C3_note <- if (is.na(deriv_at_ip))
        "IP not available → C3=0"
    else
        sprintf("|R'(IP)|=%.4f → tanh(%.4f) → C3=%.4f",
                deriv_at_ip, deriv_at_ip, C3)
    
    # --- C4: GAM ribbon width penalty ---
    amplitude   <- diff(range(r_hat, na.rm = TRUE))
    mean_ribbon <- mean(r_upper - r_lower, na.rm = TRUE)
    C4          <- exp(-mean_ribbon / (amplitude + 1e-6))
    
    C4_note <- sprintf(
        "mean_ribbon=%.4f / amplitude=%.4f → exp(-ratio=%.4f) → C4=%.4f",
        mean_ribbon, amplitude, mean_ribbon / (amplitude + 1e-6), C4)
    
    # --- C5: R(x) smoothness + post-IP sign fidelity ---
    tv <- sum(abs(diff(r_hat)), na.rm = TRUE)
    
    if (!is.na(ip_value) && nrow(residual_data) > 0L) {
        post_cells    <- residual_data$residual[residual_data$CARM1P1 >= ip_value]
        n_post        <- length(post_cells)
        neg_post_frac <- if (n_post > 0L) sum(post_cells < 0) / n_post else 0
    } else {
        n_post        <- 0L
        neg_post_frac <- 0
    }
    sign_fidelity <- (1 - neg_post_frac)^2
    
    if (amplitude < 1e-6) {
        C5      <- 0
        C5_note <- sprintf(
            "amplitude=%.6f < 1e-6 (flat R(x)) → C5=0", amplitude)
    } else {
        excess  <- tv / (amplitude + 1e-6) - 1
        wave_sc <- exp(-max(0, excess))
        C5      <- wave_sc * sign_fidelity
        C5_note <- sprintf(
            "TV=%.4f/amp=%.4f → wave=%.4f | neg_post=%.3f → sign=(1-%.3f)^2=%.4f | C5=%.4f*%.4f=%.4f",
            tv, amplitude, wave_sc, neg_post_frac, neg_post_frac,
            sign_fidelity, wave_sc, sign_fidelity, C5)
    }
    
    IPRS <- (C1 + C2 + C3 + C4 + C5) / 5
    
    # Hard rule: any component = 0 → check for structural incompleteness
    any_zero <- any(c(C1, C2, C3, C4, C5) == 0)
    
    IPRS_tier <- if (IPRS >= IPRS_VERY_STRONG) "Very Strong" else
        if (IPRS >= IPRS_STRONG)      "Strong"      else
            if (IPRS >= IPRS_MODERATE)    "Moderate"    else
                if (IPRS >= IPRS_WEAK)        "Weak"        else
                    "Failed"
    
    IPRS_flag <- if (any_zero) "Visual review needed [C=0]" else
        if (IPRS_tier %in% c("Very Strong", "Strong")) "None" else
            if (IPRS_tier %in% c("Moderate", "Weak"))  "Visual review needed" else
                "No IP - do not proceed"
    
    cat(sprintf("  IPRS components: C1=%.4f [%s]\n", C1, C1_note))
    cat(sprintf("                   C2=%.4f [%s]\n", C2, C2_note))
    cat(sprintf("                   C3=%.4f [%s]\n", C3, C3_note))
    cat(sprintf("                   C4=%.4f [%s]\n", C4, C4_note))
    cat(sprintf("                   C5=%.4f [%s]\n", C5, C5_note))
    cat(sprintf("  IPRS = (%.4f+%.4f+%.4f+%.4f+%.4f)/5 = %.4f  [%s]  %s\n",
                C1, C2, C3, C4, C5, IPRS, IPRS_tier,
                if (IPRS_flag == "None") "" else paste0("*** ", IPRS_flag, " ***")))
    
    list(C1          = C1,        C2          = C2,        C3          = C3,
         C4          = C4,        C5          = C5,
         IPRS        = IPRS,      IPRS_tier   = IPRS_tier, IPRS_flag   = IPRS_flag,
         C1_note     = C1_note,   C2_note     = C2_note,   C3_note     = C3_note,
         C4_note     = C4_note,   C5_note     = C5_note,
         deriv_at_ip = deriv_at_ip,
         ci_width    = ci_width,
         carm1p1_range = carm1p1_range,
         amplitude   = amplitude,
         mean_ribbon = mean_ribbon,
         tv          = tv)
}

# ---- 3d. Core IP detection function ----
detect_ip <- function(sample_data, gam_model, dev_explained_pct) {
    
    # Step 1: Purple-cell classification (top EP cells)
    null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data))
    model_fitted <- fitted(gam_model)
    
    # Guard against division by zero
    null_sq   <- (sample_data$Expression - null_fitted)^2
    model_sq  <- (sample_data$Expression - model_fitted)^2
    expl_power <- ifelse(null_sq > 0, 1 - model_sq / null_sq, 0)
    
    n_purple  <- round(nrow(sample_data) * dev_explained_pct)
    n_purple  <- max(1L, min(n_purple, nrow(sample_data)))
    is_purple <- rep(FALSE, nrow(sample_data))
    is_purple[order(expl_power, decreasing = TRUE)[seq_len(n_purple)]] <- TRUE
    
    # Step 2: Signed residuals for purple cells
    residual      <- sample_data$Expression -
        predict(gam_model, newdata = sample_data)
    purple_idx    <- which(is_purple)
    residual_data <- data.frame(
        CARM1P1  = sample_data$CARM1P1[purple_idx],
        residual = residual[purple_idx]
    )
    
    cat("  Purple cells:", nrow(residual_data),
        "| Above:", sum(residual_data$residual > 0),
        "| Below:", sum(residual_data$residual < 0), "\n")
    
    # FAILED_list: returned when IP detection cannot proceed
    FAILED_list <- function(status) {
        iprs_fail <- list(C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0,
                          IPRS = 0, IPRS_tier = "Failed",
                          IPRS_flag = "No IP - do not proceed",
                          C1_note = "IP detection failed",
                          C2_note = "IP detection failed",
                          C3_note = "IP detection failed",
                          C4_note = "IP detection failed",
                          C5_note = "IP detection failed",
                          deriv_at_ip   = NA_real_,
                          ci_width      = NA_real_,
                          carm1p1_range = NA_real_,
                          amplitude     = NA_real_,
                          mean_ribbon   = NA_real_,
                          tv            = NA_real_)
        list(ip_value = NA_real_, manual_flag = TRUE, status = status,
             stat_test = NA_character_, stat_p_raw = NA_real_,
             stat_p_adj = NA_real_, both_normal = NA,
             ip_ci_lower = NA_real_, ip_ci_upper = NA_real_,
             ci_lower_method = NA_character_, ci_upper_method = NA_character_,
             neg_pos_count = 0L,
             iprs = iprs_fail,
             residual_data = residual_data,
             residual_pred = data.frame(CARM1P1 = numeric(0),
                                        r_hat   = numeric(0),
                                        se      = numeric(0),
                                        r_upper = numeric(0),
                                        r_lower = numeric(0)))
    }
    
    if (nrow(residual_data) < MIN_PURPLE_CELLS ||
        length(unique(residual_data$CARM1P1)) < 3L) {
        cat("  *** INSUFFICIENT PURPLE CELLS - skipping IP detection ***\n")
        return(FAILED_list("insufficient_purple_cells"))
    }
    
    # Step 3: Residual GAM R(x) - unweighted (purple cells only)
    set.seed(123)
    n_unique     <- length(unique(residual_data$CARM1P1))
    k_resid      <- min(20L, n_unique - 1L)
    residual_gam <- tryCatch(
        gam(residual ~ s(CARM1P1, bs = "tp", k = k_resid),
            data = residual_data, method = "REML"),
        error = function(e) {
            cat("  Fallback k=8:", conditionMessage(e), "\n")
            gam(residual ~ s(CARM1P1, bs = "tp",
                             k = min(8L, n_unique - 1L)),
                data = residual_data, method = "REML")
        }
    )
    gam_sum <- summary(residual_gam)
    cat("  R(x) GAM: edf =", sprintf("%.2f", gam_sum$edf),
        "| Dev =", sprintf("%.1f%%", gam_sum$dev.expl * 100),
        "| p =", sprintf("%.2e", gam_sum$s.table[4]), "\n")
    
    # Step 4: 2000-pt grid prediction
    x_seq <- seq(min(residual_data$CARM1P1), max(residual_data$CARM1P1),
                 length.out = 2000L)
    pred  <- predict(residual_gam,
                     newdata = data.frame(CARM1P1 = x_seq), se.fit = TRUE)
    r_hat <- as.numeric(pred$fit)
    r_se  <- as.numeric(pred$se.fit)
    
    r_upper <- r_hat + 1.96 * r_se
    r_lower <- r_hat - 1.96 * r_se
    
    # Step 5: First neg→pos zero-crossing (linear interpolation)
    sc_all <- which(diff(sign(r_hat)) != 0L)
    sc_np  <- sc_all[r_hat[sc_all] < 0 & r_hat[sc_all + 1L] > 0]
    
    ip_value    <- NA_real_
    manual_flag <- FALSE
    
    if (length(sc_np) > 0L) {
        i_cross  <- sc_np[1L]
        x0 <- x_seq[i_cross];      y0 <- r_hat[i_cross]
        x1 <- x_seq[i_cross + 1L]; y1 <- r_hat[i_cross + 1L]
        ip_value <- x0 - y0 * (x1 - x0) / (y1 - y0)
        cat("  IP (R(x) zero-crossing):", sprintf("%.6f", ip_value),
            "| neg\u2192pos crossings:", length(sc_np), "\n")
    } else {
        manual_flag <- TRUE
        cat("  *** R(x) has no neg\u2192pos zero-crossing",
            "- R(x) is entirely",
            ifelse(all(r_hat >= 0), "positive", "negative"), "***\n")
    }
    
    # Step 6: Brent 95% CI - with proximity fallback for missing bounds
    ci_anchor <- if (!is.na(ip_value)) ip_value else median(x_seq)
    
    ip_ci_lower     <- tryCatch(
        find_ci_brent(residual_gam, x_seq, r_hat, r_se, +1.96, ci_anchor),
        error = function(e) NA_real_)
    ci_lower_method <- if (!is.na(ip_ci_lower)) "brent" else NA_character_
    
    ip_ci_upper     <- tryCatch(
        find_ci_brent(residual_gam, x_seq, r_hat, r_se, -1.96, ci_anchor),
        error = function(e) NA_real_)
    ci_upper_method <- if (!is.na(ip_ci_upper)) "brent" else NA_character_
    
    if (is.na(ip_ci_lower) && !is.na(ip_value)) {
        fb_lower <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                               +1.96, ip_value, side = "lower")
        if (!is.na(fb_lower)) {
            ip_ci_lower     <- fb_lower
            ci_lower_method <- "proximity_fallback"
            cat("  CI lower: Brent NA → proximity fallback =",
                sprintf("%.4f", ip_ci_lower), "[band closest-to-zero pre-IP]\n")
        }
    }
    if (is.na(ip_ci_upper) && !is.na(ip_value)) {
        fb_upper <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                               -1.96, ip_value, side = "upper")
        if (!is.na(fb_upper)) {
            ip_ci_upper     <- fb_upper
            ci_upper_method <- "proximity_fallback"
            cat("  CI upper: Brent NA → proximity fallback =",
                sprintf("%.4f", ip_ci_upper), "[band closest-to-zero post-IP]\n")
        }
    }
    
    if (!is.na(ip_ci_lower) && !is.na(ip_ci_upper) &&
        ip_ci_lower > ip_ci_upper) {
        tmp              <- ip_ci_lower
        ip_ci_lower      <- ip_ci_upper
        ip_ci_upper      <- tmp
        tmp_m            <- ci_lower_method
        ci_lower_method  <- ci_upper_method
        ci_upper_method  <- tmp_m
    }
    
    cat("  95% CI: [",
        if (is.na(ip_ci_lower)) "NA" else sprintf("%.4f [%s]", ip_ci_lower, ci_lower_method),
        ",",
        if (is.na(ip_ci_upper)) "NA" else sprintf("%.4f [%s]", ip_ci_upper, ci_upper_method),
        "] | width:",
        if (is.na(ip_ci_lower) || is.na(ip_ci_upper)) "NA"
        else sprintf("%.4f", abs(ip_ci_upper - ip_ci_lower)), "\n")
    
    # Step 7: Stat test pre vs post zero-crossing
    if (!is.na(ip_value)) {
        pre_r  <- residual_data$residual[residual_data$CARM1P1 <  ip_value]
        post_r <- residual_data$residual[residual_data$CARM1P1 >= ip_value]
        st     <- stat_test_zc(pre_r, post_r)
        cat("  Stat:", st$stat_test,
            sprintf("| p_adj: %s",
                    if (is.na(st$stat_p_adj)) "NA"
                    else sprintf("%.3e", st$stat_p_adj)), "\n")
    } else {
        st <- list(stat_test = NA_character_, stat_p_raw = NA_real_,
                   stat_p_adj = NA_real_, both_normal = NA)
    }
    
    # Step 8: Compute IPRS (C1-C5)
    iprs <- compute_iprs(ip_value    = ip_value,
                         ip_ci_lower = ip_ci_lower,
                         ip_ci_upper = ip_ci_upper,
                         stat_p_adj  = st$stat_p_adj,
                         x_seq       = x_seq,
                         r_hat       = r_hat,
                         r_upper     = r_upper,
                         r_lower     = r_lower,
                         residual_data = residual_data)
    
    list(
        ip_value        = ip_value,
        manual_flag     = manual_flag,
        status          = "ok",
        stat_test       = st$stat_test,
        stat_p_raw      = st$stat_p_raw,
        stat_p_adj      = st$stat_p_adj,
        both_normal     = st$both_normal,
        ip_ci_lower     = ip_ci_lower,
        ip_ci_upper     = ip_ci_upper,
        ci_lower_method = ci_lower_method,
        ci_upper_method = ci_upper_method,
        neg_pos_count   = length(sc_np),
        iprs            = iprs,
        residual_data   = residual_data,
        residual_pred   = data.frame(
            CARM1P1 = x_seq,
            r_hat   = r_hat,
            se      = r_se,
            r_upper = r_upper,
            r_lower = r_lower
        )
    )
}

# ------------------------------------------------------------------
# SECTION 4: Run IP detection across all donor × gene_set fits
# ------------------------------------------------------------------
set.seed(123)
cat("\n=============================================\n")
cat(sprintf("IP DETECTION | Zero-crossing + Brent CI + IPRS | STAT_ALPHA=%.2f\n",
            STAT_ALPHA))
cat("Dataset  : Allen MTG - CARM1P1 Analysis\n")
cat("Donors   :", paste(donors, collapse = ", "), "\n")
cat("Gene sets: DFG, HKG\n")
cat("=============================================\n\n")

ip_results <- list()

for (key in names(valid_results)) {
    parts          <- strsplit(key, "\\|")[[1]]
    donor_id       <- parts[1]
    gene_set_label <- parts[2]
    
    cat("\n--- Processing:", donor_id, "|",
        GS_LABELS[[gene_set_label]], "---\n")
    
    res         <- valid_results[[key]]
    sample_data <- res$data
    gam_model   <- res$best_model
    dev_expl    <- summary(gam_model)$dev.expl
    
    ip_results[[key]] <- tryCatch(
        detect_ip(sample_data, gam_model, dev_expl),
        error = function(e) {
            cat("  ERROR in detect_ip:", conditionMessage(e), "\n")
            NULL
        }
    )
}

# ------------------------------------------------------------------
# SECTION 5: Console summary
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("IP SUMMARY - IP Reliability Score (IPRS)\n")
cat("=============================================\n\n")
cat(sprintf("%-30s %8s %8s %8s %6s %6s %6s %6s %6s %6s  %-14s  %s\n",
            "Fit (Donor|GS)", "IP", "CI_Lo", "CI_Hi",
            "C1", "C2", "C3", "C4", "C5", "IPRS", "Tier", "Flag"))
cat(strrep("-", 140), "\n")

for (key in names(ip_results)) {
    r <- ip_results[[key]]; if (is.null(r)) next
    iprs <- r$iprs
    cat(sprintf("%-30s %8s %8s %8s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  %-14s  %s\n",
                key,
                if (is.na(r$ip_value))    "NA" else sprintf("%.4f", r$ip_value),
                if (is.na(r$ip_ci_lower)) "NA" else sprintf("%.4f", r$ip_ci_lower),
                if (is.na(r$ip_ci_upper)) "NA" else sprintf("%.4f", r$ip_ci_upper),
                iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5, iprs$IPRS,
                iprs$IPRS_tier, iprs$IPRS_flag))
}

cat("\n")
for (key in names(ip_results)) {
    r <- ip_results[[key]]; if (is.null(r)) next
    if (r$iprs$IPRS_flag != "None")
        cat(sprintf("  [FLAG] %s | IPRS=%.3f (%s) → %s\n",
                    key, r$iprs$IPRS, r$iprs$IPRS_tier, r$iprs$IPRS_flag))
}

# ------------------------------------------------------------------
# SECTION 6: Build ip_table and extract environment vectors
# ------------------------------------------------------------------
ip_table <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]
    parts <- strsplit(key, "\\|")[[1]]
    if (is.null(r)) {
        return(data.frame(
            Donor = parts[1], Gene_Set = parts[2],
            IP_Value = NA_real_, IP_CI_Lower = NA_real_, IP_CI_Upper = NA_real_,
            CI_Width = NA_real_,
            C1 = NA_real_, C2 = NA_real_, C3 = NA_real_,
            C4 = NA_real_, C5 = NA_real_,
            IPRS = NA_real_, IPRS_Tier = NA_character_, IPRS_Flag = NA_character_,
            stringsAsFactors = FALSE))
    }
    iprs <- r$iprs
    data.frame(
        Donor       = parts[1],
        Gene_Set    = parts[2],
        IP_Value    = r$ip_value,
        IP_CI_Lower = r$ip_ci_lower,
        IP_CI_Upper = r$ip_ci_upper,
        CI_Width    = iprs$ci_width,
        C1          = iprs$C1,
        C2          = iprs$C2,
        C3          = iprs$C3,
        C4          = iprs$C4,
        C5          = iprs$C5,
        IPRS        = iprs$IPRS,
        IPRS_Tier   = iprs$IPRS_tier,
        IPRS_Flag   = iprs$IPRS_flag,
        stringsAsFactors = FALSE
    )
}))
rownames(ip_table) <- NULL

# ---- Environment vectors for Part 5 ----
inflection_points <- setNames(
    ip_table$IP_Value,
    paste0(ip_table$Donor, "|", ip_table$Gene_Set)
)

# Per-gene-set convenience vectors
inflection_points_DFG <- setNames(
    ip_table$IP_Value[ip_table$Gene_Set == "DFG"],
    ip_table$Donor[ip_table$Gene_Set == "DFG"]
)

inflection_points_HKG <- setNames(
    ip_table$IP_Value[ip_table$Gene_Set == "HKG"],
    ip_table$Donor[ip_table$Gene_Set == "HKG"]
)

# Mirrors inflection_points_by_patient from GBM template (DFG-keyed by donor)
inflection_points_by_donor <- setNames(
    sapply(donors, function(d) inflection_points_DFG[d]),
    donors
)

cat("\n=== Environment vectors created ===\n")
cat("inflection_points           (", length(inflection_points),
    "entries, keyed 'donor|gene_set') - used by Part 5\n")
cat("inflection_points_DFG       (", length(inflection_points_DFG),
    "entries, keyed by donor)\n")
cat("inflection_points_HKG       (", length(inflection_points_HKG),
    "entries, keyed by donor)\n")
cat("inflection_points_by_donor  (", length(inflection_points_by_donor),
    "entries, DFG only, keyed by donor)\n\n")

cat("Final IP values with IPRS:\n")
print(ip_table[, c("Donor", "Gene_Set", "IP_Value",
                   "IP_CI_Lower", "IP_CI_Upper",
                   "C1", "C2", "C3", "C4", "C5", "IPRS", "IPRS_Tier", "IPRS_Flag")])

# ------------------------------------------------------------------
# SECTION 7: Excel export (3 sheets)
# ------------------------------------------------------------------
cat("\nExporting results...\n")

ip_sheet <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]
    parts <- strsplit(key, "\\|")[[1]]
    if (is.null(r)) return(NULL)
    iprs  <- r$iprs
    gs    <- parts[2]
    data.frame(
        Donor            = parts[1],
        Group            = "Allen_MTG",
        Gene_Set         = gs,
        Gene_Set_Label   = GS_LABELS[[gs]],
        Gene_Members     = GS_MEMBERS[[gs]],
        IP_Value         = r$ip_value,
        CI_Lower_95      = r$ip_ci_lower,
        CI_Upper_95      = r$ip_ci_upper,
        CI_Width         = iprs$ci_width,
        CARM1P1_Range    = iprs$carm1p1_range,
        CI_Lower_Method  = if (is.null(r$ci_lower_method)) NA_character_
        else r$ci_lower_method,
        CI_Upper_Method  = if (is.null(r$ci_upper_method)) NA_character_
        else r$ci_upper_method,
        NegPos_Crossings = r$neg_pos_count,
        N_Purple_Cells   = nrow(r$residual_data),
        Stat_Test        = r$stat_test,
        Stat_P_Raw       = r$stat_p_raw,
        Stat_P_Adj       = r$stat_p_adj,
        C1_Precision     = iprs$C1,
        C1_Note          = iprs$C1_note,
        C2_Statistical   = iprs$C2,
        C2_Note          = iprs$C2_note,
        C3_Abruptness    = iprs$C3,
        C3_Note          = iprs$C3_note,
        C4_GAM_Ribbon    = iprs$C4,
        C4_Note          = iprs$C4_note,
        C5_Smoothness    = iprs$C5,
        C5_Note          = iprs$C5_note,
        Rx_Amplitude     = iprs$amplitude,
        Rx_MeanRibbon    = iprs$mean_ribbon,
        Rx_TotalVar      = iprs$tv,
        Deriv_at_IP      = iprs$deriv_at_ip,
        IPRS             = iprs$IPRS,
        IPRS_Tier        = iprs$IPRS_tier,
        IPRS_Flag        = iprs$IPRS_flag,
        Method           = "ZeroCrossing|BrentCI|IPRS=(C1+C2+C3+C4+C5)/5|C2=tanh(-log10p)|C5=waviness*sign_fidelity|SW_welch_t_OR_mannwhitney_u",
        stringsAsFactors = FALSE
    )
}))
rownames(ip_sheet) <- NULL

resid_pred_df <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]; if (is.null(r)) return(NULL)
    if (nrow(r$residual_pred) == 0L) return(NULL)
    parts <- strsplit(key, "\\|")[[1]]
    # Thin to every 10th row for file size
    p     <- r$residual_pred[seq(1L, nrow(r$residual_pred), by = 10L), ]
    data.frame(Donor = parts[1], Gene_Set = parts[2], p,
               stringsAsFactors = FALSE)
}))
if (is.null(resid_pred_df) || length(resid_pred_df) == 0)
    resid_pred_df <- data.frame()
rownames(resid_pred_df) <- NULL

method_notes <- data.frame(
    Step = paste0("Step ", 1:8),
    Description = c(
        "Cell classification: EP = 1 - model_residual² / null_residual²; top (dev_expl × n) cells are 'purple'. Null model fitted unweighted (no dynamic zero-weighting in Allen MTG Part 1).",
        "Signed residuals: residual = observed Expression − GAM_fitted, for each purple cell.",
        "Residual GAM R(x): thin-plate spline (REML) fitted to (CARM1P1, residual) for purple cells; k = min(20, n_unique−1). Unweighted - all purple cells contribute equally to R(x).",
        "Fine grid: R(x) and SE evaluated at 2000 equally-spaced CARM1P1 values.",
        paste0("IP = first neg→pos zero-crossing of R(x) GAM (linear interpolation). ",
               "If R(x) has no neg→pos crossing, IP=NA and all IPRS components default to 0."),
        paste0("95% CI primary (Brent): CI_Lower = zero of [R_hat + 1.96×SE]; CI_Upper = zero of [R_hat - 1.96×SE]; anchored near IP. ",
               "Fallback (proximity): if Brent returns NA for one bound, ",
               "fallback finds x where band is closest to zero in pre-IP (lower) or post-IP (upper) region. ",
               "CI_Lower_Method / CI_Upper_Method column records 'brent' or 'proximity_fallback'."),
        paste0("Stat test: Shapiro-Wilk (α=", SHAPIRO_ALPHA,
               ") per group; both normal → Welch t (one-sided post>pre); ",
               "either non-normal → Mann-Whitney U; BH FDR applied."),
        paste0("IPRS = (C1 + C2 + C3 + C4 + C5) / 5  [range 0-1]. ",
               "C1 (Precision): both CI bounds → max(0, 1-CI_Width/CARM1P1_range); ",
               "one bound NA → 0.5×max(0, 1-half_width/CARM1P1_range); both NA → 0. ",
               "C2 (Statistical) = tanh(scale_c2 × (-log10(p_adj))); scale_c2 = atanh(0.5)/(-log10(0.05)); C2=0.5 at p=0.05. ",
               "C3 (Abruptness) = tanh(|R'(IP)|) via central difference on 2000-pt grid. ",
               "C4 (GAM Ribbon) = exp(-mean_ribbon_width / amplitude). ",
               "C5 (Smoothness) = exp(-(TV/amplitude - 1)) × (1-neg_post_frac)²; C5=0 if amplitude < 1e-6. ",
               "Tiers: ≥0.80 Very Strong | 0.60-0.79 Strong | 0.40-0.59 Moderate [visual review] | ",
               "0.20-0.39 Weak [visual review] | <0.20 Failed [no IP - do not proceed]. ",
               "IP/CI lines plotted in green (#2E8B57) regardless of IPRS tier.")
    ),
    stringsAsFactors = FALSE
)

hdr_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#1F3864",
                         textDecoration = "bold", fontSize = 10,
                         halign = "center", valign = "center",
                         border = "TopBottomLeftRight",
                         borderStyle = "thin", borderColour = "#AAAAAA",
                         wrapText = TRUE)

wb <- createWorkbook()

addWorksheet(wb, "IP_Results")
writeData(wb, "IP_Results", ip_sheet)
addStyle(wb, "IP_Results", hdr_style,
         rows = 1, cols = seq_len(ncol(ip_sheet)), gridExpand = TRUE)
setColWidths(wb, "IP_Results", cols = 1:2,   widths = 16)
setColWidths(wb, "IP_Results", cols = 3:5,   widths = 20)
setColWidths(wb, "IP_Results", cols = 6:17,  widths = 16)
setColWidths(wb, "IP_Results", cols = 18:31, widths = 16)
setColWidths(wb, "IP_Results", cols = 32,    widths = 90)
setRowHeights(wb, "IP_Results", rows = 1, heights = 28)
freezePane(wb, "IP_Results", firstActiveRow = 2L, firstActiveCol = 3L)

addWorksheet(wb, "Residual_Predictions")
writeData(wb, "Residual_Predictions", resid_pred_df)
addStyle(wb, "Residual_Predictions", hdr_style,
         rows = 1, cols = seq_len(ncol(resid_pred_df)), gridExpand = TRUE)
setColWidths(wb, "Residual_Predictions",
             cols = seq_len(ncol(resid_pred_df)), widths = 16)
setRowHeights(wb, "Residual_Predictions", rows = 1, heights = 28)

addWorksheet(wb, "Method_Steps")
writeData(wb, "Method_Steps", method_notes)
addStyle(wb, "Method_Steps", hdr_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Method_Steps", cols = 1, widths = 10)
setColWidths(wb, "Method_Steps", cols = 2, widths = 130)
setRowHeights(wb, "Method_Steps", rows = 1, heights = 28)

saveWorkbook(wb, "AllenMTG_CARM1P1_GAM_Results_4.xlsx", overwrite = TRUE)
cat("Exported: AllenMTG_CARM1P1_GAM_Results_4.xlsx\n")

# ------------------------------------------------------------------
# SECTION 8: Plot colour constants
# ------------------------------------------------------------------
ip_line_color <- "#2E8B57"
ip_ci_color   <- "#2E8B57"
manual_color  <- "#CC0000"

iprs_subtitle <- function(r) {
    iprs <- r$iprs
    sprintf("IPRS=%.3f (%s)  |  C1=%.3f  C2=%.3f  C3=%.3f  C4=%.3f  C5=%.3f  |  %s",
            iprs$IPRS, iprs$IPRS_tier,
            iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5,
            if (iprs$IPRS_flag == "None") "OK" else iprs$IPRS_flag)
}

# ------------------------------------------------------------------
# SECTION 9: Scatter plot (CARM1P1 vs gene-set Expression)
# Purple = dev-explained cells; gray = non-dev-explained
# ------------------------------------------------------------------
create_scatter <- function(key) {
    r <- ip_results[[key]]; if (is.null(r)) return(NULL)
    parts      <- strsplit(key, "\\|")[[1]]
    donor_id   <- parts[1]
    gs         <- parts[2]
    gs_label   <- GS_LABELS[[gs]]
    
    res          <- valid_results[[key]]
    sample_data  <- res$data
    gam_model    <- res$best_model
    dev_expl     <- summary(gam_model)$dev.expl
    
    # Null model (unweighted) for EP classification
    null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data))
    model_fitted <- fitted(gam_model)
    null_sq      <- (sample_data$Expression - null_fitted)^2
    model_sq     <- (sample_data$Expression - model_fitted)^2
    expl_power   <- ifelse(null_sq > 0, 1 - model_sq / null_sq, 0)
    
    n_dev      <- round(nrow(sample_data) * dev_expl)
    dev_cells  <- rownames(sample_data)[
        order(expl_power, decreasing = TRUE)[seq_len(n_dev)]]
    
    plot_data <- data.frame(
        CARM1P1    = sample_data$CARM1P1,
        Expression = sample_data$Expression,
        Group      = ifelse(rownames(sample_data) %in% dev_cells,
                            "Dev explained", "Non-dev explained")
    )
    plot_data <- plot_data[order(plot_data$Group == "Dev explained"), ]
    
    pred_data        <- data.frame(CARM1P1 = seq(min(plot_data$CARM1P1),
                                                 max(plot_data$CARM1P1),
                                                 length.out = 100L))
    pg               <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit    <- pg$fit
    pred_data$se.fit <- pg$se.fit
    
    lc        <- ip_line_color
    title_col <- if (r$manual_flag) manual_color else "black"
    
    title <- if (r$manual_flag || is.na(r$ip_value)) {
        paste0(donor_id, " | ", gs_label, "  *** NO IP DETECTED ***")
    } else {
        paste0(donor_id, " | ", gs_label,
               "  IP: ", sprintf("%.4f", r$ip_value),
               "  95% CI [",
               ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
               ", ",
               ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
               "]")
    }
    
    p <- ggplot()
    
    # Soft green CI band (mirrors create_rx) - drawn first so it sits behind dots
    if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
        p <- p + annotate("rect",
                          xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                          ymin = -Inf, ymax = Inf,
                          fill = lc, alpha = 0.08)
    
    p <- p +
        geom_point(data = plot_data,
                   aes(x = CARM1P1, y = Expression, color = Group),
                   size = 1.8, alpha = 0.50) +
        geom_ribbon(data = pred_data,
                    aes(x = CARM1P1, ymin = fit - 1.96 * se.fit,
                        ymax = fit + 1.96 * se.fit),
                    fill = "#FFCC99", alpha = 0.20) +
        geom_line(data = pred_data, aes(x = CARM1P1, y = fit),
                  color = "#FFCC99", linewidth = 1.2)
    
    if (!is.na(r$ip_value))
        p <- p + geom_vline(xintercept = r$ip_value,
                            linetype = "dashed", color = lc, linewidth = 0.8)
    if (!is.na(r$ip_ci_lower))
        p <- p + geom_vline(xintercept = r$ip_ci_lower,
                            linetype = "dotdash", color = lc, linewidth = 0.6)
    if (!is.na(r$ip_ci_upper))
        p <- p + geom_vline(xintercept = r$ip_ci_upper,
                            linetype = "dotdash", color = lc, linewidth = 0.6)
    
    p + scale_color_manual(
        values = c("Dev explained"     = "#4B0082",
                   "Non-dev explained" = "#C0C0C0")) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position  = "none",
            plot.title       = element_text(size = 10, face = "bold",
                                            hjust = 0.5, color = title_col),
            plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                            color = "black"),
            panel.border     = element_rect(color = "black", fill = NA,
                                            linewidth = 0.5),
            axis.line        = element_blank()
        ) +
        labs(title    = title,
             subtitle = paste0("Dev explained: ", round(dev_expl * 100, 2),
                               "%\n", iprs_subtitle(r)),
             x = "CARM1P1 Expression (log2)",
             y = paste0(gs_label, " Mean Expression (log2)")) +
        scale_x_continuous(
            limits = range(sample_data$CARM1P1, na.rm = TRUE) +
                c(-1, 1) * 0.05 * diff(range(sample_data$CARM1P1, na.rm = TRUE))
        )
}

# ------------------------------------------------------------------
# SECTION 10: R(x) diagnostic plot
# ------------------------------------------------------------------
create_rx <- function(key) {
    r <- ip_results[[key]]; if (is.null(r)) return(NULL)
    
    parts      <- strsplit(key, "\\|")[[1]]
    donor_id   <- parts[1]
    gs         <- parts[2]
    gs_label   <- GS_LABELS[[gs]]
    
    pred      <- r$residual_pred
    rd        <- r$residual_data
    iprs      <- r$iprs
    lc        <- ip_line_color
    title_col <- if (r$manual_flag) manual_color else "black"
    
    n_bins_rx <- 100L
    p_carm    <- rd$CARM1P1
    
    if (nrow(rd) > 0L) {
        if (length(unique(p_carm)) > n_bins_rx) {
            breaks_rx  <- seq(min(p_carm), max(p_carm),
                              length.out = n_bins_rx + 1L)
            bin_ids_rx <- cut(p_carm, breaks = breaks_rx,
                              include.lowest = TRUE, labels = FALSE)
            bin_ctr_rx <- (breaks_rx[-length(breaks_rx)] + breaks_rx[-1L]) / 2
            bin_w_rx   <- breaks_rx[2L] - breaks_rx[1L]
        } else {
            uq_rx      <- sort(unique(p_carm))
            bin_ids_rx <- match(p_carm, uq_rx)
            bin_ctr_rx <- uq_rx
            bin_w_rx   <- if (length(uq_rx) > 1L) min(diff(uq_rx)) else 0.1
        }
        
        bs_rx <- data.frame(
            carm_center   = bin_ctr_rx,
            mean_residual = sapply(seq_along(bin_ctr_rx), function(b) {
                idx <- which(bin_ids_rx == b)
                if (length(idx) > 0L) mean(rd$residual[idx]) else NA_real_
            })
        )
        bs_rx    <- bs_rx[!is.na(bs_rx$mean_residual), ]
        bar_w_rx <- if (nrow(bs_rx) > 1L)
            diff(range(bs_rx$carm_center)) / nrow(bs_rx) * 0.8
        else bin_w_rx * 0.8
    } else {
        bs_rx    <- data.frame(carm_center   = numeric(0),
                               mean_residual = numeric(0))
        bar_w_rx <- 0.1
    }
    
    stat_txt <- if (!is.na(r$stat_test)) {
        sprintf("%s | BH FDR q = %.3e", r$stat_test,
                ifelse(is.na(r$stat_p_adj), NA, r$stat_p_adj))
    } else "Stat: NA"
    
    subtitle <- if (r$status == "insufficient_purple_cells") {
        paste0("*** FAILED - R(x) GAM not fitted (",
               nrow(rd), " purple cells, ",
               length(unique(rd$CARM1P1)), " unique CARM1P1 values",
               " - minimum ", MIN_PURPLE_CELLS,
               " cells / 3 unique values required) ***")
    } else if (r$manual_flag) {
        paste0("*** NO IP - R(x) has no neg\u2192pos crossing ***\n",
               iprs_subtitle(r))
    } else {
        paste0("IP = ", sprintf("%.4f", r$ip_value),
               "  95% CI [",
               ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
               ", ",
               ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
               "]  |  ", stat_txt, "\n", iprs_subtitle(r))
    }
    
    p <- ggplot()
    
    if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
        p <- p + annotate("rect",
                          xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                          ymin = -Inf, ymax = Inf,
                          fill = lc, alpha = 0.08)
    
    if (nrow(bs_rx) > 0L)
        p <- p + geom_bar(data = bs_rx,
                          aes(x = carm_center, y = mean_residual),
                          stat = "identity", fill = "#4B0082", color = "#3D0065",
                          alpha = 0.40, width = bar_w_rx, linewidth = 0.4)
    
    if (nrow(pred) > 0L) {
        p <- p +
            geom_ribbon(data = pred,
                        aes(x = CARM1P1, ymin = r_lower, ymax = r_upper),
                        fill = "#FFCC99", alpha = 0.35) +
            geom_line(data = pred, aes(x = CARM1P1, y = r_hat),
                      color = "#FFCC99", linewidth = 1.5)
    }
    
    p <- p + geom_hline(yintercept = 0, linetype = "dotted",
                        color = "gray40", linewidth = 0.5)
    
    if (!is.na(r$ip_ci_lower))
        p <- p + geom_vline(xintercept = r$ip_ci_lower,
                            linetype = "dotdash", color = lc, linewidth = 0.7)
    if (!is.na(r$ip_ci_upper))
        p <- p + geom_vline(xintercept = r$ip_ci_upper,
                            linetype = "dotdash", color = lc, linewidth = 0.7)
    if (!is.na(r$ip_value))
        p <- p + geom_vline(xintercept = r$ip_value,
                            linetype = "dashed", color = lc, linewidth = 1.0)
    
    p + theme_minimal() +
        theme(
            panel.border  = element_rect(color = "black", fill = NA,
                                         linewidth = 0.5),
            axis.line     = element_blank(),
            plot.title    = element_text(size = 10, face = "bold",
                                         hjust = 0.5, color = title_col),
            plot.subtitle = element_text(size = 7.5, hjust = 0.5,
                                         color = "black"),
            axis.text     = element_text(size = 7)
        ) +
        labs(
            title    = paste0(donor_id, " | ", gs_label,
                              " \u2014 R(x) = E[residual | purple, CARM1P1=x]",
                              "  [Zero-crossing | Brent CI | IPRS]"),
            subtitle = subtitle,
            x = "CARM1P1 Expression (log2)",
            y = "Expected Signed Residual R(x)"
        )
}

# ------------------------------------------------------------------
# SECTION 11: Generate and print all plots
# ------------------------------------------------------------------
cat("\n\nGenerating scatter plots...\n")
for (key in names(ip_results)) {
    tryCatch({ p <- create_scatter(key); if (!is.null(p)) print(p) },
             error = function(e) cat("Scatter error:", key,
                                     conditionMessage(e), "\n"))
}

cat("\nGenerating R(x) diagnostic plots...\n")
for (key in names(ip_results)) {
    tryCatch({ p <- create_rx(key); if (!is.null(p)) print(p) },
             error = function(e) cat("R(x) error:", key,
                                     conditionMessage(e), "\n"))
}

# ------------------------------------------------------------------
# SECTION 12: IPRS tier summary
# ------------------------------------------------------------------
for (gs in gene_sets) {
    gs_keys     <- names(ip_results)[grepl(paste0("\\|", gs, "$"), names(ip_results))]
    tier_counts <- table(sapply(ip_results[gs_keys],
                                function(r) if (!is.null(r)) r$iprs$IPRS_tier else "NULL"))
    cat(sprintf("\nIPRS tier summary - %s (%d donors):\n",
                GS_LABELS[[gs]], length(gs_keys)))
    for (tier in c("Very Strong", "Strong", "Moderate", "Weak", "Failed"))
        cat(sprintf("  %-12s : %d\n", tier,
                    if (tier %in% names(tier_counts)) tier_counts[[tier]] else 0L))
}

cat("\n========== Part 4 complete ==========\n")
cat("\nOutput   : AllenMTG_CARM1P1_GAM_Results_4.xlsx\n")
cat("Env ready: inflection_points / inflection_points_DFG / inflection_points_HKG /",
    "inflection_points_by_donor\n")
cat("Part 5 will find inflection_points keyed as 'donor|gene_set'",
    "(e.g. inflection_points[\"H200.1030|DFG\"]).\n")


##############################################################################
# Allen MTG CARM1P1 - Part 5: CEP-IP GAM Scatter Plots
#                              Dev-explained cells (purple) vs rest (gray)
##############################################################################
library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data      <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results   <- rds_data$all_results
donors        <- rds_data$donors
valid_results <- Filter(Negate(is.null), all_results)
gene_sets     <- c("DFG", "HKG")

cat("Valid fits loaded:", length(valid_results), "\n")

# ------------------------------------------------------------------
# SECTION 2: Check for inflection points in environment
# (mirrors the ip_values guard in the original script)
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 4 (IP detection), which creates a named\n",
        "numeric vector 'inflection_points' keyed by \"donor|gene_set\".\n",
        "Example: inflection_points[\"H200.1030|DFG\"]"
    ))
}

# ------------------------------------------------------------------
# SECTION 3: Scatter plot function
# ------------------------------------------------------------------
create_carm1p1_scatter_plot <- function(donor, gene_set_label, result) {
    cat("\nCreating scatter plot for:", donor, "|", gene_set_label, "\n")

    gam_model           <- result$best_model
    sample_data         <- result$data
    model_dev_explained <- summary(gam_model)$dev.expl

    # ---- EP calculation (same formula as Part 3) ----
    null_model      <- gam(Expression ~ 1, data = sample_data)
    null_fitted     <- fitted(null_model)
    model_fitted    <- fitted(gam_model)

    null_residuals  <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted

    null_sq_diff  <- null_residuals^2
    model_sq_diff <- model_residuals^2

    # Guard against zero null_sq_diff (cell already at null mean)
    explanatory_power <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )

    # ---- Classify deviance vs non-deviance cells ----
    sorted_indices    <- order(explanatory_power, decreasing = TRUE)
    sorted_cell_names <- rownames(sample_data)[sorted_indices]

    target_cells <- round(nrow(sample_data) * model_dev_explained)
    target_cells <- max(1L, min(target_cells, nrow(sample_data)))  # bounds guard

    deviance_cells     <- sorted_cell_names[seq_len(target_cells)]
    non_deviance_cells <- sorted_cell_names[(target_cells + 1):nrow(sample_data)]

    # ---- Dynamic axis limits with 5% margin ----
    x_range  <- range(sample_data$CARM1P1, na.rm = TRUE)
    y_range  <- range(sample_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)

    # ---- plot_data: draw deviance cells on top (draw_order = 2) ----
    plot_data <- data.frame(
        CARM1P1   = sample_data$CARM1P1,
        Expression = sample_data$Expression,
        Group      = ifelse(rownames(sample_data) %in% deviance_cells,
                            "Dev explained", "Non-dev explained"),
        stringsAsFactors = FALSE
    )
    plot_data$draw_order <- ifelse(plot_data$Group == "Dev explained", 2, 1)
    plot_data <- plot_data[order(plot_data$draw_order), ]

    # ---- 1000-point GAM prediction with SE ----
    pred_data       <- data.frame(
        CARM1P1 = seq(x_range[1], x_range[2], length.out = 1000)
    )
    pred            <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit   <- pred$fit
    pred_data$se.fit <- pred$se.fit

    # ---- Inflection point for this donor (NA if not available) ----
    ip_val <- inflection_points[paste0(donor, "|", gene_set_label)]

    dev_pct <- round(model_dev_explained * 100, 2)

    # ---- Build plot ----
    p <- ggplot() +
        geom_point(data = plot_data,
                   aes(x = CARM1P1, y = Expression, color = Group),
                   size = 1.8, alpha = 0.6) +
        geom_line(data = pred_data,
                  aes(x = CARM1P1, y = fit),
                  color     = "#FFCC99",
                  linewidth = 1.2) +
        geom_ribbon(data = pred_data,
                    aes(x    = CARM1P1,
                        ymin = fit - 1.96 * se.fit,
                        ymax = fit + 1.96 * se.fit),
                    fill  = "#FFCC99",
                    alpha = 0.2) +

        # IP vertical dashed line (only drawn when IP is available)
        { if (!is.na(ip_val))
            geom_vline(xintercept = ip_val,
                       linetype   = "dashed",
                       color      = "black",
                       linewidth  = 0.6)
        } +

        scale_color_manual(
            values = c("Dev explained"     = "#4B0082",
                       "Non-dev explained" = "#C0C0C0")) +

        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position  = "none",
            plot.title       = element_text(size = 11, face = "bold",  hjust = 0.5),
            plot.subtitle    = element_text(size = 10, hjust = 0.5),
            axis.title       = element_text(size = 9),
            axis.text        = element_text(size = 8),
            panel.border     = element_rect(color = "black", fill = NA,
                                            linewidth = 0.5),
            axis.line        = element_blank()
        ) +

        labs(
            title    = paste0("CARM1P1 vs ", gene_set_label,
                              " Expression \u2014 ", donor),
            subtitle = paste0(
                if (!is.na(ip_val))
                    paste0("IP: ", sprintf("%.3f", ip_val), ";  ")
                else "",
                "Dev explained: ", dev_pct, "%"
            ),
            x = "CARM1P1 Expression (log2)",
            y = "Gene Set Expression (log2)"
        ) +

        scale_x_continuous(
            limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(
            limits = c(y_range[1] - y_margin, y_range[2] + y_margin))

    print(p)
    invisible(p)
}

# ------------------------------------------------------------------
# SECTION 4: Generate plots for all valid fits
# ------------------------------------------------------------------
cat("\n========== Generating CEP-IP scatter plots ==========\n")
cat("Fits to plot:", length(valid_results), "\n\n")

all_plots <- list()

for (key in names(valid_results)) {
    parts          <- strsplit(key, "\\|")[[1]]
    donor_id       <- parts[1]
    gene_set_label <- parts[2]

    tryCatch({
        p <- create_carm1p1_scatter_plot(
                 donor_id, gene_set_label, valid_results[[key]])
        all_plots[[key]] <- p
    }, error = function(e) {
        cat("  Error for", key, ":", conditionMessage(e), "\n")
    })
}

cat("\n========== Completed ==========\n")
cat("Plots generated:", length(all_plots), "of", length(valid_results), "\n")
cat("All plots displayed in RStudio Plots pane (no files saved).\n")
cat("IP values sourced from: inflection_points (produced by Part 4).\n")
cat("EP and deviance data: AllenMTG_CARM1P1_GAM_Results_3.xlsx\n")



##############################################################################
# Allen MTG CARM1P1 - Part 6: DEG Analysis
#                              Pre/Post-IP TREP cells (purple) vs non-TREP (gray)
#
# Requires:
#   - AllenMTG_CARM1P1_GAM_models.rds  (from Part 1)
#   - inflection_points                 (named numeric vector from Part 4/5,
#                                        keyed by "donor|DFG", e.g.
#                                        inflection_points["H200.1030|DFG"])
#   - seurat_obj                        (Seurat object with scRNA-seq data;
#                                        cell barcodes must match rownames in
#                                        result$data from the GAM models)
##############################################################################

library(Seurat)
library(writexl)
library(dplyr)
library(mgcv)
library(future)
library(parallelly)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data    <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results <- rds_data$all_results

# Restrict to DFG fits only for the three donors of interest
donors_of_interest <- c("H200.1023", "H200.1025", "H200.1030")
gene_set_label     <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)
valid_results <- valid_results[
    grepl(paste0("\\|", gene_set_label, "$"), names(valid_results)) &
    sapply(strsplit(names(valid_results), "\\|"), `[`, 1) %in% donors_of_interest
]

cat("DFG fits retained for analysis:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n")

# ------------------------------------------------------------------
# SECTION 2: Check prerequisites
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 4/5 (IP detection) which creates a named numeric\n",
        "vector 'inflection_points' keyed by \"donor|DFG\".\n",
        "Example: inflection_points[\"H200.1030|DFG\"]"
    ))
}

if (!exists("seurat_obj")) {
    stop(paste0(
        "seurat_obj not found in environment.\n",
        "Load your Seurat object before running Part 6.\n",
        "Cell barcodes in seurat_obj must match rownames(result$data)."
    ))
}

# ------------------------------------------------------------------
# SECTION 3: Parallelization
# ------------------------------------------------------------------
n_cores         <- max(1L, parallelly::availableCores() - 1L)
ram_limit_bytes <- 0.90 * 48 * 1024^3

plan(multisession, workers = n_cores)
options(future.globals.maxSize = ram_limit_bytes)

cat(sprintf("Parallelization: %d workers | RAM limit: %.1f GB\n",
            n_cores, ram_limit_bytes / 1024^3))

# ------------------------------------------------------------------
# SECTION 4: Configuration
# ------------------------------------------------------------------

# Default DEG thresholds
default_deg_params <- list(
    min_pct         = 0.10,
    logfc_threshold = 0.1,
    min_cells       = 3,
    p_cutoff        = 0.01,
    logfc_cutoff    = 0.2,
    max_deg         = 500
)

donor_deg_overrides <- list()

# ------------------------------------------------------------------
# SECTION 5: Helper - classify_cells
# ------------------------------------------------------------------
classify_cells <- function(donor, result, ip_value) {

    gam_model   <- result$best_model
    sample_data <- result$data     # data.frame: CARM1P1, Expression;
                                   # rownames = cell barcodes

    dev_explained <- summary(gam_model)$dev.expl

    # Per-cell explanatory power (identical to Part 5 scatter plot logic)
    null_model   <- gam(Expression ~ 1, data = sample_data)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(gam_model)

    null_sq_diff  <- (sample_data$Expression - null_fitted)^2
    model_sq_diff <- (sample_data$Expression - model_fitted)^2

    # Guard against zero null_sq_diff (cell already at null mean)
    explanatory_power <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )

    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    target_cells   <- round(nrow(sample_data) * dev_explained)
    target_cells   <- max(1L, min(target_cells, nrow(sample_data)))  # bounds guard

    is_trep <- rep(FALSE, nrow(sample_data))
    is_trep[sorted_indices[seq_len(target_cells)]] <- TRUE

    cell_barcodes  <- rownames(sample_data)
    carm1p1_values <- sample_data$CARM1P1

    classification <- data.frame(
        barcode  = cell_barcodes,
        CARM1P1  = carm1p1_values,
        is_trep  = is_trep,
        region   = ifelse(carm1p1_values < ip_value, "Pre_IP", "Post_IP"),
        stringsAsFactors = FALSE
    )

    classification$subpopulation <- paste0(
        classification$region, "_",
        ifelse(classification$is_trep, "TREP", "nonTREP")
    )

    return(classification)
}

# ------------------------------------------------------------------
# SECTION 6: Helper - run_deg_analysis
# ------------------------------------------------------------------
run_deg_analysis <- function(seurat_obj, trep_cells, nontrep_cells,
                             deg_params, region_label, donor_label) {

    cat("    DEG analysis:", region_label, "-",
        length(trep_cells), "TREP vs", length(nontrep_cells), "non-TREP\n")

    MIN_RELIABLE_CELLS <- 3L
    if (length(trep_cells)    < MIN_RELIABLE_CELLS ||
        length(nontrep_cells) < MIN_RELIABLE_CELLS) {
        cat("    Skipping: group too small (<", MIN_RELIABLE_CELLS,
            "cells) - DEGs statistically unreliable\n")
        return(NULL)
    }

    if (length(trep_cells)    < deg_params$min_cells ||
        length(nontrep_cells) < deg_params$min_cells) {
        cat("    Skipping: insufficient cells (<", deg_params$min_cells,
            "per group)\n")
        return(NULL)
    }

    Idents(seurat_obj, cells = trep_cells)    <- "TREP"
    Idents(seurat_obj, cells = nontrep_cells) <- "nonTREP"

    markers <- tryCatch({
        FindMarkers(
            seurat_obj,
            ident.1         = "TREP",
            ident.2         = "nonTREP",
            test.use        = "wilcox",
            min.pct         = deg_params$min_pct,
            logfc.threshold = deg_params$logfc_threshold,
            min.cells.group = deg_params$min_cells
        )
    }, error = function(e) {
        cat("    Error in FindMarkers:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(markers) || nrow(markers) == 0) {
        cat("    No DEGs found\n")
        return(NULL)
    }

    markers$gene <- rownames(markers)

    up_in_trep <- markers %>%
        filter(p_val < deg_params$p_cutoff &
               avg_log2FC > deg_params$logfc_cutoff) %>%
        arrange(p_val)

    down_in_trep <- markers %>%
        filter(p_val < deg_params$p_cutoff &
               avg_log2FC < -deg_params$logfc_cutoff) %>%
        arrange(p_val)

    up_in_trep   <- head(up_in_trep,   deg_params$max_deg)
    down_in_trep <- head(down_in_trep, deg_params$max_deg)

    cat("    Up in TREP:", nrow(up_in_trep),
        "| Down in TREP (up in non-TREP):", nrow(down_in_trep), "\n")

    return(list(
        all_markers  = markers,
        up_in_trep   = up_in_trep,
        down_in_trep = down_in_trep
    ))
}

# ------------------------------------------------------------------
# SECTION 7: Helper - safe sheet name (max 31 chars, no collision)
# ------------------------------------------------------------------
make_sheet_name <- function(name, existing_names) {
    truncated <- substr(name, 1L, 31L)
    if (!(truncated %in% existing_names)) return(truncated)
    for (i in 2L:99L) {
        suffix    <- paste0("_", i)
        candidate <- paste0(substr(name, 1L, 31L - nchar(suffix)), suffix)
        if (!(candidate %in% existing_names)) return(candidate)
    }
    stop("Cannot generate unique sheet name for: ", name)
}

# ------------------------------------------------------------------
# SECTION 8: Main loop
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("DEG ANALYSIS - Pre/Post-IP TREP vs non-TREP\n")
cat("Allen MTG CARM1P1 | DFG gene set\n")
cat("Donors:", paste(donors_of_interest, collapse = ", "), "\n")
cat("=============================================\n\n")

all_deg_results <- list()

for (donor in donors_of_interest) {

    key <- paste0(donor, "|", gene_set_label)
    cat("\n========== Processing:", donor, "(key:", key, ") ==========\n")

    # Retrieve GAM fit for this donor
    result <- valid_results[[key]]
    if (is.null(result)) {
        cat("  Skipping:", donor, "- no valid DFG fit in valid_results\n")
        next
    }

    # Retrieve IP value
    ip_value <- inflection_points[key]
    if (is.na(ip_value) || !is.numeric(ip_value)) {
        cat("  Skipping:", donor, "- no valid IP value in inflection_points\n")
        next
    }
    cat("  IP:", sprintf("%.4f", ip_value), "\n")

    # Resolve per-donor parameters (use overrides if defined, else defaults)
    deg_params <- default_deg_params
    if (!is.null(donor_deg_overrides[[donor]])) {
        deg_params <- modifyList(deg_params, donor_deg_overrides[[donor]])
        cat("  Using custom DEG thresholds for", donor, "\n")
    }

    # Classify cells into TREP / non-TREP x Pre-IP / Post-IP
    cell_class <- classify_cells(donor, result, ip_value)

    cat("  Cell counts per subpopulation:\n")
    print(table(cell_class$subpopulation))

    for (region in c("Pre_IP", "Post_IP")) {

        region_label <- gsub("_", "-", region)  # "Pre-IP" or "Post-IP"
        cat("\n  ---", region_label, "region ---\n")

        trep_cells    <- cell_class$barcode[
            cell_class$subpopulation == paste0(region, "_TREP")]
        nontrep_cells <- cell_class$barcode[
            cell_class$subpopulation == paste0(region, "_nonTREP")]

        deg <- run_deg_analysis(seurat_obj, trep_cells, nontrep_cells,
                                deg_params, region_label, donor)

        if (is.null(deg)) next

        key_deg                    <- paste0(donor, "_", region)
        all_deg_results[[key_deg]] <- deg
    }

    cat("  [", donor, "complete]\n")
}

# ------------------------------------------------------------------
# SECTION 9: Export DEG lists - per-donor-region sheets
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("EXPORTING DEG RESULTS\n")
cat("=============================================\n\n")

# ── Upregulated ──────────────────────────────────────────────────
deg_list_up <- list()

for (k in names(all_deg_results)) {
    deg <- all_deg_results[[k]]

    if (nrow(deg$up_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_upTREP"), names(deg_list_up))
        deg_list_up[[sn]] <- deg$up_in_trep
    }

    if (nrow(deg$down_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_upNonTREP"), names(deg_list_up))
        deg_list_up[[sn]] <- deg$down_in_trep
    }
}

if (length(deg_list_up) > 0) {
    write_xlsx(deg_list_up,
               path = "DEG_CARM1P1_PrePost_IP_TREP_vs_nonTREP.xlsx")
    cat("Exported", length(deg_list_up),
        "sheets (upregulated DEGs) → DEG_CARM1P1_PrePost_IP_TREP_vs_nonTREP.xlsx\n")
} else {
    cat("No upregulated DEG results found.\n")
}

# ── Downregulated ────────────────────────────────────────────────
deg_list_down <- list()

for (k in names(all_deg_results)) {
    deg <- all_deg_results[[k]]

    if (nrow(deg$down_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_downTREP"), names(deg_list_down))
        deg_list_down[[sn]] <- deg$down_in_trep
    }

    if (nrow(deg$up_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_downNonTREP"), names(deg_list_down))
        deg_list_down[[sn]] <- deg$up_in_trep
    }
}

if (length(deg_list_down) > 0) {
    write_xlsx(deg_list_down,
               path = "DEG_CARM1P1_Down_PrePost_IP_TREP_vs_nonTREP.xlsx")
    cat("Exported", length(deg_list_down),
        "sheets (downregulated DEGs) → DEG_CARM1P1_Down_PrePost_IP_TREP_vs_nonTREP.xlsx\n")
} else {
    cat("No downregulated DEG results found.\n")
}

# ------------------------------------------------------------------
# SECTION 10: Final summary
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================\n\n")

cat("Dataset:    Allen MTG | CARM1P1 vs DFG gene set\n")
cat("Donors:     H200.1023, H200.1025, H200.1030\n")
cat("IP source:  inflection_points[\"donor|DFG\"] (from Part 4/5)\n\n")

cat("Files generated:\n")
cat("  • DEG_CARM1P1_PrePost_IP_TREP_vs_nonTREP.xlsx        (upregulated DEGs)\n")
cat("  • DEG_CARM1P1_Down_PrePost_IP_TREP_vs_nonTREP.xlsx   (downregulated DEGs)\n\n")

cat("DEG criteria (default):\n")
cat("  p < 0.01, |log2FC| > 0.2, min.pct = 0.10,\n")
cat("  min 3 cells per group, max 500 DEGs\n")
cat("  Reliability guard: both groups must have >= 3 cells\n\n")

cat("Sheet naming: donor_region_upTREP / upNonTREP / downTREP / downNonTREP\n\n")

cat(sprintf("Parallelization: %d workers | RAM limit: %.1f GB\n\n",
            n_cores, ram_limit_bytes / 1024^3))
            

##############################################################################
# Allen MTG CARM1P1 - Part 7: Monocle3 Trajectory
#                              (TREP and non-TREP Cells, Quantitative Analysis)
#
# Requires:
#   - AllenMTG_CARM1P1_GAM_models.rds  (from Part 1)
#   - inflection_points                 (named numeric vector from Part 4/5,
#                                        keyed by "donor|DFG")
#   - seurat_obj                        (Seurat object with scRNA-seq data)
##############################################################################

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggridges)
library(openxlsx)
library(cluster)
library(FNN)
library(broom)
library(mgcv)

# ------------------------------------------------------------------
# SECTION 1: Load Part 1 models
# ------------------------------------------------------------------
cat("Loading Part 1 models...\n")
rds_data    <- readRDS("AllenMTG_CARM1P1_GAM_models.rds")
all_results <- rds_data$all_results

donors_of_interest <- c("H200.1023", "H200.1025", "H200.1030")
gene_set_label     <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)
valid_results <- valid_results[
    grepl(paste0("\\|", gene_set_label, "$"), names(valid_results)) &
        sapply(strsplit(names(valid_results), "\\|"), `[`, 1) %in% donors_of_interest
]

cat("DFG fits retained for analysis:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n")

# ------------------------------------------------------------------
# SECTION 2: Check prerequisites
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 4/5 (IP detection) which creates a named numeric\n",
        "vector 'inflection_points' keyed by \"donor|DFG\".\n",
        "Example: inflection_points[\"H200.1030|DFG\"]"
    ))
}

if (!exists("seurat_obj")) {
    stop(paste0(
        "seurat_obj not found in environment.\n",
        "Load your Seurat object before running Part 7.\n",
        "Cell barcodes in seurat_obj must match rownames(result$data)."
    ))
}

# ------------------------------------------------------------------
# SECTION 3: Monocle3 trajectory per donor
# ------------------------------------------------------------------
create_monocle3_trajectory <- function(donor) {
    key <- paste0(donor, "|", gene_set_label)
    cat("=== CREATING MONOCLE3 TRAJECTORY FOR", donor, "(key:", key, ") ===\n")
    
    result <- valid_results[[key]]
    if (is.null(result)) {
        cat("  Skipping:", donor, "- no valid DFG fit\n")
        return(NULL)
    }
    
    # ── GAM data ──────────────────────────────────────────────────
    original_gam_data <- result$data   # data.frame: CARM1P1, Expression;
    # rownames = cell barcodes
    if (is.null(original_gam_data)) {
        stop("GAM data not found for donor ", donor)
    }
    cat("GAM data dimensions:", nrow(original_gam_data), "x", ncol(original_gam_data), "\n")
    
    # ── IP value ──────────────────────────────────────────────────
    ip_value <- inflection_points[key]
    if (is.na(ip_value)) {
        stop("No IP value found in inflection_points for key ", key,
             ". Ensure Part 4/5 has been run.")
    }
    cat("IP value for", donor, ":", ip_value, "\n")
    
    # ── Cell classification (identical logic to Part 5 / Part 6) ──
    gam_model     <- result$best_model
    dev_explained <- summary(gam_model)$dev.expl
    
    null_model   <- gam(Expression ~ 1, data = original_gam_data)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(gam_model)
    
    null_sq_diff  <- (original_gam_data$Expression - null_fitted)^2
    model_sq_diff <- (original_gam_data$Expression - model_fitted)^2
    
    explanatory_power <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )
    
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    target_cells   <- round(nrow(original_gam_data) * dev_explained)
    target_cells   <- max(1L, min(target_cells, nrow(original_gam_data)))
    
    cat("Target TREP cells:", target_cells, "out of", nrow(original_gam_data), "total cells\n")
    
    trep_cells <- rownames(original_gam_data)[sorted_indices[seq_len(target_cells)]]
    
    original_gam_data$is_trep <- ifelse(
        rownames(original_gam_data) %in% trep_cells, "TREP", "non-TREP")
    original_gam_data$timing  <- ifelse(
        original_gam_data$CARM1P1 < ip_value, "Pre-IP", "Post-IP")
    original_gam_data$cell_group <- paste0(
        original_gam_data$timing, "_", original_gam_data$is_trep)
    
    # ── Match cells with Seurat object ────────────────────────────
    common_cells <- intersect(colnames(seurat_obj), rownames(original_gam_data))
    cat("Common cells found:", length(common_cells), "\n")
    
    if (length(common_cells) == 0) {
        stop("No common cells found between seurat_obj and GAM data for ", donor)
    }
    
    sample_seurat       <- subset(seurat_obj, cells = common_cells)
    sample_data_matched <- original_gam_data[common_cells, ]
    
    cat("Final processing:", nrow(sample_data_matched), "cells\n")
    cat("Cell group distribution:\n")
    print(table(sample_data_matched$cell_group))
    
    # ── Build Monocle3 CDS ────────────────────────────────────────
    cat("Converting to Monocle3 format...\n")
    
    count_matrix <- GetAssayData(sample_seurat, slot = "counts", assay = "RNA")
    
    cell_metadata <- data.frame(
        cell_id           = colnames(sample_seurat),
        CARM1P1           = sample_data_matched$CARM1P1,
        Expression        = sample_data_matched$Expression,   # DFG score
        cell_group        = sample_data_matched$cell_group,
        is_trep           = sample_data_matched$is_trep,
        timing            = sample_data_matched$timing,
        explanatory_power = explanatory_power[
            match(colnames(sample_seurat), rownames(original_gam_data))],
        stringsAsFactors  = FALSE
    )
    rownames(cell_metadata) <- cell_metadata$cell_id
    
    gene_metadata <- data.frame(
        gene_id         = rownames(count_matrix),
        gene_short_name = rownames(count_matrix),
        stringsAsFactors = FALSE
    )
    rownames(gene_metadata) <- gene_metadata$gene_id
    
    cds <- new_cell_data_set(
        expression_data = count_matrix,
        cell_metadata   = cell_metadata,
        gene_metadata   = gene_metadata
    )
    
    # ── Preprocessing & trajectory ────────────────────────────────
    cat("Preprocessing data...\n")
    cds <- preprocess_cds(cds, num_dim = 200, verbose = FALSE)
    
    cat("Reducing dimensions...\n")
    cds <- reduce_dimension(cds, verbose = FALSE)
    
    cat("Clustering cells...\n")
    cds <- cluster_cells(cds, verbose = FALSE, resolution = 0.0005, k = 250)
    
    cat("Learning trajectory...\n")
    cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE, close_loop = FALSE)
    
    # Root = cell with highest CARM1P1 expression
    root_cell <- colnames(cds)[which.max(colData(cds)$CARM1P1)]
    cds <- order_cells(cds, root_cells = root_cell, verbose = FALSE)
    
    # ── Colors ────────────────────────────────────────────────────
    colors <- c(
        "Pre-IP_non-TREP"  = "#CCCCCC",
        "Post-IP_non-TREP" = "#666666",
        "Pre-IP_TREP"      = "#DDA0DD",
        "Post-IP_TREP"     = "#6666FF"
    )
    
    # ── Trajectory plot ───────────────────────────────────────────
    cat("Creating trajectory plot...\n")
    
    trajectory_plot <- plot_cells(cds,
                                  color_cells_by                = "cell_group",
                                  label_cell_groups             = FALSE,
                                  label_leaves                  = FALSE,
                                  label_branch_points           = FALSE,
                                  label_roots                   = FALSE,
                                  show_trajectory_graph         = TRUE,
                                  graph_label_size              = 3,
                                  cell_size                     = 1.5,
                                  cell_stroke                   = 0.5,
                                  trajectory_graph_color        = "#333333",
                                  trajectory_graph_segment_size = 1) +
        scale_color_manual(values = colors) +
        labs(
            title    = paste("Monocle3 Trajectory -", donor),
            subtitle = paste("TREP cells (purple/blue) vs non-TREP cells (gray),",
                             "Pre/Post IP =", round(ip_value, 4)),
            color    = "Cell Group"
        ) +
        theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title  = element_text(size = 12, face = "bold"),
            legend.text   = element_text(size = 10), 
            axis.text.x   = element_text(size = 6),
            axis.text.y   = element_text(size = 6)  
        )
    
    # ── Pseudotime plot ───────────────────────────────────────────
    cat("Creating pseudotime plot...\n")
    
    pseudotime_plot <- plot_cells(cds,
                                  color_cells_by                = "pseudotime",
                                  label_cell_groups             = FALSE,
                                  label_leaves                  = FALSE,
                                  label_branch_points           = FALSE,
                                  label_roots                   = FALSE,
                                  show_trajectory_graph         = TRUE,
                                  graph_label_size              = 3,
                                  cell_size                     = 0.5,
                                  cell_stroke                   = 0.5,
                                  trajectory_graph_color        = "#333333",
                                  trajectory_graph_segment_size = 1) +
        scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
        labs(
            title    = paste("Monocle3 Pseudotime -", donor),
            subtitle = paste("Lighter colors = higher pseudotime, Pre/Post IP =",
                             round(ip_value, 4))
        ) +
        theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title  = element_text(size = 12, face = "bold"),
            legend.text   = element_text(size = 10)
        )
    
    print(trajectory_plot)
    print(pseudotime_plot)
    
    # ── Summary ───────────────────────────────────────────────────
    cat("\n=== TRAJECTORY SUMMARY ===\n")
    cat("Donor:", donor, "\n")
    cat("Total cells:", ncol(cds), "\n")
    cat("Inflection point:", ip_value, "\n")
    cat("Deviance explained:", round(dev_explained * 100, 1), "%\n")
    cat("Root cell:", root_cell, "\n")
    
    cat("\nCell group counts:\n")
    group_counts <- table(cell_metadata$cell_group)
    for (g in names(group_counts)) {
        cat(" ", g, ":", group_counts[g],
            "(", round(group_counts[g] / sum(group_counts) * 100, 1), "%)\n")
    }
    
    return(cds)
}

# ------------------------------------------------------------------
# SECTION 4: Run trajectory for all donors
# ------------------------------------------------------------------
cat("\n\nGenerating Monocle3 trajectory plots for all donors...\n")
cat(strrep("=", 60), "\n")

all_trajectories <- list()

for (donor in donors_of_interest) {
    tryCatch({
        cds_result <- create_monocle3_trajectory(donor)
        if (!is.null(cds_result)) {
            all_trajectories[[donor]] <- cds_result
        }
    }, error = function(e) {
        cat("Error creating trajectory for", donor, ":", conditionMessage(e), "\n")
    })
}

cat("\nCompleted generating Monocle3 trajectory plots for all donors.\n")


#######################################################################################
# Quantitative Analysis of Pre-IP and Post-IP TREP Cells in Monocle3 Trajectories
#######################################################################################

analyze_cell_type_clustering <- function(cds, donor_id, colors) {
    cat("\n=== ANALYZING CELL TYPE CLUSTERING FOR", donor_id, "===\n")
    
    umap_coords   <- reducedDims(cds)$UMAP
    cell_metadata <- colData(cds)
    
    analysis_data <- data.frame(
        cell_id    = rownames(umap_coords),
        UMAP1      = umap_coords[, 1],
        UMAP2      = umap_coords[, 2],
        cell_group = cell_metadata$cell_group,
        is_trep    = cell_metadata$is_trep,
        timing     = cell_metadata$timing,
        CARM1P1    = cell_metadata$CARM1P1,
        Expression = cell_metadata$Expression,
        stringsAsFactors = FALSE
    )
    
    dist_matrix <- dist(umap_coords)
    
    # ── 1. Within / between group distances ───────────────────────
    cat("Calculating within-group and between-group distances...\n")
    
    within_distances  <- list()
    between_distances <- list()
    cell_groups       <- unique(analysis_data$cell_group)
    
    for (group in cell_groups) {
        group_cells <- which(analysis_data$cell_group == group)
        if (length(group_cells) > 1) {
            gd <- as.matrix(dist_matrix)[group_cells, group_cells]
            within_distances[[group]] <- gd[upper.tri(gd)]
        }
    }
    
    for (i in seq_len(length(cell_groups) - 1)) {
        for (j in (i + 1):length(cell_groups)) {
            g1 <- cell_groups[i]; g2 <- cell_groups[j]
            idx1 <- which(analysis_data$cell_group == g1)
            idx2 <- which(analysis_data$cell_group == g2)
            if (length(idx1) > 0 && length(idx2) > 0) {
                bd <- as.matrix(dist_matrix)[idx1, idx2]
                between_distances[[paste(g1, "vs", g2, sep = "_")]] <- as.vector(bd)
            }
        }
    }
    
    # ── 2. Statistical testing (UMAP1) ────────────────────────────
    cat("Performing statistical tests on UMAP1 distributions...\n")
    
    statistical_results <- data.frame()
    
    post_ip_trep_umap1 <- analysis_data$UMAP1[analysis_data$cell_group == "Post-IP_TREP"]
    pre_ip_trep_umap1  <- analysis_data$UMAP1[analysis_data$cell_group == "Pre-IP_TREP"]
    
    if (length(post_ip_trep_umap1) > 0 && length(pre_ip_trep_umap1) > 0) {
        
        post_normal <- ifelse(
            length(post_ip_trep_umap1) >= 3 && length(post_ip_trep_umap1) <= 5000,
            shapiro.test(post_ip_trep_umap1)$p.value > 0.05, FALSE)
        pre_normal <- ifelse(
            length(pre_ip_trep_umap1) >= 3 && length(pre_ip_trep_umap1) <= 5000,
            shapiro.test(pre_ip_trep_umap1)$p.value > 0.05, FALSE)
        
        both_normal <- post_normal && pre_normal
        test_type   <- ifelse(both_normal, "t-test", "Mann-Whitney")
        
        if (both_normal) {
            tres        <- t.test(post_ip_trep_umap1, pre_ip_trep_umap1)
            p_value     <- tres$p.value
            effect_size <- abs(mean(post_ip_trep_umap1) - mean(pre_ip_trep_umap1)) /
                sqrt(((length(post_ip_trep_umap1) - 1) * var(post_ip_trep_umap1) +
                          (length(pre_ip_trep_umap1)  - 1) * var(pre_ip_trep_umap1)) /
                         (length(post_ip_trep_umap1) + length(pre_ip_trep_umap1) - 2))
        } else {
            wres        <- wilcox.test(post_ip_trep_umap1, pre_ip_trep_umap1)
            p_value     <- wres$p.value
            n1          <- length(post_ip_trep_umap1)
            n2          <- length(pre_ip_trep_umap1)
            effect_size <- (2 * wres$statistic) / (n1 * n2) - 1
        }
        
        result_row <- data.frame(
            Sample              = donor_id,
            Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
            Test_Type           = test_type,
            Post_IP_TREP_Median = median(post_ip_trep_umap1),
            Post_IP_TREP_IQR    = IQR(post_ip_trep_umap1),
            Post_IP_TREP_Mean   = mean(post_ip_trep_umap1),
            Post_IP_TREP_N      = length(post_ip_trep_umap1),
            Post_IP_TREP_Normal = post_normal,
            Pre_IP_TREP_Median  = median(pre_ip_trep_umap1),
            Pre_IP_TREP_IQR     = IQR(pre_ip_trep_umap1),
            Pre_IP_TREP_Mean    = mean(pre_ip_trep_umap1),
            Pre_IP_TREP_N       = length(pre_ip_trep_umap1),
            Pre_IP_TREP_Normal  = pre_normal,
            P_Value             = p_value,
            Effect_Size         = effect_size,
            stringsAsFactors    = FALSE
        )
        
        statistical_results <- rbind(statistical_results, result_row)
        
        cat("Post-IP_TREP: n =", length(post_ip_trep_umap1), ", Normal =", post_normal, "\n")
        cat("Pre-IP_TREP:  n =", length(pre_ip_trep_umap1),  ", Normal =", pre_normal,  "\n")
        cat("Test used:", test_type, ", p-value =", p_value, "\n")
    }
    
    # ── 3. Silhouette analysis ─────────────────────────────────────
    cat("Calculating silhouette scores for UMAP separation...\n")
    
    binary_labels <- ifelse(analysis_data$cell_group == "Post-IP_TREP", 1,
                            ifelse(analysis_data$cell_group == "Pre-IP_TREP",  2, NA))
    valid_idx    <- !is.na(binary_labels)
    valid_labels <- binary_labels[valid_idx]
    valid_coords <- umap_coords[valid_idx, ]
    
    avg_silhouette <- NA
    if (length(unique(valid_labels)) == 2 && nrow(valid_coords) > 2) {
        sil_scores     <- silhouette(valid_labels, dist(valid_coords))
        avg_silhouette <- mean(sil_scores[, 3])
        
        silhouette_row <- data.frame(
            Sample              = donor_id,
            Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
            Test_Type           = "Silhouette_Score_UMAP_Separation",
            Post_IP_TREP_Median = NA, Post_IP_TREP_IQR = NA,
            Post_IP_TREP_Mean   = NA, Post_IP_TREP_N   = NA,
            Post_IP_TREP_Normal = NA,
            Pre_IP_TREP_Median  = NA, Pre_IP_TREP_IQR  = NA,
            Pre_IP_TREP_Mean    = NA, Pre_IP_TREP_N    = NA,
            Pre_IP_TREP_Normal  = NA,
            P_Value             = NA,
            Effect_Size         = avg_silhouette,
            stringsAsFactors    = FALSE
        )
        statistical_results <- rbind(statistical_results, silhouette_row)
    }
    
    # ── 4. BH correction ──────────────────────────────────────────
    p_vals <- statistical_results$P_Value[!is.na(statistical_results$P_Value)]
    if (length(p_vals) > 0) {
        adj_p <- p.adjust(p_vals, method = "BH")
        statistical_results$BH_Adjusted_P_Value <- NA
        statistical_results$BH_Adjusted_P_Value[!is.na(statistical_results$P_Value)] <- adj_p
    }
    
    # ── 5. Ridgeline plot ─────────────────────────────────────────
    cat("Preparing ridgeline plot...\n")
    
    ridgeline_data <- data.frame(
        UMAP1      = analysis_data$UMAP1,
        Cell_Group = factor(analysis_data$cell_group,
                            levels = c("Pre-IP_TREP", "Pre-IP_non-TREP",
                                       "Post-IP_TREP", "Post-IP_non-TREP")),
        stringsAsFactors = FALSE
    )
    
    ridgeline_plot <- ggplot(ridgeline_data,
                             aes(x = UMAP1, y = Cell_Group, fill = Cell_Group)) +
        geom_density_ridges(alpha = 0.5, scale = 1.5, rel_min_height = 0.01) +
        scale_fill_manual(values = colors) +
        scale_y_discrete(
            limits = rev(c("Pre-IP_TREP", "Pre-IP_non-TREP",
                           "Post-IP_TREP", "Post-IP_non-TREP")),
            expand = expansion(mult = c(0, 0.2))) +
        coord_cartesian(clip = "off") +
        labs(
            title    = paste("Cell Type Clustering Analysis -", donor_id),
            subtitle = "UMAP1 distributions for each cell type",
            x        = "UMAP1",
            y        = "Cell Type",
            caption  = paste("Overall Silhouette Score:", round(avg_silhouette, 3),
                             "| Lower spread = tighter clustering")
        ) +
        theme_ridges() +
        theme(
            plot.title         = element_text(size = 14, face = "bold"),
            plot.subtitle      = element_text(size = 12),
            legend.title       = element_text(size = 12, face = "bold"),
            legend.text        = element_text(size = 10),
            axis.text.y        = element_text(size = 8),
            axis.text.x        = element_text(size = 8, margin = margin(t = 2)),
            plot.caption       = element_text(size = 9, hjust = 0),
            panel.grid.major.x = element_line(linetype = "dotted", color = "gray70"),
            panel.grid.minor.x = element_blank(),
            plot.margin        = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")
        ) +
        guides(fill = "none")
    
    print(ridgeline_plot)
    
    # ── 6. Summary stats table ─────────────────────────────────────
    summary_stats <- statistical_results %>%
        dplyr::filter(!grepl("Silhouette", Test_Type)) %>%
        dplyr::select(Cell_Group, Post_IP_TREP_Mean, Pre_IP_TREP_Mean,
                      P_Value, BH_Adjusted_P_Value, Effect_Size, Test_Type) %>%
        dplyr::mutate(
            Significant_BH     = ifelse(BH_Adjusted_P_Value < 0.05, "Yes", "No"),
            Clustering_Quality = dplyr::case_when(
                abs(Effect_Size) > 0.5 ~ "Strong",
                abs(Effect_Size) > 0.3 ~ "Moderate",
                abs(Effect_Size) > 0.1 ~ "Weak",
                TRUE                   ~ "Very Weak"
            )
        )
    
    cat("\nClustering Summary:\n")
    print(summary_stats)
    
    return(list(
        statistical_results = statistical_results,
        ridgeline_data      = ridgeline_data,
        ridgeline_plot      = ridgeline_plot,
        summary_stats       = summary_stats,
        overall_silhouette  = avg_silhouette
    ))
}

# ------------------------------------------------------------------
# SECTION 5: Process all donors & export
# ------------------------------------------------------------------
process_all_donors_clustering <- function(all_trajectories,
                                          output_file = "monocle3_clustering_analysis_CARM1P1_DFG.xlsx") {
    
    cat("\n=== PROCESSING ALL DONORS FOR CLUSTERING ANALYSIS ===\n")
    
    colors <- c(
        "Pre-IP_non-TREP"  = "#CCCCCC",
        "Post-IP_non-TREP" = "#666666",
        "Pre-IP_TREP"      = "#DDA0DD",
        "Post-IP_TREP"     = "#6666FF"
    )
    
    all_statistical_results <- data.frame()
    
    for (donor_name in names(all_trajectories)) {
        cat("\nProcessing", donor_name, "...\n")
        
        tryCatch({
            res <- analyze_cell_type_clustering(
                all_trajectories[[donor_name]], donor_name, colors)
            all_statistical_results <- rbind(all_statistical_results,
                                             res$statistical_results)
        }, error = function(e) {
            cat("Error analyzing clustering for", donor_name, ":", conditionMessage(e), "\n")
        })
    }
    
    # BH correction across all donors
    cat("Applying Benjamini-Hochberg correction across all donors...\n")
    
    if (nrow(all_statistical_results) > 0) {
        p_rows <- !is.na(all_statistical_results$P_Value)
        if (sum(p_rows) > 0) {
            adj_p <- p.adjust(all_statistical_results$P_Value[p_rows], method = "BH")
            all_statistical_results$BH_Adjusted_P_Value <- NA
            all_statistical_results$BH_Adjusted_P_Value[p_rows] <- adj_p
            cat("BH correction applied to", sum(p_rows), "p-values.\n")
        }
    }
    
    # Export to Excel
    cat("\nExporting results to:", output_file, "\n")
    
    wb <- createWorkbook()
    
    addWorksheet(wb, "Statistical_Results")
    writeData(wb, "Statistical_Results", all_statistical_results)
    
    metadata <- data.frame(
        Description = c(
            "Dataset: Allen MTG CARM1P1 vs DFG gene set",
            "Donors: H200.1023, H200.1025, H200.1030",
            "Post_IP_TREP_Median: Median UMAP1 value for Post-IP_TREP cells",
            "Post_IP_TREP_IQR: Interquartile range of UMAP1 for Post-IP_TREP cells",
            "Post_IP_TREP_Mean: Mean UMAP1 value for Post-IP_TREP cells",
            "Post_IP_TREP_N: Number of Post-IP_TREP cells",
            "Post_IP_TREP_Normal: Shapiro-Wilk normality (p > 0.05) for Post-IP_TREP",
            "Pre_IP_TREP_Median: Median UMAP1 value for Pre-IP_TREP cells",
            "Pre_IP_TREP_IQR: Interquartile range of UMAP1 for Pre-IP_TREP cells",
            "Pre_IP_TREP_Mean: Mean UMAP1 value for Pre-IP_TREP cells",
            "Pre_IP_TREP_N: Number of Pre-IP_TREP cells",
            "Pre_IP_TREP_Normal: Shapiro-Wilk normality (p > 0.05) for Pre-IP_TREP",
            "P_Value: Raw p-value from t-test (normal) or Mann-Whitney (non-normal)",
            "BH_Adjusted_P_Value: Benjamini-Hochberg corrected p-value (across all 3 donors)",
            "Effect_Size: Cohen's d (t-test) or Cliff's delta (Mann-Whitney)",
            "Test_Type: Statistical test used",
            "IP source: inflection_points[\"donor|DFG\"] from Part 4/5"
        )
    )
    
    addWorksheet(wb, "Metadata")
    writeData(wb, "Metadata", metadata)
    
    saveWorkbook(wb, output_file, overwrite = TRUE)
    cat("Analysis complete! Results saved to:", output_file, "\n")
    
    # Overall summary
    cat("\n=== OVERALL CLUSTERING SUMMARY ===\n")
    
    if (nrow(all_statistical_results) > 0) {
        comp <- all_statistical_results %>%
            dplyr::filter(Cell_Group == "Post-IP_TREP_vs_Pre-IP_TREP" &
                              !grepl("Silhouette", Test_Type))
        
        if (nrow(comp) > 0) {
            overall_summary <- comp %>%
                dplyr::summarise(
                    Total_Donors            = dplyr::n(),
                    Mean_Effect_Size        = mean(Effect_Size, na.rm = TRUE),
                    Mean_Post_IP_TREP_UMAP1 = mean(Post_IP_TREP_Mean, na.rm = TRUE),
                    Mean_Pre_IP_TREP_UMAP1  = mean(Pre_IP_TREP_Mean, na.rm = TRUE),
                    Significant_Raw         = sum(P_Value < 0.05, na.rm = TRUE),
                    Significant_BH          = sum(BH_Adjusted_P_Value < 0.05, na.rm = TRUE),
                    T_Test_Used             = sum(Test_Type == "t-test", na.rm = TRUE),
                    Mann_Whitney_Used       = sum(Test_Type == "Mann-Whitney", na.rm = TRUE),
                    .groups = 'drop'
                )
            print(overall_summary)
            
            cat("\nRaw vs BH-adjusted significance per donor:\n")
            print(comp %>%
                      dplyr::select(Sample, P_Value, BH_Adjusted_P_Value) %>%
                      dplyr::mutate(
                          Raw_Significant = P_Value < 0.05,
                          BH_Significant  = BH_Adjusted_P_Value < 0.05
                      ))
        } else {
            cat("No UMAP1 comparison results found.\n")
        }
    }
    
    return(list(statistical_results = all_statistical_results))
}

# ------------------------------------------------------------------
# SECTION 6: Run clustering analysis
# ------------------------------------------------------------------
cat("\n\nStarting clustering analysis for all donors...\n")
clustering_analysis_results <- process_all_donors_clustering(all_trajectories)

cat("\nPart 7 complete.\n")
cat("Output: monocle3_clustering_analysis_CARM1P1_DFG.xlsx\n")
