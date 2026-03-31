##########################################################################
# Validation of CEP-IP Framework: Neftel Glioblastoma Multiforme (GBM) Dataset
# Dataset: Neftel et al. 2019 (GSE131928) doi: https://doi.org/10.1016/j.cell.2019.06.024
# Part 1: Metadata Extraction, Sample Selection
#                        & Dominant Cell State Assignment
##########################################################################

library(data.table)
library(dplyr)
library(openxlsx)

# ---- File paths ----
base_dir    <- "C:/Users/...#directory"
meta_path   <- file.path(base_dir, "IDHwt.GBM.Metadata.SS2.txt")
expr_path   <- file.path(base_dir, "IDHwtGBM.processed.SS2.logTPM.txt")
output_path <- file.path(base_dir, "Neftel_SS2_Part1_Metadata.xlsx")

# ============================================================
# SECTION 1: Hardcoded adult/pediatric classification
# (from GSE131928 supplementary - avoids external file dependency)
# ============================================================
adult_samples <- c(
  "MGH101", "MGH100", "MGH102", "MGH104", "MGH105", "MGH106",
  "MGH110", "MGH113", "MGH115", "MGH121", "MGH122", "MGH124",
  "MGH125", "MGH128", "MGH129", "MGH136", "MGH143", "MGH151",
  "MGH152", "MGH66",  "MGH114", "MGH118", "MGH126"
)

pediatric_samples <- c(
  "BT749", "BT771", "BT830", "MGH85", "BT85",  # BT85 = MGH85 in SS2 naming
  "BT1160", "BT1187", "BT786", "BT920"
)

# ============================================================
# SECTION 2: Load and parse SS2 metadata
# The file has a TYPE descriptor row (row 2 after header) - skip it
# ============================================================
cat("Loading SS2 metadata...\n")

meta_raw <- fread(meta_path, header = TRUE, sep = "\t", data.table = FALSE)

# The TYPE row contains literal strings "group" or "numeric" in each column
# Identify and remove it by checking if the first column value is "TYPE"
meta <- meta_raw[meta_raw[[1]] != "TYPE", ]

# Rename first column to Cell_ID
colnames(meta)[1] <- "Cell_ID"

cat("Metadata rows after TYPE-row removal:", nrow(meta), "\n")
cat("Columns:", paste(colnames(meta), collapse = ", "), "\n")

# ---- Coerce numeric score columns ----
score_cols <- c("MESlike2", "MESlike1", "AClike", "OPClike",
                "NPClike1", "NPClike2", "G1S", "G2M")
meta[score_cols]    <- lapply(meta[score_cols], as.numeric)
meta$GenesExpressed <- as.numeric(meta$GenesExpressed)

# ============================================================
# SECTION 3: Annotate age group
# ============================================================
meta <- meta %>%
  mutate(
    Age_Group = case_when(
      Sample %in% adult_samples     ~ "Adult",
      Sample %in% pediatric_samples ~ "Pediatric",
      TRUE                          ~ "Unknown"
    )
  )

# Warn if any samples are unclassified
unknown_samples <- unique(meta$Sample[meta$Age_Group == "Unknown"])
if (length(unknown_samples) > 0) {
  cat("WARNING - Unclassified samples:", paste(unknown_samples, collapse = ", "), "\n")
} else {
  cat("All", n_distinct(meta$Sample), "samples successfully classified by age group.\n")
}

cat("\nCellAssignment categories:", paste(unique(meta$CellAssignment), collapse = ", "), "\n")

# ============================================================
# SECTION 4: Filter subsets
# ============================================================
meta_adult_all       <- meta %>% filter(Age_Group == "Adult")
meta_adult_malignant <- meta %>% filter(Age_Group == "Adult",
                                        CellAssignment == "Malignant")

cat("\nTotal SS2 cells (all)            :", nrow(meta), "\n")
cat("Adult cells (all assignments)    :", nrow(meta_adult_all), "\n")
cat("Adult malignant cells            :", nrow(meta_adult_malignant), "\n")

# ============================================================
# SECTION 5: Dominant Cell State Assignment
# ============================================================
# Four canonical GBM metacell states (Neftel et al.):
#   MES-like  = max(MESlike1, MESlike2)
#   AC-like   = AClike
#   OPC-like  = OPClike
#   NPC-like  = max(NPClike1, NPClike2)
#
# Cell-level dominant state = state with highest representative score.
# Co-dominant (mixed) = top-2 score gap < mixed_threshold (default 0.3).
# After running, inspect State_Confidence distribution:
#   If >40% "Low (mixed)" → increase threshold to 0.5
#   If very few "Low (mixed)" → decrease threshold to 0.2
# ============================================================

assign_dominant_state <- function(df, mixed_threshold = 0.3) {

  df <- df %>%
    mutate(
      MES_score = pmax(MESlike1, MESlike2, na.rm = TRUE),
      NPC_score = pmax(NPClike1, NPClike2, na.rm = TRUE),
      AC_score  = AClike,
      OPC_score = OPClike
    )

  score_mat <- as.matrix(df[, c("MES_score", "NPC_score", "AC_score", "OPC_score")])
  storage.mode(score_mat) <- "double"   # ensure pure numeric, no list coercion
  state_names <- c("MES", "NPC", "AC", "OPC")

  # Explicit anonymous functions - which.max() returns a named list when called
  # via apply() on matrix rows unless input is a clean numeric vector
  top_idx      <- apply(score_mat, 1, function(x) which.max(replace(x, is.na(x), -Inf)))
  top_state    <- state_names[as.integer(top_idx)]
  top_score    <- apply(score_mat, 1, function(x) max(x, na.rm = TRUE))

  second_idx   <- apply(score_mat, 1, function(x) {
    x[is.na(x)] <- -Inf
    order(x, decreasing = TRUE)[2]
  })
  second_state <- state_names[as.integer(second_idx)]
  second_score <- apply(score_mat, 1, function(x) {
    x[is.na(x)] <- -Inf
    sort(x, decreasing = TRUE)[2]
  })

  score_gap   <- top_score - second_score
  is_mixed    <- score_gap < mixed_threshold

  df <- df %>%
    mutate(
      Top_State        = top_state,
      Second_State     = second_state,
      Score_Gap        = round(score_gap, 4),
      Dominant_State   = ifelse(is_mixed,
                                paste0(top_state, "+", second_state),
                                top_state),
      State_Confidence = case_when(
        score_gap >= 1.0 ~ "High",
        score_gap >= 0.3 ~ "Medium",
        TRUE             ~ "Low (mixed)"
      )
    )
  return(df)
}

meta_adult_malignant <- assign_dominant_state(meta_adult_malignant,
                                              mixed_threshold = 0.3)

# ============================================================
# SECTION 6: Per-tumor summary
# ============================================================
tumor_summary <- meta_adult_malignant %>%
  group_by(Sample) %>%
  summarise(
    Platform              = "Smart-seq2",
    Age_Group             = first(Age_Group),
    GBMType               = first(GBMType),
    Total_Malignant_Cells = n(),
    Mean_GenesExpressed   = round(mean(GenesExpressed, na.rm = TRUE), 1),

    # Mean metacell scores
    Mean_MES = round(mean(MES_score, na.rm = TRUE), 4),
    Mean_NPC = round(mean(NPC_score, na.rm = TRUE), 4),
    Mean_AC  = round(mean(AC_score,  na.rm = TRUE), 4),
    Mean_OPC = round(mean(OPC_score, na.rm = TRUE), 4),

    # Cell state composition (%)
    Pct_MES   = round(100 * sum(Top_State == "MES") / n(), 1),
    Pct_NPC   = round(100 * sum(Top_State == "NPC") / n(), 1),
    Pct_AC    = round(100 * sum(Top_State == "AC")  / n(), 1),
    Pct_OPC   = round(100 * sum(Top_State == "OPC") / n(), 1),
    Pct_Mixed = round(100 * sum(grepl("\\+", Dominant_State)) / n(), 1),

    # Majority-vote dominant state (cell-level Top_State)
    Dominant_State_MajorityVote = names(sort(table(Top_State), decreasing = TRUE))[1],

    # Cell cycle
    Mean_G1S = round(mean(G1S, na.rm = TRUE), 4),
    Mean_G2M = round(mean(G2M, na.rm = TRUE), 4),

    Meets_100cell_threshold = ifelse(n() >= 100, "Yes", "No"),
    .groups = "drop"
  ) %>%
  # Tumor-level dominant state by highest mean score
  mutate(
    Dominant_State_MeanScore = case_when(
      pmax(Mean_MES, Mean_NPC, Mean_AC, Mean_OPC) == Mean_MES ~ "MES",
      pmax(Mean_MES, Mean_NPC, Mean_AC, Mean_OPC) == Mean_NPC ~ "NPC",
      pmax(Mean_MES, Mean_NPC, Mean_AC, Mean_OPC) == Mean_AC  ~ "AC",
      TRUE                                                     ~ "OPC"
    )
  ) %>%
  arrange(desc(Total_Malignant_Cells))

# ============================================================
# SECTION 7: Summary statistics
# ============================================================
cat("\n===== SS2 ADULT MALIGNANT DATASET SUMMARY =====\n")
cat("Total adult malignant cells  :", nrow(meta_adult_malignant), "\n")
cat("Total tumor samples          :", n_distinct(meta_adult_malignant$Sample), "\n")
cat("Samples >= 100 cells         :", sum(tumor_summary$Meets_100cell_threshold == "Yes"), "\n")
cat("\nDominant state distribution (cell-level Top_State):\n")
print(sort(table(meta_adult_malignant$Top_State), decreasing = TRUE))
cat("\nDominant state distribution (including mixed):\n")
print(sort(table(meta_adult_malignant$Dominant_State), decreasing = TRUE))
cat("\nState_Confidence distribution:\n")
print(table(meta_adult_malignant$State_Confidence))
cat("================================================\n")

cat("\nPer-tumor summary:\n")
print(
  tumor_summary %>%
    select(Sample, GBMType, Total_Malignant_Cells,
           Dominant_State_MeanScore, Dominant_State_MajorityVote,
           Pct_MES, Pct_NPC, Pct_AC, Pct_OPC, Pct_Mixed,
           Meets_100cell_threshold)
)

# ============================================================
# SECTION 8: Export to Excel (4 sheets)
# ============================================================
wb <- createWorkbook()

addWorksheet(wb, "Tumor_Summary")
writeData(wb, "Tumor_Summary", tumor_summary)

addWorksheet(wb, "Adult_Malignant_Cells")
writeData(wb, "Adult_Malignant_Cells", meta_adult_malignant)

addWorksheet(wb, "Adult_All_Cells")
writeData(wb, "Adult_All_Cells", meta_adult_all)

addWorksheet(wb, "Full_Metadata")
writeData(wb, "Full_Metadata", meta)

saveWorkbook(wb, output_path, overwrite = TRUE)
cat("\nPart 1 complete. Excel saved to:\n", output_path, "\n")


#############################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 2: Seurat Processing, UMAP & Gene Panel Analysis
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
# Adult malignant cells only (from Part 1 output)
#############################################################

# =========================================
# USER SETTINGS
# =========================================
GOI           <- "CX3CR1"         # gene of interest for UMAP feature plot
OUTPUT_PREFIX <- "NeftelGBM_SS2_AdultMalignant"

# Gene panel for box+jitter plots (Section 6)
GENE_PANEL <- c("CX3CR1", "AXL", "HAVCR2", "CSF1R",
                "TGFBI",  "CXCR4", "CLEC4A",
                "ALDOC",  "AQP4",
                "EGFR",   "STAT3",
                "MERTK",  "CD44")

# =========================================
# File paths
# =========================================
base_dir     <- "C:/Users/...#directory"
expr_path    <- file.path(base_dir, "IDHwtGBM.processed.SS2.logTPM.txt")
meta_path    <- file.path(base_dir, "Neftel_SS2_Part1_Metadata.xlsx")  # Part 1 output

# =========================================
# Libraries
# =========================================
library(data.table)
library(dplyr)
library(future)
library(ggplot2)
library(openxlsx)
library(parallel)
library(scales)
library(Seurat)
library(viridis)
library(rstatix)
library(writexl)
library(Matrix)

# =========================================
# 1. Load Part 1 metadata (adult malignant cells)
# =========================================
cat("Loading Part 1 adult malignant metadata...\n")

# Read the Adult_Malignant_Cells sheet from Part 1 Excel output
meta_malignant <- read.xlsx(meta_path, sheet = "Adult_Malignant_Cells")

cat("Adult malignant cells loaded:", nrow(meta_malignant), "\n")
cat("Tumors:", n_distinct(meta_malignant$Sample), "\n")
cat("Columns:", paste(colnames(meta_malignant), collapse = ", "), "\n")

# =========================================
# 2. Load SS2 logTPM expression matrix
# =========================================
cat("\nLoading SS2 logTPM expression matrix...\n")

expr_raw  <- fread(expr_path, header = TRUE, sep = "\t", data.table = FALSE)

# Row 2 in this file is also a TYPE descriptor row - remove it
expr_raw  <- expr_raw[expr_raw[[1]] != "TYPE", ]

# First column = gene names; remaining columns = cell IDs
gene_names <- expr_raw[[1]]
expr_mat   <- as.matrix(expr_raw[, -1])
rownames(expr_mat) <- make.unique(gene_names)

cat("Expression matrix dimensions (genes x cells):",
    nrow(expr_mat), "x", ncol(expr_mat), "\n")

# =========================================
# 3. Subset to adult malignant cells only
# =========================================
cat("\nSubsetting to adult malignant cells...\n")

# Cell IDs in expression matrix columns must match Cell_ID in metadata
malignant_ids   <- meta_malignant$Cell_ID
cells_in_matrix <- colnames(expr_mat)
cells_keep      <- intersect(malignant_ids, cells_in_matrix)

cat("Adult malignant cells in metadata  :", length(malignant_ids), "\n")
cat("Cells found in expression matrix   :", length(cells_keep), "\n")
missing_cells <- setdiff(malignant_ids, cells_in_matrix)
if (length(missing_cells) > 0) {
  cat("WARNING:", length(missing_cells), "metadata cells not found in matrix.\n")
  cat("First few missing:", paste(head(missing_cells, 5), collapse = ", "), "\n")
}

expr_subset <- expr_mat[, cells_keep]
meta_subset <- meta_malignant[match(cells_keep, meta_malignant$Cell_ID), ]
rownames(meta_subset) <- meta_subset$Cell_ID

cat("Final subset:", ncol(expr_subset), "cells x", nrow(expr_subset), "genes\n")

# =========================================
# 4. Build Seurat object
# =========================================
cat("\nBuilding Seurat object...\n")

# SS2 data is already log2(TPM/10 + 1) - store in 'data' slot directly.
# We still need a counts slot for Seurat; we back-transform to approximate
# counts by reversing the log transformation: counts = round(2^x * 10 - 10)
# This is imperfect but necessary for Seurat's architecture; all
# downstream analyses use the pre-normalised 'data' slot.
approx_counts <- round(pmax(2^expr_subset * 10 - 10, 0))
count_sparse  <- as(approx_counts, "sparseMatrix")

seurat_meta <- data.frame(
  row.names     = cells_keep,
  orig.ident    = meta_subset$Sample,
  GBMType       = meta_subset$GBMType,
  CellAssignment = meta_subset$CellAssignment,
  CrossSection  = meta_subset$CrossSection,
  GenesExpressed = meta_subset$GenesExpressed,
  GeneticSubclone = meta_subset$GeneticSubclone,
  MES_score     = meta_subset$MES_score,
  NPC_score     = meta_subset$NPC_score,
  AC_score      = meta_subset$AC_score,
  OPC_score     = meta_subset$OPC_score,
  Top_State     = meta_subset$Top_State,
  Dominant_State = meta_subset$Dominant_State,
  State_Confidence = meta_subset$State_Confidence,
  G1S_score     = meta_subset$G1S,
  G2M_score     = meta_subset$G2M,
  stringsAsFactors = FALSE
)

seurat_obj <- CreateSeuratObject(
  counts    = count_sparse,
  meta.data = seurat_meta,
  project   = OUTPUT_PREFIX
)

# Inject the pre-computed logTPM values into the data slot directly
seurat_obj <- SetAssayData(seurat_obj,
                           assay   = "RNA",
                           layer   = "data",
                           new.data = expr_subset[, cells_keep])

cat("Seurat object created:", ncol(seurat_obj), "cells,",
    nrow(seurat_obj), "genes\n")

# =========================================
# 5. Cell cycle regression using pre-computed scores
# =========================================
# SS2 metadata already provides Neftel's own G1S and G2M scores -
# use these directly rather than recomputing with Seurat's gene sets.
cat("\nCell cycle regression using pre-computed Neftel G1S/G2M scores...\n")
cat("G1S score summary:\n"); print(summary(seurat_obj$G1S_score))
cat("G2M score summary:\n"); print(summary(seurat_obj$G2M_score))

# ScaleData requires zero NAs in regression variables - impute with 0.
# These scores are mean-centred so 0 = population average, the safest
# imputation for the 43 cells where scores were not computed.
n_na_g1s <- sum(is.na(seurat_obj$G1S_score))
n_na_g2m <- sum(is.na(seurat_obj$G2M_score))
if (n_na_g1s > 0 || n_na_g2m > 0) {
  cat("Imputing", n_na_g1s, "NA G1S and", n_na_g2m,
      "NA G2M scores with 0 (mean-centred baseline).\n")
  seurat_obj$G1S_score[is.na(seurat_obj$G1S_score)] <- 0
  seurat_obj$G2M_score[is.na(seurat_obj$G2M_score)] <- 0
}

seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj,
                        vars.to.regress = c("G1S_score", "G2M_score"),
                        verbose = FALSE)

cat("Cell cycle regression complete.\n")

# =========================================
# 6. Batch correction (per tumor)
# =========================================
# Note: No ribosomal/MT/doublet filtering for SS2 -
# plate-based full-length protocol, already QC'd by Neftel et al.

correct_batch_effects_ss2 <- function(seurat_obj) {

  cat("Metadata columns:\n"); print(colnames(seurat_obj@meta.data))
  cat("Tumors (orig.ident):\n"); print(unique(seurat_obj$orig.ident))

  # Pre-integration UMAP for reference
  seurat_pre <- FindVariableFeatures(seurat_obj, nfeatures = 500)
  seurat_pre <- RunPCA(seurat_pre, verbose = FALSE)
  seurat_pre <- RunUMAP(seurat_pre, dims = 1:8, verbose = FALSE)

  plasma_cols <- viridis(n_distinct(seurat_pre$orig.ident), option = "plasma")

  p1 <- DimPlot(seurat_pre, group.by = "orig.ident", pt.size = 0.4,
                cols = plasma_cols) +
    ggtitle("UMAP Before Integration (by Tumor)")
  print(p1)

  p2 <- DimPlot(seurat_pre, group.by = "Dominant_State", pt.size = 0.4) +
    ggtitle("UMAP Before Integration (by Dominant Cell State)")
  print(p2)

  p3 <- DimPlot(seurat_pre, group.by = "GBMType", pt.size = 0.4) +
    ggtitle("UMAP Before Integration (by TCGA Subtype)")
  print(p3)

  # SCTransform per tumor then integrate
  sample_list <- SplitObject(seurat_obj, split.by = "orig.ident")
  sample_list <- lapply(sample_list, function(x) {
    x <- SCTransform(x, verbose = FALSE,
                     variable.features.n = 500,
                     vst.flavor = "v2")
    x
  })

  anchors <- FindIntegrationAnchors(object.list = sample_list,
                                    dims = 1:10, verbose = FALSE)

  tumor_counts <- table(seurat_obj$orig.ident)
  cat("Cells per tumor:\n"); print(tumor_counts)
  k_weight <- max(10, min(100, floor(min(tumor_counts) * 0.9)))
  cat("Dynamic k.weight:", k_weight, "\n")

  seurat_int <- IntegrateData(anchorset = anchors,
                              dims = 1:10,
                              k.weight = k_weight,
                              verbose = FALSE)
  seurat_int <- ScaleData(seurat_int, verbose = FALSE)
  seurat_int <- RunPCA(seurat_int, verbose = FALSE)
  seurat_int <- FindNeighbors(seurat_int, dims = 1:10)
  seurat_int <- FindClusters(seurat_int, resolution = 0.5)
  seurat_int <- RunUMAP(seurat_int, dims = 1:10, verbose = FALSE,
                        umap.method = "uwot", metric = "cosine")

  # Post-integration UMAPs
  plasma_cols2 <- viridis(n_distinct(seurat_int$orig.ident), option = "plasma")

  p4 <- DimPlot(seurat_int, group.by = "orig.ident", pt.size = 0.4,
                cols = plasma_cols2) +
    ggtitle("UMAP After Integration (by Tumor)")
  print(p4)

  p5 <- DimPlot(seurat_int, group.by = "Dominant_State", pt.size = 0.4) +
    ggtitle("UMAP After Integration (by Dominant Cell State)")
  print(p5)

  p6 <- DimPlot(seurat_int, group.by = "GBMType", pt.size = 0.4) +
    ggtitle("UMAP After Integration (by TCGA Subtype)")
  print(p6)

  p7 <- DimPlot(seurat_int, group.by = "Top_State", pt.size = 0.4,
                label = TRUE, repel = TRUE) +
    ggtitle("UMAP After Integration (by Top Cell State, labeled)")
  print(p7)

  return(seurat_int)
}

options(future.globals.maxSize = Inf)
cat("\nRunning batch correction...\n")
seurat_integrated <- correct_batch_effects_ss2(seurat_obj)

# =========================================
# 7. Elbow plot
# =========================================
generate_elbow_plot <- function(seurat_int, output_prefix) {
  seurat_int <- RunPCA(seurat_int, verbose = FALSE)
  ep <- ElbowPlot(seurat_int, ndims = 50) +
    labs(title = paste("Elbow Plot -", output_prefix),
         x = "Principal Components", y = "Standard Deviation")
  ggsave(paste0(output_prefix, "_ElbowPlot.pdf"), ep,
         width = 6.83, height = 6.41)
  print(ep)
}

generate_elbow_plot(seurat_integrated, OUTPUT_PREFIX)

# =========================================
# 8. Downstream analyses (UMAP clusters + feature plot)
# =========================================
downstream_analyses_ss2 <- function(seurat_int, gene_of_interest,
                                    output_prefix, dims = 15) {
  set.seed(10)
  DefaultAssay(seurat_int) <- "RNA"
  seurat_int <- FindVariableFeatures(seurat_int, nfeatures = 2000)
  seurat_int <- ScaleData(seurat_int, verbose = FALSE)
  seurat_int <- RunPCA(seurat_int, verbose = FALSE)
  seurat_int <- RunUMAP(seurat_int, dims = 1:dims, verbose = FALSE)
  set.seed(11); seurat_int <- FindNeighbors(seurat_int, dims = 1:dims)
  set.seed(12); seurat_int <- FindClusters(seurat_int, resolution = 0.5)
  seurat_int <- JoinLayers(seurat_int)

  set.seed(13)
  cluster_markers <- FindAllMarkers(seurat_int,
                                    only.pos = TRUE,
                                    min.pct = 0.1,
                                    logfc.threshold = 0.25)
  cat("Cluster markers found:", nrow(cluster_markers), "\n")

  # ---- UMAP helpers ----
  umap_data             <- as.data.frame(Embeddings(seurat_int, "umap"))
  umap_data$cluster_id  <- Idents(seurat_int)
  umap_data_mean        <- aggregate(. ~ cluster_id, data = umap_data, FUN = mean)

  cluster_order  <- sort(unique(as.numeric(as.character(Idents(seurat_int)))))
  n_cl           <- length(cluster_order)
  plasma_cols    <- colorRampPalette(viridis(100, direction = -1, option = "plasma"))(
                      round(n_cl / 0.8))

  # Cluster UMAP (labelled)
  set.seed(14)
  p_cl_lab <- ggplot(umap_data,
                     aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = plasma_cols) +
    geom_text(data = umap_data_mean,
              aes(label = cluster_id, x = umap_1, y = umap_2),
              color = "black", size = 3, fontface = "bold",
              check_overlap = TRUE) +
    theme_classic() +
    labs(title = "UMAP - Louvain clusters (labelled)",
         x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p_cl_lab)

  # Cluster UMAP (no labels)
  set.seed(15)
  p_cl_nol <- ggplot(umap_data,
                     aes(x = umap_1, y = umap_2, color = as.factor(cluster_id))) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = plasma_cols) +
    theme_classic() +
    labs(title = "UMAP - Louvain clusters (no labels)",
         x = "UMAP_1", y = "UMAP_2", color = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p_cl_nol)

  # UMAP by TCGA subtype
  umap_data$GBMType <- seurat_int@meta.data[rownames(umap_data), "GBMType"]
  subtype_cols <- c("Adult"       = "#4A90D9",
                    "Proneural"   = "#27AE60",
                    "Mesenchymal" = "#C0392B",
                    "Classical"   = "#F39C12",
                    "Mixed"       = "#8E44AD")
  p_sub <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = GBMType)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = subtype_cols, na.value = "grey80") +
    theme_classic() +
    labs(title = "UMAP - TCGA subtype",
         x = "UMAP_1", y = "UMAP_2", color = "TCGA Subtype") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p_sub)

  # UMAP by dominant cell state
  umap_data$Dominant_State <- seurat_int@meta.data[rownames(umap_data), "Dominant_State"]
  p_dom <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = Dominant_State)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_classic() +
    labs(title = "UMAP - Dominant cell state",
         x = "UMAP_1", y = "UMAP_2", color = "Dominant State") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p_dom)

  # UMAP by Top_State (4 canonical states only)
  umap_data$Top_State <- seurat_int@meta.data[rownames(umap_data), "Top_State"]
  state4_cols <- c("MES" = "#C0392B", "AC" = "#F39C12",
                   "NPC" = "#27AE60", "OPC" = "#4A90D9")
  p_top <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = Top_State)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = state4_cols) +
    theme_classic() +
    labs(title = "UMAP - Top canonical cell state",
         x = "UMAP_1", y = "UMAP_2", color = "Top State") +
    guides(color = guide_legend(override.aes = list(size = 3)))
  print(p_top)

  # Feature plot for GOI
  if (gene_of_interest %in% rownames(seurat_int)) {
    gene_cols <- c(scales::alpha("lightgray", 0.85),
                   scales::alpha("lightpink", 0.85),
                   scales::alpha("#FF6666", 0.85),
                   scales::alpha("#BC2727", 0.85),
                   scales::alpha("#660000", 0.85))
    set.seed(16)
    fp <- FeaturePlot(seurat_int, features = gene_of_interest,
                      min.cutoff = "q10", max.cutoff = "q90",
                      pt.size = 0.2, cols = gene_cols) +
      theme_classic() +
      labs(title = paste("UMAP -", gene_of_interest, "expression"),
           x = "UMAP_1", y = "UMAP_2")
    print(fp)
  } else {
    cat("WARNING: Gene", gene_of_interest, "not found in matrix.\n")
    cat("First 20 genes:", paste(head(rownames(seurat_int), 20), collapse = ", "), "\n")
  }

  # Write top markers
  if (nrow(cluster_markers) > 0) {
    top_markers <- cluster_markers %>%
      group_by(cluster) %>%
      top_n(n = 50, wt = avg_log2FC)
  } else {
    top_markers <- data.frame()
    warning("No cluster markers found - check cluster assignment.")
  }

  write.table(top_markers,
              file  = paste0(output_prefix, "_top_markers_per_cluster.tsv"),
              sep   = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  return(list(seurat_obj    = seurat_int,
              cluster_markers = cluster_markers,
              top_markers   = top_markers))
}

set.seed(42)
cat("\nRunning downstream analyses...\n")
gbm_results <- downstream_analyses_ss2(seurat_integrated, GOI,
                                       OUTPUT_PREFIX, dims = 20)

# Save workspace
save.image(file = paste0(OUTPUT_PREFIX, "_Seurat_processed.RData"))


##############################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 3: Gene Panel Plots, Spearman Correlation & Statistics
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#
# Prerequisite: gbm_results and seurat_obj must be in memory
# from Part 2 (Neftel_SS2_Part2_Seurat.R)
#
# Groupings:
#   (a) Louvain cluster
#   (b) Tumor (20 adult samples)
#   (c) Dominant cell state (all combinations e.g. AC+MES)
#   (d) Top canonical state (MES / AC / NPC / OPC)
#
# Also: Patient vs patient whole-transcriptome Spearman
#       correlation matrix across 20 adult tumors
##############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(writexl)
library(rstatix)
library(Seurat)
library(Matrix)

# =========================================================
# CONFIG - add or remove genes as needed
# =========================================================
GENE_PANEL_P3    <- c("RFX4", "ANXA1", "HOPX", "CHI3L1", "AXL", "HAVCR2")
OUTPUT_PREFIX_P3 <- "NeftelGBM_SS2_AdultMalignant_Part3"

# =========================================================
# STEP 1: Verify object and build colour palettes
# =========================================================
seurat_p3 <- gbm_results$seurat_obj   # final processed object from Part 2

cat(sprintf("SS2 adult malignant object: %d cells, %d clusters\n",
            ncol(seurat_p3),
            nlevels(Idents(seurat_p3))))
cat(sprintf("Tumors: %s\n",
            paste(sort(unique(seurat_p3$orig.ident)), collapse = ", ")))
# Cluster colours
cluster_order     <- sort(unique(as.numeric(as.character(Idents(seurat_p3)))))
n_cl              <- length(cluster_order)
cluster_cols      <- rev(viridis(round(n_cl / 0.8), option = "plasma"))[seq_len(n_cl)]
named_cluster_cols <- setNames(cluster_cols, as.character(cluster_order))

# Tumor colours (20 adult tumors)
tumor_order       <- sort(unique(seurat_p3$orig.ident))
tumor_cols        <- rev(viridis(length(tumor_order), option = "plasma"))
named_tumor_cols  <- setNames(tumor_cols, tumor_order)

# Top canonical state colours
state4_cols <- c(
  "MES" = "#C0392B",
  "AC"  = "#F39C12",
  "NPC" = "#27AE60",
  "OPC" = "#4A90D9"
)

# =========================================================
# STEP 2: Expression extractor
# =========================================================
get_expr_p3 <- function(gene, seurat_obj_p3) {
  DefaultAssay(seurat_obj_p3) <- "RNA"
  if (!gene %in% rownames(seurat_obj_p3)) {
    warning(sprintf("Gene not found: %s", gene))
    return(NULL)
  }
  cells <- Cells(seurat_obj_p3)
  meta  <- seurat_obj_p3@meta.data[cells, ]
  data.frame(
    gene           = gene,
    cluster        = factor(Idents(seurat_obj_p3),
                            levels = as.character(cluster_order)),
    expression     = as.numeric(
                       GetAssayData(seurat_obj_p3, layer = "data")[gene, cells]),
    tumor          = factor(meta$orig.ident,    levels = tumor_order),
    Dominant_State = factor(meta$Dominant_State),
    Top_State      = factor(meta$Top_State,
                            levels = c("MES", "AC", "NPC", "OPC")),
    stringsAsFactors = FALSE
  )
}

# =========================================================
# STEP 3: Build combined expression data frame
# =========================================================
cat("\nExtracting expression for Part 3 gene panel...\n")

panel_expr_p3  <- do.call(rbind, Filter(Negate(is.null),
  lapply(GENE_PANEL_P3, get_expr_p3, seurat_obj_p3 = seurat_p3)
))
genes_found_p3 <- unique(panel_expr_p3$gene)
missing_p3     <- setdiff(GENE_PANEL_P3, genes_found_p3)

cat(sprintf("Genes found    : %s\n", paste(genes_found_p3, collapse = ", ")))
if (length(missing_p3) > 0)
  cat(sprintf("Genes NOT found: %s\n", paste(missing_p3, collapse = ", ")))
cat(sprintf("Expression df  : %d rows (%d genes x %d cells)\n",
            nrow(panel_expr_p3), length(genes_found_p3), ncol(seurat_p3)))

# =========================================================
# STEP 4: Box + jitter plots - 5 groupings
# =========================================================

# ---- (a) By Louvain cluster ----
cat("\nPlotting by Louvain cluster...\n")
for (g in genes_found_p3) {
  df <- panel_expr_p3[panel_expr_p3$gene == g, ]
  p  <- ggplot(df, aes(x = cluster, y = expression, fill = cluster)) +
    geom_jitter(aes(color = cluster), alpha = 0.4, width = 0.3, height = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    scale_fill_manual(values  = named_cluster_cols) +
    scale_color_manual(values = named_cluster_cols) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none") +
    labs(x     = "Louvain Cluster",
         y     = paste("Expression of", g),
         title = paste(g, "- SS2 adult malignant GBM by Louvain cluster"))
  print(p)
}

# ---- (b) By tumor ----
cat("Plotting by tumor...\n")
for (g in genes_found_p3) {
  df <- panel_expr_p3[panel_expr_p3$gene == g, ]
  p  <- ggplot(df, aes(x = tumor, y = expression, fill = tumor)) +
    geom_jitter(aes(color = tumor), alpha = 0.4, width = 0.3, height = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    scale_fill_manual(values  = named_tumor_cols) +
    scale_color_manual(values = named_tumor_cols) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    labs(x     = "Tumor",
         y     = paste("Expression of", g),
         title = paste(g, "- SS2 adult malignant GBM by Tumor"))
  print(p)
}

# ---- (c) By dominant cell state (all combinations) ----
cat("Plotting by dominant cell state...\n")
for (g in genes_found_p3) {
  df <- panel_expr_p3[panel_expr_p3$gene == g &
                        !is.na(panel_expr_p3$Dominant_State), ]
  p  <- ggplot(df, aes(x = Dominant_State, y = expression,
                        fill = Dominant_State)) +
    geom_jitter(aes(color = Dominant_State), alpha = 0.4,
                width = 0.3, height = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    labs(x     = "Dominant Cell State",
         y     = paste("Expression of", g),
         title = paste(g, "- SS2 adult malignant GBM by Dominant Cell State"))
  print(p)
}

# ---- (d) By top canonical state (MES / AC / NPC / OPC) ----
cat("Plotting by top canonical state...\n")
for (g in genes_found_p3) {
  df <- panel_expr_p3[panel_expr_p3$gene == g &
                        !is.na(panel_expr_p3$Top_State), ]
  p  <- ggplot(df, aes(x = Top_State, y = expression, fill = Top_State)) +
    geom_jitter(aes(color = Top_State), alpha = 0.4, width = 0.3, height = 0) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    scale_fill_manual(values  = state4_cols) +
    scale_color_manual(values = state4_cols) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 30, vjust = 1, hjust = 1),
          legend.position = "none") +
    labs(x     = "Top Canonical State",
         y     = paste("Expression of", g),
         title = paste(g, "- SS2 adult malignant GBM: MES vs AC vs NPC vs OPC"))
  print(p)
}

# =========================================================
# STEP 5: Feature plots on UMAP
# =========================================================
cat("\nGenerating feature plots on SS2 UMAP...\n")
gene_colours <- c(scales::alpha("lightgray", 0.85),
                  scales::alpha("lightpink", 0.85),
                  scales::alpha("#FF6666",   0.85),
                  scales::alpha("#BC2727",   0.85),
                  scales::alpha("#660000",   0.85))

for (g in genes_found_p3) {
  if (!g %in% rownames(seurat_p3)) next
  fp <- FeaturePlot(seurat_p3, features = g,
                    min.cutoff = "q10", max.cutoff = "q90",
                    pt.size = 0.3, cols = gene_colours) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black")) +
    labs(title = sprintf("SS2 adult malignant GBM: %s", g),
         x = "UMAP_1", y = "UMAP_2")
  print(fp)
}

# =========================================================
# STEP 6: Patient vs patient whole-transcriptome Spearman
#         correlation matrix (20 adult tumors)
#
# For each tumor, compute mean normalised expression across
# all malignant cells → gene × tumor matrix.
# Spearman correlation between every tumor pair across the
# full expressed transcriptome (≥1% cells expressing).
# =========================================================
cat("\nComputing patient vs patient Spearman correlation (20 tumors)...\n")

DefaultAssay(seurat_p3) <- "RNA"
expr_data   <- GetAssayData(seurat_p3, layer = "data")
pct_nonzero <- rowMeans(expr_data > 0)
expr_filt   <- expr_data[pct_nonzero >= 0.01, ]

patients     <- sort(unique(seurat_p3$orig.ident))
pat_mean_mat <- do.call(cbind, lapply(patients, function(pt) {
  cells_pt <- which(seurat_p3$orig.ident == pt)
  Matrix::rowMeans(expr_filt[, cells_pt, drop = FALSE])
}))
colnames(pat_mean_mat) <- patients

cat(sprintf("  Genes used: %d  |  Patients: %d\n",
            nrow(pat_mean_mat), ncol(pat_mean_mat)))

pat_spearman <- cor(pat_mean_mat, method = "spearman")
cat("\nSpearman correlation matrix (20 adult tumors):\n")
print(round(pat_spearman, 4))

# Heatmap
pat_cor_df   <- as.data.frame(pat_spearman)
pat_cor_df$patient_row <- rownames(pat_cor_df)
pat_cor_long <- pivot_longer(pat_cor_df,
                             cols      = -patient_row,
                             names_to  = "patient_col",
                             values_to = "rho")
pat_cor_long$patient_row <- factor(pat_cor_long$patient_row, levels = patients)
pat_cor_long$patient_col <- factor(pat_cor_long$patient_col, levels = rev(patients))

rho_min <- floor(min(pat_spearman[lower.tri(pat_spearman)]) * 20) / 20

p_heatmap <- ggplot(pat_cor_long,
                    aes(x = patient_row, y = patient_col, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", rho)),
            size = 2.5, colour = "black") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
    limits  = c(rho_min, 1),
    name    = "Spearman ρ") +
  theme_minimal(base_size = 11) +
  theme(panel.border    = element_rect(fill = NA, colour = "black"),
        panel.grid      = element_blank(),
        axis.text       = element_text(colour = "black", size = 8),
        axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "right") +
  labs(title    = "Patient vs patient whole-transcriptome Spearman correlation",
       subtitle = sprintf(
         "SS2 adult malignant GBM | %d genes (≥1%% expressed) | mean per-tumor expression",
         nrow(pat_mean_mat)),
       x = NULL, y = NULL)
print(p_heatmap)

# =========================================================
# STEP 7: Statistics - 5 groupings
# (a) KW + Dunn - Louvain clusters
# (b) KW + Dunn - tumors
# (d) KW + Dunn - dominant cell states (all combinations)
# (e) KW + Dunn - top canonical states (MES/AC/NPC/OPC)
# =========================================================
cat("\nRunning statistics on Part 3 gene panel...\n")

run_stats_p3 <- function(gene, panel_df) {
  df <- panel_df[panel_df$gene == gene, ]
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df$cluster        <- droplevels(df$cluster)
  df$tumor          <- droplevels(df$tumor)
  df$Dominant_State <- droplevels(df$Dominant_State)
  df$Top_State      <- droplevels(df$Top_State)

  safe_kw   <- function(f, d) tryCatch(
    as.data.frame(kruskal_test(f, data = d)),
    error = function(e) data.frame(note = e$message))
  safe_dunn <- function(f, d) tryCatch(
    as.data.frame(dunn_test(f, data = d, p.adjust.method = "BH")),
    error = function(e) data.frame(note = e$message))
  desc <- function(grp_col, d) {
    d %>% group_by(.data[[grp_col]]) %>%
      summarise(n      = n(),
                median = median(expression),
                Q1     = quantile(expression, 0.25),
                Q3     = quantile(expression, 0.75),
                IQR    = IQR(expression),
                .groups = "drop")
  }

  df_dom <- df[!is.na(df$Dominant_State), ]
  df_top <- df[!is.na(df$Top_State), ]

  list(
    kw_clusters    = safe_kw(expression ~ cluster,        df),
    dunn_clusters  = safe_dunn(expression ~ cluster,      df),
    stats_clusters = desc("cluster",                      df),
    kw_tumors      = safe_kw(expression ~ tumor,          df),
    dunn_tumors    = safe_dunn(expression ~ tumor,        df),
    stats_tumors   = desc("tumor",                        df),
    kw_domstate    = safe_kw(expression ~ Dominant_State, df_dom),
    dunn_domstate  = safe_dunn(expression ~ Dominant_State, df_dom),
    stats_domstate = desc("Dominant_State",               df_dom),
    kw_topstate    = safe_kw(expression ~ Top_State,      df_top),
    dunn_topstate  = safe_dunn(expression ~ Top_State,    df_top),
    stats_topstate = desc("Top_State",                    df_top)
  )
}

all_stats_p3 <- lapply(setNames(genes_found_p3, genes_found_p3),
                       run_stats_p3,
                       panel_df = panel_expr_p3)

# =========================================================
# STEP 8: Export to Excel
# =========================================================
sheet_list_p3 <- list()
for (g in genes_found_p3) {
  s <- all_stats_p3[[g]]; if (is.null(s)) next
  add <- function(suffix, obj)
    sheet_list_p3[[substr(paste0(g, suffix), 1, 31)]] <<- as.data.frame(obj)
  add("_KW_Clusters",    s$kw_clusters)
  add("_Dunn_Clusters",  s$dunn_clusters)
  add("_Stats_Clusters", s$stats_clusters)
  add("_KW_Tumors",      s$kw_tumors)
  add("_Dunn_Tumors",    s$dunn_tumors)
  add("_Stats_Tumors",   s$stats_tumors)
  add("_KW_DomState",    s$kw_domstate)
  add("_Dunn_DomState",  s$dunn_domstate)
  add("_Stats_DomState", s$stats_domstate)
  add("_KW_TopState",    s$kw_topstate)
  add("_Dunn_TopState",  s$dunn_topstate)
  add("_Stats_TopState", s$stats_topstate)
}
sheet_list_p3[["Spearman_tumor_vs_tumor"]] <-
  as.data.frame(round(pat_spearman, 4))

out_xlsx_p3 <- paste0(OUTPUT_PREFIX_P3, "_stats.xlsx")
write_xlsx(sheet_list_p3, out_xlsx_p3)
cat(sprintf("\nStatistics written to: %s\n", out_xlsx_p3))

# =========================================================
# STEP 9: Summary
# =========================================================
cat("\n===== Part 3 complete =====\n")
cat(sprintf("  Genes plotted  : %s\n", paste(genes_found_p3, collapse = ", ")))
cat(sprintf("  Cells          : %d\n", ncol(seurat_p3)))
cat(sprintf("  Clusters       : %d\n", nlevels(Idents(seurat_p3))))
cat(sprintf("  Tumors         : %d\n", length(patients)))
cat(sprintf("  Output xlsx    : %s\n", out_xlsx_p3))
cat("\nReady for Part 4 (CEP-IP inflection point analysis)\n")


###############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4 (Revised): State-Stratified Dual-Filter Correlation Analysis
#                   Automated loop over full gene longlist
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#          Adult malignant cells only
#
# Prerequisite: gbm_results must be in memory from Part 2, OR load workspace:
#   load("NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData")
###############################################################################

setwd("C:/Users/...#directory")
load("NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData") # Load RData if haven't already
# RAM Cleanup - run before starting Part 4 (State-Stratified) if needed
'''
# Frees memory from previous interrupted run
# ---- 1. Stop any lingering parallel clusters ----
library(parallel)
library(doParallel)

# Kill all registered clusters gracefully
tryCatch({
  registered <- showConnections(all = FALSE)
  if (nrow(registered) > 0) {
    cat("Open connections found - closing...\n")
    closeAllConnections()
    cat("Done.\n")
  } else {
    cat("No open connections found.\n")
  }
}, error = function(e) cat("Connection check skipped:", e$message, "\n"))

tryCatch({
  stopImplicitCluster()
  cat("Implicit cluster stopped.\n")
}, error = function(e) cat("No implicit cluster to stop.\n"))

# Explicitly stop any cl_global left from previous run
if (exists("cl_global")) {
  tryCatch({
    stopCluster(cl_global)
    rm(cl_global)
    cat("cl_global stopped and removed.\n")
  }, error = function(e) cat("cl_global already stopped.\n"))
}

# ---- 2. Remove large objects from previous Part 4 run ----
large_objects <- c(
  "full_expr_mat",   # full gene x cell matrix - largest object
  "pool_mat",        # per-state subset matrix
  "pool_mat_kd",     # Kendall candidate subset
  "sp_df", "kd_df",  # per-state correlation results
  "sp_sheets", "kd_sheets", "df_sheets",  # per-GOI workbook lists
  "all_spearman", "all_kendall",          # per-patient results (old approach)
  "sp_vals", "kd_vals",                   # raw correlation vectors
  "state_cell_idx",                       # cell index lists
  "master_rows", "master_df",             # master summary
  "dual_filter_all", "dual_filter_per_patient",
  "gene_frequency", "dfg_ranking",
  "panel_expr", "panel_expr_p3",          # gene panel expression dfs
  "d", "data_list"                        # temp objects
)

removed <- c()
for (obj in large_objects) {
  if (exists(obj)) {
    size_mb <- round(object.size(get(obj)) / 1024^2, 1)
    rm(list = obj, envir = .GlobalEnv)
    removed <- c(removed, sprintf("%s (%.1f MB)", obj, size_mb))
  }
}

if (length(removed) > 0) {
  cat("\nRemoved objects:\n")
  for (r in removed) cat(" ", r, "\n")
} else {
  cat("\nNo large Part 4 objects found in environment.\n")
}

# ---- 3. Force garbage collection - multiple passes ----
cat("\nRunning garbage collection...\n")
before_mb <- round(sum(gc(reset = TRUE)[, 2]) * 8 / 1024, 1)  # Ncells x 8 bytes
gc(); gc(); gc()  # three passes to fully release
after_mb  <- round(sum(gc()[, 2]) * 8 / 1024, 1)
cat(sprintf("GC complete. Approximate R heap: %.1f MB\n", after_mb))

# ---- 4. Report current environment ----
cat("\nRemaining objects in global environment:\n")
env_objects <- ls(envir = .GlobalEnv)
if (length(env_objects) == 0) {
  cat("  (empty)\n")
} else {
  sizes <- sapply(env_objects, function(x) {
    tryCatch(object.size(get(x)), error = function(e) 0)
  })
  obj_df <- data.frame(
    Object  = env_objects,
    Size_MB = round(as.numeric(sizes) / 1024^2, 2)
  )
  obj_df <- obj_df[order(-obj_df$Size_MB), ]
  print(obj_df, row.names = FALSE)
}

cat("\nRAM cleanup complete. Safe to run Part 4 (State-Stratified).\n")
'''

library(Seurat)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(writexl)

# ---------------------------------------------------------------------------
# Global parallelisation - register once, reuse throughout entire run.
# Uses all cores minus 1 to keep the OS responsive.
# ---------------------------------------------------------------------------
NUM_CORES <- max(1L, detectCores() - 1L)
cl_global <- makeCluster(NUM_CORES)
registerDoParallel(cl_global)
cat(sprintf("Parallel backend registered: %d cores\n\n", NUM_CORES))

# ===========================================================================
# CONFIG
# ===========================================================================
# Total: 100 genes
# Gene longlist: 659 upregulated DEGs in GBM vs normal brain
# Source: Prasad B, Tian Y, Li X. Large-Scale Analysis Reveals Gene Signature
#   for Survival Prediction in Primary Glioblastoma.
#   Mol Neurobiol. 2020 Dec;57(12):5235-5246. doi:10.1007/s12035-020-02088-w
# Derived from Supplementary Table S2 (Additional file 2):
#   maxFC > 0 AND Pval_Bonf < 0.05 → 686 genes upregulated in GBM
#   minus immune genes (10 HLA + CD74 + CD14/CD163/CD68 +
#   8 complement + FCER1G/FCGR2A/LY96) → 659 final genes
GENE_LONGLIST <- c(
    # ── Prasad et al. 2020 - 659 GBM upregulated DEGs ──────────────────────
    "A2M", "ABCA1", "ABCC3", "ACKR3", "ACTA2", "ACTL6A", "ACTN1", "ADAM12", "ADAM9", "ADAMTS9",
    "ADGRL4", "ADM", "AEBP1", "AFAP1L1", "AGTRAP", "AJUBA", "AKT3", "ANGPT2", "ANGPTL2", "ANKRD10-IT1",
    "ANO6", "ANTXR2", "ANXA1", "ANXA2", "ANXA2P2", "ANXA5", "APLN", "APOC1", "AQP1", "ARHGAP18",
    "ARSJ", "ASCL1", "ASPM", "ATP10D", "AURKA", "B3GNT5", "BACE2", "BARD1", "BAZ1A", "BCAN",
    "BCAT1", "BCHE", "BCL6", "BEST3", "BGN", "BICD1", "BIRC5", "BRI3", "BST2", "BTG1",
    "BTG3", "BTN3A2", "BUB1", "BUB1B", "C11orf96", "C1orf162", "C1orf226", "C1orf61", "C21orf62", "CA12",
    "CA3", "CALD1", "CALU", "CAPG", "CARD16", "CASK", "CASP1", "CAV1", "CAV2", "CBX3",
    "CCDC102B", "CCDC18-AS1", "CCDC50", "CCDC80", "CCN1", "CCN2", "CCN4", "CCNA2", "CCNB1", "CCNB2",
    "CCND2", "CD151", "CD276", "CD44", "CD58", "CD63", "CD93", "CD99", "CDC20", "CDCA7",
    "CDCA7L", "CDH11", "CDK1", "CDK2", "CDK4", "CDK6", "CDKN1A", "CDKN2C", "CDKN3", "CEBPD",
    "CENPF", "CENPK", "CENPU", "CEP55", "CFI", "CHEK1", "CHEK2", "CHI3L1", "CHI3L2", "CHRNB4",
    "CKS2", "CLEC2B", "CLEC5A", "CLIC1", "CLIC4", "CMTM3", "CMTM6", "CNIH4", "CNN3", "COL15A1",
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "COL8A1", "CP", "CPVL", "CPXM1", "CRISPLD1", "CRNDE", "CSRP2", "CSTA", "CTHRC1", "CTSC",
    "CTSS", "CXCL10", "CXCR4", "DAB2", "DBF4", "DBI", "DDR2", "DDX39A", "DEPDC1B", "DEPP1",
    "DESI2", "DLEU2", "DLGAP5", "DMRTA2", "DPH3", "DPY19L1", "DPYD", "DPYSL3", "DRAM1", "DSE",
    "DTL", "DTNA", "DTX3L", "DUSP4", "DUSP6", "DYNLT1", "E2F7", "ECM2", "ECSCR", "ECT2",
    "EDNRA", "EFEMP1", "EFHC1", "EGFEM1P", "EGFR", "EIF4E2", "EIF4EBP1", "ELN", "EMILIN1", "EMILIN2",
    "EMP1", "EMP3", "ENPEP", "ERAP2", "ERI1", "ESM1", "ETS1", "ETV1", "EXOC4", "EZH2",
    "F13A1", "F2R", "FAM111A", "FAM20C", "FANCD2", "FANCI", "FBXO32", "FBXO5", "FIGN", "FILIP1L",
    "FJX1", "FKBP10", "FLNA", "FN1", "FNDC3B", "FOXD1", "FOXM1", "FPR3", "FREM2", "FRMD3",
    "FSTL1", "FTX", "FZD7", "GABPB1-AS1", "GADD45A", "GART", "GAS1", "GAS2L3", "GBE1", "GBP1",
    "GBP2", "GBP3", "GDF15", "GEM", "GFAP", "GIMAP2", "GINS1", "GJC1", "GLIPR1", "GLIPR2",
    "GLIS2", "GLIS3", "GMNN", "GNA13", "GNAI3", "GNAS", "GNG12", "GNG5", "GNS", "GPNMB",
    "GPR65", "GPX3", "GPX7", "GPX8", "GRN", "GUSB", "H2AFV", "HAS2", "HAT1", "HAUS1",
    "HELLS", "HES1", "HIF1A", "HILPDA", "HK2", "HMG20B", "HMGB2", "HMGN1", "HMOX1", "HOTAIRM1",
    "HOXB3", "HOXC6", "HS3ST3B1", "HSPA6", "HSPB1", "HSPG2", "IBSP", "ID3", "ID4", "IDH1",
    "IFI16", "IFI30", "IFI44", "IFRD1", "IGDCC4", "IGF2BP2", "IGF2BP3", "IGFBP2", "IGFBP3", "IGFBP4",
    "IGFBP5", "IGFBP7", "IL13RA2", "IL1RAP", "IQGAP1", "IQGAP2", "IRAK1", "ITGA7", "ITGB1", "ITGB2",
    "ITGB3BP", "ITPRIPL2", "JAG1", "JAM2", "KDELR2", "KIAA0040", "KIAA0754", "KIF11", "KIF14", "KIF15",
    "KIF1B", "KIF20A", "KIF23", "KIF2C", "KIF4A", "KLHDC8A", "KLHL7", "KNL1", "LAMA4", "LAMB1",
    "LAMB2", "LAMC1", "LAP3", "LAPTM5", "LATS2", "LDLRAD3", "LGALS1", "LGALS3", "LGALS3BP", "LIMA1",
    "LINC01088", "LINC01279", "LMNB1", "LOC101929876", "LOC154761", "LOX", "LOXL1", "LOXL2", "LPL", "LRR1",
    "LRRC17", "LSM5", "LSM8", "LTBP3", "LTF", "LUM", "LYPD1", "LYZ", "MAD2L1", "MAFB",
    "MAFF", "MALAT1", "MAML2", "MAP3K7CL", "MATN2", "MBD6", "MCM2", "MCM3", "MCM7", "MCUB",
    "MDK", "MDM2", "MELK", "MEOX2", "MEST", "METTL7B", "MGP", "MIR570HG", "MMP2", "MMP7",
    "MMP9", "MRC2", "MS4A4A", "MS4A6A", "MS4A7", "MSI2", "MSN", "MSR1", "MTHFD2", "MVB12A",
    "MXRA5", "MYC", "MYL12A", "MYO1B", "MYOF", "N4BP2", "N4BP2L2", "NAMPT", "NAPSB", "NASP",
    "NCAPG", "NDC80", "NEDD1", "NEDD9", "NEMP1", "NES", "NFKBIZ", "NIBAN1", "NID1", "NID2",
    "NLRC5", "NMB", "NMI", "NNMT", "NOTCH2NLA", "NOX4", "NPC2", "NPL", "NPNT", "NRP1",
    "NRP2", "NT5DC2", "NTN1", "NUF2", "NUP160", "NUPR1", "NUSAP1", "ODC1", "P4HB", "PABPC1L",
    "PALLD", "PARP9", "PBK", "PBX3", "PCDHB16", "PCDHB7", "PCLAF", "PCNA", "PCOLCE", "PCOLCE2",
    "PDIA4", "PDIA5", "PDLIM1", "PDLIM3", "PDPN", "PECAM1", "PFN1", "PGGT1B", "PGM2", "PHF14",
    "PHLDA1", "PIMREG", "PLA2G2A", "PLA2G5", "PLAT", "PLAU", "PLEKHA4", "PLEKHA8P1", "PLIN2", "PLOD1",
    "PLOD2", "PLP2", "PLSCR1", "PLTP", "POGLUT3", "POSTN", "PPIC", "PPP1R14B", "PPP1R18", "PPP1R3B",
    "PRC1", "PRDX4", "PRKD3", "PROS1", "PRR11", "PRRX1", "PRSS23", "PSMB8", "PSMB9", "PSPH",
    "PTBP1", "PTPN12", "PTPRC", "PTPRZ1", "PTTG1", "PTTG3P", "PTX3", "PXDN", "PYGL", "QKI",
    "RAB13", "RAD51AP1", "RBBP8", "RBM47", "RBMS1", "RBP1", "RCAN1", "RCN1", "RDH10", "RELL1",
    "RFC4", "RFX4", "RGMA", "RGS1", "RHOJ", "RHPN2", "RIT1", "RND3", "RNF135", "RNF180",
    "RNF213", "RP2", "RPL10", "RPL18A", "RPLP0", "RPN2", "RPS13", "RPS2P32", "RPS8", "RRM2",
    "S100A10", "S100A11", "S100A4", "S100A6", "S100A8", "S1PR3", "SALL3", "SAMD9L", "SAP30", "SCARA3",
    "SEC24D", "SEC61A1", "SEC61G", "SERPINA1", "SERPINA3", "SERPINE1", "SERPING1", "SERPINH1", "SFRP4", "SH3PXD2B",
    "SHC1", "SHCBP1", "SHMT2", "SHOX2", "SINHCAF", "SLC16A1", "SLC16A3", "SLC16A4", "SLC2A10", "SLC39A14",
    "SLC39A8", "SLC40A1", "SLC43A3", "SLC47A2", "SLPI", "SMAD1", "SMC4", "SMC5", "SMIM3", "SNAI2",
    "SNCAIP", "SNRPG", "SOAT1", "SOCS2", "SOCS3", "SOX11", "SOX2", "SOX4", "SOX9", "SPARC",
    "SPOCD1", "SPPL2A", "SPRY1", "SPRY4", "SRGAP1", "SRGN", "SRPX", "SRPX2", "STAB1", "STC1",
    "STEAP3", "STK17A", "STMP1", "STON1", "SULF2", "SYNE2", "SYNPO", "TACC3", "TAGLN", "TAGLN2",
    "TANC1", "TAP1", "TCF12", "TCIM", "TEAD2", "TENT2", "TENT5A", "TFAP2A", "TFPI", "TGFB1I1",
    "TGFBI", "TGFBR1", "TGIF1", "THOC2", "TIMELESS", "TIMP1", "TIMP4", "TMEM176B", "TMEM248", "TMEM255A",
    "TMEM45A", "TMSB15A", "TMX1", "TNC", "TNFRSF10B", "TNFRSF12A", "TNFRSF19", "TNFRSF1A", "TNFSF13B", "TNPO1",
    "TOP2A", "TP53", "TP53I3", "TP53INP1", "TPM2", "TPM4", "TPST1", "TPX2", "TREM1", "TRIB2",
    "TRIM14", "TRIM22", "TRIM47", "TRIM5", "TRIP6", "TRPM8", "TSPAN12", "TSPAN6", "TTK", "TUBB6",
    "TWSG1", "TXNDC17", "TYMS", "TYROBP", "UBE2C", "UHRF1", "VAMP8", "VANGL2", "VCAM1", "VCAN",
    "VCL", "VEGFA", "VIM", "VMP1", "VOPP1", "VSIG4", "WEE1", "WLS", "WNT5A", "WSB1",
    "WTAP", "WWTR1", "XPR1", "YAP1", "YBX1", "YBX3", "ZBTB20", "ZBTB42", "ZC3HAV1", "ZFP36L1",
    "ZFP36L2", "ZIC1", "ZNF217", "ZNF300", "ZNF367", "ZNF521", "ZNF548", "ZWINT", "ZYX"
    # ── Additional GBM-related genes absent from Prasad et al. 2020 ─────────────────────
    # Surface / Membrane Antigens & Receptors
    #"EPHA2", "ERBB2", "CD70", "PROM1", "PDGFRA",
    # Proliferation / Cell Cycle
    #"MKI67",
    # Transcription Factors & Regulators
    #"OLIG2", "HOXA3", "EN1", "FOXD3", "STAT3",
    # Other GBM-Enriched / Specific
    #"ARHGAP9", "ARHGAP30", "CLEC7A", "MAN2B1", "ARPC1B",
    #"PLB1", "SOD2", "CALM3", "SIGLEC9", "MYO1F",
    #"AGBL2", "BLOC1S6", "MAP1A", "ZSWIM5",
    # Cancer-Testis Antigens
    #"MAGEA1", "MAGEA3", "MAGEA4", "MAGEA12",
    #"CTAG1B", "FBXO39", "CTCFL", "CASC5", "DDX43", "ACTL8",
    #"OIP5", "XAGE3", "RGS22"
)

TOP_STATES         <- c("MES", "AC", "NPC", "OPC")

# Quality filter: exclude patient-state pairs with fewer than this many cells
MIN_PT_STATE_CELLS <- 5

# GOI expression gate within each state pool
GOI_PCT_THRESHOLD  <- 0.25   # GOI must be expressed in >= 25% of pooled cells
GOI_EXPRESSOR_MEAN <- 0.30   # mean expression among GOI-expressing cells >= 0.3

# Partner gene filter within each state pool
MIN_CELLS          <- 10     # partner gene expressed in >= 10 cells in pool

# Dual-filter thresholds
SP_THRESHOLD       <- 0.6
KD_THRESHOLD       <- 0.5

OUTPUT_PREFIX      <- "NeftelGBM_SS2_StateStratified"

# ===========================================================================
# SECTION 1: Load object, extract full matrix once, build state cell indices
# ===========================================================================
seurat_p4 <- gbm_results$seurat_obj
DefaultAssay(seurat_p4) <- "RNA"

cat("=================================================================\n")
cat("STATE-STRATIFIED DUAL-FILTER CORRELATION ANALYSIS\n")
cat("SS2 Adult Malignant GBM - Neftel et al. 2019\n")
cat("Genes in longlist      :", length(GENE_LONGLIST), "\n")
cat("States analysed        :", paste(TOP_STATES, collapse = ", "), "\n")
cat("Total cells            :", ncol(seurat_p4), "\n")
cat("Min cells per pt-state :", MIN_PT_STATE_CELLS, "\n")
cat("GOI pct threshold      :", GOI_PCT_THRESHOLD * 100, "% (>=25% expressing to avoid tie inflation)\n")
cat("GOI expressor mean     :", GOI_EXPRESSOR_MEAN, "\n")
cat("Min partner cells      :", MIN_CELLS, "\n")
cat("Spearman threshold     :", SP_THRESHOLD, "\n")
cat("Kendall threshold      :", KD_THRESHOLD, "\n")
cat("=================================================================\n\n")

# Extract full normalised matrix once
cat("Extracting full normalised expression matrix (once)...\n")
full_expr_mat <- as.matrix(GetAssayData(seurat_p4, layer = "data",
                                        assay = "RNA"))
cat("Matrix:", nrow(full_expr_mat), "genes x", ncol(full_expr_mat),
    "cells\n\n")

# Per-patient cell counts per state - for quality filter reporting
pt_state_counts <- as.data.frame(table(
    Patient   = seurat_p4$orig.ident,
    Top_State = seurat_p4$Top_State
))
cat("Patient-state cell counts (all):\n")
print(pt_state_counts %>%
          filter(Freq > 0) %>%
          tidyr::pivot_wider(names_from = Top_State, values_from = Freq,
                             values_fill = 0) %>%
          arrange(Patient),
      row.names = FALSE)

# Build per-state pooled cell indices AFTER applying MIN_PT_STATE_CELLS filter
cat(sprintf("\nApplying min patient-state filter (>= %d cells)...\n",
            MIN_PT_STATE_CELLS))

state_cell_idx <- list()
for (st in TOP_STATES) {
    # patients that have >= MIN_PT_STATE_CELLS cells in this state
    eligible_pts <- pt_state_counts %>%
        filter(Top_State == st, Freq >= MIN_PT_STATE_CELLS) %>%
        pull(Patient)
    
    # cell indices belonging to eligible patients AND this state
    idx <- which(seurat_p4$Top_State == st &
                     seurat_p4$orig.ident %in% eligible_pts)
    
    n_excluded <- sum(seurat_p4$Top_State == st) - length(idx)
    cat(sprintf("  %s: %d cells total | %d excluded (tiny pt-state pairs) | %d retained | from %d patients\n",
                st, sum(seurat_p4$Top_State == st),
                n_excluded, length(idx), length(eligible_pts)))
    
    state_cell_idx[[st]] <- idx
}
cat("\n")

# Check which longlist genes are in the matrix
genes_in_matrix <- intersect(GENE_LONGLIST, rownames(full_expr_mat))
genes_missing   <- setdiff(GENE_LONGLIST, rownames(full_expr_mat))
cat("Longlist genes found in matrix:", length(genes_in_matrix), "\n")
if (length(genes_missing) > 0)
    cat("NOT in matrix:", paste(genes_missing, collapse = ", "), "\n")
cat("\n")

# ===========================================================================
# SECTION 2: Helper functions
# ===========================================================================

# ---- 2a. GOI expression gate within a state pool ----
check_goi_gate <- function(goi_vec, state_label) {
    n_cells     <- length(goi_vec)
    n_expr      <- sum(goi_vec > 0)
    pct_expr    <- n_expr / n_cells
    mean_expr   <- if (n_expr > 0) mean(goi_vec[goi_vec > 0]) else 0
    
    passes <- pct_expr >= GOI_PCT_THRESHOLD & mean_expr >= GOI_EXPRESSOR_MEAN
    
    cat(sprintf("    Gate [%s]: cells=%d | expressing=%d (%.1f%%) | mean_expr(>0)=%.3f | %s\n",
                state_label, n_cells, n_expr, pct_expr * 100, mean_expr,
                if (passes) "PASS" else "SKIP"))
    return(passes)
}

# ---- 2b. Parallelised Spearman ----
# Splits genes into NUM_CORES chunks; each worker computes cor() on its chunk.
# Faster than the single-threaded vectorised cor() for large gene x cell matrices.
run_spearman_state <- function(state_mat, goi_vec, goi, state_label) {
    n_genes <- nrow(state_mat)
    cat(sprintf("    [Spearman] %s: %d genes x %d cells, %d cores ... ",
                state_label, n_genes, ncol(state_mat), NUM_CORES))
    
    gene_rank  <- rank(goi_vec)
    gene_chunks <- split(seq_len(n_genes),
                         cut(seq_len(n_genes), breaks = NUM_CORES, labels = FALSE))
    
    sp_vals <- unlist(foreach(
        chunk    = gene_chunks,
        .combine = "c",
        .export  = c("gene_rank", "state_mat")
    ) %dopar% {
        cor(gene_rank, t(state_mat[chunk, , drop = FALSE]),
            method = "spearman")[1, ]
    })
    names(sp_vals) <- rownames(state_mat)
    
    df <- data.frame(
        Gene                 = names(sp_vals),
        Spearman_Correlation = as.numeric(sp_vals),
        State                = state_label,
        stringsAsFactors     = FALSE
    ) %>%
        filter(Gene != goi) %>%
        arrange(desc(Spearman_Correlation))
    
    cat(sprintf("done | Pass(>=0.6): %d\n",
                sum(df$Spearman_Correlation >= 0.6, na.rm = TRUE)))
    return(df)
}

# ---- 2c. Parallelised Kendall (uses global cluster cl_global) ----
run_kendall_state <- function(state_mat, goi_vec, goi, state_label) {
    n_genes <- nrow(state_mat)
    cat(sprintf("    [Kendall]  %s: %d genes x %d cells, %d cores ... ",
                state_label, n_genes, ncol(state_mat), NUM_CORES))
    
    kd_vals <- foreach(
        i        = seq_len(n_genes),
        .combine = "c",
        .export  = c("goi_vec", "state_mat")
    ) %dopar% {
        cor(goi_vec, state_mat[i, ], method = "kendall")
    }
    
    df <- data.frame(
        Gene                = rownames(state_mat),
        Kendall_Correlation = kd_vals,
        State               = state_label,
        stringsAsFactors    = FALSE
    ) %>%
        filter(Gene != goi) %>%
        arrange(desc(Kendall_Correlation))
    
    cat(sprintf("done | Pass(>=0.5): %d\n",
                sum(df$Kendall_Correlation >= 0.5, na.rm = TRUE)))
    return(df)
}

# ---- 2d. Dual-filter for one state ----
dual_filter_state <- function(sp_df, kd_df, state_label) {
    hits <- inner_join(
        sp_df %>% select(Gene, Spearman_Correlation),
        kd_df %>% select(Gene, Kendall_Correlation),
        by = "Gene"
    ) %>%
        filter(Spearman_Correlation >= SP_THRESHOLD,
               Kendall_Correlation  >= KD_THRESHOLD) %>%
        mutate(State = state_label) %>%
        arrange(desc(Spearman_Correlation))
    
    cat(sprintf("    DualFilter [%s]: %d DFGs (Sp>=%.1f & Kd>=%.1f)\n",
                state_label, nrow(hits), SP_THRESHOLD, KD_THRESHOLD))
    return(hits)
}

# ===========================================================================
# SECTION 3: Master loop - all longlist genes x all states
# ===========================================================================
# Storage for master summary
master_rows <- list()

for (goi in genes_in_matrix) {
    
    cat("################################################################\n")
    cat(sprintf("GOI: %s\n", goi))
    cat("################################################################\n")
    
    goi_row <- full_expr_mat[goi, ]
    
    # Per-state storage
    sp_sheets  <- list()   # Spearman results per state
    kd_sheets  <- list()   # Kendall results per state
    df_sheets  <- list()   # DualFilter results per state
    dfg_counts <- setNames(rep(NA_integer_, length(TOP_STATES)), TOP_STATES)
    
    for (st in TOP_STATES) {
        
        cat(sprintf("  -- State: %s (%d pooled cells) --\n",
                    st, length(state_cell_idx[[st]])))
        
        idx     <- state_cell_idx[[st]]
        goi_vec <- goi_row[idx]
        
        # GOI expression gate
        passes <- check_goi_gate(goi_vec, st)
        if (!passes) {
            sp_sheets[[st]] <- data.frame(Gene = character(), Spearman_Correlation = numeric(),
                                          State = character(), Note = character())
            kd_sheets[[st]] <- data.frame(Gene = character(), Kendall_Correlation = numeric(),
                                          State = character(), Note = character())
            df_sheets[[st]] <- data.frame(Gene = character(), Spearman_Correlation = numeric(),
                                          Kendall_Correlation = numeric(), State = character(),
                                          Note = character())
            dfg_counts[st]  <- 0L
            next
        }
        
        # Subset expression matrix to this state pool
        pool_mat <- full_expr_mat[, idx, drop = FALSE]
        
        # Filter partner genes expressed in >= MIN_CELLS cells in this pool
        genes_keep <- rowSums(pool_mat > 0, na.rm = TRUE) >= MIN_CELLS
        pool_mat   <- pool_mat[genes_keep, , drop = FALSE]
        cat(sprintf("    Partner genes after filter (>=%d cells): %d\n",
                    MIN_CELLS, nrow(pool_mat)))
        
        # Spearman
        sp_df <- run_spearman_state(pool_mat, goi_vec, goi, st)
        
        # Gate: only run Kendall if at least 1 gene clears Spearman T2 (>= SP_THRESHOLD)
        # Rationale: dual-filter requires Sp >= SP_THRESHOLD, so zero T2 Spearman hits
        # means zero possible DFGs regardless of Kendall - skip expensive parallelisation
        sp_t2_count <- sum(sp_df$Spearman_Correlation >= SP_THRESHOLD, na.rm = TRUE)
        
        if (sp_t2_count == 0) {
            cat(sprintf("    [Kendall]  %s: SKIPPED - 0 Spearman T2 hits, no DFGs possible\n", st))
            kd_df <- data.frame(Gene = character(), Kendall_Correlation = numeric(),
                                State = character(), stringsAsFactors = FALSE)
            hits  <- data.frame(Gene = character(), Spearman_Correlation = numeric(),
                                Kendall_Correlation = numeric(), State = character(),
                                stringsAsFactors = FALSE)
        } else {
            cat(sprintf("    [Kendall]  proceeding - %d Spearman T2 candidate(s)\n", sp_t2_count))
            # Run Kendall only on the subset of genes passing Spearman T2 -
            # further reduces Kendall compute when T2 count is small
            sp_t2_genes <- sp_df$Gene[sp_df$Spearman_Correlation >= SP_THRESHOLD]
            pool_mat_kd <- pool_mat[rownames(pool_mat) %in% sp_t2_genes, , drop = FALSE]
            kd_df <- run_kendall_state(pool_mat_kd, goi_vec, goi, st)
            hits  <- dual_filter_state(sp_df, kd_df, st)
        }
        
        # Print top 20 DFGs
        if (nrow(hits) > 0) {
            cat(sprintf("    Top DFGs [%s]:\n", st))
            cat(sprintf("    %-20s %12s %12s\n", "Gene", "Spearman", "Kendall"))
            cat("    ", strrep("-", 48), "\n")
            top20 <- head(hits, 20)
            for (i in seq_len(nrow(top20))) {
                r <- top20[i, ]
                cat(sprintf("    %-20s %12.4f %12.4f\n",
                            r$Gene, r$Spearman_Correlation, r$Kendall_Correlation))
            }
        }
        cat("\n")
        
        sp_sheets[[st]] <- sp_df
        kd_sheets[[st]] <- kd_df
        df_sheets[[st]] <- hits
        dfg_counts[st]  <- nrow(hits)
    }
    
    # ---- Export per-GOI workbooks ----
    
    # Spearman workbook - 4 sheets (one per state)
    sp_file <- paste0(OUTPUT_PREFIX, "_", goi, "_Spearman.xlsx")
    write_xlsx(lapply(sp_sheets, as.data.frame), path = sp_file)
    cat(sprintf("  Spearman saved : %s\n", sp_file))
    
    # Kendall workbook - 4 sheets (one per state)
    kd_file <- paste0(OUTPUT_PREFIX, "_", goi, "_Kendall.xlsx")
    write_xlsx(lapply(kd_sheets, as.data.frame), path = kd_file)
    cat(sprintf("  Kendall saved  : %s\n", kd_file))
    
    # DualFilter summary workbook - 4 state sheets + combined
    all_hits <- bind_rows(df_sheets)
    dfg_summary_sheets <- c(
        list(All_States_Combined = as.data.frame(all_hits)),
        lapply(df_sheets, as.data.frame)
    )
    df_file <- paste0(OUTPUT_PREFIX, "_", goi, "_DualFilter.xlsx")
    write_xlsx(dfg_summary_sheets, path = df_file)
    cat(sprintf("  DualFilter saved: %s\n", df_file))
    
    # ---- Store master row ----
    master_rows[[goi]] <- data.frame(
        GOI          = goi,
        DFG_MES      = dfg_counts["MES"],
        DFG_AC       = dfg_counts["AC"],
        DFG_NPC      = dfg_counts["NPC"],
        DFG_OPC      = dfg_counts["OPC"],
        DFG_Total    = sum(dfg_counts, na.rm = TRUE),
        DFG_MaxState = TOP_STATES[which.max(dfg_counts)],
        DFG_MaxCount = max(dfg_counts, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    
    cat(sprintf("  DFG summary - MES:%d | AC:%d | NPC:%d | OPC:%d | Total:%d\n\n",
                dfg_counts["MES"], dfg_counts["AC"],
                dfg_counts["NPC"], dfg_counts["OPC"],
                sum(dfg_counts, na.rm = TRUE)))
}

# ===========================================================================
# SECTION 4: Master DFG ranking - all GOIs ranked by total DFG count
# ===========================================================================
master_df <- bind_rows(master_rows) %>%
    arrange(desc(DFG_Total), desc(DFG_MaxCount)) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, everything())

cat("\n#################################################################\n")
cat("MASTER DFG RANKING - state-stratified, all longlist genes\n")
cat(sprintf("Dual-filter thresholds: Spearman >= %.1f AND Kendall >= %.1f\n",
            SP_THRESHOLD, KD_THRESHOLD))
cat("#################################################################\n\n")

cat(sprintf("%-5s %-12s %8s %8s %8s %8s %10s  %s\n",
            "Rank", "GOI", "MES", "AC", "NPC", "OPC",
            "Total", "MaxState"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(master_df))) {
    r <- master_df[i, ]
    cat(sprintf("%-5d %-12s %8s %8s %8s %8s %10s  %s\n",
                r$Rank, r$GOI,
                ifelse(is.na(r$DFG_MES), "SKIP", r$DFG_MES),
                ifelse(is.na(r$DFG_AC),  "SKIP", r$DFG_AC),
                ifelse(is.na(r$DFG_NPC), "SKIP", r$DFG_NPC),
                ifelse(is.na(r$DFG_OPC), "SKIP", r$DFG_OPC),
                r$DFG_Total, r$DFG_MaxState))
}
cat(strrep("-", 70), "\n")

# Top gene per state
cat("\nTop GOI per state:\n")
for (st in TOP_STATES) {
    col <- paste0("DFG_", st)
    top <- master_df %>% arrange(desc(.data[[col]])) %>% slice(1)
    cat(sprintf("  %s: %s (%d DFGs)\n", st, top$GOI, top[[col]]))
}

# Export master ranking
master_file <- paste0(OUTPUT_PREFIX, "_MasterDFG_Ranking.xlsx")
write_xlsx(list(Master_DFG_Ranking = as.data.frame(master_df)),
           path = master_file)
cat(sprintf("\nMaster ranking saved to: %s\n", master_file))

cat("\n=================================================================\n")
cat("Part 4 (State-Stratified) complete.\n")
cat("Per-GOI files: *_Spearman.xlsx | *_Kendall.xlsx | *_DualFilter.xlsx\n")
cat("Master summary:", master_file, "\n")
cat("=================================================================\n")
cat("Ready for Part 5 (CEP-IP inflection point analysis)\n")

# Clean up global parallel cluster
stopCluster(cl_global)
cat("Parallel cluster stopped.\n")


##################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4b: DFG Diverging Bar - FOXM1 as sole MES-viable candidate in GBM
##################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(scales)

# ----------------------------------------------------------
# DATA
# ----------------------------------------------------------
dfg <- data.frame(
    GOI = c("FOXM1","TOP2A","ZWINT","BIRC5","TYMS","MELK","MCM3","PCNA","ATAD2",
            "BCAN","EGFR","OLIG2","PTPRZ1","IL13RA2","EPHA2","CD276","ERBB2",
            "CD70","PROM1","PDGFRA","CD44","MMP9","TNC","CHI3L1","POSTN","FN1",
            "LGALS3","SPARC","MMP2","SOX2","NES","VIM","AURKA","CDC20","UBE2C",
            "MKI67","CCNB1","CDK4","RFC4","HOXA3","EN1","ZIC1","FOXD3","MYC",
            "YAP1","STAT3","SLC16A3","ARHGAP9","ARHGAP30","CLEC7A","MAN2B1",
            "ARPC1B","PLB1","GADD45A","FSTL1","SOD2","AEBP1","CALM3","SIGLEC9",
            "MYO1F","LAPTM5","AGBL2","BLOC1S6","MAP1A","ZSWIM5","PBK","CEP55",
            "MAGEA1","MAGEA3","MAGEA4","MAGEA12","CTAG1B","FBXO39","CTCFL",
            "CASC5","DDX43","ACTL8","OIP5","XAGE3","RGS22","ANXA1","CCNA2",
            "CDC6","DCAF12","E2F7","FABP7","FOSL1","HIF1A","HOPX","IGFBP2",
            "LAMC1","MMP1","NUF2","PRAME","RFX4","SEMA6A","TIMP1","TTK",
            "TWIST1","VEGFA"),
    DFG_MES = c(7, rep(0, 99)),
    DFG_AC  = rep(0, 100),
    DFG_NPC = c(2,0,0,0,3,0,0,0,0,2,1,1,1, rep(0, 87)),
    DFG_OPC = c(4,48,33,26,21,5,4,2,2,0,0,0,0, rep(0, 87))
)

# ----------------------------------------------------------
# COLOURS
# ----------------------------------------------------------
col_mes    <- "#A0445A"
col_nonmes <- "gray82"

# ----------------------------------------------------------
# DERIVED COLUMNS & GENE ORDER
# ----------------------------------------------------------
dfg <- dfg %>%
    mutate(
        total   = DFG_MES + DFG_AC + DFG_NPC + DFG_OPC,
        non_mes = DFG_AC + DFG_NPC + DFG_OPC
    )

gene_order <- c(
    dfg %>% filter(GOI != "FOXM1") %>% arrange(desc(total)) %>% pull(GOI),
    "FOXM1"
)
dfg$GOI <- factor(dfg$GOI, levels = gene_order)

# ----------------------------------------------------------
# AXIS TRUNCATION (y-axis, downward / non-MES direction)
# ----------------------------------------------------------
k1 <- 0.28
k2 <- 0.18
k3 <- 0.28
gap_w <- 0.4

d_at_8  <- 8  * k1
d_at_20 <- d_at_8  + gap_w
d_at_35 <- d_at_20 + 15 * k2
d_at_45 <- d_at_35 + gap_w
d_at_50 <- d_at_45 + 5  * k3

compress_nonmes <- function(x) {
    dplyr::case_when(
        x <=  8 ~ x * k1,
        x <= 20 ~ d_at_8  + gap_w * (x -  8) / 12,
        x <= 35 ~ d_at_20 + (x - 20) * k2,
        x <= 45 ~ d_at_35 + gap_w * (x - 35) / 10,
        TRUE    ~ d_at_45 + (x - 45) * k3
    )
}

dfg <- dfg %>%
    mutate(
        non_mes_disp = -compress_nonmes(non_mes),
        mes_disp     =  DFG_MES
    )

dfg_div <- dfg %>%
    select(GOI, non_mes_disp, mes_disp) %>%
    pivot_longer(cols = c(mes_disp, non_mes_disp),
                 names_to = "direction", values_to = "val")

# ----------------------------------------------------------
# Y-AXIS BREAKS & LABELS
# ----------------------------------------------------------
raw_down  <- c(0, 5, 8, 20, 25, 35, 45, 50)
disp_down <- -compress_nonmes(raw_down)

raw_up    <- c(2, 4, 6, 7)
disp_up   <- as.numeric(raw_up)

all_breaks <- c(disp_down, disp_up)
all_labels <- c(as.character(raw_down), as.character(raw_up))

# ----------------------------------------------------------
# PLOT
# ----------------------------------------------------------
pA <- ggplot(dfg_div,
             aes(x = GOI, y = val,
                 fill = direction, alpha = direction)) +
    
    geom_col(width = 0.7) +
    
    geom_col(
        data = filter(dfg_div, GOI == "FOXM1", direction == "mes_disp"),
        aes(x = GOI, y = val),
        fill      = col_mes,
        color     = "#7A2C40",
        linewidth = 0.55,
        width     = 0.7,
        inherit.aes = FALSE
    ) +
    
    geom_hline(yintercept = 0, color = "gray20", linewidth = 0.5) +
    
    scale_fill_manual(
        values = c(mes_disp = col_mes, non_mes_disp = col_nonmes),
        labels = c(mes_disp = "DFG_MES", non_mes_disp = "Non-MES total"),
        name   = NULL
    ) +
    scale_alpha_manual(
        values = c(mes_disp = 1, non_mes_disp = 0.65),
        guide  = "none"
    ) +
    scale_x_discrete(limits = gene_order) +
    scale_y_continuous(
        breaks = all_breaks,
        labels = all_labels,
        name   = "Non-MES DFG count \u2193  |  \u2191 MES DFG count",
        expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    labs(
        title    = "Diverging bar: MES signal vs non-MES signal",
        subtitle = "Only FOXM1 rises above zero \u2014 every other gene\u2019s signal is entirely non-MES",
        x        = NULL
    ) +
    
    theme_minimal(base_size = 9) +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size  = 3,
            color = ifelse(gene_order == "FOXM1", col_mes, "gray40"),
            face  = ifelse(gene_order == "FOXM1", "bold.italic", "italic")
        ),
        axis.text.y        = element_text(size = 5.5),   # reduced from 7
        axis.title.y       = element_text(size = 6, color = "gray40"),
        axis.ticks.y       = element_line(color = "gray50", linewidth = 0.4),
        axis.ticks.length  = unit(3, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
        legend.position    = "top",
        legend.text        = element_text(size = 6),
        plot.title         = element_text(face = "bold", size = 11),
        plot.subtitle      = element_text(color = "gray45", size = 8)
    )

print(pA)

###############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4c: Gene Panel Box/Jitter Plots & Statistics
# GOIs: FOXM1, BIRC5, TYMS, ZWINT, TOP2A
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#          Adult malignant cells only
#
# Prerequisite: gbm_results and seurat_obj must be in memory from Part 2
# Groupings:
#   (a) Louvain cluster
#   (b) Tumor
#   (c) Dominant cell state (all combinations)
#   (d) Top canonical state (MES / AC / NPC / OPC)
###############################################################################

library(ggplot2)
library(dplyr)
library(viridis)
library(scales)
library(writexl)
library(rstatix)
library(Seurat)

# ===========================================================================
# CONFIG
# ===========================================================================
GENE_PANEL    <- c("FOXM1", "BIRC5", "MKI67")
OUTPUT_PREFIX <- "NeftelGBM_SS2_AdultMalignant_Part4b"

# ===========================================================================
# STEP 1: Verify object and build colour palettes
# ===========================================================================
seurat_p <- gbm_results$seurat_obj

cat(sprintf("SS2 adult malignant object: %d cells, %d clusters\n",
            ncol(seurat_p), nlevels(Idents(seurat_p))))
cat(sprintf("Tumors: %s\n",
            paste(sort(unique(seurat_p$orig.ident)), collapse = ", ")))

# Cluster colours - full plasma spectrum reversed
cluster_order      <- sort(unique(as.numeric(as.character(Idents(seurat_p)))))
n_cl               <- length(cluster_order)
plasma_cols        <- rev(viridis(n_cl, option = "plasma", begin = 0, end = 1))
named_cluster_cols <- setNames(plasma_cols, as.character(cluster_order))

# Tumor colours - turbo palette
tumor_order      <- sort(unique(seurat_p$orig.ident))
tumor_cols       <- viridis(length(tumor_order), option = "turbo",
                            begin = 0, end = 1)
named_tumor_cols <- setNames(tumor_cols, tumor_order)

# Top canonical state colours
state4_cols <- c("MES" = "#C0392B", "AC" = "#F39C12",
                 "NPC" = "#27AE60", "OPC" = "#4A90D9")

# ===========================================================================
# STEP 2: Feature plots on UMAP
# ===========================================================================
cat("\nGenerating feature plots on SS2 UMAP...\n")
gene_colours <- c(scales::alpha("lightgray", 0.85),
                  scales::alpha("lightpink", 0.85),
                  scales::alpha("#FF6666",   0.85),
                  scales::alpha("#BC2727",   0.85),
                  scales::alpha("#660000",   0.85))

for (g in GENE_PANEL) {
    if (!g %in% rownames(seurat_p)) {
        cat("Gene", g, "not found in matrix. Skipping feature plot.\n")
        next
    }
    fp <- FeaturePlot(seurat_p, features = g,
                      min.cutoff = "q10", max.cutoff = "q90",
                      pt.size = 0.2, cols = gene_colours) +
        theme_classic() +
        labs(title = paste("UMAP -", g, "expression"),
             x = "UMAP_1", y = "UMAP_2")
    print(fp)
}

# ===========================================================================
# STEP 3: Expression extractor
# ===========================================================================
get_expr_one <- function(gene, seurat_res, seurat_raw) {
    if (!gene %in% rownames(seurat_raw)) {
        warning(paste("Gene not found, skipping:", gene)); return(NULL)
    }
    cells <- Cells(seurat_res$seurat_obj)
    meta  <- seurat_res$seurat_obj@meta.data[cells, ]
    data.frame(
        gene           = gene,
        cluster        = factor(Idents(seurat_res$seurat_obj),
                                levels = as.character(cluster_order)),
        expression     = as.numeric(
            GetAssayData(seurat_raw, layer = "data")[gene, cells]),
        tumor          = factor(meta$orig.ident,     levels = tumor_order),
        Dominant_State = factor(meta$Dominant_State),
        Top_State      = factor(meta$Top_State,
                                levels = c("MES", "AC", "NPC", "OPC")),
        stringsAsFactors = FALSE
    )
}

cat("\nExtracting expression for gene panel...\n")
panel_expr  <- do.call(rbind,
                       lapply(GENE_PANEL, get_expr_one,
                              seurat_res = gbm_results,
                              seurat_raw = seurat_obj))
genes_found <- unique(panel_expr$gene)
cat("Genes found:", paste(genes_found, collapse = ", "), "\n")
missing_genes <- setdiff(GENE_PANEL, genes_found)
if (length(missing_genes) > 0)
    cat("Genes NOT found:", paste(missing_genes, collapse = ", "), "\n")

# ===========================================================================
# STEP 4: Box + jitter plots - 4 groupings
# ===========================================================================

# ---- (a) By Louvain cluster ----
cat("\nPlotting by Louvain cluster...\n")
for (g in genes_found) {
    df <- panel_expr[panel_expr$gene == g, ]
    p  <- ggplot(df, aes(x = cluster, y = expression, fill = cluster)) +
        geom_jitter(aes(color = cluster), alpha = 0.4, width = 0.3, height = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.2) +
        scale_fill_manual(values  = named_cluster_cols) +
        scale_color_manual(values = named_cluster_cols) +
        theme_minimal() +
        theme(axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none") +
        labs(x     = "Louvain Cluster",
             y     = paste("Expression of", g),
             title = paste(g, "- SS2 adult malignant GBM by Louvain cluster"))
    print(p)
}

# ---- (b) By tumor ----
cat("Plotting by tumor...\n")
for (g in genes_found) {
    df <- panel_expr[panel_expr$gene == g, ]
    p  <- ggplot(df, aes(x = tumor, y = expression, fill = tumor)) +
        geom_jitter(aes(color = tumor), alpha = 0.4, width = 0.3, height = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.2) +
        scale_fill_manual(values  = named_tumor_cols) +
        scale_color_manual(values = named_tumor_cols) +
        theme_minimal() +
        theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "none") +
        labs(x     = "Tumor",
             y     = paste("Expression of", g),
             title = paste(g, "- SS2 adult malignant GBM by Tumor"))
    print(p)
}

# ---- (c) By dominant cell state (all combinations) ----
cat("Plotting by dominant cell state...\n")
for (g in genes_found) {
    df <- panel_expr[panel_expr$gene == g &
                         !is.na(panel_expr$Dominant_State), ]
    p  <- ggplot(df, aes(x = Dominant_State, y = expression,
                         fill = Dominant_State)) +
        geom_jitter(aes(color = Dominant_State), alpha = 0.4,
                    width = 0.3, height = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.2) +
        theme_minimal() +
        theme(axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "none") +
        labs(x     = "Dominant Cell State",
             y     = paste("Expression of", g),
             title = paste(g, "- SS2 adult malignant GBM by Dominant Cell State"))
    print(p)
}

# ---- (d) By top canonical state (MES / AC / NPC / OPC) ----
cat("Plotting by top canonical state...\n")
for (g in genes_found) {
    df <- panel_expr[panel_expr$gene == g &
                         !is.na(panel_expr$Top_State), ]
    p  <- ggplot(df, aes(x = Top_State, y = expression, fill = Top_State)) +
        geom_jitter(aes(color = Top_State), alpha = 0.4, width = 0.3, height = 0) +
        geom_boxplot(outlier.shape = NA, alpha = 0.2) +
        scale_fill_manual(values  = state4_cols) +
        scale_color_manual(values = state4_cols) +
        theme_minimal() +
        theme(axis.text.x     = element_text(angle = 30, vjust = 1, hjust = 1),
              legend.position = "none") +
        labs(x     = "Top Canonical State",
             y     = paste("Expression of", g),
             title = paste(g, "- MES vs AC vs NPC vs OPC"))
    print(p)
}

# ===========================================================================
# STEP 5: Statistics
# (a) KW + Dunn across Louvain clusters
# (b) KW + Dunn across tumors
# (c) KW + Dunn across dominant cell states (all combinations)
# (d) KW + Dunn across top canonical states (MES/AC/NPC/OPC)
# ===========================================================================
cat("\nRunning statistics...\n")

run_stats_one <- function(gene, panel_df) {
    df <- panel_df[panel_df$gene == gene, ]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df$cluster        <- droplevels(df$cluster)
    df$tumor          <- droplevels(df$tumor)
    df$Dominant_State <- droplevels(df$Dominant_State)
    df$Top_State      <- droplevels(df$Top_State)
    
    safe_kw   <- function(f, d) tryCatch(
        as.data.frame(kruskal_test(f, data = d)),
        error = function(e) data.frame(note = e$message))
    safe_dunn <- function(f, d) tryCatch(
        as.data.frame(dunn_test(f, data = d, p.adjust.method = "BH")),
        error = function(e) data.frame(note = e$message))
    desc_stat <- function(grp_col, d) {
        d %>% group_by(.data[[grp_col]]) %>%
            summarise(n      = n(),
                      median = median(expression),
                      Q1     = quantile(expression, 0.25),
                      Q3     = quantile(expression, 0.75),
                      IQR    = IQR(expression),
                      .groups = "drop")
    }
    
    df_dom <- df[!is.na(df$Dominant_State), ]
    df_top <- df[!is.na(df$Top_State), ]
    
    list(
        kw_clusters    = safe_kw(expression ~ cluster,        df),
        dunn_clusters  = safe_dunn(expression ~ cluster,      df),
        stats_clusters = desc_stat("cluster",                 df),
        kw_tumors      = safe_kw(expression ~ tumor,          df),
        dunn_tumors    = safe_dunn(expression ~ tumor,        df),
        stats_tumors   = desc_stat("tumor",                   df),
        kw_domstate    = safe_kw(expression ~ Dominant_State, df_dom),
        dunn_domstate  = safe_dunn(expression ~ Dominant_State, df_dom),
        stats_domstate = desc_stat("Dominant_State",          df_dom),
        kw_topstate    = safe_kw(expression ~ Top_State,      df_top),
        dunn_topstate  = safe_dunn(expression ~ Top_State,    df_top),
        stats_topstate = desc_stat("Top_State",               df_top)
    )
}

all_stats <- lapply(setNames(genes_found, genes_found),
                    run_stats_one, panel_df = panel_expr)

# ===========================================================================
# STEP 6: Export to Excel
# ===========================================================================
sheet_list <- list()
for (g in genes_found) {
    s <- all_stats[[g]]; if (is.null(s)) next
    add <- function(suffix, obj)
        sheet_list[[substr(paste0(g, suffix), 1, 31)]] <<- as.data.frame(obj)
    add("_KW_Clusters",    s$kw_clusters)
    add("_Dunn_Clusters",  s$dunn_clusters)
    add("_Stats_Clusters", s$stats_clusters)
    add("_KW_Tumors",      s$kw_tumors)
    add("_Dunn_Tumors",    s$dunn_tumors)
    add("_Stats_Tumors",   s$stats_tumors)
    add("_KW_DomState",    s$kw_domstate)
    add("_Dunn_DomState",  s$dunn_domstate)
    add("_Stats_DomState", s$stats_domstate)
    add("_KW_TopState",    s$kw_topstate)
    add("_Dunn_TopState",  s$dunn_topstate)
    add("_Stats_TopState", s$stats_topstate)
}

out_xlsx <- paste0(OUTPUT_PREFIX, "_stats.xlsx")
write_xlsx(sheet_list, path = out_xlsx)
cat("Statistics written to:", out_xlsx, "\n")

# ===========================================================================
# STEP 7: Summary
# ===========================================================================
cat("\n===== Part 4b complete =====\n")
cat(sprintf("  Genes plotted : %s\n", paste(genes_found, collapse = ", ")))
cat(sprintf("  Cells         : %d\n", ncol(seurat_p)))
cat(sprintf("  Clusters      : %d\n", nlevels(Idents(seurat_p))))
cat(sprintf("  Output xlsx   : %s\n", out_xlsx))
cat("\nReady for Part 5 (CEP-IP inflection point analysis)\n")


###############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4d: Heatmap - FOXM1, DFG Members & Housekeeping Genes - By Top Cell States
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#          Adult malignant cells only
#
# Prerequisite: gbm_results must be in memory from Part 2
###############################################################################

library(Seurat)
library(ComplexHeatmap)
library(viridis)

# ===========================================================================
# SECTION 1: Define gene sets
# ===========================================================================
GOI             <- "FOXM1"
DFG_MEMBERS     <- c("BIRC5","MKI67", "CENPF", "TOP2A", "PBK", "TROAP", "NUSAP1")
HKG_MEMBERS     <- c("ACTB", "PGK1", "PPIA", "RPL13A", "SDHA")

# DFG weights - proportional to number of cell states each gene was identified in:
#   BIRC5, MKI67  : 3 states (MES, OPC, NPC)  → weight = 3/13
#   CENPF, TOP2A  : 2 states (OPC, NPC)        → weight = 2/13
#   PBK, TROAP,
#   NUSAP1        : 1 state  (MES only)         → weight = 1/13
DFG_WEIGHTS <- c(
    BIRC5   = 3,
    MKI67   = 3,
    CENPF   = 2,
    TOP2A   = 2,
    PBK     = 1,
    TROAP   = 1,
    NUSAP1  = 1
)
DFG_WEIGHTS <- DFG_WEIGHTS / sum(DFG_WEIGHTS)   # normalise → sum = 1

# Column grouping - cell states in this display order
STATE_ORDER <- c("AC", "NPC", "OPC", "MES")

# ===========================================================================
# SECTION 2: Subset Seurat object to all 4 cell states
# ===========================================================================
seurat_hm <- gbm_results$seurat_obj
DefaultAssay(seurat_hm) <- "RNA"

# Use Top_State metadata for grouping (not Louvain clusters)
all_cells   <- Cells(seurat_hm)
state_labels <- seurat_hm$Top_State[all_cells]

# Exclude cells with NA Top_State (43 cells from Part 1)
valid_cells  <- all_cells[!is.na(state_labels)]
state_labels <- state_labels[!is.na(state_labels)]

cat(sprintf("Total cells after NA removal: %d\n", length(valid_cells)))
cat("Cells per state:\n")
print(table(state_labels))

# ===========================================================================
# SECTION 3: Compute weighted DFG and simple-average HKG expression
# ===========================================================================

# Weighted average helper - weights must be a named numeric vector
compute_weighted_avg_expr <- function(members, weights, seurat_obj, cells) {
    members_found <- members[members %in% rownames(seurat_obj)]
    missing       <- setdiff(members, members_found)
    if (length(missing) > 0)
        cat("  WARNING - genes not found:", paste(missing, collapse = ", "), "\n")
    
    # Re-normalise weights in case any gene was missing from the dataset
    w <- weights[members_found]
    w <- w / sum(w)
    
    cat("  Effective weights used:\n")
    print(round(w, 6))
    
    mat <- GetAssayData(seurat_obj, layer = "data")[members_found, cells,
                                                    drop = FALSE]
    # Weighted sum across genes per cell
    as.numeric(Matrix::t(mat) %*% w)
}

# Simple average helper (unchanged for HKG)
compute_avg_expr <- function(members, seurat_obj, cells) {
    members_found <- members[members %in% rownames(seurat_obj)]
    missing       <- setdiff(members, members_found)
    if (length(missing) > 0)
        cat("  WARNING - HKG not found:", paste(missing, collapse = ", "), "\n")
    mat <- GetAssayData(seurat_obj, layer = "data")[members_found, cells,
                                                    drop = FALSE]
    Matrix::colMeans(mat)
}

cat("\nComputing weighted DFG composite expression...\n")
cat("  Members:", paste(DFG_MEMBERS, collapse = ", "), "\n")
dfg_avg <- compute_weighted_avg_expr(DFG_MEMBERS, DFG_WEIGHTS, seurat_hm, valid_cells)

cat("Computing HKG averaged expression...\n")
cat("  Members:", paste(HKG_MEMBERS, collapse = ", "), "\n")
hkg_avg <- compute_avg_expr(HKG_MEMBERS, seurat_hm, valid_cells)

# ===========================================================================
# SECTION 4: Assemble expression matrix
# ===========================================================================
extract_indiv <- function(genes, seurat_obj, cells) {
    found   <- genes[genes %in% rownames(seurat_obj)]
    missing <- setdiff(genes, found)
    if (length(missing) > 0)
        cat("  WARNING - genes not found:", paste(missing, collapse = ", "), "\n")
    as.matrix(GetAssayData(seurat_obj, layer = "data")[found, cells,
                                                       drop = FALSE])
}

mat_goi <- extract_indiv(GOI,         seurat_hm, valid_cells)
mat_dfg <- extract_indiv(DFG_MEMBERS, seurat_hm, valid_cells)
mat_hkg <- extract_indiv(HKG_MEMBERS, seurat_hm, valid_cells)

expression_matrix <- rbind(
    mat_goi,
    DFG = as.numeric(dfg_avg),
    mat_dfg,
    HKG = as.numeric(hkg_avg),
    mat_hkg
)
rownames(expression_matrix) <- c(
    GOI,
    "DFG",
    rownames(mat_dfg),
    "HKG",
    rownames(mat_hkg)
)

cat(sprintf("\nExpression matrix: %d rows x %d cells\n",
            nrow(expression_matrix), ncol(expression_matrix)))
cat("Row order:", paste(rownames(expression_matrix), collapse = ", "), "\n")

# ===========================================================================
# SECTION 5: Log2 transform and z-score scaling (cap -1 to 1.5)
# ===========================================================================
expression_matrix_log <- log2(expression_matrix + 1)

scale_to_zscore <- function(x) {
    z <- scale(x)
    pmin(pmax(z, -1), 2)
}
expression_matrix_scaled <- t(apply(expression_matrix_log, 1,
                                    scale_to_zscore))

# ===========================================================================
# SECTION 6: Order cells by state (AC | NPC | OPC | MES) then
#            by FOXM1 expression ascending within each state
# ===========================================================================
state_factor  <- factor(state_labels, levels = STATE_ORDER)
foxm1_expr    <- GetAssayData(seurat_hm, layer = "data")[GOI, valid_cells]

cell_order               <- order(state_factor, foxm1_expr)
expression_matrix_log    <- expression_matrix_log[,  cell_order, drop = FALSE]
expression_matrix_scaled <- t(apply(expression_matrix_log, 1,
                                    scale_to_zscore))
state_factor             <- state_factor[cell_order]
foxm1_expr               <- foxm1_expr[cell_order]

# ===========================================================================
# SECTION 7: Column annotation (state labels with cell counts)
# ===========================================================================
# Four distinct colours - one per canonical state, matching Part 4b palette
state_colours <- c(
    "AC"  = "#F39C12",
    "NPC" = "#27AE60",
    "OPC" = "#4A90D9",
    "MES" = "#C0392B"
)

cell_counts      <- table(state_factor)
cell_percentages <- round(prop.table(cell_counts) * 100, 1)
group_labels     <- paste0(names(cell_counts), "\n",
                           cell_counts, " cells\n(",
                           cell_percentages, "%)")

column_annotation <- HeatmapAnnotation(
    State = anno_block(
        gp        = gpar(fill = state_colours[levels(state_factor)]),
        labels    = group_labels,
        labels_gp = gpar(col = "white", fontsize = 9, fontface = "bold")),
    show_legend          = FALSE,
    show_annotation_name = FALSE
)

# ===========================================================================
# SECTION 8: Row annotation - distinguish GOI / DFG / HKG blocks
# ===========================================================================
# row_block_labels must map every row to one of the colour keys
row_block_labels <- c(
    "GOI",                                  # FOXM1
    "DFG_avg",                              # DFG composite (immediately below GOI)
    rep("DFG", length(rownames(mat_dfg))),  # individual DFG genes
    "HKG_avg",                              # HKG composite
    rep("HKG", length(rownames(mat_hkg)))   # individual HKG genes
)

row_colours <- c(
    "GOI"     = "#8E44AD",   # purple  - FOXM1
    "DFG"     = "#2980B9",   # blue    - individual DFG members
    "DFG_avg" = "#1A5276",   # dark blue - DFG composite
    "HKG_avg" = "#566573",   # dark grey - HKG composite
    "HKG"     = "#7F8C8D"    # grey    - individual HKG members
)

row_annotation <- rowAnnotation(
    Block = anno_simple(
        row_block_labels,
        col       = row_colours,
        width     = unit(3, "mm")),
    show_annotation_name = FALSE,
    show_legend          = FALSE
)

# ===========================================================================
# SECTION 9: Draw heatmap
# ===========================================================================
cat("\nDrawing heatmap...\n")

heatmap_obj <- Heatmap(
    expression_matrix_scaled,
    name               = "Z-score",
    column_title       = sprintf("%s & DFG Members - SS2 Adult Malignant GBM by Cell State",
                                 GOI),
    column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    row_title          = NULL,
    show_row_names     = TRUE,
    show_column_names  = FALSE,
    cluster_rows       = FALSE,
    cluster_columns    = FALSE,
    top_annotation     = column_annotation,
    left_annotation    = row_annotation,
    col                = viridis(100, option = "mako"),
    row_names_gp       = gpar(fontsize = 10),
    row_names_side     = "right",
    column_split       = state_factor,
    row_gap            = unit(1,   "mm"),
    column_gap         = unit(1.5, "mm"),
    border             = TRUE,
    use_raster         = TRUE,
    raster_quality     = 2,
    heatmap_legend_param = list(
        title          = "Z-score\n(capped -1 to 2)",
        title_gp       = gpar(fontsize = 9),
        labels_gp      = gpar(fontsize = 8),
        legend_height  = unit(3, "cm")
    )
)

draw(heatmap_obj)

# ===========================================================================
# SECTION 10: Summary
# ===========================================================================
cat("\n=== Part 4d complete ===\n")
cat(sprintf("  GOI          : %s\n", GOI))
cat(sprintf("  DFG members  : %s\n", paste(rownames(mat_dfg), collapse = ", ")))
cat(sprintf("  DFG weights  : %s\n",
            paste(sprintf("%s=%.4f", names(DFG_WEIGHTS), DFG_WEIGHTS), collapse = ", ")))
cat(sprintf("  HKG members  : %s\n", paste(HKG_MEMBERS, collapse = ", ")))
cat(sprintf("  Cell states  : %s\n", paste(STATE_ORDER, collapse = " | ")))
cat(sprintf("  Cells plotted: %d\n", ncol(expression_matrix_scaled)))
cat(sprintf("  Rows plotted : %s\n",
            paste(rownames(expression_matrix_scaled), collapse = ", ")))


###############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4d (FOXM1): Heatmap - FOXM1, DFG Members & Housekeeping Genes - By Patient ID
###############################################################################
library(Seurat)
library(ComplexHeatmap)
library(viridis)

# ===========================================================================
# SECTION 1: Define gene sets and patients
# ===========================================================================
GOI         <- "FOXM1"
DFG_MEMBERS <- c("BIRC5", "MKI67", "CENPF", "TOP2A", "PBK", "TROAP", "NUSAP1")
HKG_MEMBERS <- c("ACTB", "PGK1", "PPIA", "RPL13A", "SDHA")

# DFG weights - proportional to number of cell states each gene was identified in:
#   BIRC5, MKI67: 3 states (MES, OPC, NPC)  → weight = 3/13
#   CENPF, TOP2A : 2 states (OPC, NPC)        → weight = 2/13
#   PBK, TROAP, NUSAP1: 1 state  (MES only)   → weight = 1/13
DFG_WEIGHTS <- c(
    BIRC5   = 3,
    MKI67   = 3,
    CENPF   = 2,
    TOP2A   = 2,
    PBK     = 1,
    TROAP   = 1,
    NUSAP1  = 1
)
DFG_WEIGHTS <- DFG_WEIGHTS / sum(DFG_WEIGHTS)   # normalise → sum = 1

# Patient display order (left → right in heatmap)
PATIENT_ORDER <- c("MGH152", "MGH104", "MGH100", "MGH110", "MGH136",
                   "MGH66",  "MGH124", "MGH128", "MGH143", "MGH121",
                   "MGH125", "MGH122", "MGH102", "MGH113", "MGH106",
                   "MGH129", "MGH101", "MGH115", "MGH105", "MGH151")

patients <- PATIENT_ORDER

# ===========================================================================
# SECTION 2: Subset Seurat object to selected patients
# ===========================================================================
seurat_hm <- gbm_results$seurat_obj
DefaultAssay(seurat_hm) <- "RNA"

all_cells      <- Cells(seurat_hm)
patient_labels <- seurat_hm$orig.ident[all_cells]

# Keep only cells belonging to the selected patients
valid_cells    <- all_cells[patient_labels %in% PATIENT_ORDER]
patient_labels <- patient_labels[patient_labels %in% PATIENT_ORDER]

cat(sprintf("Total cells after patient filter: %d\n", length(valid_cells)))
cat("Cells per patient:\n")
print(table(patient_labels))

# ===========================================================================
# SECTION 3: Compute composite averaged expression
# ===========================================================================

# Weighted average helper - for DFG
compute_weighted_avg_expr <- function(members, weights, seurat_obj, cells) {
    members_found <- members[members %in% rownames(seurat_obj)]
    missing       <- setdiff(members, members_found)
    if (length(missing) > 0)
        cat("  WARNING - genes not found:", paste(missing, collapse = ", "), "\n")
    w <- weights[members_found]
    w <- w / sum(w)
    cat("  Effective weights used:\n")
    print(round(w, 6))
    mat <- GetAssayData(seurat_obj, layer = "data")[members_found, cells,
                                                    drop = FALSE]
    as.numeric(Matrix::t(mat) %*% w)
}

# Simple average helper - for HKG
compute_avg_expr <- function(members, seurat_obj, cells) {
    members_found <- members[members %in% rownames(seurat_obj)]
    missing       <- setdiff(members, members_found)
    if (length(missing) > 0)
        cat("  WARNING - genes not found:", paste(missing, collapse = ", "), "\n")
    mat <- GetAssayData(seurat_obj, layer = "data")[members_found, cells,
                                                    drop = FALSE]
    Matrix::colMeans(mat)
}

cat("\nComputing weighted DFG composite expression...\n")
cat("  Members:", paste(DFG_MEMBERS, collapse = ", "), "\n")
dfg_avg <- compute_weighted_avg_expr(DFG_MEMBERS, DFG_WEIGHTS, seurat_hm, valid_cells)

cat("Computing HKG averaged expression...\n")
cat("  Members:", paste(HKG_MEMBERS, collapse = ", "), "\n")
hkg_avg <- compute_avg_expr(HKG_MEMBERS, seurat_hm, valid_cells)

# ===========================================================================
# SECTION 4: Assemble expression matrix
# ===========================================================================
extract_indiv <- function(genes, seurat_obj, cells) {
    found   <- genes[genes %in% rownames(seurat_obj)]
    missing <- setdiff(genes, found)
    if (length(missing) > 0)
        cat("  WARNING - genes not found:", paste(missing, collapse = ", "), "\n")
    as.matrix(GetAssayData(seurat_obj, layer = "data")[found, cells,
                                                       drop = FALSE])
}

mat_goi <- extract_indiv(GOI,         seurat_hm, valid_cells)
mat_dfg <- extract_indiv(DFG_MEMBERS, seurat_hm, valid_cells)
mat_hkg <- extract_indiv(HKG_MEMBERS, seurat_hm, valid_cells)

expression_matrix <- rbind(
    mat_goi,
    DFG = as.numeric(dfg_avg),
    mat_dfg,
    HKG = as.numeric(hkg_avg),
    mat_hkg
)
rownames(expression_matrix) <- c(
    GOI,
    "DFG",
    rownames(mat_dfg),
    "HKG",
    rownames(mat_hkg)
)

cat(sprintf("\nExpression matrix: %d rows x %d cells\n",
            nrow(expression_matrix), ncol(expression_matrix)))
cat("Row order:", paste(rownames(expression_matrix), collapse = ", "), "\n")

# ===========================================================================
# SECTION 5: Log2 transform and z-score scaling (cap -1 to 2)
# ===========================================================================
expression_matrix_log <- log2(expression_matrix + 1)

scale_to_zscore <- function(x) {
    z <- scale(x)
    pmin(pmax(z, -1), 2)
}
expression_matrix_scaled <- t(apply(expression_matrix_log, 1,
                                    scale_to_zscore))

# ===========================================================================
# SECTION 6: Compute FOXM1-positive % per patient, reorder patients
#            (highest → lowest % positive), then order cells by patient
#            then FOXM1 expression ascending within each patient
# ===========================================================================
FOXM1_expr_all <- GetAssayData(seurat_hm, layer = "data")[GOI, valid_cells]

# FOXM1-positive = log2(count+1) > 0, i.e. raw count > 0
foxm1_pos_pct <- sapply(PATIENT_ORDER, function(p) {
    idx     <- which(patient_labels == p)
    pos_pct <- sum(FOXM1_expr_all[idx] > 0) / length(idx) * 100
    round(pos_pct, 1)
})

# Reorder patients from highest to lowest FOXM1-positive %
PATIENT_ORDER <- names(sort(foxm1_pos_pct, decreasing = TRUE))
foxm1_pos_pct <- foxm1_pos_pct[PATIENT_ORDER]

cat("\nFOXM1-positive % per patient (high → low):\n")
print(foxm1_pos_pct)

patient_factor <- factor(patient_labels, levels = PATIENT_ORDER)
FOXM1_expr     <- FOXM1_expr_all

cell_order               <- order(patient_factor, FOXM1_expr)
expression_matrix_log    <- expression_matrix_log[,  cell_order, drop = FALSE]
expression_matrix_scaled <- t(apply(expression_matrix_log, 1,
                                    scale_to_zscore))
patient_factor           <- patient_factor[cell_order]
FOXM1_expr               <- FOXM1_expr[cell_order]

# ===========================================================================
# SECTION 7: Column annotation - continuous red → pink → violet → navy
#            palette; per-patient font sizes for cramped boxes
# ===========================================================================

# Extended palette now routes through pink and violet between blue and navy,
# giving a fuller spectral arc:
# dark red → red → orange → yellow → green → teal → blue → violet → pink → navy
hot_to_cold_palette <- colorRampPalette(c(
    "#7B0000",   # darkest red
    "#C0392B",   # red
    "#E67E22",   # orange
    "#F1C40F",   # yellow
    "#27AE60",   # green
    "#1ABC9C",   # teal
    "#2980B9",   # blue
    "#1A2F6E",   # darkest navy
    "#8E44AD",   # violet
    "#E91E8C"    # pink
))(length(PATIENT_ORDER))

# First patient (highest FOXM1+%) gets darkest red; last gets darkest navy
patient_colours <- setNames(hot_to_cold_palette, PATIENT_ORDER)

cell_counts <- table(patient_factor)

# Extract single dominant (modal) Top_State per patient
top_state_meta <- seurat_hm@meta.data[valid_cells, "Top_State"]
names(top_state_meta) <- valid_cells

dominant_state_per_patient <- sapply(names(cell_counts), function(p) {
    idx    <- which(patient_labels == p)
    states <- top_state_meta[valid_cells[idx]]
    states <- states[!is.na(states) & states != ""]
    if (length(states) == 0) return("NA")
    # Modal state - single label only (AC, MES, NPC, or OPC)
    names(sort(table(states), decreasing = TRUE))[1]
})

# Label: patient ID, cell count, FOXM1-positive %, modal dominant state
group_labels <- paste0(
    names(cell_counts), "\n",
    cell_counts, " cells\n",
    foxm1_pos_pct[names(cell_counts)], "%\n",
    dominant_state_per_patient[names(cell_counts)]
)

# Per-patient font sizes - reduce for patients with narrow boxes
# (MGH128, MGH106, MGH101, MGH115, MGH151 identified as cramped)
CRAMPED_PATIENTS <- c("MGH128", "MGH106", "MGH101", "MGH115", "MGH151")
label_fontsizes  <- ifelse(PATIENT_ORDER %in% CRAMPED_PATIENTS, 5.0, 6.5)
label_fontsizes[PATIENT_ORDER == "MGH151"] <- 4.2   # adjust value

column_annotation <- HeatmapAnnotation(
    Patient = anno_block(
        gp        = gpar(fill = patient_colours[levels(patient_factor)]),
        labels    = group_labels,
        labels_gp = gpar(col      = "white",
                         fontsize = label_fontsizes,
                         fontface = "bold")),
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    annotation_height    = unit(1.35, "cm")   # ~15% reduction from default 2 cm
)

# ===========================================================================
# SECTION 8: Row annotation - distinguish GOI / DFG / HKG blocks
# ===========================================================================
row_block_labels <- c(
    "GOI",                                  # FOXM1
    "DFG_avg",                              # DFG composite
    rep("DFG", length(rownames(mat_dfg))),  # individual DFG members
    "HKG_avg",                              # HKG composite
    rep("HKG", length(rownames(mat_hkg)))   # individual HKG members
)

row_colours <- c(
    "GOI"     = "#8E44AD",   # purple    - FOXM1
    "DFG"     = "#2980B9",   # blue      - individual DFG members
    "DFG_avg" = "#1A5276",   # dark blue - DFG composite
    "HKG_avg" = "#566573",   # dark grey - HKG composite
    "HKG"     = "#7F8C8D"    # grey      - individual HKG members
)

row_annotation <- rowAnnotation(
    Block = anno_simple(
        row_block_labels,
        col       = row_colours,
        width     = unit(3, "mm")),
    show_annotation_name = FALSE,
    show_legend          = FALSE
)

# ===========================================================================
# SECTION 9: Draw heatmap
# ===========================================================================
cat("\nDrawing heatmap...\n")

heatmap_obj <- Heatmap(
    expression_matrix_scaled,
    name               = "Z-score",
    column_title       = sprintf("%s & DFG Members (Weighted) - SS2 Adult Malignant GBM by Patient",
                                 GOI),
    column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    row_title          = NULL,
    show_row_names     = TRUE,
    show_column_names  = FALSE,
    cluster_rows       = FALSE,
    cluster_columns    = FALSE,
    top_annotation     = column_annotation,
    left_annotation    = row_annotation,
    col                = viridis(100, option = "mako"),
    row_names_gp       = gpar(fontsize = 10),
    row_names_side     = "right",
    column_split       = patient_factor,
    row_gap            = unit(1,   "mm"),
    column_gap         = unit(0.6, "mm"),   # reduced by 60% from original 1.5 mm
    border             = TRUE,
    use_raster         = TRUE,
    raster_quality     = 2,
    heatmap_legend_param = list(
        title          = "Z-score\n(capped -1 to 2)",
        title_gp       = gpar(fontsize = 9),
        labels_gp      = gpar(fontsize = 8),
        legend_height  = unit(3, "cm")
    )
)

draw(heatmap_obj)

# ===========================================================================
# SECTION 10: Summary
# ===========================================================================
cat("\n=== Part 4d (FOXM1) complete ===\n")
cat(sprintf("  GOI          : %s\n", GOI))
cat(sprintf("  DFG members  : %s\n", paste(rownames(mat_dfg), collapse = ", ")))
cat(sprintf("  DFG weights  : %s\n",
            paste(sprintf("%s=%.4f", names(DFG_WEIGHTS), DFG_WEIGHTS), collapse = ", ")))
cat(sprintf("  HKG members  : %s\n", paste(HKG_MEMBERS, collapse = ", ")))
cat(sprintf("  Patients     : %s\n", paste(PATIENT_ORDER, collapse = " | ")))
cat(sprintf("  Cells plotted: %d\n", ncol(expression_matrix_scaled)))
cat(sprintf("  Rows plotted : %s\n",
            paste(rownames(expression_matrix_scaled), collapse = ", ")))




##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 4e: Within-Positive Monotonicity Diagnostics
#          FOXM1 vs DFG composite & individual DFG genes
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#          Adult malignant cells only
#
# Loads:  NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds  (from Part 5)
#
# Aims:
#   Determine whether the high Spearman/Kendall correlations between FOXM1
#   and DFG genes reflect genuine within-positive monotonicity or are driven
#   by zero-concordance (cells jointly zero dominating the rank statistic).
#
# Prerequisite: gbm_results in memory AND
#               NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds on disk
##############################################################################
library(mgcv)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)

set.seed(123)

# ===========================================================================
# SECTION 1: Load data
# ===========================================================================
cat("Loading Part 5 models...\n")
rds_data     <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
patient_data <- rds_data$patient_data
patients     <- rds_data$patients

PLOT_PATIENTS <- c("MGH152", "MGH104", "MGH100", "MGH110", "MGH136",
              "MGH66",  "MGH124", "MGH128", "MGH143", "MGH121",
              "MGH125", "MGH122", "MGH102", "MGH113", "MGH106",
              "MGH129", "MGH101", "MGH115", "MGH105", "MGH151")
              
DFG_GENES   <- c("BIRC5", "MKI67", "CENPF", "TOP2A", "PBK", "TROAP", "NUSAP1")
DFG_WEIGHTS <- c(BIRC5=3, MKI67=3, CENPF=2, TOP2A=2, PBK=1, TROAP=1, NUSAP1=1)
DFG_WEIGHTS <- DFG_WEIGHTS / sum(DFG_WEIGHTS)

cat("Patients loaded:", length(patients), "\n")
cat("Plot patients  :", paste(PLOT_PATIENTS, collapse = ", "), "\n\n")

# ===========================================================================
# SECTION 2: Shared theme
# ===========================================================================
theme_diag <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1,
                                            hjust = 0.5),
            plot.subtitle    = element_text(size = base_size - 1, hjust = 0.5,
                                            color = "gray40"),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            panel.border     = element_rect(color = "black", fill = NA,
                                            linewidth = 0.4),
            panel.grid.major = element_line(color = "gray92"),
            panel.grid.minor = element_line(color = "gray96"),
            plot.background  = element_rect(fill = "white", color = NA)
        )
}

# ===========================================================================
# SECTION 3: Build per-patient expression data frames
# ===========================================================================
# For each patient extract: foxm1, dfg_comp, individual DFGs, phase (if any)
# Uses the full mat (all cells) stored in patient_data from Part 5.
# ===========================================================================

seurat_p5 <- gbm_results$seurat_obj
has_phase  <- "Phase" %in% colnames(seurat_p5@meta.data)
cat("Cell cycle Phase metadata available:", has_phase, "\n\n")

build_patient_df <- function(patient) {
    dd <- patient_data[[patient]]
    if (is.null(dd)) return(NULL)

    mat   <- dd$mat          # genes x cells, ALL cells
    cells <- colnames(mat)

    foxm1 <- log2(mat["FOXM1", ] + 1)

    # Weighted DFG composite (all cells)
    dfg_present <- DFG_GENES[DFG_GENES %in% rownames(mat)]
    w           <- DFG_WEIGHTS[dfg_present]
    w           <- w / sum(w)
    dfg_comp    <- log2(as.numeric(t(mat[dfg_present, , drop = FALSE]) %*% w) + 1)

    # Individual DFG genes
    indiv <- lapply(dfg_present, function(g) log2(mat[g, ] + 1))
    names(indiv) <- dfg_present

    # Cell cycle phase (if available)
    phase <- if (has_phase) seurat_p5$Phase[cells] else
                             rep(NA_character_, length(cells))

    df <- data.frame(
        Patient  = patient,
        Cell     = cells,
        FOXM1    = foxm1,
        DFG_comp = dfg_comp,
        Phase    = phase,
        stringsAsFactors = FALSE
    )
    for (g in dfg_present) df[[g]] <- indiv[[g]]
    df
}

cat("Building per-patient expression data frames...\n")
patient_dfs       <- lapply(patients, build_patient_df)
names(patient_dfs) <- patients
cat("Done.\n\n")

# ===========================================================================
# SECTION 4: D1 - Within-positive Spearman rho table
# ===========================================================================
cat("========== D1: Within-positive dual-filter (Spearman + Kendall) ==========\n\n")

spearman_rows <- list()

for (patient in patients) {
    df <- patient_dfs[[patient]]
    if (is.null(df)) next

    is_pos <- df$FOXM1 > 0
    n_all  <- nrow(df)
    n_pos  <- sum(is_pos)

    genes_to_test <- c("DFG_comp", intersect(DFG_GENES, colnames(df)))

    for (gene in genes_to_test) {
        y_all <- df[[gene]]
        x_all <- df$FOXM1
        x_pos <- x_all[is_pos]
        y_pos <- y_all[is_pos]

        # --- Spearman rho ---
        rho_all <- if (length(unique(x_all)) > 2 && length(unique(y_all)) > 2)
            suppressWarnings(cor(x_all, y_all, method = "spearman"))
        else NA_real_

        rho_pos <- if (n_pos > 5 && length(unique(x_pos)) > 2 &&
                       length(unique(y_pos)) > 2)
            suppressWarnings(cor(x_pos, y_pos, method = "spearman"))
        else NA_real_

        p_rho_pos <- if (!is.na(rho_pos) && n_pos > 5)
            suppressWarnings(
                cor.test(x_pos, y_pos, method = "spearman")$p.value)
        else NA_real_

        rho_delta <- if (!is.na(rho_all) && !is.na(rho_pos))
            rho_pos - rho_all else NA_real_

        # --- Kendall tau ---
        tau_all <- if (length(unique(x_all)) > 2 && length(unique(y_all)) > 2)
            suppressWarnings(cor(x_all, y_all, method = "kendall"))
        else NA_real_

        tau_pos <- if (n_pos > 5 && length(unique(x_pos)) > 2 &&
                       length(unique(y_pos)) > 2)
            suppressWarnings(cor(x_pos, y_pos, method = "kendall"))
        else NA_real_

        p_tau_pos <- if (!is.na(tau_pos) && n_pos > 5)
            suppressWarnings(
                cor.test(x_pos, y_pos, method = "kendall")$p.value)
        else NA_real_

        tau_delta <- if (!is.na(tau_all) && !is.na(tau_pos))
            tau_pos - tau_all else NA_real_

        # --- Dual-filter interpretation ---
        # Mirrors original DFG selection: Spearman >= 0.50 AND Kendall >= 0.40
        # within FOXM1-positive cells only.
        # n_pos < 50 is evaluated first regardless of correlation values -
        # too few positive cells for reliable within-positive estimation.
        pass_rho <- !is.na(rho_pos) && rho_pos >= 0.50
        pass_tau <- !is.na(tau_pos) && tau_pos >= 0.40

        interpretation <-
            if (is.na(rho_pos) || is.na(tau_pos)) "insufficient data" else
            if (n_pos < 50)                        "insufficient FOXM1+ cells (n<50)" else
            if (pass_rho && pass_tau)              "dual-filter pass: genuine within-positive monotonicity" else
            if (pass_rho || pass_tau)              "partial pass: one filter met" else
            if (rho_pos >= 0.20 || tau_pos >= 0.15) "weak within-positive signal" else
                                                   "zero-concordance dominated"

        spearman_rows[[length(spearman_rows) + 1]] <- data.frame(
            Patient        = patient,
            Gene           = gene,
            N_all_cells    = n_all,
            N_pos_cells    = n_pos,
            Pct_pos        = round(100 * n_pos / n_all, 1),
            Rho_all_cells  = round(rho_all, 4),
            Rho_pos_only   = round(rho_pos, 4),
            Rho_delta      = round(rho_delta, 4),
            P_rho_pos      = p_rho_pos,
            Tau_all_cells  = round(tau_all, 4),
            Tau_pos_only   = round(tau_pos, 4),
            Tau_delta      = round(tau_delta, 4),
            P_tau_pos      = p_tau_pos,
            Pass_Rho       = pass_rho,
            Pass_Tau       = pass_tau,
            Dual_Filter_Pass = pass_rho && pass_tau,
            Interpretation = interpretation,
            stringsAsFactors = FALSE
        )

        cat(sprintf("  %-10s | %-10s | n_pos=%3d | rho_pos=%6.3f (pass=%s) | tau_pos=%6.3f (pass=%s) | %s\n",
                    patient, gene, n_pos,
                    ifelse(is.na(rho_pos), NA, rho_pos),
                    ifelse(pass_rho, "Y", "N"),
                    ifelse(is.na(tau_pos), NA, tau_pos),
                    ifelse(pass_tau, "Y", "N"),
                    interpretation))
    }
}

sheet_spearman <- do.call(rbind, spearman_rows)
rownames(sheet_spearman) <- NULL

# Median summary across all patients (Spearman + Kendall)
cat("\n--- Median within-positive Spearman rho & Kendall tau across all 20 patients ---\n")
cat(sprintf("  %-10s  %-26s  %-26s  %s\n",
            "Gene", "rho_all | rho_pos | delta",
            "tau_all | tau_pos | delta", "Dual-pass"))
for (gene in c("DFG_comp", DFG_GENES)) {
    sub <- sheet_spearman[sheet_spearman$Gene == gene, ]
    if (nrow(sub) == 0) next
    n_dual <- sum(sub$Dual_Filter_Pass, na.rm = TRUE)
    cat(sprintf("  %-10s : rho %+.3f | %+.3f | %+.3f   tau %+.3f | %+.3f | %+.3f   dual-pass=%d/%d\n",
                gene,
                median(sub$Rho_all_cells, na.rm = TRUE),
                median(sub$Rho_pos_only,  na.rm = TRUE),
                median(sub$Rho_delta,     na.rm = TRUE),
                median(sub$Tau_all_cells, na.rm = TRUE),
                median(sub$Tau_pos_only,  na.rm = TRUE),
                median(sub$Tau_delta,     na.rm = TRUE),
                n_dual, nrow(sub)))
}

# ===========================================================================
# SECTION 5: D2 - Pseudobulk decile plot
# ===========================================================================
cat("\n\n========== D2: Pseudobulk decile plots ==========\n")

decile_rows <- list()

for (patient in PLOT_PATIENTS) {
    df <- patient_dfs[[patient]]
    if (is.null(df)) next

    df_pos <- df[df$FOXM1 > 0, ]
    if (nrow(df_pos) < 10) {
        cat("  Skipping", patient, "- fewer than 10 positive cells\n"); next
    }

    n_bins    <- min(10L, floor(nrow(df_pos) / 3))
    probs_seq <- seq(0, 1, length.out = n_bins + 1)
    breaks    <- unique(quantile(df_pos$FOXM1, probs = probs_seq, na.rm = TRUE))
    n_bins    <- length(breaks) - 1L

    df_pos$decile <- as.integer(cut(df_pos$FOXM1, breaks = breaks,
                                    include.lowest = TRUE,
                                    labels = seq_len(n_bins)))

    genes_avail <- intersect(DFG_GENES, colnames(df_pos))

    decile_summary <- df_pos %>%
        group_by(decile) %>%
        summarise(
            n_cells    = n(),
            mean_foxm1 = mean(FOXM1,    na.rm = TRUE),
            mean_dfg   = mean(DFG_comp, na.rm = TRUE),
            across(all_of(genes_avail),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "mean_{.col}"),
            .groups = "drop"
        ) %>%
        mutate(Patient = patient)

    decile_rows[[patient]] <- as.data.frame(decile_summary)

    # Long format for individual genes
    gene_cols  <- paste0("mean_", genes_avail)
    long_genes <- decile_summary %>%
        select(decile, mean_foxm1, all_of(gene_cols)) %>%
        pivot_longer(cols      = all_of(gene_cols),
                     names_to  = "Gene",
                     values_to = "Mean_Expression") %>%
        mutate(Gene = sub("^mean_", "", Gene))

    gene_pal <- setNames(
        viridis::viridis(length(genes_avail), option = "turbo"),
        genes_avail)

    p_dec <- ggplot() +
        geom_line(data = long_genes,
                  aes(x = decile, y = Mean_Expression,
                      color = Gene, group = Gene),
                  linewidth = 0.7, linetype = "dashed", alpha = 0.8) +
        geom_point(data = long_genes,
                   aes(x = decile, y = Mean_Expression, color = Gene),
                   size = 2.0, alpha = 0.8) +
        geom_line(data = decile_summary,
                  aes(x = decile, y = mean_dfg),
                  color = "#3A3A3A", linewidth = 1.4) +
        geom_point(data = decile_summary,
                   aes(x = decile, y = mean_dfg),
                   color = "#3A3A3A", size = 3.5, shape = 18) +
        scale_color_manual(values = gene_pal, name = "DFG gene") +
        scale_x_continuous(breaks = seq_len(n_bins)) +
        labs(title    = paste0(patient, " - Pseudobulk Decile Plot"),
             subtitle = paste0("FOXM1-positive cells binned into ", n_bins,
                               " deciles by FOXM1 expression\n",
                               "Dark gray diamond = DFG composite | dashed = individual genes (turbo palette)"),
             x = paste0("FOXM1 Decile  (1 = lowest, ", n_bins, " = highest)"),
             y = "Mean log2(Expression + 1)") +
        theme_diag() +
        theme(legend.position = "right")
    print(p_dec)
    cat("  Plotted D2:", patient, "| deciles:", n_bins,
        "| n_pos:", nrow(df_pos), "\n")
}

sheet_decile <- if (length(decile_rows) > 0)
    do.call(rbind, decile_rows) else data.frame()
rownames(sheet_decile) <- NULL

# ===========================================================================
# SECTION 7: D3 - Cell cycle phase overlay (within positives)
# ===========================================================================
cat("\n\n========== D3: Cell cycle phase overlay ==========\n")

if (!has_phase) {
    cat("  Phase metadata not found in seurat_p5@meta.data.\n")
    cat("  To enable D5, run Seurat::CellCycleScoring() before Part 4e.\n")
    cat("  Example:\n")
    cat("    cc_genes <- readRDS(system.file('extdata','cc.genes.updated.2019.rds',\n")
    cat("                        package='Seurat'))\n")
    cat("    seurat_p5 <- CellCycleScoring(seurat_p5,\n")
    cat("                    s.features   = cc_genes$s.genes,\n")
    cat("                    g2m.features = cc_genes$g2m.genes)\n")
    cat("  Then re-run Part 4e.\n")
} else {
    phase_colours <- c("G1" = "#95A5A6", "S" = "#3498DB", "G2M" = "#E74C3C")

    for (patient in PLOT_PATIENTS) {
        df <- patient_dfs[[patient]]
        if (is.null(df)) next

        df_pos <- df[df$FOXM1 > 0 & !is.na(df$Phase), ]
        if (nrow(df_pos) < 10) {
            cat("  Skipping", patient, "- fewer than 10 positive cells\n"); next
        }

        df_pos$rank_foxm1 <- rank(df_pos$FOXM1,    ties.method = "average")
        df_pos$rank_dfg   <- rank(df_pos$DFG_comp,  ties.method = "average")

        phase_rho <- df_pos %>%
            group_by(Phase) %>%
            summarise(
                n   = n(),
                rho = ifelse(n() > 5,
                             round(suppressWarnings(
                                 cor(FOXM1, DFG_comp,
                                     method = "spearman")), 3),
                             NA_real_),
                .groups = "drop")

        phase_label <- phase_rho %>%
            filter(!is.na(rho)) %>%
            mutate(lbl = sprintf("%s: rho=%.3f (n=%d)", Phase, rho, n)) %>%
            pull(lbl) %>%
            paste(collapse = " | ")

        p_phase <- ggplot(df_pos,
                          aes(x = rank_foxm1, y = rank_dfg, color = Phase)) +
            geom_point(alpha = 0.55, size = 1.8) +
            geom_smooth(aes(group = Phase), method = "lm", se = FALSE,
                        linewidth = 0.8, linetype = "dashed") +
            scale_color_manual(values = phase_colours,
                               name   = "Cell cycle phase") +
            annotate("text", x = -Inf, y = Inf,
                     label  = phase_label,
                     hjust  = -0.05, vjust = 1.5,
                     size   = 2.8, color = "gray30") +
            labs(title    = paste0(patient, " - Cell Cycle Phase Overlay"),
                 subtitle = paste0("Rank-rank scatter coloured by cell cycle phase\n",
                                   "Tests whether phase-offset confounds ",
                                   "within-positive monotonicity"),
                 x = "FOXM1 Rank (within positive cells)",
                 y = "DFG Composite Rank (within positive cells)") +
            theme_diag() +
            theme(legend.position = "right")
        print(p_phase)
        cat("  Plotted D3:", patient, "\n")
        cat("  Phase breakdown:", phase_label, "\n")
    }
}

# ===========================================================================
# SECTION 8: Export to Excel
# ===========================================================================
cat("\n\n========== Exporting ==========\n")

export_list <- list(
    DualFilter_Summary = sheet_spearman,
    Decile_Values    = sheet_decile
)

out_file <- "NeftelGBM_SS2_FOXM1_DFG_Monotonicity.xlsx"
tryCatch({
    write_xlsx(export_list, path = out_file)
    cat("Saved:", out_file, "\n")
    cat("  DualFilter_Summary:", nrow(sheet_spearman), "rows (Spearman + Kendall per patient x gene)\n")
    cat("  Decile_Values    :", nrow(sheet_decile),   "rows\n")
}, error = function(e) cat("Export error:", conditionMessage(e), "\n"))

# ===========================================================================
# SECTION 9: Interpretation guide
# ===========================================================================
cat("\n========== Interpretation Guide ==========\n")
cat("\nD1 Dual-filter (Spearman + Kendall) - key diagnostic:\n")
cat("  rho_pos / tau_pos : within FOXM1-positive cells ONLY (the critical numbers)\n")
cat("  rho_delta / tau_delta : drop after removing zeros (negative = inflation)\n")
cat("  Dual-filter mirrors original DFG selection criteria (evaluated in order):\n")
cat("    Insufficient: n_pos < 50 -> insufficient FOXM1+ cells for reliable estimation\n")
cat("    Pass        : rho_pos >= 0.50 AND tau_pos >= 0.40 -> genuine within-positive monotonicity\n")
cat("    Partial     : one filter met (rho>=0.50 OR tau>=0.40) -> asymmetric signal\n")
cat("    Weak        : rho_pos >= 0.20 OR tau_pos >= 0.15 -> weak but real signal\n")
cat("    Dominated   : both below weak thresholds -> zero-concordance dominated\n")
cat("\nD2 Decile plot - what to look for:\n")
cat("  Rising DFG composite across deciles = population-level monotonicity\n")
cat("  Flat or non-monotonic = no graded relationship within positives\n")
cat("  Individual gene divergence across deciles = phase-offset confounding\n")
cat("\nD3 Phase overlay - what to look for:\n")
cat("  G2M cells should show highest rho if FOXM1-DFG is phase-dependent\n")
cat("  Similar rho across phases = phase is not the main confounder\n")

cat("\n========== Part 4e complete ==========\n")
cat("All plots displayed in RStudio Plots pane (no files saved).\n")
cat("Tabular results: NeftelGBM_SS2_FOXM1_DFG_Monotonicity.xlsx\n")




###############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 5: GAM-REML-PRSS Analysis - FOXM1 vs DFG composite & HKG composite
# FOXM1-positive adult malignant cells only
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#
# x-axis : FOXM1 expression
# y-axis : DFG composite (frequency-weighted mean of BIRC5, MKI67, CENPF,
#                         TOP2A, PBK, TROAP, NUSAP1)
#        | HKG composite (simple mean of ACTB, PGK1, PPIA, RPL13A, SDHA)
#
# Patients : 20 patients in specified order (MGH152 first ... MGH151 last)
#
# Zero-inflation mitigation:
#   GAM fitting restricted to FOXM1-positive cells only
#
# Prerequisite: gbm_results must be in memory from Part 2, OR load workspace:
#   load("NeftelGBM_SS2_AdultMalignant_Seurat_processed.RData")
###############################################################################
library(mgcv)
library(Seurat)
library(dplyr)
library(writexl)
set.seed(123)

# ===========================================================================
# SECTION 1: Per-patient data preparation
# ===========================================================================

patients <- c("MGH152", "MGH104", "MGH100", "MGH110", "MGH136",
              "MGH66",  "MGH124", "MGH128", "MGH143", "MGH121",
              "MGH125", "MGH122", "MGH102", "MGH113", "MGH106",
              "MGH129", "MGH101", "MGH115", "MGH105", "MGH151")

DFG_GENES <- c("BIRC5", "MKI67", "CENPF", "TOP2A", "PBK", "TROAP", "NUSAP1")
HKG_GENES <- c("ACTB", "PGK1", "PPIA", "RPL13A", "SDHA")

# --- DFG frequency weights (proportional to cell-state appearances) ---
# BIRC5, MKI67: 3 states (MES, OPC, NPC) → 3/13
# CENPF, TOP2A: 2 states (OPC, NPC)       → 2/13
# PBK, TROAP, NUSAP1: 1 state (MES)       → 1/13
DFG_WEIGHTS <- c(BIRC5=3, MKI67=3, CENPF=2, TOP2A=2, PBK=1, TROAP=1, NUSAP1=1)
DFG_WEIGHTS <- DFG_WEIGHTS / sum(DFG_WEIGHTS)   # normalise → sum = 1

# Minimum FOXM1-positive cells to consider a fit adequately powered
MIN_POS_CELLS <- 50L

# Within-positive Spearman rho_pos for DFG composite per patient
# Source: Part 4e diagnostics (D1 Spearman table, Gene = DFG_comp)
# Used for relationship classifier in build_gam_row.
# Patients with n_pos < MIN_POS_CELLS are flagged Underpowered regardless of rho.
RHO_POS_DFG <- c(
    MGH152 = 0.6362, MGH104 = 0.3921, MGH100 = 0.5569, MGH110 = 0.6282,
    MGH136 = 0.3718, MGH66  = 0.6593, MGH124 = 0.4462, MGH128 = 0.6011,
    MGH143 = 0.4553, MGH121 = 0.5364, MGH125 = 0.2627, MGH122 = 0.6422,
    MGH102 = 0.5390, MGH113 = 0.5597, MGH106 = 0.8164, MGH129 = 0.6670,
    MGH101 = 0.0688, MGH115 = 0.5515, MGH105 = 0.5320, MGH151 = 0.3187
)

seurat_p5 <- gbm_results$seurat_obj
DefaultAssay(seurat_p5) <- "RNA"

all_cells <- Cells(seurat_p5)

cat("=================================================================\n")
cat("GAM-REML-PRSS: FOXM1 vs DFG & HKG  -  SS2 Adult Malignant GBM\n")
cat("Patients         :", length(patients), "\n")
cat("DFG genes        :", paste(DFG_GENES, collapse = ", "), "\n")
cat("DFG weights      :", paste(sprintf("%s=%.4f", names(DFG_WEIGHTS),
                                        DFG_WEIGHTS), collapse = ", "), "\n")
cat("HKG genes        :", paste(HKG_GENES, collapse = ", "), "\n")
cat("Total cells      :", length(all_cells), "\n")
cat("Zero-inflation  : FOXM1-zero cells excluded from GAM fitting entirely\n")
cat("Min pos cells   : patients with <", MIN_POS_CELLS, "FOXM1-pos cells flagged Underpowered\n")
cat("=================================================================\n\n")

cat("Extracting full normalised matrix from SS2 object...\n")
full_mat <- GetAssayData(seurat_p5, layer = "data")
cat(sprintf("  Matrix: %d genes x %d cells\n", nrow(full_mat), ncol(full_mat)))

for (g in c("FOXM1", DFG_GENES, HKG_GENES)) {
    if (!g %in% rownames(full_mat))
        warning(sprintf("Gene not found in matrix: %s", g))
}

cat("\nBuilding per-patient data...\n")

patient_data <- lapply(patients, function(patient) {
    
    cells <- all_cells[seurat_p5$orig.ident == patient]
    
    if (length(cells) < 10) {
        cat("  Patient", patient, ": only", length(cells),
            "cells -- skipping (min 10 required)\n")
        return(NULL)
    }
    
    mat_pt <- as.matrix(full_mat[, cells, drop = FALSE])
    
    dfg_present <- DFG_GENES[DFG_GENES %in% rownames(mat_pt)]
    hkg_present <- HKG_GENES[HKG_GENES %in% rownames(mat_pt)]
    
    foxm1_all <- log2(mat_pt["FOXM1", ] + 1)
    
    # --- Restrict to FOXM1-positive cells only for GAM fitting ---
    # "Positive" = log2(normalised + 1) > 0, consistent across all patients.
    # Zero-FOXM1 cells are excluded entirely from the data frames passed to
    # fit_gam_prss. This eliminates zero-inflation distortion of the GAM curve
    # and IP placement without requiring any weighting scheme.
    # Full mat_pt (all cells) is retained for LOO composite calculation.
    is_pos    <- foxm1_all > 0
    p_pos     <- mean(is_pos)
    mat_pos   <- mat_pt[, is_pos, drop = FALSE]
    foxm1     <- foxm1_all[is_pos]
    
    if (sum(is_pos) < 10) {
        cat("  Patient", patient, ": only", sum(is_pos),
            "FOXM1-positive cells -- skipping (min 10 required)\n")
        return(NULL)
    }
    
    # --- Weighted DFG composite (positive cells only) ---
    dfg_avg <- if (length(dfg_present) > 0) {
        w <- DFG_WEIGHTS[dfg_present]
        w <- w / sum(w)   # re-normalise in case any gene absent
        log2(as.numeric(t(mat_pos[dfg_present, , drop = FALSE]) %*% w) + 1)
    } else {
        cat("  WARNING - no DFG genes found for patient", patient, "\n")
        rep(NA_real_, sum(is_pos))
    }
    
    # --- Simple-mean HKG composite (positive cells only) ---
    hkg_avg <- if (length(hkg_present) > 0)
        log2(colMeans(mat_pos[hkg_present, , drop = FALSE]) + 1)
    else {
        cat("  WARNING - no HKG genes found for patient", patient, "\n")
        rep(NA_real_, sum(is_pos))
    }
    
    # --- All fitted cells receive equal weight = 1 (no zeros present) ---
    cell_weights  <- rep(1.0, sum(is_pos))
    n_pos         <- sum(is_pos)
    underpowered  <- n_pos < MIN_POS_CELLS
    
    cat("  Patient", patient,
        "| cells (total):", length(cells),
        "| FOXM1-pos (fitted):", n_pos,
        sprintf("(%.1f%%)", 100 * p_pos),
        if (underpowered) "*** UNDERPOWERED ***" else "",
        "| DFG genes:", length(dfg_present),
        "| HKG genes:", length(hkg_present), "\n")
    
    list(
        dfg          = data.frame(FOXM1 = foxm1, Expression = dfg_avg),
        hkg          = data.frame(FOXM1 = foxm1, Expression = hkg_avg),
        cell_weights = cell_weights,   # all 1.0 - kept for API consistency
        p_pos        = p_pos,
        n_pos        = n_pos,
        underpowered = underpowered,
        foxm1        = foxm1,          # positive cells only
        foxm1_all    = foxm1_all,      # all cells - retained for reference
        mat          = mat_pt,         # all cells
        mat_pos      = mat_pos,        # positive cells only - used in LOO
        dfg_present  = dfg_present,
        hkg_present  = hkg_present
    )
})
names(patient_data) <- patients
cat("Patients with sufficient cells:",
    sum(!sapply(patient_data, is.null)), "/", length(patients), "\n")

# ===========================================================================
# SECTION 2: Helper functions  
# ===========================================================================

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
    pt <- coef(model)[!grepl("s\\(", names(coef(model)))]
    st <- coef(model)[ grepl("s\\(", names(coef(model)))]
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

# ===========================================================================
# SECTION 3: GAM fitting with two-phase PRSS optimisation
#
# Only change from original: gam() call now passes weights = cell_weights
# (dynamic zero-weighting). All PRSS / REML / convergence logic unchanged.
# ===========================================================================
fit_gam_prss <- function(data, patient_id, gene_set_label,
                         cell_weights, num_iterations = 100) {
    set.seed(123)
    unique_x <- length(unique(data$FOXM1))
    max_k    <- min(10, unique_x - 1)
    cat("  [", patient_id, "|", gene_set_label, "]",
        "cells:", nrow(data), "| unique FOXM1:", unique_x,
        "| max k:", max_k, "\n")
    
    best_prss <- Inf; best_model <- NULL; best_k <- NULL
    best_iter <- 0;   prss_values <- numeric(num_iterations)
    model_params <- list(); reml_iters_list <- list()
    early_stop_counter <- 0; early_stop_triggered <- FALSE
    
    run_one <- function(i, k) {
        set.seed(123 + i)
        # --- only change vs original: weights argument added ---
        gam_model <- gam(Expression ~ FOXM1 + s(FOXM1, bs = "tp", k = k),
                         data    = data,
                         weights = cell_weights,
                         method  = "REML",
                         select  = TRUE,
                         gamma   = 1.5)
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
    
    k_range <- if (best_k == max_k) c(max_k - 1, max_k) else
        if (best_k == 3)    c(3, 4) else
            c(best_k - 1, best_k, best_k + 1)
    cat("    best k after phase 1:", best_k,
        "| PRSS:", round(best_prss, 4), "\n")
    
    # ---- Phase 2: adaptive refinement (iterations 9-100) ----
    for (i in 9:num_iterations) {
        if (early_stop_counter >= 20) {
            cat("    early stopping at iter", i, "\n")
            early_stop_triggered <- TRUE; break
        }
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
         patient        = patient_id,
         gene_set       = gene_set_label,
         data           = data,
         cell_weights   = cell_weights)   # stored for downstream use
}

# ===========================================================================
# SECTION 4: Run 40 fits (20 patients x 2 gene sets)
# ===========================================================================
cat("\n========== Running GAM fits: 20 patients x 2 gene sets ==========\n")

all_results <- list()

for (patient in patients) {
    dd <- patient_data[[patient]]
    if (is.null(dd)) {
        cat("Skipping patient", patient, "(insufficient cells)\n")
        next
    }
    cat("\n--- Patient:", patient,
        sprintf("| p_pos = %.3f | FOXM1-pos cells only ---\n",
                dd$p_pos))
    
    cat(" Fitting DFG (FOXM1 composite)...\n")
    key_dfg <- paste0(patient, "|DFG")
    all_results[[key_dfg]] <- tryCatch(
        fit_gam_prss(dd$dfg, patient, "DFG", dd$cell_weights),
        error = function(e) {
            cat("  FAILED:", conditionMessage(e), "\n"); NULL })
    
    cat(" Fitting HKG...\n")
    key_hkg <- paste0(patient, "|HKG")
    all_results[[key_hkg]] <- tryCatch(
        fit_gam_prss(dd$hkg, patient, "HKG", dd$cell_weights),
        error = function(e) {
            cat("  FAILED:", conditionMessage(e), "\n"); NULL })
}

n_success <- sum(!sapply(all_results, is.null))
cat("\nCompleted:", n_success, "/ 40 fits\n")

# ===========================================================================
# SECTION 5: Sheet-building helpers  
# ===========================================================================

build_convergence_row <- function(result) {
    cv <- extract_convergence_details(result$best_model)
    data.frame(
        Patient                 = result$patient,
        Group                   = "Neftel_GBM_SS2",
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

build_prss_rows <- function(result) {
    n     <- result$actual_iters
    iters <- seq_len(n)
    get_p <- function(i, key) {
        x <- result$model_params[[i]]
        if (is.null(x)) NA else x[[key]]
    }
    data.frame(
        Patient             = result$patient,
        Group               = "Neftel_GBM_SS2",
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
        Patient        = result$patient,
        Group          = "Neftel_GBM_SS2",
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

build_gam_row <- function(result) {
    m  <- result$best_model
    pc <- calculate_prss(m, result$data)
    sm <- summary(m)
    
    p_val <- sm$s.table[4]
    
    smooth_names  <- grepl("s\\(", names(coef(m)))
    parametric_df <- sum(!smooth_names)
    total_edf     <- sum(sm$edf)
    
    smooth_coefs      <- coef(m)[smooth_names]
    nonzero_basis     <- sum(abs(smooth_coefs) > 1e-8)
    total_basis       <- length(smooth_coefs)
    smooth_coef_count <- length(smooth_coefs)
    
    p_pos        <- if (!is.null(result$cell_weights))
        mean(result$cell_weights == 1.0) else NA_real_
    n_pos        <- nrow(result$data)
    underpowered <- isTRUE(n_pos < MIN_POS_CELLS)
    
    # Relationship classifier - based on within-positive rho_pos (Part 4e D1)
    # Only applied to DFG gene set; HKG classifier set to NA
    rho_pos <- if (result$gene_set == "DFG")
        RHO_POS_DFG[result$patient] else NA_real_
    
    relationship_type <- if (result$gene_set != "DFG") NA_character_ else
        if (underpowered)           "Insufficient FOXM1+ cells (n<50)" else
            if (is.na(rho_pos))         "No rho_pos data" else
                if (rho_pos >= 0.50)        "Graded" else
                    if (rho_pos >= 0.20)        "Weak" else
                        "Binary"
    
    data.frame(
        Patient                    = result$patient,
        Group                      = "Neftel_GBM_SS2",
        Gene_Set                   = result$gene_set,
        FOXM1_pos_proportion       = p_pos,
        N_cells_fitted             = n_pos,
        Underpowered               = underpowered,
        Rho_pos_DFG                = rho_pos,
        Relationship_Type          = relationship_type,
        PRSS                       = pc$PRSS,
        RSS                        = pc$RSS,
        Lambda_REML                = pc$lambda,
        f_double_prime_int         = pc$f_double_prime_integral,
        k                          = result$best_k,
        Model_deviance             = deviance(m),
        Null_deviance              = deviance(gam(Expression ~ 1,
                                                  data    = result$data,
                                                  weights = result$cell_weights)),
        DE                         = sm$dev.expl,
        Equation                   = extract_equation(m),
        Smooth_EDF                 = sm$edf[length(sm$edf)],
        Smooth_p_value             = p_val,
        Smooth_p_BH_FDR            = NA_real_,
        Total_EDF                  = total_edf,
        Parametric_DF              = parametric_df,
        Is_Significantly_Nonlinear = NA,
        Nonzero_Basis_Functions    = nonzero_basis,
        Total_Basis_Functions      = total_basis,
        Smooth_Coef_Count          = smooth_coef_count,
        stringsAsFactors = FALSE)
}

build_obs_pred_rows <- function(result) {
    data.frame(
        Patient         = result$patient,
        Gene_Set        = result$gene_set,
        FOXM1           = result$data$FOXM1,
        Observed_value  = result$data$Expression,
        Predicted_value = as.numeric(predict(result$best_model,
                                             newdata = result$data)),
        Cell_Weight     = result$cell_weights,
        stringsAsFactors = FALSE)
}

build_basis_rows <- function(result) {
    m           <- result$best_model
    X           <- predict(m, type = "lpmatrix")
    smooth_cols <- grep("s\\(FOXM1\\)", colnames(X))
    bm          <- X[, smooth_cols, drop = FALSE]
    colnames(bm) <- paste0("phi_", seq_len(ncol(bm)))
    data.frame(Patient  = result$patient,
               Gene_Set = result$gene_set,
               FOXM1    = result$data$FOXM1,
               as.data.frame(bm), stringsAsFactors = FALSE)
}

# ===========================================================================
# SECTION 6: Collect all results into sheets 
# ===========================================================================
valid_results <- Filter(Negate(is.null), all_results)

sheet_conv      <- do.call(rbind, lapply(valid_results, build_convergence_row))
sheet_prss      <- do.call(rbind, lapply(valid_results, build_prss_rows))
sheet_reml_best <- do.call(rbind, lapply(valid_results, build_reml_best_rows))

dfg_results <- Filter(function(r) r$gene_set == "DFG", valid_results)
hkg_results <- Filter(function(r) r$gene_set == "HKG", valid_results)

sheet_perf_dfg <- do.call(rbind, lapply(dfg_results, build_gam_row))
sheet_perf_hkg <- do.call(rbind, lapply(hkg_results, build_gam_row))

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

sheet_obs_pred <- do.call(rbind, lapply(valid_results, build_obs_pred_rows))

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

# ===========================================================================
# SECTION 6b: Leave-one-out gene contribution analysis  
# ===========================================================================
fit_gam_loo <- function(data, best_k, cell_weights, goi_col = "FOXM1") {
    set.seed(123)
    k_use <- max(3L, min(best_k, length(unique(data[[goi_col]])) - 1L))
    tryCatch({
        fml <- as.formula(paste0("Expression ~ ", goi_col,
                                 " + s(", goi_col, ", bs='tp', k=k_use)"))
        m <- gam(fml, data = data, weights = cell_weights,
                 method = "REML", select = TRUE, gamma = 1.5)
        summary(m)$dev.expl
    }, error = function(e) NA_real_)
}

build_gene_contribution_sheet <- function(patient, dd, full_results) {
    rows <- list()
    for (gs in c("DFG", "HKG")) {
        key      <- paste0(patient, "|", gs)
        full_res <- full_results[[key]]
        if (is.null(full_res) || is.null(full_res$best_model)) next
        de_full <- summary(full_res$best_model)$dev.expl
        best_k  <- full_res$best_k
        genes   <- if (gs == "DFG") dd$dfg_present else dd$hkg_present
        mat     <- dd$mat_pos      # positive cells only - consistent with fitted GAM
        foxm1   <- dd$foxm1        # positive cells only
        cw      <- dd$cell_weights # all 1.0
        
        if (length(genes) < 2) {
            rows[[length(rows) + 1]] <- data.frame(
                Patient          = patient,
                Group            = "Neftel_GBM_SS2",
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
            
            # Weighted composite for LOO DFG; simple mean for LOO HKG
            if (gs == "DFG") {
                w <- DFG_WEIGHTS[genes_in]
                w <- w / sum(w)
                loo_avg <- log2(as.numeric(t(mat[genes_in, , drop = FALSE]) %*% w) + 1)
            } else {
                loo_avg <- log2(colMeans(mat[genes_in, , drop = FALSE]) + 1)
            }
            
            loo_data <- data.frame(FOXM1 = foxm1, Expression = loo_avg)
            de_loo   <- fit_gam_loo(loo_data, best_k, cw, goi_col = "FOXM1")
            delta    <- de_full - de_loo
            pct_contrib <- if (!is.na(de_full) && de_full > 0)
                100 * delta / de_full else NA_real_
            cat("    LOO [", patient, "|", gs, "] remove", gene_out,
                "| DE_full:", round(de_full, 4),
                "| DE_loo:", round(de_loo, 4),
                "| delta:", round(delta, 4), "\n")
            rows[[length(rows) + 1]] <- data.frame(
                Patient          = patient,
                Group            = "Neftel_GBM_SS2",
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
cat("Fitting 20 patients x (DFG 7 genes + HKG 5 genes) = up to 240 LOO GAMs\n\n")

loo_list <- lapply(patients, function(patient) {
    dd <- patient_data[[patient]]
    if (is.null(dd)) return(NULL)
    cat("  LOO patient:", patient, "\n")
    tryCatch(
        build_gene_contribution_sheet(patient, dd, all_results),
        error = function(e) {
            cat("  LOO FAILED for patient", patient, ":", conditionMessage(e), "\n")
            NULL
        })
})
sheet_loo <- do.call(rbind, Filter(Negate(is.null), loo_list))
cat("\nLOO sheet rows:", nrow(sheet_loo), "\n")

# ===========================================================================
# SECTION 7: Export  
# ===========================================================================
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

write_xlsx(export_list, path = "NeftelGBM_SS2_FOXM1_DFG_GAM_Results.xlsx")
cat("\nSaved: NeftelGBM_SS2_FOXM1_DFG_GAM_Results.xlsx\n")
cat("Sheets:\n")
cat("  REML_Convergence          ", nrow(sheet_conv),      "rows\n")
cat("  All_Models_PRSS           ", nrow(sheet_prss),      "rows\n")
cat("  REML_Iterations_BestModel ", nrow(sheet_reml_best), "rows\n")
cat("  GAM_Performance_DFG       ", nrow(sheet_perf_dfg),
    "rows | columns: p_pos, N_fitted, Underpowered, Rho_pos_DFG, Relationship_Type,",
    "PRSS, RSS, Lambda, f'', k, DE, EDF, p, BH-FDR, Nonlinear\n")
cat("  GAM_Performance_HKG       ", nrow(sheet_perf_hkg),  "rows (same columns as DFG)\n")
cat("  Observed_Predicted        ", nrow(sheet_obs_pred),  "rows (includes Cell_Weight col)\n")
cat("  Basis_Functions           ", nrow(sheet_basis),     "rows\n")
cat("  Gene_Contribution_LOO     ", nrow(sheet_loo),
    "rows (20 patients x DFG 7 genes + HKG 5 genes)\n")

saveRDS(list(all_results   = all_results,
             patient_data  = patient_data,
             patients      = patients,
             sheet_loo     = sheet_loo),
        "NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
cat("Models saved to: NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds\n")

# Relationship classifier and underpowered summary (DFG only)
if (!is.null(sheet_perf_dfg) && nrow(sheet_perf_dfg) > 0) {
    cat("\n=== DFG Relationship Type Summary (based on within-positive rho_pos) ===\n")
    cat("  Classifier: n<50 = Insufficient | Graded = rho_pos>=0.50 | Weak = 0.20-0.49 | Binary = <0.20\n")
    rt_counts <- table(sheet_perf_dfg$Relationship_Type)
    for (rt in c("Graded", "Weak", "Binary", "Insufficient FOXM1+ cells (n<50)", "No rho_pos data"))
        cat(sprintf("  %-18s : %d patients\n", rt,
                    if (rt %in% names(rt_counts)) rt_counts[[rt]] else 0L))
    cat("  Per-patient:\n")
    for (i in seq_len(nrow(sheet_perf_dfg)))
        cat(sprintf("    %-10s  n_pos=%3d  rho_pos=%s  DE=%.3f  -> %s%s\n",
                    sheet_perf_dfg$Patient[i],
                    sheet_perf_dfg$N_cells_fitted[i],
                    ifelse(is.na(sheet_perf_dfg$Rho_pos_DFG[i]), "  NA ",
                           sprintf("%.3f", sheet_perf_dfg$Rho_pos_DFG[i])),
                    sheet_perf_dfg$DE[i],
                    sheet_perf_dfg$Relationship_Type[i],
                    ifelse(isTRUE(sheet_perf_dfg$Underpowered[i]),
                           " *** UNDERPOWERED ***", "")))
}
cat("\n=== GBM FOXM1 vs DFG/HKG GAM-REML-PRSS Part 5 complete ===\n")



##################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 5b: Plots and Console Summary Statistics
#          FOXM1 vs DFG composite & HKG composite
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
# FOXM1-positive adult malignant cells only
#
# Loads:  NeftelGBM_SS2_FOXM1_DFG_GAM_models.rds  (from Part 5)
# Output: all tabular results are in
#         NeftelGBM_SS2_FOXM1_DFG_GAM_Results.xlsx (Part 5).
##################################################################
library(mgcv)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(viridis)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Load Part 5 models
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data     <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results  <- rds_data$all_results
patient_data <- rds_data$patient_data
patients     <- rds_data$patients
sheet_loo    <- rds_data$sheet_loo

valid_results <- Filter(function(x) !is.null(x) && !is.null(x$best_model), all_results)
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
#
# PRSS plot: first 28 iterations shown (convergence is always complete
#   by iteration 28 given the early-stop at 20 no-improvement rule).
#   Best iteration highlighted with red diamond + dashed vertical line.
#
# REML plot: iterations of the best PRSS model only.
#   Final converged point highlighted with red diamond.
#
# Layout: 4 fits per page (2 rows × n cols), PRSS row on top.
# ------------------------------------------------------------------
cat("========== Building PRSS + REML plots ==========\n")

PRSS_PLOT_MAX_ITER <- 28L   # show only first 28 PRSS iterations

# ---- 3a. PRSS plot builder ----
make_prss_plot <- function(result, compact = FALSE) {
    bs <- if (compact) 7 else 9
    
    # Step 1: hard-cap at PRSS_PLOT_MAX_ITER
    n_cap  <- min(PRSS_PLOT_MAX_ITER, result$actual_iters)
    iters  <- seq_len(n_cap)
    prss_v <- result$prss_values[iters]
    
    df <- data.frame(
        Iteration = iters,
        PRSS      = prss_v,
        Is_Best   = iters == result$best_iter
    )
    df <- df[!is.na(df$PRSS), ]
    
    # Step 2: trim near-zero PRSS tail
    prss_max    <- max(df$PRSS, na.rm = TRUE)
    zero_thresh <- max(prss_max * 0.01, 1e-6)
    first_zero  <- which(df$PRSS < zero_thresh)[1]
    
    if (!is.na(first_zero) && first_zero > 1) {
        df <- df[seq_len(first_zero - 1), ]
    }
    
    best_row  <- df[df$Is_Best, ]
    beyond_28 <- result$best_iter > PRSS_PLOT_MAX_ITER
    
    if (nrow(best_row) == 0 && nrow(df) > 0) {
        best_row  <- df[which.min(df$PRSS), ]
        beyond_28 <- TRUE
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
        labs(title    = paste0(result$patient, " | ", result$gene_set, " - PRSS"),
             subtitle = paste0("Best iter: ", result$best_iter,
                               "  k = ", result$best_k,
                               "  PRSS = ", round(result$best_prss, 3),
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
                   labs(title = paste0(result$patient, " | ", result$gene_set,
                                       " - REML")) +
                   theme_void())
    
    best_ri <- all_ri[all_ri$is_best_model == TRUE, ]
    if (nrow(best_ri) == 0)
        best_ri <- all_ri[all_ri$prss_iteration == result$best_iter, ]
    if (nrow(best_ri) == 0) best_ri <- all_ri
    
    best_ri   <- best_ri[order(best_ri$iteration), ]
    final_row <- best_ri[nrow(best_ri), ]
    
    ggplot(best_ri, aes(x = iteration, y = score)) +
        geom_line(linewidth = 0.6, color = "gray30") +
        geom_point(size = 2, color = "gray20") +
        geom_point(data = final_row,
                   aes(x = iteration, y = score),
                   color = "firebrick", size = 3.5, shape = 18,
                   inherit.aes = FALSE) +
        labs(title    = paste0(result$patient, " | ", result$gene_set, " - REML"),
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
    
    n_cols     <- length(page_keys)
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
# SECTION 4: GAM curve plots - one per patient, DFG + HKG overlaid
# ------------------------------------------------------------------
cat("\n========== GAM curve plots ==========\n")

# DFG composite: BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1
# HKG composite: ACTB, PGK1, PPIA, RPL13A, SDHA
gene_set_colors <- c("DFG" = "#4B0082", "HKG" = "#4F75DE")
gene_set_fills  <- c("DFG" = "#4B0082", "HKG" = "#4F75DE")

make_gam_curve_plot <- function(patient, vres) {
    plot_data <- data.frame()
    for (gs in c("DFG", "HKG")) {
        res <- vres[[paste0(patient, "|", gs)]]
        if (is.null(res) || is.null(res$best_model)) next
        m    <- res$best_model
        data <- res$data
        # Observed points - x column is FOXM1
        pts  <- data.frame(
            FOXM1      = data$FOXM1,
            Expression = data$Expression,
            Gene_Set   = gs,
            se.fit     = NA_real_,
            type       = "data",
            stringsAsFactors = FALSE)
        set.seed(123)
        x_seq <- seq(min(data$FOXM1), max(data$FOXM1), length.out = 1000)
        pred  <- predict(m, newdata = data.frame(FOXM1 = x_seq), se.fit = TRUE)
        fit_df <- data.frame(
            FOXM1      = x_seq,
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
                            label = paste("No data for patient:", patient)) +
                   theme_void())
    
    de_labels <- sapply(c("DFG", "HKG"), function(gs) {
        res <- vres[[paste0(patient, "|", gs)]]
        if (is.null(res) || is.null(res$best_model)) return(NA_character_)
        gs_label <- if (gs == "DFG") "DFG composite" else "HKG"
        sprintf("%s DE=%.1f%%", gs_label, summary(res$best_model)$dev.expl * 100)
    })
    de_labels <- de_labels[!is.na(de_labels)]
    
    set.seed(123)
    ggplot(plot_data, aes(x = FOXM1, y = Expression, color = Gene_Set)) +
        geom_point(data = subset(plot_data, type == "data"),
                   position = position_jitter(width = 0.01, height = 0, seed = 123),
                   alpha = 0.25, size = 0.8) +
        geom_line(data  = subset(plot_data, type == "fit"), linewidth = 1.2) +
        geom_ribbon(data = subset(plot_data, type == "fit"),
                    aes(ymin = Expression - 1.96 * se.fit,
                        ymax = Expression + 1.96 * se.fit,
                        fill = Gene_Set),
                    alpha = 0.20, color = NA) +
        scale_color_manual(values = gene_set_colors,
                           labels = c("DFG" = "DFG composite", "HKG" = "HKG"),
                           name   = "Gene Set") +
        scale_fill_manual( values = gene_set_fills,
                           labels = c("DFG" = "DFG composite", "HKG" = "HKG"),
                           name   = "Gene Set") +
        labs(title    = paste("GAM: Gene Set Expression vs FOXM1 -", patient),
             subtitle = paste(de_labels, collapse = "  |  "),
             x        = "FOXM1 Expression (log2)",
             y        = "Gene Set Mean Expression (log2)") +
        scale_x_continuous(labels = number_format(accuracy = 0.1)) +
        scale_y_continuous(labels = number_format(accuracy = 0.1)) +
        theme_gam() +
        theme(legend.position = "right")
}

for (patient in patients) {
    has_result <- any(sapply(c("DFG", "HKG"), function(gs) {
        res <- valid_results[[paste0(patient, "|", gs)]]
        !is.null(res) && !is.null(res$best_model)
    }))
    if (!has_result) { cat("  Skipping", patient, "(no valid fits)\n"); next }
    cat(" ", patient, "\n")
    p <- tryCatch(make_gam_curve_plot(patient, valid_results),
                  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(p)) print(p)
}

# ------------------------------------------------------------------
# SECTION 5: LOO gene contribution plots
#
# Both DFG (7 genes: BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)
# and HKG (5 genes: ACTB, PGK1, PPIA, RPL13A, SDHA) have full LOO.
# Heatmaps and bar-plot grids are produced for each gene set separately.
# ------------------------------------------------------------------
cat("\n========== LOO gene contribution plots (DFG & HKG) ==========\n")

if (!is.null(sheet_loo) && nrow(sheet_loo) > 0) {
    
    for (gs in c("DFG", "HKG")) {
        
        loo_gs <- sheet_loo[
            sheet_loo$Gene_Set == gs &
                !is.na(sheet_loo$Delta_DE) &
                sheet_loo$Removed_Gene != "(only 1 gene present - LOO not applicable)", ]
        
        if (nrow(loo_gs) == 0) {
            cat("  No valid LOO data for", gs, "- skipping\n")
            next
        }
        
        gs_title <- if (gs == "DFG")
            "DFG composite (BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)"
        else
            "HKG composite (ACTB, PGK1, PPIA, RPL13A, SDHA)"
        
        cat("  Plotting LOO for:", gs_title, "\n")
        
        # Enforce patient display order (MGH152 first ... MGH151 last)
        loo_gs$Patient <- factor(loo_gs$Patient, levels = rev(patients))
        
        # ---- 5a. Heatmap: all patients × genes ----
        p_heat <- ggplot(loo_gs,
                         aes(x = Removed_Gene, y = Patient,
                             fill = Delta_DE * 100)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = sprintf("%.1f", Delta_DE * 100)),
                      size = 2.6, color = "black") +
            scale_fill_distiller(palette = "Spectral",
                                 name   = expression(Delta*"DE (%)"),
                                 labels = number_format(accuracy = 0.1)) +
            labs(title    = paste0("Leave-One-Out Gene Contribution - ", gs, " | FOXM1 GAM"),
                 subtitle = paste0("\u0394DE = DE(full ", gs, ") \u2212 DE(LOO)  |  ",
                                   "red = gene helps fit,  blue = gene suppresses fit\n",
                                   gs_title),
                 x = paste("Removed", gs, "Gene"), y = "Patient") +
            theme_gam() +
            theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                  strip.text  = element_text(face = "bold", size = 10))
        print(p_heat)
        
        # ---- 5b. Per-patient bar plots (2 columns) ----
        bar_plots <- list()
        for (patient in patients) {
            dd_loo <- loo_gs[as.character(loo_gs$Patient) == patient, ]
            if (nrow(dd_loo) == 0) next
            gene_order <- dd_loo %>%
                arrange(desc(Delta_DE)) %>%
                pull(Removed_Gene)
            dd_loo$Removed_Gene <- factor(dd_loo$Removed_Gene, levels = gene_order)
            bar_color <- if (gs == "DFG") "#4B0082" else "#4F75DE"
            p_bar <- ggplot(dd_loo,
                            aes(x = Removed_Gene, y = Delta_DE * 100)) +
                geom_col(width = 0.6, fill = bar_color) +
                geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
                geom_text(aes(label = sprintf("%.1f%%", Delta_DE * 100),
                              vjust = ifelse(Delta_DE >= 0, -0.4, 1.2)),
                          size = 2.2, color = "gray20") +
                labs(title = patient, x = NULL,
                     y = expression(Delta*"DE (%)")) +
                theme_gam(base_size = 7) +
                theme(axis.text.x = element_text(angle = 30, hjust = 1))
            bar_plots[[patient]] <- p_bar
        }
        
        if (length(bar_plots) > 0) {
            cat("  Displaying", gs, "LOO bar plots grid (",
                length(bar_plots), "patients)...\n")
            grid_bars <- tryCatch(
                gridExtra::grid.arrange(
                    grobs = bar_plots,
                    ncol  = 2,
                    top   = grid::textGrob(
                        paste0("LOO ", gs, " Gene Contribution - \u0394DE per Gene per Patient\n",
                               "\u0394DE = DE(full ", gs, ") \u2212 DE(LOO);  positive = gene contributes\n",
                               gs_title),
                        gp = grid::gpar(fontface = "bold", fontsize = 10))),
                error = function(e) {
                    cat("  grid.arrange error:", conditionMessage(e), "\n"); NULL })
            if (!is.null(grid_bars)) print(grid_bars)
        }
    }
    
} else {
    cat("  No LOO data found - skipping\n")
}

# ------------------------------------------------------------------
# SECTION 6: Console summary statistics
# ------------------------------------------------------------------
cat("\n========== Summary Statistics ==========\n")
cat(sprintf("Dataset  : Neftel GBM SS2 - Adult malignant cells\n"))
cat(sprintf("Patients : %s\n", paste(patients, collapse = ", ")))
cat(sprintf("DFG genes: %s\n", paste(c("BIRC5","MKI67","CENPF","TOP2A","PBK","TROAP","NUSAP1"),
                                     collapse = ", ")))
cat(sprintf("HKG genes: %s\n", paste(c("ACTB","PGK1","PPIA","RPL13A","SDHA"),
                                     collapse = ", ")))

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
    
    gs_label <- if (gs == "DFG")
        "DFG (BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)"
    else
        "HKG (ACTB, PGK1, PPIA, RPL13A, SDHA)"
    
    cat("\nGene Set:", gs_label, "(", length(gs_res), "patients)\n")
    cat("  Mean best k:       ", round(mean(k_vals,    na.rm=TRUE), 2), "\n")
    cat("  Mean best PRSS:    ", round(mean(prss_vals, na.rm=TRUE), 4), "\n")
    cat("  Mean DE:           ", round(mean(de_vals,   na.rm=TRUE)*100, 2), "%\n")
    cat("  Mean smooth EDF:   ", round(mean(edf_vals,  na.rm=TRUE), 4), "\n")
    cat("  Mean best iter:    ", round(mean(iter_vals, na.rm=TRUE), 1), "\n")
    
    # Per-patient breakdown (in specified patient order)
    cat("  Per-patient DE:\n")
    for (patient in patients) {
        nm <- paste0(patient, "|", gs)
        r  <- gs_res[[nm]]
        if (is.null(r)) next
        de <- if (!is.null(r$best_model)) summary(r$best_model)$dev.expl else NA_real_
        cat(sprintf("    %-10s  DE = %.2f%%  k = %d  PRSS = %.4f\n",
                    r$patient,
                    ifelse(is.na(de), NA, de * 100),
                    r$best_k,
                    r$best_prss))
    }
}

cat("\n========== Part 5b complete ==========\n")
cat("All plots displayed in RStudio Plots pane (no files saved).\n")
cat("All tabular results are in: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_FOXM1-pos.xlsx\n")




##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 6: GAM Deviance Components & Explanatory Power
#         FOXM1 vs DFG composite & HKG composite
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
# FOXM1-positive adult malignant cells only
#
# Loads:  NeftelGBM_SS2_FOXM1_DFG_GAM_models.rds  (from Part 5)
# Output: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_FOXM1-pos.xlsx  (single file, multiple sheets)
#
#   Sheet 1 - Deviance_Summary      : ND / MD / DE and all component metrics
#                                     (one row per patient × gene set)
#   Sheet 2 - Deviance_Explanations : formula reference for every metric
#   Sheet 3 - Deviance_Detailed     : long-format per-metric breakdown
#   Sheet 4 - EP_All_Combined       : cell-level EP values, all patients × gene sets
#   Sheet 5 - EP_Summary            : patient × gene set summary statistics
#   Sheets 6+ - EP_<Patient>_<GS>  : per-fit cell-level EP (one sheet per fit)
##############################################################################
library(mgcv)
library(writexl)
set.seed(123)

# ------------------------------------------------------------------
# SECTION 1: Load Part 5 models
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data     <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results  <- rds_data$all_results
patients     <- rds_data$patients

valid_results <- Filter(Negate(is.null), all_results)
gene_sets     <- c("DFG", "HKG")

cat("Valid fits loaded:", length(valid_results), "\n")
cat("Patients:", paste(patients, collapse = ", "), "\n\n")

# ------------------------------------------------------------------
# SECTION 2: Deviance component functions
# ------------------------------------------------------------------

# ---- 2a. Extract all deviance components from one model ----
# (w_zero = p_pos for FOXM1-zero cells, 1.0 for positive cells).
# The null model must use identical weights for deviance to be comparable.
extract_deviance_components <- function(model, data, cell_weights) {
    set.seed(123)
    n             <- nrow(data)
    response_var  <- "Expression"
    y             <- data[[response_var]]
    
    fitted_values      <- fitted(model)
    residuals_response <- residuals(model, type = "response")
    
    edf_total   <- sum(model$edf)
    df_residual <- model$df.residual
    df_model    <- n - df_residual
    
    scale_parameter <- model$scale
    
    # Null model: intercept only - same weights as the fitted model
    null_model <- gam(Expression ~ 1, data = data,
                      weights = cell_weights, family = gaussian())
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

# ---- 2b. Run extract_deviance_components across all 40 fits ----
analyze_all_deviances <- function(valid_res) {
    deviance_data <- list()
    for (key in names(valid_res)) {
        res            <- valid_res[[key]]
        parts          <- strsplit(key, "\\|")[[1]]
        patient_id     <- parts[1]
        gene_set_label <- parts[2]
        cat("Processing:", key, "\n")
        tryCatch({
            # Pass cell_weights stored in result for consistent null model
            dc <- extract_deviance_components(res$best_model, res$data,
                                              res$cell_weights)
            dc$patient  <- patient_id
            dc$gene_set <- gene_set_label
            dc$group    <- "Neftel_GBM_SS2"
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
            Patient               = dc$patient,
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
            "mgcv null deviance calculation (weighted)",
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
            "Null deviance from mgcv (weighted, same weights as fitted model)",
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
        extra <- list(
            y_min      = dc$y_range[1],
            y_max      = dc$y_range[2],
            fitted_min = dc$fitted_range[1],
            fitted_max = dc$fitted_range[2]
        )
        vals <- c(dc[scalar_fields], extra)
        data.frame(
            Patient  = dc$patient,
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
# cell_weights: passed from result$cell_weights so the null model is
# fitted on the same weighted scale as the GAM, keeping EP consistent.
extract_ep_for_fit <- function(result, patient_id, gene_set_label) {
    cat("  Extracting EP for:", patient_id, "|", gene_set_label, "\n")
    
    model        <- result$best_model
    data         <- result$data
    cell_weights <- result$cell_weights   # dynamic zero-weights from Part 5
    n_cells      <- nrow(data)
    
    # Overall model deviance explained
    model_dev_explained <- summary(model)$dev.expl
    
    # Null model: same weights as fitted model
    null_model   <- gam(Expression ~ 1, data = data, weights = cell_weights)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(model)
    
    null_residuals  <- data$Expression - null_fitted
    model_residuals <- data$Expression - model_fitted
    
    null_sq_diff  <- null_residuals^2
    model_sq_diff <- model_residuals^2
    
    # Individual cell explanatory power
    # EP_i = 1 - (model_sq_diff_i / null_sq_diff_i)
    # Cells where null_sq_diff == 0 get EP = 0
    ep_value <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )
    
    # Sort cells by EP (highest first)
    sorted_idx        <- order(ep_value, decreasing = TRUE)
    sorted_cell_names <- rownames(data)[sorted_idx]
    
    # Deviance cells: top (model_dev_explained × n) cells by EP
    target_cells       <- round(n_cells * model_dev_explained)
    target_cells       <- max(1L, min(target_cells, n_cells))
    deviance_cells     <- sorted_cell_names[seq_len(target_cells)]
    non_deviance_cells <- sorted_cell_names[(target_cells + 1):n_cells]
    
    output_df <- data.frame(
        Patient              = patient_id,
        Group                = "Neftel_GBM_SS2",
        Gene_Set             = gene_set_label,
        Cell_Index           = rownames(data),
        EP_Value             = ep_value,
        Sorted_Index         = match(rownames(data), sorted_cell_names),
        Deviance_Cell        = ifelse(rownames(data) %in% deviance_cells,
                                      "YES", "NO"),
        Non_Deviance_Cell    = ifelse(rownames(data) %in% non_deviance_cells,
                                      "YES", "NO"),
        FOXM1_Expression     = data$FOXM1,
        Expression           = data$Expression,
        Cell_Weight          = cell_weights,   # recorded for transparency
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

all_ep_list <- list()
per_fit_ep  <- list()

for (key in names(valid_results)) {
    parts          <- strsplit(key, "\\|")[[1]]
    patient_id     <- parts[1]
    gene_set_label <- parts[2]
    cat("Processing EP:", key, "\n")
    tryCatch({
        # cell_weights retrieved from result$cell_weights inside the function
        ep_df <- extract_ep_for_fit(valid_results[[key]], patient_id, gene_set_label)
        all_ep_list[[key]] <- ep_df
        sheet_nm <- paste0("EP_", patient_id, "_", gene_set_label)
        sheet_nm <- substr(sheet_nm, 1, 31)
        per_fit_ep[[sheet_nm]] <- ep_df
    }, error = function(e) {
        cat("  Error:", conditionMessage(e), "\n")
    })
}

ep_combined <- if (length(all_ep_list) > 0)
    do.call(rbind, all_ep_list) else data.frame()
rownames(ep_combined) <- NULL

ep_summary <- if (length(all_ep_list) > 0) {
    rows <- list()
    for (patient in patients) {
        for (gs in gene_sets) {
            key <- paste0(patient, "|", gs)
            df  <- all_ep_list[[key]]
            if (is.null(df)) next
            rows[[key]] <- data.frame(
                Patient              = df$Patient[1],
                Group                = df$Group[1],
                Gene_Set             = df$Gene_Set[1],
                Total_Cells          = nrow(df),
                Mean_EP_Value        = round(mean(df$EP_Value,  na.rm = TRUE), 4),
                Max_EP_Value         = round(max(df$EP_Value,   na.rm = TRUE), 4),
                Min_EP_Value         = round(min(df$EP_Value,   na.rm = TRUE), 4),
                Deviance_Cells_Count = sum(df$Deviance_Cell == "YES"),
                Model_Dev_Explained  = round(unique(df$Model_Dev_Explained), 4),
                stringsAsFactors = FALSE
            )
        }
    }
    do.call(rbind, rows)
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
    per_fit_ep
)

out_file <- "NeftelGBM_SS2_FOXM1_DFG_GAM_Results_FOXM1-pos.xlsx"

tryCatch({
    write_xlsx(export_list, path = out_file)
    cat("\nSaved:", out_file, "\n")
    cat("Sheets:\n")
    cat("  Deviance_Summary      :", nrow(deviance_summary),
        "rows (one per patient x gene set)\n")
    cat("  Deviance_Explanations :", nrow(deviance_expl),
        "rows (formula reference)\n")
    cat("  Deviance_Detailed     :", nrow(deviance_detail),
        "rows (long-format per-metric breakdown)\n")
    cat("  EP_All_Combined       :", nrow(ep_combined),
        "rows (cell-level EP, all fits)\n")
    cat("  EP_Summary            :", nrow(ep_summary),
        "rows (one per patient x gene set)\n")
    for (nm in names(per_fit_ep))
        cat(" ", nm, ":", nrow(per_fit_ep[[nm]]), "rows\n")
}, error = function(e) {
    cat("Error writing Excel:", conditionMessage(e), "\n")
})

# ------------------------------------------------------------------
# SECTION 7: Console summary statistics  
# ------------------------------------------------------------------
cat("\n========== Summary Statistics ==========\n")
cat(sprintf("Dataset  : Neftel GBM SS2 - Adult malignant cells\n"))
cat(sprintf("Patients : %s\n", paste(patients, collapse = ", ")))
cat(sprintf("DFG genes: %s\n",
            paste(c("BIRC5","MKI67","CENPF","TOP2A","PBK","TROAP","NUSAP1"),
                  collapse = ", ")))
cat(sprintf("HKG genes: %s\n",
            paste(c("ACTB","PGK1","PPIA","RPL13A","SDHA"),
                  collapse = ", ")))

if (nrow(deviance_summary) > 0) {
    cat("\nScale parameter (phi):\n")
    cat("  Min: ", round(min(deviance_summary$Scale_Parameter,  na.rm=TRUE), 6), "\n")
    cat("  Max: ", round(max(deviance_summary$Scale_Parameter,  na.rm=TRUE), 6), "\n")
    cat("  Mean:", round(mean(deviance_summary$Scale_Parameter, na.rm=TRUE), 6), "\n")
    
    cat("\nDeviance explained (Dev_Explained_Summary):\n")
    cat("  Min: ", round(min(deviance_summary$Dev_Explained_Summary,  na.rm=TRUE), 4), "\n")
    cat("  Max: ", round(max(deviance_summary$Dev_Explained_Summary,  na.rm=TRUE), 4), "\n")
    cat("  Mean:", round(mean(deviance_summary$Dev_Explained_Summary, na.rm=TRUE), 4), "\n")
    
    cat("\nTotal EDF:\n")
    cat("  Min: ", round(min(deviance_summary$EDF_Total,  na.rm=TRUE), 3), "\n")
    cat("  Max: ", round(max(deviance_summary$EDF_Total,  na.rm=TRUE), 3), "\n")
    cat("  Mean:", round(mean(deviance_summary$EDF_Total, na.rm=TRUE), 3), "\n")
    
    for (gs in gene_sets) {
        sub_df   <- deviance_summary[deviance_summary$Gene_Set == gs, ]
        sub_df   <- sub_df[match(patients, sub_df$Patient), ]
        sub_df   <- sub_df[!is.na(sub_df$Patient), ]
        gs_label <- if (gs == "DFG")
            "DFG (BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)"
        else
            "HKG (ACTB, PGK1, PPIA, RPL13A, SDHA)"
        cat(sprintf("\nGene Set: %s (%d patients)\n", gs_label, nrow(sub_df)))
        cat("  Mean Dev_Explained_Summary:",
            round(mean(sub_df$Dev_Explained_Summary, na.rm=TRUE), 4), "\n")
        cat("  Mean Scale_Parameter:      ",
            round(mean(sub_df$Scale_Parameter, na.rm=TRUE), 6), "\n")
        cat("  Mean EDF_Total:            ",
            round(mean(sub_df$EDF_Total, na.rm=TRUE), 3), "\n")
        
        cat("  Per-patient Dev_Explained_Summary:\n")
        for (i in seq_len(nrow(sub_df))) {
            cat(sprintf("    %-10s  DE = %.4f  EDF = %.3f  phi = %.6f\n",
                        sub_df$Patient[i],
                        sub_df$Dev_Explained_Summary[i],
                        sub_df$EDF_Total[i],
                        sub_df$Scale_Parameter[i]))
        }
    }
}

cat("\n========== Part 6 complete ==========\n")
cat("Output: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_FOXM1-pos.xlsx\n")

##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 7: Inflection Point (IP) Detection - DFG composite only
#         FOXM1 as GOI (x-axis)
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
# FOXM1-positive adult malignant cells only
#
# Method: R(x) GAM zero-crossing with Brent 95% CI
#         + IP Reliability Score (IPRS) replacing categorical tiers
#
# Loads:  NeftelGBM_SS2_FOXM1_DFG_GAM_models.rds  (from Part 5)
# Output: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_5_FOXM1-pos.xlsx  (3 sheets)
#         inflection_points            - named vector keyed "patient|DFG"
#         inflection_points_DFG        - named vector keyed by patient
#         inflection_points_by_patient - named vector keyed by patient
#
# Gene set : DFG composite (BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)
#            HKG excluded
##############################################################################
library(ggplot2)
library(dplyr)
library(mgcv)
library(openxlsx)

# ------------------------------------------------------------------
# SECTION 1: Load Part 5 models
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data    <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results <- rds_data$all_results
patients    <- rds_data$patients   # full 20-patient list (kept for reference)

# All patients
TARGET_PATIENTS <- c("MGH152", "MGH129", "MGH66", "MGH122", "MGH128", "MGH110",  "MGH100", "MGH113", "MGH121", "MGH143", "MGH105", "MGH102", "MGH124", "MGH106", "MGH115", "MGH151", "MGH101", "MGH104", "MGH125", "MGH136")
TARGET_GENE_SET <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)

valid_results <- valid_results[
    names(valid_results) %in%
        paste0(TARGET_PATIENTS, "|", TARGET_GENE_SET)
]
cat("Restricted to:", paste(names(valid_results), collapse = ", "), "\n")
cat("Valid fits:", length(valid_results), "\n\n")

DFG_LABEL   <- "DFG composite"
DFG_MEMBERS <- "BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1"

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

# ---- 3a. SW-gated two-group stat test   ----
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

# ---- 3b. Brent CI on R(x) GAM confidence bands   ----
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
                             newdata = data.frame(FOXM1 = x), se.fit = TRUE)
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

# ---- 3b2. Proximity fallback CI bound   ----
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

# ---- 3c. IPRS component calculation   ----
compute_iprs <- function(ip_value, ip_ci_lower, ip_ci_upper,
                         stat_p_adj, x_seq, r_hat, r_upper, r_lower,
                         residual_data) {
    
    # --- C1: Precision ---
    foxm1_range <- max(x_seq) - min(x_seq)
    
    both_ci  <- !is.na(ip_ci_lower) && !is.na(ip_ci_upper)
    lower_ok <- !is.na(ip_ci_lower)
    upper_ok <- !is.na(ip_ci_upper)
    
    ci_width <- if (both_ci) abs(ip_ci_upper - ip_ci_lower) else NA_real_
    
    C1 <- if (both_ci && foxm1_range > 0) {
        max(0, 1 - ci_width / foxm1_range)
    } else if ((lower_ok || upper_ok) && !is.na(ip_value) && foxm1_range > 0) {
        half_width <- if (lower_ok) abs(ip_value - ip_ci_lower)
        else          abs(ip_ci_upper - ip_value)
        0.5 * max(0, 1 - half_width / foxm1_range)
    } else {
        0
    }
    
    C1_note <- if (both_ci)
        sprintf("Both CI bounds: CI_width=%.4f / FOXM1_range=%.4f → C1=%.4f",
                ci_width, foxm1_range, C1)
    else if (lower_ok || upper_ok) {
        known_bound <- if (lower_ok) ip_ci_lower else ip_ci_upper
        bound_label <- if (lower_ok) "lower" else "upper"
        half_w      <- abs(ip_value - known_bound)
        sprintf("Partial CI (%s only): half_width=%.4f / FOXM1_range=%.4f → 0.5×max(0,1-ratio)=%.4f",
                bound_label, half_w, foxm1_range, C1)
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
        post_cells    <- residual_data$residual[residual_data$FOXM1 >= ip_value]
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
    
    # Hard rule: any component = 0 → Failed (structurally incomplete IP)
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
         foxm1_range = foxm1_range,
         amplitude   = amplitude,
         mean_ribbon = mean_ribbon,
         tv          = tv)
}

# ---- 3d. Core IP detection function ----
# cell_weights: passed from result$cell_weights so the null model used
# for EP/purple-cell classification is on the same weighted scale as the
# fitted GAM. The residual GAM R(x) itself is unweighted (it operates on
# purple-cell residuals only, where all cells have weight 1.0 or near so).
detect_ip <- function(sample_data, gam_model, dev_explained_pct,
                      cell_weights) {
    
    # Step 1: Purple-cell classification (top EP cells)
    # Null model uses same weights as fitted GAM for consistent EP
    null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data,
                               weights = cell_weights))
    model_fitted <- fitted(gam_model)
    expl_power   <- 1 - ((sample_data$Expression - model_fitted)^2 /
                             (sample_data$Expression - null_fitted)^2)
    n_purple  <- round(nrow(sample_data) * dev_explained_pct)
    is_purple <- rep(FALSE, nrow(sample_data))
    is_purple[order(expl_power, decreasing = TRUE)[seq_len(n_purple)]] <- TRUE
    
    # Step 2: Signed residuals for purple cells
    residual      <- sample_data$Expression -
        predict(gam_model, newdata = sample_data)
    purple_idx    <- which(is_purple)
    residual_data <- data.frame(
        FOXM1    = sample_data$FOXM1[purple_idx],
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
                          deriv_at_ip = NA_real_,
                          ci_width    = NA_real_,
                          foxm1_range = NA_real_,
                          amplitude   = NA_real_,
                          mean_ribbon = NA_real_,
                          tv          = NA_real_)
        list(ip_value = NA_real_, manual_flag = TRUE, status = status,
             stat_test = NA_character_, stat_p_raw = NA_real_,
             stat_p_adj = NA_real_, both_normal = NA,
             ip_ci_lower = NA_real_, ip_ci_upper = NA_real_,
             neg_pos_count = 0L,
             iprs = iprs_fail,
             residual_data = residual_data,
             residual_pred = data.frame(FOXM1 = numeric(0), r_hat = numeric(0),
                                        se = numeric(0), r_upper = numeric(0),
                                        r_lower = numeric(0)))
    }
    
    if (nrow(residual_data) < MIN_PURPLE_CELLS ||
        length(unique(residual_data$FOXM1)) < 3L) {
        cat("  *** INSUFFICIENT PURPLE CELLS - skipping IP detection ***\n")
        return(FAILED_list("insufficient_purple_cells"))
    }
    
    # Step 3: Residual GAM R(x) - unweighted (purple cells only)
    set.seed(123)
    n_unique     <- length(unique(residual_data$FOXM1))
    k_resid      <- min(20L, n_unique - 1L)
    residual_gam <- tryCatch(
        gam(residual ~ s(FOXM1, bs = "tp", k = k_resid),
            data = residual_data, method = "REML"),
        error = function(e) {
            cat("  Fallback k=8:", conditionMessage(e), "\n")
            gam(residual ~ s(FOXM1, bs = "tp",
                             k = min(8L, n_unique - 1L)),
                data = residual_data, method = "REML")
        }
    )
    gam_sum <- summary(residual_gam)
    cat("  R(x) GAM: edf =", sprintf("%.2f", gam_sum$edf),
        "| Dev =", sprintf("%.1f%%", gam_sum$dev.expl * 100),
        "| p =", sprintf("%.2e", gam_sum$s.table[4]), "\n")
    
    # Step 4: 2000-pt grid prediction
    x_seq <- seq(min(residual_data$FOXM1), max(residual_data$FOXM1),
                 length.out = 2000L)
    pred  <- predict(residual_gam,
                     newdata = data.frame(FOXM1 = x_seq), se.fit = TRUE)
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
        pre_r  <- residual_data$residual[residual_data$FOXM1 <  ip_value]
        post_r <- residual_data$residual[residual_data$FOXM1 >= ip_value]
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
    iprs <- compute_iprs(ip_value       = ip_value,
                         ip_ci_lower    = ip_ci_lower,
                         ip_ci_upper    = ip_ci_upper,
                         stat_p_adj     = st$stat_p_adj,
                         x_seq          = x_seq,
                         r_hat          = r_hat,
                         r_upper        = r_upper,
                         r_lower        = r_lower,
                         residual_data  = residual_data)
    
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
            FOXM1   = x_seq,
            r_hat   = r_hat,
            se      = r_se,
            r_upper = r_upper,
            r_lower = r_lower
        )
    )
}

# ------------------------------------------------------------------
# SECTION 4: Run IP detection
# ------------------------------------------------------------------
set.seed(123)
cat("\n=============================================\n")
cat(sprintf("IP DETECTION | Zero-crossing + Brent CI + IPRS | STAT_ALPHA=%.2f\n",
            STAT_ALPHA))
cat("Dataset  : Neftel GBM SS2 - Adult malignant cells\n")
cat("Patients :", paste(TARGET_PATIENTS, collapse = ", "), "\n")
cat("Gene set : DFG composite (", DFG_MEMBERS, ")\n")
cat("=============================================\n\n")

ip_results <- list()

for (key in paste0(TARGET_PATIENTS, "|", TARGET_GENE_SET)[
    paste0(TARGET_PATIENTS, "|", TARGET_GENE_SET) %in% names(valid_results)]) {
    parts          <- strsplit(key, "\\|")[[1]]
    patient_id     <- parts[1]
    gene_set_label <- parts[2]
    cat("\n--- Processing:", patient_id, "|", DFG_LABEL, "---\n")
    
    res         <- valid_results[[key]]
    sample_data <- res$data
    gam_model   <- res$best_model
    dev_expl    <- summary(gam_model)$dev.expl
    
    # Pass cell_weights from result for consistent weighted null model
    ip_results[[key]] <- detect_ip(sample_data, gam_model, dev_expl,
                                   res$cell_weights)
}

# ------------------------------------------------------------------
# SECTION 5: Console summary  
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("IP SUMMARY - IP Reliability Score (IPRS)\n")
cat("=============================================\n\n")
cat(sprintf("%-22s %8s %8s %8s %6s %6s %6s %6s %6s %6s  %-14s  %s\n",
            "Fit (Patient|GS)", "IP", "CI_Lo", "CI_Hi",
            "C1", "C2", "C3", "C4", "C5", "IPRS", "Tier", "Flag"))
cat(strrep("-", 130), "\n")

for (key in names(ip_results)) {
    r <- ip_results[[key]]; if (is.null(r)) next
    iprs <- r$iprs
    cat(sprintf("%-22s %8s %8s %8s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  %-14s  %s\n",
                key,
                if (is.na(r$ip_value))    "NA" else sprintf("%.4f", r$ip_value),
                if (is.na(r$ip_ci_lower)) "NA" else sprintf("%.4f", r$ip_ci_lower),
                if (is.na(r$ip_ci_upper)) "NA" else sprintf("%.4f", r$ip_ci_upper),
                iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5, iprs$IPRS,
                iprs$IPRS_tier,
                iprs$IPRS_flag))
}

cat("\n")
for (key in names(ip_results)) {
    r <- ip_results[[key]]; if (is.null(r)) next
    if (r$iprs$IPRS_flag != "None")
        cat(sprintf("  [FLAG] %s | IPRS=%.3f (%s) → %s\n",
                    key, r$iprs$IPRS, r$iprs$IPRS_tier, r$iprs$IPRS_flag))
}

# ------------------------------------------------------------------
# SECTION 6: Excel export (3 sheets)  
# ------------------------------------------------------------------
cat("\nExporting results...\n")

ip_sheet <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]
    parts <- strsplit(key, "\\|")[[1]]
    iprs  <- r$iprs
    data.frame(
        Patient          = parts[1],
        Gene_Set         = parts[2],
        Gene_Set_Label   = DFG_LABEL,
        Gene_Members     = DFG_MEMBERS,
        IP_Value         = r$ip_value,
        CI_Lower_95      = r$ip_ci_lower,
        CI_Upper_95      = r$ip_ci_upper,
        CI_Width         = iprs$ci_width,
        FOXM1_Range      = iprs$foxm1_range,
        CI_Lower_Method  = if (is.null(r$ci_lower_method)) NA_character_ else r$ci_lower_method,
        CI_Upper_Method  = if (is.null(r$ci_upper_method)) NA_character_ else r$ci_upper_method,
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
        Method = "ZeroCrossing|BrentCI|IPRS=(C1+C2+C3+C4+C5)/5|C2=tanh(-log10p)|C5=waviness*sign_fidelity|SW_welch_t_OR_mannwhitney_u",
        stringsAsFactors = FALSE
    )
}))
rownames(ip_sheet) <- NULL

resid_pred_df <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]; if (is.null(r)) return(NULL)
    if (nrow(r$residual_pred) == 0L) return(NULL)
    parts <- strsplit(key, "\\|")[[1]]
    p     <- r$residual_pred[seq(1L, nrow(r$residual_pred), by = 10L), ]
    data.frame(Patient = parts[1], Gene_Set = parts[2], p,
               stringsAsFactors = FALSE)
}))
if (is.null(resid_pred_df)) resid_pred_df <- data.frame()
rownames(resid_pred_df) <- NULL

method_notes <- data.frame(
    Step = paste0("Step ", 1:8),
    Description = c(
        "Cell classification: EP = 1 - model_residual\u00b2 / null_residual\u00b2; top (dev_expl \u00d7 n) cells are 'purple'. Null model uses same cell weights as fitted GAM (dynamic zero-weighting from Part 5) for consistent EP.",
        "Signed residuals: residual = observed Expression \u2212 GAM_fitted, for each purple cell.",
        "Residual GAM R(x): thin-plate spline (REML) fitted to (FOXM1, residual) for purple cells; k = min(20, n_unique\u22121). Unweighted - all purple cells contribute equally to R(x).",
        "Fine grid: R(x) and SE evaluated at 2000 equally-spaced FOXM1 values.",
        paste0("IP = first neg\u2192pos zero-crossing of R(x) GAM (linear interpolation). ",
               "If R(x) has no neg\u2192pos crossing, IP=NA and all IPRS components default to 0."),
        paste0("95% CI primary (Brent): CI_Lower = zero of [R_hat + 1.96\u00d7SE]; CI_Upper = zero of [R_hat - 1.96\u00d7SE]; anchored near IP. ",
               "Fallback (proximity): if Brent returns NA for one bound (confidence band does not cross zero on that side due to data sparsity), ",
               "the fallback finds the x-value where the band is closest to zero in the pre-IP (lower) or post-IP (upper) region. ",
               "Fallback activates only when: (1) Brent=NA, (2) IP exists, (3) search region has \u22653 grid points, (4) band truly has no zero-crossing in that region. ",
               "CI_Lower_Method / CI_Upper_Method column records 'brent' or 'proximity_fallback' for each bound."),
        paste0("Stat test: Shapiro-Wilk (\u03b1=", SHAPIRO_ALPHA,
               ") per group; both normal \u2192 Welch t (one-sided post>pre); ",
               "either non-normal \u2192 Mann-Whitney U; BH FDR applied."),
        paste0("IPRS = (C1 + C2 + C3 + C4 + C5) / 5  [range 0-1]. ",
               "C1 (Precision): both CI bounds \u2192 max(0, 1-CI_Width/FOXM1_range); ",
               "one bound NA \u2192 0.5\u00d7max(0, 1-half_width/FOXM1_range) [partial CI penalty]; both bounds NA \u2192 0. ",
               "C2 (Statistical) = tanh(scale_c2 * (-log10(p_adj))); scale_c2 = atanh(0.5)/(-log10(1e-5)) \u2248 0.110; C2=0.5 at p=0.05 (standard significance anchor); parallel structure to C3=tanh(|R\u2019(IP)|). ",
               "C3 (Abruptness) = tanh(|R'(IP)|) via central difference on 2000-pt grid. ",
               "C4 (GAM Ribbon) = exp(-mean_ribbon_width / amplitude); ",
               "mean_ribbon_width = mean pointwise width of R(x) 95% CI band; amplitude = diff(range(r_hat)); ",
               "lambda=1 by first principles (C4=exp(-1)\u22480.368 when ribbon width equals amplitude). ",
               "C5 (Smoothness) = exp(-(TV/amplitude - 1)); ",
               "TV = total variation = sum(|diff(r_hat)|); ",
               "mu=1 by first principles (C5=1 for monotone curve where TV=amplitude; ",
               "decays exp(-1) per unit excess waviness); ",
               "C5=0 if amplitude < 1e-6 (flat line). ",
               "Equal weights assigned a priori (principle of insufficient reason): no empirical basis ",
               "to prioritise any single dimension of IP quality over others. ",
               "Tiers: \u22650.80 Very Strong | 0.60-0.79 Strong | 0.40-0.59 Moderate [visual review] | ",
               "0.20-0.39 Weak [visual review] | <0.20 Failed [no IP - do not proceed]. ",
               "All IP/CI lines plotted in green (#2E8B57) regardless of IPRS tier.")
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
setColWidths(wb, "IP_Results", cols = 1:2,   widths = 14)
setColWidths(wb, "IP_Results", cols = 3:4,   widths = 20)
setColWidths(wb, "IP_Results", cols = 5:16,  widths = 16)
setColWidths(wb, "IP_Results", cols = 17:30, widths = 16)
setColWidths(wb, "IP_Results", cols = 31,    widths = 90)
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

saveWorkbook(wb, "NeftelGBM_SS2_FOXM1_DFG_GAM_Results_5_FOXM1-pos.xlsx", overwrite = TRUE)
cat("Exported: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_5_FOXM1-pos.xlsx\n")

# ------------------------------------------------------------------
# SECTION 7: Extract IP values into environment  
# ------------------------------------------------------------------
ip_table <- do.call(rbind, lapply(names(ip_results), function(key) {
    r     <- ip_results[[key]]
    parts <- strsplit(key, "\\|")[[1]]
    data.frame(
        Patient     = parts[1],
        Gene_Set    = parts[2],
        IP_Value    = if (!is.null(r)) r$ip_value          else NA_real_,
        IP_CI_Lower = if (!is.null(r)) r$ip_ci_lower       else NA_real_,
        IP_CI_Upper = if (!is.null(r)) r$ip_ci_upper       else NA_real_,
        CI_Width    = if (!is.null(r)) r$iprs$ci_width     else NA_real_,
        C1          = if (!is.null(r)) r$iprs$C1           else NA_real_,
        C2          = if (!is.null(r)) r$iprs$C2           else NA_real_,
        C3          = if (!is.null(r)) r$iprs$C3           else NA_real_,
        C4          = if (!is.null(r)) r$iprs$C4           else NA_real_,
        C5          = if (!is.null(r)) r$iprs$C5           else NA_real_,
        IPRS        = if (!is.null(r)) r$iprs$IPRS         else NA_real_,
        IPRS_Tier   = if (!is.null(r)) r$iprs$IPRS_tier    else NA_character_,
        IPRS_Flag   = if (!is.null(r)) r$iprs$IPRS_flag    else NA_character_,
        stringsAsFactors = FALSE
    )
}))
rownames(ip_table) <- NULL

inflection_points <- setNames(ip_table$IP_Value,
                              paste0(ip_table$Patient, "|", ip_table$Gene_Set))

inflection_points_DFG <- setNames(
    ip_table$IP_Value[ip_table$Gene_Set == "DFG"],
    ip_table$Patient[ip_table$Gene_Set == "DFG"]
)

inflection_points_by_patient <- setNames(
    sapply(TARGET_PATIENTS, function(p) inflection_points_DFG[p]),
    TARGET_PATIENTS
)

cat("\nFinal IP values with IPRS:\n")
print(ip_table[, c("Patient", "Gene_Set", "IP_Value",
                   "IP_CI_Lower", "IP_CI_Upper",
                   "C1", "C2", "C3", "C4", "C5", "IPRS", "IPRS_Tier", "IPRS_Flag")])

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
# SECTION 9: Scatter plot
# Purple-cell colouring uses weighted null model (consistent with GAM fit)
# ------------------------------------------------------------------
create_scatter <- function(key) {
    r <- ip_results[[key]]; if (is.null(r)) return(NULL)
    parts      <- strsplit(key, "\\|")[[1]]
    patient_id <- parts[1]
    
    res          <- valid_results[[key]]
    sample_data  <- res$data
    gam_model    <- res$best_model
    cell_weights <- res$cell_weights   # use stored weights for null model
    dev_expl     <- summary(gam_model)$dev.expl
    
    # Weighted null model - consistent with detect_ip purple classification
    null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data,
                               weights = cell_weights))
    model_fitted <- fitted(gam_model)
    expl_power   <- 1 - ((sample_data$Expression - model_fitted)^2 /
                             (sample_data$Expression - null_fitted)^2)
    dev_cells    <- rownames(sample_data)[
        order(expl_power, decreasing = TRUE)[
            seq_len(round(nrow(sample_data) * dev_expl))]]
    
    plot_data <- data.frame(
        FOXM1      = sample_data$FOXM1,
        Expression = sample_data$Expression,
        Group      = ifelse(rownames(sample_data) %in% dev_cells,
                            "Dev explained", "Non-dev explained")
    )
    plot_data <- plot_data[order(plot_data$Group == "Dev explained"), ]
    
    pred_data        <- data.frame(FOXM1 = seq(min(plot_data$FOXM1),
                                               max(plot_data$FOXM1),
                                               length.out = 100L))
    pg               <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit    <- pg$fit
    pred_data$se.fit <- pg$se.fit
    
    lc        <- ip_line_color
    title_col <- if (r$manual_flag) manual_color else "black"
    
    title <- if (r$manual_flag || is.na(r$ip_value)) {
        paste0(patient_id, " | ", DFG_LABEL, "  *** NO IP DETECTED ***")
    } else {
        paste0(patient_id, " | ", DFG_LABEL,
               "  IP: ", sprintf("%.4f", r$ip_value),
               "  95% CI [",
               ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
               ", ",
               ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
               "]")
    }
    
    p <- ggplot() +
        geom_point(data = plot_data,
                   aes(x = FOXM1, y = Expression, color = Group),
                   size = 1.8, alpha = 0.65) +
        geom_ribbon(data = pred_data,
                    aes(x = FOXM1, ymin = fit - 1.96 * se.fit,
                        ymax = fit + 1.96 * se.fit),
                    fill = "#FFCC99", alpha = 0.20) +
        geom_line(data = pred_data, aes(x = FOXM1, y = fit),
                  color = "#FFCC99", linewidth = 1.2)
    
    if (!is.na(r$ip_value))
        p <- p + geom_vline(xintercept = r$ip_value,
                            linetype = "dashed", color = "black", linewidth = 0.5)
    
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
             x = "FOXM1 Expression (log2)",
             y = paste0(DFG_LABEL, " Mean Expression (log2)")) +
        scale_x_continuous(
            limits = range(sample_data$FOXM1, na.rm = TRUE) +
                c(-1, 1) * 0.05 * diff(range(sample_data$FOXM1, na.rm = TRUE))
        )
}

# ------------------------------------------------------------------
# SECTION 10: R(x) diagnostic plot  
# ------------------------------------------------------------------
create_rx <- function(key) {
    r <- ip_results[[key]]; if (is.null(r)) return(NULL)
    
    parts      <- strsplit(key, "\\|")[[1]]
    patient_id <- parts[1]
    
    pred      <- r$residual_pred
    rd        <- r$residual_data
    iprs      <- r$iprs
    lc        <- ip_line_color
    title_col <- if (r$manual_flag) manual_color else "black"
    
    n_bins_rx <- 100L
    p_foxm1   <- rd$FOXM1
    
    if (nrow(rd) > 0L) {
        if (length(unique(p_foxm1)) > n_bins_rx) {
            breaks_rx  <- seq(min(p_foxm1), max(p_foxm1),
                              length.out = n_bins_rx + 1L)
            bin_ids_rx <- cut(p_foxm1, breaks = breaks_rx,
                              include.lowest = TRUE, labels = FALSE)
            bin_ctr_rx <- (breaks_rx[-length(breaks_rx)] + breaks_rx[-1L]) / 2
            bin_w_rx   <- breaks_rx[2L] - breaks_rx[1L]
        } else {
            uq_rx      <- sort(unique(p_foxm1))
            bin_ids_rx <- match(p_foxm1, uq_rx)
            bin_ctr_rx <- uq_rx
            bin_w_rx   <- if (length(uq_rx) > 1L) min(diff(uq_rx)) else 0.1
        }
        
        bs_rx <- data.frame(
            foxm1_center  = bin_ctr_rx,
            mean_residual = sapply(seq_along(bin_ctr_rx), function(b) {
                idx <- which(bin_ids_rx == b)
                if (length(idx) > 0L) mean(rd$residual[idx]) else NA_real_
            })
        )
        bs_rx    <- bs_rx[!is.na(bs_rx$mean_residual), ]
        bar_w_rx <- if (nrow(bs_rx) > 1L)
            diff(range(bs_rx$foxm1_center)) / nrow(bs_rx) * 0.8
        else bin_w_rx * 0.8
    } else {
        bs_rx    <- data.frame(foxm1_center  = numeric(0),
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
               length(unique(rd$FOXM1)), " unique FOXM1 values",
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
                          aes(x = foxm1_center, y = mean_residual),
                          stat = "identity", fill = "#4B0082", color = "#3D0065",
                          alpha = 0.40, width = bar_w_rx, linewidth = 0.4)
    
    if (nrow(pred) > 0L) {
        p <- p +
            geom_ribbon(data = pred,
                        aes(x = FOXM1, ymin = r_lower, ymax = r_upper),
                        fill = "#FFCC99", alpha = 0.35) +
            geom_line(data = pred, aes(x = FOXM1, y = r_hat),
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
            title    = paste0(patient_id, " | ", DFG_LABEL,
                              " \u2014 R(x) = E[residual | purple, FOXM1=x]",
                              "  [Zero-crossing | Brent CI | IPRS]"),
            subtitle = subtitle,
            x = "FOXM1 Expression (log2)",
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

tier_counts <- table(sapply(ip_results, function(r) r$iprs$IPRS_tier))
cat(sprintf("\nIPRS tier summary - %s (%d Tier1+Tier2 patients):\n",
            DFG_LABEL, length(ip_results)))
for (tier in c("Very Strong", "Strong", "Moderate", "Weak", "Failed"))
    cat(sprintf("  %-12s : %d\n", tier,
                if (tier %in% names(tier_counts)) tier_counts[[tier]] else 0L))

cat("\n========== Part 7 complete ==========\n")
cat("\nOutput: NeftelGBM_SS2_FOXM1_DFG_GAM_Results_FOXM1-pos.xlsx\n")
cat("Environment: inflection_points / inflection_points_DFG / inflection_points_by_patient\n")

##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 8: DEG Analysis - Pre/Post-IP TREP cells (purple) vs non-TREP (gray)
#         FOXM1 as GOI (x-axis) | DFG composite (y-axis)
# Dataset: Neftel et al. 2019 - GBM Smart-seq2 (SS2)
#          FOXM1-positive adult malignant cells only
#
# Requires:
#   - NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds  (from Part 5)
#   - inflection_points                 (named numeric vector from Part 7,
#                                        keyed by "patient|DFG", e.g.
#                                        inflection_points["MGH152|DFG"])
#   - gbm_results$seurat_obj            (Seurat object; cell barcodes must
#                                        match rownames in result$data)
#
# Output files:
#   DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Up.xlsx    (upregulated in TREP)
#   DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Down.xlsx  (downregulated in TREP)
##############################################################################

library(Seurat)
library(writexl)
library(dplyr)
library(mgcv)
library(future)
library(parallelly)

# ------------------------------------------------------------------
# SECTION 1: Load Part 5 models and Seurat object
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data    <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results <- rds_data$all_results
patients    <- rds_data$patients

# Five patients that passed dual-filter or partial-pass correlation criteria
# (Spearman >= 0.50 AND/OR Kendall >= 0.40 within FOXM1-positive cells,
#  n_pos >= 50; source: Part 4e / Part 5c D1 dual-filter table)
TARGET_PATIENTS <- c("MGH152", "MGH110", "MGH66", "MGH100", "MGH121")
TARGET_GENE_SET <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)
valid_results <- valid_results[
    names(valid_results) %in%
        paste0(TARGET_PATIENTS, "|", TARGET_GENE_SET)
]

cat("DFG fits retained for analysis:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n")

# Seurat object - FOXM1-positive adult malignant GBM cells
seurat_obj <- gbm_results$seurat_obj

# ------------------------------------------------------------------
# SECTION 2: Check prerequisites
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 7 (IP detection) which creates a named numeric\n",
        "vector 'inflection_points' keyed by \"patient|DFG\".\n",
        "Example: inflection_points[\"MGH152|DFG\"]"
    ))
}

# Validate that IP values exist for the target patients
missing_ips <- TARGET_PATIENTS[
    is.na(inflection_points[paste0(TARGET_PATIENTS, "|", TARGET_GENE_SET)])
]
if (length(missing_ips) > 0) {
    warning(paste0(
        "No valid IP values found for: ",
        paste(missing_ips, collapse = ", "),
        "\nThese patients will be skipped in DEG analysis."
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

# Per-patient relaxed thresholds hook (extend as needed).
# Example: patient_deg_overrides[["MGH100"]] <- list(p_cutoff = 0.05,
#                                                     logfc_cutoff = 0.1)
patient_deg_overrides <- list()

# ------------------------------------------------------------------
# SECTION 5: Helper - classify_cells
#
# Classifies FOXM1-positive cells into TREP / non-TREP using the GAM
# model stored in valid_results[[key]].
# TREP = deviance-explained cells (top EP cells, purple in Part 7).
# EP formula is identical to Part 7 scatter plot logic.
# cell_weights (all 1.0 for FOXM1-positive cells) passed for consistency.
#
# Returns a data.frame with columns:
#   barcode | FOXM1 | is_trep | region | subpopulation
# ------------------------------------------------------------------
classify_cells <- function(patient, result, ip_value) {
    
    gam_model    <- result$best_model
    sample_data  <- result$data        # data.frame: FOXM1, Expression;
    # rownames = cell barcodes
    cell_weights <- result$cell_weights  # all 1.0 - FOXM1-positive only
    
    dev_explained <- summary(gam_model)$dev.expl
    
    # Per-cell explanatory power - null model uses same weights as GAM
    # (all 1.0 for FOXM1-positive cells; consistent with Parts 6 & 7)
    null_model   <- gam(Expression ~ 1, data = sample_data,
                        weights = cell_weights)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(gam_model)
    
    null_sq_diff  <- (sample_data$Expression - null_fitted)^2
    model_sq_diff <- (sample_data$Expression - model_fitted)^2
    
    explanatory_power <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff),
        0
    )
    
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    target_cells   <- round(nrow(sample_data) * dev_explained)
    target_cells   <- max(1L, min(target_cells, nrow(sample_data)))
    
    is_trep <- rep(FALSE, nrow(sample_data))
    is_trep[sorted_indices[seq_len(target_cells)]] <- TRUE
    
    cell_barcodes <- rownames(sample_data)
    foxm1_values  <- sample_data$FOXM1
    
    classification <- data.frame(
        barcode = cell_barcodes,
        FOXM1   = foxm1_values,
        is_trep = is_trep,
        region  = ifelse(foxm1_values < ip_value, "Pre_IP", "Post_IP"),
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
#
# Wilcoxon rank-sum via Seurat FindMarkers, TREP vs non-TREP.
# Minimum cell guard (< 3 cells per group) triggers early exit.
# ------------------------------------------------------------------
run_deg_analysis <- function(seurat_obj, trep_cells, nontrep_cells,
                             deg_params, region_label, patient_label) {
    
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
    
    # Cross-check barcodes against Seurat object
    valid_trep    <- trep_cells[trep_cells %in% colnames(seurat_obj)]
    valid_nontrep <- nontrep_cells[nontrep_cells %in% colnames(seurat_obj)]
    
    if (length(valid_trep) < MIN_RELIABLE_CELLS ||
        length(valid_nontrep) < MIN_RELIABLE_CELLS) {
        cat("    Skipping: too few barcodes matched in seurat_obj\n")
        cat("    (TREP matched:", length(valid_trep),
            "| non-TREP matched:", length(valid_nontrep), ")\n")
        return(NULL)
    }
    
    if (length(valid_trep)    != length(trep_cells) ||
        length(valid_nontrep) != length(nontrep_cells)) {
        cat("    Warning: some barcodes not found in seurat_obj\n")
        cat("    (TREP:", length(trep_cells), "->", length(valid_trep),
            "| non-TREP:", length(nontrep_cells), "->", length(valid_nontrep), ")\n")
    }
    
    Idents(seurat_obj, cells = valid_trep)    <- "TREP"
    Idents(seurat_obj, cells = valid_nontrep) <- "nonTREP"
    
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
cat("Neftel GBM SS2 | FOXM1 vs DFG composite\n")
cat("FOXM1-positive adult malignant cells only\n")
cat("Patients:", paste(TARGET_PATIENTS, collapse = ", "), "\n")
cat("Filter criteria: dual-filter pass OR partial pass\n")
cat("  (Spearman >= 0.50 AND/OR Kendall >= 0.40, n_pos >= 50)\n")
cat("=============================================\n\n")

all_deg_results <- list()

for (patient in TARGET_PATIENTS) {
    
    key <- paste0(patient, "|", TARGET_GENE_SET)
    cat("\n========== Processing:", patient, "(key:", key, ") ==========\n")
    
    result <- valid_results[[key]]
    if (is.null(result)) {
        cat("  Skipping:", patient, "- no valid DFG fit in valid_results\n")
        next
    }
    
    ip_value <- inflection_points[key]
    if (is.na(ip_value) || !is.numeric(ip_value)) {
        cat("  Skipping:", patient, "- no valid IP value in inflection_points\n")
        next
    }
    cat("  IP:", sprintf("%.4f", ip_value), "\n")
    
    # Resolve per-patient parameters
    deg_params <- default_deg_params
    if (!is.null(patient_deg_overrides[[patient]])) {
        deg_params <- modifyList(deg_params, patient_deg_overrides[[patient]])
        cat("  Using custom DEG thresholds for", patient, "\n")
    }
    
    # Classify FOXM1-positive cells into TREP / non-TREP × Pre/Post-IP
    cell_class <- classify_cells(patient, result, ip_value)
    
    cat("  Cell counts per subpopulation:\n")
    print(table(cell_class$subpopulation))
    
    for (region in c("Pre_IP", "Post_IP")) {
        
        region_label <- gsub("_", "-", region)   # "Pre-IP" or "Post-IP"
        cat("\n  ---", region_label, "region ---\n")
        
        trep_cells    <- cell_class$barcode[
            cell_class$subpopulation == paste0(region, "_TREP")]
        nontrep_cells <- cell_class$barcode[
            cell_class$subpopulation == paste0(region, "_nonTREP")]
        
        deg <- run_deg_analysis(seurat_obj, trep_cells, nontrep_cells,
                                deg_params, region_label, patient)
        
        if (is.null(deg)) next
        
        key_deg                    <- paste0(patient, "_", region)
        all_deg_results[[key_deg]] <- deg
    }
    
    cat("  [", patient, "complete]\n")
}

# ------------------------------------------------------------------
# SECTION 9: Export DEG lists - per-patient-region sheets
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("EXPORTING DEG RESULTS\n")
cat("=============================================\n\n")

# ── Upregulated in TREP ──────────────────────────────────────────
deg_list_up <- list()

for (k in names(all_deg_results)) {
    deg <- all_deg_results[[k]]
    
    if (!is.null(deg$up_in_trep) && nrow(deg$up_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_upTREP"), names(deg_list_up))
        deg_list_up[[sn]] <- deg$up_in_trep
    }
    
    if (!is.null(deg$down_in_trep) && nrow(deg$down_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_upNonTREP"), names(deg_list_up))
        deg_list_up[[sn]] <- deg$down_in_trep
    }
}

if (length(deg_list_up) > 0) {
    write_xlsx(deg_list_up,
               path = "DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Up.xlsx")
    cat("Exported", length(deg_list_up),
        "sheets → DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Up.xlsx\n")
} else {
    cat("No upregulated DEG results found.\n")
}

# ── Downregulated in TREP ────────────────────────────────────────
deg_list_down <- list()

for (k in names(all_deg_results)) {
    deg <- all_deg_results[[k]]
    
    if (!is.null(deg$down_in_trep) && nrow(deg$down_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_downTREP"), names(deg_list_down))
        deg_list_down[[sn]] <- deg$down_in_trep
    }
    
    if (!is.null(deg$up_in_trep) && nrow(deg$up_in_trep) > 0) {
        sn <- make_sheet_name(paste0(k, "_downNonTREP"), names(deg_list_down))
        deg_list_down[[sn]] <- deg$up_in_trep
    }
}

if (length(deg_list_down) > 0) {
    write_xlsx(deg_list_down,
               path = "DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Down.xlsx")
    cat("Exported", length(deg_list_down),
        "sheets → DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Down.xlsx\n")
} else {
    cat("No downregulated DEG results found.\n")
}

# ------------------------------------------------------------------
# SECTION 10: Final summary
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================\n\n")

cat("Dataset    : Neftel GBM SS2 - FOXM1-positive adult malignant cells\n")
cat("GOI        : FOXM1 (x-axis)\n")
cat("Gene set   : DFG composite (BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1)\n")
cat("Patients   :", paste(TARGET_PATIENTS, collapse = ", "), "\n")
cat("Filter     : dual-filter pass (MGH152, MGH66) or partial pass (MGH110, MGH100, MGH121)\n")
cat("IP source  : inflection_points[\"patient|DFG\"] (from Part 7)\n\n")

cat("Files generated:\n")
cat("  • DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Up.xlsx\n")
cat("    (upregulated in TREP - sheets: patient_region_upTREP / upNonTREP)\n")
cat("  • DEG_FOXM1_DFG_PrePost_IP_TREP_vs_nonTREP_Down.xlsx\n")
cat("    (downregulated in TREP - sheets: patient_region_downTREP / downNonTREP)\n\n")

cat("DEG criteria (default):\n")
cat("  p < 0.01, |log2FC| > 0.2, min.pct = 0.10,\n")
cat("  min 3 cells per group, max 500 DEGs\n")
cat("  Reliability guard: both groups must have >= 3 cells\n")
cat("  Barcode validation: cross-checked against seurat_obj before FindMarkers\n\n")

cat(sprintf("Parallelization: %d workers | RAM limit: %.1f GB\n\n",
            n_cores, ram_limit_bytes / 1024^3))

cat("========== Part 8 complete ==========\n")



##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 9: Monocle3 Trajectory
#         (TREP and non-TREP Cells, Quantitative Analysis)
#         FOXM1-positive adult malignant GBM cells
#
# Requires:
#   - NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds  (from Part 5)
#   - inflection_points          (named numeric vector from Part 7,
#                                 keyed by "patient|DFG")
#   - gbm_results$seurat_obj     (Seurat object with scRNA-seq data;
#                                 cell barcodes must match rownames
#                                 in result$data from the GAM models)
#
# Patients analysed (DFG gene set only - dual/partial filter passes):
#   MGH152, MGH110, MGH66, MGH100, MGH121
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
# SECTION 1: Load Part 5 models
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data    <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results <- rds_data$all_results

TARGET_PATIENTS <- c("MGH152", "MGH66", "MGH110", "MGH100", "MGH104", "MGH121")
gene_set_label  <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)
valid_results <- valid_results[
    names(valid_results) %in%
        paste0(TARGET_PATIENTS, "|", gene_set_label)
]

cat("DFG fits retained for analysis:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n")

# Seurat object - FOXM1-positive adult malignant GBM cells
seurat_obj <- gbm_results$seurat_obj

# ------------------------------------------------------------------
# SECTION 2: Check prerequisites
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 7 (IP detection) which creates a named numeric\n",
        "vector 'inflection_points' keyed by \"patient|DFG\".\n",
        "Example: inflection_points[\"MGH152|DFG\"]"
    ))
}

# ------------------------------------------------------------------
# SECTION 3: Monocle3 trajectory per patient
# ------------------------------------------------------------------
create_monocle3_trajectory <- function(patient) {
    key <- paste0(patient, "|", gene_set_label)
    cat("=== CREATING MONOCLE3 TRAJECTORY FOR", patient, "(key:", key, ") ===\n")

    result <- valid_results[[key]]
    if (is.null(result)) {
        cat("  Skipping:", patient, "- no valid DFG fit\n")
        return(NULL)
    }

    # ── GAM data ──────────────────────────────────────────────────
    original_gam_data <- result$data   # data.frame: FOXM1, Expression;
                                       # rownames = cell barcodes
    if (is.null(original_gam_data)) {
        stop("GAM data not found for patient ", patient)
    }
    cat("GAM data dimensions:", nrow(original_gam_data), "x",
        ncol(original_gam_data), "\n")

    # ── IP value ──────────────────────────────────────────────────
    ip_value <- inflection_points[key]
    if (is.na(ip_value)) {
        stop("No IP value found in inflection_points for key ", key,
             ". Ensure Part 7 has been run.")
    }
    cat("IP value for", patient, ":", ip_value, "\n")

    # ── Cell classification (identical logic to Parts 6, 7, 8) ───
    gam_model    <- result$best_model
    cell_weights <- result$cell_weights   # all 1.0 - FOXM1-positive only
    dev_explained <- summary(gam_model)$dev.expl

    # Null model uses same weights as fitted GAM for consistent EP
    null_model   <- gam(Expression ~ 1, data = original_gam_data,
                        weights = cell_weights)
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

    cat("Target TREP cells:", target_cells, "out of",
        nrow(original_gam_data), "total cells\n")

    trep_cells <- rownames(original_gam_data)[
        sorted_indices[seq_len(target_cells)]]

    original_gam_data$is_trep <- ifelse(
        rownames(original_gam_data) %in% trep_cells, "TREP", "non-TREP")
    original_gam_data$timing  <- ifelse(
        original_gam_data$FOXM1 < ip_value, "Pre-IP", "Post-IP")
    original_gam_data$cell_group <- paste0(
        original_gam_data$timing, "_", original_gam_data$is_trep)

    # ── Match cells with Seurat object ────────────────────────────
    common_cells <- intersect(colnames(seurat_obj), rownames(original_gam_data))
    cat("Common cells found:", length(common_cells), "\n")

    if (length(common_cells) == 0) {
        stop("No common cells found between seurat_obj and GAM data for ",
             patient)
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
        FOXM1             = sample_data_matched$FOXM1,
        Expression        = sample_data_matched$Expression,   # DFG composite
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
    cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE,
                       close_loop = FALSE)

    # Root = cell with highest FOXM1 expression
    root_cell <- colnames(cds)[which.max(colData(cds)$FOXM1)]
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
            title    = paste("Monocle3 Trajectory -", patient,
                             "| FOXM1-positive GBM cells"),
            subtitle = paste("TREP (purple/blue) vs non-TREP (gray) |",
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
            title    = paste("Monocle3 Pseudotime -", patient,
                             "| FOXM1-positive GBM cells"),
            subtitle = paste("Lighter = higher pseudotime | Pre/Post IP =",
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
    cat("Patient:", patient, "\n")
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
# SECTION 4: Run trajectory for all patients
# ------------------------------------------------------------------
cat("\n\nGenerating Monocle3 trajectory plots for all patients...\n")
cat(strrep("=", 60), "\n")

all_trajectories <- list()

for (patient in TARGET_PATIENTS) {
    tryCatch({
        cds_result <- create_monocle3_trajectory(patient)
        if (!is.null(cds_result)) {
            all_trajectories[[patient]] <- cds_result
        }
    }, error = function(e) {
        cat("Error creating trajectory for", patient, ":",
            conditionMessage(e), "\n")
    })
}

cat("\nCompleted generating Monocle3 trajectory plots for all patients.\n")


##############################################################################
# Quantitative Analysis of Pre-IP and Post-IP TREP Cells
##############################################################################

analyze_cell_type_clustering <- function(cds, patient_id, colors) {
    cat("\n=== ANALYZING CELL TYPE CLUSTERING FOR", patient_id, "===\n")

    umap_coords   <- reducedDims(cds)$UMAP
    cell_metadata <- colData(cds)

    analysis_data <- data.frame(
        cell_id    = rownames(umap_coords),
        UMAP1      = umap_coords[, 1],
        UMAP2      = umap_coords[, 2],
        cell_group = cell_metadata$cell_group,
        is_trep    = cell_metadata$is_trep,
        timing     = cell_metadata$timing,
        FOXM1      = cell_metadata$FOXM1,
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
                between_distances[[paste(g1, "vs", g2, sep = "_")]] <-
                    as.vector(bd)
            }
        }
    }

    # ── 2. Statistical testing (UMAP1) ────────────────────────────
    cat("Performing statistical tests on UMAP1 distributions...\n")

    statistical_results <- data.frame()

    post_ip_trep_umap1 <- analysis_data$UMAP1[
        analysis_data$cell_group == "Post-IP_TREP"]
    pre_ip_trep_umap1  <- analysis_data$UMAP1[
        analysis_data$cell_group == "Pre-IP_TREP"]

    if (length(post_ip_trep_umap1) > 0 && length(pre_ip_trep_umap1) > 0) {

        post_normal <- ifelse(
            length(post_ip_trep_umap1) >= 3 &&
                length(post_ip_trep_umap1) <= 5000,
            shapiro.test(post_ip_trep_umap1)$p.value > 0.05, FALSE)
        pre_normal <- ifelse(
            length(pre_ip_trep_umap1) >= 3 &&
                length(pre_ip_trep_umap1) <= 5000,
            shapiro.test(pre_ip_trep_umap1)$p.value > 0.05, FALSE)

        both_normal <- post_normal && pre_normal
        test_type   <- ifelse(both_normal, "t-test", "Mann-Whitney")

        if (both_normal) {
            tres        <- t.test(post_ip_trep_umap1, pre_ip_trep_umap1)
            p_value     <- tres$p.value
            effect_size <- abs(mean(post_ip_trep_umap1) -
                                    mean(pre_ip_trep_umap1)) /
                sqrt(((length(post_ip_trep_umap1) - 1) *
                          var(post_ip_trep_umap1) +
                          (length(pre_ip_trep_umap1) - 1) *
                          var(pre_ip_trep_umap1)) /
                         (length(post_ip_trep_umap1) +
                              length(pre_ip_trep_umap1) - 2))
        } else {
            wres        <- wilcox.test(post_ip_trep_umap1, pre_ip_trep_umap1)
            p_value     <- wres$p.value
            n1          <- length(post_ip_trep_umap1)
            n2          <- length(pre_ip_trep_umap1)
            effect_size <- (2 * wres$statistic) / (n1 * n2) - 1
        }

        result_row <- data.frame(
            Sample              = patient_id,
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

        cat("Post-IP_TREP: n =", length(post_ip_trep_umap1),
            ", Normal =", post_normal, "\n")
        cat("Pre-IP_TREP:  n =", length(pre_ip_trep_umap1),
            ", Normal =", pre_normal, "\n")
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
            Sample              = patient_id,
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
        statistical_results$BH_Adjusted_P_Value[
            !is.na(statistical_results$P_Value)] <- adj_p
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
                             aes(x = UMAP1, y = Cell_Group,
                                 fill = Cell_Group)) +
        geom_density_ridges(alpha = 0.5, scale = 1.5,
                            rel_min_height = 0.01) +
        scale_fill_manual(values = colors) +
        scale_y_discrete(
            limits = rev(c("Pre-IP_TREP", "Pre-IP_non-TREP",
                           "Post-IP_TREP", "Post-IP_non-TREP")),
            expand = expansion(mult = c(0, 0.2))) +
        coord_cartesian(clip = "off") +
        labs(
            title    = paste("Cell Type Clustering -", patient_id,
                             "| FOXM1-positive GBM cells"),
            subtitle = "UMAP1 distributions for each cell type",
            x        = "UMAP1",
            y        = "Cell Type",
            caption  = paste("Overall Silhouette Score:",
                             round(avg_silhouette, 3),
                             "| Lower spread = tighter clustering")
        ) +
        theme_ridges() +
        theme(
            plot.title         = element_text(size = 14, face = "bold"),
            plot.subtitle      = element_text(size = 12),
            legend.title       = element_text(size = 12, face = "bold"),
            legend.text        = element_text(size = 10),
            axis.text.y        = element_text(size = 8),
            axis.text.x        = element_text(size = 8,
                                              margin = margin(t = 2)),
            plot.caption       = element_text(size = 9, hjust = 0),
            panel.grid.major.x = element_line(linetype = "dotted",
                                              color = "gray70"),
            panel.grid.minor.x = element_blank(),
            plot.margin        = margin(t = 20, r = 5, b = 5, l = 5,
                                        unit = "pt")
        ) +
        guides(fill = "none")

    print(ridgeline_plot)

    # ── 6. Summary stats table ─────────────────────────────────────
    summary_stats <- statistical_results %>%
        dplyr::filter(!grepl("Silhouette", Test_Type)) %>%
        dplyr::select(Cell_Group, Post_IP_TREP_Mean, Pre_IP_TREP_Mean,
                      P_Value, BH_Adjusted_P_Value, Effect_Size,
                      Test_Type) %>%
        dplyr::mutate(
            Significant_BH     = ifelse(BH_Adjusted_P_Value < 0.05,
                                        "Yes", "No"),
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
# SECTION 5: Process all patients & export
# ------------------------------------------------------------------
process_all_patients_clustering <- function(
        all_trajectories,
        output_file = "monocle3_clustering_FOXM1_DFG_FOXM1-pos.xlsx") {

    cat("\n=== PROCESSING ALL PATIENTS FOR CLUSTERING ANALYSIS ===\n")

    colors <- c(
        "Pre-IP_non-TREP"  = "#CCCCCC",
        "Post-IP_non-TREP" = "#666666",
        "Pre-IP_TREP"      = "#DDA0DD",
        "Post-IP_TREP"     = "#6666FF"
    )

    all_statistical_results <- data.frame()

    for (patient_name in names(all_trajectories)) {
        cat("\nProcessing", patient_name, "...\n")

        tryCatch({
            res <- analyze_cell_type_clustering(
                all_trajectories[[patient_name]], patient_name, colors)
            all_statistical_results <- rbind(all_statistical_results,
                                             res$statistical_results)
        }, error = function(e) {
            cat("Error analyzing clustering for", patient_name, ":",
                conditionMessage(e), "\n")
        })
    }

    # BH correction across all patients
    cat("Applying Benjamini-Hochberg correction across all patients...\n")

    if (nrow(all_statistical_results) > 0) {
        p_rows <- !is.na(all_statistical_results$P_Value)
        if (sum(p_rows) > 0) {
            adj_p <- p.adjust(all_statistical_results$P_Value[p_rows],
                              method = "BH")
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
            "Dataset: Neftel GBM SS2 - FOXM1-positive adult malignant cells",
            "GOI: FOXM1 (x-axis) | Gene set: DFG composite",
            "DFG members: BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1",
            "Patients: MGH152 (dual-pass), MGH110 (partial), MGH66 (dual-pass), MGH100 (partial), MGH121 (partial)",
            "Filter: dual-filter pass (rho>=0.50 AND tau>=0.40) or partial pass, n_pos>=50",
            "IP source: inflection_points[\"patient|DFG\"] from Part 7",
            "TREP classification: top (dev_expl x n) cells by EP = 1 - model_res^2 / null_res^2",
            "Null model uses same weights as fitted GAM (all 1.0 for FOXM1-positive cells)",
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
            "BH_Adjusted_P_Value: Benjamini-Hochberg corrected p-value (across all 5 patients)",
            "Effect_Size: Cohen's d (t-test) or Cliff's delta (Mann-Whitney)",
            "Test_Type: Statistical test used (t-test / Mann-Whitney / Silhouette_Score_UMAP_Separation)"
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
                    Total_Patients          = dplyr::n(),
                    Mean_Effect_Size        = mean(Effect_Size, na.rm = TRUE),
                    Mean_Post_IP_TREP_UMAP1 = mean(Post_IP_TREP_Mean,
                                                    na.rm = TRUE),
                    Mean_Pre_IP_TREP_UMAP1  = mean(Pre_IP_TREP_Mean,
                                                    na.rm = TRUE),
                    Significant_Raw         = sum(P_Value < 0.05,
                                                  na.rm = TRUE),
                    Significant_BH          = sum(BH_Adjusted_P_Value < 0.05,
                                                  na.rm = TRUE),
                    T_Test_Used             = sum(Test_Type == "t-test",
                                                  na.rm = TRUE),
                    Mann_Whitney_Used       = sum(Test_Type == "Mann-Whitney",
                                                  na.rm = TRUE),
                    .groups = "drop"
                )
            print(overall_summary)

            cat("\nRaw vs BH-adjusted significance per patient:\n")
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
cat("\n\nStarting clustering analysis for all patients...\n")
clustering_analysis_results <- process_all_patients_clustering(all_trajectories)

cat("\n========== Part 9 complete ==========\n")
cat("Output: monocle3_clustering_FOXM1_DFG_FOXM1-pos.xlsx\n")


##############################################################################
# CEP-IP Framework - Glioblastoma Analysis
# Part 9: Monocle3 Trajectory
#         (TREP and non-TREP Cells, Quantitative Analysis)
#         FOXM1-positive adult malignant GBM cells
#
# Requires:
#   - NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds  (from Part 5)
#   - inflection_points          (named numeric vector from Part 7,
#                                 keyed by "patient|DFG")
#   - gbm_results$seurat_obj     (Seurat object with scRNA-seq data;
#                                 cell barcodes must match rownames
#                                 in result$data from the GAM models)
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
library(plotly)
library(vegan)    # PERMANOVA (adonis2) and PERMDISP (betadisper)

# ------------------------------------------------------------------
# SECTION 1: Load Part 5 models
# ------------------------------------------------------------------
cat("Loading Part 5 models...\n")
rds_data    <- readRDS("NeftelGBM_SS2_FOXM1_DFG_GAM_models_FOXM1-pos.rds")
all_results <- rds_data$all_results

TARGET_PATIENTS <- c("MGH152", "MGH110", "MGH66", "MGH100", "MGH121", "MGH104")
gene_set_label  <- "DFG"

valid_results <- Filter(Negate(is.null), all_results)
valid_results <- valid_results[
    names(valid_results) %in%
        paste0(TARGET_PATIENTS, "|", gene_set_label)
]

cat("DFG fits retained for analysis:", length(valid_results), "\n")
cat("Keys:", paste(names(valid_results), collapse = ", "), "\n")

seurat_obj <- gbm_results$seurat_obj

# ------------------------------------------------------------------
# SECTION 2: Check prerequisites
# ------------------------------------------------------------------
if (!exists("inflection_points")) {
    stop(paste0(
        "inflection_points not found in environment.\n",
        "Run Part 7 (IP detection) to create inflection_points.\n",
        "Example key: inflection_points[\"MGH152|DFG\"]"
    ))
}

# ------------------------------------------------------------------
# SECTION 3: Monocle3 trajectory per patient
# ------------------------------------------------------------------
create_monocle3_trajectory <- function(patient) {
    key <- paste0(patient, "|", gene_set_label)
    cat("=== CREATING MONOCLE3 TRAJECTORY FOR", patient,
        "(key:", key, ") ===\n")

    result <- valid_results[[key]]
    if (is.null(result)) {
        cat("  Skipping:", patient, "- no valid DFG fit\n")
        return(NULL)
    }

    original_gam_data <- result$data
    if (is.null(original_gam_data))
        stop("GAM data not found for patient ", patient)

    cat("GAM data dimensions:", nrow(original_gam_data), "x",
        ncol(original_gam_data), "\n")

    ip_value <- inflection_points[key]
    if (is.na(ip_value))
        stop("No IP value found for key ", key, ". Run Part 7 first.")
    cat("IP value for", patient, ":", ip_value, "\n")

    # ---- Cell classification (identical to Parts 6, 7, 8) --------
    gam_model     <- result$best_model
    cell_weights  <- result$cell_weights
    dev_explained <- summary(gam_model)$dev.expl

    null_model   <- gam(Expression ~ 1, data = original_gam_data,
                        weights = cell_weights)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(gam_model)

    null_sq_diff  <- (original_gam_data$Expression - null_fitted)^2
    model_sq_diff <- (original_gam_data$Expression - model_fitted)^2

    explanatory_power <- ifelse(
        null_sq_diff > 0,
        1 - (model_sq_diff / null_sq_diff), 0)

    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    target_cells   <- max(1L, min(
        round(nrow(original_gam_data) * dev_explained),
        nrow(original_gam_data)))

    cat("Target TREP cells:", target_cells, "out of",
        nrow(original_gam_data), "\n")

    trep_cells <- rownames(original_gam_data)[
        sorted_indices[seq_len(target_cells)]]

    original_gam_data$is_trep <- ifelse(
        rownames(original_gam_data) %in% trep_cells, "TREP", "non-TREP")
    original_gam_data$timing  <- ifelse(
        original_gam_data$FOXM1 < ip_value, "Pre-IP", "Post-IP")
    original_gam_data$cell_group <- paste0(
        original_gam_data$timing, "_", original_gam_data$is_trep)

    # ---- Match with Seurat ---------------------------------------
    common_cells <- intersect(colnames(seurat_obj),
                              rownames(original_gam_data))
    cat("Common cells found:", length(common_cells), "\n")
    if (length(common_cells) == 0)
        stop("No common cells for patient ", patient)

    sample_seurat       <- subset(seurat_obj, cells = common_cells)
    sample_data_matched <- original_gam_data[common_cells, ]

    cat("Final processing:", nrow(sample_data_matched), "cells\n")
    cat("Cell group distribution:\n")
    print(table(sample_data_matched$cell_group))

    # ---- Build CDS -----------------------------------------------
    cat("Converting to Monocle3 format...\n")
    count_matrix <- GetAssayData(sample_seurat,
                                 slot = "counts", assay = "RNA")

    cell_metadata <- data.frame(
        cell_id           = colnames(sample_seurat),
        FOXM1             = sample_data_matched$FOXM1,
        Expression        = sample_data_matched$Expression,
        cell_group        = sample_data_matched$cell_group,
        is_trep           = sample_data_matched$is_trep,
        timing            = sample_data_matched$timing,
        explanatory_power = explanatory_power[
            match(colnames(sample_seurat),
                  rownames(original_gam_data))],
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

    cat("Preprocessing data...\n")
    cds <- preprocess_cds(cds, num_dim = 200, verbose = FALSE)

    # ---- 2D UMAP -------------------------------------------------
    cat("Reducing dimensions (2D UMAP)...\n")
    cds <- reduce_dimension(cds, max_components = 2L, verbose = FALSE)

    cat("Clustering cells...\n")
    k_clust <- min(250L, ncol(cds) - 2L)
    cds <- cluster_cells(cds, verbose = FALSE,
                         resolution = 0.0005, k = k_clust)

    cat("Learning trajectory...\n")
    cds <- learn_graph(cds, verbose = FALSE,
                       use_partition = FALSE, close_loop = FALSE)

    root_cell <- colnames(cds)[which.max(colData(cds)$FOXM1)]
    cds <- order_cells(cds, root_cells = root_cell, verbose = FALSE)

    # ---- 3D UMAP (max_components = 3L) ---------------------------
    # v2 FIX: correct argument is max_components.
    cat("Reducing dimensions (3D UMAP)...\n")
    cds_3d <- tryCatch({
        reduce_dimension(cds, max_components = 3L, verbose = FALSE)
    }, error = function(e) {
        cat("  monocle3 3D UMAP failed:", conditionMessage(e), "\n")
        cat("  Falling back to uwot::umap...\n")
        if (!requireNamespace("uwot", quietly = TRUE)) {
            cat("  uwot not available - 3D UMAP skipped.\n")
            return(NULL)
        }
        pca_embed  <- reducedDims(cds)$PCA
        umap3d_mat <- uwot::umap(pca_embed, n_components = 3L,
                                 n_neighbors = 15L, min_dist = 0.1,
                                 metric = "cosine")
        colnames(umap3d_mat) <- c("UMAP_1", "UMAP_2", "UMAP_3")
        cds_copy <- cds
        reducedDims(cds_copy)[["UMAP"]] <- umap3d_mat
        cds_copy
    })

    # ---- Colors --------------------------------------------------
    colors <- c(
        "Pre-IP_non-TREP"  = "#CCCCCC",
        "Post-IP_non-TREP" = "#666666",
        "Pre-IP_TREP"      = "#DDA0DD",
        "Post-IP_TREP"     = "#6666FF"
    )

    # ---- 2D trajectory plot --------------------------------------
    cat("Creating trajectory plot...\n")
    trajectory_plot <- plot_cells(cds,
        color_cells_by = "cell_group", label_cell_groups = FALSE,
        label_leaves = FALSE, label_branch_points = FALSE,
        label_roots = FALSE, show_trajectory_graph = TRUE,
        graph_label_size = 3, cell_size = 1.5, cell_stroke = 0.5,
        trajectory_graph_color = "#333333",
        trajectory_graph_segment_size = 1) +
        scale_color_manual(values = colors) +
        labs(
            title    = paste("Monocle3 Trajectory -", patient,
                             "| FOXM1-positive GBM cells"),
            subtitle = paste("TREP (purple/blue) vs non-TREP (gray) |",
                             "Pre/Post IP =", round(ip_value, 4)),
            color    = "Cell Group") +
        theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title  = element_text(size = 12, face = "bold"),
            legend.text   = element_text(size = 10),
            axis.text.x   = element_text(size = 6),
            axis.text.y   = element_text(size = 6))

    # ---- 2D pseudotime plot --------------------------------------
    cat("Creating pseudotime plot...\n")
    pseudotime_plot <- plot_cells(cds,
        color_cells_by = "pseudotime", label_cell_groups = FALSE,
        label_leaves = FALSE, label_branch_points = FALSE,
        label_roots = FALSE, show_trajectory_graph = TRUE,
        graph_label_size = 3, cell_size = 0.5, cell_stroke = 0.5,
        trajectory_graph_color = "#333333",
        trajectory_graph_segment_size = 1) +
        scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
        labs(
            title    = paste("Monocle3 Pseudotime -", patient,
                             "| FOXM1-positive GBM cells"),
            subtitle = paste("Lighter = higher pseudotime | Pre/Post IP =",
                             round(ip_value, 4))) +
        theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title  = element_text(size = 12, face = "bold"),
            legend.text   = element_text(size = 10))

    suppressWarnings(print(trajectory_plot))
    suppressWarnings(print(pseudotime_plot))

    # ---- 3D interactive plots ------------------------------------
    if (!is.null(cds_3d)) {
        cat("Creating 3D UMAP plots...\n")
        tryCatch({
            umap3d    <- reducedDims(cds_3d)$UMAP
            umap3d_df <- data.frame(
                UMAP1      = umap3d[, 1],
                UMAP2      = umap3d[, 2],
                UMAP3      = umap3d[, 3],
                cell_group = cell_metadata$cell_group,
                FOXM1      = cell_metadata$FOXM1,
                pseudotime = tryCatch(pseudotime(cds_3d),
                    error = function(e) rep(NA_real_, ncol(cds_3d))),
                stringsAsFactors = FALSE)

            color_map_3d <- c(
                "Pre-IP_non-TREP"  = "#CCCCCC",
                "Post-IP_non-TREP" = "#666666",
                "Pre-IP_TREP"      = "#DDA0DD",
                "Post-IP_TREP"     = "#6666FF")

            fig_group <- plot_ly(
                data = umap3d_df,
                x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
                color = ~cell_group, colors = color_map_3d,
                type = "scatter3d", mode = "markers",
                marker = list(size = 4, opacity = 0.8),
                text = ~paste0("Group: ", cell_group,
                               "<br>FOXM1: ", round(FOXM1, 3)),
                hoverinfo = "text") %>%
            layout(
                title = paste0("3D UMAP - ", patient,
                               " | FOXM1-positive GBM cells<br>",
                               "IP = ", round(ip_value, 4),
                               " | Rotate to inspect continuous trajectory"),
                scene = list(
                    xaxis = list(title = "UMAP1"),
                    yaxis = list(title = "UMAP2"),
                    zaxis = list(title = "UMAP3")))
            print(fig_group)

            if (!all(is.na(umap3d_df$pseudotime))) {
                fig_pseudo <- plot_ly(
                    data = umap3d_df,
                    x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
                    color = ~pseudotime,
                    colors = viridis::viridis(100, option = "plasma"),
                    type = "scatter3d", mode = "markers",
                    marker = list(size = 4, opacity = 0.8),
                    text = ~paste0("Pseudotime: ",
                                   round(pseudotime, 3),
                                   "<br>FOXM1: ", round(FOXM1, 3)),
                    hoverinfo = "text") %>%
                layout(
                    title = paste0("3D UMAP Pseudotime - ", patient,
                                   " | FOXM1-positive GBM cells"),
                    scene = list(
                        xaxis = list(title = "UMAP1"),
                        yaxis = list(title = "UMAP2"),
                        zaxis = list(title = "UMAP3")))
                print(fig_pseudo)
            } else {
                cat("  Pseudotime unavailable for cds_3d (uwot fallback)",
                    "- skipping pseudotime 3D plot.\n")
            }

            cat("  3D UMAP rendered in RStudio Viewer pane.\n")
        }, error = function(e) {
            cat("  3D UMAP plot failed:", conditionMessage(e), "\n")
        })
    } else {
        cat("  3D UMAP skipped (cds_3d is NULL).\n")
    }

    # ---- Summary -------------------------------------------------
    cat("\n=== TRAJECTORY SUMMARY ===\n")
    cat("Patient:", patient, "| Cells:", ncol(cds), "\n")
    cat("IP:", ip_value, "| Dev.expl:",
        round(dev_explained * 100, 1), "%\n")
    cat("Root cell:", root_cell, "\n")
    group_counts <- table(cell_metadata$cell_group)
    for (g in names(group_counts))
        cat(" ", g, ":", group_counts[g],
            "(", round(group_counts[g] / sum(group_counts) * 100, 1),
            "%)\n")

    return(list(cds_2d = cds, cds_3d = cds_3d))
}

# ------------------------------------------------------------------
# SECTION 4: Run trajectory for all patients
# ------------------------------------------------------------------
cat("\n\nGenerating Monocle3 trajectory plots...\n")
cat(strrep("=", 60), "\n")

all_trajectories    <- list()   # 2D CDS objects
all_trajectories_3d <- list()   # 3D CDS objects (NEW v3)

for (patient in TARGET_PATIENTS) {
    tryCatch({
        cds_result <- create_monocle3_trajectory(patient)
        if (!is.null(cds_result)) {
            all_trajectories[[patient]]    <- cds_result$cds_2d
            all_trajectories_3d[[patient]] <- cds_result$cds_3d
        }
    }, error = function(e) {
        cat("Error creating trajectory for", patient, ":",
            conditionMessage(e), "\n")
    })
}

cat("\nCompleted trajectory generation.\n")


##############################################################################
# HELPER: single-axis Shapiro -> t-test / Mann-Whitney
# Returns a one-row data.frame or NULL.
##############################################################################
test_one_axis <- function(post_vals, pre_vals, axis_label) {
    if (length(post_vals) == 0 || length(pre_vals) == 0) return(NULL)

    post_normal <- tryCatch(
        length(post_vals) >= 3 && length(post_vals) <= 5000 &&
            shapiro.test(post_vals)$p.value > 0.05,
        error = function(e) FALSE)
    pre_normal <- tryCatch(
        length(pre_vals) >= 3 && length(pre_vals) <= 5000 &&
            shapiro.test(pre_vals)$p.value > 0.05,
        error = function(e) FALSE)

    both_normal <- post_normal && pre_normal
    test_type   <- ifelse(both_normal, "t-test", "Mann-Whitney")

    if (both_normal) {
        tres        <- t.test(post_vals, pre_vals)
        p_value     <- tres$p.value
        effect_size <- abs(mean(post_vals) - mean(pre_vals)) /
            sqrt(((length(post_vals) - 1) * var(post_vals) +
                      (length(pre_vals) - 1) * var(pre_vals)) /
                     (length(post_vals) + length(pre_vals) - 2))
    } else {
        wres        <- wilcox.test(post_vals, pre_vals)
        p_value     <- wres$p.value
        n1 <- length(post_vals); n2 <- length(pre_vals)
        effect_size <- (2 * wres$statistic) / (n1 * n2) - 1
    }

    data.frame(
        Axis                = axis_label,
        Test_Type           = test_type,
        Post_IP_TREP_N      = length(post_vals),
        Post_IP_TREP_Mean   = mean(post_vals),
        Post_IP_TREP_Median = median(post_vals),
        Post_IP_TREP_IQR    = IQR(post_vals),
        Post_IP_TREP_Normal = post_normal,
        Pre_IP_TREP_N       = length(pre_vals),
        Pre_IP_TREP_Mean    = mean(pre_vals),
        Pre_IP_TREP_Median  = median(pre_vals),
        Pre_IP_TREP_IQR     = IQR(pre_vals),
        Pre_IP_TREP_Normal  = pre_normal,
        P_Value             = p_value,
        Effect_Size         = effect_size,
        stringsAsFactors    = FALSE
    )
}


##############################################################################
# 2D CLUSTERING ANALYSIS
##############################################################################
analyze_cell_type_clustering <- function(cds, patient_id, colors) {
    cat("\n=== ANALYZING 2D CLUSTERING FOR", patient_id, "===\n")

    umap_coords   <- reducedDims(cds)$UMAP
    cell_metadata <- colData(cds)

    analysis_data <- data.frame(
        cell_id    = rownames(umap_coords),
        UMAP1      = umap_coords[, 1],
        UMAP2      = umap_coords[, 2],
        cell_group = cell_metadata$cell_group,
        is_trep    = cell_metadata$is_trep,
        timing     = cell_metadata$timing,
        FOXM1      = cell_metadata$FOXM1,
        Expression = cell_metadata$Expression,
        stringsAsFactors = FALSE
    )

    # UMAP1 test
    cat("Performing UMAP1 statistical test...\n")
    post_vals <- analysis_data$UMAP1[
        analysis_data$cell_group == "Post-IP_TREP"]
    pre_vals  <- analysis_data$UMAP1[
        analysis_data$cell_group == "Pre-IP_TREP"]

    statistical_results <- data.frame()

    if (length(post_vals) > 0 && length(pre_vals) > 0) {
        row <- test_one_axis(post_vals, pre_vals, "UMAP1")
        if (!is.null(row)) {
            result_row <- data.frame(
                Sample              = patient_id,
                Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
                Test_Type           = row$Test_Type,
                Post_IP_TREP_Median = row$Post_IP_TREP_Median,
                Post_IP_TREP_IQR    = row$Post_IP_TREP_IQR,
                Post_IP_TREP_Mean   = row$Post_IP_TREP_Mean,
                Post_IP_TREP_N      = row$Post_IP_TREP_N,
                Post_IP_TREP_Normal = row$Post_IP_TREP_Normal,
                Pre_IP_TREP_Median  = row$Pre_IP_TREP_Median,
                Pre_IP_TREP_IQR     = row$Pre_IP_TREP_IQR,
                Pre_IP_TREP_Mean    = row$Pre_IP_TREP_Mean,
                Pre_IP_TREP_N       = row$Pre_IP_TREP_N,
                Pre_IP_TREP_Normal  = row$Pre_IP_TREP_Normal,
                P_Value             = row$P_Value,
                Effect_Size         = row$Effect_Size,
                stringsAsFactors    = FALSE)
            statistical_results <- rbind(statistical_results,
                                         result_row)
            cat("Post-IP_TREP: n =", row$Post_IP_TREP_N,
                ", Normal =", row$Post_IP_TREP_Normal, "\n")
            cat("Pre-IP_TREP:  n =", row$Pre_IP_TREP_N,
                ", Normal =", row$Pre_IP_TREP_Normal, "\n")
            cat("Test:", row$Test_Type,
                "| p =", row$P_Value, "\n")
        }
    }

    # 2D Silhouette
    cat("Calculating 2D silhouette score...\n")
    dist_matrix   <- dist(umap_coords)
    binary_labels <- ifelse(analysis_data$cell_group == "Post-IP_TREP",
                            1,
                     ifelse(analysis_data$cell_group == "Pre-IP_TREP",
                            2, NA))
    valid_idx    <- !is.na(binary_labels)
    valid_labels <- binary_labels[valid_idx]
    valid_coords <- umap_coords[valid_idx, ]
    avg_silhouette <- NA

    if (length(unique(valid_labels)) == 2 && nrow(valid_coords) > 2) {
        sil_scores     <- silhouette(valid_labels, dist(valid_coords))
        avg_silhouette <- mean(sil_scores[, 3])
        sil_row <- data.frame(
            Sample              = patient_id,
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
            stringsAsFactors    = FALSE)
        statistical_results <- rbind(statistical_results, sil_row)
    }

    # Within-patient BH
    p_rows <- !is.na(statistical_results$P_Value)
    if (sum(p_rows) > 0) {
        statistical_results$BH_Adjusted_P_Value <- NA
        statistical_results$BH_Adjusted_P_Value[p_rows] <-
            p.adjust(statistical_results$P_Value[p_rows], method = "BH")
    }

    # 2D UMAP1 ridgeline
    cat("Preparing 2D ridgeline plot (UMAP1)...\n")
    ridgeline_data <- data.frame(
        UMAP1      = analysis_data$UMAP1,
        Cell_Group = factor(analysis_data$cell_group,
                            levels = c("Pre-IP_TREP", "Pre-IP_non-TREP",
                                       "Post-IP_TREP",
                                       "Post-IP_non-TREP")),
        stringsAsFactors = FALSE)

    ridgeline_plot <- ggplot(ridgeline_data,
                             aes(x = UMAP1, y = Cell_Group,
                                 fill = Cell_Group)) +
        geom_density_ridges(alpha = 0.5, scale = 1.5,
                            rel_min_height = 0.01) +
        scale_fill_manual(values = colors) +
        scale_y_discrete(
            limits = rev(c("Pre-IP_TREP", "Pre-IP_non-TREP",
                           "Post-IP_TREP", "Post-IP_non-TREP")),
            expand = expansion(mult = c(0, 0.2))) +
        coord_cartesian(clip = "off") +
        labs(
            title    = paste("Cell Type Clustering -", patient_id,
                             "| FOXM1-positive GBM cells"),
            subtitle = "UMAP1 distributions per cell group (2D)",
            x = "UMAP1", y = "Cell Type",
            caption  = paste("2D Silhouette Score:",
                             round(avg_silhouette, 3))) +
        theme_ridges() +
        theme(
            plot.title         = element_text(size = 14, face = "bold"),
            plot.subtitle      = element_text(size = 12),
            legend.title       = element_text(size = 12, face = "bold"),
            legend.text        = element_text(size = 10),
            axis.text.y        = element_text(size = 8),
            axis.text.x        = element_text(size = 8,
                                              margin = margin(t = 2)),
            plot.caption       = element_text(size = 9, hjust = 0),
            panel.grid.major.x = element_line(linetype = "dotted",
                                              color = "gray70"),
            panel.grid.minor.x = element_blank(),
            plot.margin        = margin(20, 5, 5, 5, "pt")) +
        guides(fill = "none")

    print(ridgeline_plot)

    summary_stats <- statistical_results %>%
        dplyr::filter(!grepl("Silhouette", Test_Type)) %>%
        dplyr::select(Cell_Group, Post_IP_TREP_Mean, Pre_IP_TREP_Mean,
                      P_Value, BH_Adjusted_P_Value, Effect_Size,
                      Test_Type) %>%
        dplyr::mutate(
            Significant_BH = ifelse(BH_Adjusted_P_Value < 0.05,
                                    "Yes", "No"),
            Clustering_Quality = dplyr::case_when(
                abs(Effect_Size) > 0.5 ~ "Strong",
                abs(Effect_Size) > 0.3 ~ "Moderate",
                abs(Effect_Size) > 0.1 ~ "Weak",
                TRUE                   ~ "Very Weak"))

    cat("\n2D Clustering Summary:\n")
    print(summary_stats)

    return(list(
        statistical_results = statistical_results,
        ridgeline_data      = ridgeline_data,
        ridgeline_plot      = ridgeline_plot,
        summary_stats       = summary_stats,
        overall_silhouette  = avg_silhouette
    ))
}


##############################################################################
# 3D CLUSTERING ANALYSIS
##############################################################################
analyze_cell_type_clustering_3d <- function(cds_3d, patient_id,
                                             colors) {
    if (is.null(cds_3d)) {
        cat("  Skipping 3D analysis for", patient_id,
            "- cds_3d is NULL.\n")
        return(NULL)
    }

    cat("\n=== ANALYZING 3D CLUSTERING FOR", patient_id, "===\n")

    umap3d        <- reducedDims(cds_3d)$UMAP
    cell_metadata <- colData(cds_3d)

    if (ncol(umap3d) < 3) {
        cat("  3D UMAP has fewer than 3 columns for", patient_id,
            "- skipping.\n")
        return(NULL)
    }

    analysis_data <- data.frame(
        cell_id    = rownames(umap3d),
        UMAP1      = umap3d[, 1],
        UMAP2      = umap3d[, 2],
        UMAP3      = umap3d[, 3],
        cell_group = as.character(cell_metadata$cell_group),
        is_trep    = as.character(cell_metadata$is_trep),
        timing     = as.character(cell_metadata$timing),
        FOXM1      = cell_metadata$FOXM1,
        stringsAsFactors = FALSE
    )

    post_idx <- analysis_data$cell_group == "Post-IP_TREP"
    pre_idx  <- analysis_data$cell_group == "Pre-IP_TREP"
    cat("Post-IP_TREP:", sum(post_idx), "| Pre-IP_TREP:",
        sum(pre_idx), "\n")

    stat_rows <- list()

    # ---- 1. Per-axis tests (UMAP1, UMAP2, UMAP3) -----------------
    # UMAP3 is the key axis for loop/helix structures: if Post-IP TREP
    # separates from Pre-IP TREP only in 3D (as in MGH152), UMAP3 will
    # show a significant result even when UMAP1 and UMAP2 do not.
    cat("Running per-axis tests (UMAP1, UMAP2, UMAP3)...\n")
    for (ax in c("UMAP1", "UMAP2", "UMAP3")) {
        row <- test_one_axis(
            analysis_data[[ax]][post_idx],
            analysis_data[[ax]][pre_idx],
            ax)
        if (!is.null(row)) {
            row$Sample     <- patient_id
            row$Cell_Group <- "Post-IP_TREP_vs_Pre-IP_TREP"
            stat_rows[[ax]] <- row
        }
    }

    # ---- 2. 3D centroid-distance test ----------------------------
    # Rotation-invariant spatial separation test.
    # Each group's cells are measured by their distance to their own
    # group centroid in 3D UMAP space.  This captures helix/loop
    # structures regardless of which single axis carries the signal.
    # Effect_Size is set to the Euclidean centroid separation distance
    # (more interpretable than Cliff's delta for this test).
    cat("Running 3D centroid-distance test...\n")
    tryCatch({
        coords_3d     <- as.matrix(
            analysis_data[, c("UMAP1", "UMAP2", "UMAP3")])
        post_cells_3d <- coords_3d[post_idx, , drop = FALSE]
        pre_cells_3d  <- coords_3d[pre_idx,  , drop = FALSE]

        if (nrow(post_cells_3d) > 0 && nrow(pre_cells_3d) > 0) {
            post_centroid <- colMeans(post_cells_3d)
            pre_centroid  <- colMeans(pre_cells_3d)

            dist_post_to_post <- sqrt(rowSums(
                sweep(post_cells_3d, 2, post_centroid)^2))
            dist_pre_to_pre   <- sqrt(rowSums(
                sweep(pre_cells_3d,  2, pre_centroid)^2))

            centroid_sep <- sqrt(sum(
                (post_centroid - pre_centroid)^2))

            row_cd <- test_one_axis(dist_post_to_post,
                                    dist_pre_to_pre,
                                    "3D_Centroid_Distance")
            if (!is.null(row_cd)) {
                row_cd$Sample     <- patient_id
                row_cd$Cell_Group <- "Post-IP_TREP_vs_Pre-IP_TREP"
                row_cd$Effect_Size <- centroid_sep
                row_cd$Test_Type   <- paste0(
                    row_cd$Test_Type,
                    " [CentroidSep=", round(centroid_sep, 4), "]")
                stat_rows[["3D_Centroid_Distance"]] <- row_cd
            }
            cat("  Post centroid:", round(post_centroid, 3), "\n")
            cat("  Pre centroid: ", round(pre_centroid,  3), "\n")
            cat("  3D centroid separation:", round(centroid_sep, 4),
                "\n")
        }
    }, error = function(e) {
        cat("  3D centroid test failed:", conditionMessage(e), "\n")
    })

    # ---- 3. 3D silhouette ----------------------------------------
    cat("Calculating 3D silhouette score...\n")
    tryCatch({
        binary_labels <- ifelse(post_idx, 1,
                         ifelse(pre_idx,  2, NA))
        valid_idx    <- !is.na(binary_labels)
        valid_labels <- binary_labels[valid_idx]
        valid_coords <- umap3d[valid_idx, ]
        avg_sil_3d   <- NA

        if (length(unique(valid_labels)) == 2 &&
                nrow(valid_coords) > 2) {
            sil_scores <- silhouette(valid_labels, dist(valid_coords))
            avg_sil_3d <- mean(sil_scores[, 3])
            cat("  3D Silhouette:", round(avg_sil_3d, 4), "\n")

            sil_row <- data.frame(
                Sample              = patient_id,
                Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
                Axis                = "3D_Silhouette",
                Test_Type           = "Silhouette_Score_3D_UMAP",
                Post_IP_TREP_N      = sum(post_idx),
                Post_IP_TREP_Mean   = NA, Post_IP_TREP_Median = NA,
                Post_IP_TREP_IQR    = NA, Post_IP_TREP_Normal = NA,
                Pre_IP_TREP_N       = sum(pre_idx),
                Pre_IP_TREP_Mean    = NA, Pre_IP_TREP_Median = NA,
                Pre_IP_TREP_IQR     = NA, Pre_IP_TREP_Normal = NA,
                P_Value             = NA,
                Effect_Size         = avg_sil_3d,
                stringsAsFactors    = FALSE)
            stat_rows[["3D_Silhouette"]] <- sil_row
        }
    }, error = function(e) {
        cat("  3D silhouette failed:", conditionMessage(e), "\n")
    })

    # ---- 4. Hotelling's T² (multivariate, all 3 axes jointly) -----
    cat("Running Hotelling's T2 test (joint UMAP1+2+3)...\n")
    tryCatch({
        coords_3d     <- as.matrix(
            analysis_data[, c("UMAP1", "UMAP2", "UMAP3")])
        post_cells_3d <- coords_3d[post_idx, , drop = FALSE]
        pre_cells_3d  <- coords_3d[pre_idx,  , drop = FALSE]
        n1 <- nrow(post_cells_3d); n2 <- nrow(pre_cells_3d); p <- 3L

        if (n1 >= p + 1L && n2 >= p + 1L) {
            mu1 <- colMeans(post_cells_3d)
            mu2 <- colMeans(pre_cells_3d)
            S1  <- cov(post_cells_3d)
            S2  <- cov(pre_cells_3d)
            # Pooled within-group covariance
            S_pool <- ((n1 - 1) * S1 + (n2 - 1) * S2) /
                (n1 + n2 - 2)
            diff   <- mu1 - mu2
            # T2 statistic
            T2 <- (n1 * n2) / (n1 + n2) *
                t(diff) %*% solve(S_pool) %*% diff
            T2 <- as.numeric(T2)
            # F approximation: F = T2 * (n1+n2-p-1) / (p*(n1+n2-2))
            df1 <- p; df2 <- n1 + n2 - p - 1
            F_stat <- T2 * df2 / (p * (n1 + n2 - 2))
            p_hot  <- pf(F_stat, df1, df2, lower.tail = FALSE)
            # Mahalanobis D² as standardized effect size
            mahal_d2 <- as.numeric(
                t(diff) %*% solve(S_pool) %*% diff)

            hot_row <- data.frame(
                Axis                = "Hotelling_T2_3D",
                Test_Type           = paste0("Hotelling_T2 F(",
                                             df1, ",", df2, ")=",
                                             round(F_stat, 3)),
                Post_IP_TREP_N      = n1,
                Post_IP_TREP_Mean   = NA, Post_IP_TREP_Median = NA,
                Post_IP_TREP_IQR    = NA, Post_IP_TREP_Normal = NA,
                Pre_IP_TREP_N       = n2,
                Pre_IP_TREP_Mean    = NA, Pre_IP_TREP_Median = NA,
                Pre_IP_TREP_IQR     = NA, Pre_IP_TREP_Normal = NA,
                P_Value             = p_hot,
                Effect_Size         = mahal_d2,
                Sample              = patient_id,
                Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
                stringsAsFactors    = FALSE)
            stat_rows[["Hotelling_T2"]] <- hot_row
            cat("  Hotelling T2:", round(T2, 3),
                "| F =", round(F_stat, 3),
                "| p =", round(p_hot, 6),
                "| Mahal D2 =", round(mahal_d2, 3), "\n")
        } else {
            cat("  Hotelling T2 skipped: insufficient cells",
                "(need n >= p+1 = 4 per group).\n")
        }
    }, error = function(e) {
        cat("  Hotelling T2 failed:", conditionMessage(e), "\n")
    })

    # ---- 5. Nearest-Neighbour Purity (NNP) -----------------------
    cat("Running Nearest-Neighbour Purity (k=", k_nnp, ")...\n",
        sep = "")
    tryCatch({
        coords_3d  <- as.matrix(
            analysis_data[, c("UMAP1", "UMAP2", "UMAP3")])
        trep_idx   <- post_idx | pre_idx
        trep_coords <- coords_3d[trep_idx, , drop = FALSE]
        trep_labels <- ifelse(
            analysis_data$cell_group[trep_idx] == "Post-IP_TREP",
            "Post", "Pre")

        n_trep <- nrow(trep_coords)
        k_use  <- min(k_nnp, n_trep - 1L)

        if (n_trep >= k_use + 1L && k_use >= 1L) {
            nn_idx <- FNN::knnx.index(trep_coords, trep_coords,
                                      k = k_use + 1L)
            # Drop first column (self)
            nn_idx <- nn_idx[, -1, drop = FALSE]

            purity <- vapply(seq_len(n_trep), function(i) {
                neigh_labels <- trep_labels[nn_idx[i, ]]
                mean(neigh_labels == trep_labels[i])
            }, numeric(1))

            post_purity <- purity[trep_labels == "Post"]
            pre_purity  <- purity[trep_labels == "Pre"]

            # Compare purity distributions: Post vs Pre
            nnp_row <- test_one_axis(post_purity, pre_purity,
                                     "NNP_3D")
            if (!is.null(nnp_row)) {
                nnp_row$Sample     <- patient_id
                nnp_row$Cell_Group <- "Post-IP_TREP_vs_Pre-IP_TREP"
                nnp_row$Test_Type  <- paste0(nnp_row$Test_Type,
                                             " [k=", k_use, "]")
                stat_rows[["NNP_3D"]] <- nnp_row
            }

            # One-sample t-test vs null purity = 0.5
            post_vs_null <- t.test(post_purity, mu = 0.5)
            pre_vs_null  <- t.test(pre_purity,  mu = 0.5)
            cat("  Post-IP_TREP mean purity:", round(mean(post_purity), 3),
                "| vs-null p =", round(post_vs_null$p.value, 4), "\n")
            cat("  Pre-IP_TREP  mean purity:", round(mean(pre_purity),  3),
                "| vs-null p =", round(pre_vs_null$p.value,  4), "\n")
            cat("  Post vs Pre purity p =",
                round(nnp_row$P_Value, 6), "\n")
        } else {
            cat("  NNP skipped: insufficient TREP cells",
                "(n_trep =", n_trep, ").\n")
        }
    }, error = function(e) {
        cat("  NNP test failed:", conditionMessage(e), "\n")
    })

    # ---- 6. PERMANOVA + PERMDISP (vegan) -------------------------
    cat("Running PERMANOVA + PERMDISP...\n")
    tryCatch({
        if (!requireNamespace("vegan", quietly = TRUE))
            stop("vegan package required: install.packages('vegan')")

        coords_3d   <- as.matrix(
            analysis_data[, c("UMAP1", "UMAP2", "UMAP3")])
        trep_idx    <- post_idx | pre_idx
        trep_coords <- coords_3d[trep_idx, , drop = FALSE]
        trep_groups <- analysis_data$cell_group[trep_idx]

        n_trep <- nrow(trep_coords)
        n_perm <- 999L

        if (n_trep >= 4L &&
                length(unique(trep_groups)) == 2L) {

            euc_dist <- dist(trep_coords, method = "euclidean")
            grp_fac  <- factor(trep_groups)

            # PERMANOVA
            set.seed(42L)
            perm_res <- vegan::adonis2(
                euc_dist ~ grp_fac,
                permutations = n_perm,
                method       = "euclidean")

            perm_p  <- perm_res$`Pr(>F)`[1]
            perm_F  <- perm_res$F[1]
            perm_R2 <- perm_res$R2[1]

            perm_row <- data.frame(
                Axis                = "PERMANOVA_3D",
                Test_Type           = paste0(
                    "PERMANOVA_adonis2 F=", round(perm_F, 3),
                    " R2=", round(perm_R2, 4),
                    " nperm=", n_perm),
                Post_IP_TREP_N      = sum(post_idx),
                Post_IP_TREP_Mean   = NA, Post_IP_TREP_Median = NA,
                Post_IP_TREP_IQR    = NA, Post_IP_TREP_Normal = NA,
                Pre_IP_TREP_N       = sum(pre_idx),
                Pre_IP_TREP_Mean    = NA, Pre_IP_TREP_Median = NA,
                Pre_IP_TREP_IQR     = NA, Pre_IP_TREP_Normal = NA,
                P_Value             = perm_p,
                Effect_Size         = perm_R2,
                Sample              = patient_id,
                Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
                stringsAsFactors    = FALSE)
            stat_rows[["PERMANOVA_3D"]] <- perm_row

            cat("  PERMANOVA: F =", round(perm_F, 3),
                "| R2 =", round(perm_R2, 4),
                "| p =", perm_p, "\n")

            # PERMDISP
            set.seed(42L)
            disp_res  <- vegan::betadisper(euc_dist, grp_fac,
                                           type = "centroid")
            disp_perm <- vegan::permutest(disp_res,
                                          permutations = n_perm)

            disp_p <- disp_perm$tab$`Pr(>F)`[1]
            disp_F <- disp_perm$tab$F[1]

            # Average distance-to-centroid per group
            disp_means     <- tapply(disp_res$distances,
                                     disp_res$group, mean)
            post_disp_mean <- disp_means[
                grep("Post-IP_TREP", names(disp_means))]
            pre_disp_mean  <- disp_means[
                grep("Pre-IP_TREP",  names(disp_means))]

            # Dispersion ratio: Post / Pre (>1 = Post more dispersed)
            disp_ratio <- if (!is.null(post_disp_mean) &&
                                   !is.null(pre_disp_mean) &&
                                   pre_disp_mean > 0)
                as.numeric(post_disp_mean / pre_disp_mean) else NA

            disp_row <- data.frame(
                Axis                = "PERMDISP_3D",
                Test_Type           = paste0(
                    "PERMDISP_betadisper F=", round(disp_F, 3),
                    " PostDisp=", round(post_disp_mean, 3),
                    " PreDisp=",  round(pre_disp_mean,  3),
                    " nperm=", n_perm),
                Post_IP_TREP_N      = sum(post_idx),
                Post_IP_TREP_Mean   = as.numeric(post_disp_mean),
                Post_IP_TREP_Median = NA, Post_IP_TREP_IQR = NA,
                Post_IP_TREP_Normal = NA,
                Pre_IP_TREP_N       = sum(pre_idx),
                Pre_IP_TREP_Mean    = as.numeric(pre_disp_mean),
                Pre_IP_TREP_Median  = NA, Pre_IP_TREP_IQR  = NA,
                Pre_IP_TREP_Normal  = NA,
                P_Value             = disp_p,
                # Effect_Size = dispersion ratio (Post/Pre)
                # >1: Post-IP TREP more dispersed (spread out)
                # <1: Pre-IP TREP more dispersed
                # ~1: homogeneous dispersions (clean PERMANOVA result)
                Effect_Size         = disp_ratio,
                Sample              = patient_id,
                Cell_Group          = "Post-IP_TREP_vs_Pre-IP_TREP",
                stringsAsFactors    = FALSE)
            stat_rows[["PERMDISP_3D"]] <- disp_row

            cat("  PERMDISP: F =", round(disp_F, 3),
                "| p =", disp_p,
                "| DispRatio(Post/Pre) =", round(disp_ratio, 3), "\n")

            # Interpretation guidance
            cat("  Interpretation: ")
            if (perm_p < 0.05 && disp_p >= 0.05) {
                cat("PERMANOVA sig + PERMDISP not sig",
                    "-> clear centroid separation.\n")
            } else if (perm_p < 0.05 && disp_p < 0.05) {
                cat("BOTH significant -> location AND spread differ;",
                    "interpret PERMANOVA cautiously.\n")
            } else {
                cat("PERMANOVA not significant",
                    "-> no detectable 3D separation.\n")
            }

        } else {
            cat("  PERMANOVA/PERMDISP skipped: need >= 4 TREP cells",
                "with 2 groups (n_trep =", n_trep, ").\n")
        }
    }, error = function(e) {
        cat("  PERMANOVA/PERMDISP failed:", conditionMessage(e), "\n")
    })

    # ---- Assemble results ----------------------------------------
    if (length(stat_rows) == 0) {
        cat("  No 3D stats generated for", patient_id, "\n")
        return(NULL)
    }
    all_cols  <- unique(unlist(lapply(stat_rows, names)))
    stat_rows <- lapply(stat_rows, function(r) {
        missing <- setdiff(all_cols, names(r))
        if (length(missing) > 0) r[missing] <- NA
        r[all_cols]
    })
    stat_df <- do.call(rbind, stat_rows)
    rownames(stat_df) <- NULL

    # Within-patient BH (excluding silhouette rows)
    p_rows <- !is.na(stat_df$P_Value)
    if (sum(p_rows) > 0) {
        stat_df$BH_Adjusted_P_Value <- NA
        stat_df$BH_Adjusted_P_Value[p_rows] <-
            p.adjust(stat_df$P_Value[p_rows], method = "BH")
    }

    cat("\n3D Stats for", patient_id, ":\n")
    print(stat_df %>%
              dplyr::select(Axis, Test_Type, P_Value,
                            BH_Adjusted_P_Value, Effect_Size))

    # ---- UMAP3 ridgeline -----------------------------------------
    cat("Creating UMAP3 ridgeline...\n")
    tryCatch({
        umap3_row <- stat_df[stat_df$Axis == "UMAP3" &
                                 !is.na(stat_df$Axis), ]
        umap3_cap <- if (nrow(umap3_row) > 0) {
            paste0("UMAP3: p=",
                   formatC(umap3_row$P_Value[1],
                            format = "e", digits = 2),
                   " (BH=",
                   formatC(umap3_row$BH_Adjusted_P_Value[1],
                            format = "e", digits = 2),
                   ")  Effect=",
                   round(umap3_row$Effect_Size[1], 3))
        } else "UMAP3 test not available"

        ridge3_data <- data.frame(
            UMAP3      = analysis_data$UMAP3,
            Cell_Group = factor(
                analysis_data$cell_group,
                levels = c("Pre-IP_TREP", "Pre-IP_non-TREP",
                           "Post-IP_TREP", "Post-IP_non-TREP")),
            stringsAsFactors = FALSE)

        ridge3_plot <- ggplot(ridge3_data,
                              aes(x = UMAP3, y = Cell_Group,
                                  fill = Cell_Group)) +
            geom_density_ridges(alpha = 0.5, scale = 1.5,
                                rel_min_height = 0.01) +
            scale_fill_manual(values = colors) +
            scale_y_discrete(
                limits = rev(c("Pre-IP_TREP", "Pre-IP_non-TREP",
                               "Post-IP_TREP", "Post-IP_non-TREP")),
                expand = expansion(mult = c(0, 0.2))) +
            coord_cartesian(clip = "off") +
            labs(
                title    = paste("3D UMAP3 Distribution -",
                                 patient_id,
                                 "| FOXM1-positive GBM cells"),
                subtitle = paste("UMAP3 axis distributions -",
                                 "key separator for loop/helix",
                                 "structures"),
                x = "UMAP3", y = "Cell Type",
                caption = umap3_cap) +
            theme_ridges() +
            theme(
                plot.title         = element_text(size = 14,
                                                  face = "bold"),
                plot.subtitle      = element_text(size = 12),
                axis.text.y        = element_text(size = 8),
                axis.text.x        = element_text(
                    size = 8, margin = margin(t = 2)),
                plot.caption       = element_text(size = 9, hjust = 0),
                panel.grid.major.x = element_line(linetype = "dotted",
                                                  color = "gray70"),
                panel.grid.minor.x = element_blank(),
                plot.margin        = margin(20, 5, 5, 5, "pt")) +
            guides(fill = "none")

        print(ridge3_plot)
    }, error = function(e) {
        cat("  UMAP3 ridgeline failed:", conditionMessage(e), "\n")
    })

    # ---- UMAP2 ridgeline -----------------------------------------
    cat("Creating UMAP2 ridgeline...\n")
    tryCatch({
        umap2_row <- stat_df[stat_df$Axis == "UMAP2" &
                                 !is.na(stat_df$Axis), ]
        umap2_cap <- if (nrow(umap2_row) > 0) {
            paste0("UMAP2: p=",
                   formatC(umap2_row$P_Value[1],
                            format = "e", digits = 2),
                   " (BH=",
                   formatC(umap2_row$BH_Adjusted_P_Value[1],
                            format = "e", digits = 2),
                   ")  Effect=",
                   round(umap2_row$Effect_Size[1], 3))
        } else "UMAP2 test not available"

        ridge2_data <- data.frame(
            UMAP2      = analysis_data$UMAP2,
            Cell_Group = factor(
                analysis_data$cell_group,
                levels = c("Pre-IP_TREP", "Pre-IP_non-TREP",
                           "Post-IP_TREP", "Post-IP_non-TREP")),
            stringsAsFactors = FALSE)

        ridge2_plot <- ggplot(ridge2_data,
                              aes(x = UMAP2, y = Cell_Group,
                                  fill = Cell_Group)) +
            geom_density_ridges(alpha = 0.5, scale = 1.5,
                                rel_min_height = 0.01) +
            scale_fill_manual(values = colors) +
            scale_y_discrete(
                limits = rev(c("Pre-IP_TREP", "Pre-IP_non-TREP",
                               "Post-IP_TREP", "Post-IP_non-TREP")),
                expand = expansion(mult = c(0, 0.2))) +
            coord_cartesian(clip = "off") +
            labs(
                title    = paste("3D UMAP2 Distribution -",
                                 patient_id,
                                 "| FOXM1-positive GBM cells"),
                subtitle = "UMAP2 axis distributions from 3D embedding",
                x = "UMAP2", y = "Cell Type",
                caption = umap2_cap) +
            theme_ridges() +
            theme(
                plot.title         = element_text(size = 14,
                                                  face = "bold"),
                plot.subtitle      = element_text(size = 12),
                axis.text.y        = element_text(size = 8),
                axis.text.x        = element_text(
                    size = 8, margin = margin(t = 2)),
                plot.caption       = element_text(size = 9, hjust = 0),
                panel.grid.major.x = element_line(linetype = "dotted",
                                                  color = "gray70"),
                panel.grid.minor.x = element_blank(),
                plot.margin        = margin(20, 5, 5, 5, "pt")) +
            guides(fill = "none")

        print(ridge2_plot)
    }, error = function(e) {
        cat("  UMAP2 ridgeline failed:", conditionMessage(e), "\n")
    })

    # ---- 3D centroid-distance ridgeline --------------------------
    cat("Creating 3D centroid-distance ridgeline...\n")
    tryCatch({
        coords_3d     <- as.matrix(
            analysis_data[, c("UMAP1", "UMAP2", "UMAP3")])
        post_cells_3d <- coords_3d[post_idx, , drop = FALSE]
        pre_cells_3d  <- coords_3d[pre_idx,  , drop = FALSE]

        if (nrow(post_cells_3d) > 0 && nrow(pre_cells_3d) > 0) {
            post_centroid <- colMeans(post_cells_3d)
            pre_centroid  <- colMeans(pre_cells_3d)

            trep_idx    <- post_idx | pre_idx
            trep_coords <- coords_3d[trep_idx, , drop = FALSE]
            trep_groups <- analysis_data$cell_group[trep_idx]

            dist_to_post <- sqrt(rowSums(
                sweep(trep_coords, 2, post_centroid)^2))
            dist_to_pre  <- sqrt(rowSums(
                sweep(trep_coords, 2, pre_centroid)^2))

            cd_data <- data.frame(
                Dist_Own_Centroid = ifelse(
                    trep_groups == "Post-IP_TREP",
                    dist_to_post, dist_to_pre),
                Cell_Group = factor(
                    trep_groups,
                    levels = c("Pre-IP_TREP", "Post-IP_TREP")),
                stringsAsFactors = FALSE)

            cd_colors <- c("Pre-IP_TREP"  = "#DDA0DD",
                           "Post-IP_TREP" = "#6666FF")

            cd_row <- stat_df[grepl("3D_Centroid",
                                    stat_df$Axis) &
                                   !is.na(stat_df$Axis), ]
            cd_cap <- if (nrow(cd_row) > 0) {
                paste0("Centroid sep=",
                       round(cd_row$Effect_Size[1], 3),
                       "  p=",
                       formatC(cd_row$P_Value[1],
                               format = "e", digits = 2),
                       " (BH=",
                       formatC(cd_row$BH_Adjusted_P_Value[1],
                               format = "e", digits = 2), ")")
            } else "Centroid distance test not available"

            cd_plot <- ggplot(cd_data,
                              aes(x = Dist_Own_Centroid,
                                  y = Cell_Group,
                                  fill = Cell_Group)) +
                geom_density_ridges(alpha = 0.5, scale = 1.5,
                                    rel_min_height = 0.01) +
                scale_fill_manual(values = cd_colors) +
                scale_y_discrete(
                    limits = rev(c("Pre-IP_TREP", "Post-IP_TREP")),
                    expand = expansion(mult = c(0, 0.2))) +
                coord_cartesian(clip = "off") +
                labs(
                    title    = paste("3D Centroid Distance -",
                                     patient_id, "| TREP cells"),
                    subtitle = paste(
                        "Euclidean distance of each TREP cell to",
                        "its own group centroid in 3D UMAP space\n",
                        "Separated peaks = distinct 3D zones;",
                        "merged peak = spatially interleaved"),
                    x = "3D distance to own centroid",
                    y = "Cell Type",
                    caption = cd_cap) +
                theme_ridges() +
                theme(
                    plot.title         = element_text(size = 14,
                                                      face = "bold"),
                    plot.subtitle      = element_text(size = 11),
                    axis.text.y        = element_text(size = 8),
                    axis.text.x        = element_text(
                        size = 8, margin = margin(t = 2)),
                    plot.caption       = element_text(size = 9,
                                                      hjust = 0),
                    panel.grid.major.x = element_line(
                        linetype = "dotted", color = "gray70"),
                    panel.grid.minor.x = element_blank(),
                    plot.margin        = margin(20, 5, 5, 5, "pt")) +
                guides(fill = "none")

            print(cd_plot)
        }
    }, error = function(e) {
        cat("  3D centroid ridgeline failed:", conditionMessage(e),
            "\n")
    })

    return(stat_df)
}


##############################################################################
# SECTION 5: Process all patients - 2D + 3D + export
##############################################################################
process_all_patients_clustering <- function(
        all_trajectories,
        all_trajectories_3d,
        output_file = "monocle3_clustering_FOXM1_DFG_FOXM1-pos.xlsx") {

    cat("\n=== PROCESSING ALL PATIENTS ===\n")

    colors <- c(
        "Pre-IP_non-TREP"  = "#CCCCCC",
        "Post-IP_non-TREP" = "#666666",
        "Pre-IP_TREP"      = "#DDA0DD",
        "Post-IP_TREP"     = "#6666FF"
    )

    all_2d_results <- data.frame()
    all_3d_results <- data.frame()

    for (patient_name in names(all_trajectories)) {
        cat("\nProcessing", patient_name, "...\n")

        # 2D
        tryCatch({
            res_2d <- analyze_cell_type_clustering(
                all_trajectories[[patient_name]],
                patient_name, colors)
            if (!is.null(res_2d))
                all_2d_results <- rbind(all_2d_results,
                                        res_2d$statistical_results)
        }, error = function(e) {
            cat("2D error for", patient_name, ":",
                conditionMessage(e), "\n")
        })

        # 3D
        tryCatch({
            res_3d <- analyze_cell_type_clustering_3d(
                all_trajectories_3d[[patient_name]],
                patient_name, colors)
            if (!is.null(res_3d))
                all_3d_results <- rbind(all_3d_results, res_3d)
        }, error = function(e) {
            cat("3D error for", patient_name, ":",
                conditionMessage(e), "\n")
        })
    }

    # ---- Global BH - 2D -----------------------------------------
    cat("\nApplying global BH correction (2D)...\n")
    if (nrow(all_2d_results) > 0) {
        p_rows <- !is.na(all_2d_results$P_Value)
        if (sum(p_rows) > 0) {
            all_2d_results$BH_Adjusted_P_Value <- NA
            all_2d_results$BH_Adjusted_P_Value[p_rows] <-
                p.adjust(all_2d_results$P_Value[p_rows], method = "BH")
            cat("2D BH applied to", sum(p_rows), "p-values.\n")
        }
    }

    # ---- Global BH - 3D -----------------------------------------
    cat("Applying global BH correction (3D)...\n")
    if (nrow(all_3d_results) > 0 &&
            "P_Value" %in% names(all_3d_results)) {
        p_rows <- !is.na(all_3d_results$P_Value)
        if (sum(p_rows) > 0) {
            all_3d_results$BH_Adjusted_P_Value_Global <- NA
            all_3d_results$BH_Adjusted_P_Value_Global[p_rows] <-
                p.adjust(all_3d_results$P_Value[p_rows], method = "BH")
            cat("3D global BH applied to", sum(p_rows), "p-values.\n")
        }
    }

    # ---- Export --------------------------------------------------
    cat("\nExporting to:", output_file, "\n")
    wb <- createWorkbook()

    addWorksheet(wb, "2D_Statistical_Results")
    writeData(wb, "2D_Statistical_Results", all_2d_results)

    if (nrow(all_3d_results) > 0) {
        addWorksheet(wb, "3D_Statistical_Results")
        writeData(wb, "3D_Statistical_Results", all_3d_results)
        cat("3D sheet:", nrow(all_3d_results), "rows written.\n")
    } else {
        cat("Warning: no 3D results to write.\n")
    }

    metadata <- data.frame(Description = c(
        "Dataset: Neftel GBM SS2 - FOXM1-positive adult malignant cells",
        "GOI: FOXM1 (x-axis) | Gene set: DFG composite",
        "DFG members: BIRC5, MKI67, CENPF, TOP2A, PBK, TROAP, NUSAP1",
        "Patients: MGH152 (dual), MGH110 (partial), MGH66 (dual), MGH100 (partial), MGH121 (partial)",
        "Filter: dual-filter (rho>=0.50 AND tau>=0.40) or partial, n_pos>=50",
        "IP source: inflection_points['patient|DFG'] from Part 7",
        "TREP: top (dev_expl x n) cells by EP = 1 - model_res^2 / null_res^2",
        "Null model: same weights as fitted GAM (all 1.0 for FOXM1-positive cells)",
        "=== 2D_Statistical_Results sheet ===",
        "UMAP1 axis: Shapiro-Wilk -> t-test (Cohen's d) or Mann-Whitney (Cliff's delta)",
        "Silhouette_Score_UMAP_Separation: 2D silhouette for Post-IP_TREP vs Pre-IP_TREP",
        "BH_Adjusted_P_Value: global BH correction across all 5 patients",
        "=== 3D_Statistical_Results sheet ===",
        "Axis=UMAP1/UMAP2/UMAP3: per-axis Shapiro->t/MW tests, Post-IP_TREP vs Pre-IP_TREP",
        "  UMAP3 is the key separating axis for loop/helix structures (e.g. MGH152)",
        "  MGH152: UMAP1 p=0.487 (non-significant in 2D) but UMAP3 captures 3D separation",
        "Axis=3D_Centroid_Distance: rotation-invariant test",
        "  Each group's cells measured by distance to own 3D centroid",
        "  Effect_Size = 3D Euclidean centroid separation (not Cliff's delta)",
        "  Captures helix/loop separation regardless of axis orientation",
        "Axis=3D_Silhouette: silhouette score on full 3D UMAP distance matrix",
        "Axis=Hotelling_T2_3D: multivariate test across UMAP1+2+3 simultaneously",
        "  Accounts for axis correlations; more powerful than 3 separate univariate tests",
        "  P_Value from F-approximation of T2 statistic (exact for two-group case)",
        "  Effect_Size = Mahalanobis D2 (group centroid separation in covariance units)",
        "  Significant Hotelling + non-significant univariate = oblique 3D separation",
        "Axis=NNP_3D: Nearest-Neighbour Purity (k=5) in 3D UMAP space",
        "  For each TREP cell: fraction of k=5 3D neighbours in same group (range 0-1)",
        "  Arc/helix-aware: does NOT assume spherical clusters unlike centroid distance",
        "  NNP > 0.5 = spatial segregation; NNP = 0.5 = random mixing (null hypothesis)",
        "  Effect_Size = Cohen's d or Cliff's delta comparing Post vs Pre purity scores",
        "  High NNP + low centroid sep = arc-shaped helix; High NNP + high sep = spherical",
        "Axis=PERMANOVA_3D: adonis2 on Euclidean distance matrix of 3D UMAP (999 perms)",
        "  Non-parametric multivariate test; no normality assumption; permutation p-value",
        "  Effect_Size = R2 (proportion of total variance explained by group membership)",
        "  R2 > 0.1 is generally meaningful for single-cell UMAP data",
        "  Preferred over Hotelling T2 when UMAP coordinates are non-normal (most cases)",
        "Axis=PERMDISP_3D: betadisper on same Euclidean distance matrix (999 perms)",
        "  Tests homogeneity of multivariate dispersions (within-group variance equality)",
        "  Post_IP_TREP_Mean = mean distance-to-centroid for Post-IP_TREP group",
        "  Pre_IP_TREP_Mean  = mean distance-to-centroid for Pre-IP_TREP group",
        "  Effect_Size = dispersion ratio Post/Pre (>1: Post more dispersed, <1: Pre more)",
        "  INTERPRETATION PAIRING:",
        "    PERMANOVA sig + PERMDISP NOT sig -> clean centroid separation (best case)",
        "    PERMANOVA sig + PERMDISP sig     -> location AND spread differ; interpret with caution",
        "    PERMANOVA NOT sig                -> no detectable 3D separation",
        "BH_Adjusted_P_Value: within-patient BH correction",
        "BH_Adjusted_P_Value_Global: global BH correction across all patients x 3D tests",
        "=== Ridgeline plots (per patient) ===",
        "1. 2D UMAP1 ridgeline: standard 2D UMAP1 axis distributions",
        "2. 3D UMAP3 ridgeline: key separator for loop/helix structures",
        "3. 3D UMAP2 ridgeline: often carries strongest per-axis signal",
        "   e.g. MGH152: UMAP2 Cliff delta=-0.773 p=5.76e-09 vs UMAP3 p=3.22e-02",
        "4. 3D centroid-distance ridgeline: own-centroid distance per TREP cell",
        "   Separated peaks = distinct 3D zones; merged peaks = spatially interleaved"
    ))
    addWorksheet(wb, "Metadata")
    writeData(wb, "Metadata", metadata)

    saveWorkbook(wb, output_file, overwrite = TRUE)
    cat("Saved to:", output_file, "\n")

    # ---- Summary printout ----------------------------------------
    cat("\n=== 2D CLUSTERING SUMMARY ===\n")
    if (nrow(all_2d_results) > 0) {
        comp_2d <- all_2d_results %>%
            dplyr::filter(Cell_Group == "Post-IP_TREP_vs_Pre-IP_TREP" &
                              !grepl("Silhouette", Test_Type))
        if (nrow(comp_2d) > 0) {
            print(comp_2d %>%
                      dplyr::select(Sample, P_Value,
                                    BH_Adjusted_P_Value) %>%
                      dplyr::mutate(
                          Raw_Sig = P_Value < 0.05,
                          BH_Sig  = BH_Adjusted_P_Value < 0.05))
        }
    }

    cat("\n=== 3D CLUSTERING SUMMARY ===\n")
    if (nrow(all_3d_results) > 0) {
        axis_sum <- all_3d_results %>%
            dplyr::filter(!is.na(P_Value)) %>%
            dplyr::select(Sample, Axis, Test_Type, P_Value,
                          BH_Adjusted_P_Value,
                          BH_Adjusted_P_Value_Global,
                          Effect_Size) %>%
            dplyr::mutate(
                Raw_Sig    = P_Value < 0.05,
                Global_Sig = BH_Adjusted_P_Value_Global < 0.05)
        print(axis_sum)

        sil_3d <- all_3d_results %>%
            dplyr::filter(grepl("Silhouette", Test_Type)) %>%
            dplyr::select(Sample, Axis, Effect_Size)
        if (nrow(sil_3d) > 0) {
            cat("\n3D Silhouette scores:\n")
            print(sil_3d)
        }
    }

    return(list(results_2d = all_2d_results,
                results_3d = all_3d_results))
}


# ------------------------------------------------------------------
# SECTION 6: Run clustering analysis
# ------------------------------------------------------------------
cat("\n\nStarting clustering analysis...\n")
clustering_analysis_results <- process_all_patients_clustering(
    all_trajectories,
    all_trajectories_3d
)

cat("\n========== Part 9 complete ==========\n")
cat("Output: monocle3_clustering_FOXM1_DFG_FOXM1-pos.xlsx\n")
cat("  Sheet 1: 2D_Statistical_Results\n")
cat("  Sheet 2: 3D_Statistical_Results\n")
cat("    - Axis: UMAP1 / UMAP2 / UMAP3 (per-axis Shapiro->t/MW)\n")
cat("    - Axis: 3D_Centroid_Distance (rotation-invariant,\n")
cat("            Effect_Size = Euclidean centroid separation)\n")
cat("    - Axis: 3D_Silhouette (full 3D distance matrix)\n")
cat("    - Axis: Hotelling_T2_3D (joint UMAP1+2+3 multivariate test,\n")
cat("            Effect_Size = Mahalanobis D2)\n")
cat("    - Axis: NNP_3D (Nearest-Neighbour Purity k=5, arc-aware)\n")
cat("    - Axis: PERMANOVA_3D (adonis2, Euclidean, 999 perms,\n")
cat("            Effect_Size = R2; non-parametric gold standard)\n")
cat("    - Axis: PERMDISP_3D (betadisper, 999 perms,\n")
cat("            Effect_Size = Post/Pre dispersion ratio)\n")
cat("    Interpretation: PERMANOVA sig + PERMDISP not sig = clean separation\n")
cat("                    PERMANOVA sig + PERMDISP sig     = location + spread differ\n")
cat("    - Columns: BH_Adjusted_P_Value (within-patient) +\n")
cat("               BH_Adjusted_P_Value_Global (all patients)\n")
cat("  Sheet 3: Metadata\n")
cat("Ridgeline plots (per patient):\n")
cat("  - 2D UMAP1 ridgeline\n")
cat("  - 3D UMAP3 ridgeline  [key for loop structures]\n")
cat("  - 3D UMAP2 ridgeline  [often strongest per-axis signal]\n")
cat("  - 3D centroid-distance ridgeline\n")
cat("3D interactive UMAP: RStudio Viewer (not saved to disk)\n")


# ================================================
# R SCRIPT: Create a FULLY STANDALONE combined HTML
# Embeds all 6 Monocle3 3D plots using srcdoc (no external files needed)
# ================================================

library(htmltools)

# Set the folder containing 6 original HTML files
folder <- "C: #directory"

# List of filenames (update if names change)
filenames <- c(
    "3D_MGH152_v2.html",
    "3D_MGH66_v2.html",
    "3D_MGH110_v2.html",
    "3D_MGH100_v2.html",
    "3D_MGH121_v2.html",
    "3D_MGH104_v2.html"
)

# Nice sample labels
samples <- c("MGH152", "MGH66", "MGH110", "MGH100", "MGH121", "MGH104")

# Read the full HTML content of each file
html_contents <- lapply(file.path(folder, filenames), function(f) {
    paste(readLines(f, warn = FALSE), collapse = "\n")
})

# Build the combined standalone page
combined_html <- tagList(
    
    tags$head(
        tags$style(HTML("
      body {
        font-family: Arial, Helvetica, sans-serif;
        margin: 40px;
        background-color: #f8f9fa;
        line-height: 1.6;
      }
      h1 {
        text-align: center;
        color: #2c3e50;
      }
      h2 {
        color: #2980b9;
        border-bottom: 3px solid #3498db;
        padding-bottom: 10px;
      }
      iframe {
        border: 3px solid #bdc3c7;
        border-radius: 10px;
        box-shadow: 0 6px 20px rgba(0,0,0,0.15);
        display: block;
        margin: 25px auto;
        background: white;
      }
      .note {
        text-align: center;
        font-size: 16px;
        color: #7f8c8d;
        margin-bottom: 40px;
      }
      .separator {
        margin: 70px 0;
        border-top: 2px dashed #95a5a6;
      }
    "))
    ),
    
    h1("3D Monocle3 Plots of CEP-IP Framework (Neftel GBM Dataset; FOXM1+ Cells)"),
    p("Each 3D plot is fully interactive (rotate, zoom, hover). Scroll down to switch between samples.", 
      class = "note"),
    
    # Embed each plot using srcdoc
    lapply(seq_along(html_contents), function(i) {
        tagList(
            h2(paste0("Sample: ", samples[i])),
            tags$iframe(
                srcdoc = html_contents[[i]],   # <-- full HTML content embedded here
                width = "100%",
                height = "950px",              # Increase if plots feel cramped
                frameborder = "0",
                allowfullscreen = "true",
                sandbox = "allow-scripts allow-same-origin"  # Needed for 3D interactivity
            ),
            tags$div(class = "separator")
        )
    })
)

# Save as a single standalone HTML file
output_file <- file.path(folder, "3D_Monocle3_Combined_Standalone.html")

save_html(
    combined_html,
    file = output_file,
    libdir = NULL   # No extra folder created
)

# Open it automatically
browseURL(output_file)

cat("✅ SUCCESS!\n")
cat("Standalone file created:\n")
cat(basename(output_file), "\n")
