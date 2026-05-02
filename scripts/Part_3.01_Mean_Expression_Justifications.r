#####################################################################
# Part 3.01: Justifications to Average Genes Expression for Modeling
#####################################################################
# Load required libraries
library(Seurat)
library(psych) 
library(openxlsx)
library(dplyr)

# Function to validate and analyze gene sets (simplified)
validate_gene_set_simple <- function(sample_name, gene_set_name, genes, analysis_type = "Tumor") {
    cat("\n", rep("=", 80), "\n")
    cat(gene_set_name, "GENE RELIABILITY ANALYSIS FOR:", sample_name, "(", analysis_type, ")")
    cat("\n", rep("=", 80), "\n")
    
    # Extract expression data based on analysis type
    if(analysis_type == "BP") {
        # Extract cluster 3 cells from NonCa samples using non_cancerous objects
        sample_cells <- WhichCells(non_cancerous_results$seurat_obj, 
                                   cells = grep(sample_name, colnames(non_cancerous_results$seurat_obj), value = TRUE))
        # Filter for cluster 3 cells
        cluster_info <- non_cancerous_results$seurat_obj@meta.data[sample_cells, "seurat_clusters"]
        cluster_3_cells <- sample_cells[cluster_info == 3]
        analysis_cells <- cluster_3_cells
        
        # Use non_cancerous_seurat for gene expression data
        seurat_obj <- non_cancerous_seurat
        
    } else {
        # Default tumor analysis (all cells from sample) using prostate objects
        analysis_cells <- WhichCells(prostate_results$seurat_obj, 
                                     cells = grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE))
        # Use prostate_ca_seurat for gene expression data
        seurat_obj <- prostate_ca_seurat
    }
    
    # Check if we have enough cells
    if(length(analysis_cells) < 10) {
        cat("Warning: Only", length(analysis_cells), "cells found. Skipping analysis.\n")
        return(NULL)
    }
    
    # Get expression matrix for gene set from the appropriate Seurat object
    # Ensure we're using the RNA assay
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Subset the Seurat object to the selected cells
    cluster_obj <- subset(seurat_obj, cells = analysis_cells)
    
    # Get normalized data
    gene_expression <- GetAssayData(cluster_obj, slot = "data", assay = "RNA")
    
    # Check if all genes are present
    genes_present <- genes[genes %in% rownames(gene_expression)]
    if(length(genes_present) < 3) {
        cat("Warning: Only", length(genes_present), "genes found in expression data. Skipping analysis.\n")
        return(NULL)
    }
    
    # Calculate average expression for the gene set
    avg_expression <- colMeans(gene_expression[genes_present, ])
    
    # Create data frame with individual gene expressions for reliability analysis
    gene_df <- as.data.frame(t(gene_expression[genes_present, ]))
    
    cat("Sample:", sample_name, "\n")
    cat("Analysis Type:", analysis_type, "\n")
    cat("Gene Set:", gene_set_name, "\n")
    cat("Number of cells:", nrow(gene_df), "\n")
    cat("Genes analyzed:", paste(genes_present, collapse = ", "), "\n\n")
    
    # === 1. DESCRIPTIVE STATISTICS ===
    cat("=== 1. DESCRIPTIVE STATISTICS ===\n")
    desc_stats <- data.frame(
        Gene = genes_present,
        Mean = round(apply(gene_df, 2, mean), 3),
        SD = round(apply(gene_df, 2, sd), 3),
        Min = round(apply(gene_df, 2, min), 3),
        Max = round(apply(gene_df, 2, max), 3),
        CV = round(apply(gene_df, 2, function(x) sd(x)/mean(x) * 100), 2)
    )
    print(desc_stats)
    
    # === 2. RELIABILITY ANALYSIS ===
    cat("\n=== 2. RELIABILITY ANALYSIS ===\n")
    
    # Cronbach's Alpha
    cronbach_result <- psych::alpha(gene_df)
    cronbach_alpha <- cronbach_result$total$std.alpha
    cat("Cronbach's Alpha:", round(cronbach_alpha, 3), "\n")
    
    # McDonald's Omega
    omega_result <- psych::omega(gene_df, nfactors = 1, plot = FALSE)
    mcdonald_omega <- omega_result$omega_h
    cat("McDonald's Omega (ω):", round(mcdonald_omega, 3), "\n")
    
    # === 3. KAISER-MEYER-OLKIN (KMO) TEST ===
    cat("\n=== 3. KAISER-MEYER-OLKIN (KMO) MEASURE ===\n")
    
    kmo_result <- psych::KMO(gene_df)
    kmo_score <- kmo_result$MSA
    cat("Overall KMO:", round(kmo_score, 3), "\n")
    
    # === 4. CONSOLIDATED RELIABILITY SCORE ===
    cat("\n=== 4. CONSOLIDATED RELIABILITY SCORE ===\n")
    
    # Calculate consolidated score (average of the three metrics)
    consolidated_score <- mean(c(cronbach_alpha, mcdonald_omega, kmo_score), na.rm = TRUE)
    cat("Consolidated Score (Average):", round(consolidated_score, 3), "\n")
    
    # Component interpretations
    cat("\nComponent Scores:\n")
    cat("  Cronbach's α:", round(cronbach_alpha, 3), 
        ifelse(cronbach_alpha >= 0.8, " (Good)", ifelse(cronbach_alpha >= 0.5, " (Acceptable)", " (Poor)")), "\n")
    cat("  McDonald's ω:", round(mcdonald_omega, 3), 
        ifelse(mcdonald_omega >= 0.8, " (Good)", ifelse(mcdonald_omega >= 0.5, " (Acceptable)", " (Poor)")), "\n")
    cat("  KMO:", round(kmo_score, 3), 
        ifelse(kmo_score >= 0.8, " (Excellent)", ifelse(kmo_score >= 0.6, " (Good)", ifelse(kmo_score >= 0.5, " (Adequate)", " (Poor)"))), "\n")
    
    # Return results
    results <- list(
        sample = sample_name,
        analysis_type = analysis_type,
        gene_set = gene_set_name,
        n_cells = nrow(gene_df),
        n_genes = length(genes_present),
        genes = genes_present,
        descriptive_stats = desc_stats,
        cronbach_alpha = cronbach_alpha,
        mcdonald_omega = mcdonald_omega,
        kmo_score = kmo_score,
        consolidated_score = consolidated_score
    )
    
    return(results)
}

# Function to run validation across all samples and gene sets
validate_all_samples_simple <- function() {
    cat("\n", rep("=", 100), "\n")
    cat("RUNNING GENE SET RELIABILITY ANALYSIS FOR RIBO AND AR")
    cat("\n", rep("=", 100), "\n")
    
    # Define sample groups
    tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                       "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                       "HYW_5755_Tumor")
    
    non_ca_samples <- c("HYW_4701_Benign", "HYW_4847_Benign", "HYW_4880_Benign", 
                        "HYW_4881_Benign", "HYW_5386_Benign", "HYW_5742_Benign", 
                        "HYW_5755_Benign")
    
    # Define gene sets - Only Ribo and AR
    gene_sets <- list(
        Ribo = c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26"),
        AR = c("KLK4", "KLK2", "KLK3", "PDLIM5", "ABHD2", "ALDH1A3", "SORD")
    )
    
    all_validation_results <- list()
    
    # 1. Run analysis for tumor samples
    cat("\n=== ANALYZING TUMOR SAMPLES ===\n")
    for(gene_set_name in names(gene_sets)) {
        for(sample_name in tumor_samples) {
            tryCatch({
                validation_result <- validate_gene_set_simple(sample_name, gene_set_name, gene_sets[[gene_set_name]], "Tumor")
                if(!is.null(validation_result)) {
                    result_key <- paste("Tumor", gene_set_name, sample_name, sep = "_")
                    all_validation_results[[result_key]] <- validation_result
                    cat("Completed tumor validation for:", gene_set_name, "-", sample_name, "\n")
                }
            }, error = function(e) {
                cat("Error in tumor validation for", gene_set_name, "-", sample_name, ":", conditionMessage(e), "\n")
            })
        }
    }
    
    # 2. Run analysis for BP (cluster 3 in NonCa samples)
    cat("\n=== ANALYZING BP (CLUSTER 3 IN NONCA SAMPLES) ===\n")
    for(gene_set_name in names(gene_sets)) {
        for(sample_name in non_ca_samples) {
            tryCatch({
                validation_result <- validate_gene_set_simple(sample_name, gene_set_name, gene_sets[[gene_set_name]], "BP")
                if(!is.null(validation_result)) {
                    result_key <- paste("BP", gene_set_name, sample_name, sep = "_")
                    all_validation_results[[result_key]] <- validation_result
                    cat("Completed BP validation for:", gene_set_name, "-", sample_name, "\n")
                }
            }, error = function(e) {
                cat("Error in BP validation for", gene_set_name, "-", sample_name, ":", conditionMessage(e), "\n")
            })
        }
    }
    
    # Create comprehensive Excel workbook
    wb <- createWorkbook()
    
    # === SHEET 1: SUMMARY TABLE ===
    cat("\n", rep("=", 80), "\n")
    cat("CREATING SUMMARY TABLE")
    cat("\n", rep("=", 80), "\n")
    
    summary_df <- data.frame(
        Analysis_Type = character(),
        Gene_Set = character(),
        Sample = character(),
        N_Cells = numeric(),
        N_Genes = numeric(),
        Cronbach_Alpha = numeric(),
        McDonald_Omega = numeric(),
        KMO_Score = numeric(),
        Consolidated_Score = numeric(),
        stringsAsFactors = FALSE
    )
    
    for(result_key in names(all_validation_results)) {
        result <- all_validation_results[[result_key]]
        
        summary_df <- rbind(summary_df, data.frame(
            Analysis_Type = result$analysis_type,
            Gene_Set = result$gene_set,
            Sample = result$sample,
            N_Cells = result$n_cells,
            N_Genes = result$n_genes,
            Cronbach_Alpha = round(result$cronbach_alpha, 3),
            McDonald_Omega = round(result$mcdonald_omega, 3),
            KMO_Score = round(result$kmo_score, 3),
            Consolidated_Score = round(result$consolidated_score, 3),
            stringsAsFactors = FALSE
        ))
    }
    
    # Add summary sheet
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", summary_df)
    
    # Format summary sheet
    addStyle(wb, "Summary", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(summary_df))
    setColWidths(wb, "Summary", cols = 1:ncol(summary_df), widths = "auto")
    
    print(summary_df)
    
    # === SHEET 2: DESCRIPTIVE STATISTICS ===
    cat("\nCreating descriptive statistics sheet...\n")
    
    # Compile all descriptive statistics
    all_desc_stats <- data.frame()
    
    for(result_key in names(all_validation_results)) {
        result <- all_validation_results[[result_key]]
        desc_with_ids <- result$descriptive_stats
        desc_with_ids$Analysis_Type <- result$analysis_type
        desc_with_ids$Gene_Set <- result$gene_set
        desc_with_ids$Sample <- result$sample
        desc_with_ids <- desc_with_ids[, c("Analysis_Type", "Gene_Set", "Sample", "Gene", "Mean", "SD", "Min", "Max", "CV")]
        all_desc_stats <- rbind(all_desc_stats, desc_with_ids)
    }
    
    addWorksheet(wb, "Descriptive_Stats")
    writeData(wb, "Descriptive_Stats", all_desc_stats)
    
    # Format descriptive stats sheet
    addStyle(wb, "Descriptive_Stats", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(all_desc_stats))
    setColWidths(wb, "Descriptive_Stats", cols = 1:ncol(all_desc_stats), widths = "auto")
    
    # === SHEET 3: GENE SET OVERVIEW ===
    cat("Creating gene set overview sheet...\n")
    
    gene_set_overview <- data.frame(
        Gene_Set = names(gene_sets),
        Genes = sapply(gene_sets, function(x) paste(x, collapse = ", ")),
        N_Genes = sapply(gene_sets, length),
        stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Gene_Set_Overview")
    writeData(wb, "Gene_Set_Overview", gene_set_overview)
    
    # Format gene set overview sheet
    addStyle(wb, "Gene_Set_Overview", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(gene_set_overview))
    setColWidths(wb, "Gene_Set_Overview", cols = 1:ncol(gene_set_overview), widths = c(12, 80, 10))
    
    # === SHEET 4: ANALYSIS BREAKDOWN BY TYPE AND GENE SET ===
    cat("Creating analysis breakdown sheet...\n")
    
    # Calculate summary statistics by analysis type and gene set
    analysis_summary <- summary_df %>%
        group_by(Analysis_Type, Gene_Set) %>%
        summarise(
            N_Samples = n(),
            Mean_Cronbach = round(mean(Cronbach_Alpha, na.rm = TRUE), 3),
            Mean_McDonald = round(mean(McDonald_Omega, na.rm = TRUE), 3),
            Mean_KMO = round(mean(KMO_Score, na.rm = TRUE), 3),
            Mean_Consolidated = round(mean(Consolidated_Score, na.rm = TRUE), 3),
            .groups = 'drop'
        ) %>%
        as.data.frame()
    
    addWorksheet(wb, "Analysis_Summary")
    writeData(wb, "Analysis_Summary", analysis_summary)
    
    # Format analysis summary sheet
    addStyle(wb, "Analysis_Summary", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(analysis_summary))
    setColWidths(wb, "Analysis_Summary", cols = 1:ncol(analysis_summary), widths = "auto")
    
    # === SHEET 5: SAMPLE OVERVIEW ===
    cat("Creating sample overview sheet...\n")
    
    sample_overview <- data.frame(
        Analysis_Type = c(rep("Tumor", length(tumor_samples)), 
                         rep("BP", length(non_ca_samples))),
        Sample = c(tumor_samples, non_ca_samples),
        Description = c(rep("Tumor samples (all cells)", length(tumor_samples)),
                       rep("NonCa samples - Cluster 3 (BP)", length(non_ca_samples))),
        Seurat_Object_Used = c(rep("prostate_ca_seurat", length(tumor_samples)),
                              rep("non_cancerous_seurat", length(non_ca_samples))),
        stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Sample_Overview")
    writeData(wb, "Sample_Overview", sample_overview)
    
    # Format sample overview sheet
    addStyle(wb, "Sample_Overview", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(sample_overview))
    setColWidths(wb, "Sample_Overview", cols = 1:ncol(sample_overview), widths = c(15, 20, 40, 25))
    
    # Save Excel file
    excel_filename <- "Gene_Set_Reliability_Analysis_Ribo_AR.xlsx"
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
    cat("\nAll results exported to:", excel_filename, "\n")
    
    cat("\nExcel file contains 5 sheets:\n")
    cat("1. Summary: Main results table for Tumor and BP analysis types\n")
    cat("2. Descriptive_Stats: Gene-level descriptive statistics\n") 
    cat("3. Gene_Set_Overview: List of genes in Ribo and AR sets\n")
    cat("4. Analysis_Summary: Summary statistics by analysis type and gene set\n")
    cat("5. Sample_Overview: Description of analysis types, samples, and Seurat objects used\n")
    
    # Print summary statistics
    cat("\n", rep("=", 80), "\n")
    cat("OVERALL ANALYSIS SUMMARY")
    cat("\n", rep("=", 80), "\n")
    
    print(analysis_summary)
    
    cat("\nTotal analyses performed:", nrow(summary_df), "\n")
    cat("  Tumor analyses:", sum(summary_df$Analysis_Type == "Tumor"), "\n")
    cat("  BP analyses:", sum(summary_df$Analysis_Type == "BP"), "\n")
    
    return(list(
        individual_results = all_validation_results,
        summary = summary_df,
        analysis_summary = analysis_summary,
        gene_sets = gene_sets,
        sample_overview = sample_overview
    ))
}

# Load dplyr for summary statistics (if not already loaded)
if(!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
}

# Usage instructions
cat("\nGene set reliability analysis for Ribo and AR is ready!")
cat("\nAnalyzes: Tumor samples and BP (cluster 3 in NonCa samples)")
cat("\nUse: validation_results <- validate_all_samples_simple()")
cat("\n\nConsolidated Score = Average of (Cronbach's α + McDonald's ω + KMO)")

# Run complete analysis
validation_results <- validate_all_samples_simple()
