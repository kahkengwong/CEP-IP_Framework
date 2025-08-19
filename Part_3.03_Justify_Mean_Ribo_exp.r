#####################################################################
# Part 3.01: Justifications to Average Genes Expression for Modeling
#####################################################################
# Load required libraries
library(ggplot2)
library(psych)
library(openxlsx)

# Function to validate and analyze gene sets (simplified)
validate_gene_set_simple <- function(sample_name, gene_set_name, genes) {
    cat("\n", rep("=", 80), "\n")
    cat(gene_set_name, "GENE RELIABILITY ANALYSIS FOR:", sample_name)
    cat("\n", rep("=", 80), "\n")
    
    # Extract expression data for the gene set
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE))
    
    # Get expression matrix for gene set
    gene_expression <- GetAssayData(prostate_results$seurat_obj, assay = "RNA", slot = "data")[genes, sample_cells]
    gene_expression <- as.matrix(gene_expression)
    gene_df <- as.data.frame(t(gene_expression))
    
    cat("Sample:", sample_name, "\n")
    cat("Gene Set:", gene_set_name, "\n")
    cat("Number of cells:", nrow(gene_df), "\n")
    cat("Genes analyzed:", paste(genes, collapse = ", "), "\n\n")
    
    # === 1. DESCRIPTIVE STATISTICS ===
    cat("=== 1. DESCRIPTIVE STATISTICS ===\n")
    desc_stats <- data.frame(
        Gene = genes,
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
        gene_set = gene_set_name,
        n_cells = nrow(gene_df),
        n_genes = length(genes),
        genes = genes,
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
    cat("RUNNING SIMPLIFIED GENE SET RELIABILITY ANALYSIS")
    cat("\n", rep("=", 100), "\n")
    
    # Define gene sets
    gene_sets <- list(
        Ribo = c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26"),
        AR = c("KLK4", "KLK2", "KLK3", "PDLIM5", "ABHD2", "ALDH1A3", "SORD"),
        PI3K_AKT = c("PIK3CA", "AKT1", "PTEN", "MTOR", "TSC2", "FOXO3", "GSK3B"),
        mTOR = c("MTOR", "RPTOR", "RICTOR", "MLST8", "AKT1S1", "DEPTOR", "PRR5"),
        GSK3B = c("GSK3B", "CTNNB1", "AXIN1", "CSNK1A1", "APC", "FZD1", "LRP6"),
        NFKB = c("NFKB1", "RELA", "TNFAIP3", "NFKBIA", "IKBKB", "TRAF6", "NFKB2"),
        WNT = c("CTNNB1", "APC", "AXIN1", "GSK3B", "LEF1", "TCF7", "FZD1")
    )
    
    all_validation_results <- list()
    
    # Run analysis for each gene set and sample
    for(gene_set_name in names(gene_sets)) {
        for(sample_name in tumor_samples) {
            tryCatch({
                validation_result <- validate_gene_set_simple(sample_name, gene_set_name, gene_sets[[gene_set_name]])
                if(!is.null(validation_result)) {
                    result_key <- paste(gene_set_name, sample_name, sep = "_")
                    all_validation_results[[result_key]] <- validation_result
                    cat("Completed validation for:", gene_set_name, "-", sample_name, "\n")
                }
            }, error = function(e) {
                cat("Error in validation for", gene_set_name, "-", sample_name, ":", conditionMessage(e), "\n")
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
        desc_with_ids$Gene_Set <- result$gene_set
        desc_with_ids$Sample <- result$sample
        desc_with_ids <- desc_with_ids[, c("Gene_Set", "Sample", "Gene", "Mean", "SD", "Min", "Max", "CV")]
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
    
    # === SHEET 4: RELIABILITY BREAKDOWN BY GENE SET ===
    cat("Creating reliability breakdown sheet...\n")
    
    # Calculate summary statistics by gene set
    reliability_summary <- summary_df %>%
        group_by(Gene_Set) %>%
        summarise(
            N_Samples = n(),
            Mean_Cronbach = round(mean(Cronbach_Alpha, na.rm = TRUE), 3),
            Mean_McDonald = round(mean(McDonald_Omega, na.rm = TRUE), 3),
            Mean_KMO = round(mean(KMO_Score, na.rm = TRUE), 3),
            Mean_Consolidated = round(mean(Consolidated_Score, na.rm = TRUE), 3),
            .groups = 'drop'
        ) %>%
        as.data.frame()
    
    addWorksheet(wb, "Reliability_Summary")
    writeData(wb, "Reliability_Summary", reliability_summary)
    
    # Format reliability summary sheet
    addStyle(wb, "Reliability_Summary", 
             style = createStyle(textDecoration = "bold", fgFill = "#D3D3D3"), 
             rows = 1, cols = 1:ncol(reliability_summary))
    setColWidths(wb, "Reliability_Summary", cols = 1:ncol(reliability_summary), widths = "auto")
    
    # Save Excel file
    excel_filename <- "Gene_Set_Reliability_Analysis.xlsx"
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
    cat("\nAll results exported to:", excel_filename, "\n")
    
    cat("\nExcel file contains 4 sheets:\n")
    cat("1. Summary: Main results table with reliability scores\n")
    cat("2. Descriptive_Stats: Gene-level descriptive statistics\n") 
    cat("3. Gene_Set_Overview: List of genes in each set\n")
    cat("4. Reliability_Summary: Summary statistics by gene set\n")
    
    # Print summary statistics
    cat("\n", rep("=", 80), "\n")
    cat("OVERALL RELIABILITY SUMMARY")
    cat("\n", rep("=", 80), "\n")
    
    print(reliability_summary)
    
    cat("\nTotal analyses performed:", nrow(summary_df), "\n")
    
    return(list(
        individual_results = all_validation_results,
        summary = summary_df,
        reliability_summary = reliability_summary,
        gene_sets = gene_sets
    ))
}

# Load dplyr for summary statistics (if not already loaded)
if(!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
}

# Usage instructions
cat("\nSimplified gene set reliability analysis code is ready!")
cat("\nUse: validation_results <- validate_all_samples_simple()")
cat("\nOr for individual gene sets: validate_gene_set_simple('HYW_4701_Tumor', 'Ribo', c('RPL10', 'RPL27', 'RPL28', 'RPS2', 'RPS8', 'RPS12', 'RPS26'))")
cat("\n\nConsolidated Score = Average of (Cronbach's α + McDonald's ω + KMO)")

# Run complete analysis
validation_results <- validate_all_samples_simple()
