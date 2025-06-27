##################################################################
# Part 3.01: Justifications to Average 7 Ribosomal Genes for Modeling with GAM 
##################################################################
# Load required libraries
library(corrplot)
library(psych)
library(factoextra)
library(FactoMineR)
library(GGally)
library(pheatmap)
library(RColorBrewer)
library(energy)

# Function to capture console output
capture_console_output <- function(expr) {
    output <- capture.output(expr)
    return(paste(output, collapse = "\n"))
}

# Function to validate and analyze 7 ribosomal genes
validate_ribosomal_genes <- function(sample_name, ribosomal_genes = c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26")) {
    cat("\n", rep("=", 80), "\n")
    cat("RIBOSOMAL GENE VALIDATION ANALYSIS FOR:", sample_name)
    cat("\n", rep("=", 80), "\n")
    
    # Extract expression data for the 7 ribosomal genes
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE))
    
    # Get expression matrix for ribosomal genes
    ribo_expression <- GetAssayData(prostate_results$seurat_obj, assay = "RNA", slot = "data")[ribosomal_genes, sample_cells]
    ribo_expression <- as.matrix(ribo_expression)
    ribo_df <- as.data.frame(t(ribo_expression))
    
    cat("Sample:", sample_name, "\n")
    cat("Number of cells:", nrow(ribo_df), "\n")
    cat("Ribosomal genes analyzed:", paste(ribosomal_genes, collapse = ", "), "\n\n")
    
    # === 1. DESCRIPTIVE STATISTICS ===
    cat("=== 1. DESCRIPTIVE STATISTICS ===\n")
    desc_stats <- data.frame(
        Gene = ribosomal_genes,
        Mean = round(apply(ribo_df, 2, mean), 3),
        SD = round(apply(ribo_df, 2, sd), 3),
        Min = round(apply(ribo_df, 2, min), 3),
        Max = round(apply(ribo_df, 2, max), 3),
        CV = round(apply(ribo_df, 2, function(x) sd(x)/mean(x) * 100), 2)
    )
    print(desc_stats)
    
    # === 2. DISTANCE CORRELATION ANALYSIS ===
    cat("\n=== 2. DISTANCE CORRELATION ANALYSIS ===\n")
    
    # Distance correlations (captures linear and non-linear relationships)
    dcor_matrix <- matrix(0, nrow = length(ribosomal_genes), ncol = length(ribosomal_genes))
    rownames(dcor_matrix) <- ribosomal_genes
    colnames(dcor_matrix) <- ribosomal_genes
    
    for(i in 1:length(ribosomal_genes)) {
        for(j in 1:length(ribosomal_genes)) {
            dcor_matrix[i, j] <- dcor(ribo_df[, i], ribo_df[, j])
        }
    }
    
    cat("Distance Correlation Matrix:\n")
    print(round(dcor_matrix, 3))
    
    # Correlation statistics
    upper_tri <- upper.tri(dcor_matrix)
    dcor_values <- dcor_matrix[upper_tri]
    
    cat("\nDistance Correlation Summary Statistics:\n")
    cat("  Mean:", round(mean(dcor_values), 3), "\n")
    cat("  Median:", round(median(dcor_values), 3), "\n")
    cat("  Range:", round(min(dcor_values), 3), "to", round(max(dcor_values), 3), "\n")
    cat("  Proportion > 0.5:", round(sum(dcor_values > 0.5) / length(dcor_values), 3), "\n")
    cat("  Proportion > 0.6:", round(sum(dcor_values > 0.6) / length(dcor_values), 3), "\n")
    cat("  Proportion > 0.7:", round(sum(dcor_values > 0.7) / length(dcor_values), 3), "\n")
    
    # === 3. PRINCIPAL COMPONENT ANALYSIS ===
    cat("\n=== 3. PRINCIPAL COMPONENT ANALYSIS ===\n")
    
    # Standardize data for PCA
    ribo_scaled <- scale(ribo_df)
    
    # Perform PCA
    pca_result <- prcomp(ribo_scaled, center = FALSE, scale. = FALSE)
    
    # PCA summary
    pca_summary <- summary(pca_result)
    cat("PCA Summary:\n")
    print(pca_summary$importance)
    
    # Variance explained
    var_explained <- pca_summary$importance[2, ] * 100
    cat("\nVariance Explained by Components:\n")
    for(i in 1:min(4, length(var_explained))) {
        cat("  PC", i, ":", round(var_explained[i], 2), "%\n")
    }
    
    cat("\nCumulative Variance Explained:\n")
    cum_var <- cumsum(var_explained)
    for(i in 1:min(4, length(cum_var))) {
        cat("  PC1-PC", i, ":", round(cum_var[i], 2), "%\n")
    }
    
    # PC1 loadings (gene contributions)
    pc1_loadings <- pca_result$rotation[, 1]
    cat("\nPC1 Loadings (Gene Contributions):\n")
    pc1_df <- data.frame(
        Gene = names(pc1_loadings),
        Loading = round(pc1_loadings, 3),
        Abs_Loading = round(abs(pc1_loadings), 3)
    )
    pc1_df <- pc1_df[order(pc1_df$Abs_Loading, decreasing = TRUE), ]
    print(pc1_df)
    
    # Check if genes load similarly on PC1 (similar signs and magnitudes)
    pc1_range <- max(pc1_loadings) - min(pc1_loadings)
    pc1_signs_same <- all(pc1_loadings > 0) || all(pc1_loadings < 0)
    
    cat("\nPC1 Analysis:\n")
    cat("  All genes have same sign on PC1:", pc1_signs_same, "\n")
    cat("  Range of PC1 loadings:", round(pc1_range, 3), "\n")
    cat("  Standard deviation of |PC1 loadings|:", round(sd(abs(pc1_loadings)), 3), "\n")
    
    # === 4. RELIABILITY ANALYSIS ===
    cat("\n=== 4. RELIABILITY ANALYSIS ===\n")
    
    # Cronbach's Alpha
    cronbach_result <- psych::alpha(ribo_df)
    cat("Cronbach's Alpha:", round(cronbach_result$total$std.alpha, 3), "\n")
    
    # McDonald's Omega
    omega_result <- psych::omega(ribo_df, nfactors = 1, plot = FALSE)
    mcdonald_omega <- omega_result$omega_h
    cat("McDonald's Omega (ω):", round(mcdonald_omega, 3), "\n")
    
    # Average Inter-Item Correlation (AIIC)
    cor_matrix <- cor(ribo_df)
    upper_tri_cor <- upper.tri(cor_matrix)
    aiic <- mean(cor_matrix[upper_tri_cor])
    cat("Average Inter-Item Correlation (AIIC):", round(aiic, 3), "\n")
    
    cat("\nReliability Interpretation:\n")
    if(cronbach_result$total$std.alpha >= 0.9) {
        cat("  Cronbach's α: Excellent internal consistency (α ≥ 0.9)\n")
    } else if(cronbach_result$total$std.alpha >= 0.8) {
        cat("  Cronbach's α: Good internal consistency (α ≥ 0.8)\n")
    } else if(cronbach_result$total$std.alpha >= 0.7) {
        cat("  Cronbach's α: Acceptable internal consistency (α ≥ 0.7)\n")
    } else {
        cat("  Cronbach's α: Poor internal consistency (α < 0.7)\n")
    }
    
    if(mcdonald_omega >= 0.9) {
        cat("  McDonald's ω: Excellent reliability (ω ≥ 0.9)\n")
    } else if(mcdonald_omega >= 0.8) {
        cat("  McDonald's ω: Good reliability (ω ≥ 0.8)\n")
    } else if(mcdonald_omega >= 0.7) {
        cat("  McDonald's ω: Acceptable reliability (ω ≥ 0.7)\n")
    } else {
        cat("  McDonald's ω: Poor reliability (ω < 0.7)\n")
    }
    
    if(aiic >= 0.5) {
        cat("  AIIC: Strong inter-item correlations (r ≥ 0.5)\n")
    } else if(aiic >= 0.3) {
        cat("  AIIC: Moderate inter-item correlations (r ≥ 0.3)\n")
    } else if(aiic >= 0.1) {
        cat("  AIIC: Weak inter-item correlations (r ≥ 0.1)\n")
    } else {
        cat("  AIIC: Very weak inter-item correlations (r < 0.1)\n")
    }
    
    # === 5. KAISER-MEYER-OLKIN (KMO) TEST ===
    cat("\n=== 5. KAISER-MEYER-OLKIN (KMO) MEASURE ===\n")
    
    kmo_result <- psych::KMO(ribo_df)
    cat("Overall KMO:", round(kmo_result$MSA, 3), "\n")
    cat("Individual KMO values:\n")
    for(i in 1:length(kmo_result$MSAi)) {
        cat("  ", names(kmo_result$MSAi)[i], ":", round(kmo_result$MSAi[i], 3), "\n")
    }
    
    cat("KMO Interpretation:\n")
    if(kmo_result$MSA >= 0.8) {
        cat("  Excellent sampling adequacy (KMO ≥ 0.8)\n")
    } else if(kmo_result$MSA >= 0.7) {
        cat("  Good sampling adequacy (KMO ≥ 0.7)\n")
    } else if(kmo_result$MSA >= 0.6) {
        cat("  Adequate sampling adequacy (KMO ≥ 0.6)\n")
    } else {
        cat("  Poor sampling adequacy (KMO < 0.6)\n")
    }
    
    # === 6. AVERAGE EXPRESSION VALIDATION ===
    cat("\n=== 6. AVERAGE EXPRESSION VALIDATION ===\n")
    
    # Calculate average expression
    ribo_avg <- rowMeans(ribo_df)
    
    # Distance correlations between average and individual genes
    avg_dcor_correlations <- sapply(ribosomal_genes, function(gene) {
        dcor(ribo_avg, ribo_df[, gene])
    })
    
    cat("Distance Correlations between Average Ribo and Individual Genes:\n")
    for(i in 1:length(avg_dcor_correlations)) {
        cat("  ", names(avg_dcor_correlations)[i], ":", round(avg_dcor_correlations[i], 3), "\n")
    }
    
    cat("\nDistance Correlation Average Summary:\n")
    cat("  Mean correlation with average:", round(mean(avg_dcor_correlations), 3), "\n")
    cat("  Min correlation with average:", round(min(avg_dcor_correlations), 3), "\n")
    cat("  All correlations > 0.7:", all(avg_dcor_correlations > 0.7), "\n")
    cat("  All correlations > 0.8:", all(avg_dcor_correlations > 0.8), "\n")
    
    # === 7. VISUALIZATION PLOTS ===
    cat("\n=== 7. CREATING VISUALIZATION PLOTS ===\n")
    
    # Scree plot only
    p3 <- fviz_eig(pca_result, 
                   addlabels = TRUE,
                   title = paste("PCA Scree Plot -", sample_name))
    
    # Export Scree plot as PDF (4x4 inches)
    pdf_filename <- paste0("PCA_Scree_Plot_", gsub("HYW_", "", sample_name), ".pdf")
    pdf(pdf_filename, width = 4, height = 4)
    print(p3)
    dev.off()
    cat("Scree plot exported as:", pdf_filename, "\n")
    
    # Display Scree plot in console
    print(p3)
    
    # === 8. OVERALL ASSESSMENT ===
    cat("\n=== 8. OVERALL ASSESSMENT FOR AVERAGING ===\n")
    
    assessment_score <- 0
    max_score <- 7
    
    # Store assessment details for Excel export
    assessment_details <- c()
    
    # Criteria for averaging suitability
    if(mean(dcor_values) > 0.6) {
        detail_line <- paste0("✓ High mean distance correlation ( ", round(mean(dcor_values), 3), " > 0.6)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ Low mean distance correlation ( ", round(mean(dcor_values), 3), " ≤ 0.6)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(var_explained[1] > 60) {
        detail_line <- paste0("✓ PC1 explains substantial variance ( ", round(var_explained[1], 1), "% > 60%)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ PC1 explains insufficient variance ( ", round(var_explained[1], 1), "% ≤ 60%)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(pc1_signs_same) {
        detail_line <- "✓ All genes load with same sign on PC1 (unidirectional)"
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- "✗ Genes load with different signs on PC1 (bidirectional)"
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(cronbach_result$total$std.alpha > 0.8) {
        detail_line <- paste0("✓ High internal consistency (α = ", round(cronbach_result$total$std.alpha, 3), " > 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ Low internal consistency (α = ", round(cronbach_result$total$std.alpha, 3), " ≤ 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(mcdonald_omega > 0.8) {
        detail_line <- paste0("✓ High reliability (ω = ", round(mcdonald_omega, 3), " > 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ Low reliability (ω = ", round(mcdonald_omega, 3), " ≤ 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(kmo_result$MSA > 0.7) {
        detail_line <- paste0("✓ Good sampling adequacy (KMO = ", round(kmo_result$MSA, 3), " > 0.7)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ Poor sampling adequacy (KMO = ", round(kmo_result$MSA, 3), " ≤ 0.7)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    if(mean(avg_dcor_correlations) > 0.8) {
        detail_line <- paste0("✓ High distance correlation with average ( ", round(mean(avg_dcor_correlations), 3), " > 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
        assessment_score <- assessment_score + 1
    } else {
        detail_line <- paste0("✗ Low distance correlation with average ( ", round(mean(avg_dcor_correlations), 3), " ≤ 0.8)")
        cat(detail_line, "\n")
        assessment_details <- c(assessment_details, detail_line)
    }
    
    cat("\nOVERALL ASSESSMENT SCORE:", assessment_score, "/", max_score, "\n")
    
    if(assessment_score >= 6) {
        recommendation <- "✓ EXCELLENT - Averaging is highly appropriate"
        cat("RECOMMENDATION:", recommendation, "\n")
    } else if(assessment_score >= 4) {
        recommendation <- "✓ GOOD - Averaging is appropriate"
        cat("RECOMMENDATION:", recommendation, "\n")
    } else if(assessment_score >= 2) {
        recommendation <- "⚠ CAUTION - Averaging may introduce bias"
        cat("RECOMMENDATION:", recommendation, "\n")
    } else {
        recommendation <- "✗ POOR - Averaging not recommended"
        cat("RECOMMENDATION:", recommendation, "\n")
    }
    
    # Return results
    results <- list(
        sample = sample_name,
        descriptive_stats = desc_stats,
        correlations = list(
            distance = dcor_matrix
        ),
        pca = pca_result,
        variance_explained = var_explained,
        pc1_loadings = pc1_loadings,
        cronbach_alpha = cronbach_result$total$std.alpha,
        mcdonald_omega = mcdonald_omega,
        aiic = aiic,
        kmo = kmo_result$MSA,
        avg_dcor_correlations = avg_dcor_correlations,
        assessment_score = assessment_score,
        max_score = max_score,
        assessment_details = assessment_details,
        recommendation = recommendation
    )
    
    return(results)
}

# Function to export console output to Excel
export_console_to_excel <- function(sample_name, results, ribosomal_genes) {
    # Create a comprehensive text output
    console_text <- paste0(
        "RIBOSOMAL GENE VALIDATION ANALYSIS FOR: ", sample_name, "\n",
        "================================================================================\n\n",
        
        "Sample: ", sample_name, "\n",
        "Number of cells: ", nrow(as.data.frame(t(GetAssayData(prostate_results$seurat_obj, assay = "RNA", slot = "data")[ribosomal_genes, grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE)]))), "\n",
        "Ribosomal genes analyzed: ", paste(ribosomal_genes, collapse = ", "), "\n\n",
        
        "=== 1. DESCRIPTIVE STATISTICS ===\n",
        paste(capture.output(print(results$descriptive_stats)), collapse = "\n"), "\n\n",
        
        "=== 2. DISTANCE CORRELATION ANALYSIS ===\n",
        "Distance Correlation Matrix:\n",
        paste(capture.output(print(round(results$correlations$distance, 3))), collapse = "\n"), "\n\n",
        
        "=== 3. PRINCIPAL COMPONENT ANALYSIS ===\n",
        "PCA Summary:\n",
        paste(capture.output(print(summary(results$pca)$importance)), collapse = "\n"), "\n\n",
        
        "Variance Explained by Components:\n"
    )
    
    # Add variance explained details
    for(i in 1:min(4, length(results$variance_explained))) {
        console_text <- paste0(console_text, "  PC", i, ": ", round(results$variance_explained[i], 2), "%\n")
    }
    
    console_text <- paste0(console_text, "\nPC1 Loadings (Gene Contributions):\n")
    pc1_df <- data.frame(
        Gene = names(results$pc1_loadings),
        Loading = round(results$pc1_loadings, 3),
        Abs_Loading = round(abs(results$pc1_loadings), 3)
    )
    pc1_df <- pc1_df[order(pc1_df$Abs_Loading, decreasing = TRUE), ]
    console_text <- paste0(console_text, paste(capture.output(print(pc1_df)), collapse = "\n"), "\n\n")
    
    console_text <- paste0(console_text,
        "=== 4. RELIABILITY ANALYSIS ===\n",
        "Cronbach's Alpha: ", round(results$cronbach_alpha, 3), "\n",
        "McDonald's Omega (ω): ", round(results$mcdonald_omega, 3), "\n",
        "Average Inter-Item Correlation (AIIC): ", round(results$aiic, 3), "\n\n",
        
        "=== 5. KAISER-MEYER-OLKIN (KMO) MEASURE ===\n",
        "Overall KMO: ", round(results$kmo, 3), "\n\n",
        
        "=== 6. AVERAGE EXPRESSION VALIDATION ===\n",
        "Distance Correlations between Average Ribo and Individual Genes:\n"
    )
    
    for(i in 1:length(results$avg_dcor_correlations)) {
        console_text <- paste0(console_text, "  ", names(results$avg_dcor_correlations)[i], ": ", round(results$avg_dcor_correlations[i], 3), "\n")
    }
    
    console_text <- paste0(console_text, "\n",
        "Mean distance correlation with average: ", round(mean(results$avg_dcor_correlations), 3), "\n",
        "Min distance correlation with average: ", round(min(results$avg_dcor_correlations), 3), "\n\n",
        
        "=== 8. OVERALL ASSESSMENT FOR AVERAGING ===\n"
    )
    
    # Add detailed assessment lines
    for(detail in results$assessment_details) {
        console_text <- paste0(console_text, detail, "\n")
    }
    
    console_text <- paste0(console_text, "\n",
        "OVERALL ASSESSMENT SCORE: ", results$assessment_score, "/", results$max_score, "\n",
        "RECOMMENDATION: ", results$recommendation, "\n"
    )
    
    # Split text into lines for Excel
    text_lines <- strsplit(console_text, "\n")[[1]]
    
    # Create data frame for Excel
    excel_data <- data.frame(
        Console_Output = text_lines,
        stringsAsFactors = FALSE
    )
    
    return(excel_data)
}

# Function to run validation across all samples
validate_all_samples <- function() {
    cat("\n", rep("=", 100), "\n")
    cat("RUNNING RIBOSOMAL GENE VALIDATION ACROSS ALL SAMPLES")
    cat("\n", rep("=", 100), "\n")
    
    ribosomal_genes <- c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26")
    all_validation_results <- list()
    
    # Create Excel workbook for console outputs
    wb_console <- createWorkbook()
    
    for(sample_name in tumor_samples) {
        tryCatch({
            validation_result <- validate_ribosomal_genes(sample_name, ribosomal_genes)
            if(!is.null(validation_result)) {
                all_validation_results[[sample_name]] <- validation_result
                
                # Export console output to Excel
                console_data <- export_console_to_excel(sample_name, validation_result, ribosomal_genes)
                sheet_name <- gsub("HYW_", "", sample_name)  # Remove HYW_ prefix
                addWorksheet(wb_console, sheet_name)
                writeData(wb_console, sheet_name, console_data)
                
                # Auto-adjust column width
                setColWidths(wb_console, sheet_name, cols = 1, widths = 100)
                
                cat("Completed validation for:", sample_name, "\n")
            }
        }, error = function(e) {
            cat("Error in validation for", sample_name, ":", conditionMessage(e), "\n")
        })
    }
    
    # Save console outputs Excel file
    console_filename <- "Ribosomal_Gene_Validation_Console_Outputs.xlsx"
    saveWorkbook(wb_console, console_filename, overwrite = TRUE)
    cat("\nConsole outputs exported to:", console_filename, "\n")
    
    # Create summary table
    cat("\n", rep("=", 80), "\n")
    cat("SUMMARY TABLE ACROSS ALL SAMPLES")
    cat("\n", rep("=", 80), "\n")
    
    summary_df <- data.frame(
        Sample = character(),
        Mean_Distance_Cor = numeric(),
        PC1_Variance_Pct = numeric(),
        Cronbach_Alpha = numeric(),
        McDonald_Omega = numeric(),
        AIIC = numeric(),
        KMO = numeric(),
        Mean_Distance_Avg_Correlation = numeric(),
        Min_Distance_Avg_Correlation = numeric(),
        Assessment_Score = character(),
        Recommendation = character(),
        stringsAsFactors = FALSE
    )
    
    for(sample_name in names(all_validation_results)) {
        result <- all_validation_results[[sample_name]]
        
        # Calculate mean correlations
        dcor_matrix <- result$correlations$distance
        upper_tri <- upper.tri(dcor_matrix)
        
        mean_dcor <- mean(dcor_matrix[upper_tri])
        
        # Get recommendation
        score <- result$assessment_score
        max_score <- result$max_score
        score_text <- paste(score, "/", max_score)
        
        # Clean recommendation text
        recommendation_clean <- gsub("✓ |⚠ |✗ ", "", result$recommendation)
        recommendation_clean <- gsub(" - .*", "", recommendation_clean)
        
        summary_df <- rbind(summary_df, data.frame(
            Sample = sample_name,
            Mean_Distance_Cor = round(mean_dcor, 3),
            PC1_Variance_Pct = round(result$variance_explained[1], 1),
            Cronbach_Alpha = round(result$cronbach_alpha, 3),
            McDonald_Omega = round(result$mcdonald_omega, 3),
            AIIC = round(result$aiic, 3),
            KMO = round(result$kmo, 3),
            Mean_Distance_Avg_Correlation = round(mean(result$avg_dcor_correlations), 3),
            Min_Distance_Avg_Correlation = round(min(result$avg_dcor_correlations), 3),
            Assessment_Score = score_text,
            Recommendation = recommendation_clean,
            stringsAsFactors = FALSE
        ))
    }
    
    print(summary_df)
    
    # Export summary table to Excel
    wb_summary <- createWorkbook()
    addWorksheet(wb_summary, "Summary")
    writeData(wb_summary, "Summary", summary_df)
    
    # Format summary sheet
    addStyle(wb_summary, "Summary", 
             style = createStyle(textDecoration = "bold"), 
             rows = 1, cols = 1:ncol(summary_df))
    setColWidths(wb_summary, "Summary", cols = 1:ncol(summary_df), widths = "auto")
    
    summary_filename <- "Ribosomal_Gene_Validation_Summary.xlsx"
    saveWorkbook(wb_summary, summary_filename, overwrite = TRUE)
    cat("\nSummary table exported to:", summary_filename, "\n")
    
    cat("\nFiles created:\n")
    cat("1. Individual PCA Scree plots: PCA_Scree_Plot_[sample].pdf (4x4 inches each)\n")
    cat("2. Console outputs:", console_filename, "(one sheet per sample)\n")
    cat("3. Summary table:", summary_filename, "\n")
    
    return(list(
        individual_results = all_validation_results,
        summary = summary_df
    ))
}

# Run complete validation across all samples
validation_results <- validate_all_samples()

# Run the validation
cat("\nStarting comprehensive ribosomal gene validation analysis...\n")

# Usage instructions
cat("\nValidation analysis code is ready to run!")
cat("\nUse: validation_results <- validate_all_samples()")
cat("\nOr for individual samples: validate_ribosomal_genes('HYW_4701_Tumor')")