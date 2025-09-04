# ========================================================
# Family Distribution Analysis with gam.check() and AIC/BIC
# ========================================================
library(mgcv)
library(writexl)
library(dplyr)

# Function to create improved gam.check plots and save as PDF
create_improved_gam_plots <- function(gam_model, sample_name, gene_set_name, output_dir = ".") {
    if (is.null(gam_model)) {
        return(NULL)
    }
    
    # Create filename for PDF
    safe_sample_name <- gsub("[^A-Za-z0-9_]", "_", sample_name)
    safe_gene_set <- gsub("[^A-Za-z0-9_]", "_", gene_set_name)
    pdf_filename <- file.path(output_dir, paste0("GAM_Diagnostics_", safe_sample_name, "_", safe_gene_set, ".pdf"))
    
    # Open PDF device
    pdf(pdf_filename, width = 12, height = 10)
    
    # Set up 2x2 layout
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 3, 0))
    
    # Extract data for plots
    residuals_dev <- residuals(gam_model, type = "deviance")
    fitted_vals <- fitted(gam_model)
    response_vals <- gam_model$model[[1]]
    
    # 1. QQ Plot for residuals
    qqnorm(residuals_dev, main = "Q-Q Plot of Deviance Residuals", 
           pch = 16, col = "#303030", cex = 0.8)
    qqline(residuals_dev, col = "#6666FF", lwd = 2)
    
    # 2. Histogram of residuals
    hist(residuals_dev, main = "Histogram of Deviance Residuals", 
         xlab = "Deviance Residuals", col = "#C0C0C0", border = "#000000",
         breaks = 20, freq = FALSE)
    # Add normal curve overlay
    x_seq <- seq(min(residuals_dev), max(residuals_dev), length.out = 100)
    lines(x_seq, dnorm(x_seq, mean(residuals_dev), sd(residuals_dev)), 
          col = "#6666FF", lwd = 2)
    
    # 3. Residuals vs Fitted Values
    plot(fitted_vals, residuals_dev, 
         main = "Residuals vs Fitted Values",
         xlab = "Fitted Values", ylab = "Deviance Residuals",
         pch = 16, col = "#303030", cex = 0.8)
    abline(h = 0, col = "#1EE17F", lwd = 2, lty = 2)
    # Add smoothed line
    lines(lowess(fitted_vals, residuals_dev), col = "#6666FF", lwd = 2)
    
    # 4. Response vs Fitted Values
    plot(fitted_vals, response_vals,
         main = "Response vs Fitted Values",
         xlab = "Fitted Values", ylab = "Response",
         pch = 16, col = "#303030", cex = 0.8)
    # Add 1:1 line
    abline(a = 0, b = 1, col = "#1EE17F", lwd = 2, lty = 2)
    # Add smoothed line
    lines(lowess(fitted_vals, response_vals), col = "#6666FF", lwd = 2)
    
    # Add overall title
    main_title <- paste("GAM Diagnostic Plots:", sample_name, "-", gene_set_name)
    mtext(main_title, outer = TRUE, cex = 1.3, font = 2)
    
    # Close PDF device
    dev.off()
    
    cat("Diagnostic plots saved to:", pdf_filename, "\n")
    return(pdf_filename)
}

# Function to test alternative families
test_alternative_families <- function(gam_data, sample_name, gene_set_name, original_model) {
    
    families_to_test <- list(
        Gaussian = gaussian(),
        QuasiPoisson = quasipoisson(),
        Gamma = Gamma(link = "log"),
        Inverse_Gaussian = inverse.gaussian(link = "1/mu^2"),
        Negative_Binomial = nb()
    )
    
    family_results <- list()
    
    for (family_name in names(families_to_test)) {
        cat("Testing", family_name, "family for", sample_name, "-", gene_set_name, "...\n")
        
        tryCatch({
            # Prepare data for specific families
            test_data <- gam_data
            
            if (family_name %in% c("Gamma", "QuasiPoisson") && any(test_data$Expression <= 0)) {
                test_data$Expression <- pmax(test_data$Expression, 0.001)
            }
            
            # For Inverse Gaussian, ensure positive values
            if (family_name == "Inverse_Gaussian" && any(test_data$Expression <= 0)) {
                test_data$Expression <- pmax(test_data$Expression, 0.001)
            }
            
            # Use the same model structure as the original
            original_k <- original_model$smooth[[1]]$bs.dim
            
            # Fit model
            test_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = original_k), 
                              family = families_to_test[[family_name]], 
                              data = test_data, 
                              method = "REML")
            
            # Extract metrics
            aic_val <- AIC(test_model)
            bic_val <- BIC(test_model)
            dev_expl <- summary(test_model)$dev.expl * 100
            
            family_results[[family_name]] <- data.frame(
                Sample = sample_name,
                Gene_Set = gene_set_name,
                Family = family_name,
                AIC = round(aic_val, 2),
                BIC = round(bic_val, 2),
                Dev_Explained_Percent = round(dev_expl, 2),
                stringsAsFactors = FALSE
            )
            
        }, error = function(e) {
            cat(" Error with", family_name, ":", conditionMessage(e), "\n")
            family_results[[family_name]] <- data.frame(
                Sample = sample_name,
                Gene_Set = gene_set_name,
                Family = family_name,
                AIC = NA,
                BIC = NA,
                Dev_Explained_Percent = NA,
                stringsAsFactors = FALSE
            )
        })
    }
    
    return(do.call(rbind, family_results))
}

# Function to identify best family for each model
identify_best_families <- function(family_comparison_df) {
    if (nrow(family_comparison_df) == 0) return(list(AIC = data.frame(), BIC = data.frame()))
    
    # Filter for models with valid AIC and BIC
    valid_models <- family_comparison_df[!is.na(family_comparison_df$AIC) & 
                                         !is.na(family_comparison_df$BIC), ]
    
    if (nrow(valid_models) == 0) return(list(AIC = data.frame(), BIC = data.frame()))
    
    # Find best family by AIC
    best_families_aic <- valid_models %>%
        group_by(Sample, Gene_Set, Sample_Type) %>%
        slice_min(AIC, n = 1) %>%
        ungroup()
    
    # Find best family by BIC
    best_families_bic <- valid_models %>%
        group_by(Sample, Gene_Set, Sample_Type) %>%
        slice_min(BIC, n = 1) %>%
        ungroup()
    
    # Create AIC-based dataframe
    result_aic <- data.frame(
        Sample = best_families_aic$Sample,
        Gene_Set = best_families_aic$Gene_Set,
        Sample_Type = best_families_aic$Sample_Type,
        Best_Family = best_families_aic$Family,
        Best_AIC = best_families_aic$AIC,
        Best_Dev_Explained = best_families_aic$Dev_Explained_Percent,
        stringsAsFactors = FALSE
    )
    
    # Create BIC-based dataframe
    result_bic <- data.frame(
        Sample = best_families_bic$Sample,
        Gene_Set = best_families_bic$Gene_Set,
        Sample_Type = best_families_bic$Sample_Type,
        Best_Family = best_families_bic$Family,
        Best_BIC = best_families_bic$BIC,
        Best_Dev_Explained = best_families_bic$Dev_Explained_Percent,
        stringsAsFactors = FALSE
    )
    
    return(list(AIC = result_aic, BIC = result_bic))
}

# Main function to run comprehensive diagnostics on all models
run_comprehensive_gam_diagnostics <- function(results_object, sample_type = "PCa", create_plots = TRUE, plot_dir = "GAM_Diagnostic_Plots") {
    
    all_family_comparisons <- list()
    plot_files <- list()
    
    # Create plot directory if requested
    if (create_plots && !dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
        cat("Created directory for plots:", plot_dir, "\n")
    }
    
    for (sample_name in names(results_object)) {
        sample_results <- results_object[[sample_name]]
        
        if (is.null(sample_results)) next
        
        for (gene_set_name in names(sample_results)) {
            gene_set_result <- sample_results[[gene_set_name]]
            
            if (is.null(gene_set_result) || is.null(gene_set_result$best_model)) next
            
            cat("\nRunning diagnostics for", sample_name, "-", gene_set_name, "\n")
            
            # Create diagnostic plots if requested
            if (create_plots) {
                plot_file <- create_improved_gam_plots(
                    gene_set_result$best_model, 
                    sample_name, 
                    gene_set_name, 
                    plot_dir
                )
                if (!is.null(plot_file)) {
                    plot_files[[paste(sample_name, gene_set_name, sep="_")]] <- plot_file
                }
            }
            
            # Test alternative families if data is available
            if (!is.null(gene_set_result$gam_data)) {
                family_comparison <- test_alternative_families(
                    gene_set_result$gam_data, 
                    sample_name, 
                    gene_set_name, 
                    gene_set_result$best_model
                )
                family_comparison$Sample_Type <- sample_type
                all_family_comparisons[[paste(sample_name, gene_set_name, sep="_")]] <- family_comparison
            }
        }
    }
    
    # Combine results
    family_comparison_df <- if(length(all_family_comparisons) > 0) {
        do.call(rbind, all_family_comparisons)
    } else {
        data.frame()
    }
    
    return(list(
        family_comparisons = family_comparison_df,
        plot_files = plot_files
    ))
}

# ========================================================
# Usage
# ========================================================

# Run comprehensive diagnostics on the results
cat("Running comprehensive GAM diagnostics with gam.check()...\n")

# Analyze PCa results
pca_diagnostics <- run_comprehensive_gam_diagnostics(pca_results, "PCa", create_plots = TRUE)

# Analyze Non-Ca results  
non_ca_diagnostics <- run_comprehensive_gam_diagnostics(non_ca_results, "Non-Ca", create_plots = TRUE)

# Combine all results
all_family_comparisons <- rbind(pca_diagnostics$family_comparisons, non_ca_diagnostics$family_comparisons)
all_plot_files <- c(pca_diagnostics$plot_files, non_ca_diagnostics$plot_files)

# Identify best families
best_families <- identify_best_families(all_family_comparisons)
best_families_aic <- best_families$AIC
best_families_bic <- best_families$BIC

# Export comprehensive results to Excel
export_list <- list(
    Family_Comparisons = all_family_comparisons[, c("Sample", "Gene_Set", "Family", "AIC", "BIC", "Dev_Explained_Percent")],
    Best_Families_by_AIC = best_families_aic,
    Best_Family_of_Data_Distribution_by_BIC = best_families_bic
)

# Add plot files information to export
plot_files_df <- data.frame(
    Sample_Gene_Set = names(all_plot_files),
    Plot_File_Path = unlist(all_plot_files),
    stringsAsFactors = FALSE
)
export_list$Diagnostic_Plot_Files <- plot_files_df

# Export to Excel
tryCatch({
    write_xlsx(export_list, path = "GAM_Family_Analysis_to_Determine_if_Gaussian_or_Other_Distribution_is_Appropriate.xlsx")
    cat("\nGAM family analysis exported to: GAM_Family_Analysis_to_Determine_if_Gaussian_or_Other_Distribution_is_Appropriate.xlsx\n")
    cat("Sheets included:\n")
    cat("- Family_Comparisons: Performance comparison across different families\n")
    cat("- Best_Families_by_AIC: Best-performing family for each model by AIC\n")
    cat("- Best_Family_of_Data_Distribution_by_BIC: Best-performing family for each model by BIC\n")
    cat("- Diagnostic_Plot_Files: List of generated PDF plot files\n")
    
    cat("\nDiagnostic plots saved to folder: GAM_Diagnostic_Plots/\n")
    cat("Each model has its own PDF file with 4 diagnostic plots.\n")
}, error = function(e) {
    cat("Error exporting results:", conditionMessage(e), "\n")
})

cat("\nDiagnostic analysis complete!\n")