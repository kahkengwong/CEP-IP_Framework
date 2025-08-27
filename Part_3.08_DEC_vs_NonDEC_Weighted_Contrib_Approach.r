# ====================================================================
# Advanced Statistical Analysis with Proper Statistical Tests
# ====================================================================
library(dplyr)
library(openxlsx)
library(ggplot2)
library(gridExtra)

# Function to calculate comprehensive metrics for a sample (unchanged)
calculate_sample_metrics <- function(current_sample) {
    cat("Processing sample:", current_sample, "\n")
    
    # Get PCa cells (same as original code)
    pca_clusters <- c(6, 9, 11, 14, 19) 
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    # Get GAM model results
    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    model_dev_explained <- summary(sample_gam_model)$dev.expl
    
    # Match cell IDs
    sample_subset <- subset(prostate_ca_seurat, cells = selected_cells)
    common_cells <- intersect(rownames(sample_data), colnames(sample_subset))
    
    if (length(common_cells) == 0) {
        return(NULL)
    }
    
    # Filter data
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    # Calculate models and explanatory power
    null_model <- gam(Expression ~ 1, data = sample_data)
    null_fitted <- fitted(null_model)
    model_fitted <- fitted(sample_gam_model)
    null_residuals <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted
    null_sq_diff <- null_residuals^2
    model_sq_diff <- model_residuals^2
    explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
    
    # Calculate individual contributions to total deviance
    individual_contributions <- null_sq_diff - model_sq_diff
    
    return(list(
        sample = current_sample,
        model_dev_explained = model_dev_explained,
        null_sq_diff = null_sq_diff,
        model_sq_diff = model_sq_diff,
        explanatory_power = explanatory_power,
        individual_contributions = individual_contributions,
        sample_data = sample_data
    ))
}

# Process all samples (unchanged)
tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                   "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                   "HYW_5755_Tumor")

all_sample_data <- list()
for (sample_name in tumor_samples) {
    tryCatch({
        metrics <- calculate_sample_metrics(sample_name)
        if (!is.null(metrics)) {
            all_sample_data[[sample_name]] <- metrics
        }
    }, error = function(e) {
        cat("Error processing", sample_name, ":", conditionMessage(e), "\n")
    })
}

# ====================================================================
# 1. PROPERLY CORRECTED DISTRIBUTIONAL ANALYSIS
# ====================================================================
cat("\n=== PROPERLY CORRECTED DISTRIBUTIONAL ANALYSIS ===\n")

distributional_results <- list()

for (sample_name in names(all_sample_data)) {
    data <- all_sample_data[[sample_name]]
    
    # Calculate weighted average of individual EPs
    weights <- data$null_sq_diff / sum(data$null_sq_diff)
    weighted_mean_ep <- sum(data$explanatory_power * weights)
    
    # Simple average of individual EPs
    simple_mean_ep <- mean(data$explanatory_power)
    
    # Median of individual EPs (what's being plotted)
    median_ep <- median(data$explanatory_power)
    
    # Actual GAM DE
    actual_de <- data$model_dev_explained
    
    # Test statistics - comparing what's actually being plotted
    median_diff <- median_ep - actual_de
    
    # Theoretical reconstruction test
    theoretical_de <- sum(data$individual_contributions) / sum(data$null_sq_diff)
    reconstruction_error <- abs(theoretical_de - actual_de)
    
    # THE REAL ISSUE: We should test if the OBSERVED median differs from EXPECTED median under null
    # The current test is asking: "Do the individual EP values have a median significantly different from GAM DE?"
    # But this is NOT the right question! We should ask: "Is the difference between median EP and GAM DE significant?"
    
    # PROPER Statistical Test: Two-proportion test or chi-square goodness of fit
    # Alternative: Bootstrap confidence interval for the median difference
    ep_values <- data$explanatory_power
    n_cells <- length(ep_values)
    
    # Check if distribution is normal for diagnostic purposes
    if (n_cells >= 3) {
        shapiro_test <- shapiro.test(ep_values)
        is_normal <- shapiro_test$p.value > 0.05
        shapiro_p <- shapiro_test$p.value
    } else {
        is_normal <- FALSE
        shapiro_p <- NA
    }
    
    # CORRECT APPROACH: Bootstrap confidence interval for median difference
    set.seed(123)  # For reproducibility
    n_bootstrap <- 1000
    bootstrap_medians <- replicate(n_bootstrap, {
        boot_sample <- sample(ep_values, replace = TRUE)
        median(boot_sample)
    })
    
    # Calculate 95% confidence interval for median
    median_ci <- quantile(bootstrap_medians, c(0.025, 0.975))
    
    # Test if GAM DE falls within the 95% CI of the median
    gam_in_ci <- actual_de >= median_ci[1] && actual_de <= median_ci[2]
    
    # Calculate p-value: proportion of bootstrap medians more extreme than observed difference
    observed_diff <- abs(median_ep - actual_de)
    bootstrap_diffs <- abs(bootstrap_medians - actual_de)
    p_value_bootstrap <- mean(bootstrap_diffs >= observed_diff)
    
    # Alternative: Simple z-test for median (approximation for large samples)
    # Standard error of median ≈ 1.253 * sd(x) / sqrt(n) for normal distributions
    if (n_cells >= 30) {
        se_median <- 1.253 * sd(ep_values) / sqrt(n_cells)
        z_score <- (median_ep - actual_de) / se_median
        p_value_z <- 2 * (1 - pnorm(abs(z_score)))  # Two-tailed test
        test_method <- "Z-test for median difference"
        test_statistic <- z_score
        p_value <- p_value_z
    } else {
        # For small samples, use bootstrap p-value
        test_method <- "Bootstrap test for median difference"
        test_statistic <- observed_diff
        p_value <- p_value_bootstrap
    }
    
    # Calculate effect size (standardized difference)
    effect_size <- abs(median_ep - actual_de) / sd(ep_values)
    
    distributional_results[[sample_name]] <- data.frame(
        Sample = sample_name,
        N_Cells = n_cells,
        GAM_DE = actual_de,
        Median_EP = median_ep,
        Mean_EP = simple_mean_ep,
        Weighted_Mean_EP = weighted_mean_ep,
        Median_Difference = median_diff,
        Median_CI_Lower = median_ci[1],
        Median_CI_Upper = median_ci[2],
        GAM_in_CI = gam_in_ci,
        Theoretical_DE = theoretical_de,
        Reconstruction_Error = reconstruction_error,
        Distribution_Type = ifelse(is_normal, "Normal", "Non-normal"),
        Shapiro_P_Value = shapiro_p,
        Test_Method = test_method,
        Test_Statistic = as.numeric(test_statistic),
        P_Value = p_value,
        P_Value_Bootstrap = p_value_bootstrap,
        Effect_Size = effect_size,
        stringsAsFactors = FALSE
    )
}

distributional_df <- do.call(rbind, distributional_results)
distributional_df$FDR_Adjusted_P <- p.adjust(distributional_df$P_Value, method = "fdr")

# ====================================================================
# 2. EP CLASSIFICATION VALIDITY ANALYSIS
# ====================================================================
cat("\n=== EP CLASSIFICATION VALIDITY ANALYSIS ===\n")

validity_results <- list()

for (sample_name in names(all_sample_data)) {
    data <- all_sample_data[[sample_name]]
    
    # Get the original classification (using the same logic as main analysis)
    sorted_indices <- order(data$explanatory_power, decreasing = TRUE)
    target_cells <- round(length(data$explanatory_power) * data$model_dev_explained)
    selected_cells <- sorted_indices[1:target_cells]
    
    # Compare selected vs non-selected cells
    selected_ep <- data$explanatory_power[selected_cells]
    non_selected_ep <- data$explanatory_power[-selected_cells]
    selected_contrib <- data$individual_contributions[selected_cells]
    non_selected_contrib <- data$individual_contributions[-selected_cells]
    
    # EP differences
    ep_diff_test <- wilcox.test(selected_ep, non_selected_ep, alternative = "greater")
    ep_effect_size <- (mean(selected_ep) - mean(non_selected_ep)) / sd(data$explanatory_power)
    
    # Contribution differences  
    contrib_diff_test <- wilcox.test(selected_contrib, non_selected_contrib, alternative = "greater")
    contrib_effect_size <- (mean(selected_contrib) - mean(non_selected_contrib)) / sd(data$individual_contributions)
    
    # Proportion of total improvement captured by selected cells
    total_improvement <- sum(data$individual_contributions)
    selected_improvement <- sum(selected_contrib)
    improvement_proportion <- selected_improvement / total_improvement
    
    # Expected vs observed improvement (should match GAM DE)
    improvement_error <- abs(improvement_proportion - data$model_dev_explained)
    
    # Calculate overlap between dev-explained and non-dev-explained groups
    dev_min <- min(selected_ep)
    dev_max <- max(selected_ep)
    non_dev_min <- min(non_selected_ep)
    non_dev_max <- max(non_selected_ep)
    
    # Calculate overlap range
    overlap_min <- max(dev_min, non_dev_min)
    overlap_max <- min(dev_max, non_dev_max)
    
    # If overlap_max > overlap_min, there is overlap
    if (overlap_max > overlap_min) {
        overlap_range <- overlap_max - overlap_min
        total_range <- max(dev_max, non_dev_max) - min(dev_min, non_dev_min)
        overlap_percentage <- (overlap_range / total_range) * 100
    } else {
        overlap_percentage <- 0
    }
    
    validity_results[[sample_name]] <- data.frame(
        Sample = sample_name,
        Selected_Cells = length(selected_cells),
        Non_Selected_Cells = length(non_selected_ep),
        Selection_Proportion = length(selected_cells) / length(data$explanatory_power),
        Selected_EP_Mean = mean(selected_ep),
        Non_Selected_EP_Mean = mean(non_selected_ep),
        EP_Difference = mean(selected_ep) - mean(non_selected_ep),
        EP_Test_P = ep_diff_test$p.value,
        EP_Effect_Size = ep_effect_size,
        Selected_Contrib_Mean = mean(selected_contrib),
        Non_Selected_Contrib_Mean = mean(non_selected_contrib),
        Contrib_Difference = mean(selected_contrib) - mean(non_selected_contrib),
        Contrib_Test_P = contrib_diff_test$p.value,
        Contrib_Effect_Size = contrib_effect_size,
        Improvement_Proportion = improvement_proportion,
        Expected_Proportion = data$model_dev_explained,
        Improvement_Error = improvement_error,
        Overlap_Percentage = overlap_percentage,
        stringsAsFactors = FALSE
    )
}

validity_df <- do.call(rbind, validity_results)
validity_df$EP_Test_FDR <- p.adjust(validity_df$EP_Test_P, method = "fdr")
validity_df$Contrib_Test_FDR <- p.adjust(validity_df$Contrib_Test_P, method = "fdr")

# ====================================================================
# CORRECTED VISUALIZATION
# ====================================================================
cat("\n=== CREATING CORRECTED VISUALIZATIONS ===\n")

# Clean sample names for plotting
distributional_df$Sample_Clean <- gsub("HYW_", "", gsub("_Tumor", "", distributional_df$Sample))
validity_df$Sample_Clean <- gsub("HYW_", "", gsub("_Tumor", "", validity_df$Sample))

# Plot 1: GAM DE vs Median Individual EP
dist_plot_data <- data.frame(
    Sample = distributional_df$Sample_Clean,
    GAM_DE = distributional_df$GAM_DE * 100,  # Convert to percentage
    Median_EP = distributional_df$Median_EP * 100,  # Convert to percentage
    Reconstruction_Error = distributional_df$Reconstruction_Error * 100,
    Significance = ifelse(distributional_df$FDR_Adjusted_P < 0.001, "***",
                          ifelse(distributional_df$FDR_Adjusted_P < 0.01, "**",
                                 ifelse(distributional_df$FDR_Adjusted_P < 0.05, "*", "ns"))),
    N_Cells = distributional_df$N_Cells,
    Test_Method = distributional_df$Test_Method,
    GAM_in_CI = distributional_df$GAM_in_CI
)

p1 <- ggplot(dist_plot_data, aes(x = Sample)) +
    geom_col(aes(y = GAM_DE), fill = "#2CA02C", alpha = 0.7, width = 0.4, position = position_nudge(x = -0.2)) +
    geom_col(aes(y = Median_EP), fill = "#4B0082", alpha = 0.7, width = 0.4, position = position_nudge(x = 0.2)) +
    geom_text(aes(y = pmax(GAM_DE, Median_EP) + 5, label = Significance), 
              size = 4, fontface = "bold") +
    labs(title = "GAM Deviance Explained vs Median Individual EP",
         subtitle = "Purple = GAM DE, Green = Median EP | CORRECTED: Bootstrap/Z-test for median difference | ** = FDR < 0.01",
         x = "Sample", y = "Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

print(p1) 

# Plot 2: EP Classification Validity with FDR < 0.01 cutoff
validity_plot_data <- data.frame()
for (sample_name in names(all_sample_data)) {
    data <- all_sample_data[[sample_name]]
    sample_clean <- gsub("HYW_", "", gsub("_Tumor", "", sample_name))
    
    # Get classification (same logic as main analysis)
    sorted_indices <- order(data$explanatory_power, decreasing = TRUE)
    target_cells <- round(length(data$explanatory_power) * data$model_dev_explained)
    selected_cells <- sorted_indices[1:target_cells]
    
    # Create data for selected cells
    selected_data <- data.frame(
        Sample = sample_clean,
        EP_Value = data$explanatory_power[selected_cells],
        Group = "Selected (Dev-explained)",
        GAM_DE = data$model_dev_explained * 100
    )
    
    # Create data for non-selected cells  
    non_selected_data <- data.frame(
        Sample = sample_clean,
        EP_Value = data$explanatory_power[-selected_cells],
        Group = "Non-selected",
        GAM_DE = data$model_dev_explained * 100
    )
    
    validity_plot_data <- rbind(validity_plot_data, selected_data, non_selected_data)
}

# Get significance annotations with FDR < 0.01 cutoff
val_sig_values <- ifelse(validity_df$EP_Test_FDR < 0.001, "***",
                         ifelse(validity_df$EP_Test_FDR < 0.01, "**",
                                ifelse(validity_df$EP_Test_FDR < 0.05, "*", "ns")))

p2 <- ggplot(validity_plot_data, aes(x = Sample, y = EP_Value)) +
    geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("Selected (Dev-explained)" = "#4B0082", 
                                 "Non-selected" = "#C0C0C0")) +
    geom_text(data = data.frame(Sample = validity_df$Sample_Clean, 
                                EP_Value = rep(1.2, nrow(validity_df)),
                                sig = val_sig_values),
              aes(x = Sample, y = EP_Value, label = sig), 
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    labs(title = "EP Classification Validity",
         subtitle = "Selected (Dev-explained) vs Non-selected cells | Tests if classification captures high-EP cells | ** = FDR < 0.01",
         x = "Sample", y = "Explanatory Power", fill = "Classification") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom") +
    scale_y_continuous(limits = c(-1, 1.4))

print(p2)

# Plot 3: Individual EP values with FDR < 0.01 cutoff
all_ep_plot_data <- data.frame()
for (sample_name in names(all_sample_data)) {
    data <- all_sample_data[[sample_name]]
    sample_clean <- gsub("HYW_", "", gsub("_Tumor", "", sample_name))
    
    # Get classification (same logic as main analysis)
    sorted_indices <- order(data$explanatory_power, decreasing = TRUE)
    target_cells <- round(length(data$explanatory_power) * data$model_dev_explained)
    selected_cells <- sorted_indices[1:target_cells]
    
    temp_data <- data.frame(
        Sample = sample_clean,
        EP_Value = data$explanatory_power,
        Classification = ifelse(seq_along(data$explanatory_power) %in% selected_cells, 
                                "Dev-explained", "Non-dev-explained"),
        GAM_DE = round(data$model_dev_explained * 100, 1)
    )
    all_ep_plot_data <- rbind(all_ep_plot_data, temp_data)
}

# Create sample labels with GAM DE percentages
all_ep_plot_data$Sample_Label <- paste0(all_ep_plot_data$Sample, "\n(DE: ", 
                                        all_ep_plot_data$GAM_DE, "%)")

p3 <- ggplot(all_ep_plot_data, aes(x = Sample, y = EP_Value)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6, color = "black", fill = "lightgray") +
    geom_jitter(aes(color = Classification), width = 0.3, alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("Dev-explained" = "#4B0082", 
                                  "Non-dev-explained" = "#C0C0C0")) +
    labs(title = "Individual Cell EP Values by GAM Classification",
         subtitle = "Purple = Dev-explained cells, Gray = Non-dev-explained cells",
         x = "Sample (GAM Deviance Explained %)", 
         y = "Explanatory Power", 
         color = "Classification") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom") +
    scale_y_continuous(limits = c(-10, 1.5))

print(p3)

# ====================================================================
# EXPORT TO EXCEL
# ====================================================================
cat("\n=== EXPORTING TO EXCEL ===\n")

wb <- createWorkbook()

# Sheet 1: Corrected Distributional Analysis
addWorksheet(wb, "1_Distributional_Analysis")
writeData(wb, "1_Distributional_Analysis", distributional_df)

# Sheet 2: EP Classification Validity Analysis
addWorksheet(wb, "2_EP_Classification_Validity")
writeData(wb, "2_EP_Classification_Validity", validity_df)

# Save Excel file
excel_filename <- "Cleaned_Advanced_GAM_Analysis.xlsx"
saveWorkbook(wb, excel_filename, overwrite = TRUE)

# Print summary results
cat("\n", rep("=", 80), "\n")
cat("CLEANED ADVANCED STATISTICAL ANALYSIS SUMMARY (FDR < 0.01)\n")
cat(rep("=", 80), "\n")

cat("\n1. DISTRIBUTIONAL ANALYSIS (CORRECTED):\n")
cat("Samples with significant Median EP vs GAM DE differences (FDR < 0.01):", sum(distributional_df$FDR_Adjusted_P < 0.01), "/", nrow(distributional_df), "\n")
cat("Mean reconstruction error:", round(mean(distributional_df$Reconstruction_Error), 4), "\n")
cat("Sample sizes range:", min(distributional_df$N_Cells), "-", max(distributional_df$N_Cells), "cells\n")

# Show the actual values being compared for clarity
cat("\nDetailed comparison (GAM DE vs Median EP):\n")
for(i in 1:nrow(distributional_df)) {
    significance_status <- ifelse(distributional_df$FDR_Adjusted_P[i] < 0.01, "Significant (FDR < 0.01)", "Non-significant")
    cat(sprintf("%-10s: GAM DE = %5.1f%%, Median EP = %5.1f%%, Diff = %5.1f%%, P = %.4f (%s)\n",
                distributional_df$Sample_Clean[i],
                distributional_df$GAM_DE[i] * 100,
                distributional_df$Median_EP[i] * 100,
                distributional_df$Median_Difference[i] * 100,
                distributional_df$P_Value[i],
                significance_status))
}

cat("\n2. EP CLASSIFICATION VALIDITY ANALYSIS:\n")
cat("Samples with significant EP differences (selected vs non-selected, FDR < 0.01):", sum(validity_df$EP_Test_FDR < 0.01), "/", nrow(validity_df), "\n")
cat("Samples with significant contribution differences (FDR < 0.01):", sum(validity_df$Contrib_Test_FDR < 0.01), "/", nrow(validity_df), "\n")
cat("Mean improvement proportion error:", round(mean(validity_df$Improvement_Error), 4), "\n")
cat("Mean overlap percentage:", round(mean(validity_df$Overlap_Percentage), 1), "%\n")

cat("\nCleaned Excel file exported as:", excel_filename, "\n")
cat("File contains 2 analysis sheets (cell heterogeneity analysis removed).\n")
cat("Three visualization plots have been displayed with FDR < 0.01 significance cutoff.\n")

# Additional diagnostic information
cat("\n=== DIAGNOSTIC INFORMATION ===\n")
cat("Statistical test explanation for GAM DE vs Median EP comparison:\n")
cat("- SIGNIFICANCE (FDR < 0.01): The difference between GAM Deviance Explained and\n")
cat("  the median of individual cell explanatory powers is statistically significant.\n")
cat("  This suggests the GAM summary statistic may not represent the typical cell experience.\n")
cat("- NON-SIGNIFICANCE: The difference is not statistically significant, suggesting\n")
cat("  the GAM DE is a reasonable representation of the median individual cell EP.\n")

# Export all plots as SVG 375 x 345
# Filename: EP1_GAM_vs_EP


# ====================================================================
# 3.09 DEC (purple) vs Non-DEC (gray) Cells Scatter Plot: Weighted Contribution Approach
# ====================================================================
library(dplyr)
library(tibble)

# Function to create scatter plot for a sample
create_sample_scatter_plot <- function(current_sample) {
  cat("\nCreating scatter plot for sample:", current_sample, "\n")
  
  # Get PCa cells
  pca_clusters <- c(6, 9, 11, 14, 19) 
  
  # Identify cluster cells and sample cells
  cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
  sample_cells <- WhichCells(prostate_results$seurat_obj, 
                             cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
  selected_cells <- intersect(cluster_cells, sample_cells)
  
  # Subset the data for PCa cells in the current sample
  sample_subset <- subset(prostate_ca_seurat, cells = selected_cells)
  
  # Get the sample's GAM model results
  sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
  sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
  
  # Get the overall deviance explained from the model
  model_dev_explained <- summary(sample_gam_model)$dev.expl
  
  # Match cell IDs between sample_data and sample_subset
  common_cells <- intersect(rownames(sample_data), colnames(sample_subset))
  
  if (length(common_cells) == 0) {
    cat("No common cells found between GAM data and Seurat object for", current_sample, "\n")
    return(NULL)
  }
  
  # Filter sample_data to include only cells present in the Seurat object
  sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
  
  # Calculate the null model (intercept only)
  null_model <- gam(Expression ~ 1, data = sample_data)
  
  # Calculate fitted values and residuals for both models
  null_fitted <- fitted(null_model)
  model_fitted <- fitted(sample_gam_model)
  null_residuals <- sample_data$Expression - null_fitted
  model_residuals <- sample_data$Expression - model_fitted
  
  # Calculate squared residuals (individual contributions to deviance)
  null_sq_diff <- null_residuals^2
  model_sq_diff <- model_residuals^2
  
  # Define individual explanatory power
  explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
  
  # CORRECTED APPROACH: Calculate weighted contributions
  # Each cell's weight in the overall deviance explained calculation
  null_weights <- null_sq_diff / sum(null_sq_diff)
  
  # Each cell's weighted contribution to deviance explained
  weighted_contribution <- explanatory_power * null_weights
  
  # Sort cells by their weighted contribution (not raw explanatory power)
  sorted_indices <- order(weighted_contribution, decreasing = TRUE)
  sorted_contributions <- weighted_contribution[sorted_indices]
  sorted_cell_names <- rownames(sample_data)[sorted_indices]
  
  # Find cells that collectively explain the observed deviance
  cumulative_contribution <- cumsum(sorted_contributions)
  
  # Find the threshold index where cumulative contribution reaches model_dev_explained
  threshold_index <- which(cumulative_contribution >= model_dev_explained)[1]
  
  # Handle edge case where no cells reach the threshold
  if (is.na(threshold_index)) {
    threshold_index <- length(sorted_cell_names)
  }
  
  # Select cells based on cumulative weighted contribution
  deviance_cells <- sorted_cell_names[1:threshold_index]
  non_deviance_cells <- sorted_cell_names[(threshold_index+1):length(sorted_cell_names)]
  
  # Calculate dynamic axis limits with a small margin
  x_range <- range(sample_data$TRPM4, na.rm = TRUE)
  y_range <- range(sample_data$Expression, na.rm = TRUE)
  x_margin <- 0.05 * diff(x_range)
  y_margin <- 0.05 * diff(y_range)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(
    TRPM4 = sample_data$TRPM4,
    Expression = sample_data$Expression,
    Group = ifelse(rownames(sample_data) %in% deviance_cells, "Dev explained", "Non-dev explained"),
    explanatory_power = explanatory_power,
    weighted_contribution = weighted_contribution
  )
  
  # Add a drawing order column to control which points appear on top
  plot_data$draw_order <- ifelse(plot_data$Group == "Dev explained", 2, 1)
  
  # Sort the data frame by the draw order
  plot_data <- plot_data[order(plot_data$draw_order), ]
  
  # Get the GAM model from existing analysis
  gam_model <- sample_gam_model
  
  # Create prediction data for the GAM line
  pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
  pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
  pred_data$fit <- pred$fit
  pred_data$se.fit <- pred$se.fit
  
  # Calculate percentage for subtitle
  dev_explained_percentage <- round(model_dev_explained * 100, 2)
  
  # Calculate actual percentage of cells classified as "dev explained"
  actual_dev_cells_pct <- round(length(deviance_cells) / nrow(sample_data) * 100, 1)
  
  # Create the plot with ggplot2
  p <- ggplot() +
    # Add points with colors and transparency
    geom_point(data = plot_data, 
               aes(x = TRPM4, y = Expression, color = Group),
               size = 1.8, alpha = 0.3) +
    # Add the GAM line
    geom_line(data = pred_data,
              aes(x = TRPM4, y = fit),
              color = "#FFCC99",  # Orange
              size = 1.2) +
    # Add confidence interval ribbon
    geom_ribbon(data = pred_data,
                aes(x = TRPM4, 
                    ymin = fit - 1.96 * se.fit, 
                    ymax = fit + 1.96 * se.fit),
                fill = "#FFCC99", 
                alpha = 0.2) +
    # Set colors for the groups
    scale_color_manual(values = c("Dev explained" = "#4B0082",  # Purple
                                 "Non-dev explained" = "#C0C0C0")) +  # Gray
    # Styling
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_line(color = "#F5F5F5"),
      legend.position = "none", 
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      # Add thin black border around the plot area
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_blank()
    ) +
    # Set title and axis labels with sample name and dev_explained percentage
    labs(
      title = paste("TRPM4 vs Ribo Expression (", current_sample, ")", sep = ""),
      subtitle = paste("Dev explained: ", dev_explained_percentage, "% | Purple cells: ", actual_dev_cells_pct, "%", sep = ""),
      x = "TRPM4 Expression",
      y = "Ribo Expression"
    ) +
    # Dynamically set axis limits with small margins
    scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
    scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
  
  # Print diagnostic information
  cat("Sample:", current_sample, "\n")
  cat("Model deviance explained:", round(model_dev_explained * 100, 2), "%\n")
  cat("Cells classified as 'dev explained':", length(deviance_cells), "out of", nrow(sample_data), 
      "(", actual_dev_cells_pct, "%)\n")
  cat("Cumulative contribution of selected cells:", round(sum(weighted_contribution[rownames(sample_data) %in% deviance_cells]) * 100, 2), "%\n\n")
  
  # Display the plot
  print(p)
  
  # Return additional diagnostic information
  return(list(
    plot = p,
    dev_cells_count = length(deviance_cells),
    total_cells = nrow(sample_data),
    dev_cells_percentage = actual_dev_cells_pct,
    model_dev_explained = model_dev_explained * 100,
    cumulative_contribution_check = sum(weighted_contribution[rownames(sample_data) %in% deviance_cells])
  ))
}

# List of samples to process
tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                   "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                   "HYW_5755_Tumor")

# Generate scatter plots for all samples
cat("\n\nGenerating scatter plots for all samples...\n")
cat(rep("=", 60), "\n")

# Create a list to store all plots and diagnostics
all_plots <- list()
all_diagnostics <- list()

# Process each sample for plotting
for (sample_name in tumor_samples) {
  tryCatch({
    result <- create_sample_scatter_plot(sample_name)
    if (!is.null(result)) {
      all_plots[[sample_name]] <- result$plot
      all_diagnostics[[sample_name]] <- result[names(result) != "plot"]
    }
  }, error = function(e) {
    cat("Error creating plot for", sample_name, ":", conditionMessage(e), "\n")
  })
}

cat("\nCompleted generating scatter plots for all samples.\n")

# Print summary diagnostics
cat("\n", rep("=", 60), "\n")
cat("SUMMARY DIAGNOSTICS\n")
cat(rep("=", 60), "\n")
for (sample_name in names(all_diagnostics)) {
  diag <- all_diagnostics[[sample_name]]
  cat(sprintf("%-20s | Model DE: %5.1f%% | Purple cells: %5.1f%% | Cumulative check: %5.1f%%\n", 
              sample_name, 
              diag$model_dev_explained, 
              diag$dev_cells_percentage, 
              diag$cumulative_contribution_check * 100))

}
