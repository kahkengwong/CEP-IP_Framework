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