# =========================================================================================
# Part 3.10: TRDE (purple) vs non-TRDE (gray) GAM Plots, and Explanatory Power Values
# =========================================================================================
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
  
  # Sort cells by their explanatory power (from highest to lowest)
  sorted_indices <- order(explanatory_power, decreasing = TRUE)
  sorted_cell_names <- rownames(sample_data)[sorted_indices]
  
  # Target number of cells (model_dev_explained of total cells, rounded)
  target_cells <- round(nrow(sample_data) * model_dev_explained)
  
  # Take the top cells by explanatory power
  deviance_cells <- sorted_cell_names[1:target_cells]
  non_deviance_cells <- sorted_cell_names[(target_cells+1):length(sorted_cell_names)]
  
  # Calculate dynamic axis limits with a small margin
  x_range <- range(sample_data$TRPM4, na.rm = TRUE)
  y_range <- range(sample_data$Expression, na.rm = TRUE)
  x_margin <- 0.05 * diff(x_range)
  y_margin <- 0.05 * diff(y_range)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(
    TRPM4 = sample_data$TRPM4,
    Expression = sample_data$Expression,
    Group = ifelse(rownames(sample_data) %in% deviance_cells, "Dev explained", "Non-dev explained")
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
      subtitle = paste("Dev explained: ", dev_explained_percentage, "%", sep = ""),
      x = "TRPM4 Expression",
      y = "Ribo Expression"
    ) +
    # Dynamically set axis limits with small margins
    scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
    scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
  
  # Display the plot
  print(p)
  
  return(p)
}

# List of samples to process
tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                   "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                   "HYW_5755_Tumor")

# Generate scatter plots for all samples
cat("\n\nGenerating scatter plots for all samples...\n")
cat(rep("=", 60), "\n")

# Create a list to store all plots
all_plots <- list()

# Process each sample for plotting
for (sample_name in tumor_samples) {
  tryCatch({
    plot_obj <- create_sample_scatter_plot(sample_name)
    if (!is.null(plot_obj)) {
      all_plots[[sample_name]] <- plot_obj
    }
  }, error = function(e) {
    cat("Error creating plot for", sample_name, ":", conditionMessage(e), "\n")
  })
}

cat("\nCompleted generating scatter plots for all samples.\n")


# Function to extract and compile cell-level data including TRDE and non-TRDE explanatory power values
extract_cell_data_for_export <- function(current_sample) {
    cat("\nExtracting data for sample:", current_sample, "\n")
    
    # Get sample data for current sample
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    
    # Calculate position vs inflection for current sample
    ip_value <- inflection_points[current_sample]
    position_vs_inflection <- ifelse(sample_data$TRPM4 < ip_value, "Before", "After")
    
    # Get PCa cells using the correct method 
    pca_clusters <- c(6, 9, 11, 14, 19) 
    
    # Correctly identify cluster cells and sample cells
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    # Subset the data for PCa cells in the current sample
    sample_subset <- subset(prostate_ca_seurat, cells = selected_cells)
    
    # Get the sample's GAM model results
    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    
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
    
    # Sort cells by their explanatory power (from highest to lowest)
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    sorted_cell_names <- rownames(sample_data)[sorted_indices]
    
    # Target number of cells (model_dev_explained of total cells, rounded)
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    
    # Take the top cells by explanatory power
    deviance_cells <- sorted_cell_names[1:target_cells]
    non_deviance_cells <- sorted_cell_names[(target_cells+1):length(sorted_cell_names)]
    
    # Extract patient ID from sample name (assuming format like "HYW_4701_Tumor")
    patient_id <- gsub("_Tumor$", "", current_sample)
    
    # Create the output dataframe
    output_df <- data.frame(
        Patient_ID = patient_id,
        Cell_Index = rownames(sample_data),  # This should contain the full barcode like "HYW_4701_Tumor_ATCGGGGATTT..."
        EP_Value = explanatory_power,
        Sorted_Index = match(rownames(sample_data), sorted_cell_names),  # Position in sorted order
        Sorted_Cell_Names = rownames(sample_data),  # Same as Cell_Index but renamed for clarity
        Deviance_Cells = ifelse(rownames(sample_data) %in% deviance_cells, "YES", "NO"),
        Non_Deviance_Cells = ifelse(rownames(sample_data) %in% non_deviance_cells, "YES", "NO"),
        TRPM4_Expression = sample_data$TRPM4,
        Ribo_Expression = sample_data$Expression,
        Model_Dev_Explained = model_dev_explained,
        stringsAsFactors = FALSE
    )
    
    # Sort by explanatory power (highest first)
    output_df <- output_df[order(output_df$EP_Value, decreasing = TRUE), ]
    
    cat("Extracted data for", nrow(output_df), "cells from sample", current_sample, "\n")
    
    return(output_df)
}

# Extract data for all samples
cat("\n\nExtracting cell-level data for Excel export...\n")
cat(rep("=", 60), "\n")

# Create a list to store all sample data
all_sample_data <- list()

# Process each sample
for (sample_name in tumor_samples) {
    tryCatch({
        sample_df <- extract_cell_data_for_export(sample_name)
        if (!is.null(sample_df)) {
            all_sample_data[[sample_name]] <- sample_df
        }
    }, error = function(e) {
        cat("Error extracting data for", sample_name, ":", conditionMessage(e), "\n")
    })
}

# Combine all sample data into one dataframe
if (length(all_sample_data) > 0) {
    combined_data <- do.call(rbind, all_sample_data)
    rownames(combined_data) <- NULL  # Reset row names
    
    # Create Excel file with multiple sheets
    wb <- createWorkbook()
    
    # Add a combined sheet with all samples
    addWorksheet(wb, "All_Samples_Combined")
    writeData(wb, "All_Samples_Combined", combined_data)
    
    # Add individual sheets for each sample
    for (sample_name in names(all_sample_data)) {
        sheet_name <- gsub("HYW_", "", sample_name)  # Shorten sheet name
        sheet_name <- gsub("_Tumor", "", sheet_name)
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, all_sample_data[[sample_name]])
    }
    
    # Add a summary sheet
    summary_data <- data.frame(
        Sample = names(all_sample_data),
        Patient_ID = sapply(names(all_sample_data), function(x) gsub("_Tumor$", "", x)),
        Total_Cells = sapply(all_sample_data, nrow),
        Mean_EP_Value = sapply(all_sample_data, function(x) round(mean(x$EP_Value, na.rm = TRUE), 4)),
        Max_EP_Value = sapply(all_sample_data, function(x) round(max(x$EP_Value, na.rm = TRUE), 4)),
        Min_EP_Value = sapply(all_sample_data, function(x) round(min(x$EP_Value, na.rm = TRUE), 4)),
        Deviance_Cells_Count = sapply(all_sample_data, function(x) sum(x$Deviance_Cells == "YES")),
        Model_Dev_Explained = sapply(all_sample_data, function(x) round(unique(x$Model_Dev_Explained), 4)),
        stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", summary_data)
    
    # Create filename with timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("single_cell_analysis_results_", timestamp, ".xlsx")
    
    # Save the Excel file
    saveWorkbook(wb, filename, overwrite = TRUE)
    
    cat("\n")
    cat(rep("=", 60), "\n")
    cat("Excel file exported successfully!\n")
    cat("Filename:", filename, "\n")
    cat("Total samples processed:", length(all_sample_data), "\n")
    cat("Total cells exported:", nrow(combined_data), "\n")
    cat("\nFile contains the following sheets:\n")
    cat("- All_Samples_Combined: All cells from all samples\n")
    for (sample_name in names(all_sample_data)) {
        short_name <- gsub("HYW_", "", gsub("_Tumor", "", sample_name))
        cat("-", short_name, ": Individual sample data\n")
    }
    cat("- Summary: Sample-level statistics\n")
    cat(rep("=", 60), "\n")
    
} else {
    cat("No data was successfully extracted. Please check for errors above.\n")
}

# Optional: Display first few rows of the combined data to verify
if (exists("combined_data") && nrow(combined_data) > 0) {
    cat("\nFirst 10 rows of combined data (preview):\n")
    print(head(combined_data, 10))
    
    cat("\nColumn summary:\n")
    cat("- Patient_ID: Patient identifier (e.g., HYW_4701)\n")
    cat("- Cell_Index: Full cell barcode from scRNA-seq dataset\n") 
    cat("- EP_Value: Explanatory power value for each cell\n")
    cat("- Sorted_Index: Position in sorted order (by EP_Value)\n")
    cat("- Sorted_Cell_Names: Cell names in sorted order\n")
    cat("- Deviance_Cells: YES if cell contributes to deviance explained\n")
    cat("- Non_Deviance_Cells: YES if cell does not contribute to deviance explained\n")
    cat("- TRPM4_Expression: TRPM4 gene expression level\n")
    cat("- Ribo_Expression: Ribosomal gene expression level\n")
    cat("- Model_Dev_Explained: Overall model deviance explained for the sample\n")
}
