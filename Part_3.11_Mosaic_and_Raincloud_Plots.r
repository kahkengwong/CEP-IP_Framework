###################################################################################
# Part 3.11: Mosaic and Raincloud Plots of Ribo Expression Pre-IP and Post-IP
###################################################################################

# ====================================================================
# Extract Contingency Table Data for Mosaic Plots using Direct Inflection Points
# ====================================================================
# Initialize lists to store contingency tables
all_gray_tables <- list()
all_purple_tables <- list()

# Define patient IDs and their corresponding inflection points
patient_ids <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                 "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                 "HYW_5755_Tumor")

# Define inflection points for each patient
inflection_points <- c(3.800, 2.214, 3.179, 3.306, 2.636, 3.465, 3.476)
names(inflection_points) <- patient_ids

cat("Using the following inflection points:\n")
for (i in 1:length(patient_ids)) {
  cat("Pt.", i, "(", patient_ids[i], "):", inflection_points[i], "\n")
}

# Function to extract contingency data for a sample using direct inflection points
extract_contingency_data <- function(current_sample) {
  cat("\nExtracting contingency data for sample:", current_sample, "\n")
  
  # Get PCa cells (same as in scatter plot and DEG analysis code)
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
  
  # Define individual explanatory power (same as DEG analysis)
  explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
  
  # Sort cells by their explanatory power (from highest to lowest)
  sorted_indices <- order(explanatory_power, decreasing = TRUE)
  sorted_cell_names <- rownames(sample_data)[sorted_indices]
  
  # Target number of cells (model_dev_explained of total cells, rounded)
  target_cells <- round(nrow(sample_data) * model_dev_explained)
  
  # Take the top cells by explanatory power
  deviance_cells <- sorted_cell_names[1:target_cells]
  non_deviance_cells <- sorted_cell_names[(target_cells+1):length(sorted_cell_names)]
  
  # Create cell type classification
  sample_data$cell_type <- ifelse(rownames(sample_data) %in% deviance_cells, "purple", "gray")
  
  # Use the direct inflection point for this sample
  ip_value <- inflection_points[current_sample]
  cat("Using inflection point:", ip_value, "for sample:", current_sample, "\n")
  
  # Classify cells as before/after inflection point (same logic as DEG analysis)
  sample_data$timing <- ifelse(sample_data$TRPM4 < ip_value, "Before", "After")
  
  # For the "Below/Above" classification, use the GAM fitted line as threshold
  # This makes biological sense - cells below the fitted curve vs above the fitted curve
  fitted_values <- fitted(sample_gam_model)
  sample_data$location <- ifelse(sample_data$Expression <= fitted_values, "Below", "Above")
  
  # Create contingency tables
  # Gray cells (non-dev explained)
  gray_data <- sample_data[sample_data$cell_type == "gray", ]
  if (nrow(gray_data) > 0) {
    gray_table <- table(gray_data$timing, gray_data$location)
  } else {
    gray_table <- table(factor(c(), levels = c("Before", "After")), 
                        factor(c(), levels = c("Below", "Above")))
  }
  
  # Purple cells (dev explained) 
  purple_data <- sample_data[sample_data$cell_type == "purple", ]
  if (nrow(purple_data) > 0) {
    purple_table <- table(purple_data$timing, purple_data$location)
  } else {
    purple_table <- table(factor(c(), levels = c("Before", "After")), 
                          factor(c(), levels = c("Below", "Above")))
  }
  
  # Ensure tables have the correct structure (2x2 with proper row/column names)
  expected_rows <- c("Before", "After")
  expected_cols <- c("Below", "Above")
  
  # Create properly formatted tables
  gray_formatted <- matrix(0, nrow = 2, ncol = 2, 
                           dimnames = list(expected_rows, expected_cols))
  purple_formatted <- matrix(0, nrow = 2, ncol = 2, 
                             dimnames = list(expected_rows, expected_cols))
  
  # Fill in the values
  for (i in expected_rows) {
    for (j in expected_cols) {
      if (i %in% rownames(gray_table) && j %in% colnames(gray_table)) {
        gray_formatted[i, j] <- gray_table[i, j]
      }
      if (i %in% rownames(purple_table) && j %in% colnames(purple_table)) {
        purple_formatted[i, j] <- purple_table[i, j]
      }
    }
  }
  
  cat("Gray cells distribution:\n")
  cat("  Before IP - Below fitted line:", gray_formatted["Before", "Below"], "\n")
  cat("  Before IP - Above fitted line:", gray_formatted["Before", "Above"], "\n")
  cat("  After IP  - Below fitted line:", gray_formatted["After", "Below"], "\n")
  cat("  After IP  - Above fitted line:", gray_formatted["After", "Above"], "\n")
  
  cat("Purple cells distribution:\n")
  cat("  Before IP - Below fitted line:", purple_formatted["Before", "Below"], "\n")
  cat("  Before IP - Above fitted line:", purple_formatted["Before", "Above"], "\n")
  cat("  After IP  - Below fitted line:", purple_formatted["After", "Below"], "\n")
  cat("  After IP  - Above fitted line:", purple_formatted["After", "Above"], "\n")
  
  # Calculate odds ratios for reporting
  gray_or <- tryCatch({
    (gray_formatted[1,1] * gray_formatted[2,2]) / (gray_formatted[1,2] * gray_formatted[2,1])
  }, error = function(e) NA)
  
  purple_or <- tryCatch({
    (purple_formatted[1,1] * purple_formatted[2,2]) / (purple_formatted[1,2] * purple_formatted[2,1])
  }, error = function(e) NA)
  
  cat("Gray cells odds ratio:", round(gray_or, 3), "\n")
  cat("Purple cells odds ratio:", round(purple_or, 3), "\n")
  
  return(list(gray = gray_formatted, purple = purple_formatted))
}

# Extract contingency data for all patients
cat("\nExtracting contingency tables for all patients...\n")
cat(rep("=", 60), "\n")

for (i in 1:length(patient_ids)) {
  patient_id <- patient_ids[i]
  
  tryCatch({
    contingency_data <- extract_contingency_data(patient_id)
    
    if (!is.null(contingency_data)) {
      all_gray_tables[[i]] <- contingency_data$gray
      all_purple_tables[[i]] <- contingency_data$purple
    } else {
      cat("Warning: No contingency data for", patient_id, ". Creating empty tables.\n")
      # Create empty tables if extraction failed
      all_gray_tables[[i]] <- matrix(c(0, 0, 0, 0), nrow = 2, 
                                     dimnames = list(c("Before", "After"), c("Below", "Above")))
      all_purple_tables[[i]] <- matrix(c(0, 0, 0, 0), nrow = 2, 
                                       dimnames = list(c("Before", "After"), c("Below", "Above")))
    }
  }, error = function(e) {
    cat("Error extracting contingency data for", patient_id, ":", conditionMessage(e), "\n")
    # Create empty tables on error
    all_gray_tables[[i]] <- matrix(c(0, 0, 0, 0), nrow = 2, 
                                   dimnames = list(c("Before", "After"), c("Below", "Above")))
    all_purple_tables[[i]] <- matrix(c(0, 0, 0, 0), nrow = 2, 
                                     dimnames = list(c("Before", "After"), c("Below", "Above")))
  })
}

cat("\nContingency table extraction complete!\n")
cat("Gray tables created:", length(all_gray_tables), "\n")
cat("Purple tables created:", length(all_purple_tables), "\n")

# Verify the tables are properly created
cat("\nVerification - showing table dimensions:\n")
for (i in 1:length(patient_ids)) {
  if (i <= length(all_gray_tables) && i <= length(all_purple_tables)) {
    cat("Patient", i, "(", patient_ids[i], "):\n")
    cat("  Gray table dimensions:", dim(all_gray_tables[[i]]), "\n")
    cat("  Purple table dimensions:", dim(all_purple_tables[[i]]), "\n")
    cat("  Gray table sum:", sum(all_gray_tables[[i]]), "\n")
    cat("  Purple table sum:", sum(all_purple_tables[[i]]), "\n")
  } else {
    cat("Patient", i, "(", patient_ids[i], "): No tables created\n")
  }
}

cat("\nReady for mosaic plot generation\n")


# ============
# Mosaic Plots 
# ============
library(tidyverse)  
library(ggmosaic)   # For mosaic plots
library(svglite)    # For SVG export

# Define the patient IDs
patient_ids <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                "HYW_5755_Tumor")

#------------------------------------------------
# CUSTOM MOSAIC PLOT FUNCTION
#------------------------------------------------

# Function to create a single patient's horizontal mosaic plot
create_horizontal_mosaic <- function(patient_index) {
    pid <- patient_ids[patient_index]
    
    # Get contingency tables
    gray_table <- all_gray_tables[[patient_index]]
    purple_table <- all_purple_tables[[patient_index]]
    
    # Calculate odds ratios
    gray_or <- (gray_table[1,1] * gray_table[2,2]) / (gray_table[1,2] * gray_table[2,1])
    purple_or <- (purple_table[1,1] * purple_table[2,2]) / (purple_table[1,2] * purple_table[2,1])
    
    # Extract counts
    gray_before_below <- gray_table[1,1]
    gray_before_above <- gray_table[1,2]
    gray_after_below <- gray_table[2,1]
    gray_after_above <- gray_table[2,2]
    
    purple_before_below <- purple_table[1,1]
    purple_before_above <- purple_table[1,2]
    purple_after_below <- purple_table[2,1]
    purple_after_above <- purple_table[2,2]
    
    # Calculate totals
    gray_total <- sum(gray_table)
    purple_total <- sum(purple_table)
    
    # Create data frames with proper structure for horizontal mosaic plots
    
    # Gray cells data
    gray_df <- data.frame(
        variable = rep(c("Before", "After"), c(gray_before_below + gray_before_above, 
                                               gray_after_below + gray_after_above)),
        location = c(rep("Below", gray_before_below), 
                     rep("Above", gray_before_above),
                     rep("Below", gray_after_below), 
                     rep("Above", gray_after_above))
    )
    
    # Purple cells data
    purple_df <- data.frame(
        variable = rep(c("Before", "After"), c(purple_before_below + purple_before_above, 
                                               purple_after_below + purple_after_above)),
        location = c(rep("Below", purple_before_below), 
                     rep("Above", purple_before_above),
                     rep("Below", purple_after_below), 
                     rep("Above", purple_after_above))
    )
    
    # Ensure proper factor levels for horizontal layout
    gray_df$variable <- factor(gray_df$variable, levels = c("Before", "After"))
    gray_df$location <- factor(gray_df$location, levels = c("Below", "Above"))
    
    purple_df$variable <- factor(purple_df$variable, levels = c("Before", "After"))
    purple_df$location <- factor(purple_df$location, levels = c("Below", "Above"))
    
    # Calculate totals for each timing group
    gray_before_total <- gray_before_below + gray_before_above
    gray_after_total <- gray_after_below + gray_after_above
    
    gray_labels <- data.frame(
        variable = c("Before", "Before", "After", "After"),
        location = c("Below", "Above", "Below", "Above"),
        count = c(gray_before_below, gray_before_above, gray_after_below, gray_after_above),
        percentage = c(ifelse(gray_before_total > 0, round(gray_before_below/gray_before_total*100), 0),
                       ifelse(gray_before_total > 0, round(gray_before_above/gray_before_total*100), 0),
                       ifelse(gray_after_total > 0, round(gray_after_below/gray_after_total*100), 0),
                       ifelse(gray_after_total > 0, round(gray_after_above/gray_after_total*100), 0))
    )
    
    # Calculate totals for each timing group
    purple_before_total <- purple_before_below + purple_before_above
    purple_after_total <- purple_after_below + purple_after_above
    
    purple_labels <- data.frame(
        variable = c("Before", "Before", "After", "After"),
        location = c("Below", "Above", "Below", "Above"),
        count = c(purple_before_below, purple_before_above, purple_after_below, purple_after_above),
        percentage = c(ifelse(purple_before_total > 0, round(purple_before_below/purple_before_total*100), 0),
                       ifelse(purple_before_total > 0, round(purple_before_above/purple_before_total*100), 0),
                       ifelse(purple_after_total > 0, round(purple_after_below/purple_after_total*100), 0),
                       ifelse(purple_after_total > 0, round(purple_after_above/purple_after_total*100), 0))
    )
    
    # Make sure variable and location are factors with the correct levels
    gray_labels$variable <- factor(gray_labels$variable, levels = c("Before", "After"))
    gray_labels$location <- factor(gray_labels$location, levels = c("Below", "Above"))
    
    purple_labels$variable <- factor(purple_labels$variable, levels = c("Before", "After"))
    purple_labels$location <- factor(purple_labels$location, levels = c("Below", "Above"))
    
    # Create gray mosaic plot with horizontal layout
    gray_plot <- ggplot(data = gray_df) +
        geom_mosaic(aes(x = product(variable), fill = location), offset = 0.01) +
        geom_mosaic_text(aes(x = product(variable), fill = location), 
                         na.rm = TRUE, show.legend = FALSE) +
        scale_fill_manual(values = c("Below" = "#E6E6E6", "Above" = "#808080")) +
        labs(
            title = paste(pid, "- Gray Cells"),
            subtitle = paste("OR:", round(gray_or, 2)),
            x = "Timing",
            fill = "Location"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        )
    
    # Create purple mosaic plot with horizontal layout
    purple_plot <- ggplot(data = purple_df) +
        geom_mosaic(aes(x = product(variable), fill = location), offset = 0.01) +
        geom_mosaic_text(aes(x = product(variable), fill = location), 
                         na.rm = TRUE, show.legend = FALSE) +
        scale_fill_manual(values = c("Below" = "#E8E1F2", "Above" = "#9B72CF")) +
        labs(
            title = paste(pid, "- Purple Cells"),
            subtitle = paste("OR:", round(purple_or, 2)),
            x = "Timing",
            fill = "Location"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        )
    
    # Add custom labels for counts and percentages
    gray_plot_with_labels <- gray_plot +
        geom_text(data = gray_labels,
                  aes(x = ifelse(variable == "Before", 0.25, 0.75),
                      y = ifelse(location == "Below", 0.25, 0.75),
                      label = paste0(count, "\n(", percentage, "%)")),
                  size = 4, fontface = "bold")
    
    purple_plot_with_labels <- purple_plot +
        geom_text(data = purple_labels,
                  aes(x = ifelse(variable == "Before", 0.25, 0.75),
                      y = ifelse(location == "Below", 0.25, 0.75),
                      label = paste0(count, "\n(", percentage, "%)")),
                  size = 4, fontface = "bold")
    
    # Update axis labels
    gray_plot_final <- gray_plot_with_labels +
        scale_x_continuous(
            breaks = c(0.25, 0.75),
            labels = c("Before Inflection", "After Inflection")
        ) +
        scale_y_continuous(breaks = NULL)
    
    purple_plot_final <- purple_plot_with_labels +
        scale_x_continuous(
            breaks = c(0.25, 0.75),
            labels = c("Before Inflection", "After Inflection")
        ) +
        scale_y_continuous(breaks = NULL)
    
    # Create a combined plot with gray on top and purple on bottom
    combined_plot <- cowplot::plot_grid(
        gray_plot_final, 
        purple_plot_final,
        ncol = 1,
        align = "v",
        rel_heights = c(1, 1)
    )
    
    return(list(
        gray_plot = gray_plot_final,
        purple_plot = purple_plot_final,
        combined_plot = combined_plot
    ))
}

#------------------------------------------------
# ALTERNATIVE APPROACH USING MANUAL POSITIONING
#------------------------------------------------

# Function to create a properly positioned horizontal mosaic
create_manual_horizontal_mosaic <- function(patient_index) {
    pid <- patient_ids[patient_index]
    
    # Get contingency tables
    gray_table <- all_gray_tables[[patient_index]]
    purple_table <- all_purple_tables[[patient_index]]
    
    # Calculate odds ratios
    gray_or <- (gray_table[1,1] * gray_table[2,2]) / (gray_table[1,2] * gray_table[2,1])
    purple_or <- (purple_table[1,1] * purple_table[2,2]) / (purple_table[1,2] * purple_table[2,1])
    
    # Extract counts
    gray_before_below <- gray_table[1,1]
    gray_before_above <- gray_table[1,2]
    gray_after_below <- gray_table[2,1]
    gray_after_above <- gray_table[2,2]
    
    purple_before_below <- purple_table[1,1]
    purple_before_above <- purple_table[1,2]
    purple_after_below <- purple_table[2,1]
    purple_after_above <- purple_table[2,2]
    
    # Calculate totals
    gray_total <- sum(gray_table)
    purple_total <- sum(purple_table)
    
    # Calculate proportions for horizontal mosaic
    # For Gray cells
    gray_before_total <- gray_before_below + gray_before_above
    gray_after_total <- gray_after_below + gray_after_above
    
    # Width proportions (horizontal axis for Before/After)
    gray_before_width <- 0.5
    gray_after_width <- 0.5
    
    # Height proportions within each section
    if (gray_before_total > 0) {
        gray_before_below_height <- gray_before_below / gray_before_total
        gray_before_above_height <- gray_before_above / gray_before_total
    } else {
        gray_before_below_height <- 0.5
        gray_before_above_height <- 0.5
    }
    
    if (gray_after_total > 0) {
        gray_after_below_height <- gray_after_below / gray_after_total
        gray_after_above_height <- gray_after_above / gray_after_total
    } else {
        gray_after_below_height <- 0.5
        gray_after_above_height <- 0.5
    }
    
    # For Purple cells
    purple_before_total <- purple_before_below + purple_before_above
    purple_after_total <- purple_after_below + purple_after_above
    
    # Width proportions (horizontal axis for Before/After)
    purple_before_width <- 0.5
    purple_after_width <- 0.5
    
    # Height proportions within each section
    if (purple_before_total > 0) {
        purple_before_below_height <- purple_before_below / purple_before_total
        purple_before_above_height <- purple_before_above / purple_before_total
    } else {
        purple_before_below_height <- 0.5
        purple_before_above_height <- 0.5
    }
    
    if (purple_after_total > 0) {
        purple_after_below_height <- purple_after_below / purple_after_total
        purple_after_above_height <- purple_after_above / purple_after_total
    } else {
        purple_after_below_height <- 0.5
        purple_after_above_height <- 0.5
    }
    
    # Create rectangle data for gray mosaic
    gray_rects <- data.frame(
        timing = c("Before", "Before", "After", "After"),
        location = c("Below", "Above", "Below", "Above"),
        xmin = c(0, 0, gray_before_width, gray_before_width),
        xmax = c(gray_before_width, gray_before_width, 1, 1),
        ymin = c(0, gray_before_below_height, 0, gray_after_below_height),
        ymax = c(gray_before_below_height, 1, gray_after_below_height, 1),
        count = c(gray_before_below, gray_before_above, gray_after_below, gray_after_above),
        percentage = c(
            ifelse(gray_before_total > 0, round(gray_before_below/gray_before_total*100), 0),
            ifelse(gray_before_total > 0, round(gray_before_above/gray_before_total*100), 0),
            ifelse(gray_after_total > 0, round(gray_after_below/gray_after_total*100), 0),
            ifelse(gray_after_total > 0, round(gray_after_above/gray_after_total*100), 0)
        ),
        label = c(
            paste0(gray_before_below, "\n(", ifelse(gray_before_total > 0, round(gray_before_below/gray_before_total*100), 0), "%)"),
            paste0(gray_before_above, "\n(", ifelse(gray_before_total > 0, round(gray_before_above/gray_before_total*100), 0), "%)"),
            paste0(gray_after_below, "\n(", ifelse(gray_after_total > 0, round(gray_after_below/gray_after_total*100), 0), "%)"),
            paste0(gray_after_above, "\n(", ifelse(gray_after_total > 0, round(gray_after_above/gray_after_total*100), 0), "%)")
        )
    )
    # Create rectangle data for purple mosaic
    purple_rects <- data.frame(
        timing = c("Before", "Before", "After", "After"),
        location = c("Below", "Above", "Below", "Above"),
        xmin = c(0, 0, purple_before_width, purple_before_width),
        xmax = c(purple_before_width, purple_before_width, 1, 1),
        ymin = c(0, purple_before_below_height, 0, purple_after_below_height),
        ymax = c(purple_before_below_height, 1, purple_after_below_height, 1),
        count = c(purple_before_below, purple_before_above, purple_after_below, purple_after_above),
		percentage = c(
            ifelse(purple_before_total > 0, round(purple_before_below/purple_before_total*100), 0),
            ifelse(purple_before_total > 0, round(purple_before_above/purple_before_total*100), 0),
            ifelse(purple_after_total > 0, round(purple_after_below/purple_after_total*100), 0),
            ifelse(purple_after_total > 0, round(purple_after_above/purple_after_total*100), 0)
        ),
        label = c(
            paste0(purple_before_below, "\n(", ifelse(purple_before_total > 0, round(purple_before_below/purple_before_total*100), 0), "%)"),
            paste0(purple_before_above, "\n(", ifelse(purple_before_total > 0, round(purple_before_above/purple_before_total*100), 0), "%)"),
            paste0(purple_after_below, "\n(", ifelse(purple_after_total > 0, round(purple_after_below/purple_after_total*100), 0), "%)"),
            paste0(purple_after_above, "\n(", ifelse(purple_after_total > 0, round(purple_after_above/purple_after_total*100), 0), "%)")
        )
    )
    # Calculate positions for labels
    gray_rects$label_x <- (gray_rects$xmin + gray_rects$xmax) / 2
    gray_rects$label_y <- (gray_rects$ymin + gray_rects$ymax) / 2
    
    purple_rects$label_x <- (purple_rects$xmin + purple_rects$xmax) / 2
    purple_rects$label_y <- (purple_rects$ymin + purple_rects$ymax) / 2
    
    # Create the gray mosaic plot
    gray_plot <- ggplot() +
        # Add rectangles
        geom_rect(data = gray_rects,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = location),
                  color = "white", size = 0.5) +
        # Add labels
        geom_text(data = gray_rects,
                  aes(x = label_x, y = label_y, label = label),
                  size = 4, fontface = "bold") +
        # Set colors
        scale_fill_manual(values = c("Below" = "#E6E6E6", "Above" = "#808080")) +
        # Add titles
        labs(
            title = paste(pid, "- Gray Cells"),
            subtitle = paste("OR:", round(gray_or, 2))
        ) +
        # Set axis scales
        scale_x_continuous(
            breaks = c(0.25, 0.75),
            labels = c("Before Inflection", "After Inflection"),
            expand = c(0, 0)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        # Theme
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        )
    
    # Create the purple mosaic plot
    purple_plot <- ggplot() +
        # Add rectangles
        geom_rect(data = purple_rects,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = location),
                  color = "white", size = 0.5) +
        # Add labels
        geom_text(data = purple_rects,
                  aes(x = label_x, y = label_y, label = label),
                  size = 4, fontface = "bold") +
        # Set colors
        scale_fill_manual(values = c("Below" = "#E8E1F2", "Above" = "#9B72CF")) +
        # Add titles
        labs(
            title = paste(pid, "- Purple Cells"),
            subtitle = paste("OR:", round(purple_or, 2))
        ) +
        # Set axis scales
        scale_x_continuous(
            breaks = c(0.25, 0.75),
            labels = c("Before Inflection", "After Inflection"),
            expand = c(0, 0)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        # Theme
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()
        )
    
    # Create a combined plot with gray on top and purple on bottom
    combined_plot <- cowplot::plot_grid(
        gray_plot, 
        purple_plot,
        ncol = 1,
        align = "v",
        rel_heights = c(1, 1)
    )
    
    return(list(
        gray_plot = gray_plot,
        purple_plot = purple_plot,
        combined_plot = combined_plot
    ))
}

#------------------------------------------------
# GENERATE AND SAVE PLOTS
#------------------------------------------------

# Generate all plots using the manual method
all_patient_plots <- lapply(1:7, create_manual_horizontal_mosaic)

# Display plots in R console
for (i in 1:7) {
    cat(paste("\nDisplaying plot for", patient_ids[i], ":\n"))
    print(all_patient_plots[[i]]$combined_plot)
}

# Save plots as png
for (i in 1:7) {
    pid <- patient_ids[i]
    filename <- paste0(pid, "_Mosaic_Plot.png")
    
    # Save png
    ggsave(
        filename = filename,
        plot = all_patient_plots[[i]]$combined_plot,
        device = "png",
        width = 10,
        height = 8,
        dpi = 300
    )
    
    cat(paste("Saved PNG file:", filename, "\n"))
}

# Create a combined visualization of all patients
create_all_patients_mosaic <- function() {
    # Placeholder for results
    all_plots <- list()
    
    # For each patient, get both gray and purple plots
    for (i in 1:7) {
        all_plots[[paste0("gray_", i)]] <- all_patient_plots[[i]]$gray_plot
        all_plots[[paste0("purple_", i)]] <- all_patient_plots[[i]]$purple_plot
    }
    
    # Create a gray cells mosaic with all patients
    combined_gray_mosaic <- cowplot::plot_grid(
        all_plots$gray_1, all_plots$gray_2, all_plots$gray_3,
        all_plots$gray_4, all_plots$gray_5, all_plots$gray_6,
        all_plots$gray_7,
        ncol = 1,
        align = "v"
    )
    
    # Create a purple cells mosaic with all patients
    combined_purple_mosaic <- cowplot::plot_grid(
        all_plots$purple_1, all_plots$purple_2, all_plots$purple_3,
        all_plots$purple_4, all_plots$purple_5, all_plots$purple_6,
        all_plots$purple_7,
        ncol = 1,
        align = "v"
    )
    
    # Add titles to the combined mosaics
    combined_gray_titled <- cowplot::ggdraw(combined_gray_mosaic) +
        cowplot::draw_label("Gray Cells Distribution - All Patients", 
                            x = 0.5, y = 0.99, size = 18, fontface = "bold")
    
    combined_purple_titled <- cowplot::ggdraw(combined_purple_mosaic) +
        cowplot::draw_label("Purple Cells Distribution - All Patients", 
                            x = 0.5, y = 0.99, size = 18, fontface = "bold")
    
    return(list(
        gray_mosaic = combined_gray_titled,
        purple_mosaic = combined_purple_titled
    ))
}


# =================
# Raincloud Plots
# =================
# Load required libraries
library(readxl)
library(overlapping)
library(ggplot2)
library(dplyr)
library(introdataviz)

# Read the data
data <- read_excel("C:.../Pt1_to_Pt7_Pre-Post-IP.xlsx") # Ribo expression values in pre-IP and post-IP for each patient

# Extract data for all patients and convert to numeric
patients_data <- list()
for(i in 1:7) {
    pre_col <- paste0("Pt", i, "_Pre-IP")
    post_col <- paste0("Pt", i, "_Post-IP")
    
    # Extract and clean data
    pre_data <- as.numeric(data[[pre_col]])
    post_data <- as.numeric(data[[post_col]])
    
    # Remove NA values
    pre_data <- pre_data[!is.na(pre_data)]
    post_data <- post_data[!is.na(post_data)]
    
    patients_data[[i]] <- list(
        pre = pre_data,
        post = post_data,
        patient_name = paste0("Patient ", i)
    )
}

# Function to calculate overlapping coefficient
calc_ovl <- function(x, y) {
    ovl_result <- overlap(list(x, y))
    return(ovl_result$OV)
}

# NEW FUNCTION: Calculate Extremity-Weighted OVL (EW-OVL) with improved tail weighting
calc_ew_ovl <- function(x, y) {
    # Calculate standard OVL
    standard_ovl <- calc_ovl(x, y)
    
    # Define extreme tails as values in one distribution that have no overlap with the other
    x_min <- min(x)
    x_max <- max(x)
    y_min <- min(y)
    y_max <- max(y)
    
    # Identify extreme tail cells as specified:
    # - pre-IP cells with values lower than the lowest value in post-IP
    # - post-IP cells with values higher than the highest value in pre-IP
    x_extreme_low <- x[x < y_min]  # Assuming x is pre-IP
    y_extreme_high <- y[y > x_max]  # Assuming y is post-IP
    
    # Calculate proportion of cells in extreme tails
    prop_x_extreme <- length(x_extreme_low) / length(x)
    prop_y_extreme <- length(y_extreme_high) / length(y)
    
    # Calculate mean distance of extreme values from the boundary
    if(length(x_extreme_low) > 0) {
        mean_x_dist <- mean(y_min - x_extreme_low)
    } else {
        mean_x_dist <- 0
    }
    
    if(length(y_extreme_high) > 0) {
        mean_y_dist <- mean(y_extreme_high - x_max)
    } else {
        mean_y_dist <- 0
    }
    
    # Normalize distances by range of all data
    total_range <- max(y_max, x_max) - min(y_min, x_min)
    norm_x_dist <- mean_x_dist / total_range
    norm_y_dist <- mean_y_dist / total_range
    
    # Calculate tail extremity factor that accounts for both:
    # 1. Proportion of cells in extreme tails
    # 2. How far these extreme cells are from the boundary
    
    # Scale by square of proportion to make it more sensitive to larger proportions
    x_factor <- norm_x_dist * (prop_x_extreme^2) * 4
    y_factor <- norm_y_dist * (prop_y_extreme^2) * 4
    
    # Combined extremity factor
    extremity_factor <- x_factor + y_factor
    
    # Apply a stronger penalty using square root to enhance sensitivity for smaller differences
    penalty_factor <- sqrt(extremity_factor) * 0.7
    
    # Cap penalty to avoid negative values
    penalty_factor <- min(penalty_factor, 0.9)
    
    # Calculate Extremity-Weighted OVL
    ew_ovl <- standard_ovl * (1 - penalty_factor)
    
    # Return diagnostic info
    attr(ew_ovl, "diagnostics") <- list(
        prop_x_extreme = prop_x_extreme,
        prop_y_extreme = prop_y_extreme,
        norm_x_dist = norm_x_dist, 
        norm_y_dist = norm_y_dist,
        penalty_factor = penalty_factor
    )
    
    return(ew_ovl)
}

# Calculate statistics for all patients
results <- list()
p_values <- c()

for(i in 1:7) {
    pre_data <- patients_data[[i]]$pre
    post_data <- patients_data[[i]]$post
    
    # Calculate OVL and EW-OVL
    ovl <- calc_ovl(pre_data, post_data)
    ew_ovl <- calc_ew_ovl(pre_data, post_data)
    
    # Calculate IQR
    iqr_pre <- quantile(pre_data, c(0.25, 0.75))
    iqr_post <- quantile(post_data, c(0.25, 0.75))
    
    # Perform Kolmogorov-Smirnov test
    ks_test <- ks.test(pre_data, post_data)
    p_values <- c(p_values, ks_test$p.value)
    
    # Store results
    results[[i]] <- list(
        patient = i,
        ovl = ovl,
        ew_ovl = ew_ovl,
        iqr_pre = iqr_pre,
        iqr_post = iqr_post,
        ks_test = ks_test,
        ew_diag = attr(ew_ovl, "diagnostics")
    )
    
    # Print summary for each patient
    cat("=== PATIENT", i, "RESULTS ===\n")
    cat("Data summary:\n")
    cat("Pre-IP: n =", length(pre_data), "range:", round(range(pre_data), 2), "\n")
    cat("Post-IP: n =", length(post_data), "range:", round(range(post_data), 2), "\n")
    cat("Overlapping Coefficient:", round(ovl, 4), "\n")
    cat("Extremity-Weighted OVL:", round(ew_ovl, 4), "\n")
    cat("KS Test p-value:", round(ks_test$p.value, 4), "\n")
    cat("KS Test D-statistic:", round(ks_test$statistic, 4), "\n\n")
}

# Apply Benjamini-Hochberg correction for multiple testing
p_adjusted <- p.adjust(p_values, method = "BH")

# Print adjusted p-values
cat("=== BENJAMINI-HOCHBERG CORRECTED P-VALUES ===\n")
for(i in 1:7) {
    cat("Patient", i, "BH-adjusted p-value:", round(p_adjusted[i], 4), "\n")
}

# Create data frame for plotting
plot_data <- data.frame()

for(i in 1:7) {
    pre_data <- patients_data[[i]]$pre
    post_data <- patients_data[[i]]$post
    patient_name <- patients_data[[i]]$patient_name
    
    patient_df <- data.frame(
        value = c(pre_data, post_data),
        patient = rep(patient_name, length(pre_data) + length(post_data)),
        condition = c(rep("Pre-IP", length(pre_data)), rep("Post-IP", length(post_data)))
    )
    
    plot_data <- rbind(plot_data, patient_df)
}

# Create facet labels with OVL, EW-OVL, and IQR information
facet_labels <- c()
for(i in 1:7) {
    res <- results[[i]]
    label <- paste0("Pt.", i, " (OVL = ", round(res$ovl * 100, 1), "%, EW-OVL = ", round(res$ew_ovl * 100, 1), 
                    "%;\nIQR pre = ", round(res$iqr_pre[1], 1), "-", round(res$iqr_pre[2], 1),
                    ", post = ", round(res$iqr_post[1], 1), "-", round(res$iqr_post[2], 1), ")")
    facet_labels <- c(facet_labels, label)
}

names(facet_labels) <- paste0("Patient ", 1:7)

# Set raincloud parameters
rain_height <- 0.08

# Create individual plots for each patient
individual_plots <- list()

for(i in 1:7) {
    # Filter data for current patient
    patient_data <- plot_data[plot_data$patient == paste0("Patient ", i), ]
    
    # Create individual plot
    individual_plot <- ggplot(patient_data, aes(x = "", y = value, fill = condition)) +
        # Density clouds
        geom_flat_violin(trim = FALSE, alpha = 0.6,
                         position = position_nudge(x = rain_height + 0.05)) +
        # Raw data points (rain)
        geom_point(aes(colour = condition), size = 1.2, alpha = 0.5, 
                   show.legend = FALSE,
                   position = position_jitter(width = rain_height, height = 0)) +
        # Boxplots
        geom_boxplot(width = rain_height, alpha = 0.6, 
                     show.legend = FALSE, outlier.shape = NA,
                     position = position_nudge(x = -rain_height * 2)) +
        # Mean and 95% CI
        stat_summary(fun.data = mean_cl_normal, 
                     mapping = aes(color = condition), 
                     show.legend = FALSE,
                     position = position_nudge(x = rain_height * 3)) +
        # Layout adjustments
        scale_x_discrete(name = "", expand = c(rain_height * 3, 0, 0, 0.7)) +
        scale_y_continuous(name = "Expression Value") +
        coord_flip() +
        # Custom colors
        scale_fill_manual(values = c("Pre-IP" = "#000066", "Post-IP" = "#21908CFF")) +
        scale_colour_manual(values = c("Pre-IP" = "#000066", "Post-IP" = "#21908CFF")) +
        labs(title = facet_labels[paste0("Patient ", i)],
             subtitle = paste("Pre-IP vs Post-IP Expression Distribution")) +
        theme_minimal() +
        theme(panel.grid.major.y = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold", size = 14),
              plot.subtitle = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    # Store the plot
    individual_plots[[i]] <- individual_plot
    
    # Display the plot
    print(individual_plot)
    
    # Save individual plot
    filename <- paste0("Raincloud_Patient_", i, ".png")
    ggsave(filename, plot = individual_plot, 
           width = 8, height = 6, dpi = 300)
    
    cat("Saved plot for Patient", i, "as", filename, "\n")
}

# Also create the combined plot for comparison
combined_plot <- ggplot(plot_data, aes(x = "", y = value, fill = condition)) +
    # Density clouds
    geom_flat_violin(trim = FALSE, alpha = 0.6,
                     position = position_nudge(x = rain_height + 0.05)) +
    # Raw data points (rain)
    geom_point(aes(colour = condition), size = 0.8, alpha = 0.4, 
               show.legend = FALSE,
               position = position_jitter(width = rain_height, height = 0)) +
    # Boxplots
    geom_boxplot(width = rain_height, alpha = 0.6, 
                 show.legend = FALSE, outlier.shape = NA,
                 position = position_nudge(x = -rain_height * 2)) +
    # Mean and 95% CI
    stat_summary(fun.data = mean_cl_normal, 
                 mapping = aes(color = condition), 
                 show.legend = FALSE,
                 position = position_nudge(x = rain_height * 3)) +
    # Layout adjustments
    scale_x_discrete(name = "", expand = c(rain_height * 3, 0, 0, 0.7)) +
    scale_y_continuous(name = "Expression Value") +
    coord_flip() +
    facet_wrap(~factor(patient, 
                       levels = paste0("Patient ", 1:7),
                       labels = facet_labels), 
               ncol = 2) +
    # Custom colors
    scale_fill_manual(values = c("Pre-IP" = "#000066", "Post-IP" = "#21908CFF")) +
    scale_colour_manual(values = c("Pre-IP" = "#000066", "Post-IP" = "#21908CFF")) +
    labs(title = "Raincloud Plot: Pre-IP vs Post-IP Expression Distributions (All Patients)",
         subtitle = "Comparing overlapping coefficients and extremity-weighted OVL across all 7 patients") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          strip.background = element_rect(fill = "gray95", color = "black", size = 0.25),
          strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold", size = 14),

          plot.subtitle = element_text(size = 12))

