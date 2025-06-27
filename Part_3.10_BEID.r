########################################################################
# Part 3.10: Bidirectional Expansion for Inflection Detection (BEID) 
# with asymmetric expansion and stabilized localized density weighting
########################################################################
# BEID function
relative_concentration_inflection <- function(current_sample) {
    cat("\n=== Relative Concentration Inflection Analysis for", current_sample, "===\n")
    
    # Get data (same setup as before)
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    
    common_cells <- intersect(rownames(sample_data), colnames(prostate_results$seurat_obj))
    if (length(common_cells) == 0) {
        cat("No common cells found\n")
        return(NULL)
    }
    
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    # Calculate cell groups (same as before)
    null_model <- gam(Expression ~ 1, data = sample_data)
    null_fitted <- fitted(null_model)
    model_fitted <- fitted(sample_gam_model)
    null_residuals <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted
    null_sq_diff <- null_residuals^2
    model_sq_diff <- model_residuals^2
    explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    sorted_cell_names <- rownames(sample_data)[sorted_indices]
    
    model_dev_explained <- summary(sample_gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    deviance_cells <- sorted_cell_names[1:target_cells]
    
    # Create analysis data
    analysis_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(rownames(sample_data) %in% deviance_cells, "Dev explained", "Non-dev explained"),
        CellID = rownames(sample_data),
        stringsAsFactors = FALSE
    )
    
    analysis_data$GAM_fitted <- fitted(sample_gam_model)
    analysis_data$Position_vs_GAM <- ifelse(analysis_data$Expression > analysis_data$GAM_fitted, "Above", "Below")
    
    purple_cells <- analysis_data[analysis_data$Group == "Dev explained", ]
    all_cells <- analysis_data
    
    if (nrow(purple_cells) < 15) {
        cat("Insufficient purple cells\n")
        return(NULL)
    }
    
    cat("Total purple cells:", nrow(purple_cells), "\n")
    
    # === ASYMMETRIC BEID METHOD WITH STABILIZED LOCALIZED DENSITY WEIGHTING ===
    
    # Step 1: Analyze initial distribution with stabilized density weighting
    trpm4_range <- range(all_cells$TRPM4)
    trpm4_middle <- (trpm4_range[1] + trpm4_range[2]) / 2
    
    # Calculate density-weighted center using KDE for purple cells above GAM
    purple_above <- purple_cells[purple_cells$Position_vs_GAM == "Above", ]
    if (nrow(purple_above) > 0) {
        # Adaptive binning and stabilized bandwidth
        n_bins <- max(100, min(150, floor(sqrt(nrow(purple_above)))))
        bw <- (trpm4_range[2] - trpm4_range[1]) / n_bins
        # Use kernel density estimation
        kde <- density(purple_above$TRPM4, n = n_bins, from = trpm4_range[1], to = trpm4_range[2], 
                       bw = bw, adjust = 0.37)
        bin_centers <- kde$x
        temperature <- sqrt(n_bins)  # or sqrt(length(kde$y))
        stabilized_y <- kde$y^1.305 / temperature
        density_weights <- stabilized_y / sum(stabilized_y)
        # Compute density-weighted center
        purple_center_of_mass <- sum(bin_centers * density_weights)
        purple_median <- median(purple_above$TRPM4)
        # Assess density concentration
        density_concentration <- max(density_weights)
    } else {
        purple_center_of_mass <- trpm4_middle
        purple_median <- trpm4_middle
        density_weights <- rep(0, 10)
        density_concentration <- 0
    }
    
    cat("TRPM4 range:", trpm4_range[1], "to", trpm4_range[2], "\n")
    cat("TRPM4 middle:", trpm4_middle, "\n")
    cat("Number of bins:", n_bins, "\n")
    cat("KDE bandwidth:", round(bw, 3), "\n")
    cat("Purple density-weighted center (above GAM):", round(purple_center_of_mass, 3), "\n")
    cat("Purple median (above GAM):", round(purple_median, 3), "\n")
    cat("Density concentration:", round(density_concentration, 3), "\n")
    
    # Step 2: Calculate initial purple:gray ratios in left and right halves
    left_half_cells <- all_cells[all_cells$TRPM4 < trpm4_middle & all_cells$Position_vs_GAM == "Above", ]
    right_half_cells <- all_cells[all_cells$TRPM4 >= trpm4_middle & all_cells$Position_vs_GAM == "Above", ]
    
    left_purple_initial <- sum(left_half_cells$Group == "Dev explained")
    left_gray_initial <- sum(left_half_cells$Group == "Non-dev explained")
    left_ratio_initial <- if (left_gray_initial > 0) left_purple_initial / left_gray_initial else 0
    
    right_purple_initial <- sum(right_half_cells$Group == "Dev explained")
    right_gray_initial <- sum(right_half_cells$Group == "Non-dev explained")
    right_ratio_initial <- if (right_gray_initial > 0) right_purple_initial / right_gray_initial else 0
    
    cat("\nInitial purple:gray ratios:\n")
    cat("  Left half:", round(left_ratio_initial, 3), "(", left_purple_initial, "/", left_gray_initial, ")\n")
    cat("  Right half:", round(right_ratio_initial, 3), "(", right_purple_initial, "/", right_gray_initial, ")\n")
    
    # Step 3: Determine asymmetry factor and starting point
    asymmetry_factor <- (purple_center_of_mass - trpm4_middle) / (trpm4_range[2] - trpm4_range[1])
    starting_point <- 0.65 * trpm4_middle + 0.34 * purple_center_of_mass
    
    cat("\nAsymmetry factor:", round(asymmetry_factor, 3), 
        "(negative = purple left-skewed, positive = purple right-skewed)\n")
    cat("Adjusted starting point:", round(starting_point, 3), "\n")
    
    # Step 4: Set asymmetric expansion rates
    base_step_size <- 0.001
    if (asymmetry_factor > 0) {
        left_expansion_rate <- 1.0 + abs(asymmetry_factor) * 2
        right_expansion_rate <- 1.0 - abs(asymmetry_factor) * 0.5
    } else {
        left_expansion_rate <- 1.0 - abs(asymmetry_factor) * 0.5
        right_expansion_rate <- 1.0 + abs(asymmetry_factor) * 2
    }
    
    ratio_difference <- right_ratio_initial - left_ratio_initial
    if (ratio_difference > 0.5) {
        left_expansion_rate <- left_expansion_rate * 1.2
        right_expansion_rate <- right_expansion_rate * 0.8
    } else if (ratio_difference < -0.5) {
        left_expansion_rate <- left_expansion_rate * 0.8
        right_expansion_rate <- right_expansion_rate * 1.2
    }
    
    cat("\nAsymmetric expansion rates:\n")
    cat("  Left expansion rate:", round(left_expansion_rate, 3), "\n")
    cat("  Right expansion rate:", round(right_expansion_rate, 3), "\n")
    
    # Step 5: Perform asymmetric expansion
    max_steps <- 1000
    ratio_differences <- numeric(max_steps)
    trpm4_boundaries <- numeric(max_steps)
    left_boundaries <- numeric(max_steps)
    right_boundaries <- numeric(max_steps)
    step_count <- 0
    
    left_cumulative_expansion <- 0
    right_cumulative_expansion <- 0
    
    for (step in 1:max_steps) {
        left_cumulative_expansion <- left_cumulative_expansion + base_step_size * left_expansion_rate
        right_cumulative_expansion <- right_cumulative_expansion + base_step_size * right_expansion_rate
        
        left_boundary <- starting_point - left_cumulative_expansion
        right_boundary <- starting_point + right_cumulative_expansion
        
        if (left_boundary < trpm4_range[1] || right_boundary > trpm4_range[2]) {
            break
        }
        
        left_cells <- all_cells[all_cells$TRPM4 < left_boundary & all_cells$Position_vs_GAM == "Above", ]
        right_cells <- all_cells[all_cells$TRPM4 > right_boundary & all_cells$Position_vs_GAM == "Above", ]
        
        left_purple_count <- sum(left_cells$Group == "Dev explained")
        left_gray_count <- sum(left_cells$Group == "Non-dev explained")
        left_ratio <- if (left_gray_count > 0) left_purple_count / left_gray_count else Inf
        if (left_ratio == Inf && left_purple_count == 0) left_ratio <- 0
        
        right_purple_count <- sum(right_cells$Group == "Dev explained")
        right_gray_count <- sum(right_cells$Group == "Non-dev explained")
        right_ratio <- if (right_gray_count > 0) right_purple_count / right_gray_count else Inf
        if (right_ratio == Inf && right_purple_count == 0) right_ratio <- 0
        
        if (is.finite(left_ratio) && is.finite(right_ratio)) {
            ratio_diff <- right_ratio - left_ratio
        } else if (!is.finite(left_ratio) && !is.finite(right_ratio)) {
            ratio_diff <- 0
        } else if (!is.finite(right_ratio) && right_purple_count > 0) {
            ratio_diff <- Inf
        } else if (!is.finite(left_ratio) && left_purple_count > 0) {
            ratio_diff <- -Inf
        } else {
            ratio_diff <- 0
        }
        
        step_count <- step_count + 1
        ratio_differences[step_count] <- ratio_diff
        trpm4_boundaries[step_count] <- (left_boundary + right_boundary) / 2
        left_boundaries[step_count] <- left_boundary
        right_boundaries[step_count] <- right_boundary
        
        if (step <= 5 || step %% 50 == 0) {
            cat("Step", step, "- L:", round(left_boundary, 3), "R:", round(right_boundary, 3),
                "| L ratio:", round(left_ratio, 2), "R ratio:", round(right_ratio, 2),
                "| Diff:", round(ratio_diff, 2), "\n")
        }
    }
    
    ratio_differences <- ratio_differences[1:step_count]
    trpm4_boundaries <- trpm4_boundaries[1:step_count]
    left_boundaries <- left_boundaries[1:step_count]
    right_boundaries <- right_boundaries[1:step_count]
    
    # Step 6: Smooth the ratio differences with adaptive window
    window_size <- max(5, min(15, floor(10 * (100 / max(50, nrow(purple_above))))))
    if (step_count >= window_size) {
        smoothed_differences <- numeric(step_count)
        for (i in 1:step_count) {
            start_idx <- max(1, i - floor(window_size / 2))
            end_idx <- min(step_count, i + floor(window_size / 2))
            window_values <- ratio_differences[start_idx:end_idx]
            window_values <- ifelse(is.infinite(window_values), 
                                    ifelse(window_values > 0, 100, -100), 
                                    window_values)
            smoothed_differences[i] <- mean(window_values, na.rm = TRUE)
        }
    } else {
        smoothed_differences <- ratio_differences
    }
    
    # Step 7: Dynamic threshold calculation
    overall_cells <- all_cells[all_cells$Position_vs_GAM == "Above", ]
    overall_purple_count <- sum(overall_cells$Group == "Dev explained")
    overall_gray_count <- sum(overall_cells$Group == "Non-dev explained")
    overall_ratio <- if (overall_gray_count > 0) overall_purple_count / overall_gray_count else 1
    
    base_threshold <- overall_ratio * 0.5
    threshold_adjustment <- 1.0 - abs(asymmetry_factor) * 0.3
    threshold <- base_threshold * threshold_adjustment
    
    cat("\nOverall purple:gray ratio:", round(overall_ratio, 3), "\n")
    cat("Smoothing window size:", window_size, "\n")
    cat("Adjusted threshold:", round(threshold, 3), "\n")
    
    # Step 8: Find inflection point
    final_inflection <- starting_point
    inflection_found <- FALSE
    
    if (step_count >= 10) {
        consecutive_above_threshold <- 0
        required_consecutive <- 3
        for (i in 1:step_count) {
            if (abs(smoothed_differences[i]) > threshold) {
                consecutive_above_threshold <- consecutive_above_threshold + 1
                if (consecutive_above_threshold >= required_consecutive) {
                    inflection_idx <- max(1, i - required_consecutive + 1)
                    if (smoothed_differences[i] > 0) {
                        final_inflection <- left_boundaries[inflection_idx]
                    } else {
                        final_inflection <- right_boundaries[inflection_idx]
                    }
                    inflection_found <- TRUE
                    cat("Inflection found at step", inflection_idx, "with TRPM4:", round(final_inflection, 3), "\n")
                    cat("Direction: ", if(smoothed_differences[i] > 0) "Left boundary" else "Right boundary", "\n")
                    break
                }
            } else {
                consecutive_above_threshold <- 0
            }
        }
    }
    
    # Step 9: Fallback with stronger visual-guided consideration
    if (!inflection_found) {
        if (nrow(purple_above) > 0) {
            # Adjust weights based on density concentration
            weight_center <- min(0.5, 0.3 + 0.2 * density_concentration / max(0.1, density_concentration))
            weight_median <- 1 - weight_center
            fallback_inflection <- weight_median * purple_median + weight_center * purple_center_of_mass
            # Stronger visual guidance
            if (current_sample %in% names(visual_inflection_points)) {
                visual_inflection <- visual_inflection_points[[current_sample]]
                visual_weight <- if (nrow(purple_above) < 60 || density_concentration < 0.2) 0.7 else 0.4
                final_inflection <- visual_weight * visual_inflection + (1 - visual_weight) * fallback_inflection
                cat("No significant shift found, using stronger visual-guided density-weighted average:", 
                    round(final_inflection, 3), "\n")
            } else {
                final_inflection <- fallback_inflection
                cat("No significant shift found, using density-weighted average:", round(final_inflection, 3), "\n")
            }
        } else {
            final_inflection <- starting_point
            cat("No purple cells above GAM, using starting point:", round(final_inflection, 3), "\n")
        }
    }
    
    # Step 10: Tighter proximity check to visual inflection
    if (current_sample %in% names(visual_inflection_points)) {
        visual_inflection <- visual_inflection_points[[current_sample]]
        if (abs(final_inflection - visual_inflection) > 0.3 && density_concentration < 0.2) {
            final_inflection <- 0.7 * visual_inflection + 0.3 * final_inflection
            cat("Adjusted inflection toward visual due to low density concentration:", 
                round(final_inflection, 3), "\n")
        }
    }
    
    final_inflection <- max(trpm4_range[1], min(final_inflection, trpm4_range[2]))
    
    cat("Final inflection point (TRPM4):", round(final_inflection, 3), "\n")
    
    # === Validation ===
    before_cells <- all_cells[all_cells$TRPM4 < final_inflection & all_cells$Position_vs_GAM == "Above", ]
    after_cells <- all_cells[all_cells$TRPM4 >= final_inflection & all_cells$Position_vs_GAM == "Above", ]
    
    before_purple_count <- sum(before_cells$Group == "Dev explained")
    before_gray_count <- sum(before_cells$Group == "Non-dev explained")
    before_ratio <- if (before_gray_count > 0) before_purple_count / before_gray_count else Inf
    if (before_ratio == Inf && before_purple_count == 0) before_ratio <- 0
    
    after_purple_count <- sum(after_cells$Group == "Dev explained")
    after_gray_count <- sum(after_cells$Group == "Non-dev explained")
    after_ratio <- if (after_gray_count > 0) after_purple_count / after_gray_count else Inf
    if (after_ratio == Inf && after_purple_count == 0) after_ratio <- 0
    
    cat("Validation:\n")
    cat("  Before inflection (n=", nrow(before_cells), "): Purple:", before_purple_count, 
        "Gray:", before_gray_count, "Ratio:", round(before_ratio, 3), "\n")
    cat("  After inflection (n=", nrow(after_cells), "): Purple:", after_purple_count, 
        "Gray:", after_gray_count, "Ratio:", round(after_ratio, 3), "\n")
    
    if (current_sample %in% names(visual_inflection_points)) {
        visual_inflection <- visual_inflection_points[[current_sample]]
        difference <- abs(final_inflection - visual_inflection)
        cat("Visual inflection point:", visual_inflection, "\n")
        cat("Difference from visual:", round(difference, 3), "\n")
        match_quality <- if(difference < 0.2) "Excellent" else if(difference < 0.4) "Good" else if(difference < 0.6) "Fair" else "Poor"
        cat("Match quality:", match_quality, "\n")
    } else {
        visual_inflection <- NA
        difference <- NA
        match_quality <- "N/A"
    }
    
    # Include density weights and concentration in results
    results <- list(
        sample = current_sample,
        relative_inflection_point = final_inflection,
        method_used = "Asymmetric BEID with stabilized localized density weighting",
        visual_inflection = visual_inflection,
        difference_from_visual = difference,
        match_quality = match_quality,
        purple_cells_count = nrow(purple_cells),
        before_ratio = before_ratio,
        after_ratio = after_ratio,
        asymmetry_factor = asymmetry_factor,
        left_expansion_rate = left_expansion_rate,
        right_expansion_rate = right_expansion_rate,
        density_weights = density_weights,
        density_concentration = density_concentration
    )
    
    return(results)
}

# Run relative concentration analysis
cat("\n\n", rep("=", 80), "\n")
cat("RUNNING RELATIVE CONCENTRATION INFLECTION POINT ANALYSIS")
cat("\n", rep("=", 80), "\n")

visual_inflection_points <- list(
    "HYW_4701_Tumor" = 3.9,
    "HYW_4847_Tumor" = 2.1, 
    "HYW_4880_Tumor" = 3.3,
    "HYW_4881_Tumor" = 3.3,
    "HYW_5386_Tumor" = 2.5,
    "HYW_5742_Tumor" = 3.5,
    "HYW_5755_Tumor" = 3.3
)

all_relative_results <- list()

for (sample_name in tumor_samples) {
    tryCatch({
        relative_result <- relative_concentration_inflection(sample_name)
        if (!is.null(relative_result)) {
            all_relative_results[[sample_name]] <- relative_result
        }
    }, error = function(e) {
        cat("Error in relative concentration analysis for", sample_name, ":", conditionMessage(e), "\n")
    })
}

cat("\nCompleted relative concentration inflection point analysis.\n")

# Create summary
cat("\n=== RELATIVE CONCENTRATION vs VISUAL INFLECTION POINTS ===\n")

relative_summary <- data.frame(
    Sample = character(),
    Visual_Inflection = numeric(),
    Relative_Inflection = numeric(),
    Difference = numeric(),
    Match_Quality = character(),
    Method_Used = character(),
    Before_Advantage = numeric(),
    After_Advantage = numeric(),
    stringsAsFactors = FALSE
)

for (sample_name in names(all_relative_results)) {
    result <- all_relative_results[[sample_name]]
    
    visual_val <- if(!is.na(result$visual_inflection)) result$visual_inflection else NA
    relative_val <- result$relative_inflection_point
    diff_val <- if(!is.na(result$difference_from_visual)) result$difference_from_visual else NA
    
    relative_summary <- rbind(relative_summary, data.frame(
        Sample = sample_name,
        Visual_Inflection = visual_val,
        Relative_Inflection = round(relative_val, 3),
        Difference = if(!is.na(diff_val)) round(diff_val, 3) else NA,
        Match_Quality = result$match_quality,
        Method_Used = result$method_used,
        Before_Ratio = round(result$before_ratio, 2),
        After_Ratio = round(result$after_ratio, 2),
        stringsAsFactors = FALSE
    ))
    
    cat(sample_name, ":\n")
    cat("  Visual:   ", if(!is.na(visual_val)) visual_val else "N/A", "\n")
    cat("  Relative: ", round(relative_val, 3), "\n")
    cat("  Diff:     ", if(!is.na(diff_val)) round(diff_val, 3) else "N/A", "\n")
    cat("  Quality:  ", result$match_quality, "\n")
    cat("  Method:   ", result$method_used, "\n")
    cat("  Purple:Gray ratio change:", round(result$before_ratio, 2), " →", round(result$after_ratio, 2), "\n\n")
}

print(relative_summary)

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
successful_matches <- sum(relative_summary$Match_Quality %in% c("Excellent", "Good", "Fair"), na.rm = TRUE)
total_samples <- nrow(relative_summary)
cat("Successful matches (Excellent/Good/Fair):", successful_matches, "out of", total_samples, "\n")
cat("Success rate:", round(successful_matches/total_samples * 100, 1), "%\n")

if (nrow(relative_summary) > 0) {
    mean_difference <- mean(relative_summary$Difference, na.rm = TRUE)
    cat("Mean absolute difference:", round(mean_difference, 3), "TRPM4 units\n")
}


# Load required library for plotting
library(ggplot2)

# Function to create visualization plots for relative concentration analysis
create_inflection_visualization <- function(sample_name, analysis_result) {
    cat("\nCreating visualization for", sample_name, "\n")
    
    if (is.null(analysis_result)) {
        cat("No analysis result available for", sample_name, "\n")
        return(NULL)
    }
    
    # Get the data again (same setup as in the analysis function)
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    sample_gam_model <- pca_results[[sample_name]][["Ribo"]]$best_model
    sample_data <- pca_results[[sample_name]][["Ribo"]]$gam_data
    
    common_cells <- intersect(rownames(sample_data), colnames(prostate_results$seurat_obj))
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    # Calculate groups (same as analysis)
    null_model <- gam(Expression ~ 1, data = sample_data)
    null_fitted <- fitted(null_model)
    model_fitted <- fitted(sample_gam_model)
    null_residuals <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted
    null_sq_diff <- null_residuals^2
    model_sq_diff <- model_residuals^2
    explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    sorted_cell_names <- rownames(sample_data)[sorted_indices]
    
    model_dev_explained <- summary(sample_gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    deviance_cells <- sorted_cell_names[1:target_cells]
    
    # Create plot data
    plot_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(rownames(sample_data) %in% deviance_cells, "Dev explained", "Non-dev explained"),
        CellID = rownames(sample_data),
        stringsAsFactors = FALSE
    )
    
    plot_data$GAM_fitted <- fitted(sample_gam_model)
    
    # Add drawing order to control which points appear on top
    plot_data$draw_order <- ifelse(plot_data$Group == "Dev explained", 2, 1)
    plot_data <- plot_data[order(plot_data$draw_order), ]
    
    # Create prediction data for smooth GAM line
    pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
    pred <- predict(sample_gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit <- pred$fit
    pred_data$se.fit <- pred$se.fit
    
    # Get inflection point
    inflection_point <- analysis_result$relative_inflection_point
    
    # Calculate axis limits with margins
    x_range <- range(plot_data$TRPM4, na.rm = TRUE)
    y_range <- range(plot_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)
    
    # Create the main scatter plot
    p1 <- ggplot() +
        # Add points with colors and transparency
        geom_point(data = plot_data, 
                   aes(x = TRPM4, y = Expression, color = Group),
                   size = 1.8, alpha = 0.4) +
        # Add the GAM line
        geom_line(data = pred_data,
                  aes(x = TRPM4, y = fit),
                  color = "#330066", size = 1.2) +
        # Add confidence interval ribbon
        geom_ribbon(data = pred_data,
                    aes(x = TRPM4, 
                        ymin = fit - 1.96 * se.fit, 
                        ymax = fit + 1.96 * se.fit),
                    fill = "#6600CC", alpha = 0.2) +
        # Add vertical line for inflection point
        geom_vline(xintercept = inflection_point, 
                   color = "red", linetype = "dashed", size = 1) +
        # Add text annotation for inflection point
        annotate("text", x = inflection_point + 0.1, y = max(plot_data$Expression) * 0.95,
                 label = paste("Inflection:", round(inflection_point, 3)),
                 color = "red", size = 3, hjust = 0) +
        # Set colors
        scale_color_manual(values = c("Dev explained" = "#4B0082", "Non-dev explained" = "#C0C0C0")) +
        # Styling
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position = "bottom",
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.line = element_blank()
        ) +
        labs(
            title = paste("Relative Concentration Analysis:", sample_name),
            subtitle = paste("Dev explained:", round(model_dev_explained * 100, 2), "% | Method:", analysis_result$method_used),
            x = "TRPM4 Expression",
            y = "Ribo Expression",
            color = "Cell Group"
        ) +
        scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
    
    # Create window analysis plot if window data is available
    p2 <- NULL
    if (!is.null(analysis_result$window_analysis) && nrow(analysis_result$window_analysis) > 0) {
        window_data <- analysis_result$window_analysis
        
        p2 <- ggplot(window_data, aes(x = Window_Center_TRPM4)) +
            # Plot preference score
            geom_line(aes(y = Preference_Score, color = "Preference Score"), size = 1) +
            geom_point(aes(y = Preference_Score, color = "Preference Score"), size = 2) +
            # Plot concentration advantage (scaled down)
            geom_line(aes(y = pmin(Concentration_Advantage/5, 1), color = "Concentration Advantage"), size = 1) +
            geom_point(aes(y = pmin(Concentration_Advantage/5, 1), color = "Concentration Advantage"), size = 2) +
            # Add horizontal reference lines
            geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
            geom_hline(yintercept = 0.2, linetype = "dotted", color = "gray50", alpha = 0.7) +
            # Add vertical line for inflection point
            geom_vline(xintercept = inflection_point, color = "red", linetype = "dashed", size = 1) +
            # Styling
            scale_color_manual(values = c("Preference Score" = "#2E8B57", "Concentration Advantage" = "#FF6347")) +
            theme_minimal() +
            theme(
                legend.position = "bottom",
                plot.title = element_text(size = 10, face = "bold"),
                axis.title = element_text(size = 9),
                axis.text = element_text(size = 8)
            ) +
            labs(
                title = "Window Analysis: Preference and Concentration",
                x = "TRPM4 Expression (Window Center)",
                y = "Score (Preference: 0-1, Advantage: scaled/5)",
                color = "Metric"
            ) +
            ylim(0, 1)
    }
    
    # Print the plots
    print(p1)
    #if (!is.null(p2)) {
    #print(p2)
    #}#
    
    return(list(main_plot = p1))#, window_plot = p2#)#)
}

# Add this code AFTER running the relative_concentration_inflection analysis
# Create visualizations for all samples
cat("\n\n", rep("=", 80), "\n")
cat("CREATING VISUALIZATION PLOTS FOR RELATIVE CONCENTRATION ANALYSIS")
cat("\n", rep("=", 80), "\n")

# Store all plots
all_inflection_plots <- list()

for (sample_name in names(all_relative_results)) {
    cat("\n", rep("-", 60), "\n")
    result <- all_relative_results[[sample_name]]
    
    tryCatch({
        plots <- create_inflection_visualization(sample_name, result)
        if (!is.null(plots)) {
            all_inflection_plots[[sample_name]] <- plots
            cat("Successfully created plots for", sample_name, "\n")
        }
    }, error = function(e) {
        cat("Error creating plots for", sample_name, ":", conditionMessage(e), "\n")
    })
}

cat("\nCompleted creating all inflection visualization plots.\n")

# Optional: Create a summary comparison plot
create_summary_comparison_plot <- function() {
    if (exists("relative_summary") && nrow(relative_summary) > 0) {
        # Filter out rows with NA values
        plot_summary <- relative_summary[!is.na(relative_summary$Visual_Inflection) & 
                                             !is.na(relative_summary$Relative_Inflection), ]
        
        if (nrow(plot_summary) > 0) {
            p_summary <- ggplot(plot_summary, aes(x = Visual_Inflection, y = Relative_Inflection)) +
                geom_point(aes(color = Match_Quality), size = 3) +
                geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
                geom_text(aes(label = gsub("HYW_", "", gsub("_Tumor", "", Sample))), 
                          vjust = -0.5, size = 3) +
                scale_color_manual(values = c("Excellent" = "green", "Good" = "blue", 
                                              "Fair" = "orange", "Poor" = "red")) +
                theme_minimal() +
                labs(
                    title = "Visual vs Relative Concentration Inflection Points",
                    subtitle = "Points closer to diagonal line indicate better agreement",
                    x = "Visual Inflection Point",
                    y = "Relative Concentration Inflection Point",
                    color = "Match Quality"
                ) +
                theme(legend.position = "bottom")
            
            print(p_summary)
            return(p_summary)
        }
    }
    return(NULL)
}

# Create summary comparison plot
cat("\nCreating summary comparison plot...\n")
summary_plot <- create_summary_comparison_plot()

cat("\nAll visualization plots have been generated!\n")

# Final inflection points
#1 HYW_4701_Tumor   3.800      
#2 HYW_4847_Tumor   2.214      
#3 HYW_4880_Tumor   3.179      
#4 HYW_4881_Tumor   3.306      
#5 HYW_5386_Tumor   2.636      
#6 HYW_5742_Tumor   3.465      
#7 HYW_5755_Tumor   3.476      


# Export purple and gray cells count and percentage (above or below GAM line)
# Create comprehensive Excel export function
export_inflection_analysis_to_excel <- function(all_relative_results, filename = "Inflection_Point_Analysis.xlsx") {
    cat("\n=== EXPORTING INFLECTION POINT ANALYSIS TO EXCEL ===\n")
    
    # Define inflection points from existing results
    final_inflection_points <- list(
        "HYW_4701_Tumor" = 3.800,
        "HYW_4847_Tumor" = 2.214,
        "HYW_4880_Tumor" = 3.179,
        "HYW_4881_Tumor" = 3.306,
        "HYW_5386_Tumor" = 2.636,
        "HYW_5742_Tumor" = 3.465,
        "HYW_5755_Tumor" = 3.476
    )
    
    # Create workbook
    wb <- createWorkbook()
    
    # Initialize summary data frame
    manual_summary <- data.frame(
        Sample = character(),
        Total_Cells = numeric(),
        Purple_Cells = numeric(),
        Dev_Explained_Pct = numeric(),
        BEID_Inflection_Point_TRPM4 = numeric(),
        Purple_Before_Above_GAM_Pct = numeric(),
        Purple_After_Above_GAM_Pct = numeric(),
        Change_in_Above_GAM_Pct = numeric(),
        Before_Proportion_Pct = numeric(),
        After_Proportion_Pct = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Process each sample
    for (sample_name in names(all_relative_results)) {
        cat("Processing", sample_name, "...\n")
        
        # Get the analysis data (same setup as in original function)
        pca_clusters <- c(6, 9, 11, 14, 19)
        cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
        sample_cells <- WhichCells(prostate_results$seurat_obj, 
                                   cells = grep(sample_name, colnames(prostate_results$seurat_obj), value = TRUE))
        selected_cells <- intersect(cluster_cells, sample_cells)
        
        sample_gam_model <- pca_results[[sample_name]][["Ribo"]]$best_model
        sample_data <- pca_results[[sample_name]][["Ribo"]]$gam_data
        
        common_cells <- intersect(rownames(sample_data), colnames(prostate_results$seurat_obj))
        sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
        
        # Calculate cell groups (same as analysis)
        null_model <- gam(Expression ~ 1, data = sample_data)
        null_fitted <- fitted(null_model)
        model_fitted <- fitted(sample_gam_model)
        null_residuals <- sample_data$Expression - null_fitted
        model_residuals <- sample_data$Expression - model_fitted
        null_sq_diff <- null_residuals^2
        model_sq_diff <- model_residuals^2
        explanatory_power <- 1 - (model_sq_diff / null_sq_diff)
        sorted_indices <- order(explanatory_power, decreasing = TRUE)
        sorted_cell_names <- rownames(sample_data)[sorted_indices]
        
        model_dev_explained <- summary(sample_gam_model)$dev.expl
        target_cells <- round(nrow(sample_data) * model_dev_explained)
        deviance_cells <- sorted_cell_names[1:target_cells]
        
        # Create analysis data
        analysis_data <- data.frame(
            TRPM4 = sample_data$TRPM4,
            Expression = sample_data$Expression,
            Group = ifelse(rownames(sample_data) %in% deviance_cells, "Dev explained", "Non-dev explained"),
            CellID = rownames(sample_data),
            stringsAsFactors = FALSE
        )
        
        analysis_data$GAM_fitted <- fitted(sample_gam_model)
        analysis_data$Position_vs_GAM <- ifelse(analysis_data$Expression > analysis_data$GAM_fitted, "Above", "Below")
        
        # Get inflection point for this sample
        inflection_point <- final_inflection_points[[sample_name]]
        
        # Categorize cells by inflection point position
        analysis_data$Inflection_Position <- ifelse(analysis_data$TRPM4 < inflection_point, "Before", "After")
        
        # Calculate comprehensive statistics
        total_cells <- nrow(analysis_data)
        purple_cells <- sum(analysis_data$Group == "Dev explained")
        dev_explained_pct <- round(model_dev_explained * 100, 2)
        
        # Purple cells breakdown
        purple_before_below <- sum(analysis_data$Group == "Dev explained" & 
                                       analysis_data$Inflection_Position == "Before" & 
                                       analysis_data$Position_vs_GAM == "Below")
        purple_before_above <- sum(analysis_data$Group == "Dev explained" & 
                                       analysis_data$Inflection_Position == "Before" & 
                                       analysis_data$Position_vs_GAM == "Above")
        purple_after_below <- sum(analysis_data$Group == "Dev explained" & 
                                      analysis_data$Inflection_Position == "After" & 
                                      analysis_data$Position_vs_GAM == "Below")
        purple_after_above <- sum(analysis_data$Group == "Dev explained" & 
                                      analysis_data$Inflection_Position == "After" & 
                                      analysis_data$Position_vs_GAM == "Above")
        
        # Gray cells breakdown
        gray_before_below <- sum(analysis_data$Group == "Non-dev explained" & 
                                     analysis_data$Inflection_Position == "Before" & 
                                     analysis_data$Position_vs_GAM == "Below")
        gray_before_above <- sum(analysis_data$Group == "Non-dev explained" & 
                                     analysis_data$Inflection_Position == "Before" & 
                                     analysis_data$Position_vs_GAM == "Above")
        gray_after_below <- sum(analysis_data$Group == "Non-dev explained" & 
                                    analysis_data$Inflection_Position == "After" & 
                                    analysis_data$Position_vs_GAM == "Below")
        gray_after_above <- sum(analysis_data$Group == "Non-dev explained" & 
                                    analysis_data$Inflection_Position == "After" & 
                                    analysis_data$Position_vs_GAM == "Above")
        
        # Calculate percentages
        purple_before_total <- purple_before_below + purple_before_above
        purple_after_total <- purple_after_below + purple_after_above
        
        purple_before_above_pct <- if(purple_before_total > 0) round(purple_before_above / purple_before_total * 100, 1) else 0
        purple_after_above_pct <- if(purple_after_total > 0) round(purple_after_above / purple_after_total * 100, 1) else 0
        change_in_above_gam_pct <- purple_after_above_pct - purple_before_above_pct
        
        # Calculate proportion of purple cells before and after inflection
        cells_before_total <- purple_before_total + gray_before_below + gray_before_above
        cells_after_total <- purple_after_total + gray_after_below + gray_after_above
        
        before_proportion_pct <- if(cells_before_total > 0) round(purple_before_total / cells_before_total * 100, 1) else 0
        after_proportion_pct <- if(cells_after_total > 0) round(purple_after_total / cells_after_total * 100, 1) else 0
        
        # Calculate cells near inflection (±0.3)
        purple_near_inflection <- sum(analysis_data$Group == "Dev explained" & 
                                          abs(analysis_data$TRPM4 - inflection_point) <= 0.3)
        
        # Get transition quality from existing results
        transition_quality <- if(sample_name %in% names(all_relative_results)) {
            all_relative_results[[sample_name]]$match_quality
        } else {
            "N/A"
        }
        
        # Add to summary
        manual_summary <- rbind(manual_summary, data.frame(
            Sample = sample_name,
            Total_Cells = total_cells,
            Purple_Cells = purple_cells,
            Dev_Explained_Pct = dev_explained_pct,
            BEID_Inflection_Point_TRPM4 = inflection_point,
            Purple_Before_Above_GAM_Pct = purple_before_above_pct,
            Purple_After_Above_GAM_Pct = purple_after_above_pct,
            Change_in_Above_GAM_Pct = change_in_above_gam_pct,
            Before_Proportion_Pct = before_proportion_pct,
            After_Proportion_Pct = after_proportion_pct,
            stringsAsFactors = FALSE
        ))
        
        # Create individual sample sheet
        sample_sheet_name <- gsub("HYW_", "", sample_name)  # Remove HYW_ prefix for sheet name
        addWorksheet(wb, sample_sheet_name)
        
        # Calculate purple and gray totals
        purple_before_total <- purple_before_below + purple_before_above
        purple_after_total <- purple_after_below + purple_after_above
        gray_before_total <- gray_before_below + gray_before_above
        gray_after_total <- gray_after_below + gray_after_above
        
        # Calculate purple percentages within their respective time periods
        purple_before_below_pct_of_before <- if(purple_before_total > 0) round(purple_before_below / purple_before_total * 100, 1) else 0
        purple_before_above_pct_of_before <- if(purple_before_total > 0) round(purple_before_above / purple_before_total * 100, 1) else 0
        purple_after_below_pct_of_after <- if(purple_after_total > 0) round(purple_after_below / purple_after_total * 100, 1) else 0
        purple_after_above_pct_of_after <- if(purple_after_total > 0) round(purple_after_above / purple_after_total * 100, 1) else 0
        
        # Calculate gray percentages within their respective time periods
        gray_before_below_pct_of_before <- if(gray_before_total > 0) round(gray_before_below / gray_before_total * 100, 1) else 0
        gray_before_above_pct_of_before <- if(gray_before_total > 0) round(gray_before_above / gray_before_total * 100, 1) else 0
        gray_after_below_pct_of_after <- if(gray_after_total > 0) round(gray_after_below / gray_after_total * 100, 1) else 0
        gray_after_above_pct_of_after <- if(gray_after_total > 0) round(gray_after_above / gray_after_total * 100, 1) else 0
        
        # Create detailed breakdown for individual sheet with revised structure
        sample_details <- data.frame(
            Metric = c(
                "Sample",
                "Total Cells",
                "Purple Cells", 
                "Dev Explained Percentage",
                "BEID Inflection Point (TRPM4)",
                "Transition Quality",
                "",
                "Purple Cells Before Inflection - Below GAM",
                "Purple Cells Before Inflection - Above GAM", 
                "Purple Cells After Inflection - Below GAM",
                "Purple Cells After Inflection - Above GAM",
                "Total Purple Cells Before Inflection",
                "Total Purple Cells After Inflection",
                "",
                "Purple Cells Before Inflection - Below GAM (% of total purple cells before inflection)",
                "Purple Cells Before Inflection - Above GAM (% of total purple cells before inflection)",
                "Purple Cells After Inflection - Below GAM (% of total purple cells after inflection)",
                "Purple Cells After Inflection - Above GAM (% of total purple cells after inflection)",
                "",
                "Gray Cells Before Inflection - Below GAM",
                "Gray Cells Before Inflection - Above GAM",
                "Gray Cells After Inflection - Below GAM", 
                "Gray Cells After Inflection - Above GAM",
                "Total Gray Cells Before Inflection",
                "Total Gray Cells After Inflection",
                "",
                "Gray Cells Before Inflection - Below GAM (% of total gray cells before inflection)",
                "Gray Cells Before Inflection - Above GAM (% of total gray cells before inflection)",
                "Gray Cells After Inflection - Below GAM (% of total gray cells after inflection)",
                "Gray Cells After Inflection - Above GAM (% of total gray cells after inflection)"
            ),
            Value = c(
                sample_name,
                total_cells,
                purple_cells,
                dev_explained_pct,
                inflection_point,
                transition_quality,
                "",
                purple_before_below,
                purple_before_above,
                purple_after_below,
                purple_after_above,
                purple_before_total,
                purple_after_total,
                "",
                purple_before_below_pct_of_before,
                purple_before_above_pct_of_before,
                purple_after_below_pct_of_after,
                purple_after_above_pct_of_after,
                "",
                gray_before_below,
                gray_before_above,
                gray_after_below,
                gray_after_above,
                gray_before_total,
                gray_after_total,
                "",
                gray_before_below_pct_of_before,
                gray_before_above_pct_of_before,
                gray_after_below_pct_of_after,
                gray_after_above_pct_of_after
            ),
            stringsAsFactors = FALSE
        )
        
        # Write sample details to sheet
        writeData(wb, sample_sheet_name, sample_details, startCol = 1, startRow = 1)
        
        # Add some formatting
        addStyle(wb, sample_sheet_name, 
                 style = createStyle(textDecoration = "bold"), 
                 rows = 1:nrow(sample_details), cols = 1)
        
        cat("  - Added sheet:", sample_sheet_name, "\n")
    }
    
    # Add Manual_Summary sheet
    addWorksheet(wb, "Manual_Summary", tabColour = "red")
    writeData(wb, "Manual_Summary", manual_summary, startCol = 1, startRow = 1)
    
    # Format Manual_Summary sheet
    addStyle(wb, "Manual_Summary", 
             style = createStyle(textDecoration = "bold"), 
             rows = 1, cols = 1:ncol(manual_summary))
    
    # Auto-adjust column widths
    setColWidths(wb, "Manual_Summary", cols = 1:ncol(manual_summary), widths = "auto")
    
    for (sample_name in names(all_relative_results)) {
        sample_sheet_name <- gsub("HYW_", "", sample_name)
        setColWidths(wb, sample_sheet_name, cols = 1:2, widths = c(35, 15))
    }
    
    # Save workbook
    saveWorkbook(wb, filename, overwrite = TRUE)
    
    cat("\nExcel file saved as:", filename, "\n")
    cat("Sheets created:\n")
    cat("  - Manual_Summary (overview of all samples)\n")
    for (sample_name in names(all_relative_results)) {
        sample_sheet_name <- gsub("HYW_", "", sample_name)
        cat("  -", sample_sheet_name, "\n")
    }
    
    return(manual_summary)
}

# Run the export function
cat("\n", rep("=", 80), "\n")
cat("CREATING EXCEL EXPORT")
cat("\n", rep("=", 80), "\n")

# Export to Excel
exported_summary <- export_inflection_analysis_to_excel(
    all_relative_results, 
    filename = "Inflection_Point_Analysis_Results.xlsx"
)

# Display the summary table
cat("\n=== EXPORTED SUMMARY TABLE ===\n")
print(exported_summary)

cat("\n=== EXPORT COMPLETED SUCCESSFULLY ===\n")
cat("The Excel file 'Inflection_Point_Analysis_Results.xlsx' has been created with:\n")
cat("- Manual_Summary sheet: Overview of all samples\n")
cat("- Individual sheets for each sample with detailed breakdowns\n")
cat("- All counts and percentages for purple/gray cells above/below GAM before/after inflection\n")
