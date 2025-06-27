#############################################################
# Part 3.6: TPRS Visualization and GAM Components Analysis
##############################################################

# =========================================
# 1. TPRS Visualization
# =========================================
# Final visualization function with individual basis function plots (modified)
visualize_tprs_improved <- function(data, k = 10, lambda = 0.52801317696145, 
                                    gamma = 1.5, sample_id = "HYW_4881_Tumor", 
                                    gene_set = "Ribo", predefined_coefficients = NULL) {
    
    # Ensures reproducibility of random processes.
    set.seed(123)
    
    # Generates a sequence of TRPM4 values for smooth plotting.
    x_seq <- seq(min(data$TRPM4), max(data$TRPM4), length.out = 500)
    
    # Fits a temporary GAM model using REML with specified gamma for basis function extraction.
    temp_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                      data = data, method = "REML", sp = lambda, gamma = gamma)
    
    # Assigns predefined coefficients for consistent model interpretation.
    intercept <- predefined_coefficients$intercept
    linear_term <- predefined_coefficients$linear_term
    
    # Calculates the number of meaningful basis functions based on k.
    n_basis <- k - 2
    
    # Selects the first n_basis smooth coefficients, excluding near-zero values.
    smooth_coefs <- predefined_coefficients$smooth_coefs[1:n_basis]
    
    # Creates a prediction matrix for basis function evaluation at x_seq.
    X_pred <- predict(temp_model, newdata = data.frame(TRPM4 = x_seq), type = "lpmatrix")
    
    # Identifies columns corresponding to smooth terms in the prediction matrix.
    smooth_cols <- grep("s\\(TRPM4", colnames(X_pred))[1:n_basis]
    
    # Initializes a data frame to store basis function values.
    basis_data <- data.frame(TRPM4 = x_seq)
    
    # Extracts individual basis functions into the data frame.
    for (i in 1:n_basis) {
        basis_data[[paste0("phi_", i)]] <- X_pred[, smooth_cols[i]]
    }
    
    # Computes weighted basis functions by applying coefficients.
    for (i in 1:n_basis) {
        basis_data[[paste0("weighted_phi_", i)]] <- smooth_coefs[i] * basis_data[[paste0("phi_", i)]]
    }
    
    # Calculates the linear component using the predefined linear term.
    basis_data$linear_only <- linear_term * basis_data$TRPM4
    
    # Combines intercept and linear term for plotting.
    basis_data$intercept_plus_linear <- intercept + basis_data$linear_only
    
    # Sums weighted basis functions to obtain the smooth component without intercept.
    weighted_cols <- grep("^weighted_phi_", names(basis_data))
    if (length(weighted_cols) > 0) {
        basis_data$smooth_only <- rowSums(basis_data[, weighted_cols, drop = FALSE])
    } else {
        basis_data$smooth_only <- 0
    }
    
    # Computes the smooth component with intercept included.
    basis_data$intercept_plus_smooth <- intercept + basis_data$smooth_only
    
    # Calculates the full model fit by combining all components.
    basis_data$full_fit <- intercept + basis_data$linear_only + basis_data$smooth_only
    
    # Retrieves knot positions from the fitted model.
    knots <- temp_model$smooth[[1]]$xp
    if (is.null(knots) || length(knots) == 0) {
        n_knots <- k - 2
        knots <- quantile(data$TRPM4, probs = seq(0, 1, length.out = n_knots))
    }
    
    # Defines a custom color palette using the mako scheme from viridis.
    mako_palette <- viridis::viridis(100, option = "mako")[5:95]
    basis_colors <- mako_palette[seq(1, length(mako_palette), length.out = n_basis)]
    
    # Maps colors directly to basis functions for consistent visualization.
    direct_color_mapping <- setNames(basis_colors, paste0("φ", 1:n_basis))
    
    # Specifies a darker gray color for knot lines in plots.
    knot_color <- "#333333"
    
    # Plots all TPRS basis functions with knot locations.
    p1 <- ggplot() +
        geom_vline(xintercept = knots, linetype = "dashed", color = knot_color, alpha = 0.8) +
        geom_line(data = pivot_longer(basis_data, cols = starts_with("phi_"), 
                                      names_to = "basis_function", values_to = "value"),
                  aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        scale_color_manual(values = setNames(basis_colors, paste0("phi_", 1:n_basis)),
                           labels = paste0("φ", 1:n_basis)) +
        labs(title = paste("All TPRS Basis Functions for", sample_id, "-", gene_set),
             subtitle = paste("k =", k, ", λ =", lambda, ", γ =", gamma, "with", length(knots), "knots"),
             x = "TRPM4 Expression (log2)", 
             y = "Basis Function Value") +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold"))
    
    # Visualizes weighted basis functions with knot lines.
    p2 <- ggplot() +
        geom_vline(xintercept = knots, linetype = "dashed", color = knot_color, alpha = 0.8) +
        geom_line(data = pivot_longer(basis_data, cols = starts_with("weighted_phi_"), 
                                      names_to = "basis_function", values_to = "value"),
                  aes(x = TRPM4, y = value, color = gsub("weighted_phi_", "phi_", basis_function), 
                      group = basis_function), 
                  size = 1) +
        scale_color_manual(values = setNames(basis_colors, paste0("phi_", 1:n_basis)),
                           labels = paste0("φ", 1:n_basis)) +
        labs(title = paste("Weighted Basis Functions for", sample_id, "-", gene_set),
             subtitle = "Each basis function multiplied by its coefficient",
             x = "TRPM4 Expression (log2)", 
             y = "Weighted Value") +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold"))
    
    # Prepares unscaled basis function data for individual plotting.
    basis_long <- pivot_longer(basis_data[, c("TRPM4", paste0("phi_", 1:n_basis))], 
                               cols = -TRPM4, 
                               names_to = "basis_function", 
                               values_to = "value")
    
    # Extracts basis function numbers for faceting.
    basis_long$basis_num <- as.numeric(gsub("phi_", "", basis_long$basis_function))
    
    # Creates a grid of individual unscaled basis function plots.
    p_basis_grid_unscaled <- ggplot(basis_long, aes(x = TRPM4, y = value)) +
        geom_line(aes(color = basis_function), size = 1) +
        scale_color_manual(values = setNames(basis_colors[1:n_basis], paste0("phi_", 1:n_basis)),
                           labels = paste0("φ", 1:n_basis)) +
        facet_wrap(~ basis_num, ncol = 4) +
        labs(title = "Individual TPRS Basis Functions (Unscaled)",
             x = "TRPM4 Expression (log2)",
             y = "Value") +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"),
              strip.text = element_text(size = 12, face = "bold", color = "navy"))
    
    # Prepares unscaled weighted basis function data for individual plotting.
    weighted_basis_long <- pivot_longer(basis_data[, c("TRPM4", paste0("weighted_phi_", 1:n_basis))], 
                                        cols = -TRPM4, 
                                        names_to = "basis_function", 
                                        values_to = "value")
    
    # Extracts weighted basis function numbers for faceting.
    weighted_basis_long$basis_num <- as.numeric(gsub("weighted_phi_", "", weighted_basis_long$basis_function))
    
    # Creates a grid of individual unscaled weighted basis function plots.
    p_weighted_basis_grid_unscaled <- ggplot(weighted_basis_long, aes(x = TRPM4, y = value)) +
        geom_line(aes(color = factor(basis_num)), size = 1) +
        scale_color_manual(values = basis_colors[1:n_basis],
                           labels = paste0("φ", 1:n_basis)) +
        facet_wrap(~ basis_num, ncol = 4) +
        labs(title = "Weighted Individual TPRS Basis Functions (Unscaled)",
             subtitle = "Each basis function multiplied by its coefficient",
             x = "TRPM4 Expression (log2)",
             y = "Weighted Value") +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"),
              strip.text = element_text(size = 12, face = "bold", color = "navy"))
    
    # Plots GAM components including data points and knot lines.
    p3 <- ggplot() +
        geom_point(data = data, aes(x = TRPM4, y = Expression), 
                   alpha = 0.15, color = "gray50", size = 2) +
        geom_vline(xintercept = knots, linetype = "dashed", color = knot_color, alpha = 0.5) +
        geom_line(data = basis_data, aes(x = TRPM4, y = smooth_only, 
                                         color = "Smooth Only"), 
                  size = 1) +
        geom_line(data = basis_data, aes(x = TRPM4, y = intercept_plus_linear, 
                                         color = "I+L"), 
                  size = 1) +
        geom_line(data = basis_data, aes(x = TRPM4, y = full_fit, color = "Full Fit"), 
                  size = 1.2) +
        geom_hline(yintercept = intercept, linetype = "dotted", color = "black") +
        annotate("text", x = min(basis_data$TRPM4), y = intercept * 1.02, 
                 label = paste("Intercept =", round(intercept, 2)), hjust = 0) +
        scale_color_manual(values = c("Full Fit" = "#9467BD", 
                                      "I+L" = "#FF7F0E", 
                                      "Smooth Only" = "#2CA02C")) +
        labs(title = paste("GAM Components for", sample_id, "-", gene_set),
             subtitle = "Showing individual components and how they add to form the full fit",
             x = "TRPM4 Expression (log2)", 
             y = "Expression Value") +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold"))
    
    # Initializes a list to store basis function contributions.
    basis_contributions <- list()
    for (i in 1:n_basis) {
        basis_contributions[[i]] <- basis_data[[paste0("weighted_phi_", i)]]
    }
    
    # Computes variance explained by each basis function.
    contribution_var <- sapply(basis_contributions, var)
    total_smooth_var <- var(basis_data$smooth_only)
    
    # Normalizes variance contributions relative to total smooth variance.
    contribution_importance <- contribution_var / total_smooth_var
    
    # Uses absolute coefficient values as an importance measure.
    coefficient_importance <- abs(smooth_coefs)
    
    # Calculates modified eigenvalue importance with regularization.
    S <- temp_model$smooth[[1]]$S[[1]]
    eigen_result <- eigen(S)
    eigen_values <- eigen_result$values[1:n_basis]
    epsilon <- 1e-6
    eigenvalue_importance <- 1/(eigen_values + epsilon)
    eigenvalue_importance <- eigenvalue_importance / max(eigenvalue_importance)
    
    # Computes L2 norm of each basis function as an importance metric.
    l2_norms <- numeric(n_basis)
    for (i in 1:n_basis) {
        l2_norms[i] <- sqrt(sum(basis_data[[paste0("phi_", i)]]^2))
    }
    l2_importance <- l2_norms / max(l2_norms)
    
    # Calculates weighted L2 norm combining coefficient and function energy.
    weighted_l2_importance <- coefficient_importance * l2_importance
    weighted_l2_importance <- weighted_l2_importance / max(weighted_l2_importance)
    
    # Constructs a data frame with various importance measures for basis functions.
    basis_importance <- data.frame(
        basis_function = paste0("φ", 1:n_basis),
        original_index = 1:n_basis,
        eigen_importance = eigenvalue_importance,
        coef_magnitude = coefficient_importance,
        variance_contribution = contribution_importance,
        variance_percentage = contribution_importance / sum(contribution_importance) * 100,
        l2_norm = l2_importance,
        weighted_l2 = weighted_l2_importance
    )
    
    # Orders basis functions by variance contribution and computes cumulative variance.
    basis_importance <- basis_importance[order(-basis_importance$variance_contribution), ]
    basis_importance$cum_variance <- cumsum(basis_importance$variance_percentage)
    
    # Reorders by original index for consistent plotting.
    basis_importance <- basis_importance[order(basis_importance$original_index), ]
    
    # Visualizes variance contribution of basis functions to the smooth component.
    p4a_variance <- ggplot(basis_importance) +
        geom_col(aes(x = original_index, y = variance_percentage, fill = basis_function), alpha = 0.8) +
        geom_text(aes(x = original_index, y = variance_percentage, 
                      label = sprintf("%.1f%%", variance_percentage)),
                  vjust = -0.5, size = 3.5) +
        geom_text(aes(x = original_index, y = 0, 
                      label = basis_function),
                  vjust = 1.5, size = 3, fontface = "bold") +
        scale_fill_manual(values = direct_color_mapping) +
        labs(title = "Variance Contribution of Basis Functions to Smooth Only",
             x = "Basis Function Index (φ1-φ8)", 
             y = "Variance Contribution (%)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Displays coefficient magnitudes of basis functions.
    p4a_coef <- ggplot(basis_importance) +
        geom_col(aes(x = original_index, y = coef_magnitude, fill = basis_function), alpha = 0.8) +
        geom_hline(yintercept = abs(linear_term), 
                   linetype = "dashed", color = "#FF7F0E", size = 1) +
        geom_text(aes(x = original_index, y = coef_magnitude, 
                      label = sprintf("%.3f", coef_magnitude)),
                  vjust = -0.5, size = 3.5) +
        geom_text(aes(x = original_index, y = 0, 
                      label = basis_function),
                  vjust = 1.5, size = 3, fontface = "bold") +
        scale_fill_manual(values = direct_color_mapping) +
        labs(title = "Coefficient Magnitudes of Basis Functions",
             subtitle = "Absolute coefficient values from GAM fit",
             x = "Basis Function Index (φ1-φ8)", 
             y = "Absolute Coefficient Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Shows pure L2 norm importance of basis functions.
    p4a_pure_l2 <- ggplot(basis_importance) +
        geom_col(aes(x = original_index, y = l2_norm, fill = basis_function), alpha = 0.8) +
        geom_text(aes(x = original_index, y = l2_norm, 
                      label = sprintf("%.3f", l2_norm)),
                  vjust = -0.5, size = 3.5) +
        geom_text(aes(x = original_index, y = 0, 
                      label = basis_function),
                  vjust = 1.5, size = 3, fontface = "bold") +
        scale_fill_manual(values = direct_color_mapping) +
        labs(title = "Pure L2 Norm Importance",
             subtitle = "Measures intrinsic 'energy' of each basis function (no coefficient influence)",
             x = "Basis Function Index (φ1-φ8)", 
             y = "L2 Norm (scaled)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Illustrates weighted L2 norm importance combining coefficient and energy.
    p4a_wl2 <- ggplot(basis_importance) +
        geom_col(aes(x = original_index, y = weighted_l2, fill = basis_function), alpha = 0.8) +
        geom_text(aes(x = original_index, y = weighted_l2, 
                      label = sprintf("%.3f", weighted_l2)),
                  vjust = -0.5, size = 3.5) +
        geom_text(aes(x = original_index, y = 0, 
                      label = basis_function),
                  vjust = 1.5, size = 3, fontface = "bold") +
        scale_fill_manual(values = direct_color_mapping) +
        labs(title = "Weighted L2 Norm Importance",
             subtitle = "Combines coefficient magnitude with function 'energy'",
             x = "Basis Function Index (φ1-φ8)", 
             y = "Weighted L2 Norm (scaled)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Prepares components for full model variance contribution analysis.
    full_model_components <- list()
    full_model_components[["I+L"]] <- basis_data$intercept_plus_linear - mean(basis_data$full_fit)
    for (i in 1:n_basis) {
        full_model_components[[paste0("φ", i)]] <- basis_data[[paste0("weighted_phi_", i)]]
    }
    
    # Computes variance contributions of all model components.
    full_var_contributions <- sapply(full_model_components, var)
    full_var_percentages <- full_var_contributions / sum(full_var_contributions) * 100
    
    # Creates a data frame for full model variance contributions.
    full_model_importance <- data.frame(
        component = names(full_var_percentages),
        variance_percentage = full_var_percentages
    )
    
    # Assigns plotting order indices to components.
    full_model_importance$plot_index <- 1:nrow(full_model_importance)
    
    # Defines custom colors for full model components matching previous schemes.
    full_model_colors <- c("I+L" = "#FF7F0E")
    for (i in 1:n_basis) {
        full_model_colors[paste0("φ", i)] <- direct_color_mapping[paste0("φ", i)]
    }
    
    # Visualizes variance contribution to the full model fit.
    p_full_var_contribution <- ggplot(full_model_importance) +
        geom_col(aes(x = plot_index, y = variance_percentage, fill = component), alpha = 0.8) +
        geom_text(aes(x = plot_index, y = variance_percentage, 
                      label = sprintf("%.1f%%", variance_percentage)),
                  vjust = -0.5, size = 3.5) +
        geom_text(aes(x = plot_index, y = 0, 
                      label = component),
                  vjust = 1.5, size = 3, fontface = "bold") +
        scale_fill_manual(values = full_model_colors) +
        scale_x_continuous(breaks = full_model_importance$plot_index,
                           labels = full_model_importance$component) +
        labs(title = "Variance Contribution to Full Fit",
             x = "Model Components", 
             y = "Variance Contribution (%)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Orders basis functions by variance contribution for cumulative plotting.
    ordered_basis <- basis_importance[order(-basis_importance$variance_contribution), ]
    ordered_basis$basis_order <- 1:n_basis
    
    # Plots cumulative variance explained by basis functions.
    p4b <- ggplot(ordered_basis) +
        geom_line(aes(x = basis_order, y = cum_variance), color = "#CC6666", size = 1.2) +
        geom_point(aes(x = basis_order, y = cum_variance), color = "#CC6666", size = 3) +
        geom_text(aes(x = basis_order, y = cum_variance, 
                      label = sprintf("%.1f%%", cum_variance)),
                  vjust = -0.8, hjust = 0.5, size = 3.5) +
        scale_x_continuous(breaks = 1:n_basis, 
                           labels = ordered_basis$basis_function) +
        labs(title = "Cumulative Variance Explained",
             subtitle = "Basis functions ordered by variance contribution",
             x = "Basis Functions (ordered by importance)", 
             y = "Cumulative Variance Explained (%)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combines variance contribution plots into a single page.
    p_importance_page1 <- p4a_variance / p_full_var_contribution +
        plot_layout(heights = c(1, 1)) +
        plot_annotation(theme = theme(plot.title = element_text(face = "bold"))
        )
    
    # Combines L2 norm-based importance plots into a second page.
    p_importance_page2 <- p4a_pure_l2 / p4a_wl2 +
        plot_layout(heights = c(1, 1)) +
        plot_annotation(
            title = "Basis Function Importance: Function Energy Measures",
            subtitle = "L2 norm-based importance measures",
            theme = theme(plot.title = element_text(face = "bold"))
        )
    
    # Returns a list containing all plots and associated data for further analysis.
    return(list(
        basis_functions_plot = p1,
        weighted_basis_plot = p2,
        individual_basis_plots_unscaled = p_basis_grid_unscaled,
        weighted_individual_basis_plots_unscaled = p_weighted_basis_grid_unscaled,
        components_plot = p3,
        variance_contribution_plot = p4a_variance,
        coefficient_importance_plot = p4a_coef,
        pure_l2_plot = p4a_pure_l2,
        weighted_l2_plot = p4a_wl2,
        full_variance_contribution_plot = p_full_var_contribution,
        cumulative_variance_plot = p4b,
        importance_page1 = p_importance_page1,
        importance_page2 = p_importance_page2,
        importance_data = basis_importance,
        knots = knots,
        coefficients = list(
            intercept = intercept,
            linear_term = linear_term,
            smooth_coefs = smooth_coefs
        ),
        full_model_variance_data = full_model_importance
    ))
}

# Loads the patchwork package for multi-plot layouts.
library(patchwork)

# Retrieves the data frame for analysis from pca_results.
dataframe_name <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$gam_data

# Executes the improved TPRS visualization function with specified parameters.
tprs_viz <- visualize_tprs_improved(
    data = dataframe_name, 
    k = 10, 
    lambda = 0.52801317696145,
    gamma = 1.5,
    sample_id = "HYW_4881_Tumor",
    gene_set = "Ribo", 
    predefined_coefficients = list(
        intercept = 5.154747802,
        linear_term = 0.240575063,
        smooth_coefs = c(-1.279918562, -0.654654700, 0.406200567, 0.305860784, 
                         0.258917042, 0.268765202, 0.246497970, -0.800828120)
    )
)

# Displays all generated plots sequentially.
print(tprs_viz$basis_functions_plot)
print(tprs_viz$weighted_basis_plot)
print(tprs_viz$individual_basis_plots_unscaled)
print(tprs_viz$weighted_individual_basis_plots_unscaled)
print(tprs_viz$components_plot)
print(tprs_viz$importance_page1)
print(tprs_viz$importance_page2)

# Outputs the model formula constructed from coefficients.
cat("\nModel formula:\n")
cat(paste0("f(x) = ", tprs_viz$coefficients$intercept, 
           " + ", tprs_viz$coefficients$linear_term, "*x"))
for (i in 1:length(tprs_viz$coefficients$smooth_coefs)) {
    coef_value <- tprs_viz$coefficients$smooth_coefs[i]
    if (coef_value >= 0) {
        cat(paste0(" + ", coef_value, "*φ", i, "(x)"))
    } else {
        cat(paste0(" - ", abs(coef_value), "*φ", i, "(x)"))
    }
}


# =========================================
# 2. Phi1-8 cumulative plot
# =========================================
# Revised function to plot cumulative smooth components with specified colors
plot_cumulative_components_revised <- function(data, k, smooth_coefs, 
                                               show_variance_label = FALSE) {
    # Sets the number of basis functions based on k.
    n_basis <- k - 2
    
    # Creates a sequence of TRPM4 values for visualization.
    x_seq <- seq(min(data$TRPM4), max(data$TRPM4), length.out = 500)
    
    # Fits a temporary GAM model to extract basis functions.
    temp_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                      data = data, method = "REML")
    
    # Generates a prediction matrix for basis function evaluation.
    X_pred <- predict(temp_model, newdata = data.frame(TRPM4 = x_seq), type = "lpmatrix")
    
    # Identifies smooth term columns in the prediction matrix.
    smooth_cols <- grep("s\\(TRPM4", colnames(X_pred))[1:n_basis]
    
    # Initializes a data frame for basis function storage.
    basis_data <- data.frame(TRPM4 = x_seq)
    
    # Extracts individual basis functions into the data frame.
    for (i in 1:n_basis) {
        basis_data[[paste0("phi_", i)]] <- X_pred[, smooth_cols[i]]
    }
    
    # Applies coefficients to compute weighted basis functions.
    for (i in 1:n_basis) {
        basis_data[[paste0("weighted_phi_", i)]] <- smooth_coefs[i] * basis_data[[paste0("phi_", i)]]
    }
    
    # Defines a specific order for cumulative component addition.
    specified_order <- c(1, 8, 2, 3, 6, 4, 7)
    
    # Calculates the total smooth component from all weighted basis functions.
    basis_data$smooth_component <- rowSums(basis_data[, paste0("weighted_phi_", 1:n_basis)])
    
    # Prepares a data frame for cumulative component plotting.
    components_data <- data.frame(TRPM4 = basis_data$TRPM4)
    
    # Adds initial component combining phi1 and phi8.
    components_data$component1 <- basis_data[["weighted_phi_1"]] + basis_data[["weighted_phi_8"]]
    
    # Adds phi2 to the cumulative component.
    components_data$component2 <- components_data$component1 + basis_data[["weighted_phi_2"]]
    
    # Adds phi3 to the cumulative component.
    components_data$component3 <- components_data$component2 + basis_data[["weighted_phi_3"]]
    
    # Adds phi6 to the cumulative component.
    components_data$component4 <- components_data$component3 + basis_data[["weighted_phi_6"]]
    
    # Adds phi4 to the cumulative component.
    components_data$component5 <- components_data$component4 + basis_data[["weighted_phi_4"]]
    
    # Adds phi7 to the cumulative component.
    components_data$component6 <- components_data$component5 + basis_data[["weighted_phi_7"]]
    
    # Includes all basis functions for the final component.
    components_data$all_components <- basis_data$smooth_component
    
    # Converts cumulative data to long format for plotting.
    components_long <- pivot_longer(
        components_data, 
        cols = c("component1", "component2", "component3", "component4", 
                 "component5", "component6", "all_components"),
        names_to = "component", 
        values_to = "value"
    )
    
    # Assigns descriptive labels to each cumulative component.
    label1 <- "φ1+φ8"
    label2 <- "φ1+φ8+φ2"
    label3 <- "φ1+φ8+φ2+φ3"
    label4 <- "φ1+φ8+φ2+φ3+φ6"
    label5 <- "φ1+φ8+φ2+φ3+φ6+φ4"
    label6 <- "φ1+φ8+φ2+φ3+φ6+φ4+φ7"
    label7 <- "All φ combined"
    
    # Orders components as factors for consistent plotting.
    components_long$component <- factor(
        components_long$component,
        levels = c("component1", "component2", "component3", "component4", 
                   "component5", "component6", "all_components"),
        labels = c(label1, label2, label3, label4, label5, label6, label7)
    )
    
    # Defines colors for cumulative components, with the final one distinct.
    cumulative_colors <- c(
        "#E0E0E0", "#E0E0E0", "#E0E0E0", "#E0E0E0", "#E0E0E0", "#E0E0E0", "#2CA02C"
    )
    names(cumulative_colors) <- c(label1, label2, label3, label4, label5, label6, label7)
    
    # Computes variance of the full smooth component.
    var_all <- var(components_data$all_components)
    
    # Calculates percentage variance explained by each cumulative component.
    var_percentages <- c(
        var(components_data$component1) / var_all * 100,
        var(components_data$component2) / var_all * 100,
        var(components_data$component3) / var_all * 100,
        var(components_data$component4) / var_all * 100,
        var(components_data$component5) / var_all * 100,
        var(components_data$component6) / var_all * 100,
        100
    )
    
    # Creates a cumulative components plot with specified colors.
    p_cumulative <- ggplot(components_long, aes(x = TRPM4, y = value, color = component, group = component)) +
        geom_line(size = 1.2) +
        scale_color_manual(values = cumulative_colors) +
        labs(
            title = "Cumulative Smooth Components by Specified Order",
            x = "TRPM4 Expression (log2)",
            y = "Contribution to Smooth Component"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(face = "bold")
        )
    
    # Formats variance percentages for optional annotation.
    variance_text <- sprintf(
        "Variance explained: %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%, 100%%", 
        var_percentages[1], var_percentages[2], var_percentages[3], 
        var_percentages[4], var_percentages[5], var_percentages[6]
    )
    
    # Adds variance explanation text to the plot if requested.
    if (show_variance_label) {
        p_cumulative <- p_cumulative +
            annotate(
                "text", x = min(components_data$TRPM4), y = max(components_data$all_components),
                label = variance_text,
                hjust = 0, vjust = 1, fontface = "bold"
            )
    }
    
    # Returns the plot, variance percentages, and text for further use.
    return(list(
        plot = p_cumulative,
        variance_explained = var_percentages,
        variance_text = variance_text
    ))
}

# Generates the revised cumulative plot without variance labels.
p_revised <- plot_cumulative_components_revised(
    data = dataframe_name,
    k = 10,
    smooth_coefs = c(-1.279918562, -0.654654700, 0.406200567, 0.305860784, 
                     0.258917042, 0.268765202, 0.246497970, -0.800828120),
    show_variance_label = FALSE
)

# Displays the cumulative components plot.
print(p_revised$plot)

# Outputs variance percentages for each cumulative component.
cat("\nVariance explained by each cumulative component:\n")
component_names <- c("φ1+φ8", "φ1+φ8+φ2", "φ1+φ8+φ2+φ3", "φ1+φ8+φ2+φ3+φ6", 
                     "φ1+φ8+φ2+φ3+φ6+φ4", "φ1+φ8+φ2+φ3+φ6+φ4+φ7", "All φ combined")
for (i in 1:length(p_revised$variance_explained)) {
    cat(sprintf("%s: %.1f%%\n", component_names[i], p_revised$variance_explained[i]))
}

# Prints variance text for potential later addition to the plot.
cat("\nVariance text for adding to the plot later:\n")
cat(p_revised$variance_text)


# =============================================================
# 3. Visualization of GAM components using concrete examples from the data
# =============================================================
# Function
plot_gam_with_examples <- function(data, k = 10, lambda = 0.52801317696145, 
                                   gamma = 1.5, sample_id = "HYW_4881_Tumor", 
                                   gene_set = "Ribo", predefined_coefficients = NULL) {
    
    # Ensures reproducibility of random processes.
    set.seed(123)
    
    # Extracts predefined coefficients for model consistency.
    intercept <- predefined_coefficients$intercept
    linear_term <- predefined_coefficients$linear_term
    
    # Creates a sequence of TRPM4 values for smooth visualization.
    x_seq <- seq(min(data$TRPM4), max(data$TRPM4), length.out = 500)
    
    # Fits the GAM model using REML with specified parameters.
    model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                 data = data, method = "REML", sp = lambda, gamma = gamma)
    
    # Initializes a data frame for plotting model predictions.
    plot_data <- data.frame(TRPM4 = x_seq)
    
    # Predicts the full model fit across the TRPM4 sequence.
    plot_data$full_model <- predict(model, newdata = plot_data)
    
    # Computes the intercept plus linear component using predefined coefficients.
    plot_data$intercept_plus_linear <- intercept + linear_term * plot_data$TRPM4
    
    # Derives the smooth component as the difference from the full fit.
    plot_data$smooth_only <- plot_data$full_model - plot_data$intercept_plus_linear
    
    # Selects example TRPM4 values at low, medium, and high percentiles.
    example_points <- data.frame(
        low_x = quantile(data$TRPM4, 0.05),
        mid_x = median(data$TRPM4),
        high_x = quantile(data$TRPM4, 0.95)
    )
    
    # Identifies data points closest to the selected example values.
    example_indices <- c(
        which.min(abs(data$TRPM4 - example_points$low_x)),
        which.min(abs(data$TRPM4 - example_points$mid_x)),
        which.min(abs(data$TRPM4 - example_points$high_x))
    )
    
    # Extracts the corresponding example data subset.
    example_data <- data[example_indices, ]
    
    # Calculates model components for the example points.
    example_data$intercept_plus_linear <- intercept + linear_term * example_data$TRPM4
    example_data$full_model <- predict(model, newdata = example_data)
    example_data$smooth_only <- example_data$full_model - example_data$intercept_plus_linear
    
    # Computes residuals as the difference between observed and predicted values.
    example_data$residual <- example_data$Expression - example_data$full_model
    
    # Creates the main plot showing data, components, and example points.
    p1 <- ggplot() +
        geom_point(data = data, aes(x = TRPM4, y = Expression), 
                   alpha = 0.25, color = "gray50", size = 2) +
        geom_line(data = plot_data, aes(x = TRPM4, y = intercept_plus_linear, 
                                        color = "Intercept + Linear"), 
                  size = 1) +
        geom_line(data = plot_data, aes(x = TRPM4, y = smooth_only, 
                                        color = "Smooth Component"), 
                  size = 1) +
        geom_line(data = plot_data, aes(x = TRPM4, y = full_model, 
                                        color = "Full Model Fit"), 
                  size = 1.2) +
        geom_hline(yintercept = intercept, linetype = "dotted", color = "black") +
        geom_point(data = example_data, aes(x = TRPM4, y = Expression), 
                   color = "black", size = 4, shape = 21, fill = "yellow") +
        scale_color_manual(values = c("Full Model Fit" = "#9467BD", 
                                      "Intercept + Linear" = "#FF7F0E", 
                                      "Smooth Component" = "#2CA02C")) +
        geom_text(data = example_data, 
                  aes(x = TRPM4, y = Expression, 
                      label = c("Example 1", "Example 2", "Example 3")),
                  hjust = -0.2, vjust = -0.5, size = 4) +
        labs(title = paste("GAM Components for", sample_id, "-", gene_set),
             subtitle = "True mathematical relationships with concrete examples",
             x = "TRPM4 Expression (log2)", 
             y = "Expression Value") +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank(),
              plot.title = element_text(face = "bold"))
    
    # Initializes a list to store breakdown plots for each example.
    example_breakdowns <- list()
    
    for (i in 1:nrow(example_data)) {
        # Extracts values for the current example point.
        example_x <- example_data$TRPM4[i]
        example_y <- example_data$Expression[i]
        example_full <- example_data$full_model[i]
        example_linear <- example_data$intercept_plus_linear[i]
        example_smooth <- example_data$smooth_only[i]
        
        # Prepares data for a waterfall chart of component contributions.
        breakdown_data <- data.frame(
            Component = factor(c("Intercept", "Linear Term", "Smooth Component", "Full Model"),
                               levels = c("Intercept", "Linear Term", "Smooth Component", "Full Model")),
            Value = c(intercept, linear_term * example_x, example_smooth, example_full),
            Cumulative = c(intercept, 
                           intercept + linear_term * example_x, 
                           intercept + linear_term * example_x + example_smooth,
                           example_full)
        )
        
        # Creates a breakdown plot for the current example.
        p_breakdown <- ggplot(breakdown_data, aes(x = Component, y = Value, fill = Component)) +
            geom_col() +
            geom_text(aes(label = sprintf("%.2f", Value)), vjust = ifelse(breakdown_data$Value >= 0, -0.5, 1.5)) +
            scale_fill_manual(values = c("Intercept" = "#999999", 
                                         "Linear Term" = "#FF7F0E", 
                                         "Smooth Component" = "#2CA02C",
                                         "Full Model" = "#9467BD")) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black") +
            labs(title = paste("Example", i, "Breakdown"),
                 subtitle = paste("TRPM4 =", round(example_x, 2), ", Actual Expression =", round(example_y, 2)),
                 x = "", y = "Contribution to Prediction") +
            theme_minimal() +
            theme(legend.position = "none")
        
        example_breakdowns[[i]] <- p_breakdown
    }
    
    # Prepares data for a dot plot of component values across examples.
    component_dot_data <- data.frame(
        Example = rep(paste("Example", 1:3), each = 3),
        Component = rep(c("Intercept + Linear", "Smooth Component", "Full Model"), 3),
        Value = c(
            example_data$intercept_plus_linear[1], example_data$smooth_only[1], example_data$full_model[1],
            example_data$intercept_plus_linear[2], example_data$smooth_only[2], example_data$full_model[2],
            example_data$intercept_plus_linear[3], example_data$smooth_only[3], example_data$full_model[3]
        ),
        TRPM4 = rep(example_data$TRPM4, each = 3)
    )
    
    # Orders component factors for consistent plotting.
    component_dot_data$Component <- factor(component_dot_data$Component, 
                                           levels = c("Intercept + Linear", "Smooth Component", "Full Model"))
    
    # Creates a dot plot showing component contributions for each example.
    p_dots <- ggplot(component_dot_data, aes(x = Component, y = Value, color = Component)) +
        geom_point(size = 4) +
        geom_text(aes(label = sprintf("%.2f", Value)), hjust = -0.3, size = 3.5) +
        scale_color_manual(values = c("Full Model" = "#9467BD", 
                                      "Intercept + Linear" = "#FF7F0E", 
                                      "Smooth Component" = "#2CA02C")) +
        facet_wrap(~Example, ncol = 3) +
        labs(title = "Component Values at Example Points",
             subtitle = "Shows how components combine to form the full model prediction",
             x = "", y = "Value") +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Initializes a data frame for component flow visualization.
    example_flow_data <- data.frame()
    
    for (i in 1:nrow(example_data)) {
        example_label <- paste("Example", i)
        
        # Prepares data for flow from intercept to intercept plus linear.
        flow1 <- data.frame(
            Example = example_label,
            x = c("Intercept", "Plus Linear"),
            y = c(intercept, example_data$intercept_plus_linear[i]),
            Component = "Flow 1"
        )
        
        # Prepares data for flow from intercept plus linear to full model.
        flow2 <- data.frame(
            Example = example_label,
            x = c("Plus Linear", "Full Model"),
            y = c(example_data$intercept_plus_linear[i], example_data$full_model[i]),
            Component = "Flow 2"
        )
        
        # Combines flow data into the main data frame.
        example_flow_data <- rbind(example_flow_data, flow1, flow2)
    }
    
    # Creates a flow plot showing component contributions across examples.
    p_flow <- ggplot(example_flow_data, aes(x = x, y = y, color = Component, group = interaction(Example, Component))) +
        geom_line(aes(linetype = Component), size = 1) +
        geom_point(size = 3) +
        geom_text(aes(label = sprintf("%.2f", y)), vjust = -1, size = 3) +
        scale_color_manual(values = c("Flow 1" = "#FF7F0E", "Flow 2" = "#2CA02C")) +
        scale_linetype_manual(values = c("Flow 1" = "solid", "Flow 2" = "dashed")) +
        facet_wrap(~Example, ncol = 3) +
        labs(title = "Component Flow for Each Example",
             subtitle = "Showing how the value changes as components are added",
             x = "", y = "Value") +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Constructs a detailed textual explanation of example contributions.
    explanation_text <- paste(
        "Example 1 (Low TRPM4):",
        sprintf("- TRPM4 = %.2f", example_data$TRPM4[1]),
        sprintf("- Intercept + Linear = %.2f", example_data$intercept_plus_linear[1]),
        sprintf("- Smooth Component = %.2f", example_data$smooth_only[1]),
        sprintf("- Full Model = %.2f", example_data$full_model[1]),
        sprintf("- Actual Expression = %.2f", example_data$Expression[1]),
        "",
        "Example 2 (Medium TRPM4):",
        sprintf("- TRPM4 = %.2f", example_data$TRPM4[2]),
        sprintf("- Intercept + Linear = %.2f", example_data$intercept_plus_linear[2]),
        sprintf("- Smooth Component = %.2f", example_data$smooth_only[2]),
        sprintf("- Full Model = %.2f", example_data$full_model[2]),
        sprintf("- Actual Expression = %.2f", example_data$Expression[2]),
        "",
        "Example 3 (High TRPM4):",
        sprintf("- TRPM4 = %.2f", example_data$TRPM4[3]),
        sprintf("- Intercept + Linear = %.2f", example_data$intercept_plus_linear[3]),
        sprintf("- Smooth Component = %.2f", example_data$smooth_only[3]),
        sprintf("- Full Model = %.2f", example_data$full_model[3]),
        sprintf("- Actual Expression = %.2f", example_data$Expression[3]),
        "",
        "Key Observations:",
        "1. At low TRPM4, the smooth component is negative, pulling the prediction below the linear component.",
        "2. At medium TRPM4, the smooth component becomes less negative/more positive.",
        "3. At high TRPM4, the smooth component becomes even more positive, pushing the prediction above the linear component.",
        sep = "\n"
    )
    
    # Creates a summary table of component values for all examples.
    summary_table <- data.frame(
        Example = paste("Example", 1:3),
        TRPM4 = example_data$TRPM4,
        Intercept = intercept,
        Linear_Term = linear_term * example_data$TRPM4,
        Intercept_plus_Linear = example_data$intercept_plus_linear,
        Smooth_Component = example_data$smooth_only,
        Full_Model = example_data$full_model,
        Actual_Expression = example_data$Expression
    )
    
    # Returns all plots, data, and explanations for further use.
    return(list(
        main_plot = p1,
        breakdowns = example_breakdowns,
        dot_plot = p_dots,
        flow_plot = p_flow,
        example_data = example_data,
        explanation = explanation_text,
        summary_table = summary_table
    ))
}

# Executes the GAM visualization with example points.
examples_viz <- plot_gam_with_examples(
    dataframe_name, 
    k = 10, 
    lambda = 0.52801317696145,
    gamma = 1.5,
    predefined_coefficients = list(intercept = 5.154747802, linear_term = 0.240575063, smooth_coefs = c(-1.279918562, -0.654654700, 0.406200567, 0.305860784, 0.258917042, 0.268765202, 0.246497970, -0.800828120))
)

# Displays the main GAM components plot with examples.
print(examples_viz$main_plot)

# Shows breakdown plots for each example point.
print(examples_viz$breakdowns[[1]])
print(examples_viz$breakdowns[[2]])
print(examples_viz$breakdowns[[3]])

# Displays the dot plot of component values.
print(examples_viz$dot_plot)

# Shows the flow plot of component contributions.
print(examples_viz$flow_plot)

# Outputs the textual explanation of example contributions.
cat(examples_viz$explanation)

# Prints the summary table of component values.
print(examples_viz$summary_table)


# ================================
# 4. Building TPRS's phi visualization
# ================================
# Specifies the sample data for TPRS basis function analysis.
sample_data <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$gam_data

# Visualizes the emergence of TPRS basis function phi1.
visualize_phi1_emergence <- function(sample_data, k = 10) {
    # Loads necessary libraries for plotting and modeling.
    library(ggplot2)
    library(mgcv)
    library(viridis)
    library(tidyr)
    library(gridExtra)
    
    # Ensures reproducibility of random processes.
    set.seed(123)
    
    # Creates a sequence of TRPM4 values for visualization.
    x_seq <- seq(min(sample_data$TRPM4), max(sample_data$TRPM4), length.out = 500)
    
    # Fits a GAM model with TPRS basis to extract basis functions.
    model <- gam(Expression ~ s(TRPM4, bs = "tp", k = k), 
                 data = sample_data, method = "REML")
    
    # Retrieves knot locations from the fitted model.
    knots <- model$smooth[[1]]$xp
    if (is.null(knots)) {
        knots <- quantile(sample_data$TRPM4, probs = seq(0, 1, length.out = k-2))
    }
    
    # Generates the prediction matrix containing basis functions.
    X_pred <- predict(model, newdata = data.frame(TRPM4 = x_seq), type = "lpmatrix")
    
    # Identifies columns for smooth terms in the prediction matrix.
    smooth_cols <- grep("s\\(TRPM4", colnames(X_pred))
    phi_data <- data.frame(TRPM4 = x_seq)
    
    # Extracts all basis functions into a data frame.
    for(i in 1:length(smooth_cols)) {
        phi_data[paste0("phi", i)] <- X_pred[, smooth_cols[i]]
    }
    
    # Fits a linear model to compute residuals.
    linear_model <- lm(Expression ~ TRPM4, data = sample_data)
    sample_data$linear_pred <- predict(linear_model, newdata = sample_data)
    sample_data$residuals <- sample_data$Expression - sample_data$linear_pred
    linear_pred_seq <- predict(linear_model, newdata = data.frame(TRPM4 = x_seq))
    
    # Smooths residuals using a loess fit for visualization.
    smoother <- loess(residuals ~ TRPM4, data = sample_data, span = 0.5)
    residual_values <- predict(smoother, newdata = data.frame(TRPM4 = x_seq))
    flipped_residuals <- -residual_values
    
    # Scales values to match the range of phi1 for comparison.
    scale_to_range <- function(x, target_min, target_max) {
        x_min <- min(x, na.rm = TRUE)
        x_max <- max(x, na.rm = TRUE)
        scaled <- ((x - x_min) / (x_max - x_min)) * (target_max - target_min) + target_min
        return(scaled)
    }
    
    # Applies scaling to flipped residuals based on phi1 range.
    phi1_range <- range(phi_data$phi1)
    flipped_scaled <- scale_to_range(flipped_residuals, phi1_range[1], phi1_range[2])
    
    # Plots the initial linear fit to the data.
    p1 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = Expression), alpha = 0.3, color = "gray40") +
        geom_line(data = data.frame(TRPM4 = x_seq, value = linear_pred_seq), 
                  aes(x = TRPM4, y = value), color = "#FFCC99", size = 1.2) +
        labs(title = "1. Linear Fit",
             subtitle = "Initial attempt to capture overall data shape",
             x = "TRPM4 Expression (log2)", 
             y = "Expression") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Visualizes residuals from the linear fit.
    p2 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = residuals), alpha = 0.3, color = "gray40") +
        geom_line(data = data.frame(TRPM4 = x_seq, value = residual_values), 
                  aes(x = TRPM4, y = value), color = "#F25B5B", size = 1.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", alpha = 0.5) +
        labs(title = "2. Actual Residual Plot",
             subtitle = "Residuals from linear fit",
             x = "TRPM4 Expression (log2)", 
             y = "Residual Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Shows flipped residuals with knot locations.
    p3 <- ggplot() +
        geom_line(data = data.frame(TRPM4 = x_seq, value = flipped_residuals), 
                  aes(x = TRPM4, y = value), color = "#F25B5B", size = 1.2) +
        geom_point(data = sample_data, aes(x = TRPM4, y = -residuals), 
                   alpha = 0.3, color = "gray40") +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        labs(title = "3. Initial Flipped Residual Pattern",
             subtitle = "Starting point for basis function development",
             x = "TRPM4 Expression (log2)", 
             y = "Flipped Residual Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Defines colors for basis functions using the mako palette.
    phi_colors <- viridis::mako(8, begin = 0.1, end = 0.9)
    names(phi_colors) <- paste0("phi", 1:8)
    
    # Plots initial bell-shaped basis functions centered at knots.
    transition_data_100_bell <- data.frame()
    for (i in 1:min(8, length(knots))) {
        width <- diff(range(sample_data$TRPM4)) / (1.5 * length(knots))
        curve_x <- x_seq
        bell_curve <- dnorm(curve_x, mean = knots[i], sd = width)
        bell_curve <- scale_to_range(bell_curve, 0, max(abs(flipped_residuals)) * 0.8)
        trans_df <- data.frame(
            TRPM4 = curve_x,
            value = bell_curve,
            basis_function = paste0("phi", i)
        )
        transition_data_100_bell <- rbind(transition_data_100_bell, trans_df)
    }
    
    p4 <- ggplot() +
        geom_line(data = transition_data_100_bell, aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "4. 100% Bell Shape",
             subtitle = "Initial localized basis functions",
             x = "TRPM4 Expression (log2)", 
             y = "Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Visualizes final unpenalized basis function shapes.
    transition_data_100_final <- data.frame()
    for (i in 1:min(8, length(knots))) {
        curve_x <- x_seq
        final_shape <- phi_data[, paste0("phi", i)]
        trans_df <- data.frame(
            TRPM4 = curve_x,
            value = final_shape,
            basis_function = paste0("phi", i)
        )
        transition_data_100_final <- rbind(transition_data_100_final, trans_df)
    }
    
    p5 <- ggplot() +
        geom_line(data = transition_data_100_final, aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "5. Final Unpenalized Shapes",
             subtitle = "Computed basis functions used in the model",
             x = "TRPM4 Expression (log2)", 
             y = "Basis Function Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "right",
              legend.title = element_blank())
    
    # Compiles all plots into a list for arrangement.
    all_plots <- list(p1, p2, p3, p4, p5)
    all_plots <- all_plots[!sapply(all_plots, is.null)]
    
    # Arranges all plots vertically in a single layout.
    final_plot <- do.call(arrangeGrob, c(all_plots, list(ncol = 1)))
    
    # Returns individual plots and the combined layout.
    return(list(
        linear_fit = p1,
        actual_residuals = p2,
        flipped_residual = p3,
        transition_100_bell = p4,
        basis_final = p5,
        all_plots = final_plot
    ))
}

# Executes the phi1 emergence visualization function.
results <- visualize_phi1_emergence(sample_data, k = 10)

# Displays each plot individually.
print(results$linear_fit)
print(results$actual_residuals)
print(results$flipped_residual)
print(results$transition_100_bell)
print(results$basis_final)
