###############################################################
# Part 3.7: Consequences of Different Gamma Values on GAM Fitting
###############################################################
# Specifies the dataset for analysis from the PCA results.
dataframe_name <- pca_results[["HYW_4701_Tumor"]][["Ribo"]]$gam_data

# Visualizes GAM components across multiple gamma values with free lambda selection.
visualize_gam_gammas_free <- function(data, k = 6, 
                                      gammas = c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), 
                                      sample_id = "HYW_4701_Tumor", 
                                      gene_set = "Ribo") {
    
    # Verifies required libraries are installed before proceeding.
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required but not installed")
    }
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("mgcv package is required but not installed")
    }
    if (!requireNamespace("viridis", quietly = TRUE)) {
        stop("viridis package is required but not installed")
    }
    
    # Loads necessary libraries for modeling and visualization.
    library(ggplot2)
    library(mgcv)
    library(viridis)
    library(tidyr)
    
    # Sets a seed for reproducible model fitting and predictions.
    set.seed(123)
    
    # Generates a sequence of TRPM4 values for smooth curve visualization.
    x_seq <- seq(min(data$TRPM4), max(data$TRPM4), length.out = 500)
    
    # Calculates the number of meaningful basis functions based on k.
    n_basis <- k - 2
    
    # Initializes a data frame to store predictions for all gamma values.
    all_predictions <- data.frame(TRPM4 = x_seq)
    
    # Creates a vector to store selected lambda values for each gamma.
    selected_lambdas <- numeric(length(gammas))
    names(selected_lambdas) <- paste0("gamma_", gammas)
    
    # Fits GAM models and generates predictions for each gamma value.
    for (i in seq_along(gammas)) {
        gamma <- gammas[i]
        tryCatch({
            # Fits a GAM model with the current gamma, allowing REML to select lambda.
            model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                         data = data, method = "REML", gamma = gamma)
            
            # Records the lambda value selected by REML for this gamma.
            selected_lambdas[i] <- model$sp[1]
            
            # Adds model predictions to the data frame for visualization.
            col_name <- paste0("gamma_", gamma)
            all_predictions[[col_name]] <- predict(model, newdata = data.frame(TRPM4 = x_seq))
            
            # Outputs a confirmation of successful model fitting with selected lambda.
            cat("Successfully fit model for gamma =", gamma, 
                "- Selected lambda:", selected_lambdas[i], "\n")
            
        }, error = function(e) {
            cat("Error fitting model for gamma =", gamma, ":", e$message, "\n")
        })
    }
    
    # Fits a reference GAM model with gamma = 1 to extract knot locations.
    ref_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                     data = data, method = "REML", gamma = 1)
    
    # Extracts knot positions from the reference model or approximates if unavailable.
    knots <- ref_model$smooth[[1]]$xp
    if (is.null(knots) || length(knots) == 0) {
        n_knots <- k - 2
        knots <- quantile(data$TRPM4, probs = seq(0, 1, length.out = n_knots))
    }
    
    # Verifies the number of gamma values with successful predictions.
    gamma_cols <- grep("^gamma_", names(all_predictions), value = TRUE)
    cat("Generated predictions for", length(gamma_cols), "out of", length(gammas), "gamma values\n")
    
    # Defines a custom viridis color palette with 1-99% coverage for visualization.
    viridis_palette <- viridis::viridis(100)[1:99]
    
    # Assigns evenly spaced colors from the palette to gamma values.
    gamma_colors <- viridis_palette[seq(1, length(viridis_palette), length.out = length(gammas))]
    
    # Sets gamma labels for consistent legend display.
    gamma_labels <- paste0("", gammas)
    
    # Defines a darker gray color for knot lines in the plot.
    knot_color <- "#333333"
    
    # Reshapes prediction data into long format for plotting.
    long_predictions <- pivot_longer(
        all_predictions, 
        cols = starts_with("gamma_"),
        names_to = "gamma", 
        values_to = "prediction"
    )
    
    # Extracts numeric gamma values and creates a factor for ordered plotting.
    long_predictions$gamma_value <- as.numeric(gsub("gamma_", "", long_predictions$gamma))
    long_predictions$gamma_factor <- factor(
        long_predictions$gamma_value,
        levels = gammas,
        labels = gamma_labels
    )
    
    # Orders data to draw smaller gamma values last, ensuring visibility.
    long_predictions <- long_predictions[order(-long_predictions$gamma_value),]
    
    # Creates the primary plot showing full GAM fits for different gamma values.
    p_main <- ggplot() +
        # Adds data points with reduced opacity for context.
        geom_point(data = data, aes(x = TRPM4, y = Expression), 
                   alpha = 0.15, color = "gray50", size = 2) +
        
        # Includes vertical lines to indicate knot positions.
        geom_vline(xintercept = knots, linetype = "dashed", color = knot_color, alpha = 0.5) +
        
        # Plots smooth curves for each gamma value with distinct colors.
        geom_line(data = long_predictions, 
                  aes(x = TRPM4, y = prediction, color = gamma_factor, group = gamma), 
                  size = 1.2) +
        
        # Applies the custom color palette to gamma values.
        scale_color_manual(values = setNames(gamma_colors, gamma_labels)) +
        
        # Sets plot titles and labels for clarity.
        labs(title = paste("GAM Full Fits with Different Gamma Values for", sample_id, "-", gene_set),
             subtitle = paste("k =", k, "with free lambda selection"),
             x = "TRPM4 expression (log2)", 
             y = "Expression Value",
             color = "gamma") +
        
        # Configures a minimal theme with a right-aligned legend.
        theme_minimal() +
        theme(legend.position = "right",
              plot.title = element_text(face = "bold"))
    
    # Retrieves reference predictions for gamma = 1 if available.
    ref_col <- "gamma_1"
    if (ref_col %in% names(all_predictions)) {
        ref_predictions <- all_predictions[[ref_col]]
        
        # Initializes a data frame to store differences from the reference.
        diff_predictions <- data.frame(TRPM4 = x_seq)
        
        # Calculates prediction differences for each gamma relative to gamma = 1.
        for (gamma in gammas) {
            col_name <- paste0("gamma_", gamma)
            if (col_name %in% names(all_predictions)) {
                diff_col_name <- paste0("diff_", gamma)
                diff_predictions[[diff_col_name]] <- all_predictions[[col_name]] - ref_predictions
            }
        }
        
        # Reshapes difference data into long format for plotting.
        long_diff <- pivot_longer(
            diff_predictions,
            cols = starts_with("diff_"),
            names_to = "gamma",
            values_to = "difference"
        )
        
        # Extracts gamma values and creates a factor for ordered plotting.
        long_diff$gamma_value <- as.numeric(gsub("diff_", "", long_diff$gamma))
        long_diff$gamma_factor <- factor(
            long_diff$gamma_value,
            levels = gammas,
            labels = gamma_labels
        )
        
        # Creates a plot showing differences from the gamma = 1 reference.
        p_diff <- ggplot() +
            geom_hline(yintercept = 0, linetype = "solid", color = "gray50") +
            geom_line(data = long_diff, 
                      aes(x = TRPM4, y = difference, color = gamma_factor, group = gamma),
                      size = 0.8) +
            scale_color_manual(values = setNames(gamma_colors, gamma_labels)) +
            labs(title = paste("Difference from γ = 1 for Different Gamma Values"),
                 subtitle = paste("k =", k, "with free lambda selection"),
                 x = "TRPM4 expression (log2)",
                 y = "Difference from γ = 1",
                 color = "gamma") +
            theme_minimal() +
            theme(legend.position = "right",
                  plot.title = element_text(face = "bold"))
    } else {
        p_diff <- NULL
    }
    
    # Prepares a data frame for plotting selected lambda values across gammas.
    selected_lambdas_df <- data.frame(
        gamma = as.numeric(gsub("gamma_", "", names(selected_lambdas))),
        lambda = selected_lambdas
    )
    
    # Orders the data frame by gamma value for consistent plotting.
    selected_lambdas_df <- selected_lambdas_df[order(selected_lambdas_df$gamma), ]
    
    # Adds gamma factor labels for bar plot consistency.
    selected_lambdas_df$gamma_factor <- factor(
        selected_lambdas_df$gamma,
        levels = gammas,
        labels = gamma_labels
    )
    
    # Creates a bar plot of selected lambda values for each gamma.
    p_lambda <- ggplot(selected_lambdas_df, aes(x = gamma_factor, y = lambda, fill = gamma_factor)) +
        geom_col() +
        geom_text(aes(label = sprintf("%.3f", lambda)), vjust = -0.5, size = 3) +
        scale_fill_manual(values = setNames(gamma_colors, gamma_labels)) +
        labs(x = "Gamma",
             y = "Selected lambda value") +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Returns a list containing all plots and lambda data for further use.
    return(list(
        main_plot = p_main, 
        difference_plot = p_diff, 
        lambda_plot = p_lambda,
        selected_lambdas = selected_lambdas_df
    ))
}

# Executes the visualization function with specified parameters.
result <- visualize_gam_gammas_free(
    data = dataframe_name, 
    k = 6, 
    gammas = c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
    sample_id = "HYW_4701_Tumor",
    gene_set = "Ribo"
)

# Displays the main GAM fit plot.
print(result$main_plot)

# Displays the difference plot if it exists.
if (!is.null(result$difference_plot)) {
    print(result$difference_plot)
}

# Displays the bar plot of selected lambda values.
print(result$lambda_plot)

# Outputs the selected lambda values for review.
print(result$selected_lambdas)

# Specifies dimensions for line and bar plots.
# Lineplot 500 x 370
# Barplot 390 x 370
# CV_RMSE 500 x 430

# Evaluates gamma effects on RMSE using k-fold cross-validation and boxplot visualization.
evaluate_gamma_boxplot <- function(data, k_cv = 10, k_basis = 6, 
                                   gammas = c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5),
                                   ref_gamma = 1) {
    
    # Loads required libraries for modeling, visualization, and statistical analysis.
    library(mgcv)
    library(ggplot2)
    library(dplyr)
    library(caret)
    library(tidyr)
    library(viridis)
    
    # Sets a seed for reproducible fold generation and model fitting.
    set.seed(123)
    
    # Creates k-fold splits for cross-validation based on the response variable.
    folds <- createFolds(data$Expression, k = k_cv, list = TRUE, returnTrain = FALSE)
    
    # Initializes a matrix to store RMSE values for each gamma and fold.
    cv_results <- matrix(NA, nrow = length(gammas), ncol = k_cv)
    rownames(cv_results) <- paste0("gamma_", gammas)
    
    # Performs k-fold cross-validation for each gamma value.
    for (i in seq_along(gammas)) {
        gamma_val <- gammas[i]
        
        # Iterates over folds to compute RMSE for the current gamma.
        for (j in seq_along(folds)) {
            test_indices <- folds[[j]]
            train_data <- data[-test_indices, ]
            test_data <- data[test_indices, ]
            
            # Fits a GAM model on the training data with the specified gamma.
            model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k_basis), 
                         data = train_data, method = "REML", gamma = gamma_val)
            
            # Generates predictions on the test data.
            predictions <- predict(model, newdata = test_data)
            
            # Calculates RMSE for the current fold and stores it.
            cv_results[i, j] <- sqrt(mean((test_data$Expression - predictions)^2))
        }
    }
    
    # Converts the RMSE matrix to a data frame for easier handling.
    cv_df <- as.data.frame(cv_results)
    
    # Adds gamma values as a column for reference.
    cv_df$gamma <- gammas
    
    # Prepares to test differences from the reference gamma value.
    ref_row <- which(gammas == ref_gamma)
    p_values <- numeric(length(gammas))
    
    # Conducts statistical tests comparing each gamma to the reference.
    for (i in seq_along(gammas)) {
        if (i == ref_row) {
            p_values[i] <- 1.0
            next
        }
        
        # Applies a Wilcoxon test for paired differences, with fallback to t-test.
        tryCatch({
            test_result <- wilcox.test(as.numeric(cv_results[i, ]), 
                                       as.numeric(cv_results[ref_row, ]), 
                                       paired = TRUE)
            p_values[i] <- test_result$p.value
        }, error = function(e) {
            tryCatch({
                test_result <- t.test(as.numeric(cv_results[i, ]), 
                                      as.numeric(cv_results[ref_row, ]), 
                                      paired = TRUE)
                p_values[i] <- test_result$p.value
            }, error = function(e2) {
                p_values[i] <- 0.5
            })
        })
    }
    
    # Adjusts p-values using the Benjamini-Hochberg method for multiple comparisons.
    q_values <- p.adjust(p_values, method = "BH")
    
    # Prepares labels for the plot, marking the reference gamma.
    plot_labels <- sprintf("%.3f", q_values)
    plot_labels[ref_row] <- "ref"
    
    # Reshapes RMSE data into long format for boxplot visualization.
    long_data <- tidyr::pivot_longer(
        cv_df, 
        cols = 1:k_cv,
        names_to = "fold", 
        values_to = "rmse"
    )
    
    # Adds a factor column for gamma to ensure proper ordering in the plot.
    long_data$gamma_factor <- factor(long_data$gamma)
    
    # Defines a viridis color palette for consistent visualization.
    viridis_palette <- viridis::viridis(length(gammas))
    
    # Sets a y-position for labels above the maximum RMSE value.
    y_pos_for_labels <- max(long_data$rmse, na.rm = TRUE) + 0.002
    
    # Creates a data frame for plotting q-value labels.
    label_data <- data.frame(
        gamma_factor = factor(gammas),
        q_value = q_values,
        label = plot_labels,
        y_pos = y_pos_for_labels,
        stringsAsFactors = FALSE
    )
    
    # Creates a boxplot to visualize RMSE across gamma values.
    p <- ggplot(long_data, aes(x = gamma_factor, y = rmse, fill = gamma_factor)) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#F5F5F5"),
            panel.grid.minor = element_line(color = "#FAFAFA"),
            plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 0),
            text = element_text(color = "black"),
            legend.position = "none"
        ) +
        # Draws boxplots with customized appearance.
        geom_boxplot(outlier.shape = 16, outlier.size = 2, width = 0.6,
                     color = "black", coef = 5, fatten = 2) +
        # Adds horizontal caps to boxplot whiskers.
        geom_segment(data = long_data %>% 
                       group_by(gamma_factor) %>% 
                       summarise(min = min(rmse), max = max(rmse)),
                     aes(x = as.numeric(gamma_factor) - 0.2, 
                         xend = as.numeric(gamma_factor) + 0.2,
                         y = max, yend = max), 
                     color = "black") +
        geom_segment(data = long_data %>% 
                       group_by(gamma_factor) %>% 
                       summarise(min = min(rmse), max = max(rmse)),
                     aes(x = as.numeric(gamma_factor) - 0.2, 
                         xend = as.numeric(gamma_factor) + 0.2,
                         y = min, yend = min), 
                     color = "black") +
        # Adds q-value labels above each boxplot.
        geom_text(
            data = label_data,
            aes(x = gamma_factor, y = y_pos, label = label),
            size = 3.5, vjust = -0.5
        ) +
        scale_fill_manual(values = viridis_palette) +
        labs(
            title = "Cross-Validation RMSE for Different Gamma Values",
            x = "Gamma value",
            y = "RMSE"
        )
    
    # Compiles summary statistics for RMSE across gamma values.
    summary_stats <- data.frame(
        gamma = gammas,
        mean_rmse = rowMeans(cv_results),
        se_rmse = apply(cv_results, 1, function(x) sd(x) / sqrt(length(x))),
        median_rmse = apply(cv_results, 1, function(x) median(x)),
        q1_rmse = apply(cv_results, 1, function(x) quantile(x, 0.25)),
        q3_rmse = apply(cv_results, 1, function(x) quantile(x, 0.75)),
        iqr_rmse = apply(cv_results, 1, function(x) IQR(x)),
        p_value = p_values,
        q_value = q_values,
        test_used = ifelse(gammas == ref_gamma, "N/A", "Wilcoxon"),
        significance = plot_labels,
        stringsAsFactors = FALSE
    )
    
    # Tests normality of RMSE distributions for each gamma.
    is_normal <- numeric(length(gammas))
    for (i in seq_along(gammas)) {
        tryCatch({
            if (length(unique(cv_results[i, ])) > 1) {
                shapiro_test <- shapiro.test(cv_results[i, ])
                is_normal[i] <- shapiro_test$p.value > 0.05
            } else {
                is_normal[i] <- TRUE
            }
        }, error = function(e) {
            is_normal[i] <- TRUE
        })
    }
    
    # Adds normality results to the summary statistics.
    summary_stats$is_normal <- is_normal
    summary_stats$distribution <- ifelse(is_normal, "Normal", "Non-normal")
    
    # Reorders columns in the summary for improved readability.
    summary_stats <- summary_stats[, c("gamma", "is_normal", "distribution", "mean_rmse", "se_rmse", 
                                        "median_rmse", "q1_rmse", "q3_rmse", "iqr_rmse", 
                                        "p_value", "q_value", "test_used", "significance")]
    
    # Returns a list containing the boxplot, summary, and raw RMSE data.
    return(list(
        cv_plot = p,
        cv_summary = summary_stats,
        rmse_data = as.data.frame(cv_results)
    ))
}

# Executes the cross-validation evaluation with the specified dataset.
result <- evaluate_gamma_boxplot(dataframe_name)

# Displays the RMSE boxplot for gamma comparison.
print(result$cv_plot)

# Outputs the summary statistics for review.
print(result$cv_summary)