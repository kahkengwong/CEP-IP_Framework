#######################################################################################################################
# Part 3.11: Monte Carlo Cross-Validation (MCCV) Analysis of Decomposition of Deviance Explained (DDE) Classification
#######################################################################################################################

# =================================================================================================================
# Part I: DDE Classification MCCV Analysis with Random Negative Control Classification 
# =================================================================================================================
# Cross-validation function for deviance cells
validate_deviance_cells <- function(current_sample, train_prop = 0.7, n_iterations = 20) {
    cat("\n=== Cross-Validation for", current_sample, "===\n")
    
    # Define pre-optimized parameters for each sample
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Extract and filter sample data for analysis
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    original_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    original_dev_explained <- summary(original_model)$dev.expl
    
    if (nrow(sample_data) < 50) {
        cat("Insufficient data points for meaningful CV (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    validation_results <- list()
    
    for (iter in 1:n_iterations) {
        set.seed(2 + iter)  # Ensure reproducible train-test splits
        
        # Split data into training and testing sets
        n_total <- nrow(sample_data)
        train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
        test_idx <- setdiff(1:n_total, train_idx)
        
        train_data <- sample_data[train_idx, ]
        test_data <- sample_data[test_idx, ]
        
        tryCatch({
            # Fit GAM model on training data with fixed parameters
            train_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                               data = train_data, method = "REML", select = TRUE, 
                               gamma = 1.5, sp = optimal_lambda)
            
            # Calculate null model and residuals for training data
            null_model_train <- gam(Expression ~ 1, data = train_data)
            
            # Compute fitted values and residuals
            null_fitted_train <- fitted(null_model_train)
            model_fitted_train <- fitted(train_model)
            null_residuals_train <- train_data$Expression - null_fitted_train
            model_residuals_train <- train_data$Expression - model_fitted_train
            
            # Calculate squared residuals for deviance
            null_sq_diff_train <- null_residuals_train^2
            model_sq_diff_train <- model_residuals_train^2
            
            # Compute explanatory power for training data
            explanatory_power_train <- 1 - (model_sq_diff_train / pmax(null_sq_diff_train, 1e-8))
            
            # Sort cells by explanatory power
            sorted_indices_train <- order(explanatory_power_train, decreasing = TRUE)
            
            # Calculate deviance explained and select top cells
            dev_explained_train <- summary(train_model)$dev.expl
            target_cells_train <- round(nrow(train_data) * dev_explained_train)
            
            # Identify deviance and non-deviance cells
            deviance_cells_train_idx <- sorted_indices_train[1:target_cells_train]
            non_deviance_cells_train_idx <- sorted_indices_train[(target_cells_train+1):length(sorted_indices_train)]
            
            # Predict on test data
            test_pred <- predict(train_model, newdata = test_data)
            
            # Calculate explanatory power for test data
            null_model_test <- gam(Expression ~ 1, data = test_data)
            null_fitted_test <- fitted(null_model_test)
            
            test_null_residuals <- test_data$Expression - null_fitted_test
            test_model_residuals <- test_data$Expression - test_pred
            test_null_sq_diff <- test_null_residuals^2
            test_model_sq_diff <- test_model_residuals^2
            
            # Compute explanatory power for test cells
            test_explanatory_power <- 1 - (test_model_sq_diff / pmax(test_null_sq_diff, 1e-8))
            
            # Sort and classify test cells
            test_sorted_indices <- order(test_explanatory_power, decreasing = TRUE)
            test_target_cells <- round(nrow(test_data) * dev_explained_train)
            
            test_deviance_cells_idx <- test_sorted_indices[1:test_target_cells]
            test_non_deviance_cells_idx <- test_sorted_indices[(test_target_cells+1):length(test_sorted_indices)]
            
            # Calculate RMSE for deviance and non-deviance cells
            if (length(test_deviance_cells_idx) > 0 && length(test_non_deviance_cells_idx) > 0) {
                test_residuals <- test_data$Expression - test_pred
                test_dec_rmse <- sqrt(mean(test_residuals[test_deviance_cells_idx]^2, na.rm = TRUE))
                test_non_dec_rmse <- sqrt(mean(test_residuals[test_non_deviance_cells_idx]^2, na.rm = TRUE))
                
                validation_results[[iter]] <- data.frame(
                    iteration = iter,
                    sample = current_sample,
                    train_dev_explained = dev_explained_train,
                    train_n = nrow(train_data),
                    test_n = nrow(test_data),
                    test_dec_n = length(test_deviance_cells_idx),
                    test_non_dec_n = length(test_non_deviance_cells_idx),
                    test_dec_rmse = test_dec_rmse,
                    test_non_dec_rmse = test_non_dec_rmse,
                    rmse_difference = test_non_dec_rmse - test_dec_rmse,
                    used_k = optimal_k,
                    used_lambda = optimal_lambda,
                    train_edf = summary(train_model)$edf
                )
            }
            
        }, error = function(e) {
            cat("Error in iteration", iter, ":", conditionMessage(e), "\n")
        })
    }
    
    # Combine and summarize validation results
    if (length(validation_results) > 0) {
        combined_results <- do.call(rbind, validation_results)
        
        cat("\n--- Cross-Validation Summary ---\n")
        cat("Used k =", optimal_k, ", lambda =", optimal_lambda, "\n")
        cat("Successful iterations:", nrow(combined_results), "/", n_iterations, "\n")
        cat("Mean RMSE difference (Non-TDE - DEC):", round(mean(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("SD RMSE difference:", round(sd(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("Mean training EDF:", round(mean(combined_results$train_edf, na.rm = TRUE), 2), "\n")
        
        return(combined_results)
    } else {
        cat("No successful validation iterations\n")
        return(NULL)
    }
}

# Random control validation function
validate_random_deviance_cells <- function(current_sample, train_prop = 0.7, n_iterations = 20) {
    cat("\n=== Random Negative Control Cross-Validation for", current_sample, "===\n")
    
    # Define pre-optimized parameters for each sample
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Extract and filter sample data for analysis
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    original_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    original_dev_explained <- summary(original_model)$dev.expl
    
    if (nrow(sample_data) < 50) {
        cat("Insufficient data points for meaningful CV (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    validation_results <- list()
    
    for (iter in 1:n_iterations) {
        set.seed(2 + iter)  # Ensure reproducible train-test splits
        
        # Split data into training and testing sets
        n_total <- nrow(sample_data)
        train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
        test_idx <- setdiff(1:n_total, train_idx)
        
        train_data <- sample_data[train_idx, ]
        test_data <- sample_data[test_idx, ]
        
        tryCatch({
            # Fit GAM model on training data with fixed parameters
            train_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                               data = train_data, method = "REML", select = TRUE, 
                               gamma = 1.5, sp = optimal_lambda)
            
            # Calculate deviance explained and target cells
            dev_explained_train <- summary(train_model)$dev.expl
            target_cells_train <- round(nrow(train_data) * dev_explained_train)
            
            # Predict on test data
            test_pred <- predict(train_model, newdata = test_data)
            test_residuals <- test_data$Expression - test_pred
            
            # Randomly classify test cells based on training proportion
            test_target_cells <- round(nrow(test_data) * dev_explained_train)
            random_dec_idx <- sample(1:nrow(test_data), size = test_target_cells)
            random_non_dec_idx <- setdiff(1:nrow(test_data), random_dec_idx)
            
            # Calculate RMSE for random groups
            if (length(random_dec_idx) > 0 && length(random_non_dec_idx) > 0) {
                test_random_dec_rmse <- sqrt(mean(test_residuals[random_dec_idx]^2, na.rm = TRUE))
                test_random_non_dec_rmse <- sqrt(mean(test_residuals[random_non_dec_idx]^2, na.rm = TRUE))
                
                validation_results[[iter]] <- data.frame(
                    iteration = iter,
                    sample = current_sample,
                    train_dev_explained = dev_explained_train,
                    train_n = nrow(train_data),
                    test_n = nrow(test_data),
                    test_dec_n = length(random_dec_idx),
                    test_non_dec_n = length(random_non_dec_idx),
                    test_dec_rmse = test_random_dec_rmse,
                    test_non_dec_rmse = test_random_non_dec_rmse,
                    rmse_difference = test_random_non_dec_rmse - test_random_dec_rmse,
                    used_k = optimal_k,
                    used_lambda = optimal_lambda,
                    train_edf = summary(train_model)$edf
                )
            }
            
        }, error = function(e) {
            cat("Error in iteration", iter, ":", conditionMessage(e), "\n")
        })
    }
    
    # Combine and summarize random validation results
    if (length(validation_results) > 0) {
        combined_results <- do.call(rbind, validation_results)
        
        cat("\n--- Random Cross-Validation Summary ---\n")
        cat("Used k =", optimal_k, ", lambda =", optimal_lambda, "\n")
        cat("Successful iterations:", nrow(combined_results), "/", n_iterations, "\n")
        cat("Mean RMSE difference (Random Non-TDE - Random DEC):", round(mean(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("SD RMSE difference:", round(sd(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("Mean training EDF:", round(mean(combined_results$train_edf, na.rm = TRUE), 2), "\n")
        
        return(combined_results)
    } else {
        cat("No successful random validation iterations\n")
        return(NULL)
    }
}

# Statistical testing function for RMSE differences
perform_sample_statistical_test <- function(validation_results, sample_name) {
    if (is.null(validation_results) || nrow(validation_results) == 0) {
        return(NULL)
    }
    
    rmse_differences <- validation_results$rmse_difference
    rmse_differences <- rmse_differences[!is.na(rmse_differences)]
    
    if (length(rmse_differences) < 3) {
        cat("Insufficient data Statistical testing function for RMSE differences for statistical testing in", sample_name, "\n")
        return(NULL)
    }
    
    # Test normality of RMSE differences
    normality_test <- shapiro.test(rmse_differences)
    is_normal <- normality_test$p.value > 0.05
    
    # Perform t-test or Wilcoxon test based on normality
    if (is_normal) {
        stat_test <- t.test(rmse_differences, mu = 0)
        test_method <- "One-sample t-test"
        test_statistic <- stat_test$statistic
        test_pvalue <- stat_test$p.value
        conf_int <- stat_test$conf.int
    } else {
        stat_test <- wilcox.test(rmse_differences, mu = 0)
        test_method <- "Wilcoxon signed-rank test"
        test_statistic <- stat_test$statistic
        test_pvalue <- stat_test$p.value
        conf_int <- c(NA, NA)
    }
    
    # Calculate effect size
    if (is_normal) {
        effect_size <- mean(rmse_differences) / sd(rmse_differences)
        effect_size_name <- "Cohen's d"
    } else {
        effect_size <- abs(qnorm(test_pvalue/2)) / sqrt(length(rmse_differences))
        effect_size_name <- "Effect size r"
    }
    
    cat("\n--- Statistical Test Results for", sample_name, "---\n")
    cat("Shapiro-Wilk normality test p-value:", format.pval(normality_test$p.value, digits = 3), "\n")
    cat("Data distribution:", if(is_normal) "Normal" else "Non-normal", "\n")
    cat("Statistical test used:", test_method, "\n")
    cat("Test statistic:", round(as.numeric(test_statistic), 4), "\n")
    cat("P-value:", format.pval(test_pvalue, digits = 3), "\n")
    if (!is.na(conf_int[1])) {
        cat("95% CI: [", round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]\n")
    }
    cat(effect_size_name, ":", round(effect_size, 4), "\n")
    
    return(data.frame(
        sample = sample_name,
        n_iterations = length(rmse_differences),
        mean_rmse_diff = mean(rmse_differences),
        sd_rmse_diff = sd(rmse_differences),
        shapiro_pvalue = normality_test$p.value,
        is_normal = is_normal,
        test_method = test_method,
        test_statistic = as.numeric(test_statistic),
        test_pvalue = test_pvalue,
        conf_int_lower = conf_int[1],
        conf_int_upper = conf_int[2],
        effect_size = effect_size,
        effect_size_name = effect_size_name
    ))
}

# Visualization function for validation results
create_validation_plots <- function(validation_results, sample_name) {
    if (is.null(validation_results) || nrow(validation_results) == 0) {
        cat("No validation results to plot for", sample_name, "\n")
        return(NULL)
    }
    
    # Plot RMSE differences across iterations
    p1 <- ggplot(validation_results, aes(x = iteration, y = rmse_difference)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
        labs(title = paste("RMSE Difference Across CV Iterations (", sample_name, ")"),
             subtitle = "Positive values indicate better prediction for TDE cells",
             x = "CV Iteration", y = "RMSE Difference (Non-TDE - DEC)") +
        theme_minimal()
    
    # Plot distribution of RMSE differences
    p2 <- ggplot(validation_results, aes(x = rmse_difference)) +
        geom_histogram(bins = 10, fill = "#4B0082", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        geom_vline(xintercept = mean(validation_results$rmse_difference), linetype = "solid", color = "blue") +
        labs(title = paste("Distribution of RMSE Differences (", sample_name, ")"),
             subtitle = "Red line at 0, blue line at mean",
             x = "RMSE Difference (Non-TDE - DEC)", y = "Frequency") +
        theme_minimal()
    
    print(p1)
    print(p2)
    
    return(list(iteration_plot = p1, distribution_plot = p2))
}

# Execute validation for all tumor samples
cat("\n", rep("=", 60), "\n")
cat("GAM-ALIGNED VALIDATION WITH STATISTICAL TESTING AND RANDOM CONTROL\n")
cat(rep("=", 60), "\n")

all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                       "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                       "HYW_5755_Tumor")
all_validation_results <- list()
all_random_validation_results <- list()
all_statistical_tests <- list()
all_random_statistical_tests <- list()

# Validate DDE and random control classifications for all tumor samples
for (sample_name in all_tumor_samples) {
    if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
        cat("\n", rep("=", 60), "\n")
        # Run DDE validation
        validation_result <- validate_deviance_cells(sample_name, train_prop = 0.7, n_iterations = 20)
        
        if (!is.null(validation_result)) {
            all_validation_results[[sample_name]] <- validation_result
            
            # Perform statistical test for validation results
            stat_test_result <- perform_sample_statistical_test(validation_result, sample_name)
            if (!is.null(stat_test_result)) {
                all_statistical_tests[[sample_name]] <- stat_test_result
            }
        }
        
        # Run random negative control validation
        random_validation_result <- validate_random_deviance_cells(sample_name, train_prop = 0.7, n_iterations = 20)
        
        if (!is.null(random_validation_result)) {
            all_random_validation_results[[sample_name]] <- random_validation_result
            
            # Perform statistical test for random control
            random_stat_test_result <- perform_sample_statistical_test(random_validation_result, paste0(sample_name, "_Random"))
            if (!is.null(random_stat_test_result)) {
                all_random_statistical_tests[[sample_name]] <- random_stat_test_result
            }
        }
    } else {
        cat("\n--- No results found for", sample_name, "---\n")
    }
}

# Summarize statistical tests with FDR correction
if (length(all_statistical_tests) > 0 || length(all_random_statistical_tests) > 0) {
    # Combine statistical test results
    combined_tests <- do.call(rbind, all_statistical_tests)
    if (!is.null(combined_tests)) {
        # Apply FDR correction to p-values
        combined_tests$fdr_corrected_pvalue <- p.adjust(combined_tests$test_pvalue, method = "BH")
        
        # Map sample names to patient IDs
        sample_to_patient <- c(
            "HYW_4701_Tumor" = "Pt.1",
            "HYW_4847_Tumor" = "Pt.2", 
            "HYW_4880_Tumor" = "Pt.3",
            "HYW_4881_Tumor" = "Pt.4",
            "HYW_5386_Tumor" = "Pt.5",
            "HYW_5742_Tumor" = "Pt.6",
            "HYW_5755_Tumor" = "Pt.7"
        )
        combined_tests$patient_id <- sample_to_patient[combined_tests$sample]
    } else {
        combined_tests <- data.frame()
    }
    
    # Combine random control statistical test results
    combined_random_tests <- do.call(rbind, all_random_statistical_tests)
    if (!is.null(combined_random_tests)) {
        combined_random_tests$fdr_corrected_pvalue <- p.adjust(combined_random_tests$test_pvalue, method = "BH")
        combined_random_tests$patient_id <- sample_to_patient[gsub("_Random$", "", combined_random_tests$sample)]
    } else {
        combined_random_tests <- data.frame()
    }
    
    cat("\n", rep("=", 60), "\n")
    cat("OVERALL VALIDATION SUMMARY WITH STATISTICAL TESTING\n")
    cat(rep("=", 60), "\n")
    
    # Display results for DDE classification
    if (nrow(combined_tests) > 0) {
        cat("\n--- DDE Classification Results ---\n")
        print(combined_tests[, c("patient_id", "mean_rmse_diff", "test_method", "test_pvalue", "fdr_corrected_pvalue", "effect_size")])
    }
    
    # Display results for random control classification
    if (nrow(combined_random_tests) > 0) {
        cat("\n--- Random Negative Control Classification Results ---\n")
        print(combined_random_tests[, c("patient_id", "mean_rmse_diff", "test_method", "test_pvalue", "fdr_corrected_pvalue", "effect_size")])
    }
    
    # Summarize statistics for DDE validation
    if (nrow(combined_tests) > 0) {
        cat("\n--- Overall Summary (EP) ---\n")
        cat("Number of samples tested:", nrow(combined_tests), "\n")
        cat("Samples with normal distribution:", sum(combined_tests$is_normal), "\n")
        cat("Samples with non-normal distribution:", sum(!combined_tests$is_normal), "\n")
        cat("Mean RMSE difference across all samples:", round(mean(combined_tests$mean_rmse_diff), 4), "\n")
        cat("Samples with significant results (p < 0.05):", sum(combined_tests$test_pvalue < 0.05), "\n")
        cat("Samples with significant results (FDR < 0.05):", sum(combined_tests$fdr_corrected_pvalue < 0.05), "\n")
    }
    
    # Summarize statistics for random control validation
    if (nrow(combined_random_tests) > 0) {
        cat("\n--- Overall Summary (Random Negative Control) ---\n")
        cat("Number of samples tested:", nrow(combined_random_tests), "\n")
        cat("Samples with normal distribution:", sum(combined_random_tests$is_normal), "\n")
        cat("Samples with non-normal distribution:", sum(!combined_random_tests$is_normal), "\n")
        cat("Mean RMSE difference across all samples:", round(mean(combined_random_tests$mean_rmse_diff), 4), "\n")
        cat("Samples with significant results (p < 0.05):", sum(combined_random_tests$test_pvalue < 0.05), "\n")
        cat("Samples with significant results (FDR < 0.05):", sum(combined_random_tests$fdr_corrected_pvalue < 0.05), "\n")
    }
    
    # Create visualization for validation results
    if (length(all_validation_results) > 0 || length(all_random_validation_results) > 0) {
        # Combine DDE validation results
        combined_all <- do.call(rbind, all_validation_results)
        if (!is.null(combined_all)) {
            combined_all$patient_id <- sample_to_patient[combined_all$sample]
        }
        
        # Combine random validation results
        combined_random_all <- do.call(rbind, all_random_validation_results)
        if (!is.null(combined_random_all)) {
            combined_random_all$patient_id <- sample_to_patient[combined_random_all$sample]
        }
        
        # Prepare data for DDE boxplot
        validation_plot_data <- data.frame()
        for (sample_name in names(all_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data <- all_validation_results[[sample_name]]
            for (i in 1:nrow(sample_data)) {
                validation_plot_data <- rbind(validation_plot_data, 
                                              data.frame(
                                                  Sample = patient_id,
                                                  RMSE = sample_data$test_dec_rmse[i],
                                                  Group = "DEC"
                                              ))
                validation_plot_data <- rbind(validation_plot_data, 
                                              data.frame(
                                                  Sample = patient_id,
                                                  RMSE = sample_data$test_non_dec_rmse[i],
                                                  Group = "Non-TDE"
                                              ))
            }
        }
        
        # Prepare data for random control boxplot
        random_validation_plot_data <- data.frame()
        for (sample_name in names(all_random_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data <- all_random_validation_results[[sample_name]]
            for (i in 1:nrow(sample_data)) {
                random_validation_plot_data <- rbind(random_validation_plot_data, 
                                                     data.frame(
                                                         Sample = patient_id,
                                                         RMSE = sample_data$test_dec_rmse[i],
                                                         Group = "Random DEC"
                                                     ))
                random_validation_plot_data <- rbind(random_validation_plot_data, 
                                                     data.frame(
                                                         Sample = patient_id,
                                                         RMSE = sample_data$test_non_dec_rmse[i],
                                                         Group = "Random Non-TDE"
                                                     ))
            }
        }
        
        # Add significance annotations for DDE results
        if (nrow(combined_tests) > 0) {
            sig_values <- ifelse(combined_tests$fdr_corrected_pvalue < 0.001, "***",
                                 ifelse(combined_tests$fdr_corrected_pvalue < 0.01, "**",
                                        ifelse(combined_tests$fdr_corrected_pvalue < 0.05, "ns", "ns")))
        } else {
            sig_values <- character(0)
        }
        
        # Add significance annotations for random control results
        if (nrow(combined_random_tests) > 0) {
            random_sig_values <- ifelse(combined_random_tests$fdr_corrected_pvalue < 0.001, "***",
                                        ifelse(combined_random_tests$fdr_corrected_pvalue < 0.01, "**",
                                               ifelse(combined_random_tests$fdr_corrected_pvalue < 0.05, "ns", "ns")))
        } else {
            random_sig_values <- character(0)
        }
        
        # Perform overall statistical test for DDE results
        if (!is.null(combined_all) && nrow(combined_all) > 0) {
            all_rmse_diffs <- combined_all$rmse_difference[!is.na(combined_all$rmse_difference)]
            overall_normality <- shapiro.test(all_rmse_diffs)
            
            if (overall_normality$p.value > 0.05) {
                overall_test <- t.test(all_rmse_diffs, mu = 0)
                overall_method <- "One-sample t-test"
            } else {
                overall_test <- wilcox.test(all_rmse_diffs, mu = 0)
                overall_method <- "Wilcoxon signed-rank test"
            }
            
            cat("\n--- Overall Test Across All Samples (EP) ---\n")
            cat("Overall normality test p-value:", format.pval(overall_normality$p.value, digits = 3), "\n")
            cat("Overall test method:", overall_method, "\n")
            cat("Overall test p-value:", format.pval(overall_test$p.value, digits = 3), "\n")
        }
        
        # Perform overall statistical test for random control results
        if (!is.null(combined_random_all) && nrow(combined_random_all) > 0) {
            all_random_rmse_diffs <- combined_random_all$rmse_difference[!is.na(combined_random_all$rmse_difference)]
            overall_random_normality <- shapiro.test(all_random_rmse_diffs)
            
            if (overall_random_normality$p.value > 0.05) {
                overall_random_test <- t.test(all_random_rmse_diffs, mu = 0)
                overall_random_method <- "One-sample t-test"
            } else {
                overall_random_test <- wilcox.test(all_random_rmse_diffs, mu = 0)
                overall_random_method <- "Wilcoxon signed-rank test"
            }
            
            cat("\n--- Overall Test Across All Samples (Random Negative Control) ---\n")
            cat("Overall normality test p-value:", format.pval(overall_random_normality$p.value, digits = 3), "\n")
            cat("Overall test method:", overall_random_method, "\n")
            cat("Overall test p-value:", format.pval(overall_random_test$p.value, digits = 3), "\n")
        }
        
        # Create boxplot for DDE classification
        if (nrow(validation_plot_data) > 0) {
            validation_plot_data$Group <- factor(validation_plot_data$Group, levels = c("Non-TDE", "DEC"))
            
            # Set dynamic y-axis limits
            y_max <- max(validation_plot_data$RMSE, na.rm = TRUE) + 0.1
            y_min <- min(validation_plot_data$RMSE, na.rm = TRUE) - 0.05
            
            p1 <- ggplot(validation_plot_data, aes(x = Sample, y = RMSE)) +
                geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.shape = NA, width = 0.6,
                             position = position_dodge(width = 0.8)) +
                geom_point(aes(color = Group), alpha = 0.5, size = 1.2,
                           position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15)) +
                scale_fill_manual(values = c("Non-TDE" = "#C0C0C0", "DEC" = "#9954C1")) +
                scale_color_manual(values = c("Non-TDE" = "#A0A0A0", "DEC" = "#4B0082")) +
                geom_text(data = data.frame(Sample = combined_tests$patient_id, 
                                            RMSE = rep(y_max - 0.02, nrow(combined_tests)),
                                            sig = sig_values),
                          aes(x = Sample, y = RMSE, label = sig), 
                          inherit.aes = FALSE, size = 4, fontface = "bold") +
                labs(title = "RMSE Classification Validity",
                     subtitle = "TDE vs Non-TDE cells | Tests if classification captures low RMSE cells | ** = FDR < 0.01",
                     x = "Sample", y = "RMSE", fill = "Classification", color = "Classification") +
                theme_minimal() +
                theme(
                    plot.background = element_rect(fill = "white", color = "#2E2E2E", size = 0.3),
                    panel.background = element_rect(fill = "white"),
                    panel.border = element_rect(color = "#2E2E2E", fill = NA, size = 0.3),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold"),
                    legend.position = "bottom"
                ) +
                scale_y_continuous(limits = c(y_min, y_max)) +
                guides(color = "none")
            
            print(p1)
        }
        
        # Create boxplot for random control classification
        if (nrow(random_validation_plot_data) > 0) {
            random_validation_plot_data$Group <- factor(random_validation_plot_data$Group, levels = c("Random Non-TDE", "Random DEC"))
            
            # Set dynamic y-axis limits for random plot
            y_max_random <- max(random_validation_plot_data$RMSE, na.rm = TRUE) + 0.1
            y_min_random <- min(random_validation_plot_data$RMSE, na.rm = TRUE) - 0.05
            
            p1_random <- ggplot(random_validation_plot_data, aes(x = Sample, y = RMSE)) +
                geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.shape = NA, width = 0.6,
                             position = position_dodge(width = 0.8)) +
                geom_point(aes(color = Group), alpha = 0.5, size = 1.2,
                           position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15)) +
                scale_fill_manual(values = c("Random Non-TDE" = "#C0C0C0", "Random DEC" = "#9954C1")) +
                scale_color_manual(values = c("Random Non-TDE" = "#A0A0A0", "Random DEC" = "#4B0082")) +
                geom_text(data = data.frame(Sample = combined_random_tests$patient_id, 
                                            RMSE = rep(y_max_random - 0.02, nrow(combined_random_tests)),
                                            sig = random_sig_values),
                          aes(x = Sample, y = RMSE, label = sig), 
                          inherit.aes = FALSE, size = 4, fontface = "bold") +
                labs(title = "Random Control Classification",
                     subtitle = "Random TDE vs Random Non-TDE cells | Should show no significant difference",
                     x = "Sample", y = "RMSE", fill = "Classification", color = "Classification") +
                theme_minimal() +
                theme(
                    plot.background = element_rect(fill = "white", color = "#2E2E2E", size = 0.3),
                    panel.background = element_rect(fill = "white"),
                    panel.border = element_rect(color = "#2E2E2E", fill = NA, size = 0.3),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold"),
                    legend.position = "bottom"
                ) +
                scale_y_continuous(limits = c(y_min_random, y_max_random)) +
                guides(color = "none")
            
            print(p1_random)
        }
        
        # Create effect size plot for DDE results
        if (nrow(combined_tests) > 0) {
            p2 <- ggplot(combined_tests, aes(x = reorder(patient_id, effect_size), y = effect_size, fill = fdr_corrected_pvalue < 0.01)) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkgreen"),
                                  name = "FDR < 0.05", labels = c("No", "Yes")) +
                labs(title = "Effect Sizes by Sample (EP)",
                     subtitle = "Green bars indicate FDR < 0.01",
                     x = "Sample", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            print(p2)
        }
        
        # Create effect size plot for random control results
        if (nrow(combined_random_tests) > 0) {
            p3 <- ggplot(combined_random_tests, aes(x = reorder(patient_id, effect_size), y = effect_size, fill = fdr_corrected_pvalue < 0.01)) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkgreen"),
                                  name = "FDR < 0.05", labels = c("No", "Yes")) +
                labs(title = "Effect Sizes by Sample (Random Negative Control)",
                     subtitle = "Green bars indicate FDR < 0.01",
                     x = "Sample", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            print(p3)
        }
        
        # Export results to Excel
        library(openxlsx)
        wb <- createWorkbook()
        
        # Export DDE sample results
        for (sample_name in names(all_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data_with_id <- all_validation_results[[sample_name]]
            sample_data_with_id$patient_id <- patient_id
            addWorksheet(wb, paste0(patient_id, "_TRPM4"))
            writeData(wb, paste0(patient_id, "_TRPM4"), sample_data_with_id)
        }
        
        # Export random control sample results
        for (sample_name in names(all_random_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data_with_id <- all_random_validation_results[[sample_name]]
            sample_data_with_id$patient_id <- patient_id
            addWorksheet(wb, paste0(patient_id, "_Random"))
            writeData(wb, paste0(patient_id, "_Random"), sample_data_with_id)
        }
        
        # Export DDE statistical test results
        if (nrow(combined_tests) > 0) {
            addWorksheet(wb, "Statistical_Tests_TRPM4")
            writeData(wb, "Statistical_Tests_TRPM4", combined_tests)
        }
        
        # Export random control statistical test results
        if (nrow(combined_random_tests) > 0) {
            addWorksheet(wb, "Statistical_Tests_Random")
            writeData(wb, "Statistical_Tests_Random", combined_random_tests)
        }
        
        # Export combined DDE results
        if (!is.null(combined_all)) {
            addWorksheet(wb, "Combined_Results_TRPM4")
            writeData(wb, "Combined_Results_TRPM4", combined_all)
        }
        
        # Export combined random control results
        if (!is.null(combined_random_all)) {
            addWorksheet(wb, "Combined_Results_Random")
            writeData(wb, "Combined_Results_Random", combined_random_all)
        }
        
        # Save Excel workbook
        saveWorkbook(wb, "[1]_EP_CV.xlsx", overwrite = TRUE)
    }
}
cat("\n", rep("=", 60), "\n")
cat("GAM-ALIGNED VALIDATION WITH STATISTICAL TESTING AND RANDOM CONTROL COMPLETE\n")
cat(rep("=", 60), "\n")

# Create GAM plots for DDE classification
create_gam_validation_plots <- function(current_sample) {
    cat("\n=== Creating GAM Visualization for", current_sample, "===\n")
    
    # Define pre-optimized parameters for each sample
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Extract and filter sample data
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    if (nrow(sample_data) < 10) {
        cat("Insufficient data points for plotting (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    # Fit GAM model with pre-optimized parameters
    gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                     data = sample_data, method = "REML", select = TRUE, 
                     gamma = 1.5, sp = optimal_lambda)
    
    # Calculate null model and residuals
    null_model <- gam(Expression ~ 1, data = sample_data)
    
    # Compute fitted values and residuals
    null_fitted <- fitted(null_model)
    model_fitted <- fitted(gam_model)
    null_residuals <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted
    
    # Calculate squared residuals for deviance
    null_sq_diff <- null_residuals^2
    model_sq_diff <- model_residuals^2
    
    # Compute explanatory power
    explanatory_power <- 1 - (model_sq_diff / pmax(null_sq_diff, 1e-8))
    
    # Sort cells by explanatory power
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    
    # Calculate deviance explained and select top cells
    model_dev_explained <- summary(gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    
    # Identify deviance and non-deviance cells
    deviance_cells_idx <- sorted_indices[1:target_cells]
    non_deviance_cells_idx <- sorted_indices[(target_cells+1):length(sorted_indices)]
    
    # Set dynamic axis limits
    x_range <- range(sample_data$TRPM4, na.rm = TRUE)
    y_range <- range(sample_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)
    
    # Prepare data for plotting
    plot_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(1:nrow(sample_data) %in% deviance_cells_idx, "Dev explained", "Non-dev explained"),
        ExplanatoryPower = explanatory_power
    )
    
    # Set drawing order for points
    plot_data$draw_order <- ifelse(plot_data$Group == "Dev explained", 2, 1)
    
    # Sort by draw order
    plot_data <- plot_data[order(plot_data$draw_order), ]
    
    # Create prediction data for GAM line
    pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
    pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit <- pred$fit
    pred_data$se.fit <- pred$se.fit
    
    # Calculate deviance explained percentage
    dev_explained_percentage <- round(model_dev_explained * 100, 2)
    
    # Map sample names to patient IDs
    sample_to_patient <- c(
        "HYW_4701_Tumor" = "Pt.1",
        "HYW_4847_Tumor" = "Pt.2", 
        "HYW_4880_Tumor" = "Pt.3",
        "HYW_4881_Tumor" = "Pt.4",
        "HYW_5386_Tumor" = "Pt.5",
        "HYW_5742_Tumor" = "Pt.6",
        "HYW_5755_Tumor" = "Pt.7"
    )
    
    patient_id <- sample_to_patient[current_sample]
    
    # Create scatter plot with GAM fit
    p <- ggplot() +
        geom_point(data = plot_data, 
                   aes(x = TRPM4, y = Expression, color = Group, alpha = Group),
                   size = 1.8) +
        scale_alpha_manual(values = c("Dev explained" = 0.4, "Non-dev explained" = 0.6), guide = "none") +
        geom_line(data = pred_data, aes(x = TRPM4, y = fit), color = "#FFCC99", size = 1.2) +
        geom_ribbon(data = pred_data, aes(x = TRPM4, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), 
                    fill = "#FFCC99", alpha = 0.2) +
        scale_color_manual(values = c("Dev explained" = "#4B0082", "Non-dev explained" = "#C0C0C0")) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position = "bottom",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.line = element_blank()
        ) +
        labs(
            title = paste("TRPM4 vs Ribo Expression (", patient_id, ")", sep = ""),
            subtitle = paste("Dev explained: ", dev_explained_percentage, "% | k=", optimal_k, ", EDF=", round(summary(gam_model)$edf, 2), " | Top ", target_cells, " cells", sep = ""),
            x = "TRPM4 Expression",
            y = "Ribo Expression",
            color = "Classification"
        ) +
        scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
    
    # Display the plot
    print(p)
    
    # Print classification summary
    cat("\n--- Classification Summary ---\n")
    cat("Total cells:", nrow(sample_data), "\n")
    cat("Dev explained cells:", length(deviance_cells_idx), "\n")
    cat("Non-dev explained cells:", length(non_deviance_cells_idx), "\n")
    cat("Classification completeness:", (length(deviance_cells_idx) + length(non_deviance_cells_idx)) == nrow(sample_data), "\n")
    cat("Model deviance explained:", round(model_dev_explained * 100, 2), "%\n")
    cat("Target cells (top X%):", target_cells, "out of", nrow(sample_data), "\n")
    cat("Effective degrees of freedom:", round(summary(gam_model)$edf, 2), "\n")
    
    # Summarize explanatory power
    cat("\n--- Explanatory Power Summary ---\n")
    cat("Dev explained cells - Min EP:", round(min(explanatory_power[deviance_cells_idx]), 4), 
        "| Max EP:", round(max(explanatory_power[deviance_cells_idx]), 4), "\n")
    cat("Non-dev explained cells - Min EP:", round(min(explanatory_power[non_deviance_cells_idx]), 4), 
        "| Max EP:", round(max(explanatory_power[non_deviance_cells_idx]), 4), "\n")
    
    return(list(plot = p, model = gam_model, classification_data = plot_data))
}

# Generate GAM plots for all samples
create_all_gam_validation_plots <- function() {
    cat("\n", rep("=", 60), "\n")
    cat("CREATING GAM VALIDATION PLOTS FOR ALL SAMPLES\n")
    cat(rep("=", 60), "\n")
    
    all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                           "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                           "HYW_5755_Tumor")
    
    all_plots <- list()
    
    for (sample_name in all_tumor_samples) {
        if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
            plot_result <- create_gam_validation_plots(sample_name)
            if (!is.null(plot_result)) {
                all_plots[[sample_name]] <- plot_result
            }
        } else {
            cat("\n--- No results found for", sample_name, "---\n")
        }
    }
    
    return(all_plots)
}

# Run the plotting function for all samples
all_gam_plots <- create_all_gam_validation_plots()

cat("\n", rep("=", 60), "\n")
cat("GAM VALIDATION PLOTS COMPLETE\n")
cat(rep("=", 60), "\n")

# Create GAM plots for random control classification
create_random_gam_validation_plots <- function(current_sample, random_seed = 42) {
    cat("\n=== Creating Random Control GAM Visualization for", current_sample, "===\n")
    
    # Define pre-optimized parameters for each sample
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Extract and filter sample data
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    if (nrow(sample_data) < 10) {
        cat("Insufficient data points for plotting (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    # Fit GAM model with pre-optimized parameters
    gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                     data = sample_data, method = "REML", select = TRUE, 
                     gamma = 1.5, sp = optimal_lambda)
    
    # Calculate deviance explained and target cells
    model_dev_explained <- summary(gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    
    # Perform random classification
    set.seed(random_seed)
    random_deviance_cells_idx <- sample(1:nrow(sample_data), size = target_cells)
    random_non_deviance_cells_idx <- setdiff(1:nrow(sample_data), random_deviance_cells_idx)
    
    # Calculate explanatory power for verification
    null_model <- gam(Expression ~ 1, data = sample_data)
    null_fitted <- fitted(null_model)
    model_fitted <- fitted(gam_model)
    null_residuals <- sample_data$Expression - null_fitted
    model_residuals <- sample_data$Expression - model_fitted
    null_sq_diff <- null_residuals^2
    model_sq_diff <- model_residuals^2
    explanatory_power <- 1 - (model_sq_diff / pmax(null_sq_diff, 1e-8))
    
    # Set dynamic axis limits
    x_range <- range(sample_data$TRPM4, na.rm = TRUE)
    y_range <- range(sample_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)
    
    # Prepare data for plotting
    plot_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(1:nrow(sample_data) %in% random_deviance_cells_idx, "Random DEC", "Random Non-TDE"),
        ExplanatoryPower = explanatory_power
    )
    
    # Set drawing order for points
    plot_data$draw_order <- ifelse(plot_data$Group == "Random DEC", 2, 1)
    
    # Sort by draw order
    plot_data <- plot_data[order(plot_data$draw_order), ]
    
    # Create prediction data for GAM line
    pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
    pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit <- pred$fit
    pred_data$se.fit <- pred$se.fit
    
    # Calculate deviance explained percentage
    dev_explained_percentage <- round(model_dev_explained * 100, 2)
    
    # Map sample names to patient IDs
    sample_to_patient <- c(
        "HYW_4701_Tumor" = "Pt.1",
        "HYW_4847_Tumor" = "Pt.2", 
        "HYW_4880_Tumor" = "Pt.3",
        "HYW_4881_Tumor" = "Pt.4",
        "HYW_5386_Tumor" = "Pt.5",
        "HYW_5742_Tumor" = "Pt.6",
        "HYW_5755_Tumor" = "Pt.7"
    )
    
    patient_id <- sample_to_patient[current_sample]
    
    # Create scatter plot with GAM fit for random classification
    p <- ggplot() +
        geom_point(data = plot_data, 
                   aes(x = TRPM4, y = Expression, color = Group, alpha = Group),
                   size = 1.8) +
        scale_alpha_manual(values = c("Random DEC" = 0.25, "Random Non-TDE" = 0.6), guide = "none") +
        geom_line(data = pred_data, aes(x = TRPM4, y = fit), color = "#FFCC99", size = 1.2) +
        geom_ribbon(data = pred_data, aes(x = TRPM4, ymin = fit - 1.96 * se.fit, ymax = fit + 1.96 * se.fit), 
                    fill = "#FFCC99", alpha = 0.2) +
        scale_color_manual(values = c("Random DEC" = "#4B0082", "Random Non-TDE" = "#C0C0C0")) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position = "bottom",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.line = element_blank()
        ) +
        labs(
            title = paste("Random Control: TRPM4 vs Ribo (", patient_id, ")", sep = ""),
            subtitle = paste("Random classification | Same model: ", dev_explained_percentage, "% | k=", optimal_k, ", EDF=", round(summary(gam_model)$edf, 2), " | Random ", target_cells, " cells", sep = ""),
            x = "TRPM4 Expression",
            y = "Ribo Expression",
            color = "Random Classification"
        ) +
        scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
    
    # Display the plot
    print(p)
    
    # Print random classification summary
    cat("\n--- Random Classification Summary ---\n")
    cat("Total cells:", nrow(sample_data), "\n")
    cat("Random TDE cells:", length(random_deviance_cells_idx), "\n")
    cat("Random Non-TDE cells:", length(random_non_deviance_cells_idx), "\n")
    cat("Classification completeness:", (length(random_deviance_cells_idx) + length(random_non_deviance_cells_idx)) == nrow(sample_data), "\n")
    cat("Model deviance explained:", round(model_dev_explained * 100, 2), "%\n")
    cat("Random target cells:", target_cells, "out of", nrow(sample_data), "\n")
    cat("Effective degrees of freedom:", round(summary(gam_model)$edf, 2), "\n")
    
    # Summarize explanatory power (EP) for random groups
    cat("\n--- Explanatory Power Summary (Random Groups) ---\n")
    cat("Random TDE cells - Min EP:", round(min(explanatory_power[random_deviance_cells_idx]), 4), 
        "| Max EP:", round(max(explanatory_power[random_deviance_cells_idx]), 4),
        "| Mean EP:", round(mean(explanatory_power[random_deviance_cells_idx]), 4), "\n")
    cat("Random Non-TDE cells - Min EP:", round(min(explanatory_power[random_non_deviance_cells_idx]), 4), 
        "| Max EP:", round(max(explanatory_power[random_non_deviance_cells_idx]), 4),
        "| Mean EP:", round(mean(explanatory_power[random_non_deviance_cells_idx]), 4), "\n")
    cat("EP difference between groups:", round(mean(explanatory_power[random_deviance_cells_idx]) - mean(explanatory_power[random_non_deviance_cells_idx]), 4), "\n")
    
    return(list(plot = p, model = gam_model, classification_data = plot_data))
}

# Generate random control GAM plots for all samples
create_all_random_gam_validation_plots <- function(random_seed = 42) {
    cat("\n", rep("=", 60), "\n")
    cat("CREATING RANDOM CONTROL GAM VALIDATION PLOTS FOR ALL SAMPLES\n")
    cat(rep("=", 60), "\n")
    
    all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                           "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                           "HYW_5755_Tumor")
    
    all_random_plots <- list()
    
    for (sample_name in all_tumor_samples) {
        if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
            random_plot_result <- create_random_gam_validation_plots(sample_name, random_seed = random_seed)
            if (!is.null(random_plot_result)) {
                all_random_plots[[sample_name]] <- random_plot_result
            }
        } else {
            cat("\n--- No results found for", sample_name, "---\n")
        }
    }
    
    return(all_random_plots)
}

# Run the random control plotting function
all_random_gam_plots <- create_all_random_gam_validation_plots(random_seed = 42)

cat("\n", rep("=", 60), "\n")
cat("RANDOM CONTROL GAM VALIDATION PLOTS COMPLETE\n")
cat(rep("=", 60), "\n")



# ====================================================================================
# Part II: Leverage-based Classification Monte-Carlo CV Analysis with Random Negative Control Classification 
# ====================================================================================
# Function to calculate leverage-based cell influence
calculate_cell_influence <- function(model, data) {
    tryCatch({
        # Get prediction matrix (design matrix)
        X <- predict(model, type = "lpmatrix")
        
        # Calculate leverage (hat values)
        # For large datasets, use approximate method to avoid memory issues
        if (nrow(X) > 1000) {
            # Approximate leverage using row sums method
            XtX_inv_diag <- 1 / colSums(X^2)  # Approximate diagonal of (X'X)^-1
            leverage <- rowSums(X^2 * rep(XtX_inv_diag, each = nrow(X)))
        } else {
            # Exact calculation for smaller datasets
            XtX <- crossprod(X)  # X'X
            XtX_inv <- solve(XtX + diag(1e-8, ncol(X)))  # Add small ridge for stability
            leverage <- rowSums((X %*% XtX_inv) * X)
        }
        
        # Ensure leverage values are in valid range [0,1]
        leverage <- pmax(0, pmin(1, leverage))
        
        # Get residuals and calculate standardized residuals
        residuals <- residuals(model)
        sigma <- sqrt(summary(model)$dispersion)
        std_residuals <- residuals / (sigma * sqrt(pmax(1 - leverage, 0.01)))
        
        # Calculate Cook's distance
        p <- ncol(X)  # number of parameters
        cooks_d <- (std_residuals^2 * leverage) / (p * pmax(1 - leverage, 0.01))
        
        # Calculate overall influence score (combines leverage and residual magnitude)
        influence_score <- sqrt(leverage * std_residuals^2)
        
        return(data.frame(
            leverage = leverage,
            std_residuals = std_residuals,
            cooks_distance = cooks_d,
            influence_score = influence_score
        ))
        
    }, error = function(e) {
        cat("Error calculating cell influence:", conditionMessage(e), "\n")
        # Return data frame with NA values if calculation fails
        return(data.frame(
            leverage = rep(NA, nrow(data)),
            std_residuals = rep(NA, nrow(data)),
            cooks_distance = rep(NA, nrow(data)),
            influence_score = rep(NA, nrow(data))
        ))
    })
}

# Leverage-based validation with simplified classification method
validate_deviance_cells <- function(current_sample, train_prop = 0.7, n_iterations = 20) {
    cat("\n=== Cross-Validation for", current_sample, "===\n")
    
    # Pre-optimized parameters from analyze_sample results
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Get sample data using exact same method as analyze_sample
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    original_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    original_dev_explained <- summary(original_model)$dev.expl
    
    if (nrow(sample_data) < 50) {
        cat("Insufficient data points for meaningful CV (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    validation_results <- list()
    
    for (iter in 1:n_iterations) {
        set.seed(6 + iter)  # Reproducible splits
        
        # Create train-test split
        n_total <- nrow(sample_data)
        train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
        test_idx <- setdiff(1:n_total, train_idx)
        
        train_data <- sample_data[train_idx, ]
        test_data <- sample_data[test_idx, ]
        
        tryCatch({
            # Fit GAM on training data with fixed k and lambda
            train_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                               data = train_data, method = "REML", select = TRUE, 
                               gamma = 1.5, sp = optimal_lambda)
            
            # === LEVERAGE-BASED RANKING ON TRAINING DATA ===
            
            # Calculate leverage-based influence metrics for training data
            cat("Calculating leverage-based influence metrics for iteration", iter, "...\n")
            influence_metrics_train <- calculate_cell_influence(train_model, train_data)
            
            # Check if influence calculation was successful, use fallback if needed
            if (all(is.na(influence_metrics_train$influence_score))) {
                cat("Warning: Influence calculation failed for iteration", iter, ", falling back to explanatory power method\n")
                
                # Fallback to explanatory power method
                null_model_train <- gam(Expression ~ 1, data = train_data)
                null_fitted_train <- fitted(null_model_train)
                model_fitted_train <- fitted(train_model)
                null_residuals_train <- train_data$Expression - null_fitted_train
                model_residuals_train <- train_data$Expression - model_fitted_train
                null_sq_diff_train <- null_residuals_train^2
                model_sq_diff_train <- model_residuals_train^2
                ranking_metric_train <- 1 - (model_sq_diff_train / pmax(null_sq_diff_train, 1e-8))
            } else {
                # Use influence score for ranking
                cat("Using influence score for cell ranking in iteration", iter, "\n")
                ranking_metric_train <- influence_metrics_train$influence_score
            }
            
            # === SIMPLIFIED CLASSIFICATION METHOD ===
            
            # Sort cells by their ranking metric (highest to lowest)
            sorted_indices_train <- order(ranking_metric_train, decreasing = TRUE, na.last = TRUE)
            
            # Get model deviance explained and target number of cells
            dev_explained_train <- summary(train_model)$dev.expl
            target_cells_train <- round(nrow(train_data) * dev_explained_train)
            target_cells_train <- max(1, min(target_cells_train, nrow(train_data) - 1))  # Ensure valid range
            
            # Take the top cells by ranking metric (TDE cells)
            deviance_cells_train_idx <- sorted_indices_train[1:target_cells_train]
            non_deviance_cells_train_idx <- sorted_indices_train[(target_cells_train+1):length(sorted_indices_train)]
            
            # === APPLY SAME CLASSIFICATION TO TEST DATA ===
            
            # Predict on test data
            test_pred <- predict(train_model, newdata = test_data)
            
            # Calculate ranking metric for test data using same method as training
            influence_metrics_test <- calculate_cell_influence(train_model, test_data)
            
            if (all(is.na(influence_metrics_test$influence_score))) {
                # Fallback to explanatory power method for test data
                null_model_test <- gam(Expression ~ 1, data = test_data)
                null_fitted_test <- fitted(null_model_test)
                
                test_null_residuals <- test_data$Expression - null_fitted_test
                test_model_residuals <- test_data$Expression - test_pred
                test_null_sq_diff <- test_null_residuals^2
                test_model_sq_diff <- test_model_residuals^2
                test_ranking_metric <- 1 - (test_model_sq_diff / pmax(test_null_sq_diff, 1e-8))
            } else {
                # Use influence score for test data
                test_ranking_metric <- influence_metrics_test$influence_score
            }
            
            # Sort test cells by ranking metric and apply same top-X% classification
            test_sorted_indices <- order(test_ranking_metric, decreasing = TRUE, na.last = TRUE)
            test_target_cells <- round(nrow(test_data) * dev_explained_train)  # Use same percentage as training
            test_target_cells <- max(1, min(test_target_cells, nrow(test_data) - 1))  # Ensure valid range
            
            # Classify test cells using same method as training
            test_deviance_cells_idx <- test_sorted_indices[1:test_target_cells]
            test_non_deviance_cells_idx <- test_sorted_indices[(test_target_cells+1):length(test_sorted_indices)]
            
            # Calculate test performance metrics
            if (length(test_deviance_cells_idx) > 0 && length(test_non_deviance_cells_idx) > 0) {
                # Calculate RMSE for TDE and Non-TDE groups
                test_residuals <- test_data$Expression - test_pred
                test_dec_rmse <- sqrt(mean(test_residuals[test_deviance_cells_idx]^2, na.rm = TRUE))
                test_non_dec_rmse <- sqrt(mean(test_residuals[test_non_deviance_cells_idx]^2, na.rm = TRUE))
                
                validation_results[[iter]] <- data.frame(
                    iteration = iter,
                    sample = current_sample,
                    train_dev_explained = dev_explained_train,
                    train_n = nrow(train_data),
                    test_n = nrow(test_data),
                    test_dec_n = length(test_deviance_cells_idx),
                    test_non_dec_n = length(test_non_deviance_cells_idx),
                    test_dec_rmse = test_dec_rmse,
                    test_non_dec_rmse = test_non_dec_rmse,
                    rmse_difference = test_non_dec_rmse - test_dec_rmse,
                    used_k = optimal_k,
                    used_lambda = optimal_lambda,
                    train_edf = summary(train_model)$edf
                )
            }
            
        }, error = function(e) {
            cat("Error in iteration", iter, ":", conditionMessage(e), "\n")
        })
    }
    
    # Combine results
    if (length(validation_results) > 0) {
        combined_results <- do.call(rbind, validation_results)
        
        # Summary statistics
        cat("\n--- Cross-Validation Summary ---\n")
        cat("Used k =", optimal_k, ", lambda =", optimal_lambda, "\n")
        cat("Successful iterations:", nrow(combined_results), "/", n_iterations, "\n")
        cat("Mean RMSE difference (Non-TDE - DEC):", round(mean(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("SD RMSE difference:", round(sd(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("Mean training EDF:", round(mean(combined_results$train_edf, na.rm = TRUE), 2), "\n")
        
        return(combined_results)
    } else {
        cat("No successful validation iterations\n")
        return(NULL)
    }
}

# Perform train-test validation with random cell classification
validate_random_deviance_cells <- function(current_sample, train_prop = 0.7, n_iterations = 20) {
    cat("\n=== Random Negative Control Cross-Validation for", current_sample, "===\n")
    
    # Use same parameters and data setup as validate_deviance_cells
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Get sample data using exact same method as analyze_sample
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    original_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    original_dev_explained <- summary(original_model)$dev.expl
    
    if (nrow(sample_data) < 50) {
        cat("Insufficient data points for meaningful CV (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    validation_results <- list()
    
    for (iter in 1:n_iterations) {
        set.seed(6 + iter)  # Reproducible splits
        
        # Create train-test split
        n_total <- nrow(sample_data)
        train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
        test_idx <- setdiff(1:n_total, train_idx)
        
        train_data <- sample_data[train_idx, ]
        test_data <- sample_data[test_idx, ]
        
        tryCatch({
            # Fit GAM on training data with fixed k and lambda
            train_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                               data = train_data, method = "REML", select = TRUE, 
                               gamma = 1.5, sp = optimal_lambda)
            
            # Calculate deviance explained cells on training data
            null_model_train <- gam(Expression ~ 1, data = train_data)
            null_fitted_train <- fitted(null_model_train)
            model_fitted_train <- fitted(train_model)
            null_residuals_train <- train_data$Expression - null_fitted_train
            model_residuals_train <- train_data$Expression - model_fitted_train
            
            null_sq_diff_train <- null_residuals_train^2
            model_sq_diff_train <- model_residuals_train^2
            explanatory_power_train <- 1 - (model_sq_diff_train / pmax(null_sq_diff_train, 1e-8))
            
            # Randomly select top X% cells from training (instead of EP-based)
            dev_explained_train <- summary(train_model)$dev.expl
            target_cells_train <- round(nrow(train_data) * dev_explained_train)
            top_cells_train_idx <- sample(1:nrow(train_data), size = target_cells_train)
            
            # Predict on test data
            test_pred <- predict(train_model, newdata = test_data)
            test_residuals <- test_data$Expression - test_pred
            
            # Randomly classify test cells (same proportion as training)
            test_top_like <- sample(1:nrow(test_data), size = round(nrow(test_data) * dev_explained_train))
            test_top_like <- 1:nrow(test_data) %in% test_top_like
            test_bottom_like <- !test_top_like
            
            # Calculate test performance metrics for random top vs bottom cells
            if (sum(test_top_like) > 0 && sum(test_bottom_like) > 0) {
                # Calculate RMSE for random top and bottom groups
                test_top_rmse <- sqrt(mean(test_residuals[test_top_like]^2, na.rm = TRUE))
                test_bottom_rmse <- sqrt(mean(test_residuals[test_bottom_like]^2, na.rm = TRUE))
                
                validation_results[[iter]] <- data.frame(
                    iteration = iter,
                    sample = current_sample,
                    train_dev_explained = dev_explained_train,
                    train_n = nrow(train_data),
                    test_n = nrow(test_data),
                    test_dec_n = sum(test_top_like),
                    test_non_dec_n = sum(test_bottom_like),
                    test_dec_rmse = test_top_rmse,
                    test_non_dec_rmse = test_bottom_rmse,
                    rmse_difference = test_bottom_rmse - test_top_rmse,
                    used_k = optimal_k,
                    used_lambda = optimal_lambda,
                    train_edf = summary(train_model)$edf
                )
            }
            
        }, error = function(e) {
            cat("Error in iteration", iter, ":", conditionMessage(e), "\n")
        })
    }
    
    # Combine results
    if (length(validation_results) > 0) {
        combined_results <- do.call(rbind, validation_results)
        
        # Summary statistics
        cat("\n--- Random Cross-Validation Summary ---\n")
        cat("Used k =", optimal_k, ", lambda =", optimal_lambda, "\n")
        cat("Successful iterations:", nrow(combined_results), "/", n_iterations, "\n")
        cat("Mean RMSE difference (Random Non-TDE - Random DEC):", round(mean(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("SD RMSE difference:", round(sd(combined_results$rmse_difference, na.rm = TRUE), 4), "\n")
        cat("Mean training EDF:", round(mean(combined_results$train_edf, na.rm = TRUE), 2), "\n")
        
        return(combined_results)
    } else {
        cat("No successful random validation iterations\n")
        return(NULL)
    }
}

# Function to perform statistical testing for each sample
perform_sample_statistical_test <- function(validation_results, sample_name) {
    if (is.null(validation_results) || nrow(validation_results) == 0) {
        return(NULL)
    }
    
    rmse_differences <- validation_results$rmse_difference
    rmse_differences <- rmse_differences[!is.na(rmse_differences)]
    
    if (length(rmse_differences) < 3) {
        cat("Insufficient data for statistical testing in", sample_name, "\n")
        return(NULL)
    }
    
    # Test for normality using Shapiro-Wilk
    normality_test <- shapiro.test(rmse_differences)
    is_normal <- normality_test$p.value > 0.05
    
    # Perform appropriate test
    if (is_normal) {
        # Use one-sample t-test
        stat_test <- t.test(rmse_differences, mu = 0)
        test_method <- "One-sample t-test"
        test_statistic <- stat_test$statistic
        test_pvalue <- stat_test$p.value
        conf_int <- stat_test$conf.int
    } else {
        # Use Wilcoxon signed-rank test (one-sample)
        stat_test <- wilcox.test(rmse_differences, mu = 0)
        test_method <- "Wilcoxon signed-rank test"
        test_statistic <- stat_test$statistic
        test_pvalue <- stat_test$p.value
        conf_int <- c(NA, NA)  # Wilcoxon doesn't provide confidence interval
    }
    
    # Calculate effect size (Cohen's d for t-test, r for Wilcoxon)
    if (is_normal) {
        effect_size <- mean(rmse_differences) / sd(rmse_differences)
        effect_size_name <- "Cohen's d"
    } else {
        effect_size <- abs(qnorm(test_pvalue/2)) / sqrt(length(rmse_differences))
        effect_size_name <- "Effect size r"
    }
    
    cat("\n--- Statistical Test Results for", sample_name, "---\n")
    cat("Shapiro-Wilk normality test p-value:", format.pval(normality_test$p.value, digits = 3), "\n")
    cat("Data distribution:", if(is_normal) "Normal" else "Non-normal", "\n")
    cat("Statistical test used:", test_method, "\n")
    cat("Test statistic:", round(as.numeric(test_statistic), 4), "\n")
    cat("P-value:", format.pval(test_pvalue, digits = 3), "\n")
    if (!is.na(conf_int[1])) {
        cat("95% CI: [", round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]\n")
    }
    cat(effect_size_name, ":", round(effect_size, 4), "\n")
    
    return(data.frame(
        sample = sample_name,
        n_iterations = length(rmse_differences),
        mean_rmse_diff = mean(rmse_differences),
        sd_rmse_diff = sd(rmse_differences),
        shapiro_pvalue = normality_test$p.value,
        is_normal = is_normal,
        test_method = test_method,
        test_statistic = as.numeric(test_statistic),
        test_pvalue = test_pvalue,
        conf_int_lower = conf_int[1],
        conf_int_upper = conf_int[2],
        effect_size = effect_size,
        effect_size_name = effect_size_name
    ))
}

# Function to create validation plots
create_validation_plots <- function(validation_results, sample_name) {
    if (is.null(validation_results) || nrow(validation_results) == 0) {
        cat("No validation results to plot for", sample_name, "\n")
        return(NULL)
    }
    
    # Plot 1: RMSE differences across iterations
    p1 <- ggplot(validation_results, aes(x = iteration, y = rmse_difference)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
        labs(title = paste("RMSE Difference Across CV Iterations (", sample_name, ")"),
             subtitle = "Positive values indicate better prediction for TDE cells",
             x = "CV Iteration", y = "RMSE Difference (Non-TDE - DEC)") +
        theme_minimal()
    
    # Plot 2: Distribution of RMSE differences
    p2 <- ggplot(validation_results, aes(x = rmse_difference)) +
        geom_histogram(bins = 10, fill = "#4B0082", alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        geom_vline(xintercept = mean(validation_results$rmse_difference), linetype = "solid", color = "blue") +
        labs(title = paste("Distribution of RMSE Differences (", sample_name, ")"),
             subtitle = "Red line at 0, blue line at mean",
             x = "RMSE Difference (Non-TDE - DEC)", y = "Frequency") +
        theme_minimal()
    
    print(p1)
    print(p2)
    
    return(list(iteration_plot = p1, distribution_plot = p2))
}

# Run validation on all samples
cat("\n", rep("=", 60), "\n")
cat("GAM-ALIGNED VALIDATION WITH PROPER STATISTICAL TESTING AND RANDOM CONTROL\n")
cat(rep("=", 60), "\n")

all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                       "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                       "HYW_5755_Tumor")
all_validation_results <- list()
all_random_validation_results <- list()  # Store random validation results
all_statistical_tests <- list()
all_random_statistical_tests <- list()  # Store random statistical tests

# Validate for all tumor samples
for (sample_name in all_tumor_samples) {
    if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
        cat("\n", rep("=", 60), "\n")
        # Original EP-based validation
        validation_result <- validate_deviance_cells(sample_name, train_prop = 0.7, n_iterations = 20)
        
        if (!is.null(validation_result)) {
            all_validation_results[[sample_name]] <- validation_result
            
            # Perform statistical test for this sample
            stat_test_result <- perform_sample_statistical_test(validation_result, sample_name)
            if (!is.null(stat_test_result)) {
                all_statistical_tests[[sample_name]] <- stat_test_result
            }
        }
        
        # Random negative control validation
        random_validation_result <- validate_random_deviance_cells(sample_name, train_prop = 0.7, n_iterations = 20)
        
        if (!is.null(random_validation_result)) {
            all_random_validation_results[[sample_name]] <- random_validation_result
            
            # Perform statistical test for random control
            random_stat_test_result <- perform_sample_statistical_test(random_validation_result, paste0(sample_name, "_Random"))
            if (!is.null(random_stat_test_result)) {
                all_random_statistical_tests[[sample_name]] <- random_stat_test_result
            }
        }
    } else {
        cat("\n--- No results found for", sample_name, "---\n")
    }
}

# Combined analysis with proper FDR correction
if (length(all_statistical_tests) > 0 || length(all_random_statistical_tests) > 0) {
    # Combine original statistical tests
    combined_tests <- do.call(rbind, all_statistical_tests)
    if (!is.null(combined_tests)) {
        # Apply FDR correction across samples
        combined_tests$fdr_corrected_pvalue <- p.adjust(combined_tests$test_pvalue, method = "BH")
        
        # Map sample names to patient IDs
        sample_to_patient <- c(
            "HYW_4701_Tumor" = "Pt.1",
            "HYW_4847_Tumor" = "Pt.2", 
            "HYW_4880_Tumor" = "Pt.3",
            "HYW_4881_Tumor" = "Pt.4",
            "HYW_5386_Tumor" = "Pt.5",
            "HYW_5742_Tumor" = "Pt.6",
            "HYW_5755_Tumor" = "Pt.7"
        )
        combined_tests$patient_id <- sample_to_patient[combined_tests$sample]
    } else {
        combined_tests <- data.frame()
    }
    
    # Combine random statistical tests
    combined_random_tests <- do.call(rbind, all_random_statistical_tests)
    if (!is.null(combined_random_tests)) {
        combined_random_tests$fdr_corrected_pvalue <- p.adjust(combined_random_tests$test_pvalue, method = "BH")
        combined_random_tests$patient_id <- sample_to_patient[gsub("_Random$", "", combined_random_tests$sample)]
    } else {
        combined_random_tests <- data.frame()
    }
    
    cat("\n", rep("=", 60), "\n")
    cat("OVERALL VALIDATION SUMMARY WITH PROPER STATISTICAL TESTING\n")
    cat(rep("=", 60), "\n")
    
    # Print original results
    if (nrow(combined_tests) > 0) {
        cat("\n--- EP-based Classification Results ---\n")
        print(combined_tests[, c("patient_id", "mean_rmse_diff", "test_method", "test_pvalue", "fdr_corrected_pvalue", "effect_size")])
    }
    
    # Print random control results
    if (nrow(combined_random_tests) > 0) {
        cat("\n--- Random Negative Control Classification Results ---\n")
        print(combined_random_tests[, c("patient_id", "mean_rmse_diff", "test_method", "test_pvalue", "fdr_corrected_pvalue", "effect_size")])
    }
    
    # Overall statistics for original validation
    if (nrow(combined_tests) > 0) {
        cat("\n--- Overall Summary (EP-based) ---\n")
        cat("Number of samples tested:", nrow(combined_tests), "\n")
        cat("Samples with normal distribution:", sum(combined_tests$is_normal), "\n")
        cat("Samples with non-normal distribution:", sum(!combined_tests$is_normal), "\n")
        cat("Mean RMSE difference across all samples:", round(mean(combined_tests$mean_rmse_diff), 4), "\n")
        cat("Samples with significant results (p < 0.05):", sum(combined_tests$test_pvalue < 0.05), "\n")
        cat("Samples with significant results (FDR < 0.05):", sum(combined_tests$fdr_corrected_pvalue < 0.05), "\n")
    }
    
    # Overall statistics for random control
    if (nrow(combined_random_tests) > 0) {
        cat("\n--- Overall Summary (Random Negative Control) ---\n")
        cat("Number of samples tested:", nrow(combined_random_tests), "\n")
        cat("Samples with normal distribution:", sum(combined_random_tests$is_normal), "\n")
        cat("Samples with non-normal distribution:", sum(!combined_random_tests$is_normal), "\n")
        cat("Mean RMSE difference across all samples:", round(mean(combined_random_tests$mean_rmse_diff), 4), "\n")
        cat("Samples with significant results (p < 0.05):", sum(combined_random_tests$test_pvalue < 0.05), "\n")
        cat("Samples with significant results (FDR < 0.05):", sum(combined_random_tests$fdr_corrected_pvalue < 0.05), "\n")
    }
    
    # Create overall plots
    if (length(all_validation_results) > 0 || length(all_random_validation_results) > 0) {
        # Combine original validation results
        combined_all <- do.call(rbind, all_validation_results)
        if (!is.null(combined_all)) {
            combined_all$patient_id <- sample_to_patient[combined_all$sample]
        }
        
        # Combine random validation results
        combined_random_all <- do.call(rbind, all_random_validation_results)
        if (!is.null(combined_random_all)) {
            combined_random_all$patient_id <- sample_to_patient[combined_random_all$sample]
        }
        
        # Create data for EP-based TDE vs Non-TDE boxplot comparison
        validation_plot_data <- data.frame()
        for (sample_name in names(all_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data <- all_validation_results[[sample_name]]
            for (i in 1:nrow(sample_data)) {
                validation_plot_data <- rbind(validation_plot_data, 
                                              data.frame(
                                                  Sample = patient_id,
                                                  RMSE = sample_data$test_dec_rmse[i],
                                                  Group = "DEC"
                                              ))
                validation_plot_data <- rbind(validation_plot_data, 
                                              data.frame(
                                                  Sample = patient_id,
                                                  RMSE = sample_data$test_non_dec_rmse[i],
                                                  Group = "Non-TDE"
                                              ))
            }
        }
        
        # Create data for random TDE vs Non-TDE boxplot comparison
        random_validation_plot_data <- data.frame()
        for (sample_name in names(all_random_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data <- all_random_validation_results[[sample_name]]
            for (i in 1:nrow(sample_data)) {
                random_validation_plot_data <- rbind(random_validation_plot_data, 
                                                     data.frame(
                                                         Sample = patient_id,
                                                         RMSE = sample_data$test_dec_rmse[i],
                                                         Group = "Random DEC"
                                                     ))
                random_validation_plot_data <- rbind(random_validation_plot_data, 
                                                     data.frame(
                                                         Sample = patient_id,
                                                         RMSE = sample_data$test_non_dec_rmse[i],
                                                         Group = "Random Non-TDE"
                                                     ))
            }
        }
        
        # Get significance annotations for original results
        if (nrow(combined_tests) > 0) {
            sig_values <- ifelse(combined_tests$fdr_corrected_pvalue < 0.001, "***",
                                 ifelse(combined_tests$fdr_corrected_pvalue < 0.01, "**",
                                        ifelse(combined_tests$fdr_corrected_pvalue < 0.05, "ns", "ns")))
        } else {
            sig_values <- character(0)
        }
        
        # Get significance annotations for random results
        if (nrow(combined_random_tests) > 0) {
            random_sig_values <- ifelse(combined_random_tests$fdr_corrected_pvalue < 0.001, "***",
                                        ifelse(combined_random_tests$fdr_corrected_pvalue < 0.01, "**",
                                               ifelse(combined_random_tests$fdr_corrected_pvalue < 0.05, "ns", "ns")))
        } else {
            random_sig_values <- character(0)
        }
        
        # Overall test across all samples (original)
        if (!is.null(combined_all) && nrow(combined_all) > 0) {
            all_rmse_diffs <- combined_all$rmse_difference[!is.na(combined_all$rmse_difference)]
            overall_normality <- shapiro.test(all_rmse_diffs)
            
            if (overall_normality$p.value > 0.05) {
                overall_test <- t.test(all_rmse_diffs, mu = 0)
                overall_method <- "One-sample t-test"
            } else {
                overall_test <- wilcox.test(all_rmse_diffs, mu = 0)
                overall_method <- "Wilcoxon signed-rank test"
            }
            
            cat("\n--- Overall Test Across All Samples (EP-based) ---\n")
            cat("Overall normality test p-value:", format.pval(overall_normality$p.value, digits = 3), "\n")
            cat("Overall test method:", overall_method, "\n")
            cat("Overall test p-value:", format.pval(overall_test$p.value, digits = 3), "\n")
        }
        
        # Overall test across all samples (random)
        if (!is.null(combined_random_all) && nrow(combined_random_all) > 0) {
            all_random_rmse_diffs <- combined_random_all$rmse_difference[!is.na(combined_random_all$rmse_difference)]
            overall_random_normality <- shapiro.test(all_random_rmse_diffs)
            
            if (overall_random_normality$p.value > 0.05) {
                overall_random_test <- t.test(all_random_rmse_diffs, mu = 0)
                overall_random_method <- "One-sample t-test"
            } else {
                overall_random_test <- wilcox.test(all_random_rmse_diffs, mu = 0)
                overall_random_method <- "Wilcoxon signed-rank test"
            }
            
            cat("\n--- Overall Test Across All Samples (Random Negative Control) ---\n")
            cat("Overall normality test p-value:", format.pval(overall_random_normality$p.value, digits = 3), "\n")
            cat("Overall test method:", overall_random_method, "\n")
            cat("Overall test p-value:", format.pval(overall_random_test$p.value, digits = 3), "\n")
        }
        
        # Reorder the groups so Non-TDE appears first
        if (nrow(validation_plot_data) > 0) {
            validation_plot_data$Group <- factor(validation_plot_data$Group, levels = c("Non-TDE", "DEC"))
            
            # Calculate dynamic y-axis limits for original plot
            y_max <- max(validation_plot_data$RMSE, na.rm = TRUE) + 0.1
            y_min <- min(validation_plot_data$RMSE, na.rm = TRUE) - 0.05
            
            # Main Plot with enhanced styling
            p1 <- ggplot(validation_plot_data, aes(x = Sample, y = RMSE)) +
                geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.shape = NA, width = 0.6,
                             position = position_dodge(width = 0.8)) +
                geom_point(aes(color = Group), alpha = 0.5, size = 1.2,
                           position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15)) +
                scale_fill_manual(values = c("Non-TDE" = "#C0C0C0", "DEC" = "#9954C1")) +
                scale_color_manual(values = c("Non-TDE" = "#A0A0A0", "DEC" = "#4B0082")) +
                geom_text(data = data.frame(Sample = combined_tests$patient_id, 
                                            RMSE = rep(y_max - 0.02, nrow(combined_tests)),
                                            sig = sig_values),
                          aes(x = Sample, y = RMSE, label = sig), 
                          inherit.aes = FALSE, size = 4, fontface = "bold") +
                labs(title = "RMSE Classification Validity",
                     subtitle = "TDE vs Non-TDE cells | Tests if classification captures low RMSE cells | ** = FDR < 0.01",
                     x = "Sample", y = "RMSE", fill = "Classification", color = "Classification") +
                theme_minimal() +
                theme(
                    plot.background = element_rect(fill = "white", color = "#2E2E2E", size = 0.3),
                    panel.background = element_rect(fill = "white"),
                    panel.border = element_rect(color = "#2E2E2E", fill = NA, size = 0.3),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold"),
                    legend.position = "bottom"
                ) +
                scale_y_continuous(limits = c(y_min, y_max)) +
                guides(color = "none")
            
            print(p1)
        }
        
        # Random control plot
        if (nrow(random_validation_plot_data) > 0) {
            random_validation_plot_data$Group <- factor(random_validation_plot_data$Group, levels = c("Random Non-TDE", "Random DEC"))
            
            # Calculate dynamic y-axis limits for random plot
            y_max_random <- max(random_validation_plot_data$RMSE, na.rm = TRUE) + 0.1
            y_min_random <- min(random_validation_plot_data$RMSE, na.rm = TRUE) - 0.05
            
            p1_random <- ggplot(random_validation_plot_data, aes(x = Sample, y = RMSE)) +
                geom_boxplot(aes(fill = Group), alpha = 0.7, outlier.shape = NA, width = 0.6,
                             position = position_dodge(width = 0.8)) +
                geom_point(aes(color = Group), alpha = 0.5, size = 1.2,
                           position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15)) +
                scale_fill_manual(values = c("Random Non-TDE" = "#C2E6F6", "Random DEC" = "#6161B7")) +
                scale_color_manual(values = c("Random Non-TDE" = "#A6C8E2", "Random DEC" = "#000080")) +
                geom_text(data = data.frame(Sample = combined_random_tests$patient_id, 
                                            RMSE = rep(y_max_random - 0.02, nrow(combined_random_tests)),
                                            sig = random_sig_values),
                          aes(x = Sample, y = RMSE, label = sig), 
                          inherit.aes = FALSE, size = 4, fontface = "bold") +
                labs(title = "Random Control Classification",
                     subtitle = "Random TDE vs Random Non-TDE cells | Should show no significant difference",
                     x = "Sample", y = "RMSE", fill = "Classification", color = "Classification") +
                theme_minimal() +
                theme(
                    plot.background = element_rect(fill = "white", color = "#2E2E2E", size = 0.3),
                    panel.background = element_rect(fill = "white"),
                    panel.border = element_rect(color = "#2E2E2E", fill = NA, size = 0.3),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold"),
                    legend.position = "bottom"
                ) +
                scale_y_continuous(limits = c(y_min_random, y_max_random)) +
                guides(color = "none")
            
            print(p1_random)
        }
        
        # Plot 2: Effect sizes (original)
        if (nrow(combined_tests) > 0) {
            p2 <- ggplot(combined_tests, aes(x = reorder(patient_id, effect_size), y = effect_size, fill = fdr_corrected_pvalue < 0.01)) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkgreen"),
                                  name = "FDR < 0.05", labels = c("No", "Yes")) +
                labs(title = "Effect Sizes by Sample (EP-based)",
                     subtitle = "Green bars indicate FDR < 0.01",
                     x = "Sample", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            print(p2)
        }
        
        # Plot 3: Effect sizes (random)
        if (nrow(combined_random_tests) > 0) {
            p3 <- ggplot(combined_random_tests, aes(x = reorder(patient_id, effect_size), y = effect_size, fill = fdr_corrected_pvalue < 0.01)) +
                geom_col(alpha = 0.8) +
                scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkgreen"),
                                  name = "FDR < 0.05", labels = c("No", "Yes")) +
                labs(title = "Effect Sizes by Sample (Random Negative Control)",
                     subtitle = "Green bars indicate FDR < 0.01",
                     x = "Sample", y = "Effect Size") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            print(p3)
        }
        
        # Export results to Excel
        library(openxlsx)
        wb <- createWorkbook()
        
        # Add individual sample results with patient IDs (original)
        for (sample_name in names(all_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data_with_id <- all_validation_results[[sample_name]]
            sample_data_with_id$patient_id <- patient_id
            addWorksheet(wb, paste0(patient_id, "_TRPM4"))
            writeData(wb, paste0(patient_id, "_TRPM4"), sample_data_with_id)
        }
        
        # Add individual sample results with patient IDs (random)
        for (sample_name in names(all_random_validation_results)) {
            patient_id <- sample_to_patient[sample_name]
            sample_data_with_id <- all_random_validation_results[[sample_name]]
            sample_data_with_id$patient_id <- patient_id
            addWorksheet(wb, paste0(patient_id, "_Random"))
            writeData(wb, paste0(patient_id, "_Random"), sample_data_with_id)
        }
        
        # Add statistical test results (original)
        if (nrow(combined_tests) > 0) {
            addWorksheet(wb, "Statistical_Tests_TRPM4")
            writeData(wb, "Statistical_Tests_TRPM4", combined_tests)
        }
        
        # Add statistical test results (random)
        if (nrow(combined_random_tests) > 0) {
            addWorksheet(wb, "Statistical_Tests_Random")
            writeData(wb, "Statistical_Tests_Random", combined_random_tests)
        }
        
        # Add combined raw results with patient IDs (original)
        if (!is.null(combined_all)) {
            addWorksheet(wb, "Combined_Results_TRPM4")
            writeData(wb, "Combined_Results_TRPM4", combined_all)
        }
        
        # Add combined raw results with patient IDs (random)
        if (!is.null(combined_random_all)) {
            addWorksheet(wb, "Combined_Results_Random")
            writeData(wb, "Combined_Results_Random", combined_random_all)
        }
        
        # Save workbook
        saveWorkbook(wb, "[2]_Leverage_CV.xlsx", overwrite = TRUE)
    }
}

cat("\n", rep("=", 60), "\n")
cat("GAM-ALIGNED VALIDATION WITH PROPER STATISTICAL TESTING AND RANDOM CONTROL COMPLETE\n")
cat(rep("=", 60), "\n")

# Function to create GAM plots using leverage-based method
create_leverage_gam_validation_plots <- function(current_sample) {
    cat("\n=== Creating Leverage-Based GAM Visualization for", current_sample, "===\n")
    
    # Pre-optimized parameters from analyze_sample results
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Get sample data using exact same method as validate_deviance_cells
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    if (nrow(sample_data) < 10) {
        cat("Insufficient data points for plotting (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    # Fit GAM model with pre-optimized parameters
    gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                     data = sample_data, method = "REML", select = TRUE, 
                     gamma = 1.5, sp = optimal_lambda)
    
    # === LEVERAGE-BASED RANKING ===
    
    # Calculate leverage-based influence metrics
    cat("Calculating leverage-based influence metrics...\n")
    influence_metrics <- calculate_cell_influence(gam_model, sample_data)
    
    # Check if influence calculation was successful, use fallback if needed
    if (all(is.na(influence_metrics$influence_score))) {
        cat("Warning: Influence calculation failed, falling back to explanatory power method\n")
        
        # Fallback to explanatory power method
        null_model <- gam(Expression ~ 1, data = sample_data)
        null_fitted <- fitted(null_model)
        model_fitted <- fitted(gam_model)
        null_residuals <- sample_data$Expression - null_fitted
        model_residuals <- sample_data$Expression - model_fitted
        null_sq_diff <- null_residuals^2
        model_sq_diff <- model_residuals^2
        ranking_metric <- 1 - (model_sq_diff / pmax(null_sq_diff, 1e-8))
        method_used <- "Explanatory Power (Fallback)"
    } else {
        # Use influence score for ranking
        cat("Using influence score for cell ranking\n")
        ranking_metric <- influence_metrics$influence_score
        method_used <- "Leverage-Based Influence"
    }
    
    # === CLASSIFICATION METHOD ===
    
    # Sort cells by their ranking metric (highest to lowest)
    sorted_indices <- order(ranking_metric, decreasing = TRUE, na.last = TRUE)
    
    # Get model deviance explained and target number of cells
    model_dev_explained <- summary(gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    target_cells <- max(1, min(target_cells, nrow(sample_data) - 1))  # Ensure valid range
    
    # Take the top cells by ranking metric (TDE cells)
    deviance_cells_idx <- sorted_indices[1:target_cells]
    non_deviance_cells_idx <- sorted_indices[(target_cells+1):length(sorted_indices)]
    
    # === END OF LEVERAGE-BASED METHOD ===
    
    # Calculate dynamic axis limits with a small margin
    x_range <- range(sample_data$TRPM4, na.rm = TRUE)
    y_range <- range(sample_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)
    
    # Create a data frame for ggplot - ALL cells are classified
    plot_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(1:nrow(sample_data) %in% deviance_cells_idx, "Leverage DEC", "Leverage Non-TDE"),
        RankingMetric = ranking_metric
    )
    
    # Add a drawing order column to control which points appear on top
    plot_data$draw_order <- ifelse(plot_data$Group == "Leverage DEC", 2, 1)
    
    # Sort the data frame by the draw order
    plot_data <- plot_data[order(plot_data$draw_order), ]
    
    # Create prediction data for the GAM line
    pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
    pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit <- pred$fit
    pred_data$se.fit <- pred$se.fit
    
    # Calculate percentage for subtitle
    dev_explained_percentage <- round(model_dev_explained * 100, 2)
    
    # Map sample names to patient IDs
    sample_to_patient <- c(
        "HYW_4701_Tumor" = "Pt.1",
        "HYW_4847_Tumor" = "Pt.2", 
        "HYW_4880_Tumor" = "Pt.3",
        "HYW_4881_Tumor" = "Pt.4",
        "HYW_5386_Tumor" = "Pt.5",
        "HYW_5742_Tumor" = "Pt.6",
        "HYW_5755_Tumor" = "Pt.7"
    )
    
    patient_id <- sample_to_patient[current_sample]
    
    # Create the plot with ggplot2
    p <- ggplot() +
        # Add points with colors and different transparency for purple vs gray
        geom_point(data = plot_data, 
                   aes(x = TRPM4, y = Expression, color = Group, alpha = Group),
                   size = 1.8) +
        scale_alpha_manual(values = c("Leverage DEC" = 0.4,        # More transparent purple
                                      "Leverage Non-TDE" = 0.6),   # Keep gray as is
                           guide = "none") +  # Hide alpha legend
        # Add the GAM line with orange color
        geom_line(data = pred_data,
                  aes(x = TRPM4, y = fit),
                  color = "#FFCC99",  # Orange
                  size = 1.2) +
        # Add confidence interval ribbon
        geom_ribbon(data = pred_data,
                    aes(x = TRPM4, 
                        ymin = fit - 1.96 * se.fit, 
                        ymax = fit + 1.96 * se.fit),
                    fill = "#FFCC99",  # Orange with transparency
                    alpha = 0.2) +
        # Set colors for the groups
        scale_color_manual(values = c("Leverage DEC" = "#4B0082",      # Purple
                                      "Leverage Non-TDE" = "#C0C0C0")) + # Gray
        # Styling
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position = "bottom",  # Show legend for classification
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            # Add thin black border around the plot area
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.line = element_blank()
        ) +
        # Set title and axis labels with sample name and dev explained percentage
        labs(
            title = paste("Leverage-Based: TRPM4 vs Ribo (", patient_id, ")", sep = ""),
            subtitle = paste("Method: ", method_used, " | Dev explained: ", dev_explained_percentage, "% | k=", optimal_k, ", EDF=", round(summary(gam_model)$edf, 2), " | Top ", target_cells, " cells", sep = ""),
            x = "TRPM4 Expression",
            y = "Ribo Expression",
            color = "Classification"
        ) +
        # Dynamically set axis limits with small margins
        scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
    
    # Display the plot
    print(p)
    
    # Print summary statistics
    cat("\n--- Classification Summary (Leverage-Based Method) ---\n")
    cat("Method used:", method_used, "\n")
    cat("Total cells:", nrow(sample_data), "\n")
    cat("Leverage TDE cells:", length(deviance_cells_idx), "\n")
    cat("Leverage Non-TDE cells:", length(non_deviance_cells_idx), "\n")
    cat("Classification completeness:", (length(deviance_cells_idx) + length(non_deviance_cells_idx)) == nrow(sample_data), "\n")
    cat("Model deviance explained:", round(model_dev_explained * 100, 2), "%\n")
    cat("Target cells (top X%):", target_cells, "out of", nrow(sample_data), "\n")
    cat("Effective degrees of freedom:", round(summary(gam_model)$edf, 2), "\n")
    
    # Verify ranking metric ranges
    cat("\n--- Ranking Metric Summary ---\n")
    cat("Leverage TDE cells - Min:", round(min(ranking_metric[deviance_cells_idx]), 4), 
        "| Max:", round(max(ranking_metric[deviance_cells_idx]), 4), 
        "| Mean:", round(mean(ranking_metric[deviance_cells_idx]), 4), "\n")
    cat("Leverage Non-TDE cells - Min:", round(min(ranking_metric[non_deviance_cells_idx]), 4), 
        "| Max:", round(max(ranking_metric[non_deviance_cells_idx]), 4),
        "| Mean:", round(mean(ranking_metric[non_deviance_cells_idx]), 4), "\n")
    
    return(list(plot = p, model = gam_model, classification_data = plot_data, method_used = method_used))
}

# Function to create leverage-based plots for all samples
create_all_leverage_gam_validation_plots <- function() {
    cat("\n", rep("=", 60), "\n")
    cat("CREATING LEVERAGE-BASED GAM VALIDATION PLOTS FOR ALL SAMPLES\n")
    cat(rep("=", 60), "\n")
    
    all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                           "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                           "HYW_5755_Tumor")
    
    all_leverage_plots <- list()
    
    for (sample_name in all_tumor_samples) {
        if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
            leverage_plot_result <- create_leverage_gam_validation_plots(sample_name)
            if (!is.null(leverage_plot_result)) {
                all_leverage_plots[[sample_name]] <- leverage_plot_result
            }
        } else {
            cat("\n--- No results found for", sample_name, "---\n")
        }
    }
    
    return(all_leverage_plots)
}

# Function to create leverage-based random control plots
create_leverage_random_gam_validation_plots <- function(current_sample, random_seed = 42) {
    cat("\n=== Creating Leverage-Based Random Control GAM Visualization for", current_sample, "===\n")
    
    # Pre-optimized parameters from analyze_sample results
    sample_params <- list(
        "HYW_4701_Tumor" = list(k = 6, lambda = 1.56970197914887),
        "HYW_4847_Tumor" = list(k = 3, lambda = 0.0793413231468075),
        "HYW_4880_Tumor" = list(k = 6, lambda = 0.138867472274075),
        "HYW_4881_Tumor" = list(k = 10, lambda = 0.52801317696145),
        "HYW_5386_Tumor" = list(k = 4, lambda = 0.287456401733474),
        "HYW_5742_Tumor" = list(k = 4, lambda = 5099.35774093852),
        "HYW_5755_Tumor" = list(k = 3, lambda = 1.28680173371137)
    )
    
    if (!current_sample %in% names(sample_params)) {
        cat("No pre-optimized parameters for", current_sample, "\n")
        return(NULL)
    }
    
    optimal_k <- sample_params[[current_sample]]$k
    optimal_lambda <- sample_params[[current_sample]]$lambda
    
    # Get sample data using exact same method as validate_deviance_cells
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells <- intersect(rownames(sample_data), selected_cells)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    if (nrow(sample_data) < 10) {
        cat("Insufficient data points for plotting (n =", nrow(sample_data), ")\n")
        return(NULL)
    }
    
    # Fit GAM model with pre-optimized parameters (same model as leverage method)
    gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = optimal_k), 
                     data = sample_data, method = "REML", select = TRUE, 
                     gamma = 1.5, sp = optimal_lambda)
    
    # Get model deviance explained to determine how many cells to randomly select
    model_dev_explained <- summary(gam_model)$dev.expl
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    
    # === RANDOM CLASSIFICATION (NEGATIVE CONTROL) ===
    set.seed(random_seed)  # For reproducible random classification
    
    # Randomly select the same percentage of cells as the leverage method would
    random_deviance_cells_idx <- sample(1:nrow(sample_data), size = target_cells)
    random_non_deviance_cells_idx <- setdiff(1:nrow(sample_data), random_deviance_cells_idx)
    
    # Calculate leverage-based metrics for verification (same calculation as leverage method)
    influence_metrics <- calculate_cell_influence(gam_model, sample_data)
    
    if (all(is.na(influence_metrics$influence_score))) {
        # Fallback to explanatory power method for verification
        null_model <- gam(Expression ~ 1, data = sample_data)
        null_fitted <- fitted(null_model)
        model_fitted <- fitted(gam_model)
        null_residuals <- sample_data$Expression - null_fitted
        model_residuals <- sample_data$Expression - model_fitted
        null_sq_diff <- null_residuals^2
        model_sq_diff <- model_residuals^2
        ranking_metric <- 1 - (model_sq_diff / pmax(null_sq_diff, 1e-8))
        method_used <- "Random + EP Verification"
    } else {
        ranking_metric <- influence_metrics$influence_score
        method_used <- "Random + Leverage Verification"
    }
    
    # === END OF RANDOM CLASSIFICATION ===
    
    # Calculate dynamic axis limits with a small margin
    x_range <- range(sample_data$TRPM4, na.rm = TRUE)
    y_range <- range(sample_data$Expression, na.rm = TRUE)
    x_margin <- 0.05 * diff(x_range)
    y_margin <- 0.05 * diff(y_range)
    
    # Create a data frame for ggplot - ALL cells are classified randomly
    plot_data <- data.frame(
        TRPM4 = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group = ifelse(1:nrow(sample_data) %in% random_deviance_cells_idx, "Random Leverage DEC", "Random Leverage Non-TDE"),
        RankingMetric = ranking_metric
    )
    
    # Add a drawing order column to control which points appear on top
    plot_data$draw_order <- ifelse(plot_data$Group == "Random Leverage DEC", 2, 1)
    
    # Sort the data frame by the draw order
    plot_data <- plot_data[order(plot_data$draw_order), ]
    
    # Create prediction data for the GAM line (same as leverage method)
    pred_data <- data.frame(TRPM4 = seq(min(plot_data$TRPM4), max(plot_data$TRPM4), length.out = 1000))
    pred <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit <- pred$fit
    pred_data$se.fit <- pred$se.fit
    
    # Calculate percentage for subtitle
    dev_explained_percentage <- round(model_dev_explained * 100, 2)
    
    # Map sample names to patient IDs
    sample_to_patient <- c(
        "HYW_4701_Tumor" = "Pt.1",
        "HYW_4847_Tumor" = "Pt.2", 
        "HYW_4880_Tumor" = "Pt.3",
        "HYW_4881_Tumor" = "Pt.4",
        "HYW_5386_Tumor" = "Pt.5",
        "HYW_5742_Tumor" = "Pt.6",
        "HYW_5755_Tumor" = "Pt.7"
    )
    
    patient_id <- sample_to_patient[current_sample]
    
    # Create the plot with ggplot2 (similar styling but different colors for random)
    p <- ggplot() +
        # Add points with colors and different transparency
        geom_point(data = plot_data, 
                   aes(x = TRPM4, y = Expression, color = Group, alpha = Group),
                   size = 1.8) +
        scale_alpha_manual(values = c("Random Leverage DEC" = 0.25,           # More transparent purple
                                      "Random Leverage Non-TDE" = 0.6),     # Keep gray as is
                           guide = "none") +  # Hide alpha legend
        # Add the GAM line with orange color (same model as leverage method)
        geom_line(data = pred_data,
                  aes(x = TRPM4, y = fit),
                  color = "#FFCC99",  # Orange
                  size = 1.2) +
        # Add confidence interval ribbon
        geom_ribbon(data = pred_data,
                    aes(x = TRPM4, 
                        ymin = fit - 1.96 * se.fit, 
                        ymax = fit + 1.96 * se.fit),
                    fill = "#FFCC99",  # Orange with transparency
                    alpha = 0.2) +
        # Set colors for random groups
        scale_color_manual(values = c("Random Leverage DEC" = "#4B0082",        # Purple
                                      "Random Leverage Non-TDE" = "#C0C0C0")) +  # Gray
        # Styling
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#EEEEEE"),
            panel.grid.minor = element_line(color = "#F5F5F5"),
            legend.position = "bottom",  # Show legend for classification
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            # Add thin black border around the plot area
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.line = element_blank()
        ) +
        # Set title and axis labels with sample name and dev explained percentage
        labs(
            title = paste("Random Control Leverage: TRPM4 vs Ribo (", patient_id, ")", sep = ""),
            subtitle = paste("Method: ", method_used, " | Same model: ", dev_explained_percentage, "% | k=", optimal_k, ", EDF=", round(summary(gam_model)$edf, 2), " | Random ", target_cells, " cells", sep = ""),
            x = "TRPM4 Expression",
            y = "Ribo Expression",
            color = "Random Classification"
        ) +
        # Dynamically set axis limits with small margins
        scale_x_continuous(limits = c(x_range[1] - x_margin, x_range[2] + x_margin)) +
        scale_y_continuous(limits = c(y_range[1] - y_margin, y_range[2] + y_margin))
    
    # Display the plot
    print(p)
    
    # Print summary statistics
    cat("\n--- Random Classification Summary (Leverage-Based Negative Control) ---\n")
    cat("Method used:", method_used, "\n")
    cat("Total cells:", nrow(sample_data), "\n")
    cat("Random Leverage TDE cells:", length(random_deviance_cells_idx), "\n")
    cat("Random Leverage Non-TDE cells:", length(random_non_deviance_cells_idx), "\n")
    cat("Classification completeness:", (length(random_deviance_cells_idx) + length(random_non_deviance_cells_idx)) == nrow(sample_data), "\n")
    cat("Model deviance explained (same as leverage method):", round(model_dev_explained * 100, 2), "%\n")
    cat("Random target cells:", target_cells, "out of", nrow(sample_data), "\n")
    cat("Effective degrees of freedom:", round(summary(gam_model)$edf, 2), "\n")
    
    # Compare ranking metric between random groups (should be similar)
    cat("\n--- Ranking Metric Summary (Random Groups) ---\n")
    cat("Random Leverage TDE cells - Min:", round(min(ranking_metric[random_deviance_cells_idx]), 4), 
        "| Max:", round(max(ranking_metric[random_deviance_cells_idx]), 4),
        "| Mean:", round(mean(ranking_metric[random_deviance_cells_idx]), 4), "\n")
    cat("Random Leverage Non-TDE cells - Min:", round(min(ranking_metric[random_non_deviance_cells_idx]), 4), 
        "| Max:", round(max(ranking_metric[random_non_deviance_cells_idx]), 4),
        "| Mean:", round(mean(ranking_metric[random_non_deviance_cells_idx]), 4), "\n")
    cat("Ranking metric difference between groups:", round(mean(ranking_metric[random_deviance_cells_idx]) - mean(ranking_metric[random_non_deviance_cells_idx]), 4), "(should be near zero)\n")
    
    return(list(plot = p, model = gam_model, classification_data = plot_data, method_used = method_used))
}

# Function to create leverage-based random control plots for all samples
create_all_leverage_random_gam_validation_plots <- function(random_seed = 42) {
    cat("\n", rep("=", 60), "\n")
    cat("CREATING LEVERAGE-BASED RANDOM CONTROL GAM VALIDATION PLOTS FOR ALL SAMPLES\n")
    cat(rep("=", 60), "\n")
    
    all_tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                           "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                           "HYW_5755_Tumor")
    
    all_leverage_random_plots <- list()
    
    for (sample_name in all_tumor_samples) {
        if (exists("pca_results") && !is.null(pca_results[[sample_name]][["Ribo"]])) {
            leverage_random_plot_result <- create_leverage_random_gam_validation_plots(sample_name, random_seed = random_seed)
            if (!is.null(leverage_random_plot_result)) {
                all_leverage_random_plots[[sample_name]] <- leverage_random_plot_result
            }
        } else {
            cat("\n--- No results found for", sample_name, "---\n")
        }
    }
    
    return(all_leverage_random_plots)
}

# Run the leverage-based plotting functions for all samples
all_leverage_gam_plots <- create_all_leverage_gam_validation_plots()

cat("\n", rep("=", 60), "\n")
cat("LEVERAGE-BASED GAM VALIDATION PLOTS COMPLETE\n")
cat(rep("=", 60), "\n")

# Run the leverage-based random control plotting function for all samples
all_leverage_random_gam_plots <- create_all_leverage_random_gam_validation_plots(random_seed = 42)

cat("\n", rep("=", 60), "\n")
cat("LEVERAGE-BASED RANDOM CONTROL GAM VALIDATION PLOTS COMPLETE\n")

cat(rep("=", 60), "\n")



