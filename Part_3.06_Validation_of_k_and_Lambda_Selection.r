############################################################################################
# Part 3.06: Validation of Selected k (by PRSS) and Lambda (by REML) Values by Refitting with Fixed k or Lambda
############################################################################################

# ============================================================================
# 1. Validation #1: Refitting with Different k Values
# ============================================================================
# Example model
best_model <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_model

# Calculates PRSS components including residual sum of squares and penalty terms.
calculate_prss <- function(model, data) {
    # Extracts the smooth component from the model for penalty calculation.
    smooth_comp <- model$smooth[[1]]
    
    # Retrieves the penalty matrix for smoothness penalization.
    S <- smooth_comp$S[[1]]
    
    # Obtains the smoothing parameter lambda from the model.
    lambda <- model$sp[1]
    
    # Gets the basis dimension k from the smooth term.
    k <- smooth_comp$bs.dim
    
    # Extracts coefficients corresponding to the smooth terms.
    coefs <- model$coefficients[smooth_comp$first.para:smooth_comp$last.para]
    
    # Computes the residual sum of squares from model predictions.
    fitted_values <- predict(model, newdata = data)
    residuals <- data$Expression - fitted_values
    RSS <- sum(residuals^2)
    
    # Calculates the smoothness penalty term using the penalty matrix.
    f_double_prime_integral <- as.numeric(t(coefs) %*% S %*% coefs)
    
    # Combines RSS and penalty to compute the penalized residual sum of squares.
    PRSS <- RSS + lambda * f_double_prime_integral
    
    # Returns a list of PRSS components for analysis.
    return(list(
        RSS = RSS,
        f_double_prime_integral = f_double_prime_integral,
        PRSS = PRSS,
        lambda = lambda
    ))
}

# Extracts PRSS metrics for a range of k values with optional fixed lambda.
extract_prss_metrics <- function(original_data, k_values, lambda = NULL) {
    # Initializes a data frame to store PRSS metrics for each k value.
    results <- data.frame(
        k_value = k_values,
        prss = numeric(length(k_values)),
        rss = numeric(length(k_values)),
        penalty = numeric(length(k_values)),
        lambda = numeric(length(k_values)),
        edf = numeric(length(k_values)),
        reml_score = numeric(length(k_values)),
        dev_explained = numeric(length(k_values)),
        r_squared = numeric(length(k_values))
    )
    
    # Iterates over k values to fit models and extract metrics.
    for (i in seq_along(k_values)) {
        k <- k_values[i]
        
        # Constructs a formula with the current k value for model fitting.
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k, ")")
        current_formula <- as.formula(formula_string)
        
        # Fits the GAM model with error handling for robustness.
        model <- tryCatch({
            if (is.null(lambda)) {
                gam(
                    formula = current_formula,
                    data = original_data,
                    method = "REML",
                    select = TRUE,
                    gamma = 1.5
                )
            } else {
                gam(
                    formula = current_formula,
                    data = original_data,
                    method = "REML",
                    sp = lambda,
                    select = FALSE,
                    gamma = 1.5
                )
            }
        }, error = function(e) {
            cat("Error fitting model with k =", k, ":", conditionMessage(e), "\n")
            return(NULL)
        })
        
        # Processes metrics if the model fit is successful.
        if (!is.null(model)) {
            prss_components <- calculate_prss(model, original_data)
            model_summary <- summary(model)
            
            # Stores computed metrics in the results data frame.
            results$prss[i] <- prss_components$PRSS
            results$rss[i] <- prss_components$RSS
            results$penalty[i] <- prss_components$f_double_prime_integral
            results$lambda[i] <- prss_components$lambda
            results$edf[i] <- sum(model$edf)
            results$reml_score[i] <- model$gcv.ubre
            results$dev_explained[i] <- model_summary$dev.expl * 100
            results$r_squared[i] <- model_summary$r.sq
        }
    }
    
    # Returns the compiled PRSS metrics.
    return(results)
}

# Creates a streamlined plot of PRSS versus k values.
create_prss_plot <- function(prss_metrics) {
    # Plots PRSS against basis dimension k with a plasma color scheme.
    p <- ggplot(prss_metrics, aes(x = k_value, y = prss, color = factor(k_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = -1) +
        labs(
            title = "PRSS by Basis Dimension (k)",
            x = "Basis Dimension (k)",
            y = "Penalized Residual Sum of Squares",
            color = "k Value"
        ) +
        theme_minimal()
    
    # Returns the PRSS plot for display or further use.
    return(p)
}

# Analyzes PRSS across different k values and visualizes results.
analyze_prss_for_k_values <- function(original_data, k_values, fixed_lambda = NULL) {
    # Extracts PRSS metrics for the specified k values.
    prss_metrics <- extract_prss_metrics(original_data, k_values, fixed_lambda)
    
    # Generates a plot of PRSS versus k values.
    prss_plot <- create_prss_plot(prss_metrics)
    
    # Outputs the PRSS metrics to the console.
    cat("\nPRSS Metrics by Basis Dimension (k):\n")
    print(prss_metrics)
    
    # Determines the optimal k value based on minimum PRSS.
    min_prss_idx <- which.min(prss_metrics$prss)
    optimal_k <- prss_metrics$k_value[min_prss_idx]
    
    # Reports the optimal k and its associated metrics.
    cat("\nOptimal basis dimension (k) by PRSS criterion:", optimal_k, "\n")
    cat("PRSS metrics for optimal k:\n")
    print(prss_metrics[min_prss_idx, ])
    
    # Returns a list of analysis results including metrics and plot.
    return(list(
        metrics = prss_metrics,
        plot = prss_plot,
        optimal_k = optimal_k,
        optimal_metrics = prss_metrics[min_prss_idx, ]
    ))
}

# Executes PRSS analysis with example data and displays results.
original_data <- best_model$model
k_values <- 3:10
prss_results <- analyze_prss_for_k_values(original_data, k_values)
print(prss_results$plot)


# =================================================================================
# 2. Validation #2: Refitting with Different Lambda Values
# =================================================================================
# Extracts eigenvalues and related metrics from GAM model components across lambda values.
extract_eigenvalues <- function(best_model, lambda_values) {
    # Initializes a data frame to store eigenvalue metrics for each lambda.
    results <- data.frame(
        lambda_value = lambda_values,
        reml_score = numeric(length(lambda_values)),
        min_eigen = numeric(length(lambda_values)),
        max_eigen = numeric(length(lambda_values)),
        condition_number = numeric(length(lambda_values)),
        effective_rank = numeric(length(lambda_values)),
        trace = numeric(length(lambda_values)),
        penalty_norm = numeric(length(lambda_values)),
        edf = numeric(length(lambda_values)),
        penalized_deviance = numeric(length(lambda_values))
    )
    
    # Preserves the original dataset for consistent refitting.
    original_data <- best_model$model
    
    # Retrieves the basis dimension k from the smooth term or sets a default value.
    if (length(best_model$smooth) > 0) {
        k_value <- best_model$smooth[[1]]$bs.dim
    } else {
        k_value <- 10
    }
    
    # Constructs a formula using the specified k value and dataset column names.
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k_value, ")")
        original_formula <- as.formula(formula_string)
    } else {
        original_formula <- best_model$formula
    }
    
    # Iterates over lambda values to extract metrics from refitted models.
    for (i in seq_along(lambda_values)) {
        lambda <- lambda_values[i]
        
        # Refits the GAM model with a fixed lambda using REML optimization.
        refitted_model <- gam(
            formula = original_formula,
            data = original_data,
            method = "REML",
            sp = lambda,
            select = FALSE,
            gamma = 1.5
        )
        
        # Records the REML score from the refitted model.
        reml_score <- refitted_model$gcv.ubre
        
        # Sums the effective degrees of freedom for the model.
        edf <- sum(refitted_model$edf)
        
        # Obtains the model’s deviance as a measure of fit.
        deviance <- refitted_model$deviance
        
        # Generates the model matrix for eigenvalue analysis.
        X <- model.matrix(refitted_model)
        
        # Processes the smooth component if available for penalty matrix extraction.
        if (length(refitted_model$smooth) > 0) {
            smooth <- refitted_model$smooth[[1]]
            first_idx <- smooth$first.para
            last_idx <- smooth$last.para
            S <- smooth$S[[1]]
            
            # Ensures valid smooth component indices and matrices for analysis.
            if (!is.null(S) && !is.null(first_idx) && !is.null(last_idx) && 
                first_idx <= last_idx && last_idx <= ncol(X)) {
                
                # Extracts the smooth term columns from the model matrix.
                X_s <- X[, first_idx:last_idx, drop = FALSE]
                
                # Scales the penalty matrix by the current lambda value.
                S_lambda <- lambda * S
                
                # Computes the cross-product matrix of smooth terms.
                XtX <- crossprod(X_s)
                
                # Forms the penalized normal equations matrix with the penalty term.
                penalized_matrix <- XtX + S_lambda
                
                # Calculates eigenvalues of the penalized matrix for stability analysis.
                eigvals <- eigen(penalized_matrix, symmetric = TRUE, only.values = TRUE)$values
                
                # Computes eigenvalues of the scaled penalty matrix for comparison.
                S_eigvals <- eigen(S_lambda, symmetric = TRUE, only.values = TRUE)$values
                
                # Filters out invalid eigenvalues for reliable metrics.
                eigvals <- eigvals[is.finite(eigvals)]
                
                # Processes eigenvalues if any are valid.
                if (length(eigvals) > 0) {
                    min_eigval <- min(eigvals[eigvals > 1e-10])
                    max_eigval <- max(eigvals)
                    cond_num <- max_eigval / min_eigval
                    eff_rank <- sum(eigvals > (max_eigval * 1e-10))
                    trace_val <- sum(eigvals)
                    penalty_norm <- sqrt(sum(S_eigvals^2))
                    beta_s <- refitted_model$coefficients[first_idx:last_idx]
                    penalty_value <- as.numeric(t(beta_s) %*% S_lambda %*% beta_s)
                    penalized_dev <- deviance + penalty_value
                    
                    # Stores computed metrics in the results data frame.
                    results$min_eigen[i] <- min_eigval
                    results$max_eigen[i] <- max_eigval
                    results$condition_number[i] <- cond_num
                    results$effective_rank[i] <- eff_rank
                    results$trace[i] <- trace_val
                    results$penalty_norm[i] <- penalty_norm
                    results$edf[i] <- edf
                    results$penalized_deviance[i] <- penalized_dev
                }
            }
        }
        
        # Assigns the REML score to the results.
        results$reml_score[i] <- reml_score
    }
    
    # Returns the compiled eigenvalue metrics.
    return(results)
}

# Extracts detailed eigenvalue distributions for each lambda value.
extract_detailed_eigenvalues <- function(best_model, lambda_values) {
    # Retains the original data for refitting consistency.
    original_data <- best_model$model
    
    # Determines the basis dimension k from the smooth term or uses a default.
    if (length(best_model$smooth) > 0) {
        k_value <- best_model$smooth[[1]]$bs.dim
    } else {
        k_value <- 10
    }
    
    # Builds a formula with the explicit k value for refitting.
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k_value, ")")
        original_formula <- as.formula(formula_string)
    } else {
        original_formula <- best_model$formula
    }
    
    # Initializes a list to store eigenvalue distributions.
    eigenvalues_list <- list()
    
    # Loops through lambda values to extract detailed eigenvalues.
    for (i in seq_along(lambda_values)) {
        lambda <- lambda_values[i]
        
        # Refits the GAM model with a fixed lambda using REML.
        refitted_model <- gam(
            formula = original_formula,
            data = original_data,
            method = "REML",
            sp = lambda,
            select = FALSE,
            gamma = 1.5
        )
        
        # Extracts the penalty matrix if a smooth term exists.
        if (length(refitted_model$smooth) > 0) {
            S <- refitted_model$smooth[[1]]$S[[1]]
            
            # Verifies the penalty matrix is valid before processing.
            if (!is.null(S) && nrow(S) > 0 && ncol(S) > 0) {
                S_lambda <- lambda * S
                eigvals <- eigen(S_lambda, symmetric = TRUE, only.values = TRUE)$values
                eigvals <- sort(eigvals[is.finite(eigvals)])
                eigenvalues_list[[as.character(lambda)]] <- eigvals
            }
        }
    }
    
    # Returns the list of eigenvalue distributions.
    return(eigenvalues_list)
}

# Extracts model components using an LME approach or fallback method.
extract_lme_components <- function(best_model, lambda_values) {
    # Preserves the original dataset for refitting.
    original_data <- best_model$model
    
    # Obtains the basis dimension k from the smooth term or sets a default.
    if (length(best_model$smooth) > 0) {
        k_value <- best_model$smooth[[1]]$bs.dim
    } else {
        k_value <- 10
    }
    
    # Formulates a model equation with the specified k value.
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k_value, ")")
        original_formula <- as.formula(formula_string)
    } else {
        original_formula <- best_model$formula
    }
    
    # Initializes a list to store model components for each lambda.
    components_list <- list()
    
    # Iterates over lambda values to extract components from refitted models.
    for (i in seq_along(lambda_values)) {
        lambda <- lambda_values[i]
        
        # Refits the GAM model with a fixed lambda using REML.
        refitted_model <- gam(
            formula = original_formula,
            data = original_data,
            method = "REML",
            sp = lambda,
            select = FALSE,
            gamma = 1.5
        )
        
        # Attempts to extract LME components with error handling.
        tryCatch({
            components <- extract.lme.components(refitted_model)
            components_list[[as.character(lambda)]] <- components
        }, error = function(e) {
            components <- list(
                X = model.matrix(refitted_model),
                S = if (length(refitted_model$smooth) > 0) refitted_model$smooth[[1]]$S else NULL,
                sp = refitted_model$sp,
                edf = refitted_model$edf
            )
            components_list[[as.character(lambda)]] <- components
        })
    }
    
    # Returns the list of extracted components.
    return(components_list)
}

# Creates bar plots to visualize eigenvalue metrics across lambda values.
create_eigenvalue_plots <- function(eigenvalue_data) {
    # Plots condition number for each lambda value.
    p1 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = condition_number, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Condition Number by Lambda",
            x = "Lambda (λ) Value",
            y = "Condition Number"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots effective rank for each lambda value.
    p2 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = effective_rank, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Effective Rank by Lambda",
            x = "Lambda (λ) Value",
            y = "Effective Rank"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots penalty matrix norm for each lambda value.
    p3 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = penalty_norm, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Penalty Matrix Norm by Lambda",
            x = "Lambda (λ) Value",
            y = "Penalty Norm"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots REML score for each lambda value.
    p4 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = reml_score, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "REML Score by Lambda",
            x = "Lambda (λ) Value",
            y = "REML Score"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots lambda versus REML score as an additional comparison.
    p5 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = reml_score, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Lambda vs REML Score",
            x = "Lambda (λ) Value",
            y = "REML Score"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots condition number versus REML score for further insight.
    p6 <- ggplot(eigenvalue_data, aes(x = factor(lambda_value), y = condition_number, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Condition Number vs REML Score",
            x = "Lambda (λ) Value",
            y = "Condition Number"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Arranges all bar plots in a grid layout for display.
    grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
}

# Plots detailed eigenvalue distributions across lambda values.
plot_eigenvalue_distributions <- function(eigenvalues_list) {
    # Initializes a list to store eigenvalue distribution plots.
    plot_list <- list()
    
    # Iterates over eigenvalue types for plotting.
    for (type_name in c("normal", "penalty", "penalized")) {
        if (!is.null(eigenvalues_list[[type_name]])) {
            type_list <- eigenvalues_list[[type_name]]
            eigenvalue_df <- data.frame()
            
            # Compiles eigenvalue data into a long-format data frame.
            for (lambda in names(type_list)) {
                eigvals <- type_list[[lambda]]
                
                if (length(eigvals) > 0) {
                    temp_df <- data.frame(
                        lambda = as.numeric(lambda),
                        eigenvalue_index = 1:length(eigvals),
                        eigenvalue = eigvals
                    )
                    eigenvalue_df <- rbind(eigenvalue_df, temp_df)
                }
            }
            
            # Generates plots if data is available.
            if (nrow(eigenvalue_df) > 0) {
                p1 <- ggplot(eigenvalue_df, aes(x = eigenvalue_index, y = eigenvalue, color = factor(lambda))) +
                    geom_line() +
                    geom_point() +
                    scale_color_viridis_d() +
                    labs(
                        title = paste0(type_name, " Eigenvalue Distribution by Lambda"),
                        x = "Eigenvalue Index",
                        y = "Eigenvalue",
                        color = "Lambda (λ)"
                    ) +
                    theme_minimal()
                
                p2 <- ggplot(eigenvalue_df, aes(x = eigenvalue_index, y = eigenvalue, color = factor(lambda))) +
                    geom_line() +
                    geom_point() +
                    scale_y_log10() +
                    scale_color_viridis_d() +
                    labs(
                        title = paste0(type_name, " Eigenvalue Distribution (Log Scale) by Lambda"),
                        x = "Eigenvalue Index",
                        y = "Eigenvalue (Log Scale)",
                        color = "Lambda (λ)"
                    ) +
                    theme_minimal()
                
                # Displays the linear and log-scale plots.
                print(p1)
                print(p2)
                
                plot_list[[paste0(type_name, "_linear")]] <- p1
                plot_list[[paste0(type_name, "_log")]] <- p2
                plot_list[[paste0(type_name, "_data")]] <- eigenvalue_df
            } else {
                cat("No valid eigenvalue data to plot for", type_name, "\n")
            }
        }
    }
    
    # Creates comparison plots if all eigenvalue types are available.
    if (!is.null(eigenvalues_list[["normal"]]) && !is.null(eigenvalues_list[["penalty"]]) && 
        !is.null(eigenvalues_list[["penalized"]])) {
        
        for (lambda in names(eigenvalues_list[["normal"]])) {
            comparison_df <- data.frame()
            
            normal_eigvals <- eigenvalues_list[["normal"]][[lambda]]
            if (length(normal_eigvals) > 0) {
                comparison_df <- rbind(comparison_df, data.frame(
                    type = "X'X",
                    eigenvalue_index = 1:length(normal_eigvals),
                    eigenvalue = normal_eigvals
                ))
            }
            
            penalty_eigvals <- eigenvalues_list[["penalty"]][[lambda]]
            if (length(penalty_eigvals) > 0) {
                comparison_df <- rbind(comparison_df, data.frame(
                    type = "λS",
                    eigenvalue_index = 1:length(penalty_eigvals),
                    eigenvalue = penalty_eigvals
                ))
            }
            
            penalized_eigvals <- eigenvalues_list[["penalized"]][[lambda]]
            if (length(penalized_eigvals) > 0) {
                comparison_df <- rbind(comparison_df, data.frame(
                    type = "X'X + λS",
                    eigenvalue_index = 1:length(penalized_eigvals),
                    eigenvalue = penalized_eigvals
                ))
            }
            
            # Produces comparison plots if data is present.
            if (nrow(comparison_df) > 0) {
                p_comp <- ggplot(comparison_df, aes(x = eigenvalue_index, y = eigenvalue, color = type)) +
                    geom_line() +
                    geom_point() +
                    scale_color_viridis_d() +
                    labs(
                        title = paste0("Eigenvalue Comparison (λ = ", lambda, ")"),
                        x = "Eigenvalue Index",
                        y = "Eigenvalue",
                        color = "Matrix Type"
                    ) +
                    theme_minimal()
                
                p_comp_log <- ggplot(comparison_df, aes(x = eigenvalue_index, y = eigenvalue, color = type)) +
                    geom_line() +
                    geom_point() +
                    scale_y_log10() +
                    scale_color_viridis_d() +
                    labs(
                        title = paste0("Eigenvalue Comparison (Log Scale, λ = ", lambda, ")"),
                        x = "Eigenvalue Index",
                        y = "Eigenvalue (Log Scale)",
                        color = "Matrix Type"
                    ) +
                    theme_minimal()
                
                print(p_comp)
                print(p_comp_log)
                
                plot_list[[paste0("comparison_", lambda)]] <- p_comp
                plot_list[[paste0("comparison_log_", lambda)]] <- p_comp_log
            }
        }
    }
    
    # Returns the list of eigenvalue distribution plots.
    return(plot_list)
}

# Generates additional plots to compare eigenvalue metrics.
create_additional_plots <- function(eigen_metrics) {
    # Plots condition number against REML score with a plasma color scheme.
    p1 <- ggplot(eigen_metrics, aes(x = condition_number, y = reml_score, color = factor(lambda_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, guide = guide_legend(
                          override.aes = list(
                            shape = 21,
                            fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(11),
                            color = "#111111",
                            size = 4,
                            stroke = 0.5
                            )
                            )) +
        labs(
            title = "Condition Number vs REML Score",
            x = "Condition Number",
            y = "REML Score",
            color = "Lambda (λ)"
        ) +
        theme_minimal()
    
    # Plots effective rank versus penalized deviance with labels.
    p2 <- ggplot(eigen_metrics, aes(x = effective_rank, y = penalized_deviance, color = factor(lambda_value))) +
        geom_point(size = 4) +
        geom_text(aes(label = lambda_value), hjust = -0.3, vjust = 0.3, size = 3) +
        scale_color_viridis_d() +
        labs(
            title = "Effective Rank vs Penalized Deviance",
            x = "Effective Rank",
            y = "Penalized Deviance",
            color = "Lambda (λ)"
        ) +
        theme_minimal()
    
    # Plots lambda against REML score using a plasma palette.
    p3 <- ggplot(eigen_metrics, aes(x = lambda_value, y = reml_score, color = factor(lambda_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, guide = guide_legend(
                          override.aes = list(
                            shape = 21,
                            fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(11),
                            color = "#111111",
                            size = 4,
                            stroke = 0.5
                            )
                            )) +
        labs(
            title = "Lambda vs REML Score",
            x = "Lambda (λ)",
            y = "REML Score",
            color = "Lambda (λ)"
        ) +
        theme_minimal()
    
    # Plots effective degrees of freedom against lambda.
    p4 <- ggplot(eigen_metrics, aes(x = lambda_value, y = edf)) +
        geom_point(size = 4, color = "darkgreen") +
        geom_line(color = "darkgreen") +
        labs(
            title = "Lambda vs Effective Degrees of Freedom",
            x = "Lambda (λ)",
            y = "EDF"
        ) +
        theme_minimal()
    
    # Plots bar chart of lambda versus REML score.
    p5 <- ggplot(eigen_metrics, aes(x = factor(lambda_value), y = reml_score, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Lambda vs REML Score",
            x = "Lambda (λ) Value",
            y = "REML Score"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Plots bar chart of condition number versus REML score.
    p6 <- ggplot(eigen_metrics, aes(x = factor(lambda_value), y = condition_number, fill = factor(lambda_value))) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        labs(
            title = "Condition Number vs REML Score",
            x = "Lambda (λ) Value",
            y = "Condition Number"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    # Displays all additional comparison plots.
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    
    # Returns the list of additional plots.
    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
}

# Analyzes eigenvalues across different lambda values and generates visualizations.
analyze_eigenvalues_for_lambdas <- function(best_model, lambda_values) {
    # Extracts eigenvalue metrics for the specified lambda values.
    eigen_metrics <- extract_eigenvalues(best_model, lambda_values)
    
    # Obtains detailed eigenvalue distributions for further analysis.
    detailed_eigenvalues <- extract_detailed_eigenvalues(best_model, lambda_values)
    
    # Creates bar plots of eigenvalue metrics.
    create_eigenvalue_plots(eigen_metrics)
    
    # Generates detailed eigenvalue distribution plots.
    eigen_dist_results <- plot_eigenvalue_distributions(detailed_eigenvalues)
    
    # Produces additional comparison plots for eigenvalue metrics.
    additional_plots <- create_additional_plots(eigen_metrics)
    
    # Outputs the eigenvalue metrics to the console.
    cat("\nEigenvalue Metrics:\n")
    print(eigen_metrics)
    
    # Identifies the optimal lambda based on the minimum REML score.
    min_reml_idx <- which.min(eigen_metrics$reml_score)
    optimal_lambda <- eigen_metrics$lambda_value[min_reml_idx]
    
    # Reports the optimal lambda and its associated metrics.
    cat("\nOptimal lambda by REML criterion:", optimal_lambda, "\n")
    cat("Eigenvalue metrics for optimal lambda:\n")
    print(eigen_metrics[min_reml_idx, ])
    
    # Returns a comprehensive list of analysis results.
    return(list(
        metrics = eigen_metrics,
        detailed_eigenvalues = detailed_eigenvalues,
        distribution_plots = eigen_dist_results,
        additional_plots = additional_plots,
        optimal_lambda = optimal_lambda,
        optimal_metrics = eigen_metrics[min_reml_idx, ]
    ))
}

# Executes the eigenvalue analysis with example data.
best_model <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_model
lambda_values <- c(0.264006588, 0.332627458, 0.419084336, 0.528013177, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
eigen_results <- analyze_eigenvalues_for_lambdas(best_model, lambda_values)
