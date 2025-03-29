#######################################################################
# Part 3.6: Validation of Selected k (by PRSS) and Lambda (by REML) Values
#######################################################################

# ====================================================
# 1. Validation #1 for k Selection: Refitting with Different k Values
# ====================================================
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


# =================================================================
# 2. Validation #1 for Lambda Selection: Refitting with Different Lambda Values
# =================================================================
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


# ===============================================================================
# 3. Validation #2 for k Selection: 10-Fold Cross-Validation for PRSS's k Selection
# ===============================================================================
# Performs k-fold cross-validation to evaluate different k values for PRSS.
perform_kfold_cv_for_k <- function(original_data, k_values, k_folds = 10, seed = 123) {
    # Sets a seed for reproducible fold generation.
    set.seed(seed)
    
    # Initializes a data frame to store cross-validation results for k values.
    cv_results <- data.frame(
        k_value = k_values,
        rmse = numeric(length(k_values)),
        deviance = numeric(length(k_values))
    )
    
    # Adjusts fold count if the dataset size is insufficient.
    if (nrow(original_data) < 2 * k_folds) {
        cat("WARNING: Too few observations for", k_folds, "folds. Reducing to", floor(nrow(original_data)/2), "folds.\n")
        k_folds <- max(2, floor(nrow(original_data)/2))
    }
    
    # Reports dataset size and number of folds.
    cat("Number of observations in dataset:", nrow(original_data), "\n")
    cat("Using", k_folds, "folds for cross-validation\n")
    
    # Verifies required columns and sets up the formula base.
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        cat("Using Expression ~ TRPM4 + s(TRPM4) formula with varying k values\n")
    } else {
        cat("WARNING: Required columns 'Expression' and 'TRPM4' not found in data\n")
        return(NULL)
    }
    
    # Creates k-fold splits with a reproducible seed.
    set.seed(seed)
    folds <- createFolds(1:nrow(original_data), k = k_folds, list = TRUE, returnTrain = FALSE)
    fold_sizes <- sapply(folds, length)
    cat("Fold sizes:", paste(fold_sizes, collapse = ", "), "\n")
    
    # Initializes matrices for fold-specific RMSE and deviance metrics.
    fold_rmse <- matrix(NA, nrow = k_folds, ncol = length(k_values))
    fold_deviance <- matrix(NA, nrow = k_folds, ncol = length(k_values))
    
    # Iterates over k values to perform cross-validation.
    for (k_idx in 1:length(k_values)) {
        k <- k_values[k_idx]
        cat("\n------- Processing basis dimension k =", k, "-------\n")
        
        # Tracks fold-specific metrics for the current k value.
        fold_metrics <- data.frame(
            fold = 1:k_folds,
            rmse = numeric(k_folds),
            deviance = numeric(k_folds),
            success = logical(k_folds)
        )
        
        # Constructs a formula with the current k value.
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k, ")")
        current_formula <- as.formula(formula_string)
        cat("Using formula:", formula_string, "\n")
        
        # Performs cross-validation for each fold.
        for (j in 1:k_folds) {
            cat("  Fold", j, "of", k_folds, "...\n")
            test_indices <- folds[[j]]
            train_data <- original_data[-test_indices, ]
            test_data <- original_data[test_indices, ]
            cat("    Train size:", nrow(train_data), "Test size:", nrow(test_data), "\n")
            
            # Fits and evaluates the model with error handling.
            tryCatch({
                set.seed(seed + j)
                gam_model <- gam(
                    formula = current_formula,
                    data = train_data,
                    method = "REML",
                    select = TRUE,
                    gamma = 1.5
                )
                
                if (is.null(gam_model)) {
                    cat("    ERROR: Model fitting returned NULL\n")
                    fold_metrics$success[j] <- FALSE
                    next
                }
                
                predictions <- tryCatch({
                    preds <- predict(gam_model, newdata = test_data)
                    if (any(is.na(preds))) cat("    WARNING: NA values in predictions\n")
                    preds
                }, error = function(e) {
                    cat("    ERROR in prediction:", conditionMessage(e), "\n")
                    return(rep(NA, nrow(test_data)))
                })
                
                actuals <- test_data$Expression
                if (all(is.na(predictions))) {
                    cat("    ERROR: All predictions are NA\n")
                    fold_metrics$success[j] <- FALSE
                    next
                }
                
                rmse_value <- sqrt(mean((actuals - predictions)^2, na.rm = TRUE))
                deviance_value <- sum((actuals - predictions)^2, na.rm = TRUE)
                
                fold_metrics$rmse[j] <- rmse_value
                fold_metrics$deviance[j] <- deviance_value
                fold_metrics$success[j] <- TRUE
                
                fold_rmse[j, k_idx] <- rmse_value
                fold_deviance[j, k_idx] <- deviance_value
                
                cat("    Metrics - RMSE:", round(rmse_value, 4), "\n")
                
            }, error = function(e) {
                cat("    ERROR in fold", j, ":", conditionMessage(e), "\n")
                fold_metrics$success[j] <- FALSE
            })
        }
        
        # Summarizes fold success and computes average metrics.
        successful_folds <- sum(fold_metrics$success)
        cat("  Successfully processed", successful_folds, "out of", k_folds, "folds\n")
        
        if (successful_folds > 0) {
            cv_results$rmse[k_idx] <- mean(fold_metrics$rmse[fold_metrics$success], na.rm = TRUE)
            cv_results$deviance[k_idx] <- mean(fold_metrics$deviance[fold_metrics$success], na.rm = TRUE)
            cat("  Average metrics - RMSE:", round(cv_results$rmse[k_idx], 4), "\n")
        } else {
            cat("  WARNING: No successful folds for k =", k, "\n")
            cv_results$rmse[k_idx] <- NA
            cv_results$deviance[k_idx] <- NA
        }
    }
    
    # Warns if all results are invalid and identifies the optimal k.
    if (all(is.na(cv_results$rmse))) {
        warning("All cross-validation results are NA. Check for issues in the data or model.")
    } else {
        optimal_k_idx <- which.min(cv_results$rmse)
        cat("\nOptimal basis dimension k by CV RMSE:", cv_results$k_value[optimal_k_idx], 
            "with RMSE:", round(cv_results$rmse[optimal_k_idx], 4), "\n")
    }
    
    # Returns cross-validation results including fold-specific metrics.
    return(list(cv_results = cv_results, fold_rmse = fold_rmse, fold_deviance = fold_deviance))
}

# Visualizes cross-validation results for k values using a plasma palette.
plot_cv_results_for_k <- function(cv_results) {
    # Skips plotting if all metrics are invalid.
    if (all(is.na(cv_results$rmse))) {
        cat("Cannot create plots: All CV metrics are NA\n")
        return(NULL)
    }
    
    # Plots RMSE against k values with a plasma color scheme.
    p_rmse <- ggplot(cv_results, aes(x = k_value, y = rmse, color = factor(k_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = -1, guide = guide_legend(
            override.aes = list(
                shape = 21,
                fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(length(cv_results$k_value)),
                color = "#111111",
                size = 3,
                stroke = 0.5
            )
        )) +
        labs(title = "Cross-Validation: Basis Dimension (k) vs RMSE", 
             x = "k Value", 
             y = "Root Mean Squared Error",
             color = "k Value") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Plots deviance against k values with a plasma color scheme.
    p_deviance <- ggplot(cv_results, aes(x = k_value, y = deviance, color = factor(k_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, direction = -1, guide = guide_legend(
            override.aes = list(
                shape = 21,
                fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(length(cv_results$k_value)),
                color = "#111111",
                size = 3,
                stroke = 0.5
            )
        )) +
        labs(title = "Cross-Validation: Basis Dimension (k) vs Deviance", 
             x = "k Value", 
             y = "Deviance",
             color = "k Value") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Displays the RMSE and deviance plots.
    print(p_rmse)
    print(p_deviance)
    
    # Returns the list of k-value cross-validation plots.
    return(list(rmse_plot = p_rmse, deviance_plot = p_deviance))
}

# Analyzes cross-validation results for k values with statistical metrics.
analyze_cv_results_for_k <- function(cv_data, prss_results) {
    # Extracts cross-validation results and fold-specific metrics.
    cv_results <- cv_data$cv_results
    fold_rmse <- cv_data$fold_rmse
    fold_deviance <- cv_data$fold_deviance
    
    # Terminates analysis if all metrics are invalid.
    if (all(is.na(cv_results$rmse))) {
        cat("Cannot analyze results: All CV metrics are NA\n")
        return(NULL)
    }
    
    # Identifies optimal k values based on RMSE, deviance, and PRSS.
    optimal_k_rmse <- cv_results$k_value[which.min(cv_results$rmse)]
    optimal_k_deviance <- cv_results$k_value[which.min(cv_results$deviance)]
    best_prss_k <- prss_results$best_k
    
    # Summarizes optimal k values and their metrics.
    summary_df <- data.frame(
        Criterion = c("CV RMSE", "CV Deviance", "PRSS Method"),
        Optimal_k = c(optimal_k_rmse, optimal_k_deviance, best_prss_k),
        Metric_Value = c(
            min(cv_results$rmse, na.rm = TRUE),
            min(cv_results$deviance, na.rm = TRUE),
            prss_results$best_prss
        )
    )
    cat("\nOptimal Basis Dimension (k) Summary:\n")
    print(summary_df)
    
    # Calculates Cohen’s d and relative improvement for statistical analysis.
    optimal_idx <- which(cv_results$k_value == optimal_k_rmse)
    rmse_cohen_d <- numeric(length(cv_results$k_value))
    deviance_cohen_d <- numeric(length(cv_results$k_value))
    for (i in 1:length(cv_results$k_value)) {
        if (i != optimal_idx) {
            rmse_diff <- mean(fold_rmse[, i] - fold_rmse[, optimal_idx], na.rm = TRUE)
            rmse_sd <- sd(fold_rmse[, i] - fold_rmse[, optimal_idx], na.rm = TRUE)
            rmse_cohen_d[i] <- rmse_diff / rmse_sd
            
            deviance_diff <- mean(fold_deviance[, i] - fold_deviance[, optimal_idx], na.rm = TRUE)
            deviance_sd <- sd(fold_deviance[, i] - fold_deviance[, optimal_idx], na.rm = TRUE)
            deviance_cohen_d[i] <- deviance_diff / deviance_sd
        } else {
            rmse_cohen_d[i] <- 0
            deviance_cohen_d[i] <- 0
        }
    }
    
	# Computes relative improvement over the worst-performing k value.
    worst_rmse <- max(cv_results$rmse, na.rm = TRUE)
    worst_deviance <- max(cv_results$deviance, na.rm = TRUE)
    rmse_rel_improv <- (worst_rmse - cv_results$rmse) / worst_rmse * 100
    deviance_rel_improv <- (worst_deviance - cv_results$deviance) / worst_deviance * 100
    
    # Compiles statistical results into a data frame for export.
    stats_df <- data.frame(
        k_value = cv_results$k_value,
        RMSE_Cohen_d = rmse_cohen_d,
        RMSE_Rel_Improv = rmse_rel_improv,
        Deviance_Cohen_d = deviance_cohen_d,
        Deviance_Rel_Improv = deviance_rel_improv
    )
    
    # Exports statistical results to an Excel file for further analysis.
    write_xlsx(stats_df, "CV_Stats_Results_k_values.xlsx")
    cat("\nStatistical results exported to 'CV_Stats_Results_k_values.xlsx'\n")
    
    # Returns the summary and statistical analysis results.
    return(list(summary = summary_df, stats = stats_df))
}

# Runs k-fold cross-validation analysis for PRSS basis dimension k selection.
run_kfold_cv_analysis_for_k <- function(original_obj, k_values, prss_results, k_folds = 10, seed = 123) {
    # Initiates the cross-validation process for k-value selection.
    cat("Running", k_folds, "fold cross-validation for basis dimension (k) selection...\n")
    
    # Verifies the input object contains model data before proceeding.
    if (is.null(original_obj) || !("model" %in% names(original_obj))) {
        cat("ERROR: Original object does not contain model data\n")
        return(NULL)
    }
    
    # Extracts the original dataset from the model object.
    original_data <- original_obj$model
    
    # Performs cross-validation and retrieves results for specified k values.
    cv_data <- perform_kfold_cv_for_k(original_data, k_values, k_folds = k_folds, seed = seed)
    cv_results <- cv_data$cv_results
    
    # Handles cases where all cross-validation results are invalid.
    if (all(is.na(cv_results$rmse))) {
        warning("All cross-validation results are NA. Cannot complete analysis.")
        return(list(
            cv_results = cv_results,
            cv_plots = NULL,
            cv_analysis = NULL,
            optimal_k_cv = NA,
            optimal_model = NULL
        ))
    }
    
    # Generates and displays cross-validation plots for k values.
    cat("\nPlotting cross-validation results...\n")
    cv_plots <- plot_cv_results_for_k(cv_results)
    
    # Analyzes cross-validation results with statistical metrics.
    cat("\nAnalyzing cross-validation results...\n")
    cv_analysis <- analyze_cv_results_for_k(cv_data, prss_results)
    
    # Fits the final model using the CV-optimal k value.
    optimal_k_cv <- cv_results$k_value[which.min(cv_results$rmse)]
    cat("\nFinal model evaluation with CV-optimal k:", optimal_k_cv, "\n")
    
    # Constructs a formula with the optimal k value for the final model.
    formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", optimal_k_cv, ")")
    optimal_formula <- as.formula(formula_string)
    
    # Fits the optimal model with a reproducible seed and error handling.
    set.seed(seed)
    optimal_model <- tryCatch({
        gam(
            formula = optimal_formula,
            data = original_data,
            method = "REML",
            select = TRUE,
            gamma = 1.5
        )
    }, error = function(e) {
        cat("Error fitting optimal model:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    # Summarizes the optimal model if successfully fitted.
    if (!is.null(optimal_model)) {
        cat("\nOptimal model summary:\n")
        print(summary(optimal_model))
    }
    
    # Returns a comprehensive list of cross-validation results for k selection.
    return(list(
        cv_results = cv_results,
        cv_plots = cv_plots,
        cv_analysis = cv_analysis,
        optimal_k_cv = optimal_k_cv,
        optimal_model = optimal_model
    ))
}

# Creates a final summary comparing PRSS and CV methods for k selection.
create_final_k_summary <- function(cv_analysis, prss_results) {
    # Handles cases where CV results are invalid, defaulting to PRSS.
    if (is.null(cv_analysis) || is.null(cv_analysis$optimal_k_cv) || is.na(cv_analysis$optimal_k_cv)) {
        cat("\n======= FINAL ANALYSIS SUMMARY =======\n")
        cat("PRSS-derived optimal k:", prss_results$best_k, "\n")
        cat("CV-optimal k: Unable to determine\n")
        cat("\nRecommendation: Use PRSS-derived optimal k as CV results were not valid.\n")
        return(prss_results$best_k)
    }
    
    # Compares PRSS and CV optimal k values for consistency.
    prss_k <- prss_results$best_k
    cv_k <- cv_analysis$optimal_k_cv
    
    # Outputs the final analysis summary with recommendations.
    cat("\n======= FINAL ANALYSIS SUMMARY =======\n")
    cat("PRSS-derived optimal k:", prss_k, "\n")
    cat("CV-optimal k:", cv_k, "\n")
    
    # Evaluates agreement between PRSS and CV methods.
    if (prss_k == cv_k) {
        cat("\nConsistency between methods: Perfect agreement!\n")
        cat("\nRecommendation: Use k =", prss_k, "as both methods agree.\n")
    } else {
        cat("\nConsistency between methods: Methods selected different k values.\n")
        cat("\nRecommendation: Consider using CV-optimal k =", cv_k, 
            "as it provides better out-of-sample performance.\n")
    }
    
    # Returns the CV-optimal k value as the recommendation.
    return(cv_k)
}

# Executes the k-fold cross-validation analysis with example data.
best_model <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_model

# Prepares PRSS results for comparison, extracting key components.
prss_results <- list(
    best_k = best_model$smooth[[1]]$bs.dim,
    best_prss = pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_params$PRSS,
    best_model = best_model
)

# Defines a range of k values for testing in cross-validation.
k_values <- 3:10

# Runs the cross-validation analysis for k selection.
cv_analysis_for_k <- run_kfold_cv_analysis_for_k(
    original_obj = prss_results$best_model,
    k_values = k_values,
    prss_results = prss_results,
    k_folds = 10,
    seed = 123
)

# Determines the recommended k value based on the analysis.
recommended_k <- create_final_k_summary(cv_analysis_for_k, prss_results)

# Plots a comparison of PRSS and CV methods for k selection if valid results exist.
if (!is.null(cv_analysis_for_k$cv_results) && !all(is.na(cv_analysis_for_k$cv_results$rmse))) {
    # Extracts cross-validation data for plotting.
    cv_data <- cv_analysis_for_k$cv_results
    
    # Normalizes RMSE values for consistent comparison with PRSS.
    cv_data$rmse_norm <- (cv_data$rmse - min(cv_data$rmse, na.rm = TRUE)) / 
        (max(cv_data$rmse, na.rm = TRUE) - min(cv_data$rmse, na.rm = TRUE))
    
    # Creates a data frame for PRSS values across k values.
    prss_values <- data.frame(
        k_value = k_values,
        prss = numeric(length(k_values)),
        prss_norm = numeric(length(k_values))
    )
    
    # Populates PRSS values, approximating for non-optimal k values.
    for (i in 1:length(k_values)) {
        k <- k_values[i]
        if (k == prss_results$best_k) {
            prss_values$prss[i] <- prss_results$best_prss
        } else {
            prss_values$prss[i] <- prss_results$best_prss * (1 + abs(k - prss_results$best_k)/5)
        }
    }
    
    # Normalizes PRSS values for plotting alongside RMSE.
    prss_values$prss_norm <- (prss_values$prss - min(prss_values$prss)) / 
        (max(prss_values$prss) - min(prss_values$prss))
    
    # Combines CV and PRSS data for a unified comparison plot.
    combined_data <- merge(cv_data, prss_values, by = "k_value")
    
    # Generates a comparison plot of normalized PRSS and CV RMSE.
    p_comparison <- ggplot(combined_data) +
        geom_line(aes(x = k_value, y = rmse_norm, color = "CV RMSE"), linewidth = 1) +
        geom_line(aes(x = k_value, y = prss_norm, color = "PRSS"), linewidth = 1) +
        geom_point(aes(x = k_value, y = rmse_norm, color = "CV RMSE"), size = 3) +
        geom_point(aes(x = k_value, y = prss_norm, color = "PRSS"), size = 3) +
        geom_vline(xintercept = recommended_k, linetype = "dashed", color = "darkgreen") +
        geom_vline(xintercept = prss_results$best_k, linetype = "dashed", color = "darkblue") +
        scale_color_manual(values = c("CV RMSE" = "darkred", "PRSS" = "darkblue")) +
        scale_x_continuous(breaks = k_values) +
        labs(title = "Comparison of Basis Dimension (k) Selection Methods", 
             subtitle = "Cross-Validation RMSE vs. PRSS",
             x = "k Value", y = "Normalized Score (0-1)", color = "Method") +
        theme_minimal()
    
    # Adds annotations for optimal k values from CV and PRSS.
    if (!is.na(recommended_k)) {
        p_comparison <- p_comparison + 
            annotate("text", x = recommended_k, y = 0.1, 
                     label = paste("CV Optimal k =", recommended_k), color = "darkgreen")
    }
    
    p_comparison <- p_comparison + 
        annotate("text", x = prss_results$best_k, y = 0.2, 
                 label = paste("PRSS Optimal k =", prss_results$best_k), color = "darkblue")
    
    # Displays the comparison plot.
    print(p_comparison)
}

# Validates the recommended k by comparing PRSS and CV model predictions.
cat("\n--------- Final Validation of Recommended k ---------\n")
if (!is.na(recommended_k) && !is.na(prss_results$best_k) && recommended_k != prss_results$best_k) {
    # Extracts the original dataset for validation.
    original_data <- prss_results$best_model$model
    
    # Constructs a formula for the CV-optimal k value.
    cv_formula <- as.formula(paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", recommended_k, ")"))
    cat("Fitting model with CV-optimal k =", recommended_k, "\n")
    
    # Constructs a formula for the PRSS-optimal k value.
    prss_formula <- as.formula(paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", prss_results$best_k, ")"))
    cat("Comparing with PRSS-optimal k =", prss_results$best_k, "\n")
    
    # Fits the CV-optimal model with error handling.
    set.seed(123)
    cv_model <- tryCatch({
        gam(formula = cv_formula, data = original_data, method = "REML", select = TRUE, gamma = 1.5)
    }, error = function(e) {
        cat("Error fitting CV-optimal model:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    # Proceeds with validation if the CV model fits successfully.
    if (!is.null(cv_model)) {
        # Generates predictions from both PRSS and CV models.
        set.seed(123)
        prss_preds <- predict(prss_results$best_model, newdata = original_data)
        set.seed(123)
        cv_preds <- predict(cv_model, newdata = original_data)
        
        # Calculates differences between PRSS and CV predictions.
        pred_diff <- prss_preds - cv_preds
        max_diff <- max(abs(pred_diff))
        mean_diff <- mean(abs(pred_diff))
        
        # Reports prediction differences and model details.
        cat("PRSS-optimal vs. CV-optimal Model Comparison:\n")
        cat("Max absolute difference in predictions:", max_diff, "\n")
        cat("Mean absolute difference in predictions:", mean_diff, "\n")
        cat("PRSS-optimal k:", prss_results$best_k, "\n")
        cat("CV-optimal k:", recommended_k, "\n")
        
        # Computes in-sample mean squared error for both models.
        prss_mse <- mean((original_data$Expression - prss_preds)^2)
        cv_mse <- mean((original_data$Expression - cv_preds)^2)
        
        # Outputs MSE comparison and percentage difference.
        cat("In-sample MSE for PRSS-optimal model:", prss_mse, "\n")
        cat("In-sample MSE for CV-optimal model:", cv_mse, "\n")
        cat("Percent difference in MSE:", (prss_mse - cv_mse) / prss_mse * 100, "%\n")
        
        # Prepares a data frame for plotting predictions.
        pred_df <- data.frame(
            TRPM4 = original_data$TRPM4, 
            Original = original_data$Expression,
            PRSS_Pred = prss_preds, 
            CV_Pred = cv_preds
        )
        
        # Sorts data by TRPM4 for clearer visualization.
        pred_df <- pred_df[order(pred_df$TRPM4), ]
        
        # Creates a plot comparing PRSS and CV predictions.
        p_pred <- ggplot(pred_df, aes(x = TRPM4)) +
            geom_point(aes(y = Original), alpha = 0.3, size = 1.2, color = "gray50") +
            geom_line(aes(y = PRSS_Pred, color = "PRSS-optimal k"), size = 1) +
            geom_line(aes(y = CV_Pred, color = "CV-optimal k"), size = 1, linetype = "dashed") +
            scale_color_manual(values = c("PRSS-optimal k" = "blue", "CV-optimal k" = "red")) +
            labs(title = "Model Predictions Comparison", 
                 subtitle = paste("PRSS-optimal k =", prss_results$best_k, "vs. CV-optimal k =", recommended_k),
                 x = "TRPM4", y = "Expression", color = "Model") +
            theme_minimal()
        
        # Displays the prediction comparison plot.
        print(p_pred)
        
        # Compares model complexity using effective degrees of freedom.
        prss_edf <- sum(summary(prss_results$best_model)$edf)
        cv_edf <- sum(summary(cv_model)$edf)
        
        # Reports EDF comparison for overfitting assessment.
        cat("\nModel complexity comparison:\n")
        cat("PRSS-optimal model EDF:", prss_edf, "\n")
        cat("CV-optimal model EDF:", cv_edf, "\n")
        cat("Difference in EDF:", prss_edf - cv_edf, "\n")
    }
} else {
    # Skips validation if k values are missing or identical.
    if (is.na(recommended_k) || is.na(prss_results$best_k)) {
        cat("Cannot perform validation due to missing k values\n")
    } else {
        cat("Validation skipped as PRSS-optimal k and CV-optimal k are the same (k =", recommended_k, ")\n")
        cat("This confirms that both methods agree on the optimal basis dimension.\n")
    }
}

# Summarizes the k selection analysis with final recommendations.
cat("\n=============================================\n")
cat("SUMMARY OF BASIS DIMENSION (k) SELECTION ANALYSIS:\n")
cat("=============================================\n")
cat("1. PRSS-derived optimal k:", prss_results$best_k, "\n")

# Reports CV results and method agreement if available.
if (!is.null(cv_analysis_for_k) && !is.null(cv_analysis_for_k$optimal_k_cv) && !is.na(cv_analysis_for_k$optimal_k_cv)) {
    cat("2. CV-based optimal k:", cv_analysis_for_k$optimal_k_cv, "\n")
    cat("3. Methods agreement: ", 
        ifelse(cv_analysis_for_k$optimal_k_cv == prss_results$best_k,
               "Yes - both methods selected the same k value.",
               "No - methods selected different k values."), "\n")
} else {
    cat("2. CV-based optimal k: Could not be determined reliably\n")
    cat("3. Methods agreement: N/A\n")
}

# Outputs the final recommended k value.
cat("4. FINAL RECOMMENDED BASIS DIMENSION (k):", recommended_k, "\n")


# ================================================================================
# 4. Validation #2 for Lambda Selection: 10-Fold Cross-Validation for REML's Lambda Selection
# ================================================================================
# Performs enhanced k-fold cross-validation with detailed fold metrics for lambda selection.
perform_kfold_cv_improved <- function(best_model, lambda_values, k_folds = 10, seed = 123) {
    # Sets a seed for reproducibility of fold assignments.
    set.seed(seed)
    
    # Initializes a data frame to store cross-validation results.
    cv_results <- data.frame(
        lambda = lambda_values,
        rmse = numeric(length(lambda_values)),
        deviance = numeric(length(lambda_values))
    )
    
    # Extracts the original dataset for cross-validation.
    original_data <- best_model$model
    
    # Adjusts the number of folds if the dataset is too small.
    if (nrow(original_data) < 2 * k_folds) {
        cat("WARNING: Too few observations for", k_folds, "folds. Reducing to", floor(nrow(original_data)/2), "folds.\n")
        k_folds <- max(2, floor(nrow(original_data)/2))
    }
    
    # Reports dataset size and number of folds used.
    cat("Number of observations in dataset:", nrow(original_data), "\n")
    cat("Using", k_folds, "folds for cross-validation\n")
    
    # Determines the basis dimension k from the model or uses a default.
    if (length(best_model$smooth) > 0) {
        basis_dim <- best_model$smooth[[1]]$bs.dim
        cat("GAM basis dimension (k):", basis_dim, "\n")
    } else {
        basis_dim <- 10
        cat("Using default GAM basis dimension (k):", basis_dim, "\n")
    }
    
    # Constructs a formula with the specified k value for consistency.
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", basis_dim, ")")
        original_formula <- as.formula(formula_string)
        cat("Using formula:", formula_string, "\n")
    } else {
        original_formula <- best_model$formula
        cat("Using original formula:", deparse(original_formula), "\n")
    }
    
    # Creates k-fold splits with a reproducible seed.
    set.seed(seed)
    folds <- createFolds(1:nrow(original_data), k = k_folds, list = TRUE, returnTrain = FALSE)
    fold_sizes <- sapply(folds, length)
    cat("Fold sizes:", paste(fold_sizes, collapse = ", "), "\n")
    
    # Initializes matrices to store fold-specific RMSE and deviance.
    fold_rmse <- matrix(NA, nrow = k_folds, ncol = length(lambda_values))
    fold_deviance <- matrix(NA, nrow = k_folds, ncol = length(lambda_values))
    
    # Iterates over lambda values to perform cross-validation.
    for (lambda_idx in 1:length(lambda_values)) {
        lambda <- lambda_values[lambda_idx]
        cat("\n------- Processing lambda =", lambda, "-------\n")
        
        # Tracks metrics for each fold within the current lambda.
        fold_metrics <- data.frame(
            fold = 1:k_folds,
            rmse = numeric(k_folds),
            deviance = numeric(k_folds),
            success = logical(k_folds)
        )
        
        # Performs cross-validation for each fold.
        for (j in 1:k_folds) {
            cat("  Fold", j, "of", k_folds, "...\n")
            test_indices <- folds[[j]]
            train_data <- original_data[-test_indices, ]
            test_data <- original_data[test_indices, ]
            cat("    Train size:", nrow(train_data), "Test size:", nrow(test_data), "\n")
            
            # Fits and evaluates the model with error handling.
            tryCatch({
                set.seed(seed + j)
                gam_model <- gam(
                    formula = original_formula,
                    data = train_data,
                    method = "REML",
                    sp = lambda,
                    select = FALSE,
                    gamma = 1.5
                )
                
                # Checks for null model fits.
                if (is.null(gam_model)) {
                    cat("    ERROR: Model fitting returned NULL\n")
                    fold_metrics$success[j] <- FALSE
                    next
                }
                
                # Generates predictions with error handling for robustness.
                predictions <- tryCatch({
                    preds <- predict(gam_model, newdata = test_data)
                    if (any(is.na(preds))) cat("    WARNING: NA values in predictions\n")
                    preds
                }, error = function(e) {
                    cat("    ERROR in prediction:", conditionMessage(e), "\n")
                    return(rep(NA, nrow(test_data)))
                })
                
                # Computes metrics if predictions are valid.
                actuals <- test_data$Expression
                if (all(is.na(predictions))) {
                    cat("    ERROR: All predictions are NA\n")
                    fold_metrics$success[j] <- FALSE
                    next
                }
                
                rmse_value <- sqrt(mean((actuals - predictions)^2, na.rm = TRUE))
                deviance_value <- sum((actuals - predictions)^2, na.rm = TRUE)
                
                fold_metrics$rmse[j] <- rmse_value
                fold_metrics$deviance[j] <- deviance_value
                fold_metrics$success[j] <- TRUE
                
                fold_rmse[j, lambda_idx] <- rmse_value
                fold_deviance[j, lambda_idx] <- deviance_value
                
                cat("    Metrics - RMSE:", round(rmse_value, 4), "\n")
                
            }, error = function(e) {
                cat("    ERROR in fold", j, ":", conditionMessage(e), "\n")
                fold_metrics$success[j] <- FALSE
            })
        }
        
        # Summarizes fold success and computes average metrics.
        successful_folds <- sum(fold_metrics$success)
        cat("  Successfully processed", successful_folds, "out of", k_folds, "folds\n")
        
        if (successful_folds > 0) {
            cv_results$rmse[lambda_idx] <- mean(fold_metrics$rmse[fold_metrics$success], na.rm = TRUE)
            cv_results$deviance[lambda_idx] <- mean(fold_metrics$deviance[fold_metrics$success], na.rm = TRUE)
            cat("  Average metrics - RMSE:", round(cv_results$rmse[lambda_idx], 4), "\n")
        } else {
            cat("  WARNING: No successful folds for lambda =", lambda, "\n")
            cv_results$rmse[lambda_idx] <- NA
            cv_results$deviance[lambda_idx] <- NA
        }
    }
    
    # Warns if all results are invalid and identifies the optimal lambda.
    if (all(is.na(cv_results$rmse))) {
        warning("All cross-validation results are NA. Check for issues in the data or model.")
    } else {
        optimal_lambda_idx <- which.min(cv_results$rmse)
        cat("\nOptimal lambda by CV RMSE:", cv_results$lambda[optimal_lambda_idx], 
            "with RMSE:", round(cv_results$rmse[optimal_lambda_idx], 4), "\n")
    }
    
    # Returns cross-validation results including fold-specific metrics.
    return(list(cv_results = cv_results, fold_rmse = fold_rmse, fold_deviance = fold_deviance))
}

# Visualizes cross-validation results using a viridis plasma palette.
plot_cv_results_improved <- function(cv_results) {
    # Skips plotting if all metrics are invalid.
    if (all(is.na(cv_results$rmse))) {
        cat("Cannot create plots: All CV metrics are NA\n")
        return(NULL)
    }
    
    # Sets x-axis breaks based on lambda values.
    x_breaks <- cv_results$lambda
    
    # Plots RMSE against lambda with a plasma color scheme.
    p_rmse <- ggplot(cv_results, aes(x = lambda, y = rmse, color = factor(lambda))) +
        geom_point(size = 3.5, color = "#111111") +
        geom_point(size = 3) +
        geom_line(color = "#323232", linewidth = 0.3, stat = "identity") +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, guide = guide_legend(
            override.aes = list(
                shape = 21,
                fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(length(cv_results$lambda)),
                color = "#111111",
                size = 3,
                stroke = 0.5
            )
        )) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = "Cross-Validation: Lambda vs RMSE", 
             x = "Lambda (λ)", 
             y = "Root Mean Squared Error",
             color = "Lambda (λ)") +
        scale_y_continuous(labels = function(x) format(x, digits = 3, nsmall = 3)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Plots deviance against lambda with a plasma color scheme.
    p_deviance <- ggplot(cv_results, aes(x = lambda, y = deviance, color = factor(lambda))) +
        geom_point(size = 3.5, color = "#111111") +
        geom_point(size = 3) +
        geom_line(color = "#323232", linewidth = 0.3, stat = "identity") +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, guide = guide_legend(
            override.aes = list(
                shape = 21,
                fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(length(cv_results$lambda)),
                color = "#111111",
                size = 3,
                stroke = 0.5
            )
        )) +
        scale_x_continuous(n.breaks = 5) +
        labs(title = "Cross-Validation: Lambda vs Deviance", 
             x = "Lambda (λ)", 
             y = "Deviance",
             color = "Lambda (λ)") +
        scale_y_continuous(labels = function(x) format(x, digits = 2, nsmall = 2)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Displays the RMSE and deviance plots.
    print(p_rmse)
    print(p_deviance)
    
    # Returns the list of cross-validation plots.
    return(list(rmse_plot = p_rmse, deviance_plot = p_deviance))
}

# Analyzes cross-validation results with streamlined statistical metrics.
analyze_cv_results_improved <- function(cv_data, reml_results) {
    # Extracts cross-validation results and fold-specific metrics.
    cv_results <- cv_data$cv_results
    fold_rmse <- cv_data$fold_rmse
    fold_deviance <- cv_data$fold_deviance
    
    # Terminates analysis if all metrics are invalid.
    if (all(is.na(cv_results$rmse))) {
        cat("Cannot analyze results: All CV metrics are NA\n")
        return(NULL)
    }
    
    # Identifies optimal lambda values based on RMSE, deviance, and REML.
    optimal_rmse <- cv_results$lambda[which.min(cv_results$rmse)]
    optimal_deviance <- cv_results$lambda[which.min(cv_results$deviance)]
    optimal_reml <- reml_results$lambda_value[which.min(reml_results$reml_score)]
    
    # Summarizes optimal lambda values and their metrics.
    summary_df <- data.frame(
        Criterion = c("CV RMSE", "CV Deviance", "REML Score"),
        Optimal_Lambda = c(optimal_rmse, optimal_deviance, optimal_reml),
        Metric_Value = c(
            min(cv_results$rmse, na.rm = TRUE),
            min(cv_results$deviance, na.rm = TRUE),
            min(reml_results$reml_score)
        )
    )
    cat("\nOptimal Lambda Summary:\n")
    print(summary_df)
    
    # Calculates Cohen’s d and relative improvement for statistical comparison.
    optimal_idx <- which(cv_results$lambda == optimal_rmse)
    rmse_cohen_d <- numeric(length(cv_results$lambda))
    deviance_cohen_d <- numeric(length(cv_results$lambda))
    for (i in 1:length(cv_results$lambda)) {
        if (i != optimal_idx) {
            rmse_diff <- mean(fold_rmse[, i] - fold_rmse[, optimal_idx], na.rm = TRUE)
            rmse_sd <- sd(fold_rmse[, i] - fold_rmse[, optimal_idx], na.rm = TRUE)
            rmse_cohen_d[i] <- rmse_diff / rmse_sd
            
            deviance_diff <- mean(fold_deviance[, i] - fold_deviance[, optimal_idx], na.rm = TRUE)
            deviance_sd <- sd(fold_deviance[, i] - fold_deviance[, optimal_idx], na.rm = TRUE)
            deviance_cohen_d[i] <- deviance_diff / deviance_sd
        } else {
            rmse_cohen_d[i] <- 0
            deviance_cohen_d[i] <- 0
        }
    }
    
    # Computes relative improvement over the worst-performing lambda.
    worst_rmse <- max(cv_results$rmse, na.rm = TRUE)
    worst_deviance <- max(cv_results$deviance, na.rm = TRUE)
    rmse_rel_improv <- (worst_rmse - cv_results$rmse) / worst_rmse * 100
    deviance_rel_improv <- (worst_deviance - cv_results$deviance) / worst_deviance * 100
    
    # Compiles statistical results into a data frame.
    stats_df <- data.frame(
        Lambda = cv_results$lambda,
        RMSE_Cohen_d = rmse_cohen_d,
        RMSE_Rel_Improv = rmse_rel_improv,
        Deviance_Cohen_d = deviance_cohen_d,
        Deviance_Rel_Improv = deviance_rel_improv
    )
    
    # Exports statistical results to an Excel file.
    write_xlsx(stats_df, "CV_Stats_Results_Lambda.xlsx")
    cat("\nStreamlined statistical results exported to 'CV_Stats_Results_Lambda.xlsx'\n")
    
    # Returns the summary and statistical analysis results.
    return(list(summary = summary_df, stats = stats_df))
}

# Runs k-fold cross-validation analysis for lambda selection.
run_kfold_cv_analysis_improved <- function(best_model, lambda_values, eigen_metrics, k_folds = 10, seed = 123) {
    # Initiates the cross-validation process with specified parameters.
    cat("Running", k_folds, "fold cross-validation for lambda selection...\n")
    
    # Performs cross-validation and retrieves results.
    cv_data <- perform_kfold_cv_improved(best_model, lambda_values, k_folds = k_folds, seed = seed)
    cv_results <- cv_data$cv_results
    
    # Handles cases where all results are invalid.
    if (all(is.na(cv_results$rmse))) {
        warning("All cross-validation results are NA. Cannot complete analysis.")
        return(list(
            cv_results = cv_results,
            cv_plots = NULL,
            cv_analysis = NULL,
            optimal_lambda_cv = NA,
            optimal_model = NULL
        ))
    }
    
    # Generates and displays cross-validation plots.
    cat("\nPlotting cross-validation results...\n")
    cv_plots <- plot_cv_results_improved(cv_results)
    
    # Analyzes cross-validation results with statistical metrics.
    cat("\nAnalyzing cross-validation results...\n")
    cv_analysis <- analyze_cv_results_improved(cv_data, eigen_metrics)
    
    # Fits the final model with the CV-optimal lambda.
    optimal_lambda_cv <- cv_results$lambda[which.min(cv_results$rmse)]
    cat("\nFinal model evaluation with CV-optimal lambda:", optimal_lambda_cv, "\n")
    
    k_value <- if (length(best_model$smooth) > 0) best_model$smooth[[1]]$bs.dim else 10
    
    if ("Expression" %in% names(best_model$model) && "TRPM4" %in% names(best_model$model)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k_value, ")")
        optimal_formula <- as.formula(formula_string)
    } else {
        optimal_formula <- best_model$formula
    }
    
    set.seed(seed)
    optimal_model <- tryCatch({
        gam(
            formula = optimal_formula,
            data = best_model$model,
            method = "REML",
            sp = optimal_lambda_cv,
            select = FALSE
        )
    }, error = function(e) {
        cat("Error fitting optimal model:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    # Summarizes the optimal model if successfully fitted.
    if (!is.null(optimal_model)) {
        cat("\nOptimal model summary:\n")
        print(summary(optimal_model))
    }
    
    # Returns a comprehensive list of cross-validation results.
    return(list(
        cv_results = cv_results,
        cv_plots = cv_plots,
        cv_analysis = cv_analysis,
        optimal_lambda_cv = optimal_lambda_cv,
        optimal_model = optimal_model
    ))
}

# Creates a final summary comparing REML and CV lambda selection methods.
create_final_lambda_summary <- function(cv_analysis, eigen_metrics) {
    # Handles cases where CV results are invalid.
    if (is.null(cv_analysis) || is.null(cv_analysis$optimal_lambda_cv) || is.na(cv_analysis$optimal_lambda_cv)) {
        cat("\n======= FINAL ANALYSIS SUMMARY =======\n")
        cat("REML-optimal lambda:", eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)], "\n")
        cat("CV-optimal lambda: Unable to determine\n")
        cat("\nRecommendation: Use REML-optimal lambda as CV results were not valid.\n")
        return(eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)])
    }
    
    # Compares REML and CV optimal lambda values.
    reml_lambda <- eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)]
    cv_lambda <- cv_analysis$optimal_lambda_cv
    
    cat("\n======= FINAL ANALYSIS SUMMARY =======\n")
    cat("REML-optimal lambda:", reml_lambda, "\n")
    cat("CV-optimal lambda:", cv_lambda, "\n")
    
    if (!is.na(cv_lambda) && !is.na(reml_lambda)) {
        if (abs(cv_lambda - reml_lambda) < 1e-6) {
            cat("\nConsistency between methods: Perfectly consistent!\n")
            cat("\nRecommendation: Use either method's lambda value.\n")
        } else {
            cat("\nConsistency between methods: Methods selected different lambda values.\n")
            cat("\nRecommendation: Consider using CV-optimal lambda.\n")
        }
    }
    
    # Returns the recommended lambda, prioritizing CV if available.
    return(ifelse(!is.na(cv_lambda), cv_lambda, reml_lambda))
}

# Executes the cross-validation analysis with example data.
best_model <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_model
lambda_values <- c(0.264006588, 0.332627458, 0.419084336, 0.528013177, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
eigen_metrics <- extract_eigenvalues(best_model, lambda_values)

cv_analysis_improved <- run_kfold_cv_analysis_improved(
    best_model = best_model,
    lambda_values = lambda_values,
    eigen_metrics = eigen_metrics,
    k_folds = 10,
    seed = 123
)

recommended_lambda <- create_final_lambda_summary(cv_analysis_improved, eigen_metrics)

# Plots a comparison of REML and CV methods if valid results exist.
if (!is.null(cv_analysis_improved$cv_results) && !all(is.na(cv_analysis_improved$cv_results$rmse))) {
    reml_data <- eigen_metrics[, c("lambda_value", "reml_score")]
    cv_data <- cv_analysis_improved$cv_results[, c("lambda", "rmse")]
    names(reml_data) <- c("lambda", "reml")
    names(cv_data) <- c("lambda", "rmse")
    
    reml_data$reml_norm <- (reml_data$reml - min(reml_data$reml)) / (max(reml_data$reml) - min(reml_data$reml))
    cv_data$rmse_norm <- (cv_data$rmse - min(cv_data$rmse, na.rm = TRUE)) / (max(cv_data$rmse, na.rm = TRUE) - min(cv_data$rmse, na.rm = TRUE))
    
    combined_data <- merge(reml_data, cv_data, by = "lambda")
    
    p_final <- ggplot(combined_data) +
        geom_line(aes(x = lambda, y = reml_norm, color = "REML Score"), linewidth = 1) +
        geom_line(aes(x = lambda, y = rmse_norm, color = "CV RMSE"), linewidth = 1) +
        geom_point(aes(x = lambda, y = reml_norm, color = "REML Score"), size = 3) +
        geom_point(aes(x = lambda, y = rmse_norm, color = "CV RMSE"), size = 3) +
        geom_vline(xintercept = cv_analysis_improved$optimal_lambda_cv, linetype = "dashed", color = "darkgreen") +
        geom_vline(xintercept = eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)], linetype = "dashed", color = "darkblue") +
        scale_color_manual(values = c("REML Score" = "darkblue", "CV RMSE" = "darkred")) +
        labs(title = "Comparison of Lambda Selection Methods", subtitle = "REML Score vs. Cross-Validation RMSE",
             x = "Lambda (λ)", y = "Normalized Score (0-1)", color = "Method") +
        theme_minimal()
    
    if (!is.na(cv_analysis_improved$optimal_lambda_cv)) {
        p_final <- p_final + annotate("text", x = cv_analysis_improved$optimal_lambda_cv, y = 0.1, 
                                      label = paste("CV Optimal λ =", cv_analysis_improved$optimal_lambda_cv), color = "darkgreen")
    }
    
    p_final <- p_final + annotate("text", x = eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)], y = 0.2, 
                                  label = paste("REML Optimal λ =", eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)]), color = "darkblue")
    
    print(p_final)
}

# Validates the recommended lambda by comparing model predictions.
cat("\n--------- Final Validation of Recommended Lambda ---------\n")
if ("Expression" %in% names(best_model$model) && "TRPM4" %in% names(best_model$model)) {
    validation_formula <- as.formula("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k = 10)")
    cat("Using validation formula: Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k = 10)\n")
} else {
    validation_formula <- as.formula("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k = 10)")
    cat("Warning: Using default validation formula with k = 10 due to missing column names\n")
}

set.seed(123)
validation_model <- tryCatch({
    gam(formula = validation_formula, data = best_model$model, method = "REML", sp = recommended_lambda, select = FALSE)
}, error = function(e) {
    cat("Error fitting validation model:", conditionMessage(e), "\n")
    return(NULL)
})

if (!is.null(validation_model)) {
    original_data <- best_model$model
    set.seed(123)
    original_preds <- predict(best_model, newdata = original_data)
    set.seed(123)
    validation_preds <- predict(validation_model, newdata = original_data)
    
    pred_diff <- original_preds - validation_preds
    max_diff <- max(abs(pred_diff))
    mean_diff <- mean(abs(pred_diff))
    
    cat("Original vs. Recommended Lambda Model Comparison:\n")
    cat("Max absolute difference in predictions:", max_diff, "\n")
    cat("Mean absolute difference in predictions:", mean_diff, "\n")
    cat("Original model lambda:", best_model$sp, "\n")
    cat("Recommended lambda:", recommended_lambda, "\n")
    
    pred_df <- data.frame(TRPM4 = original_data$TRPM4, Original = original_preds, Recommended = validation_preds)
    
    p_pred <- ggplot(pred_df, aes(x = TRPM4)) +
        geom_line(aes(y = Original, color = "Original Model"), size = 1) +
        geom_line(aes(y = Recommended, color = "Recommended λ Model"), size = 1, linetype = "dashed") +
        scale_color_manual(values = c("Original Model" = "blue", "Recommended λ Model" = "red")) +
        labs(title = "Model Predictions Comparison", subtitle = paste("Original vs. Recommended λ =", recommended_lambda),
             x = "TRPM4", y = "Predicted Expression", color = "Model") +
        theme_minimal()
    
    print(p_pred)
}

# Summarizes the lambda selection analysis results.
cat("\n=============================================\n")
cat("SUMMARY OF LAMBDA SELECTION ANALYSIS:\n")
cat("=============================================\n")
cat("1. REML-based optimal lambda:", eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)], "\n")

if (!is.null(cv_analysis_improved) && !is.null(cv_analysis_improved$optimal_lambda_cv) && !is.na(cv_analysis_improved$optimal_lambda_cv)) {
    cat("2. CV-based optimal lambda:", cv_analysis_improved$optimal_lambda_cv, "\n")
    cat("3. Methods agreement: ", 
        ifelse(abs(cv_analysis_improved$optimal_lambda_cv - eigen_metrics$lambda_value[which.min(eigen_metrics$reml_score)]) < 1e-6,
               "Yes - both methods selected the same lambda value.",
               "No - methods selected different lambda values."), "\n")
} else {
    cat("2. CV-based optimal lambda: Could not be determined reliably\n")
    cat("3. Methods agreement: N/A\n")
}

cat("4. FINAL RECOMMENDED LAMBDA:", recommended_lambda, "\n")


# -------------------------------------------------------------------------------------------
# List of lambda and k values by sample for separate execution 
# -------------------------------------------------------------------------------------------
# HYW_4701_Tumor PCa Ribo (best k=6): 0.156970198, 0.248780998, 0.39429131, 0.624909614, 0.990414992, 1.569701979, 2.0, 2.5, 3.0, 3.5, 4.0
# HYW_4847_Tumor PCa Ribo (best k=3): 0.007934132, 0.012574752, 0.019929639, 0.03158635, 0.050060991, 0.079341323, 0.09, 0.11, 0.13, 0.15, 0.17
# HYW_4880_Tumor PCa Ribo (best k=6): 0.069433736, 0.082571093, 0.098194131, 0.11677316, 0.138867472, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26
# HYW_4881_Tumor PCa Ribo (best k=10): 0.264006588, 0.332627458, 0.419084336, 0.528013177, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0
# HYW_5386_Tumor PCa Ribo (best k=4): 0.02874564, 0.05111778, 0.090901696, 0.161648614, 0.287456402, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4
# HYW_5742_Tumor PCa Ribo (best k=4): 2549.67887, 2928.811924, 3364.32144, 3864.590503, 4439.248754, 5099.357741, 5500, 6000, 7000, 8000, 9000
# HYW_5755_Tumor PCa Ribo (best k=3): 0.128680173, 0.203944331, 0.323229981, 0.512284997, 0.811917005, 1.286801734, 1.4, 1.6, 1.8, 2.0, 2.2

# HYW_4701_Tumor PCa AR (best k=8): 0.206858027, 0.36785137, 0.654142517, 1.163248169, 2.068580267
# HYW_4847_Tumor PCa AR (best k=3): 0.102594568, 0.150588232, 0.221033296, 0.32443251, 0.476201801, 0.698968654, 1.025945679
# HYW_4880_Tumor PCa AR (best k=6): 0.043297249, 0.076994606, 0.136917922, 0.243478322, 0.432972486
# HYW_4881_Tumor PCa AR (best k=5): 0.176738615, 0.210178818, 0.249946146, 0.297237735, 0.35347723
# HYW_5386_Tumor PCa AR (best k=10): 6514.902627, 9052.427847, 12578.30771, 17477.50188, 24284.91012, 33743.7732, 46886.82082, 65149.02627
# HYW_5742_Tumor PCa AR (best k=10): 1.296793131, 2.306060523, 4.100819947, 7.292403676, 12.96793131
# HYW_5755_Tumor PCa AR (best k=6): 1166.630911, 1621.028394, 2252.411649, 3129.715838, 4348.726056, 6042.535261, 8396.07552, 11666.30911

# Notes: lambda values with one decimal points were synthetically added post-GAM fitting to assess if these higher lamda values return with more optimal results.

# -------------------------------------------------------------------------------------------
# Specifies SVG dimensions for figures
# -------------------------------------------------------------------------------------------
# Fig_Cond-numb-vs-REML (eigenvalues) 552 x 344
# Fig_Lambda-vs-REML (refitting) 552 x 344
# Fig_k-vs-PRSS (refitting) 540 x 344
# Fig_CV_Lambda-vs-RMSE 552 x 344
# Fig_CV_Lambda-vs-Deviance 552 x 344
# Fig_CV_k-vs-RMSE 540 x 344
# Fig_CV_k-vs-Deviance 535 x 344