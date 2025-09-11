############################################################################################
# Part 3.07: Validation of Selected k (by PRSS) and Lambda (by REML) Values by Refitting with Fixed k or Lambda
############################################################################################
library(mgcv)
library(ggplot2)
library(scales)

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
# Validation #2: Refitting with Different Lambda Values (REML Score Analysis)
# =================================================================================
# Function to extract REML scores for plotting
extract_metrics <- function(best_model, lambda_values) {
    # Initialize a data frame to store REML scores for each lambda
    results <- data.frame(
        lambda_value = lambda_values,
        reml_score = numeric(length(lambda_values))
    )
    
    # Preserve the original dataset for consistent refitting
    original_data <- best_model$model
    
    # Retrieve the basis dimension k from the smooth term or set a default value
    if (length(best_model$smooth) > 0) {
        k_value <- best_model$smooth[[1]]$bs.dim
    } else {
        k_value <- 10
    }
    
    # Construct a formula using the specified k value and dataset column names
    if ("Expression" %in% names(original_data) && "TRPM4" %in% names(original_data)) {
        formula_string <- paste("Expression ~ TRPM4 + s(TRPM4, bs = 'tp', k =", k_value, ")")
        original_formula <- as.formula(formula_string)
    } else {
        original_formula <- best_model$formula
    }
    
    # Iterate over lambda values to extract REML scores from refitted models
    for (i in seq_along(lambda_values)) {
        lambda <- lambda_values[i]
        
        # Refit the GAM model with a fixed lambda using REML optimization
        refitted_model <- gam(
            formula = original_formula,
            data = original_data,
            method = "REML",
            sp = lambda,
            select = FALSE,
            gamma = 1.5
        )
        
        # Record the REML score
        results$reml_score[i] <- refitted_model$gcv.ubre
    }
    
    return(results)
}

# Create the Lambda vs REML Score plot
create_lambda_reml_plot <- function(metrics) {
    # Plot lambda against REML score using a plasma palette
    p <- ggplot(metrics, aes(x = lambda_value, y = reml_score, color = factor(lambda_value))) +
        geom_point(size = 4.5, color = "#111111") +
        geom_point(size = 4) +
        geom_line(color = "#323232", linewidth = 0.3) +
        scale_color_viridis_d(option = "plasma", begin = 0, end = 1, guide = guide_legend(
            override.aes = list(
                shape = 21,
                fill = scales::viridis_pal(option = "plasma", begin = 0, end = 1)(length(unique(metrics$lambda_value))),
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
    
    # Display the plot
    print(p)
    
    # Return the plot
    return(p)
}

# Analyze REML scores across lambda values and generate the Lambda vs REML Score plot
analyze_for_lambdas <- function(best_model, lambda_values) {
    # Extract REML scores for the specified lambda values
    metrics <- extract_metrics(best_model, lambda_values)
    
    # Generate the Lambda vs REML Score plot
    lambda_reml_plot <- create_lambda_reml_plot(metrics)
    
    # Output the metrics to the console
    cat("\nREML Scores for Lambda Values:\n")
    print(metrics)
    
    # Identify the optimal lambda based on the minimum REML score
    min_reml_idx <- which.min(metrics$reml_score)
    optimal_lambda <- metrics$lambda_value[min_reml_idx]
    
    # Report the optimal lambda and its associated REML score
    cat("\nOptimal lambda by REML criterion:", optimal_lambda, "\n")
    cat("REML score for optimal lambda:\n")
    print(metrics[min_reml_idx, ])
    
    # Return the analysis results
    return(list(
        metrics = metrics,
        lambda_reml_plot = lambda_reml_plot,
        optimal_lambda = optimal_lambda,
        optimal_metrics = metrics[min_reml_idx, ]
    ))
}

# Execute the analysis with example data
best_model <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$best_model
lambda_values <- c(0.264006588, 0.332627458, 0.419084336, 0.528013177, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
results <- analyze_for_lambdas(best_model, lambda_values)


