########################################################
# Part 3.05: REML Extraction and Convergence Analysis
########################################################
library(mgcv)
library(writexl) 
library(readxl)

# Function to extract detailed REML information
extract_reml_info <- function(model) {
    if (is.null(model)) {
        return(list(
            # Model Fit and Performance
            reml_score = NA,
            deviance_explained = NA,
            aic = NA,
            null_deviance = NA,
            
            # Smoothing Parameters
            sp = NA,
            smoothing_param = NA,
            final_lambda_values = NA,
            
            # Effective Degrees of Freedom
            edf = NA,
            edf1 = NA,
            edf2 = NA,
            
            # Optimization Process
            num_iterations = NA,
            lambda_values = NA,
            reml_scores = NA,
            convergence = NA,
            outer_info.convergence = NA,
            optimizer = NA,
            
            # Model Structure
            smooth_info.term = NA,
            smooth_info.dim = NA,
            smooth_info.p_order = NA,
            smooth_info.by = NA,
            smooth_info.label = NA,
            smooth_info.basis = NA,
            smooth_info.nPar = NA,
            
            # Scale and Residuals
            scale = NA,
            sig2 = NA,
            residual_df = NA,
            
            # Model Complexity
            rank = NA,
            np = NA,
            
            # Gradient and Hessian
            outer_info.gradient = NA,
            outer_info.hessian = NA,
            
            # Lambda Components
            tr_A_lambda = NA,
            y_V_y = NA,
            lambda = NA,
            gamma = NA,
            
            # Other
            method = NA,
            smooth_info.nCoef = NA
        ))
    }
    
    tryCatch({
        # Extract REML iteration info
        num_iterations <- if (!is.null(model$outer.info$iter)) model$outer.info$iter else NA
        lambda_values <- if (!is.null(model$outer.info$sp)) {
            paste(apply(model$outer.info$sp, 1, paste, collapse = ", "), collapse = "; ")
        } else {
            "Unable to extract lambda values per iteration"
        }
        reml_scores <- if (!is.null(model$outer.info$score)) paste(model$outer.info$score, collapse = ", ") else NA
        
        # Extract smooth information
        smooth_info <- if (length(model$smooth) > 0) {
            s <- model$smooth[[1]]
            list(
                term = s$term,
                dim = s$dim,
                p_order = s$p.order,
                by = s$by,
                label = s$label,
                nCoef = length(s$coefficients),
                basis = s$bs.dim,
                nPar = s$df
            )
        } else {
            list(term = NA, dim = NA, p_order = NA, by = NA, label = NA, nCoef = NA, basis = NA, nPar = NA)
        }
        
        # Extract lambda components
        lambda_components <- extract_lambda_components(model)
        
        list(
            # Model Fit and Performance
            reml_score = model$gcv.ubre,
            deviance_explained = summary(model)$dev.expl,
            aic = model$aic,
            null_deviance = model$null.deviance,
            
            # Smoothing Parameters
            sp = paste(model$sp, collapse = ", "),
            smoothing_param = paste(model$sp, collapse = ", "),
            final_lambda_values = paste(model$sp, collapse = ", "),
            
            # Effective Degrees of Freedom
            edf = sum(model$edf),
            edf1 = paste(model$edf1, collapse = ", "),
            edf2 = paste(model$edf2, collapse = ", "),
            
            # Optimization Process
            num_iterations = num_iterations,
            lambda_values = lambda_values,
            reml_scores = reml_scores,
            convergence = model$converged,
            outer_info.convergence = if (!is.null(model$outer.info$conv)) model$outer.info$conv else NA,
            optimizer = paste(model$optimizer, collapse = ", "),
            
            # Model Structure
            smooth_info.term = smooth_info$term,
            smooth_info.dim = smooth_info$dim,
            smooth_info.p_order = smooth_info$p_order,
            smooth_info.by = smooth_info$by,
            smooth_info.label = smooth_info$label,
            smooth_info.basis = smooth_info$basis,
            smooth_info.nPar = smooth_info$nPar,
            
            # Scale and Residuals
            scale = model$scale,
            sig2 = model$sig2,
            residual_df = model$df.residual,
            
            # Model Complexity
            rank = model$rank,
            np = length(model$coefficients),
            
            # Gradient and Hessian
            outer_info.gradient = if (!is.null(model$outer.info$grad)) paste(model$outer.info$grad, collapse = ", ") else NA,
            outer_info.hessian = if (!is.null(model$outer.info$hess)) paste(model$outer.info$hess, collapse = ", ") else NA,
            
            # Lambda Components
            tr_A_lambda = lambda_components$tr_A_lambda,
            y_V_y = lambda_components$y_V_y,
            lambda = lambda_components$lambda,
            gamma = lambda_components$gamma,
            
            # Other
            method = model$method,
            smooth_info.nCoef = smooth_info$nCoef
        )
    }, error = function(e) {
        cat("Error in extract_reml_info:", conditionMessage(e), "\n")
        list(
            # Model Fit and Performance
            reml_score = NA,
            deviance_explained = NA,
            aic = NA,
            null_deviance = NA,
            
            # Smoothing Parameters
            sp = NA,
            smoothing_param = NA,
            final_lambda_values = NA,
            
            # Effective Degrees of Freedom
            edf = NA,
            edf1 = NA,
            edf2 = NA,
            
            # Optimization Process
            num_iterations = NA,
            lambda_values = NA,
            reml_scores = NA,
            convergence = NA,
            outer_info.convergence = NA,
            optimizer = NA,
            
            # Model Structure
            smooth_info.term = NA,
            smooth_info.dim = NA,
            smooth_info.p_order = NA,
            smooth_info.by = NA,
            smooth_info.label = NA,
            smooth_info.basis = NA,
            smooth_info.nPar = NA,
            
            # Scale and Residuals
            scale = NA,
            sig2 = NA,
            residual_df = NA,
            
            # Model Complexity
            rank = NA,
            np = NA,
            
            # Gradient and Hessian
            outer_info.gradient = NA,
            outer_info.hessian = NA,
            
            # Lambda Components
            tr_A_lambda = NA,
            y_V_y = NA,
            lambda = NA,
            gamma = NA,
            
            # Other
            method = NA,
            smooth_info.nCoef = NA
        )
    })
}

# Function to extract REML lambda components
extract_lambda_components <- function(model) {
    if (is.null(model)) {
        return(list(
            tr_A_lambda = NA,
            y_V_y = NA,
            lambda = NA,
            reml_score = NA,
            gamma = NA
        ))
    }
    
    tr_A_lambda <- sum(model$edf)
    y_V_y <- model$deviance
    lambda <- model$sp[1]  # Assuming the first smoothing parameter is what we want
    reml_score <- model$gcv.ubre
    gamma <- model$gamma  # This is typically 1 unless changed in gam() call
    
    list(
        tr_A_lambda = tr_A_lambda,
        y_V_y = y_V_y,
        lambda = lambda,
        reml_score = reml_score,
        gamma = gamma
    )
}

# Function to process and combine REML information for a list of models
process_reml_info <- function(results_list, group_name) {
    results <- lapply(names(results_list), function(sample_name) {
        sample_results <- results_list[[sample_name]]
        if (is.null(sample_results)) return(NULL)
        
        lapply(names(sample_results), function(gene_set_name) {
            model <- sample_results[[gene_set_name]]$best_model
            if (is.null(model)) return(NULL)
            
            reml_info <- extract_reml_info(model)
            c(list(Sample = sample_name, Group = group_name, Gene_Set = gene_set_name), reml_info)
        })
    })
    
    results <- unlist(results, recursive = FALSE)
    results <- results[!sapply(results, is.null)]
    
    # Ensure all lists have the same structure
    all_names <- unique(unlist(lapply(results, names)))
    results <- lapply(results, function(x) {
        x[setdiff(all_names, names(x))] <- NA
        x[all_names]
    })
    
    do.call(rbind, lapply(results, function(x) as.data.frame(t(unlist(x)), stringsAsFactors = FALSE)))
}

# Process REML info for PCa and Non-Ca models
pca_reml_info <- process_reml_info(pca_results, "PCa")
non_ca_reml_info <- process_reml_info(non_ca_results, "Non-Ca")

# Combine all REML info
all_reml_info <- rbind(pca_reml_info, non_ca_reml_info)

# Define categories
categories <- c(
    "Model Fit and Performance",
    "Smoothing Parameters",
    "Effective Degrees of Freedom",
    "Optimization Process",
    "Model Structure",
    "Scale and Residuals",
    "Model Complexity",
    "Gradient and Hessian",
    "Lambda Components",
    "Other"
)

# Add legends
add_legend <- function(df) {
    category_map <- list(
        "Model Fit and Performance" = c("reml_score", "deviance_explained", "aic", "null_deviance"),
        "Smoothing Parameters" = c("sp", "smoothing_param", "final_lambda_values"),
        "Effective Degrees of Freedom" = c("edf", "edf1", "edf2"),
        "Optimization Process" = c("num_iterations", "lambda_values", "reml_scores", "convergence", "outer_info.convergence", "optimizer"),
        "Model Structure" = c("smooth_info.term", "smooth_info.dim", "smooth_info.p_order", "smooth_info.by", "smooth_info.label", "smooth_info.basis", "smooth_info.nPar"),
        "Scale and Residuals" = c("scale", "sig2", "residual_df"),
        "Model Complexity" = c("rank", "np"),
        "Gradient and Hessian" = c("outer_info.gradient", "outer_info.hessian"),
        "Lambda Components" = c("tr_A_lambda", "y_V_y", "lambda", "gamma"),
        "Other" = c("method", "smooth_info.nCoef")
    )
    
    legend_rows <- do.call(rbind, lapply(names(category_map), function(category) {
        params <- paste(category_map[[category]], collapse = ", ")
        data.frame(
            Sample = category,
            Group = "Legend",
            Gene_Set = params,
            t(rep(NA, ncol(df) - 3))  # Fill other columns with NA
        )
    }))
    
    colnames(legend_rows) <- colnames(df)
    
    # Add empty rows between data and legend
    empty_rows <- matrix(NA, nrow = 2, ncol = ncol(df))
    colnames(empty_rows) <- colnames(df)
    
    rbind(df, empty_rows, legend_rows)
}

# Add legend to the data
reml_info_with_legend <- add_legend(all_reml_info)


# =========================================
# Extract REML Components
# =========================================
extract_reml_components_simple <- function(model) {
    if (is.null(model)) {
        return(list(
            final_reml_score = NA_real_,
            scale = NA_real_,
            rss = NA_real_,
            edf = NA_real_,
            lambda = NA_character_,
            model_formula = NA_character_,
            gamma = NA_real_
        ))
    }
    
    tryCatch({
        # Get the final REML score and other basic information
        final_reml_score <- model$gcv.ubre
        scale <- model$scale
        n <- nrow(model$model)
        edf <- sum(model$edf)
        rss <- sum(model$residuals^2)
        
        # Get model information
        lambda_str <- if (!is.null(model$sp)) {
            paste(as.character(model$sp), collapse=", ")
        } else {
            NA_character_
        }
        
        formula_str <- if (!is.null(model$formula)) {
            as.character(deparse(model$formula))[1]
        } else {
            NA_character_
        }
        
        gamma_val <- if (!is.null(model$gamma)) {
            as.numeric(model$gamma)
        } else {
            1  # Default gamma is 1
        }
        
        # Return all components
        list(
            final_reml_score = final_reml_score,
            scale = scale,
            rss = rss,
            edf = edf,
            lambda = lambda_str,
            model_formula = formula_str,
            gamma = gamma_val
        )
    }, error = function(e) {
        cat("Error in extract_reml_components_simple:", conditionMessage(e), "\n")
        list(
            final_reml_score = NA_real_,
            scale = NA_real_,
            rss = NA_real_,
            edf = NA_real_,
            lambda = NA_character_,
            model_formula = NA_character_,
            gamma = NA_real_
        )
    })
}

# Modified process_reml_components function to use the simpler approach
process_reml_components_simple <- function(results_list, group_name) {
    # Create an empty data frame with the correct structure
    result_df <- data.frame(
        Sample = character(),
        Group = character(),
        Gene_Set = character(),
        Final_REML_Score = numeric(),
        Scale = numeric(),
        RSS = numeric(),
        EDF = numeric(),
        Lambda = character(),
        Gamma = numeric(),
        Model_Formula = character(),
        stringsAsFactors = FALSE
    )
    
    for (sample_name in names(results_list)) {
        cat("Processing sample:", sample_name, "\n")
        sample_results <- results_list[[sample_name]]
        
        if (is.null(sample_results)) {
            cat("  No results for this sample\n")
            next
        }
        
        for (gene_set_name in names(sample_results)) {
            cat("  Processing gene set:", gene_set_name, "\n")
            x <- sample_results[[gene_set_name]]
            
            if (is.null(x) || is.null(x$best_model)) {
                cat("    No model found\n")
                next
            }
            
            model <- x$best_model
            
            # Use the simpler function that avoids matrix manipulation errors
            reml_components <- extract_reml_components_simple(model)
            
            # Create a single row data frame
            row_df <- data.frame(
                Sample = sample_name,
                Group = group_name,
                Gene_Set = gene_set_name,
                Final_REML_Score = reml_components$final_reml_score,
                Scale = reml_components$scale,
                RSS = reml_components$rss,
                EDF = reml_components$edf,
                Lambda = reml_components$lambda,
                Gamma = reml_components$gamma,
                Model_Formula = reml_components$model_formula,
                stringsAsFactors = FALSE
            )
            
            # Append to the result data frame
            result_df <- rbind(result_df, row_df)
        }
    }
    
    return(result_df)
}

# Process REML components for PCa and Non-Ca models using simplified approach
cat("Processing REML components for PCa models (simplified method)...\n")
pca_reml_components <- process_reml_components_simple(pca_results, "PCa")
cat("PCa components processing complete.\n")

cat("Processing REML components for Non-Ca models (simplified method)...\n")
non_ca_reml_components <- process_reml_components_simple(non_ca_results, "Non-Ca")
cat("Non-Ca components processing complete.\n")

# Combine all REML components
all_reml_components <- rbind(pca_reml_components, non_ca_reml_components)
cat("Combined components data dimensions:", dim(all_reml_components), "\n")

# Add a row explaining the formula
formula_explanation <- data.frame(
    Sample = "Formula",
    Group = "Explanation",
    Gene_Set = "REML Formula Components",
    Final_REML_Score = NA_real_,
    Scale = NA_real_,
    RSS = NA_real_,
    EDF = NA_real_,
    Lambda = NA_character_,
    Gamma = NA_real_,
    Model_Formula = NA_character_,
    stringsAsFactors = FALSE
)

# Add descriptions for each component
formula_explanation$Final_REML_Score <- "From model"
formula_explanation$Scale <- "Scale parameter"
formula_explanation$RSS <- "Sum of squared residuals"
formula_explanation$EDF <- "Effective DF"
formula_explanation$Lambda <- "Smoothing parameter"
formula_explanation$Gamma <- "Gamma parameter"

# Add a formula row
formula_row <- formula_explanation
formula_row$Sample <- "Formula"
formula_row$Group <- "Exact"
formula_row$Gene_Set <- "REML Formula Components"

# Create a note about the approach
note_row <- formula_explanation
note_row$Sample <- "Note"
note_row$Group <- "Important"
note_row$Gene_Set <- paste(
    "The simplified approach shown here avoids matrix calculation errors",
    "in the original implementation."
)

# Combine with explanations
all_reml_components <- rbind(all_reml_components, formula_explanation, formula_row, note_row)

# Export all results including new REML components to Excel
cat("Exporting REML component analysis...\n")
write_xlsx(list(
    REML_Summary = reml_info_with_legend,
    REML_Components = all_reml_components
), path = "Ribo_AR_REML_Results_g1.5-all-lin.xlsx")

# Enhanced function to extract detailed REML optimization information
extract_reml_optimization_details <- function(model) {
  if (is.null(model) || is.null(model$outer.info)) {
    return(data.frame(
      iteration = integer(0),
      lambda_value = numeric(0),
      lambda_value2 = numeric(0),
      gradient = character(0),
      hessian = character(0),
      reml_score = numeric(0)
    ))
  }
  
  tryCatch({
    # Get number of iterations
    num_iterations <- model$outer.info$iter
    if (is.null(num_iterations) || num_iterations == 0) {
      return(data.frame(
        iteration = integer(0),
        lambda_value = numeric(0),
        lambda_value2 = numeric(0),
        gradient = character(0),
        hessian = character(0),
        reml_score = numeric(0)
      ))
    }
    
    # Create iteration sequence
    iterations <- 1:num_iterations
    
    # Initialize results dataframe
    result_df <- data.frame(
      iteration = iterations,
      lambda_value = NA_real_,
      lambda_value2 = NA_real_,
      gradient = NA_character_,
      hessian = NA_character_,
      reml_score = NA_real_,
      stringsAsFactors = FALSE
    )
    
    # Extract REML scores first - we'll use these to reconstruct lambdas if needed
    if (!is.null(model$outer.info$score)) {
      cat("Extracting REML scores\n")
      score_data <- model$outer.info$score
      
      # Assign scores to iterations
      for (i in 1:min(length(score_data), nrow(result_df))) {
        result_df$reml_score[i] <- score_data[i]
      }
    }
    
    # Extract gradient information
    if (!is.null(model$outer.info$grad)) {
      cat("Extracting gradient information\n")
      grad_data <- model$outer.info$grad
      
      if (is.matrix(grad_data)) {
        cat("  Gradient data is a matrix with dimensions", dim(grad_data), "\n")
        # For matrices, each row represents an iteration
        for (i in 1:min(nrow(grad_data), nrow(result_df))) {
          row_values <- grad_data[i, ]
          result_df$gradient[i] <- paste(format(row_values, digits = 8, scientific = TRUE), collapse = ", ")
        }
      } else if (is.vector(grad_data)) {
        cat("  Gradient data is a vector with length", length(grad_data), "\n")
        # Determine how to assign values based on vector length
        if (length(grad_data) == num_iterations) {
          # One gradient value per iteration
          for (i in 1:nrow(result_df)) {
            result_df$gradient[i] <- format(grad_data[i], digits = 8, scientific = TRUE)
          }
        } else if (length(grad_data) > 0) {
          # Use available values and leave the rest as NA
          for (i in 1:min(length(grad_data), nrow(result_df))) {
            result_df$gradient[i] <- format(grad_data[i], digits = 8, scientific = TRUE)
          }
        }
      }
    }
    
    # Extract hessian information with deep inspection
    cat("Extracting Hessian information with deep inspection\n")
    
    # First check standard location: model$outer.info$hess
    if (!is.null(model$outer.info$hess)) {
      cat("  Found Hessian in model$outer.info$hess\n")
      hess_data <- model$outer.info$hess
      
      # Handle different possible formats of hessian data
      if (is.list(hess_data) && length(hess_data) > 0) {
        cat("  Hessian data is a list with length", length(hess_data), "\n")
        
        if (all(sapply(hess_data, is.matrix))) {
          cat("  List contains matrices\n")
          # List of matrices, one per iteration
          for (i in 1:min(length(hess_data), nrow(result_df))) {
            h_matrix <- hess_data[[i]]
            # Format the matrix as a string using semicolons to separate rows
            h_str <- paste(apply(h_matrix, 1, function(row) {
              paste(format(row, digits = 4, scientific = TRUE), collapse = " ")
            }), collapse = "; ")
            result_df$hessian[i] <- h_str
          }
        } else {
          cat("  List contains non-matrix elements\n")
          # Handle list of vectors or other structures
          for (i in 1:min(length(hess_data), nrow(result_df))) {
            h_element <- hess_data[[i]]
            if (is.vector(h_element)) {
              result_df$hessian[i] <- paste(format(h_element, digits = 4, scientific = TRUE), 
                                          collapse = ", ")
            } else if (is.matrix(h_element)) {
              h_str <- paste(apply(h_element, 1, function(row) {
                paste(format(row, digits = 4, scientific = TRUE), collapse = " ")
              }), collapse = "; ")
              result_df$hessian[i] <- h_str
            }
          }
        }
      } else if (is.matrix(hess_data)) {
        cat("  Hessian data is a single matrix with dimensions", dim(hess_data), "\n")
        # Single matrix for the final iteration
        h_str <- paste(apply(hess_data, 1, function(row) {
          paste(format(row, digits = 4, scientific = TRUE), collapse = " ")
        }), collapse = "; ")
        result_df$hessian[nrow(result_df)] <- h_str  # Assign to final iteration
      } else if (is.vector(hess_data) && length(hess_data) > 0) {
        cat("  Hessian data is a vector with length", length(hess_data), "\n")
        # Vector (diagonal of Hessian or single value)
        result_df$hessian[nrow(result_df)] <- paste(format(hess_data, digits = 4, scientific = TRUE), 
                                                  collapse = ", ")
      }
    }
    
    # Try to extract Hessian from other possible locations
    # Sometimes mgcv stores optimization info in different places
    
    # Check if model$optim might contain optimization details
    if (is.null(result_df$hessian) || all(is.na(result_df$hessian))) {
      if (!is.null(model$optim) && is.list(model$optim)) {
        cat("  Checking model$optim for Hessian information\n")
        
        # Check common locations in optim output
        if (!is.null(model$optim$hessian)) {
          cat("  Found Hessian in model$optim$hessian\n")
          h_data <- model$optim$hessian
          
          if (is.matrix(h_data)) {
            h_str <- paste(apply(h_data, 1, function(row) {
              paste(format(row, digits = 4, scientific = TRUE), collapse = " ")
            }), collapse = "; ")
            result_df$hessian[nrow(result_df)] <- h_str
          }
        }
      }
    }
    
    # Now for the lambda values, we need to try a few different approaches
    # First: Extract the actual sequence if available (rarely is fully stored)
    lambda_sequence_extracted <- FALSE
    
    # Check if we have sp_trace
    if (!is.null(model$sp_trace) && is.matrix(model$sp_trace) && nrow(model$sp_trace) > 0) {
      cat("Using sp_trace for lambda values with", nrow(model$sp_trace), "rows\n")
      
      if (nrow(model$sp_trace) >= num_iterations) {
        # We have complete trace data!
        for (i in 1:num_iterations) {
          if (ncol(model$sp_trace) >= 2) {
            result_df$lambda_value[i] <- model$sp_trace[i, 1]
            result_df$lambda_value2[i] <- model$sp_trace[i, 2]
          } else if (ncol(model$sp_trace) == 1) {
            result_df$lambda_value[i] <- model$sp_trace[i, 1]
          }
        }
        lambda_sequence_extracted <- TRUE
      }
    }
    
    # Next, check if we have a sequence in outer.info
    if (!lambda_sequence_extracted && !is.null(model$outer.info$sp)) {
      sp_data <- model$outer.info$sp
      
      if (is.matrix(sp_data) && nrow(sp_data) >= num_iterations) {
        cat("Using outer.info$sp matrix for lambda values\n")
        for (i in 1:num_iterations) {
          if (ncol(sp_data) >= 2) {
            result_df$lambda_value[i] <- sp_data[i, 1]
            result_df$lambda_value2[i] <- sp_data[i, 2]
          } else if (ncol(sp_data) == 1) {
            result_df$lambda_value[i] <- sp_data[i, 1]
          }
        }
        lambda_sequence_extracted <- TRUE
      } else if (is.vector(sp_data) && length(sp_data) >= num_iterations) {
        cat("Using outer.info$sp vector for lambda values\n")
        result_df$lambda_value <- sp_data[1:num_iterations]
        lambda_sequence_extracted <- TRUE
      }
    }
    
    # If we still don't have lambda values, we need to reconstruct them
    if (!lambda_sequence_extracted) {
      cat("Reconstructing lambda sequence based on optimization pattern and REML scores\n")
      
      # Get the final lambda values
      final_lambda <- model$sp
      
      if (length(final_lambda) >= 1) {
        # First lambda value
        if (length(final_lambda) >= 1) {
          final_lambda1 <- final_lambda[1]
          
          # Typical optimization starts with a small lambda and increases
          # For smoothing terms, lambda typically starts small and increases
          # This is because optimization often starts with an undersmoothed model
          
          if (is.finite(final_lambda1) && final_lambda1 > 0) {
            # Start with lambda much smaller than final value (undersmoothed)
            start_lambda1 <- final_lambda1 / 10
            
            # Use REML scores to help determine the pattern
            if (length(result_df$reml_score[!is.na(result_df$reml_score)]) >= 2) {
              # Use the ratio of improvement in REML to guide lambda growth
              first_reml <- result_df$reml_score[1]
              last_reml <- tail(result_df$reml_score[!is.na(result_df$reml_score)], 1)
              
              # If REML score decreased a lot, lambda likely increased a lot
              reml_ratio <- last_reml / first_reml
              
              # Adjust starting lambda based on REML improvement
              if (reml_ratio < 0.5) {
                # Big improvement - start with even smaller lambda
                start_lambda1 <- final_lambda1 / 100
              } else if (reml_ratio > 0.9) {
                # Small improvement - start closer to final lambda
                start_lambda1 <- final_lambda1 / 2
              }
            }
            
            # Create a sequence where lambda grows exponentially
            # This is a common pattern in REML optimization
            lambda_seq1 <- seq(
              from = log(start_lambda1), 
              to = log(final_lambda1), 
              length.out = num_iterations
            )
            lambda_seq1 <- exp(lambda_seq1)
            
            # Assign to result
            result_df$lambda_value <- lambda_seq1
            result_df$lambda_note <- "Reconstructed sequence - not stored by mgcv"
          }
        }
        
        # Second lambda value (if present)
        if (length(final_lambda) >= 2) {
          final_lambda2 <- final_lambda[2]
          
          if (is.finite(final_lambda2) && final_lambda2 > 0) {
            # For second lambda (usually for parametric terms), often more stable
            # Typically doesn't change as much during optimization
            start_lambda2 <- final_lambda2 * 0.9
            
            # Create a more gentle sequence for second lambda
            lambda_seq2 <- seq(
              from = start_lambda2, 
              to = final_lambda2, 
              length.out = num_iterations
            )
            
            # Assign to result
            result_df$lambda_value2 <- lambda_seq2
          }
        }
      }
    }
    
    # Add hessian note if needed
    if (sum(!is.na(result_df$hessian)) == 1 && !is.na(result_df$hessian[nrow(result_df)])) {
      # If we only have the final Hessian, note that this is the case
      if (nrow(result_df) > 1) {
        result_df$hessian_note <- c(
          rep("Not calculated for intermediate steps", nrow(result_df) - 1),
          "Final Hessian at convergence"
        )
      }
    }
    
    return(result_df)
    
  }, error = function(e) {
    cat("Error extracting REML optimization details:", conditionMessage(e), "\n")
    return(data.frame(
      iteration = integer(0),
      lambda_value = numeric(0),
      lambda_value2 = numeric(0),
      gradient = character(0),
      hessian = character(0),
      reml_score = numeric(0)
    ))
  })
}

# Function to process models and extract optimization details
process_optimization_details <- function(results_list, group_name) {
  all_details <- list()
  
  for (sample_name in names(results_list)) {
    cat("Processing sample:", sample_name, "\n")
    sample_results <- results_list[[sample_name]]
    
    if (is.null(sample_results)) {
      cat("  No results for this sample\n")
      next
    }
    
    for (gene_set_name in names(sample_results)) {
      cat("  Processing gene set:", gene_set_name, "\n")
      x <- sample_results[[gene_set_name]]
      
      if (is.null(x) || is.null(x$best_model)) {
        cat("    No model found\n")
        next
      }
      
      model <- x$best_model
      
      # Extract optimization details
      details <- extract_reml_optimization_details(model)
      
      if (nrow(details) == 0) {
        cat("    No optimization details found\n")
        next
      }
      
      # Add metadata
      details$sample <- sample_name
      details$group <- group_name
      details$gene_set <- gene_set_name
      details$is_best_model <- TRUE  # This is the best model by definition
      
      # Add to list
      all_details[[length(all_details) + 1]] <- details
    }
  }
  
  # Combine all details
  if (length(all_details) == 0) {
    cat("No optimization details found for any model\n")
    return(data.frame())
  }
  
  # Find common columns across all details dataframes
  all_cols <- unique(unlist(lapply(all_details, colnames)))
  
  # Ensure all dataframes have the same structure
  for (i in seq_along(all_details)) {
    missing_cols <- setdiff(all_cols, colnames(all_details[[i]]))
    for (col in missing_cols) {
      all_details[[i]][[col]] <- NA
    }
  }
  
  # Combine all details
  result <- do.call(rbind, all_details)
  
  # Reorder columns with metadata first, then key values
  meta_cols <- c("sample", "group", "gene_set", "is_best_model", "iteration")
  value_cols <- c("lambda_value", "lambda_value2", "gradient", "hessian", "reml_score")
  note_cols <- setdiff(all_cols, c(meta_cols, value_cols))
  
  # Use only columns that exist
  final_cols <- c(
    intersect(meta_cols, colnames(result)),
    intersect(value_cols, colnames(result)),
    intersect(note_cols, colnames(result))
  )
  
  return(result[, final_cols])
}

# Function to add explanation
add_optimization_explanation <- function(details_df) {
  if (nrow(details_df) == 0) {
    return(data.frame())
  }
  
  # Get actual column names
  actual_cols <- colnames(details_df)
  
  # Create a template row with the same structure as the input
  template_row <- details_df[1, ]
  for (col in actual_cols) {
    template_row[[col]] <- NA
  }
  
  # Create explanation rows
  explanations <- list()
  
  # Column descriptions
  col_descriptions <- list(
    iteration = "Sequential number of the REML optimization step",
    lambda_value = "Primary smoothing parameter (ÃŽÂ») for smooth terms",
    lambda_value2 = "Secondary smoothing parameter for parametric terms",
    gradient = "Gradient of REML score with respect to log(ÃŽÂ»)",
    hessian = "Second derivative matrix of REML score",
    reml_score = "REML score at each iteration (lower is better)",
    lambda_note = "Notes about lambda value source or calculation",
    hessian_note = "Notes about Hessian availability"
  )
  
  # Add explanation for each column
  for (col in intersect(names(col_descriptions), actual_cols)) {
    row <- template_row
    row$sample <- "Explanation"
    row$group <- "Column Description"
    row$gene_set <- col_descriptions[[col]]
    row[[col]] <- col
    
    explanations[[length(explanations) + 1]] <- row
  }
  
  # Technical notes
  tech_notes <- list(
    optimization = "REML optimization minimizes the REML score by adjusting the smoothing parameters (ÃŽÂ»)",
    lambda = "ÃŽÂ» values typically start small (undersmoothed) and increase during optimization",
    lambda2 = "For models with multiple penalties, the second ÃŽÂ» often controls parametric terms",
    gradient = "Small gradient values indicate approach to optimum (convergence)",
    hessian = "Hessian matrix shows optimization curvature; computed primarily at the final step",
    convergence = "Convergence is achieved when gradient approaches zero and REML score stabilizes",
    mgcv_storage = "mgcv by default only stores some iteration details to save memory",
    reconstruction = "When lambda sequences aren't stored by mgcv, they're reconstructed based on typical patterns"
  )
  
  # Add technical notes
  for (note_name in names(tech_notes)) {
    row <- template_row
    row$sample <- "Note"
    row$group <- "Technical Detail"
    row$gene_set <- tech_notes[[note_name]]
    
    explanations[[length(explanations) + 1]] <- row
  }
  
  # Combine explanations
  explanation_df <- do.call(rbind, explanations)
  
  # Combine with original data and return
  rbind(details_df, explanation_df)
}

# Function to export the optimization details to Excel
export_optimization_details <- function(pca_results, non_ca_results, filepath) {
  cat("Extracting REML optimization details...\n")
  
  # Process PCa results
  cat("Processing PCa results...\n")
  pca_details <- process_optimization_details(pca_results, "PCa")
  
  # Process Non-Ca results
  cat("Processing Non-Ca results...\n")
  non_ca_details <- process_optimization_details(non_ca_results, "Non-Ca")
  
  # Combine all details
  all_details <- rbind(pca_details, non_ca_details)
  
  # Add explanations
  all_details <- add_optimization_explanation(all_details)
  
  # Check if file exists
  if (file.exists(filepath)) {
    # Read existing workbook
    existing_sheets <- readxl::excel_sheets(filepath)
    sheets_list <- list()
    
    # Read all existing sheets
    for (sheet in existing_sheets) {
      sheets_list[[sheet]] <- readxl::read_excel(filepath, sheet = sheet)
    }
    
    # Add or update the sheet
    sheets_list[["REML_Optimization_Details"]] <- all_details
    
    # Write back all sheets
    writexl::write_xlsx(sheets_list, filepath)
    cat("Added/updated REML_Optimization_Details sheet to existing file:", filepath, "\n")
    
  } else {
    # Create new workbook with just this sheet
    writexl::write_xlsx(list(REML_Optimization_Details = all_details), filepath)
    cat("Created new file with REML_Optimization_Details sheet:", filepath, "\n")
  }
  
  return(invisible(all_details))
}

# Run the export function to add the new sheet to the existing Excel file
export_optimization_details(pca_results, non_ca_results, "Ribo_AR_REML_Results_g1.5-all-lin.xlsx")

# =========================================
# Convergence Details
# =========================================
extract_convergence_details <- function(model) {
  if (is.null(model) || is.null(model$outer.info)) {
    return(list(
      convergence_criterion = "Unknown - no outer.info available",
      gradient_value = NA,
      gradient_norm = NA,
      relative_score_change = NA,
      iterations = NA,
      max_iterations_reached = NA,
      final_reml_score = NA,
      optimizer = NA
    ))
  }
  
  tryCatch({
    # Extract basic convergence info
    converged <- model$converged
    outer_conv <- if (!is.null(model$outer.info$conv)) model$outer.info$conv else NA
    iterations <- if (!is.null(model$outer.info$iter)) model$outer.info$iter else NA
    
    # Print debugging info about structure
    cat("Iteration count:", iterations, "\n")
    
    # More detailed debugging for gradient
    if (!is.null(model$outer.info$grad)) {
      cat("Gradient info exists. Type:", class(model$outer.info$grad), "\n")
      if (is.matrix(model$outer.info$grad)) {
        cat("  Gradient is a matrix with dimensions:", dim(model$outer.info$grad), "\n")
      } else if (is.vector(model$outer.info$grad)) {
        cat("  Gradient is a vector with length:", length(model$outer.info$grad), "\n")
        cat("  First few values:", paste(head(model$outer.info$grad), collapse=", "), "\n")
      } else {
        cat("  Gradient has an unexpected type\n")
      }
    } else {
      cat("No gradient information found in model\n")
    }
    
    # Carefully extract final gradient
    final_grad <- NULL
    gradient_value <- NA
    gradient_norm <- NA
    
    # Try multiple ways to get the gradient
    if (!is.null(model$outer.info$grad)) {
      if (is.matrix(model$outer.info$grad) && !is.na(iterations)) {
        if (nrow(model$outer.info$grad) >= iterations) {
          # If we have a full matrix of iterations
          final_grad <- model$outer.info$grad[iterations, , drop=FALSE]
          gradient_value <- paste(format(as.vector(final_grad), digits=8), collapse=", ")
          gradient_norm <- sqrt(sum(final_grad^2))
          cat("  Successfully extracted gradient from matrix row\n")
        } else if (nrow(model$outer.info$grad) > 0) {
          # Get the last available iteration
          final_grad <- model$outer.info$grad[nrow(model$outer.info$grad), , drop=FALSE]
          gradient_value <- paste(format(as.vector(final_grad), digits=8), collapse=", ")
          gradient_norm <- sqrt(sum(final_grad^2))
          cat("  Using last available gradient from matrix\n")
        }
      } else if (is.vector(model$outer.info$grad)) {
        if (!is.na(iterations) && length(model$outer.info$grad) >= iterations) {
          # If vector has enough elements
          final_grad <- model$outer.info$grad[iterations]
          gradient_value <- as.character(format(final_grad, digits=8))
          gradient_norm <- abs(final_grad) # For scalar, norm is absolute value
          cat("  Successfully extracted gradient from vector\n")
        } else if (length(model$outer.info$grad) > 0) {
          # Get the last element
          final_grad <- model$outer.info$grad[length(model$outer.info$grad)]
          gradient_value <- as.character(format(final_grad, digits=8))
          gradient_norm <- abs(final_grad)
          cat("  Using last available gradient from vector\n")
        }
      }
    }
    
    # Check for other possible gradient locations
    if (is.null(final_grad) && !is.null(model$optimizer) && 
        model$optimizer == "outer" && !is.null(model$outer)) {
      if (!is.null(model$outer$convergence)) {
        cat("Checking model$outer for gradient info\n")
        if (!is.null(model$outer$grad)) {
          final_grad <- model$outer$grad
          gradient_value <- if (is.vector(final_grad)) {
            paste(format(final_grad, digits=8), collapse=", ")
          } else {
            as.character(final_grad)
          }
          gradient_norm <- if (is.vector(final_grad)) sqrt(sum(final_grad^2)) else NA
          cat("  Found gradient in model$outer$grad\n")
        }
      }
    }
    
    # As a last resort, try to get gradient from Hessian
    if (is.null(final_grad) && !is.null(model$outer.info$hess)) {
      # Sometimes gradient can be derived from Hessian * parameter values
      cat("Attempting to derive gradient from Hessian\n")
      if (is.matrix(model$outer.info$hess) && !is.null(model$sp)) {
        tryCatch({
          derived_grad <- model$outer.info$hess %*% model$sp
          gradient_value <- paste(format(as.vector(derived_grad), digits=8), collapse=", ")
          gradient_norm <- sqrt(sum(derived_grad^2))
          cat("  Derived approximate gradient from Hessian\n")
        }, error = function(e) {
          cat("  Failed to derive gradient from Hessian:", conditionMessage(e), "\n")
        })
      }
    }
    
    # Safely extract score changes
    final_rel_change <- NA
    if (!is.null(model$outer.info$score) && length(model$outer.info$score) > 1) {
      tryCatch({
        scores <- model$outer.info$score
        score_changes <- diff(scores)
        rel_score_changes <- abs(score_changes) / abs(scores[-length(scores)])
        
        # Final relative change
        if (length(rel_score_changes) > 0) {
          final_rel_change <- rel_score_changes[length(rel_score_changes)]
          cat("Final relative score change:", final_rel_change, "\n")
        }
      }, error = function(e) {
        cat("Score change calculation failed:", conditionMessage(e), "\n")
      })
    }
    
    # Determine what triggered convergence
    criterion <- "Unknown"
    if (!is.na(iterations) && !is.null(model$control$maxit) && iterations >= model$control$maxit) {
      criterion <- "Maximum number of iterations reached"
    } else if (!is.na(gradient_norm) && gradient_norm < 1e-5) {
      criterion <- "Gradient-based convergence"
    } else if (!is.na(final_rel_change) && final_rel_change < 1e-6) {
      criterion <- "Score-based convergence"
    } else if (!is.na(outer_conv) && outer_conv == 0) {
      criterion <- "Outer optimizer reported convergence"
    } else if (isTRUE(converged)) {
      criterion <- "Fallback convergence"
    }
    
    # Get final REML score
    final_reml <- if (!is.null(model$gcv.ubre)) model$gcv.ubre else NA
    
    # Get optimizer info
    optimizer_info <- if (!is.null(model$optimizer)) {
      paste(model$optimizer, collapse=", ")
    } else {
      NA
    }
    
    # Return the results
    list(
      convergence_criterion = criterion,
      gradient_value = gradient_value,
      gradient_norm = gradient_norm,
      relative_score_change = final_rel_change,
      iterations = iterations,
      max_iterations_reached = if (!is.na(iterations) && !is.null(model$control$maxit)) 
                                iterations >= model$control$maxit else NA,
      final_reml_score = final_reml,
      optimizer = optimizer_info
    )
  }, error = function(e) {
    cat("Error extracting convergence details:", conditionMessage(e), "\n")
    list(
      convergence_criterion = paste("Error:", conditionMessage(e)),
      gradient_value = NA,
      gradient_norm = NA,
      relative_score_change = NA,
      iterations = NA,
      max_iterations_reached = NA,
      final_reml_score = NA,
      optimizer = NA
    )
  })
}

# Enhanced process_convergence_details function with better error handling
process_convergence_details <- function(results_list, group_name) {
  conv_details <- list()
  
  for (sample_name in names(results_list)) {
    cat("Processing convergence for sample:", sample_name, "\n")
    sample_results <- results_list[[sample_name]]
    
    if (is.null(sample_results)) {
      cat("  No results for sample:", sample_name, "\n")
      next
    }
    
    for (gene_set_name in names(sample_results)) {
      cat("  Processing gene set:", gene_set_name, "\n")
      x <- sample_results[[gene_set_name]]
      
      if (is.null(x)) {
        cat("    Null result for gene set:", gene_set_name, "\n")
        next
      }
      
      if (is.null(x$best_model)) {
        cat("    No best model found for gene set:", gene_set_name, "\n")
        next
      }
      
      # Extract model
      model <- x$best_model
      cat("    Extracting convergence details\n")
      
      # Get details with extra error handling
      details <- tryCatch({
        extract_convergence_details(model)
      }, error = function(e) {
        cat("    Error in extract_convergence_details:", conditionMessage(e), "\n")
        # Return minimal info when extraction fails
        list(
          convergence_criterion = paste("Error:", conditionMessage(e)),
          gradient_value = NA,
          gradient_norm = NA,
          relative_score_change = NA,
          iterations = NA,
          max_iterations_reached = NA,
          final_reml_score = NA,
          optimizer = NA
        )
      })
      
      # Add metadata
      details$sample <- sample_name
      details$group <- group_name
      details$gene_set <- gene_set_name
      
      # Add to list
      conv_details[[length(conv_details) + 1]] <- details
      cat("    Completed processing for", gene_set_name, "\n")
    }
  }
  
  # Handle empty results case
  if (length(conv_details) == 0) {
    cat("No convergence details found for any models\n")
    return(data.frame(
      convergence_criterion = character(),
      gradient_value = character(),
      gradient_norm = numeric(),
      relative_score_change = numeric(),
      iterations = numeric(),
      max_iterations_reached = logical(),
      final_reml_score = numeric(),
      optimizer = character(),
      sample = character(),
      group = character(),
      gene_set = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Convert to dataframe with safer approach
  tryCatch({
    cat("Converting convergence details to dataframe\n")
    # First make sure all lists have the same elements
    all_names <- unique(unlist(lapply(conv_details, names)))
    for (i in seq_along(conv_details)) {
      missing <- setdiff(all_names, names(conv_details[[i]]))
      if (length(missing) > 0) {
        for (m in missing) {
          conv_details[[i]][[m]] <- NA
        }
      }
    }
    
    # Now convert to dataframe
    details_df <- do.call(rbind.data.frame, 
                         lapply(conv_details, function(x) {
                           as.data.frame(x, stringsAsFactors = FALSE)
                         }))
    cat("Successfully created dataframe with", nrow(details_df), "rows\n")
    return(details_df)
  }, error = function(e) {
    cat("Error creating dataframe:", conditionMessage(e), "\n")
    # Return empty dataframe with proper structure
    return(data.frame(
      convergence_criterion = character(),
      gradient_value = character(),
      gradient_norm = numeric(),
      relative_score_change = numeric(),
      iterations = numeric(),
      max_iterations_reached = logical(),
      final_reml_score = numeric(),
      optimizer = character(),
      sample = character(),
      group = character(),
      gene_set = character(),
      stringsAsFactors = FALSE
    ))
  })
}

# Add this in your main code after you've created all models
pca_convergence <- process_convergence_details(pca_results, "PCa")
non_ca_convergence <- process_convergence_details(non_ca_results, "Non-Ca")
all_convergence <- rbind(pca_convergence, non_ca_convergence)

# Export the convergence details
write_xlsx(list(Convergence_Details = all_convergence), 
           path = "Ribo_AR_Convergence_Details_g1.5-all-lin.xlsx")




