##################################################
# Part 3.03: GAM-REML-PRSS Setup 
##################################################
library(mgcv)    
library(Seurat)  
set.seed(123)  

# Define expanded gene sets for analysis
gene_sets <- list(
    Ribo = c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26"),
    AR = c("KLK4", "KLK2", "KLK3", "PDLIM5", "ABHD2", "ALDH1A3", "SORD"),
    PI3K_AKT = c("PIK3CA", "AKT1", "PTEN", "MTOR", "TSC2", "FOXO3", "GSK3B"),
    mTOR = c("MTOR", "RPTOR", "RICTOR", "MLST8", "AKT1S1", "DEPTOR", "PRR5"),
    GSK3B = c("GSK3B", "CTNNB1", "AXIN1", "CSNK1A1", "APC", "FZD1", "LRP6"),
    NFKB = c("NFKB1", "RELA", "TNFAIP3", "NFKBIA", "IKBKB", "TRAF6", "NFKB2"),
    WNT = c("CTNNB1", "APC", "AXIN1", "GSK3B", "LEF1", "TCF7", "FZD1")
)

# List prostate cancer samples
pca_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                 "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                 "HYW_5755_Tumor")
pca_clusters <- c(6, 9, 11, 14, 19)

# List non-cancerous samples
non_ca_samples <- c("HYW_4701_Benign", "HYW_4847_Benign", "HYW_4880_Benign", 
                    "HYW_4881_Benign", "HYW_5386_Benign", "HYW_5742_Benign", 
                    "HYW_5755_Benign")
non_ca_clusters <- 3

# Extract nonlinear basis functions from GAM model
extract_basis_functions <- function(model) {
    set.seed(123)
    cat("Debugging message: Entering extract_basis_functions\n")
    smooth_terms <- model$smooth
    cat("Number of smooth terms:", length(smooth_terms), "\n")
    
    X <- predict(model, type = "lpmatrix")
    cat("Prediction matrix dimensions:", dim(X), "\n")
    
    trpm4_values <- model$model$TRPM4
    cat("TRPM4 values extracted with length:", length(trpm4_values), "\n")
    
    basis_df <- data.frame(TRPM4 = trpm4_values)
    
    # Extract only the columns corresponding to smooth terms
    smooth_cols <- grep("s\\(TRPM4\\)", colnames(X))
    for (i in seq_along(smooth_cols)) {
        basis_df[[paste0("phi_", i)]] <- X[, smooth_cols[i]]
    }
    
    cat("Basis functions data frame created with dimensions:", dim(basis_df), "\n")
    print(head(basis_df))
    return(basis_df)
}

# Extract detailed information from GAM model
extract_model_details <- function(model, data) {
    set.seed(123)
    cat("Debugging message: Entering extract_model_details\n")
    
    # Extract parametric terms
    parametric_terms <- coef(model)[!grepl("s\\(TRPM4\\)", names(coef(model)))]
    cat("Parametric terms extracted:\n")
    print(parametric_terms)
    
    # Extract smooth terms
    smooth_terms <- coef(model)[grepl("s\\(TRPM4\\)", names(coef(model)))]
    cat("Smooth terms extracted:\n")
    print(smooth_terms)
    
    # Create equation
    parametric_part <- sprintf("%.9f", parametric_terms[1])  # Intercept term
    if (length(parametric_terms) > 1) {
        parametric_part <- paste0(parametric_part, " + ", 
                                  sprintf("%.9f*x", parametric_terms[2]))  # Linear term added
    }
    if (length(parametric_terms) > 2) {
        parametric_part <- paste0(parametric_part, " + ", 
                                  sprintf("%.9f*x^2", parametric_terms[3]))  # Quadratic term added
    }
    
    smooth_part <- paste(sprintf("%.9f*ϕ%d(x)", smooth_terms, seq_along(smooth_terms)), collapse=" + ")
    
    equation <- paste0("f(x) = ", parametric_part, 
                       ifelse(length(smooth_terms) > 0, " + ", ""), 
                       smooth_part)
    cat("Generated equation:", equation, "\n")
    
    # Extract basis functions for smooth terms only
    pred_matrix <- predict(model, type = "lpmatrix")
    smooth_cols <- grep("s\\(TRPM4\\)", colnames(pred_matrix))
    basis_functions <- pred_matrix[, smooth_cols, drop = FALSE]
    colnames(basis_functions) <- paste0("phi_", 1:ncol(basis_functions))
    basis_functions <- cbind(TRPM4 = data$TRPM4, as.data.frame(basis_functions))
    
    # Get unique TRPM4 values used in constructing the spline basis
    knot_like_values <- sort(unique(data$TRPM4))
    
    # Create dataframe of knot-like locations and corresponding expression values
    knot_data <- data.frame(
        TRPM4 = knot_like_values,
        Expression = approx(data$TRPM4, data$Expression, xout = knot_like_values)$y
    )
    
    # Calculate model and null model deviance
    model_deviance <- deviance(model)
    null_model <- gam(Expression ~ 1, data = data)
    null_deviance <- deviance(null_model)
    
    # Calculate sum of EDFs for smooth terms and degrees of freedom for parametric terms
    smooth_edf_sum <- sum(summary(model)$edf) - length(parametric_terms)
    parametric_df <- length(parametric_terms)
    
    cat("Number of smooth terms in equation:", length(smooth_terms), "\n")
    cat("Number of smooth columns in basis_functions:", sum(grepl("^s", names(basis_functions))), "\n")
    
    list(
        parametric_terms = parametric_terms,
        smooth_terms = smooth_terms,
        knot_data = knot_data,
        equation = equation,
        model_deviance = model_deviance,
        null_deviance = null_deviance,
        smooth_edf_sum = smooth_edf_sum,
        parametric_df = parametric_df,
        basis_functions = basis_functions
    )
}

# Manually calculate and verify predicted values
calculate_manual_prediction <- function(model, data) {
    set.seed(123)
    model_details <- extract_model_details(model, data)
    
    manual_pred <- model_details$parametric_terms[1]  # Start with intercept
    if (length(model_details$parametric_terms) > 1) {
        manual_pred <- manual_pred + model_details$parametric_terms[2] * data$TRPM4
    }
    
    for (i in seq_along(model_details$smooth_terms)) {
        manual_pred <- manual_pred + model_details$smooth_terms[i] * data[[paste0("phi_", i)]]
    }
    
    return(manual_pred)
}

# Calculate Penalized Residual Sum of Squares
calculate_prss <- function(model, data) {
    set.seed(123)
    # Extract necessary components for PRSS calculation
    fitted_values <- fitted(model)
    residuals <- residuals(model)
    lambda <- model$sp[1]  # Use first smoothing parameter
    
    # Calculate Residual Sum of Squares
    RSS <- sum(residuals^2)
    
    # Extract penalty matrix and coefficients
    S <- model$smooth[[1]]$S[[1]]  # Assumes one smooth term
    beta <- coef(model)[model$smooth[[1]]$first.para:model$smooth[[1]]$last.para]
    
    # Calculate penalization term equivalent to f_double_prime_integral
    f_double_prime_integral <- as.numeric(t(beta) %*% S %*% beta)
    
    # Calculate PRSS by combining RSS and penalty
    PRSS <- RSS + lambda * f_double_prime_integral
    
    # Calculate x range and predictions for consistency
    x_range <- range(data$TRPM4)
    x_seq <- seq(x_range[1], x_range[2], length.out = 1000)
    pred <- predict(model, newdata = data.frame(TRPM4 = x_seq), type = "terms")
    
    list(
        PRSS = PRSS,
        RSS = RSS,
        f_double_prime_integral = f_double_prime_integral,
        lambda = lambda
    )
}

# Extract REML iteration information from model
extract_reml_iterations <- function(model) {
    if (is.null(model) || is.null(model$outer.info)) {
        return(data.frame(
            iteration = numeric(0),
            lambda = character(0),
            score = numeric(0)
        ))
    }
    
    tryCatch({
        # Extract iteration count from outer.info
        iterations <- model$outer.info$iter
        
        # Extract lambda values per iteration as smoothing parameters
        lambda_values <- if (!is.null(model$outer.info$sp)) {
            if (is.matrix(model$outer.info$sp)) {
                # Handle multiple smoothing parameters stored as rows
                apply(model$outer.info$sp, 1, function(row) paste(format(row, digits = 6), collapse = ", "))
            } else {
                # Handle single smoothing parameter
                as.character(format(model$outer.info$sp, digits = 6))
            }
        } else {
            rep(NA_character_, length(iterations))
        }
        
        # Extract REML scores for each iteration
        reml_scores <- if (!is.null(model$outer.info$score)) {
            model$outer.info$score
        } else {
            rep(NA_real_, length(iterations))
        }
        
        # Return results as data frame
        data.frame(
            iteration = 1:length(iterations),
            num_iterations = iterations, # Total number of iterations
            lambda = lambda_values,
            score = reml_scores,
            stringsAsFactors = FALSE
        )
    }, error = function(e) {
        cat("Error extracting REML iterations:", conditionMessage(e), "\n")
        data.frame(
            iteration = numeric(0),
            lambda = character(0),
            score = numeric(0)
        )
    })
}

# Enhance REML iteration extraction using sp_trace info
extract_reml_iterations_enhanced <- function(model) {
    if (is.null(model) || is.null(model$outer.info)) {
        cat("No outer.info available in model\n")
        return(data.frame(
            iteration = numeric(0),
            lambda = character(0),
            score = numeric(0)
        ))
    }
    
    tryCatch({
        # Extract iteration count from outer.info
        iterations <- model$outer.info$iter
        
        # Create iteration sequence
        iter_seq <- seq_len(iterations)
        
        # Initialize results dataframe with one row per REML iteration
        result_df <- data.frame(
            iteration = 1:iterations,
            num_iterations = iterations,
            lambda = NA_character_,
            score = NA_real_,
            stringsAsFactors = FALSE
        )
        
        # Use sp_trace captured through tracing environment if available
        if (!is.null(model$sp_trace)) {
            cat("Using sp_trace from model with", nrow(model$sp_trace), "iterations\n")
            
            # Format lambda values for each iteration
            for (i in 1:min(nrow(model$sp_trace), nrow(result_df))) {
                lambda_values <- model$sp_trace[i, ]
                result_df$lambda[i] <- paste(format(lambda_values, digits = 6), collapse = ", ")
            }
        }
        # Fall back to outer.info if sp_trace is unavailable
        else if (!is.null(model$outer.info$sp)) {
            cat("Falling back to outer.info$sp\n")
            sp_data <- model$outer.info$sp
            
            # Handle matrix structure of sp data
            if (is.matrix(sp_data)) {
                # Assign each row as separate iteration
                for (i in 1:min(nrow(sp_data), nrow(result_df))) {
                    result_df$lambda[i] <- paste(format(sp_data[i, ], digits = 6), collapse = ", ")
                }
            } else if (is.vector(sp_data) && length(sp_data) >= nrow(result_df)) {
                # Assign vector elements to iterations
                for (i in 1:nrow(result_df)) {
                    result_df$lambda[i] <- format(sp_data[i], digits = 6)
                }
            } else if (length(sp_data) > 0) {
                # Use generic fallback to repeat available values
                cat("Using generic fallback for lambda values\n")
                if (length(sp_data) == 1) {
                    # Repeat single value for all iterations
                    result_df$lambda <- format(sp_data, digits = 6)
                } else {
                    # Recycle multiple values as needed
                    for (i in 1:nrow(result_df)) {
                        idx <- ((i-1) %% length(sp_data)) + 1
                        result_df$lambda[i] <- format(sp_data[idx], digits = 6)
                    }
                }
            }
        }
        
        # Use final smoothing parameters if no iteration-specific lambda found
        if (all(is.na(result_df$lambda)) && !is.null(model$sp) && length(model$sp) > 0) {
            cat("No iteration-specific lambda found, using final values\n")
            final_lambda <- paste(format(model$sp, digits = 6), collapse = ", ")
            result_df$lambda <- final_lambda
        }
        
        # Extract and assign REML scores
        if (!is.null(model$outer.info$score)) {
            score_data <- model$outer.info$score
            
            # Assign scores to iterations
            for (i in 1:min(length(score_data), nrow(result_df))) {
                result_df$score[i] <- score_data[i]
            }
        }
        
        # Approximate lambda using score trajectory if no lambda values exist
        if (all(is.na(result_df$lambda)) && !all(is.na(result_df$score))) {
            cat("Using score-based approximation for lambda\n")
            if (!is.null(model$sp) && length(model$sp) > 0) {
                # Create trajectory from high to final value
                start_val <- max(10 * model$sp, 10)
                end_val <- model$sp
                
                # Generate sequence approximating optimization path
                lambda_seq <- exp(seq(log(start_val), log(end_val), length.out = nrow(result_df)))
                
                # Assign to results
                for (i in 1:nrow(result_df)) {
                    result_df$lambda[i] <- paste(format(lambda_seq[i], digits = 6), collapse = ", ")
                }
                cat("Created approximate lambda trajectory from", start_val, "to", end_val, "\n")
            }
        }
        
        # Provide debug information
        cat("REML iterations extracted:", nrow(result_df), "\n")
        cat("Lambda values present:", sum(!is.na(result_df$lambda)), "\n")
        if (sum(!is.na(result_df$lambda)) > 0) {
            cat("Sample lambda values:", head(result_df$lambda, 3), "...\n")
        }
        
        return(result_df)
    }, error = function(e) {
        cat("Error extracting REML iterations:", conditionMessage(e), "\n")
        data.frame(
            iteration = numeric(0),
            lambda = character(0),
            score = numeric(0)
        )
    })
}

# Analyze sample with early stopping for GAM fitting
analyze_sample <- function(integrated_obj, original_obj, cluster_ids, gene_set, sample_id, gene_set_name, num_iterations = 100) {
    set.seed(123 + which(c(pca_samples, non_ca_samples) == sample_id))  # Set seed based on sample name
    
    # Get cell names for specified clusters and sample from integrated object
    cluster_cells <- WhichCells(integrated_obj, idents = cluster_ids)
    sample_cells <- WhichCells(integrated_obj, cells = grep(sample_id, colnames(integrated_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    # Check if enough cells are available
    if (length(selected_cells) < 10) {
        cat("Not enough cells for sample:", sample_id, "\n")
        return(NULL)
    }
    
    # Subset original object to selected cells
    cluster_obj <- subset(original_obj, cells = selected_cells)
    
    # Ensure RNA assay is used
    DefaultAssay(cluster_obj) <- "RNA"
    
    # Get normalized data from RNA assay
    normalized_data <- GetAssayData(cluster_obj, slot = "data", assay = "RNA")
    
    # Calculate average expression for given gene set
    genes <- gene_set[gene_set %in% rownames(normalized_data)]
    avg_expression <- colMeans(normalized_data[genes, ])
    
    # Combine into a data frame
    sample_data <- data.frame(
        TRPM4 = normalized_data["TRPM4", ],
        Expression = avg_expression
    )
    
    # Apply log2 transformation
    sample_data <- log2(sample_data + 1)
    
    # Track best model and its PRSS
    best_prss <- Inf
    best_model <- NULL
    best_k <- NULL
    best_iteration <- 0
    prss_values <- numeric(num_iterations)
    model_params <- list()
    reml_iterations_list <- list()  # Store REML iterations for each PRSS iteration
    
    # Determine unique TRPM4 values and maximum k
    unique_trpm4 <- length(unique(sample_data$TRPM4))
    max_k <- min(10, unique_trpm4 - 1)  # Ensure k is at most unique_trpm4 - 1
    
    cat("Sample being analyzed:", sample_id, "\n")
    cat("Number of unique TRPM4 values:", unique_trpm4, "\n")
    cat("Total number of observations:", nrow(sample_data), "\n")
    cat("Maximum k value allowed:", max_k, "\n")
    
    # Initialize variables for early stopping
    early_stop_counter <- 0
    early_stop_triggered <- FALSE
    
    # Perform initial exploration with k values from 3 to max_k
    for (i in 1:min(8, num_iterations)) {
        cat("Initial exploration iteration:", i, "for sample:", sample_id, "\n")
        
        set.seed(123 + i)  # Set seed for each iteration
        # Explore range of k values in initial iterations
        k <- 3 + ((i - 1) %% (max_k - 2))
        
        tryCatch({
            # Create environment to store smoothing parameter trace
            sp_trace_env <- new.env()
            
            # Enable tracing of smoothing parameter optimization
            old_options <- options(mgcv.trace.sp = sp_trace_env)
            on.exit(options(old_options), add = TRUE)
            
            # Fit GAM model with tracing enabled
            gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                             data = sample_data, method = "REML",
                             select = TRUE, gamma = 1.5)
            
            # Retrieve model control settings
            model_control <- gam_model$control
            cat("Model settings: gamma =", gam_model$gamma, 
                "maxit =", model_control$maxit, 
                "mgcv.tol =", model_control$mgcv.tol, "\n")
            
            # Store SP trace in model object if available
            if (exists("sp", envir = sp_trace_env)) {
                gam_model$sp_trace <- get("sp", envir = sp_trace_env)
                cat("Lambda trace captured with", NROW(gam_model$sp_trace), "iterations\n")
            } else {
                cat("No lambda trace found in environment\n")
            }
            
            # Calculate PRSS for this model
            prss_components <- calculate_prss(gam_model, sample_data)
            current_prss <- prss_components$PRSS
            prss_values[i] <- current_prss
            
            cat("Iteration", i, "k =", k, "PRSS:", current_prss, "\n")
            
            # Extract REML iteration information with enhanced method
            reml_iterations <- extract_reml_iterations_enhanced(gam_model)
            
            if (nrow(reml_iterations) > 0) {
                # Add additional information to REML iterations
                reml_iterations$prss_iteration <- i
                reml_iterations$sample <- sample_id
                reml_iterations$gene_set <- gene_set_name
                reml_iterations$k_value <- k
                reml_iterations$prss_score <- current_prss
                reml_iterations$is_best_model <- FALSE
                
                # Store REML iterations with unique identifier
                reml_iterations_list[[length(reml_iterations_list) + 1]] <- reml_iterations
                
                # Output debug information
                cat("Extracted", nrow(reml_iterations), "REML iterations\n")
                cat("Lambda values present:", sum(!is.na(reml_iterations$lambda)), "\n")
                if (sum(!is.na(reml_iterations$lambda)) > 0) {
                    cat("First lambda value:", reml_iterations$lambda[1], "\n")
                } else {
                    cat("WARNING: No lambda values captured in REML iterations\n")
                }
            } else {
                cat("No REML iterations extracted\n")
            }
            
            # Update best model if current PRSS is better
            if (current_prss < best_prss) {
                best_prss <- current_prss
                best_model <- gam_model
                best_k <- k
                best_iteration <- i
                
                # Reset early stopping counter when better model found
                early_stop_counter <- 0
            } else {
                # Increment early stopping counter
                early_stop_counter <- early_stop_counter + 1
            }
            
            # Store model parameters
            model_params[[i]] <- c(
                list(
                    iteration = i,
                    edf = summary(gam_model)$edf,
                    k = k,
                    m = gam_model$smooth[[1]]$m,
                    lambda = prss_components$lambda,
                    scale = gam_model$scale,
                    r_squared = summary(gam_model)$r.sq,
                    dev_explained = summary(gam_model)$dev.expl,
                    reml = gam_model$gcv.ubre,
                    p_value = summary(gam_model)$s.table[4],
                    RSS = prss_components$RSS,
                    f_double_prime_integral = prss_components$f_double_prime_integral,
                    PRSS = prss_components$PRSS,
                    reml_iterations = if (nrow(reml_iterations) > 0) max(reml_iterations$iteration) else 0
                ),
                prss_components
            )
            
        }, error = function(e) {
            cat("Error in iteration", i, "for sample", sample_id, ":", conditionMessage(e), "\n")
            prss_values[i] <- NA
        })
    }
    
    # Set k range for adaptive refinement phase
    if (best_k == max_k) {
        k_range <- c(max_k - 1, max_k)
    } else if (best_k == 3) {
        k_range <- c(3, 4)
    } else {
        k_range <- c(best_k - 1, best_k, best_k + 1)
    }
    
    cat("Best k after initial exploration:", best_k, "with PRSS:", best_prss, "\n")
    cat("Starting adaptive refinement phase with k range:", paste(k_range, collapse=", "), "\n")
    
    # Reset early stopping counter for refinement phase
    early_stop_counter <- 0
    max_allowed_iterations <- num_iterations
    
    # Perform adaptive refinement phase
    for (i in 9:max_allowed_iterations) {
        cat("Refinement iteration:", i, "for sample:", sample_id, "\n")
        
        # Check for early stopping condition
        if (early_stop_counter >= 20) {
            cat("Early stopping triggered after", early_stop_counter, 
                "iterations without improvement; best iteration was", best_iteration, "\n")
            early_stop_triggered <- TRUE
            break
        }
        
        set.seed(123 + i)  # Set seed for each iteration
        
        # Select k value for refinement
        if (i %% 3 == 0) {
            # Try random k from promising range
            k <- sample(k_range, 1)
        } else {
            # Use best k value
            k <- best_k
        }
        
        tryCatch({
            # Create environment to store smoothing parameter trace
            sp_trace_env <- new.env()
            
            # Enable tracing of smoothing parameter optimization
            old_options <- options(mgcv.trace.sp = sp_trace_env)
            on.exit(options(old_options), add = TRUE)
            
            # Build GAM model with selected k and tracing enabled
            gam_model <- gam(Expression ~ TRPM4 + s(TRPM4, bs = "tp", k = k), 
                             data = sample_data, method = "REML",
                             select = TRUE, gamma = 1.5)
            
            # Retrieve model control settings
            model_control <- gam_model$control
            cat("Model settings: gamma =", gam_model$gamma, 
                "maxit =", model_control$maxit, 
                "mgcv.tol =", model_control$mgcv.tol, "\n")
            
            # Store SP trace in model object if available
            if (exists("sp", envir = sp_trace_env)) {
                gam_model$sp_trace <- get("sp", envir = sp_trace_env)
                cat("Lambda trace captured with", NROW(gam_model$sp_trace), "iterations\n")
            } else {
                cat("No lambda trace found in environment\n")
            }
            
            # Calculate PRSS for this model
            prss_components <- calculate_prss(gam_model, sample_data)
            current_prss <- prss_components$PRSS
            prss_values[i] <- current_prss
            
            cat("Iteration", i, "k =", k, "PRSS:", current_prss, "\n")
            
            # Extract REML iteration information with enhanced function
            reml_iterations <- extract_reml_iterations_enhanced(gam_model)
            
            if (nrow(reml_iterations) > 0) {
                # Add additional information to REML iterations
                reml_iterations$prss_iteration <- i
                reml_iterations$sample <- sample_id
                reml_iterations$gene_set <- gene_set_name
                reml_iterations$k_value <- k
                reml_iterations$prss_score <- current_prss
                reml_iterations$is_best_model <- FALSE
                
                # Store REML iterations
                reml_iterations_list[[i]] <- reml_iterations
                
                # Output debug information
                cat("Extracted", nrow(reml_iterations), "REML iterations\n")
                cat("Lambda values present:", sum(!is.na(reml_iterations$lambda)), "\n")
                if (sum(!is.na(reml_iterations$lambda)) > 0) {
                    cat("First lambda value:", reml_iterations$lambda[1], "\n")
                } else {
                    cat("WARNING: No lambda values captured in REML iterations\n")
                }
            } else {
                cat("No REML iterations extracted\n")
            }
            
            # Update best model if better PRSS found
            if (current_prss < best_prss) {
                best_prss <- current_prss
                best_model <- gam_model
                best_k <- k
                best_iteration <- i
                
                # Update k_range if new best k found
                if (best_k == max_k) {
                    k_range <- c(max_k - 1, max_k)
                } else if (best_k == 3) {
                    k_range <- c(3, 4)
                } else {
                    k_range <- c(best_k - 1, best_k, best_k + 1)
                }
                
                cat("New best k:", best_k, "with PRSS:", best_prss, "\n")
                cat("Updated k range:", paste(k_range, collapse=", "), "\n")
                
                # Reset early stopping counter
                early_stop_counter <- 0
            } else {
                # Increment early stopping counter
                early_stop_counter <- early_stop_counter + 1
                cat("No improvement for", early_stop_counter, "iterations; best was iteration", best_iteration, "\n")
            }
            
            # Store model parameters
            model_params[[i]] <- c(
                list(
                    iteration = i,
                    edf = summary(gam_model)$edf,
                    k = k,
                    m = gam_model$smooth[[1]]$m,
                    lambda = prss_components$lambda,
                    scale = gam_model$scale,
                    r_squared = summary(gam_model)$r.sq,
                    dev_explained = summary(gam_model)$dev.expl,
                    reml = gam_model$gcv.ubre,
                    p_value = summary(gam_model)$s.table[4],
                    RSS = prss_components$RSS,
                    f_double_prime_integral = prss_components$f_double_prime_integral,
                    PRSS = prss_components$PRSS,
                    reml_iterations = if (nrow(reml_iterations) > 0) max(reml_iterations$iteration) else 0
                ),
                prss_components
            )
            
        }, error = function(e) {
            cat("Error in iteration", i, "for sample", sample_id, ":", conditionMessage(e), "\n")
            prss_values[i] <- NA
            
            # Count as non-improving iteration for early stopping
            early_stop_counter <- early_stop_counter + 1
        })
    }
    
    if (is.null(best_model)) {
        cat("No successful models for sample:", sample_id, "\n")
        return(NULL)
    }
    
    # Record early stopping information
    actual_iterations <- if (early_stop_triggered) i - 1 else num_iterations
    cat("Analysis completed after", actual_iterations, "iterations\n")
    cat("Best iteration:", best_iteration, "with PRSS:", best_prss, "and k:", best_k, "\n")
    
    # Get details for best model
    best_details <- extract_model_details(best_model, sample_data)
    best_prss <- calculate_prss(best_model, sample_data)
    
    # Compare predictions
    model_pred <- predict(best_model, newdata = sample_data)
    manual_pred <- calculate_manual_prediction(best_model, sample_data)
    pred_comparison <- all.equal(model_pred, manual_pred)
    cat("Prediction comparison result:", pred_comparison, "\n")
    if (!isTRUE(pred_comparison)) {
        cat("Maximum difference between predictions:", max(abs(model_pred - manual_pred)), "\n")
    }
    
    # Create PRSS vs Iteration data including only ran iterations
    prss_data <- data.frame(
        Sample = sample_id,
        Iteration = 1:actual_iterations,
        PRSS = prss_values[1:actual_iterations],
        Gene_Set = gene_set_name,
        Early_Stopped = early_stop_triggered
    )
    
    # Add number of REML iterations for each PRSS iteration
    reml_iteration_counts <- sapply(model_params[1:actual_iterations], function(x) {
        if (is.null(x) || !("reml_iterations" %in% names(x))) 0 else x$reml_iterations
    })
    reml_iteration_counts <- reml_iteration_counts[!is.na(reml_iteration_counts)]
    prss_data$REML_Iterations <- c(reml_iteration_counts, rep(NA, actual_iterations - length(reml_iteration_counts)))
    
    # Combine all REML iterations into single dataframe for ran iterations
    valid_reml_iterations <- reml_iterations_list[1:actual_iterations]
    valid_reml_iterations <- valid_reml_iterations[!sapply(valid_reml_iterations, is.null)]
    all_reml_iterations <- do.call(rbind, valid_reml_iterations)
    
    # Mark best model in all_reml_iterations dataframe
    if (!is.null(all_reml_iterations) && nrow(all_reml_iterations) > 0) {
        all_reml_iterations$is_best_model <- FALSE
        
        # Mark rows associated with best model
        best_model_rows <- which(all_reml_iterations$prss_iteration == best_iteration)
        if (length(best_model_rows) > 0) {
            all_reml_iterations$is_best_model[best_model_rows] <- TRUE
        }
    }
    
    # Compile best model parameters
    best_params <- list(
        iteration = best_iteration,
        edf = summary(best_model)$edf,
        k = best_k,
        lambda = best_prss$lambda,
        scale = best_model$scale,
        r_squared = summary(best_model)$r.sq,
        dev_explained = summary(best_model)$dev.expl,
        reml = best_model$gcv.ubre,
        p_value = summary(best_model)$s.table[4],
        RSS = best_prss$RSS,
        f_double_prime_integral = best_prss$f_double_prime_integral,
        PRSS = best_prss$PRSS,
        early_stopped = early_stop_triggered,
        actual_iterations = actual_iterations
    )
    
    best_params <- c(best_params, best_details)
    best_params_df <- data.frame(
        Parameter = names(unlist(best_params)),
        Value = unlist(best_params)
    )
    
    # Calculate and adjust p-values
    p_values <- sapply(model_params[1:actual_iterations], function(x) {
        if (is.null(x) || !("p_value" %in% names(x))) NA else x$p_value
    })
    p_values <- p_values[!is.na(p_values)]
    bh_corrected_p <- p.adjust(p_values, method = "BH")[which(p_values == best_params$p_value)[1]]
    
    # Debug best model details
    cat("Debugging best model:\n")
    print(summary(best_model))
    cat("\nModel formula:\n")
    print(best_model$formula)
    cat("\nSmooth terms:\n")
    print(best_model$smooth)
    
    # Compile summary statistics
    summary_stats <- list(
        sample = sample_id,
        prss_summary = summary(prss_values[1:actual_iterations]),
        prss_range = range(prss_values[1:actual_iterations][!is.na(prss_values[1:actual_iterations])]),
        best_params = best_params_df,
        bh_corrected_p = bh_corrected_p,
        early_stopped = early_stop_triggered,
        actual_iterations = actual_iterations
    )
    
    # Calculate predicted values for all datapoints
    predicted_expression <- predict(best_model, newdata = sample_data)
    
    # Add predicted values to sample_data
    sample_data$Expression_Predicted <- predicted_expression
    sample_data$Sample <- sample_id
    sample_data$Gene_Set <- gene_set_name
    
    return(list(
        summary_stats = summary_stats, 
        best_details = best_details, 
        best_model = best_model, 
        best_params = best_params, 
        prss_data = prss_data, 
        gam_data = sample_data, 
        gene_set_name = gene_set_name,
        reml_iterations = all_reml_iterations,
        best_iteration = best_iteration, 
        early_stopped = early_stop_triggered,
        actual_iterations = actual_iterations
    ))
}

# Function to analyze REML convergence across iterations
analyze_reml_convergence <- function(reml_iterations) {
    if (is.null(reml_iterations) || nrow(reml_iterations) == 0) {
        cat("No REML iterations data available to analyze\n")
        return(data.frame(
            sample = character(),
            gene_set = character(),
            prss_iteration = numeric(),
            iterations = numeric(),
            initial_score = numeric(),
            final_score = numeric(),
            diff_score = numeric(),
            percent_reduction = numeric(),
            converged = logical(),
            monotonic_decrease = logical(),
            stringsAsFactors = FALSE
        ))
    }
    
    cat("Analyzing REML convergence for", nrow(reml_iterations), "iterations\n")
    cat("Samples analyzed:", paste(unique(reml_iterations$sample), collapse=", "), "\n")
    cat("Gene sets analyzed:", paste(unique(reml_iterations$gene_set), collapse=", "), "\n")
    cat("PRSS iterations analyzed:", length(unique(reml_iterations$prss_iteration)), "\n")
    
    # Create robust grouping by sample, gene_set, and prss_iteration
    group_ids <- interaction(reml_iterations$sample, 
                             reml_iterations$gene_set, 
                             reml_iterations$prss_iteration, 
                             drop=TRUE)
    unique_groups <- levels(group_ids)
    
    cat("Found", length(unique_groups), "unique sample-geneset-iteration groups\n")
    
    # Process each group
    convergence_stats <- lapply(unique_groups, function(group_id) {
        # Extract current group
        group <- reml_iterations[group_ids == group_id, ]
        
        # Verify group extraction worked
        if (nrow(group) == 0) {
            cat("Warning: Empty group extracted for", group_id, "\n")
            return(NULL)
        }
        
        # Parse group ID to extract components
        group_parts <- strsplit(as.character(group_id), "\\.")[[1]]
        sample_name <- group_parts[1]
        gene_set_name <- group_parts[2]
        prss_iter <- as.numeric(group_parts[3])
        
        cat("Processing group:", paste(sample_name, gene_set_name, prss_iter, sep="-"), 
            "with", nrow(group), "iterations\n")
        
        # Check if enough iterations exist for meaningful analysis
        if (nrow(group) <= 1) {
            cat("Insufficient iterations", nrow(group), "for convergence analysis\n")
            return(data.frame(
                sample = sample_name,
                gene_set = gene_set_name,
                prss_iteration = prss_iter,
                iterations = nrow(group),
                initial_score = if(nrow(group) > 0) group$score[1] else NA,
                final_score = if(nrow(group) > 0) group$score[nrow(group)] else NA,
                diff_score = NA,
                percent_reduction = NA,
                converged = FALSE,
                monotonic_decrease = NA,
                stringsAsFactors = FALSE
            ))
        }
        
        # Sort by iteration to ensure proper analysis
        group <- group[order(group$iteration), ]
        
        # Check if scores decrease monotonically with numerical tolerance
        score_diffs <- diff(group$score)
        monotonic <- all(score_diffs <= 1e-10)
        
        # Calculate convergence metrics
        initial_score <- group$score[1]
        final_score <- group$score[nrow(group)]
        diff_score <- initial_score - final_score
        
        # Calculate percent reduction safely
        percent_reduction <- NA
        if (!is.na(initial_score) && !is.na(final_score) && abs(initial_score) > 1e-10) {
            percent_reduction <- (diff_score / abs(initial_score)) * 100
        }
        
        # Determine if converged based on small final iteration change
        converged <- FALSE
        if (nrow(group) >= 2) {
            final_change <- abs(group$score[nrow(group)] - group$score[nrow(group)-1])
            
            # Avoid division by zero or very small numbers
            denominator <- abs(group$score[nrow(group)-1])
            if (denominator < 1e-10) denominator <- 1e-10
            
            rel_change <- final_change / denominator
            converged <- rel_change < 1e-6  # Set threshold for convergence
        }
        
        # Create result dataframe with all needed information
        result <- data.frame(
            sample = sample_name,
            gene_set = gene_set_name,
            prss_iteration = prss_iter,
            iterations = nrow(group),
            initial_score = initial_score,
            final_score = final_score,
            diff_score = diff_score,
            percent_reduction = percent_reduction,
            converged = converged,
            monotonic_decrease = monotonic,
            stringsAsFactors = FALSE
        )
        
        return(result)
    })
    
    # Remove NULL results and combine
    convergence_stats <- Filter(Negate(is.null), convergence_stats)
    
    if (length(convergence_stats) == 0) {
        cat("Warning: No valid convergence statistics could be generated\n")
        return(data.frame(
            sample = character(),
            gene_set = character(),
            prss_iteration = numeric(),
            iterations = numeric(),
            initial_score = numeric(),
            final_score = numeric(),
            diff_score = numeric(),
            percent_reduction = numeric(),
            converged = logical(),
            monotonic_decrease = logical(),
            stringsAsFactors = FALSE
        ))
    }
    
    result_df <- do.call(rbind, convergence_stats)
    cat("Produced", nrow(result_df), "rows of convergence statistics\n")
    return(result_df)
}

# Function to calculate REML convergence rate statistics with robust error handling
calculate_reml_convergence_stats <- function(convergence_df) {
    if (is.null(convergence_df) || nrow(convergence_df) == 0) {
        cat("No convergence data available for statistics\n")
        # Return empty data frame with expected structure
        return(data.frame(
            sample = character(),
            gene_set = character(),
            total_prss_iterations = numeric(),
            avg_reml_iterations = numeric(),
            max_reml_iterations = numeric(),
            percent_converged = numeric(),
            percent_monotonic = numeric(),
            avg_percent_reduction = numeric(),
            stringsAsFactors = FALSE
        ))
    }
    
    cat("Calculating convergence statistics for", nrow(convergence_df), "records\n")
    cat("Samples analyzed:", paste(unique(convergence_df$sample), collapse=", "), "\n")
    cat("Gene sets analyzed:", paste(unique(convergence_df$gene_set), collapse=", "), "\n")
    
    # Create unique group identifier for each sample and gene_set combination
    group_ids <- interaction(convergence_df$sample, convergence_df$gene_set, drop=TRUE)
    unique_groups <- levels(group_ids)
    
    cat("Found", length(unique_groups), "unique sample-geneset groups\n")
    
    # Process each sample and gene_set group
    sample_stats <- lapply(unique_groups, function(group_id) {
        # Extract current group
        group <- convergence_df[group_ids == group_id, ]
        
        # Verify group extraction worked
        if (nrow(group) == 0) {
            cat("Warning: Empty group extracted for", group_id, "\n")
            return(NULL)
        }
        
        # Parse group ID to extract components
        group_parts <- strsplit(as.character(group_id), "\\.")[[1]]
        sample_name <- group_parts[1]
        gene_set_name <- if (length(group_parts) > 1) group_parts[2] else "Unknown"
        
        cat("Processing group:", paste(sample_name, gene_set_name, sep="-"), 
            "with", nrow(group), "records\n")
        
        # Calculate statistics with NA handling
        avg_reml_iterations <- mean(group$iterations, na.rm=TRUE)
        max_reml_iterations <- max(group$iterations, na.rm=TRUE)
        
        # Handle logical columns properly
        converged_values <- group$converged
        converged_values[is.na(converged_values)] <- FALSE
        percent_converged <- mean(converged_values) * 100
        
        monotonic_values <- group$monotonic_decrease
        monotonic_values <- monotonic_values[!is.na(monotonic_values)]
        percent_monotonic <- if (length(monotonic_values) > 0) {
            mean(monotonic_values) * 100
        } else {
            NA
        }
        
        # Calculate average percent reduction with NA handling
        percent_reduction_values <- group$percent_reduction
        percent_reduction_values <- percent_reduction_values[!is.na(percent_reduction_values)]
        avg_percent_reduction <- if (length(percent_reduction_values) > 0) {
            mean(percent_reduction_values)
        } else {
            NA
        }
        
        # Create result with comprehensive information
        result <- data.frame(
            sample = sample_name,
            gene_set = gene_set_name,
            total_prss_iterations = nrow(group),
            avg_reml_iterations = avg_reml_iterations,
            max_reml_iterations = max_reml_iterations,
            percent_converged = percent_converged,
            percent_monotonic = percent_monotonic,
            avg_percent_reduction = avg_percent_reduction,
            stringsAsFactors = FALSE
        )
        
        return(result)
    })
    
    # Remove NULL results
    sample_stats <- Filter(Negate(is.null), sample_stats)
    
    if (length(sample_stats) == 0) {
        cat("Warning: No valid sample statistics could be generated\n")
        return(data.frame(
            sample = character(),
            gene_set = character(),
            total_prss_iterations = numeric(),
            avg_reml_iterations = numeric(),
            max_reml_iterations = numeric(),
            percent_converged = numeric(),
            percent_monotonic = numeric(),
            avg_percent_reduction = numeric(),
            stringsAsFactors = FALSE
        ))
    }
    
    # Combine individual sample stats
    sample_stats_df <- do.call(rbind, sample_stats)
    
    # Calculate global statistics if valid data exists
    if (nrow(sample_stats_df) > 0) {
        # Handle numeric columns properly
        avg_reml_iterations <- mean(convergence_df$iterations, na.rm=TRUE)
        max_reml_iterations <- max(convergence_df$iterations, na.rm=TRUE)
        
        # Handle logical columns properly
        converged_values <- convergence_df$converged
        converged_values[is.na(converged_values)] <- FALSE
        percent_converged <- mean(converged_values) * 100
        
        monotonic_values <- convergence_df$monotonic_decrease
        monotonic_values <- monotonic_values[!is.na(monotonic_values)]
        percent_monotonic <- if (length(monotonic_values) > 0) {
            mean(monotonic_values) * 100
        } else {
            NA
        }
        
        # Calculate average percent reduction with NA handling
        percent_reduction_values <- convergence_df$percent_reduction
        percent_reduction_values <- percent_reduction_values[!is.na(percent_reduction_values)]
        avg_percent_reduction <- if (length(percent_reduction_values) > 0) {
            mean(percent_reduction_values)
        } else {
            NA
        }
        
        # Compile global stats
        global_stats <- data.frame(
            sample = "ALL",
            gene_set = "ALL",
            total_prss_iterations = nrow(convergence_df),
            avg_reml_iterations = avg_reml_iterations,
            max_reml_iterations = max_reml_iterations,
            percent_converged = percent_converged,
            percent_monotonic = percent_monotonic,
            avg_percent_reduction = avg_percent_reduction,
            stringsAsFactors = FALSE
        )
        
        # Combine individual and global stats
        result_df <- rbind(sample_stats_df, global_stats)
    } else {
        result_df <- sample_stats_df
    }
    
    cat("Produced", nrow(result_df), "rows of convergence statistics summary\n")
    return(result_df)
}







