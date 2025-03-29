##################################################
# Part 3.2: GAM-PRSS-REML Analysis and Results Export
##################################################
# Analyze all gene sets for a given sample
analyze_multiple_gene_sets <- function(integrated_obj, original_obj, cluster_ids, gene_sets, sample) {
    set.seed(123 + which(c(pca_samples, non_ca_samples) == sample)) 
    results <- list()
    for (gene_set_name in names(gene_sets)) {
        cat("Analyzing gene set:", gene_set_name, "for sample:", sample, "\n")
        result <- analyze_sample(integrated_obj, original_obj, cluster_ids, gene_sets[[gene_set_name]], sample)
        if (!is.null(result)) {
            results[[gene_set_name]] <- result
        }
    }
    return(results)
}

# Format GAM equation for readability
format_equation <- function(parametric_terms, smooth_terms) {
    # Format parametric part of the equation
    parametric_part <- sprintf("%.4f", parametric_terms[1])
    if (length(parametric_terms) > 1) {
        for (i in 2:length(parametric_terms)) {
            term <- parametric_terms[i]
            if (abs(term) > 1e-4) {  # Include terms not close to zero
                if (term >= 0) {
                    parametric_part <- paste0(parametric_part, " + ", sprintf("%.4f", term), "x^", i-1)
                } else {
                    parametric_part <- paste0(parametric_part, " - ", sprintf("%.4f", abs(term)), "x^", i-1)
                }
            }
        }
    }
    
    # Format smooth part of the equation
    smooth_part <- paste(sprintf("%.4fφ%d(x)", smooth_terms, 
                                 seq_along(smooth_terms)), 
                         collapse = " + ")
    smooth_part <- gsub("\\+ -", "- ", smooth_part)
    
    # Combine parametric and smooth parts into full equation
    full_equation <- paste("f(x) =", parametric_part)
    if (nchar(smooth_part) > 0) {
        full_equation <- paste(full_equation, "+", smooth_part)
    }
    
    return(full_equation)
}

# Create summary dataframe from analysis results
create_summary_df <- function(results_list) {
    set.seed(123)
    summaries <- lapply(names(results_list), function(sample_name) {
        sample_results <- results_list[[sample_name]]
        if (is.null(sample_results)) {
            return(data.frame(Sample = sample_name, PRSS = NA, lambda = NA,
                              PRSS_Min = NA, PRSS_Median = NA, PRSS_Max = NA,
                              Best_EDF = NA, Best_K = NA, Best_R_Squared = NA, 
                              Best_Dev_Explained = NA, Best_P_Value = NA, BH_Corrected_P = NA,
                              Parametric_Terms = NA, Smooth_Terms = NA, Equation = NA,
                              Model_Deviance = NA, Null_Deviance = NA,
                              Smooth_EDF_Sum = NA, Parametric_DF = NA,
                              Knot_Locations = NA, f_double_prime_integral = NA, RSS = NA))
        }
        
        gene_set_summaries <- lapply(names(sample_results), function(gene_set_name) {
            x <- sample_results[[gene_set_name]]
            tryCatch({
                if (is.null(x) || is.null(x$summary_stats) || is.null(x$best_details)) {
                    return(NULL)
                }
                
                best_params <- x$summary_stats$best_params
                
                data.frame(
                    Sample = sample_name,
                    Gene_Set = gene_set_name,
                    PRSS = best_params$Value[best_params$Parameter %in% c("PRSS", "PRSS.s(TRPM4)1")],
                    RSS = best_params$Value[best_params$Parameter == "RSS"],
                    lambda = best_params$Value[best_params$Parameter %in% c("lambda", "lambda.s(TRPM4)1")],
                    f_double_prime_integral = best_params$Value[best_params$Parameter == "f_double_prime_integral"],
                    PRSS_Min = x$summary_stats$prss_range[1],
                    PRSS_Median = median(x$summary_stats$prss_summary),
                    PRSS_Max = x$summary_stats$prss_range[2],
                    Best_EDF = best_params$Value[best_params$Parameter == "edf"],
                    Best_K = best_params$Value[best_params$Parameter == "k"],
                    Best_R_Squared = best_params$Value[best_params$Parameter == "r_squared"],
                    Best_Dev_Explained = best_params$Value[best_params$Parameter == "dev_explained"],
                    Best_P_Value = best_params$Value[best_params$Parameter == "p_value"],
                    BH_Corrected_P = x$summary_stats$bh_corrected_p,
                    Parametric_Terms = paste(x$best_details$parametric_terms, collapse = ", "),
                    Smooth_Terms = paste(x$best_details$smooth_terms, collapse = ", "),
                    Equation = x$best_details$equation,
                    Model_Deviance = x$best_details$model_deviance,
                    Null_Deviance = x$best_details$null_deviance,
                    Smooth_EDF_Sum = x$best_details$smooth_edf_sum,
                    Parametric_DF = x$best_details$parametric_df,
                    Knot_Locations = paste(x$best_details$knot_data$TRPM4, collapse = ", "),
                    stringsAsFactors = FALSE
                )
            }, error = function(e) {
                cat("Error processing sample:", sample_name, "gene set:", gene_set_name, "- Error:", conditionMessage(e), "\n")
                return(NULL)
            })
        })
        
        do.call(rbind, gene_set_summaries)
    })
    
    do.call(rbind, summaries)
}

# Process prostate cancer samples in main workflow
pca_results <- lapply(pca_samples, function(sample) {
    cat("\n======== Processing PCa sample:", sample, "========\n")
    tryCatch({
        set.seed(123 + which(pca_samples == sample))
        gene_set_results <- lapply(names(gene_sets), function(gene_set_name) {
            cat("\n---- Analyzing gene set:", gene_set_name, "----\n")
            detailed_result <- analyze_sample(
                prostate_results$seurat_obj, 
                prostate_ca_seurat, 
                pca_clusters, 
                gene_sets[[gene_set_name]], 
                sample, 
                gene_set_name
            )
            
            if (is.null(detailed_result)) {
                cat("No valid data available for sample:", sample, "gene set:", gene_set_name, "\n")
                return(NULL)
            }
            
            # Validate REML iterations data
            if (is.null(detailed_result$reml_iterations) || nrow(detailed_result$reml_iterations) == 0) {
                cat("Warning: No REML iterations found for", sample, gene_set_name, "\n")
            } else {
                cat("Found", nrow(detailed_result$reml_iterations), "REML iterations for", sample, gene_set_name, "\n")
            }
            
            detailed_result
        })
        names(gene_set_results) <- names(gene_sets)
        gene_set_results
    }, error = function(e) {
        cat("Error in processing PCa sample", sample, ":", conditionMessage(e), "\n")
        return(NULL)
    })
})
names(pca_results) <- pca_samples

# Process non-cancerous samples in main workflow
non_ca_results <- lapply(non_ca_samples, function(sample) {
    tryCatch({
        set.seed(123 + which(non_ca_samples == sample))
        gene_set_results <- lapply(names(gene_sets), function(gene_set_name) {
            detailed_result <- analyze_sample(non_cancerous_results$seurat_obj, non_cancerous_seurat, non_ca_clusters, gene_sets[[gene_set_name]], sample, gene_set_name)
            if (is.null(detailed_result)) {
                cat("No valid data available for sample:", sample, "gene set:", gene_set_name, "\n")
                return(NULL)
            }
            detailed_result
        })
        names(gene_set_results) <- names(gene_sets)
        gene_set_results
    }, error = function(e) {
        cat("Error in processing non-Ca sample", sample, ":", conditionMessage(e), "\n")
        return(NULL)
    })
})
names(non_ca_results) <- non_ca_samples

# Extract best models and parameters from results
extract_best_models_and_params <- function(results) {
    lapply(results, function(sample_result) {
        if (is.null(sample_result)) return(NULL)
        lapply(sample_result, function(gene_set_result) {
            if (is.null(gene_set_result)) return(NULL)
            list(
                best_model = gene_set_result$best_model,
                best_params = gene_set_result$best_params
            )
        })
    })
}

pca_best_models_and_params <- extract_best_models_and_params(pca_results)
non_ca_best_models_and_params <- extract_best_models_and_params(non_ca_results)

# Print debugging information for sample counts
cat("Number of PCa samples:", length(pca_samples), "\n")
cat("Number of non-Ca samples:", length(non_ca_samples), "\n")
cat("Number of PCa results:", length(pca_results), "\n")
cat("Number of non-Ca results:", length(non_ca_results), "\n")

# Count non-null models and parameters
count_non_null <- function(results, key) {
    sum(sapply(results, function(sample) {
        if (is.null(sample)) return(0)
        sum(sapply(sample, function(gene_set) !is.null(gene_set[[key]])))
    }))
}

cat("Number of non-null PCa best models:", count_non_null(pca_best_models_and_params, "best_model"), "\n")
cat("Number of non-null non-Ca best models:", count_non_null(non_ca_best_models_and_params, "best_model"), "\n")
cat("Number of non-null PCa best params:", count_non_null(pca_best_models_and_params, "best_params"), "\n")
cat("Number of non-null non-Ca best params:", count_non_null(non_ca_best_models_and_params, "best_params"), "\n")

# Safely print first non-null element from list
safe_print_first <- function(list_obj, name) {
    non_null <- Filter(Negate(is.null), unlist(list_obj, recursive = FALSE))
    if (length(non_null) > 0) {
        cat("First non-null", name, ":\n")
        if (name == "best model") {
            print(non_null[[1]]$best_model)
        } else if (name == "best params") {
            print(non_null[[1]]$best_params)
        }
    } else {
        cat("No non-null", name, "found.\n")
    }
}

# Print first non-null best model and parameters for PCa and non-Ca
safe_print_first(pca_best_models_and_params, "PCa best model")
safe_print_first(pca_best_models_and_params, "PCa best params")
safe_print_first(non_ca_best_models_and_params, "non-Ca best model")
safe_print_first(non_ca_best_models_and_params, "non-Ca best params")

# Create summary dataframes for PCa and non-Ca results
pca_summary <- create_summary_df(pca_results)
non_ca_summary <- create_summary_df(non_ca_results)

# Export summary results to Excel
write_xlsx(list(PCa_Summary = pca_summary, 
                NonCa_Summary = non_ca_summary), 
           path = "Ribo_AR_All_Analysis_Summary_g1.5-all-lin.xlsx")

# Calculate predicted ribosomal values from model
calculate_predicted_ribo <- function(model, data) {
    predicted <- predict(model, newdata = data)
    return(predicted)
}

# Export detailed results with debugging information
export_detailed_results <- function(results, filename) {
    set.seed(123)
    cat("Debugging message: Entering export_detailed_results function\n")
    detailed_list <- lapply(names(results), function(sample_name) {
        cat("Debugging message: Processing sample:", sample_name, "\n")
        sample_results <- results[[sample_name]]
        if (is.null(sample_results)) {
            cat("Debugging message: Null data for sample:", sample_name, "\n")
            return(NULL)
        }
        
        gene_set_data <- lapply(names(sample_results), function(gene_set_name) {
            x <- sample_results[[gene_set_name]]
            tryCatch({
                if (is.null(x) || is.null(x$summary_stats) || is.null(x$best_details) || is.null(x$best_model)) {
                    cat("Debugging message: Incomplete data for sample:", sample_name, "gene set:", gene_set_name, "\n")
                    return(NULL)
                }
                
                basis_functions <- tryCatch({
                    bf <- extract_basis_functions(x$best_model)
                    bf$Sample <- sample_name
                    bf$Gene_Set <- gene_set_name
                    cat("Debugging message: basis_functions dimensions for", gene_set_name, ":", dim(bf), "\n")
                    bf
                }, error = function(e) {
                    cat("Error extracting basis functions for sample:", sample_name, "gene set:", gene_set_name, "- Error:", conditionMessage(e), "\n")
                    return(NULL)
                })
                
                list(datapoints = x$gam_data, basis_functions = basis_functions)
            }, error = function(e) {
                cat("Error processing sample:", sample_name, "gene set:", gene_set_name, "- Error:", conditionMessage(e), "\n")
                NULL
            })
        })
        
        gene_set_data
    })
    
    cat("Debugging message: Filtering out NULL results\n")
    detailed_list <- Filter(Negate(is.null), unlist(detailed_list, recursive = FALSE))
    if (length(detailed_list) == 0) {
        cat("No valid data to export.\n")
        return()
    }
    
    cat("Debugging message: Extracting datapoints and basis_functions\n")
    datapoints <- lapply(detailed_list, function(x) x$datapoints)
    basis_functions <- lapply(detailed_list, function(x) x$basis_functions)
    
    cat("Debugging message: Combining data\n")
    all_datapoints <- do.call(rbind, datapoints)
    
    export_list <- list(Datapoints = all_datapoints)
    
    if (!is.null(basis_functions) && length(basis_functions) > 0) {
        cat("Debugging message: Processing basis functions for export\n")
        max_phi <- max(sapply(basis_functions, function(bf) sum(grepl("^phi_", names(bf)))))
        combined_basis_functions <- do.call(rbind, lapply(basis_functions, function(bf) {
            phi_cols <- grep("^phi_", names(bf), value = TRUE)
            missing_cols <- max_phi - length(phi_cols)
            if (missing_cols > 0) {
                for (i in 1:missing_cols) {
                    bf[[paste0("phi_", length(phi_cols) + i)]] <- NA
                }
            }
            bf
        }))
        cat("combined_basis_functions dimensions:", dim(combined_basis_functions), "\n")
        export_list$Basis_Functions <- combined_basis_functions
    } else {
        cat("No valid basis functions to export.\n")
    }
    
    tryCatch({
        write_xlsx(export_list, path = filename)
        cat("Data exported successfully to", filename, "\n")
    }, error = function(e) {
        cat("Error writing to Excel:", conditionMessage(e), "\n")
        print(str(export_list))
    })
}

# Export detailed results for PCa samples
tryCatch({
    export_detailed_results(pca_results, "Ribo_AR_All_Detailed_PCa_g1.5-all-lin.xlsx")
    cat("PCa detailed results for Ribo exported successfully.\n")
}, error = function(e) {
    cat("Error exporting PCa detailed results for Ribo:", conditionMessage(e), "\n")
})

# Export detailed results for non-Ca samples
tryCatch({
    export_detailed_results(non_ca_results, "Ribo_AR_All_Detailed_NonCa_g1.5-all-lin.xlsx")
    cat("Non-Ca detailed results for Ribo exported successfully.\n")
}, error = function(e) {
    cat("Error exporting Non-Ca detailed results for Ribo:", conditionMessage(e), "\n")
})

# Export PRSS data with early stopping, k, and lambda values
export_prss_data <- function(results, filename) {
    set.seed(123)  # Ensure reproducibility
    
    cat("\n====== Exporting PRSS and REML data to", filename, "======\n")
    
    # Extract PRSS data with detailed debugging
    cat("Extracting PRSS data from results...\n")
    prss_data_list <- lapply(names(results), function(sample_name) {
        cat("Processing PRSS data for sample:", sample_name, "\n")
        sample_results <- results[[sample_name]]
        
        if (is.null(sample_results)) {
            cat("No results for this sample\n")
            return(NULL)
        }
        
        # Process each gene set within sample
        gene_set_data <- lapply(names(sample_results), function(gene_set_name) {
            cat("Processing gene set:", gene_set_name, "\n")
            prss_data <- sample_results[[gene_set_name]]$prss_data
            
            if (is.null(prss_data) || nrow(prss_data) == 0) {
                cat("No PRSS data found\n")
                return(NULL)
            }
            
            cat("Found", nrow(prss_data), "PRSS data points\n")
            
            # Ensure gene set is recorded correctly
            if (!"Gene_Set" %in% colnames(prss_data) || all(is.na(prss_data$Gene_Set))) {
                prss_data$Gene_Set <- gene_set_name
            }
            
            # Add best iteration information
            best_iteration <- sample_results[[gene_set_name]]$best_iteration
            if (!is.null(best_iteration)) {
                prss_data$Is_Best_Model <- prss_data$Iteration == best_iteration
            } else {
                prss_data$Is_Best_Model <- FALSE
            }
            
            # Add early stopping information if missing
            if (!"Early_Stopped" %in% colnames(prss_data)) {
                if (!is.null(sample_results[[gene_set_name]]$early_stopped)) {
                    prss_data$Early_Stopped <- sample_results[[gene_set_name]]$early_stopped
                } else {
                    prss_data$Early_Stopped <- FALSE
                }
            }
            
            return(prss_data)
        })
        
        # Combine data for all gene sets in sample
        gene_set_data <- Filter(Negate(is.null), gene_set_data)
        
        if (length(gene_set_data) == 0) {
            cat("No valid PRSS data for any gene set in this sample\n")
            return(NULL)
        }
        
        # Combine all gene sets for this sample
        combined_data <- do.call(rbind, gene_set_data)
        cat("Combined", nrow(combined_data), "PRSS data points for sample", sample_name, "\n")
        return(combined_data)
    })
    
    # Extract all REML iterations with improved error handling
    cat("\nExtracting REML iterations from results...\n")
    reml_iterations_list <- lapply(names(results), function(sample_name) {
        cat("Extracting REML iterations for sample:", sample_name, "\n")
        sample_results <- results[[sample_name]]
        
        if (is.null(sample_results)) {
            cat("No results for this sample\n")
            return(NULL)
        }
        
        # Process each gene set for REML iterations
        gene_set_iterations <- lapply(names(sample_results), function(gene_set_name) {
            cat("Processing gene set:", gene_set_name, "\n")
            reml_iterations <- sample_results[[gene_set_name]]$reml_iterations
            
            if (is.null(reml_iterations) || nrow(reml_iterations) == 0) {
                cat("No REML iterations found\n")
                return(NULL)
            }
            
            cat("Found", nrow(reml_iterations), "REML iterations\n")
            
            # Ensure required columns exist
            if (is.null(reml_iterations$sample) || all(is.na(reml_iterations$sample))) {
                reml_iterations$sample <- sample_name
            }
            
            if (is.null(reml_iterations$gene_set) || all(is.na(reml_iterations$gene_set))) {
                reml_iterations$gene_set <- gene_set_name
            }
            
            # Ensure PRSS score column is present
            if (is.null(reml_iterations$prss_score)) {
                cat("Adding PRSS score column\n")
                # Get PRSS scores from prss_data if available
                if (!is.null(sample_results[[gene_set_name]]$prss_data)) {
                    prss_data <- sample_results[[gene_set_name]]$prss_data
                    
                    # Map PRSS scores to REML iterations by prss_iteration
                    reml_iterations$prss_score <- sapply(reml_iterations$prss_iteration, function(iter) {
                        idx <- which(prss_data$Iteration == iter)
                        if (length(idx) > 0) prss_data$PRSS[idx[1]] else NA_real_
                    })
                } else {
                    reml_iterations$prss_score <- NA_real_
                }
            }
            
            # Ensure is_best_model column exists
            if (is.null(reml_iterations$is_best_model)) {
                cat("Adding best model flag column\n")
                # Use best iteration to flag best model
                best_iteration <- sample_results[[gene_set_name]]$best_iteration
                if (!is.null(best_iteration)) {
                    reml_iterations$is_best_model <- reml_iterations$prss_iteration == best_iteration
                } else {
                    reml_iterations$is_best_model <- FALSE
                }
            }
            
            # Ensure numeric columns are numeric
            numeric_cols <- c("iteration", "prss_iteration", "num_iterations", "score", "k_value", "prss_score")
            for (col in numeric_cols) {
                if (col %in% colnames(reml_iterations) && !is.numeric(reml_iterations[[col]])) {
                    reml_iterations[[col]] <- as.numeric(reml_iterations[[col]])
                }
            }
            
            # Ensure logical columns are logical
            if ("is_best_model" %in% colnames(reml_iterations) && !is.logical(reml_iterations$is_best_model)) {
                reml_iterations$is_best_model <- as.logical(reml_iterations$is_best_model)
            }
            
            return(reml_iterations)
        })
        
        # Combine iterations for all gene sets in sample
        gene_set_iterations <- Filter(Negate(is.null), gene_set_iterations)
        
        if (length(gene_set_iterations) == 0) {
            cat("No valid REML iterations for any gene set in this sample\n")
            return(NULL)
        }
        
        # Combine all gene sets with error handling for column mismatches
        combined_iterations <- tryCatch({
            do.call(rbind, gene_set_iterations)
        }, error = function(e) {
            cat("Error combining REML iterations:", conditionMessage(e), "\n")
            
            # Find all column names across dataframes
            all_cols <- unique(unlist(lapply(gene_set_iterations, colnames)))
            
            # Create unified dataframe with all possible columns
            unified_data <- lapply(gene_set_iterations, function(df) {
                missing_cols <- setdiff(all_cols, colnames(df))
                for (col in missing_cols) {
                    if (col == "is_best_model") {
                        df[[col]] <- FALSE
                    } else if (col == "prss_score") {
                        df[[col]] <- NA_real_
                    } else {
                        df[[col]] <- NA
                    }
                }
                df[, all_cols]  # Ensure consistent column order
            })
            
            do.call(rbind, unified_data)
        })
        
        cat("Combined", nrow(combined_iterations), "REML iterations for sample", sample_name, "\n")
        return(combined_iterations)
    })
    
    # Filter out NULL entries and combine all PRSS data
    prss_data_list <- Filter(Negate(is.null), prss_data_list)
    cat("\nFound", length(prss_data_list), "sample(s) with PRSS data\n")
    
    all_prss_data <- if (length(prss_data_list) > 0) {
        # Combine all samples with error handling for column mismatches
        combined_data <- tryCatch({
            do.call(rbind, prss_data_list)
        }, error = function(e) {
            cat("Error combining PRSS data:", conditionMessage(e), "\n")
            
            # Find all column names across dataframes
            all_cols <- unique(unlist(lapply(prss_data_list, colnames)))
            
            # Create unified dataframe with all possible columns
            unified_data <- lapply(prss_data_list, function(df) {
                missing_cols <- setdiff(all_cols, colnames(df))
                for (col in missing_cols) {
                    if (col == "Is_Best_Model") {
                        df[[col]] <- FALSE
                    } else if (col == "Early_Stopped") {
                        df[[col]] <- FALSE
                    } else {
                        df[[col]] <- NA
                    }
                }
                df[, all_cols]  # Ensure consistent column order
            })
            
            do.call(rbind, unified_data)
        })
        
        cat("Combined", nrow(combined_data), "PRSS data points across all samples\n")
        combined_data
    } else {
        cat("No PRSS data available to export\n")
        data.frame()  # Return empty dataframe instead of NULL
    }
    
    # Filter out NULL entries and combine all REML iterations
    reml_iterations_list <- Filter(Negate(is.null), reml_iterations_list)
    cat("\nFound", length(reml_iterations_list), "sample(s) with REML iterations\n")
    
    all_reml_iterations <- if (length(reml_iterations_list) > 0) {
        # Combine all samples with error handling for column mismatches
        combined_iterations <- tryCatch({
            do.call(rbind, reml_iterations_list)
        }, error = function(e) {
            cat("Error combining REML iterations:", conditionMessage(e), "\n")
            
            # Find all column names across dataframes
            all_cols <- unique(unlist(lapply(reml_iterations_list, colnames)))
            
            # Create unified dataframe with all possible columns
            unified_data <- lapply(reml_iterations_list, function(df) {
                missing_cols <- setdiff(all_cols, colnames(df))
                for (col in missing_cols) {
                    if (col == "is_best_model") {
                        df[[col]] <- FALSE
                    } else if (col == "prss_score") {
                        df[[col]] <- NA_real_
                    } else {
                        df[[col]] <- NA
                    }
                }
                df[, all_cols]  # Ensure consistent column order
            })
            
            do.call(rbind, unified_data)
        })
        
        cat("Combined", nrow(combined_iterations), "REML iterations across all samples\n")
        combined_iterations
    } else {
        cat("No REML iterations data available to export\n")
        data.frame()  # Return empty dataframe instead of NULL
    }
    
    # Prepare export list
    export_list <- list()
    
    # Add PRSS data if available
    if (nrow(all_prss_data) > 0) {
        export_list$PRSS_Data <- all_prss_data
    }
    
    # Add REML iterations data if available
    if (nrow(all_reml_iterations) > 0) {
        # Make column names readable for Excel
        names(all_reml_iterations) <- gsub("is_best_model", "Is_Best_Model", names(all_reml_iterations))
        names(all_reml_iterations) <- gsub("prss_score", "PRSS_Score", names(all_reml_iterations))
        names(all_reml_iterations) <- gsub("prss_iteration", "PRSS_Iteration", names(all_reml_iterations))
        names(all_reml_iterations) <- gsub("sample", "Sample", names(all_reml_iterations))
        names(all_reml_iterations) <- gsub("gene_set", "Gene_Set", names(all_reml_iterations))
        names(all_reml_iterations) <- gsub("k_value", "K_Value", names(all_reml_iterations))
        
        # Add Excel-friendly boolean column
        all_reml_iterations$Best_Model_Flag <- ifelse(all_reml_iterations$Is_Best_Model, "YES", "NO")
        
        export_list$REML_Iterations <- all_reml_iterations
        
        # Analyze REML convergence patterns
        cat("\nAnalyzing REML convergence patterns...\n")
        convergence_analysis <- analyze_reml_convergence(all_reml_iterations)
        
        if (nrow(convergence_analysis) > 0) {
            export_list$REML_Convergence_Analysis <- convergence_analysis
            
            # Calculate overall convergence statistics
            cat("Calculating convergence statistics...\n")
            convergence_stats <- calculate_reml_convergence_stats(convergence_analysis)
            
            if (nrow(convergence_stats) > 0) {
                export_list$REML_Convergence_Stats <- convergence_stats
            } else {
                cat("Warning: No convergence statistics produced\n")
            }
        } else {
            cat("Warning: No convergence analysis produced\n")
        }
        
        # Create summary of REML iterations by PRSS iteration
        cat("Creating REML summary by PRSS iteration...\n")
        
        # Use robust aggregation approach
        reml_summary <- tryCatch({
            # Group by sample, gene_set, and PRSS_Iteration
            if ("Sample" %in% names(all_reml_iterations) && 
                "Gene_Set" %in% names(all_reml_iterations) && 
                "PRSS_Iteration" %in% names(all_reml_iterations)) {
                
                groups <- split(all_reml_iterations, 
                                list(all_reml_iterations$Sample, 
                                     all_reml_iterations$Gene_Set, 
                                     all_reml_iterations$PRSS_Iteration))
                
                # Process each group
                summary_rows <- lapply(names(groups), function(group_name) {
                    group_data <- groups[[group_name]]
                    
                    # Parse group name
                    parts <- strsplit(group_name, "\\.")[[1]]
                    sample_name <- parts[1]
                    gene_set_name <- parts[2]
                    prss_iter <- as.numeric(parts[3])
                    
                    # Calculate statistics
                    iter_count <- nrow(group_data)
                    min_score <- min(group_data$score, na.rm=TRUE)
                    max_score <- max(group_data$score, na.rm=TRUE)
                    
                    # Get PRSS score and best model flag
                    prss_score <- if ("PRSS_Score" %in% names(group_data)) {
                        unique(group_data$PRSS_Score)[1]
                    } else {
                        NA_real_
                    }
                    
                    is_best <- if ("Is_Best_Model" %in% names(group_data)) {
                        any(group_data$Is_Best_Model)
                    } else {
                        FALSE
                    }
                    
                    # Get k value
                    k_value <- if ("K_Value" %in% names(group_data)) {
                        unique(group_data$K_Value)[1]
                    } else {
                        NA
                    }
                    
                    # Get lambda value from last iteration
                    lambda_value <- if ("lambda" %in% names(group_data) && nrow(group_data) > 0) {
                        sorted_data <- group_data[order(group_data$iteration),]
                        sorted_data$lambda[nrow(sorted_data)]
                    } else {
                        NA
                    }
                    
                    # Create summary row
                    data.frame(
                        Sample = sample_name,
                        Gene_Set = gene_set_name,
                        PRSS_Iteration = prss_iter,
                        REML_Iteration_Count = iter_count,
                        REML_Min_Score = min_score,
                        REML_Max_Score = max_score,
                        PRSS_Score = prss_score,
                        Is_Best_Model = is_best,
                        Best_Model_Flag = ifelse(is_best, "YES", "NO"),
                        K_Value = k_value,
                        Lambda_Value = lambda_value,
                        stringsAsFactors = FALSE
                    )
                })
                
                # Combine all rows
                do.call(rbind, summary_rows)
            } else {
                cat("Missing required columns for REML summary\n")
                data.frame()
            }
        }, error = function(e) {
            cat("Error creating REML summary:", conditionMessage(e), "\n")
            data.frame()
        })
        
        # Add score reduction calculations if data is valid
        if (nrow(reml_summary) > 0) {
            # Calculate REML score reduction
            reml_summary$REML_Score_Reduction <- reml_summary$REML_Max_Score - reml_summary$REML_Min_Score
            
            # Calculate percent reduction safely
            reml_summary$REML_Percent_Reduction <- ifelse(
                abs(reml_summary$REML_Max_Score) > 1e-10,
                100 * reml_summary$REML_Score_Reduction / abs(reml_summary$REML_Max_Score),
                NA
            )
            
            export_list$REML_Summary <- reml_summary
        } else {
            cat("Warning: No REML summary produced\n")
        }
    }
    
    # Create detailed explanations sheet
    cat("Creating REML explanations sheet...\n")
    reml_explanations <- data.frame(
        Topic = c(
            "REML Iterations", 
            "REML Convergence", 
            "Convergence Analysis", 
            "Score Interpretation",
            "Optimization Process",
            "PRSS Score Column",
            "Is_Best_Model Column",
            "Best_Model_Flag Column",
            "Early_Stopped Column",
            "K_Value Column",
            "Lambda_Value Column"
        ),
        Explanation = c(
            "Each PRSS iteration involves fitting a GAM model using REML. This sheet shows how many REML iterations were needed within each PRSS iteration.", 
            "REML optimization should show decreasing REML scores until convergence. Monotonic decrease indicates proper optimization behavior.",
            "The Convergence Analysis sheet shows detailed metrics for each PRSS iteration, including whether REML scores decreased monotonically and the percent reduction in REML score.",
            "Lower REML scores are better. The final REML score represents the optimized model fit with penalization.",
            "The number of REML iterations needed indicates optimization difficulty. More iterations may suggest complex relationships or ill-conditioned optimization problems.",
            "The PRSS_Score column shows the Penalized Residual Sum of Squares for the model at each PRSS iteration. Lower values indicate better models.",
            "The Is_Best_Model column is TRUE for iterations that belong to the model with the lowest PRSS score (the final selected model).",
            "The Best_Model_Flag column displays 'YES' for iterations belonging to the final selected model, making it easier to filter in Excel.",
            "The Early_Stopped column shows whether the PRSS optimization process was terminated early due to no improvement for 20 consecutive iterations.",
            "The K_Value column shows the basis dimension (k) used for the smooth term in each model fitting iteration.",
            "The Lambda_Value column shows the smoothing parameter (λ) value used in the final REML optimization for each model."
        ),
        stringsAsFactors = FALSE
    )
    
    export_list$REML_Explanations <- reml_explanations
    
    # Verify data exists for export
    if (length(export_list) == 0) {
        cat("Warning: No data to export\n")
        return(FALSE)
    }
    
    # Write data to Excel file
    cat("\nWriting data to Excel file:", filename, "\n")
    tryCatch({
        write_xlsx(export_list, path = filename)
        cat("PRSS and REML data exported successfully to", filename, "\n")
        return(TRUE)
    }, error = function(e) {
        cat("Error writing to Excel file:", conditionMessage(e), "\n")
        return(FALSE)
    })
}

# Export PRSS data for PCa samples
tryCatch({
    export_prss_data(pca_results, "Ribo_AR_All_vs_Iteration_PCa_g1.5-all-lin.xlsx")
    cat("PCa PRSS vs Iteration data for Ribo exported successfully.\n")
}, error = function(e) {
    cat("Error exporting PCa PRSS vs Iteration data for Ribo:", conditionMessage(e), "\n")
})

# Export PRSS data for non-Ca samples
tryCatch({
    export_prss_data(non_ca_results, "Ribo_AR_All_vs_Iteration_NonCa_g1.5-all-lin.xlsx")
    cat("Non-Ca PRSS vs Iteration data for Ribo exported successfully.\n")
}, error = function(e) {
    cat("Error exporting Non-Ca PRSS vs Iteration data for Ribo:", conditionMessage(e), "\n")
})

cat("Analysis complete. Updated results exported to Excel files.\n")

# Return best models for each sample (Note: This section appears incorrect in original code)
pca_best_models <- lapply(pca_results, function(x) x$detailed_result$best_model)
non_ca_best_models <- lapply(non_ca_results, function(x) x$detailed_result$best_model)

cat("Analysis complete. Updated results exported to Excel files.\n")

# Add PRSS convergence column to exported data
add_prss_convergence <- function(file_path) {
    # Read PRSS_Data sheet from Excel file
    prss_data <- read_excel(file_path, sheet = "PRSS_Data")
    
    # Split data by Sample and Gene_Set for convergence analysis
    groups <- split(prss_data, list(prss_data$Sample, prss_data$Gene_Set))
    
    # Initialize convergence column with default value
    prss_data$Converged <- "NO"
    
    # Check convergence for each group
    for(group_name in names(groups)) {
        group <- groups[[group_name]]
        iterations <- sort(unique(group$Iteration))
        
        # Get PRSS values for each iteration
        prss_values <- sapply(iterations, function(i) {
            group$PRSS[group$Iteration == i]
        })
        
        # Check if converged using last 10 iterations
        if(length(iterations) >= 20) {
            # Calculate relative change in last 10 iterations
            last_prss <- prss_values[(length(prss_values)-9):length(prss_values)]
            relative_changes <- abs(diff(last_prss)) / abs(last_prss[-length(last_prss)])
            
            # Mark as converged if all changes are below 0.1%
            if(all(relative_changes < 0.001, na.rm = TRUE)) {
                idx <- prss_data$Sample == group$Sample[1] & prss_data$Gene_Set == group$Gene_Set[1]
                prss_data$Converged[idx] <- "YES"
            }
        }
    }
    
    return(prss_data)
}

# Add REML convergence columns including Is_Best_Model to summary
add_reml_convergence_columns <- function(file_path) {
    # Check if REML_Summary sheet exists in file
    sheets <- excel_sheets(file_path)
    if(!"REML_Summary" %in% sheets) return(NULL)
    
    # Read REML_Summary sheet
    reml_summary <- read_excel(file_path, sheet = "REML_Summary")
    
    # Add new columns with default values
    reml_summary$REML_Converged <- "NO"
    reml_summary$REML_Monotonic <- "NO"
    
    # Add k and lambda columns with default NA
    reml_summary$K_Value <- NA
    reml_summary$Lambda_Value <- NA
    
    # Add Is_Best_Model column with default FALSE
    reml_summary$Is_Best_Model <- FALSE
    
    # Read REML_Iterations sheet if available for detailed analysis
    has_reml_iterations <- "REML_Iterations" %in% sheets
    if(has_reml_iterations) {
        reml_iterations <- read_excel(file_path, sheet = "REML_Iterations")
    }
    
    # Process each row in REML summary
    for(i in 1:nrow(reml_summary)) {
        # Check REML convergence based on iteration count and reduction
        if(!is.na(reml_summary$REML_Iteration_Count[i]) && 
           !is.na(reml_summary$REML_Percent_Reduction[i])) {
            
            # Mark as converged if reduction > 5% and iterations >= 3
            if(reml_summary$REML_Percent_Reduction[i] > 5 && 
               reml_summary$REML_Iteration_Count[i] >= 3) {
                reml_summary$REML_Converged[i] <- "YES"
            }
        }
        
        # Extract k, lambda, and best model flag from REML iterations
        if(has_reml_iterations) {
            subset_iterations <- reml_iterations %>%
                filter(Sample == reml_summary$Sample[i],
                       Gene_Set == reml_summary$Gene_Set[i],
                       PRSS_Iteration == reml_summary$PRSS_Iteration[i])
            
            if(nrow(subset_iterations) > 0) {
                # Get k value if present
                if("K_Value" %in% colnames(subset_iterations)) {
                    k_val <- subset_iterations$K_Value[1]
                    if(!is.na(k_val)) {
                        reml_summary$K_Value[i] <- k_val
                    }
                }
                
                # Get lambda value from last iteration
                if("lambda" %in% colnames(subset_iterations) && nrow(subset_iterations) > 0) {
                    # Sort by iteration and extract last lambda
                    subset_iterations <- subset_iterations[order(subset_iterations$iteration),]
                    last_lambda <- subset_iterations$lambda[nrow(subset_iterations)]
                    
                    # Extract first lambda value from comma-separated string
                    if(!is.na(last_lambda)) {
                        lambda_parts <- strsplit(last_lambda, ",")[[1]]
                        if(length(lambda_parts) > 0) {
                            first_lambda <- trimws(lambda_parts[1])
                            
                            # Verify numeric format
                            if(grepl("^[0-9.e+-]+$", first_lambda)) {
                                reml_summary$Lambda_Value[i] <- first_lambda
                            }
                        }
                    }
                }
                
                # Check if this is best model
                if("Is_Best_Model" %in% colnames(subset_iterations)) {
                    if(any(subset_iterations$Is_Best_Model, na.rm = TRUE)) {
                        reml_summary$Is_Best_Model[i] <- TRUE
                    }
                } else if("Best_Model_Flag" %in% colnames(subset_iterations)) {
                    if(any(subset_iterations$Best_Model_Flag == "YES", na.rm = TRUE)) {
                        reml_summary$Is_Best_Model[i] <- TRUE
                    }
                }
            }
            
            # Check for monotonic decrease in REML scores
            if(nrow(subset_iterations) >= 2) {
                if("score" %in% colnames(subset_iterations)) {
                    # Sort by iteration
                    subset_iterations <- subset_iterations[order(subset_iterations$iteration),]
                    monotonic <- all(diff(subset_iterations$score) <= 0, na.rm = TRUE)
                    reml_summary$REML_Monotonic[i] <- ifelse(monotonic, "YES", "NO")
                }
            }
        }
    }
    
    # Fallback to PRSS_Score if no best model identified
    if(!any(reml_summary$Is_Best_Model) && "PRSS_Score" %in% colnames(reml_summary)) {
        min_prss <- min(reml_summary$PRSS_Score, na.rm = TRUE)
        min_idx <- which(reml_summary$PRSS_Score == min_prss)
        if(length(min_idx) > 0) {
            reml_summary$Is_Best_Model[min_idx[1]] <- TRUE
        }
    }
    
    # Convert Is_Best_Model to YES/NO format for Excel
    reml_summary$Best_Model_Flag <- ifelse(reml_summary$Is_Best_Model, "YES", "NO")
    
    return(reml_summary)
}

# Process Excel file to update with convergence information
process_file <- function(file_path) {
    cat("Processing file:", file_path, "\n")
    
    # Read all sheets from file
    all_sheets <- excel_sheets(file_path)
    sheets_list <- list()
    
    # Process PRSS_Data sheet if present
    if("PRSS_Data" %in% all_sheets) {
        cat("Adding convergence column to PRSS_Data...\n")
        sheets_list[["PRSS_Data"]] <- add_prss_convergence(file_path)
    }
    
    # Process REML_Summary sheet with additional columns
    if("REML_Summary" %in% all_sheets) {
        cat("Adding convergence columns, k, lambda values, and Is_Best_Model column to REML_Summary...\n")
        sheets_list[["REML_Summary"]] <- add_reml_convergence_columns(file_path)
    }
    
    # Copy other sheets unchanged
    for(sheet in setdiff(all_sheets, names(sheets_list))) {
        cat("Reading sheet:", sheet, "\n")
        sheets_list[[sheet]] <- read_excel(file_path, sheet = sheet)
    }
    
    # Write updated data back to file
    cat("Writing updated data back to:", file_path, "\n")
    write_xlsx(sheets_list, file_path)
    cat("File updated successfully!\n\n")
}

# Define file paths for PCa and non-Ca data
pca_file <- "Ribo_AR_All_vs_Iteration_PCa_g1.5-all-lin.xlsx"
nonca_file <- "Ribo_AR_All_vs_Iteration_NonCa_g1.5-all-lin.xlsx"

# Process files to update with Is_Best_Model column
process_file(pca_file)
process_file(nonca_file)

cat("All files have been updated with Is_Best_Model column in the REML_Summary sheets.\n")