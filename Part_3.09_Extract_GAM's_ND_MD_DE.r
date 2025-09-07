#####################################################################################################
# Part 3.09: Extract GAM's null deviance (ND), model deviance (MD), and deviance explained (DE), and their components
#####################################################################################################
# Function to extract detailed deviance components from GAM models
extract_deviance_components <- function(model, data) {
    set.seed(123)
    cat("Extracting deviance components for GAM model\n")
    
    # Basic model information
    n <- nrow(data)
    response_var <- all.vars(model$formula)[1]  # Get response variable name
    y <- data[[response_var]]
    
    # Fitted values and residuals
    fitted_values <- fitted(model)
    residuals_response <- residuals(model, type = "response")
    
    # Model degrees of freedom
    edf_total <- sum(model$edf)  # Total effective degrees of freedom
    df_residual <- model$df.residual  # Residual degrees of freedom
    df_model <- n - df_residual  # Model degrees of freedom
    
    cat("Model info: n =", n, ", total EDF =", edf_total, ", residual DF =", df_residual, "\n")
    
    # Scale parameter (φ̂) - this is key for Gaussian family
    scale_parameter <- model$scale
    
    # Calculate null model components
    null_model <- gam(as.formula(paste(response_var, "~ 1")), data = data, family = gaussian())
    y_null <- fitted(null_model)  # This should be the mean of y
    y_mean <- mean(y)
    
    # Verify null model prediction is the mean
    cat("Null model fitted value:", unique(y_null)[1], ", Mean of y:", y_mean, "\n")
    
    # Raw sum of squares calculations
    ss_total <- sum((y - y_mean)^2)  # Total sum of squares
    ss_residual <- sum((y - fitted_values)^2)  # Residual sum of squares
    ss_model <- ss_total - ss_residual  # Model sum of squares
    
    # Null deviance calculation
    # For Gaussian family: deviance = sum((y - μ)^2) / φ
    null_deviance_raw <- sum((y - y_mean)^2)
    null_deviance_scaled <- null_deviance_raw / scale_parameter
    
    # Model deviance calculation
    model_deviance_raw <- sum((y - fitted_values)^2)
    model_deviance_scaled <- model_deviance_raw / scale_parameter
    
    # Alternative model deviance calculation using (n-df)*φ̂
    model_deviance_alt <- df_residual * scale_parameter
    
    # Extract deviances directly from model object
    model_deviance_mgcv <- deviance(model)
    null_deviance_mgcv <- null_model$null.deviance
    
    # If null.deviance is not available from null model, calculate it
    if (is.null(null_deviance_mgcv)) {
        null_deviance_mgcv <- deviance(null_model)
    }
    
    # Calculate deviance explained using different methods
    dev_explained_raw <- 1 - (model_deviance_raw / null_deviance_raw)
    dev_explained_scaled <- 1 - (model_deviance_scaled / null_deviance_scaled)
    dev_explained_mgcv <- 1 - (model_deviance_mgcv / null_deviance_mgcv)
    dev_explained_summary <- summary(model)$dev.expl
    
    # Compile results
    results <- list(
        # Basic model info
        n_observations = n,
        response_variable = response_var,
        
        # Degrees of freedom
        edf_total = edf_total,
        df_residual = df_residual,
        df_model = df_model,
        
        # Scale parameter
        scale_parameter = scale_parameter,
        
        # Means and fitted values
        y_mean = y_mean,
        y_range = range(y),
        fitted_range = range(fitted_values),
        
        # Sum of squares (raw, unscaled)
        ss_total = ss_total,
        ss_residual = ss_residual,
        ss_model = ss_model,
        
        # Null deviance calculations
        null_deviance_raw = null_deviance_raw,
        null_deviance_scaled = null_deviance_scaled,
        null_deviance_mgcv = null_deviance_mgcv,
        
        # Model deviance calculations
        model_deviance_raw = model_deviance_raw,
        model_deviance_scaled = model_deviance_scaled,
        model_deviance_alt = model_deviance_alt,
        model_deviance_mgcv = model_deviance_mgcv,
        
        # Deviance explained calculations
        dev_explained_raw = dev_explained_raw,
        dev_explained_scaled = dev_explained_scaled,
        dev_explained_mgcv = dev_explained_mgcv,
        dev_explained_summary = dev_explained_summary,
        
        # Verification checks
        scale_check = all.equal(model_deviance_raw / df_residual, scale_parameter),
        deviance_alt_check = all.equal(model_deviance_alt, model_deviance_mgcv)
    )
    
    return(results)
}

# Function to analyze deviance components across all models in the results
analyze_all_deviances <- function(results_list, data_type = "PCa") {
    set.seed(123)
    cat("Analyzing deviances for", data_type, "results\n")
    
    deviance_data <- list()
    
    for (sample_name in names(results_list)) {
        cat("Processing sample:", sample_name, "\n")
        sample_results <- results_list[[sample_name]]
        
        if (is.null(sample_results)) {
            cat("  No results for this sample\n")
            next
        }
        
        for (gene_set_name in names(sample_results)) {
            cat("  Processing gene set:", gene_set_name, "\n")
            gene_set_result <- sample_results[[gene_set_name]]
            
            if (is.null(gene_set_result) || is.null(gene_set_result$best_model) || is.null(gene_set_result$gam_data)) {
                cat("    No valid model or data for this gene set\n")
                next
            }
            
            tryCatch({
                # Extract deviance components
                deviance_components <- extract_deviance_components(
                    gene_set_result$best_model, 
                    gene_set_result$gam_data
                )
                
                # Add identifying information
                deviance_components$sample <- sample_name
                deviance_components$gene_set <- gene_set_name
                deviance_components$data_type <- data_type
                
                # Store in list
                key <- paste(sample_name, gene_set_name, sep = "_")
                deviance_data[[key]] <- deviance_components
                
                cat("    Successfully extracted deviance components\n")
                
            }, error = function(e) {
                cat("    Error extracting deviance components:", conditionMessage(e), "\n")
            })
        }
    }
    
    return(deviance_data)
}

# Function to create summary dataframe from deviance analysis
create_deviance_summary <- function(deviance_data) {
    set.seed(123)
    cat("Creating deviance summary dataframe\n")
    
    if (length(deviance_data) == 0) {
        cat("No deviance data to summarize\n")
        return(data.frame())
    }
    
    summary_rows <- lapply(names(deviance_data), function(key) {
        components <- deviance_data[[key]]
        
        data.frame(
            Sample = components$sample,
            Gene_Set = components$gene_set,
            Data_Type = components$data_type,
            N_Observations = components$n_observations,
            Response_Variable = components$response_variable,
            
            # Degrees of freedom
            EDF_Total = components$edf_total,
            DF_Residual = components$df_residual,
            DF_Model = components$df_model,
            
            # Scale parameter
            Scale_Parameter = components$scale_parameter,
            
            # Basic statistics
            Y_Mean = components$y_mean,
            Y_Min = components$y_range[1],
            Y_Max = components$y_range[2],
            
            # Sum of squares
            SS_Total = components$ss_total,
            SS_Residual = components$ss_residual,
            SS_Model = components$ss_model,
            
            # Null deviance (different calculations)
            Null_Deviance_Raw = components$null_deviance_raw,
            Null_Deviance_Scaled = components$null_deviance_scaled,
            Null_Deviance_mgcv = components$null_deviance_mgcv,
            
            # Model deviance (different calculations)
            Model_Deviance_Raw = components$model_deviance_raw,
            Model_Deviance_Scaled = components$model_deviance_scaled,
            Model_Deviance_Alt = components$model_deviance_alt,
            Model_Deviance_mgcv = components$model_deviance_mgcv,
            
            # Deviance explained (different calculations)
            Dev_Explained_Raw = components$dev_explained_raw,
            Dev_Explained_Scaled = components$dev_explained_scaled,
            Dev_Explained_mgcv = components$dev_explained_mgcv,
            Dev_Explained_Summary = components$dev_explained_summary,
            
            # Verification checks
            Scale_Check = as.character(components$scale_check),
            Deviance_Alt_Check = as.character(components$deviance_alt_check),
            
            stringsAsFactors = FALSE
        )
    })
    
    result_df <- do.call(rbind, summary_rows)
    cat("Created summary with", nrow(result_df), "rows\n")
    return(result_df)
}

# Function to create detailed explanation sheet
create_deviance_explanations <- function() {
    explanations <- data.frame(
        Component = c(
            "EDF_Total",
            "DF_Residual",
            "DF_Model",
            "Scale_Parameter",
            "Y_Mean",
            "Y_Min",
            "Y_Max",
            "SS_Total",
            "SS_Residual",
            "SS_Model",
            "Null_Deviance_Raw",
            "Null_Deviance_Scaled", 
            "Null_Deviance_mgcv",
            "Model_Deviance_Raw",
            "Model_Deviance_Scaled",
            "Model_Deviance_Alt",
            "Model_Deviance_mgcv",
            "Dev_Explained_Raw",
            "Dev_Explained_Scaled",
            "Dev_Explained_mgcv",
            "Dev_Explained_Summary",
            "Scale_Check",
            "Deviance_Alt_Check"
        ),
        Formula = c(
            "Sum of effective degrees of freedom",
            "n - df_model",
            "n - df_residual",
            "φ̂ = ∑(yi - f(xi))² / (n - df)",
            "Mean of response variable",
            "Minimum value of response variable",
            "Maximum value of response variable",
            "∑(yi - ȳ)²",
            "∑(yi - f(xi))²",
            "SS_Total - SS_Residual",
            "∑(yi - ȳ)²",
            "∑(yi - ȳ)² / φ̂",
            "mgcv's null deviance calculation",
            "∑(yi - f(xi))²", 
            "∑(yi - f(xi))² / φ̂",
            "(n - df) × φ̂",
            "mgcv's model deviance calculation",
            "1 - (Model_Deviance_Raw / Null_Deviance_Raw)",
            "1 - (Model_Deviance_Scaled / Null_Deviance_Scaled)",
            "1 - (Model_Deviance_mgcv / Null_Deviance_mgcv)",
            "mgcv's summary deviance explained",
            "Verification: φ̂ = Model_Deviance_Raw / DF_Residual",
            "Verification: Model_Deviance_Alt = Model_Deviance_mgcv"
        ),
        Description = c(
            "Total effective degrees of freedom used by smooth terms",
            "Residual degrees of freedom for error estimation",
            "Model degrees of freedom",
            "Scale parameter estimate (residual variance)",
            "Mean of the response variable values",
            "Minimum value in the response variable",
            "Maximum value in the response variable",
            "Total sum of squares (same as raw null deviance)",
            "Residual sum of squares (same as raw model deviance)",
            "Explained sum of squares by the model",
            "Raw null deviance: sum of squared deviations from mean",
            "Scaled null deviance: raw null deviance divided by scale parameter",
            "Null deviance as calculated by mgcv package",
            "Raw model deviance: sum of squared residuals from fitted model",
            "Scaled model deviance: raw model deviance divided by scale parameter", 
            "Alternative model deviance calculation: (n-df) times scale parameter",
            "Model deviance as calculated by mgcv package",
            "Deviance explained using raw (unscaled) deviances",
            "Deviance explained using scaled deviances (standard for Gaussian)",
            "Deviance explained using mgcv's deviance calculations",
            "Deviance explained from mgcv's model summary",
            "Logical check if scale parameter calculation is correct",
            "Logical check if alternative deviance calculation matches mgcv"
        ),
        stringsAsFactors = FALSE
    )
    
    return(explanations)
}

# Function to export all deviance analysis to Excel
export_deviance_analysis <- function(pca_results, non_ca_results, filename) {
    set.seed(123)
    cat("Starting comprehensive deviance analysis export\n")
    
    # Analyze deviances for both datasets
    cat("Analyzing PCa deviances...\n")
    pca_deviances <- analyze_all_deviances(pca_results, "PCa")
    
    cat("Analyzing Non-Ca deviances...\n") 
    non_ca_deviances <- analyze_all_deviances(non_ca_results, "Non-Ca")
    
    # Combine all deviance data
    all_deviances <- c(pca_deviances, non_ca_deviances)
    
    if (length(all_deviances) == 0) {
        cat("No deviance data found to export\n")
        return(FALSE)
    }
    
    # Create summary dataframe
    deviance_summary <- create_deviance_summary(all_deviances)
    
    # Create explanations
    explanations <- create_deviance_explanations()
    
    # Prepare export list
    export_list <- list(
        Deviance_Summary = deviance_summary,
        Explanations = explanations
    )
    
    # Add individual sample details if requested
    cat("Creating detailed breakdown by sample and gene set...\n")
    detailed_breakdown <- lapply(names(all_deviances), function(key) {
        components <- all_deviances[[key]]
        
        # Create a detailed breakdown for this specific model
        detail_df <- data.frame(
            Metric = names(components),
            Value = sapply(components, function(x) {
                if (is.numeric(x) && length(x) == 1) {
                    as.character(x)
                } else if (is.character(x) && length(x) == 1) {
                    x
                } else if (length(x) == 2) {  # For ranges
                    paste(x, collapse = " to ")
                } else {
                    as.character(x[1])
                }
            }),
            Sample = components$sample,
            Gene_Set = components$gene_set,
            Data_Type = components$data_type,
            stringsAsFactors = FALSE
        )
        
        return(detail_df)
    })
    
    if (length(detailed_breakdown) > 0) {
        all_details <- do.call(rbind, detailed_breakdown)
        export_list$Detailed_Breakdown <- all_details
    }
    
    # Write to Excel
    cat("Writing to Excel file:", filename, "\n")
    tryCatch({
        write_xlsx(export_list, path = filename)
        cat("Deviance analysis exported successfully to", filename, "\n")
        
        # Print summary statistics
        cat("\nSummary of exported data:\n")
        cat("- Total models analyzed:", nrow(deviance_summary), "\n")
        cat("- PCa models:", sum(deviance_summary$Data_Type == "PCa"), "\n") 
        cat("- Non-Ca models:", sum(deviance_summary$Data_Type == "Non-Ca"), "\n")
        cat("- Gene sets analyzed:", length(unique(deviance_summary$Gene_Set)), "\n")
        cat("- Samples analyzed:", length(unique(deviance_summary$Sample)), "\n")
        
        return(TRUE)
    }, error = function(e) {
        cat("Error writing to Excel:", conditionMessage(e), "\n")
        return(FALSE)
    })
}

# Usage
cat("\n======== Starting Comprehensive Deviance Analysis ========\n")

# Run the comprehensive deviance analysis and export to Excel
tryCatch({
    success <- export_deviance_analysis(
        pca_results, 
        non_ca_results, 
        "GAM_Deviance_Analysis_Complete.xlsx"
    )
    
    if (success) {
        cat("Comprehensive deviance analysis completed successfully!\n")
    } else {
        cat("Deviance analysis encountered issues.\n")
    }
}, error = function(e) {
    cat("Error in comprehensive deviance analysis:", conditionMessage(e), "\n")
})

# Detailed analysis for all samples and gene sets
cat("\n======== Detailed Analysis for All Models ========\n")

# Initialize storage for detailed individual analyses
all_detailed_components <- list()
analysis_counter <- 0

# Process all PCa samples and gene sets
cat("Processing PCa samples...\n")
for (sample_name in pca_samples) {
    cat("Processing PCa sample:", sample_name, "\n")
    
    if (!is.null(pca_results[[sample_name]])) {
        for (gene_set_name in names(gene_sets)) {
            cat("  Processing gene set:", gene_set_name, "\n")
            
            if (!is.null(pca_results[[sample_name]][[gene_set_name]])) {
                model <- pca_results[[sample_name]][[gene_set_name]]$best_model
                data <- pca_results[[sample_name]][[gene_set_name]]$gam_data
                
                if (!is.null(model) && !is.null(data)) {
                    tryCatch({
                        detailed_components <- extract_deviance_components(model, data)
                        detailed_components$sample <- sample_name
                        detailed_components$gene_set <- gene_set_name
                        detailed_components$data_type <- "PCa"
                        
                        analysis_counter <- analysis_counter + 1
                        key <- paste("PCa", sample_name, gene_set_name, sep="_")
                        all_detailed_components[[key]] <- detailed_components
                        
                        cat("    Successfully analyzed - Scale parameter:", 
                            round(detailed_components$scale_parameter, 6), 
                            "Dev explained:", round(detailed_components$dev_explained_summary, 4), "\n")
                    }, error = function(e) {
                        cat("    Error:", conditionMessage(e), "\n")
                    })
                } else {
                    cat("    No valid model or data\n")
                }
            } else {
                cat("    No results for this gene set\n")
            }
        }
    } else {
        cat("  No results for this sample\n")
    }
}

# Process all Non-Ca samples and gene sets
cat("\nProcessing Non-Ca samples...\n")
for (sample_name in non_ca_samples) {
    cat("Processing Non-Ca sample:", sample_name, "\n")
    
    if (!is.null(non_ca_results[[sample_name]])) {
        for (gene_set_name in names(gene_sets)) {
            cat("  Processing gene set:", gene_set_name, "\n")
            
            if (!is.null(non_ca_results[[sample_name]][[gene_set_name]])) {
                model <- non_ca_results[[sample_name]][[gene_set_name]]$best_model
                data <- non_ca_results[[sample_name]][[gene_set_name]]$gam_data
                
                if (!is.null(model) && !is.null(data)) {
                    tryCatch({
                        detailed_components <- extract_deviance_components(model, data)
                        detailed_components$sample <- sample_name
                        detailed_components$gene_set <- gene_set_name
                        detailed_components$data_type <- "Non-Ca"
                        
                        analysis_counter <- analysis_counter + 1
                        key <- paste("NonCa", sample_name, gene_set_name, sep="_")
                        all_detailed_components[[key]] <- detailed_components
                        
                        cat("    Successfully analyzed - Scale parameter:", 
                            round(detailed_components$scale_parameter, 6), 
                            "Dev explained:", round(detailed_components$dev_explained_summary, 4), "\n")
                    }, error = function(e) {
                        cat("    Error:", conditionMessage(e), "\n")
                    })
                } else {
                    cat("    No valid model or data\n")
                }
            } else {
                cat("    No results for this gene set\n")
            }
        }
    } else {
        cat("  No results for this sample\n")
    }
}

cat("\nCompleted detailed analysis for", analysis_counter, "models\n")

# Create comprehensive summary statistics
if (length(all_detailed_components) > 0) {
    cat("\n======== Summary Statistics Across All Models ========\n")
    
    # Extract key metrics for summary
    scale_params <- sapply(all_detailed_components, function(x) x$scale_parameter)
    dev_explained <- sapply(all_detailed_components, function(x) x$dev_explained_summary)
    edf_total <- sapply(all_detailed_components, function(x) x$edf_total)
    
    cat("Scale parameter φ̂ - Min:", round(min(scale_params, na.rm=TRUE), 6), 
        "Max:", round(max(scale_params, na.rm=TRUE), 6), 
        "Mean:", round(mean(scale_params, na.rm=TRUE), 6), "\n")
    
    cat("Deviance explained - Min:", round(min(dev_explained, na.rm=TRUE), 4), 
        "Max:", round(max(dev_explained, na.rm=TRUE), 4), 
        "Mean:", round(mean(dev_explained, na.rm=TRUE), 4), "\n")
    
    cat("Total EDF - Min:", round(min(edf_total, na.rm=TRUE), 3), 
        "Max:", round(max(edf_total, na.rm=TRUE), 3), 
        "Mean:", round(mean(edf_total, na.rm=TRUE), 3), "\n")
    
    # Count by data type
    pca_count <- sum(sapply(all_detailed_components, function(x) x$data_type == "PCa"))
    non_ca_count <- sum(sapply(all_detailed_components, function(x) x$data_type == "Non-Ca"))
    
    cat("Models analyzed - PCa:", pca_count, "Non-Ca:", non_ca_count, "Total:", length(all_detailed_components), "\n")
    
    # Count by gene set
    gene_set_counts <- table(sapply(all_detailed_components, function(x) x$gene_set))
    cat("Models per gene set:\n")
    for (gs in names(gene_set_counts)) {
        cat("  ", gs, ":", gene_set_counts[gs], "\n")
    }
} else {
    cat("No detailed components were successfully extracted\n")
}

cat("\n======== Analysis Summary ========\n")
cat("The deviance analysis extracts and compares multiple calculation methods:\n")
cat("1. Raw deviances: Direct sum of squares without scaling\n")
cat("2. Scaled deviances: Divided by scale parameter φ̂ (standard for Gaussian)\n") 
cat("3. mgcv deviances: As calculated by the mgcv package\n")
cat("4. Alternative calculations: Using (n-df)*φ̂ formula\n")
cat("\nFor Gaussian family GAMs:\n")
cat("- Raw deviance = Sum of squares\n")
cat("- Scaled deviance = Sum of squares / φ̂\n")
cat("- Scale parameter φ̂ = RSS / (n - df)\n")
cat("- Standard deviance explained uses scaled deviances\n")
cat("\nAll results exported to: GAM_Deviance_Analysis_Complete.xlsx\n")