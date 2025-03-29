#############################################################
# Part 3.4: GAM-PRSS-REML Plots and EDF Analysis
##############################################################

# =========================================
# 1. PRSS and REML Iterations Plots 
# =========================================
# Load required libraries for visualization and data processing
library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(gridExtra)

# Read and prepare data from Excel files for plotting
prepare_data <- function(file_path) {
    # Read PRSS data from specified Excel file
    prss_data <- read_excel(file_path, sheet = "PRSS_Data")
    
    # Attempt to read REML summary data with error handling
    reml_summary <- tryCatch({
        read_excel(file_path, sheet = "REML_Summary")
    }, error = function(e) {
        message("No REML_Summary sheet found in ", file_path)
        NULL
    })
    
    # Attempt to read REML iterations data with error handling
    reml_iterations <- tryCatch({
        read_excel(file_path, sheet = "REML_Iterations")
    }, error = function(e) {
        message("No REML_Iterations sheet found in ", file_path)
        NULL
    })
    
    # Return list containing all prepared data
    return(list(
        prss = prss_data,
        reml_summary = reml_summary,
        reml_iterations = reml_iterations
    ))
}

# Visualize REML optimization for a specific PRSS iteration
reml_optimization_visualization <- function(data, prss_iteration, sample_id = NULL, gene_set = NULL) {
    # Extract REML iterations data from input
    reml_iterations <- data$reml_iterations
    
    # Return blank plot if REML iterations data is unavailable
    if (is.null(reml_iterations)) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "REML iterations data not available") + 
                 theme_void())
    }
    
    # Filter data by sample and gene set if specified
    if (!is.null(sample_id) && !is.null(gene_set)) {
        reml_iterations <- reml_iterations %>% 
            filter(Sample == sample_id, Gene_Set == gene_set)
    }
    
    # Filter REML data for the specified PRSS iteration
    reml_for_iteration <- reml_iterations %>% 
        filter(PRSS_Iteration == prss_iteration)
    
    # Return blank plot if no REML data exists for this iteration
    if (nrow(reml_for_iteration) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, 
                          label = paste("No REML data for PRSS iteration", prss_iteration)) + 
                 theme_void())
    }
    
    # Sort REML data by iteration number
    reml_for_iteration <- reml_for_iteration %>% 
        arrange(iteration)
    
    # Create line plot of REML scores across iterations
    p <- ggplot(reml_for_iteration, aes(x = iteration, y = score)) +
        geom_line(size = 1, color = "gray30") +
        geom_point(size = 3, color = "gray15") +
        labs(
            title = paste("REML Optimization for PRSS Iteration", prss_iteration),
            subtitle = paste0(
                if (!is.null(sample_id)) paste0("Sample: ", sample_id, " ") else "",
                if (!is.null(gene_set)) paste0("Gene Set: ", gene_set) else ""
            ),
            x = "REML Iteration",
            y = "REML Score"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 10),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95")
        )
    
    return(p)
}

# Generate all PRSS and REML visualizations for PCa and non-Ca data
generate_visualizations <- function(pca_file, nonca_file, output_prefix = "gam_optimization") {
    # Prepare data for PCa and non-Ca samples
    pca_data <- prepare_data(pca_file)
    nonca_data <- prepare_data(nonca_file)
    
    # Extract unique sample-gene set pairs from data
    get_pairs <- function(data) {
        if (nrow(data$prss) == 0) return(data.frame(Sample = character(), Gene_Set = character()))
        unique(data$prss[, c("Sample", "Gene_Set")])
    }
    
    pca_pairs <- get_pairs(pca_data)
    nonca_pairs <- get_pairs(nonca_data)
    
    # Initialize lists to store plots
    pca_plot_list <- list()
    nonca_plot_list <- list()
    
    # Process visualizations for PCa data
    for (i in 1:nrow(pca_pairs)) {
        sample_id <- pca_pairs$Sample[i]
        gene_set <- pca_pairs$Gene_Set[i]
        
        # Filter PRSS data for current sample and gene set
        prss_data <- pca_data$prss %>% 
            filter(Sample == sample_id, Gene_Set == gene_set)
        
        # Identify best PRSS iteration based on model flag or minimum PRSS
        best_iter <- prss_data %>% 
            filter(Is_Best_Model == TRUE) %>% 
            select(Iteration) %>% 
            pull()
        
        if (length(best_iter) == 0) {
            best_iter <- which.min(prss_data$PRSS)
        }
        
        # Create PRSS line plot
        prss_line <- ggplot(prss_data, aes(x = Iteration, y = PRSS)) +
            geom_line(size = 1, color = "gray30") +
            geom_point(size = 3, color = "gray15") +
            labs(
                title = "PRSS Values Across Iterations",
                subtitle = paste0("Sample: ", sample_id, ", Gene Set: ", gene_set),
                x = "Iteration",
                y = "PRSS Value"
            ) +
            theme_minimal() +
            theme(
                legend.position = "none",
                plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 10),
                axis.title = element_text(size = 9),
                axis.text = element_text(size = 8),
                panel.grid.major = element_line(color = "gray90"),
                panel.grid.minor = element_line(color = "gray95")
            )
        
        # Generate REML plot for best PRSS iteration
        reml_plot <- reml_optimization_visualization(pca_data, best_iter, sample_id, gene_set)
        
        # Store plots in list
        pca_plot_list[[length(pca_plot_list) + 1]] <- list(
            sample_id = sample_id, 
            gene_set = gene_set,
            prss_plot = prss_line,
            reml_plot = reml_plot
        )
        
        # Display individual plots for interactive viewing
        print(prss_line)
        print(reml_plot)
        
        # Combine PRSS and REML plots for display
        grid.arrange(
            prss_line, 
            reml_plot,
            ncol = 1
        )
    }
    
    # Process visualizations for non-Ca data
    for (i in 1:nrow(nonca_pairs)) {
        sample_id <- nonca_pairs$Sample[i]
        gene_set <- nonca_pairs$Gene_Set[i]
        
        # Filter PRSS data for current sample and gene set
        prss_data <- nonca_data$prss %>% 
            filter(Sample == sample_id, Gene_Set == gene_set)
        
        # Identify best PRSS iteration based on model flag or minimum PRSS
        best_iter <- prss_data %>% 
            filter(Is_Best_Model == TRUE) %>% 
            select(Iteration) %>% 
            pull()
        
        if (length(best_iter) == 0) {
            best_iter <- which.min(prss_data$PRSS)
        }
        
        # Create PRSS line plot
        prss_line <- ggplot(prss_data, aes(x = Iteration, y = PRSS)) +
            geom_line(size = 1, color = "gray30") +
            geom_point(size = 3, color = "gray15") +
            labs(
                title = "PRSS Values Across Iterations",
                subtitle = paste0("Sample: ", sample_id, ", Gene Set: ", gene_set),
                x = "Iteration",
                y = "PRSS Value"
            ) +
            theme_minimal() +
            theme(
                legend.position = "none",
                plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 10),
                axis.title = element_text(size = 9),
                axis.text = element_text(size = 8),
                panel.grid.major = element_line(color = "gray90"),
                panel.grid.minor = element_line(color = "gray95")
            )
        
        # Generate REML plot for best PRSS iteration
        reml_plot <- reml_optimization_visualization(nonca_data, best_iter, sample_id, gene_set)
        
        # Store plots in list
        nonca_plot_list[[length(nonca_plot_list) + 1]] <- list(
            sample_id = sample_id, 
            gene_set = gene_set,
            prss_plot = prss_line,
            reml_plot = reml_plot
        )
        
        # Display individual plots for interactive viewing
        print(prss_line)
        print(reml_plot)
        
        # Combine PRSS and REML plots for display
        grid.arrange(
            prss_line, 
            reml_plot,
            ncol = 1
        )
    }
    
    # Return lists of all generated plots
    return(list(
        pca_plots = pca_plot_list,
        nonca_plots = nonca_plot_list
    ))
}

# Save individual plots to PDF and JPG files
save_plot_to_files <- function(plot, base_filename, width = 8, height = 6) {
    # Save plot as PDF file
    pdf_filename <- paste0(base_filename, ".pdf")
    pdf(pdf_filename, width = width, height = height)
    print(plot)
    dev.off()
    cat("Saved plot to PDF:", pdf_filename, "\n")
    
    # Save plot as JPG file with high resolution
    jpg_filename <- paste0(base_filename, ".jpg")
    jpeg(jpg_filename, width = width * 100, height = height * 100, res = 300, quality = 90)
    print(plot)
    dev.off()
    cat("Saved plot to JPG:", jpg_filename, "\n")
}

# Generate and save all visualizations to files
generate_and_save_visualizations <- function(pca_file, nonca_file, output_dir = ".", file_prefix = "gam_optimization") {
    # Create output directory if it does not exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
    }
    
    # Generate visualizations using existing function
    results <- generate_visualizations(pca_file, nonca_file, file_prefix)
    
    # Extract plot lists from results
    pca_plots <- results$pca_plots
    nonca_plots <- results$nonca_plots
    
    # Save PCa plots to PDF and JPG files
    for (plot_group in pca_plots) {
        sample_id <- plot_group$sample_id
        gene_set <- plot_group$gene_set
        
        # Define base filenames for each plot type
        prss_base_filename <- file.path(output_dir, paste0(file_prefix, "_pca_", sample_id, "_", gene_set, "_prss"))
        reml_base_filename <- file.path(output_dir, paste0(file_prefix, "_pca_", sample_id, "_", gene_set, "_reml"))
        combined_base_filename <- file.path(output_dir, paste0(file_prefix, "_pca_", sample_id, "_", gene_set, "_combined"))
        
        # Save PRSS line plot
        save_plot_to_files(plot_group$prss_plot, prss_base_filename)
        
        # Save REML plot
        save_plot_to_files(plot_group$reml_plot, reml_base_filename)
        
        # Save combined PRSS and REML plot as PDF
        pdf(paste0(combined_base_filename, ".pdf"), width = 8, height = 10)
        grid.arrange(plot_group$prss_plot, plot_group$reml_plot, ncol = 1)
        dev.off()
        cat("Saved combined plot to PDF:", paste0(combined_base_filename, ".pdf"), "\n")
        
        # Save combined PRSS and REML plot as JPG
        jpeg(paste0(combined_base_filename, ".jpg"), width = 8 * 100, height = 10 * 100, res = 300, quality = 90)
        grid.arrange(plot_group$prss_plot, plot_group$reml_plot, ncol = 1)
        dev.off()
        cat("Saved combined plot to JPG:", paste0(combined_base_filename, ".jpg"), "\n")
    }
    
    # Save non-Ca plots to PDF and JPG files
    for (plot_group in nonca_plots) {
        sample_id <- plot_group$sample_id
        gene_set <- plot_group$gene_set
        
        # Define base filenames for each plot type
        prss_base_filename <- file.path(output_dir, paste0(file_prefix, "_nonca_", sample_id, "_", gene_set, "_prss"))
        reml_base_filename <- file.path(output_dir, paste0(file_prefix, "_nonca_", sample_id, "_", gene_set, "_reml"))
        combined_base_filename <- file.path(output_dir, paste0(file_prefix, "_nonca_", sample_id, "_", gene_set, "_combined"))
        
        # Save PRSS line plot
        save_plot_to_files(plot_group$prss_plot, prss_base_filename)
        
        # Save REML plot
        save_plot_to_files(plot_group$reml_plot, reml_base_filename)
        
        # Save combined PRSS and REML plot as PDF
        pdf(paste0(combined_base_filename, ".pdf"), width = 8, height = 10)
        grid.arrange(plot_group$prss_plot, plot_group$reml_plot, ncol = 1)
        dev.off()
        cat("Saved combined plot to PDF:", paste0(combined_base_filename, ".pdf"), "\n")
        
        # Save combined PRSS and REML plot as JPG
        jpeg(paste0(combined_base_filename, ".jpg"), width = 8 * 100, height = 10 * 100, res = 300, quality = 90)
        grid.arrange(plot_group$prss_plot, plot_group$reml_plot, ncol = 1)
        dev.off()
        cat("Saved combined plot to JPG:", paste0(combined_base_filename, ".jpg"), "\n")
    }
    
    # Confirm completion of plot saving
    cat("All plots saved to PDF and JPG files in directory:", output_dir, "\n")
}

# Execute visualization and saving for PCa and non-Ca data
generate_and_save_visualizations(
    "Ribo_AR_All_vs_Iteration_PCa_g1.5-all-lin.xlsx", 
    "Ribo_AR_All_vs_Iteration_NonCa_g1.5-all-lin.xlsx", 
    output_dir = "output_plots"
)

# Plot lambda values across iterations from optimization data
plot_lambda_values <- function(lambda_file, output_dir = "output_plots") {
    # Load required libraries for plotting
    library(readxl)
    library(ggplot2)
    
    # Create output directory if it does not exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        message("Created output directory: ", output_dir)
    }
    
    # Read lambda values from specified Excel file
    lambda_data <- read_excel(lambda_file)
    message("Successfully loaded data from: ", lambda_file)
    
    # Identify unique sample-gene set combinations
    sample_geneset_pairs <- unique(lambda_data[, c("sample", "gene_set")])
    message("Found ", nrow(sample_geneset_pairs), " sample-gene set combinations.")
    
    # Initialize list to store lambda plots
    all_plots <- list()
    
    # Generate plot for each sample-gene set combination
    for (i in 1:nrow(sample_geneset_pairs)) {
        sample_id <- sample_geneset_pairs$sample[i]
        gene_set <- sample_geneset_pairs$gene_set[i]
        
        message("Creating plot for ", sample_id, " - ", gene_set)
        
        # Filter and sort data for current sample and gene set
        this_data <- lambda_data[lambda_data$sample == sample_id & 
                               lambda_data$gene_set == gene_set, ]
        this_data <- this_data[order(this_data$iteration), ]
        
        # Create line plot of lambda values across iterations
        p <- ggplot(this_data, aes(x = iteration, y = lambda_value)) +
            geom_line(size = 1, color = "gray30") +
            geom_point(size = 3, color = "gray15") +
            labs(
                title = "Lambda Values Across Iterations",
                subtitle = paste0("Sample: ", sample_id, ", Gene Set: ", gene_set),
                x = "Iteration",
                y = "Lambda Value"
            ) +
            theme_minimal() +
            theme(
                legend.position = "none",
                plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 10),
                axis.title = element_text(size = 9),
                axis.text = element_text(size = 8),
                panel.grid.major = element_line(color = "gray90"),
                panel.grid.minor = element_line(color = "gray95")
            )
        
        # Store plot in list
        all_plots[[paste(sample_id, gene_set, sep = "_")]] <- p
        
        # Display plot for interactive viewing
        print(p)
    }
    
    # Return list of all lambda plots
    return(all_plots)
}

# Note export dimensions for PRSS and REML plots in SVG format
# Export PRSS plot as SVG at dimension 385 (width) x 310 (height)
# Export REML plot as SVG at dimension 365 (width) x 310 (height)


# ===================
# 2. GAM Plots
# ===================
# Create GAM plot for gene set expression versus TRPM4
create_gam_plot <- function(sample_results, sample_id) {
    tryCatch({
        # Set seed for reproducibility
        set.seed(123)
        
        # Return blank plot if sample results are insufficient
        if (is.null(sample_results) || length(sample_results) == 0) {
            return(ggplot() + 
                     annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient data for", sample_id)) + 
                     theme_void())
        }
        
        # Initialize empty dataframe for plot data
        plot_data <- data.frame()
        for (gene_set in c("Ribo", "AR", "PI3K_AKT", "mTOR", "GSK3B", "NFKB", "WNT")) {
            if (is.null(sample_results[[gene_set]])) {
                next
            }
            
            # Extract GAM data and best model for current gene set
            gam_data <- sample_results[[gene_set]]$gam_data
            best_model <- sample_results[[gene_set]]$best_model
            
            if (is.null(gam_data) || is.null(best_model)) {
                next
            }
            
            # Preserve original TRPM4 values
            gam_data$TRPM4_original <- gam_data$TRPM4
            
            gam_data$Gene_Set <- gene_set
            
            # Generate prediction data across TRPM4 range
            set.seed(123 + which(c("Ribo", "AR", "PI3K_AKT", "mTOR", "GSK3B", "NFKB", "WNT") == gene_set))
            pred_data <- data.frame(TRPM4 = seq(min(gam_data$TRPM4), max(gam_data$TRPM4), length.out = 1000))
            pred <- predict(best_model, newdata = pred_data, se.fit = TRUE)
            pred_data$fit <- pred$fit
            pred_data$se.fit <- pred$se.fit
            pred_data$Gene_Set <- gene_set
            
            # Prepare dataframes for actual and fitted values
            new_data_gam <- data.frame(
                TRPM4 = gam_data$TRPM4_original,
                Expression = gam_data$Expression,
                Gene_Set = gam_data$Gene_Set,
                se.fit = NA,
                type = "data",
                stringsAsFactors = FALSE
            )
            
            new_data_pred <- data.frame(
                TRPM4 = pred_data$TRPM4,
                Expression = pred_data$fit,
                Gene_Set = pred_data$Gene_Set,
                se.fit = pred_data$se.fit,
                type = "fit",
                stringsAsFactors = FALSE
            )
            
            # Combine actual and predicted data
            plot_data <- rbind(plot_data, new_data_gam, new_data_pred)
        }
        
        # Return blank plot if no valid data is collected
        if (nrow(plot_data) == 0) {
            return(ggplot() + 
                     annotate("text", x = 0.5, y = 0.5, label = paste("No valid data for", sample_id)) + 
                     theme_void())
        }
        
        # Define color scheme for gene sets
        gene_set_colors <- c("Ribo" = "#4B0082", "AR" = "#4F75DE", # Use "#1F5BFF" when plot individually
                             "PI3K_AKT" = "#D2B48C", "mTOR" = "#C0C0C0", "GSK3B" = "#D3D3A4",
                             "NFKB" = "#B0C4DE", "WNT" = "#C0A99E")
        
        # Set seed for jitter reproducibility
        set.seed(123)
        # Create GAM plot with points, fit lines, and confidence bands
        ggplot(plot_data, aes(x = TRPM4, y = Expression, color = Gene_Set)) +
            geom_point(data = subset(plot_data, type == "data"), 
                       position = position_jitter(width = 0.01, height = 0, seed = 123),
                       alpha = 0.2, size = 0.5) + # alpha = 0.3 for individual plot
            geom_line(data = subset(plot_data, type == "fit"), size = 1) +
            geom_ribbon(data = subset(plot_data, type == "fit"), 
                        aes(ymin = Expression - 1.96 * se.fit, ymax = Expression + 1.96 * se.fit, fill = Gene_Set), 
                        alpha = 0.3, color = NA) +
            scale_color_manual(values = gene_set_colors) +
            scale_fill_manual(values = gene_set_colors) +
            labs(x = "TRPM4 Expression",
                 y = "Gene Set Expression (Log2)",
                 title = paste("Gene Set Expression vs TRPM4 in", sample_id)) +
            theme_minimal() +
            theme(legend.position = "right",
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
                  axis.title = element_text(face = "bold", size = 8),
                  axis.text = element_text(size = 6),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.35),
                  plot.background = element_rect(fill = "white", color = NA),
                  panel.background = element_rect(fill = "white", color = NA)) +
            scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
            scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
    }, error = function(e) {
        cat("Error in create_gam_plot for sample", sample_id, ":\n")
        cat("Error message:", conditionMessage(e), "\n")
        print(traceback())
        return(NULL)
    })
}

# Create and save GAM plots for all samples
create_and_save_plots <- function(results, samples, output_prefix) {
    for (sample in samples) {
        if (!is.null(results[[sample]])) {
            cat("Processing sample:", sample, "\n")
            
            # Generate GAM plot for current sample
            gam_plot <- create_gam_plot(results[[sample]], sample)
            
            if (!is.null(gam_plot)) {
                # Save plot as PDF with specified dimensions
                filename <- paste0(output_prefix, "_", sample, ".pdf")
                ggsave(filename, gam_plot, width = 4.17, height = 3.07)
                cat("Plot saved for sample:", sample, "as", filename, "\n")
            } else {
                cat("No valid plot created for sample:", sample, "\n")
            }
        } else {
            cat("No results for sample:", sample, "\n")
        }
    }
}

# Generate and save GAM plots for PCa samples
create_and_save_plots(pca_results, pca_samples, "Ribo_AR_All_Plot_PCa_g1.5-all")

# Generate and save GAM plots for non-Ca samples
create_and_save_plots(non_ca_results, non_ca_samples, "Ribo_AR_All_Plot_NonCa_g1.5-all")


# ====================================
# 3. EDF and FDR Testing for Smooth Terms
# ====================================
# Extract EDF summary and apply FDR correction for all samples
extract_edf_summary <- function(results_list) {
    # Initialize dataframe to store EDF results
    edf_summary <- data.frame(
        Sample = character(),
        Gene_Set = character(),
        K_Value = numeric(),
        Total_EDF = numeric(),
        Smooth_EDF = numeric(),
        Parametric_DF = numeric(),
        Smooth_Term_P_Value = numeric(),
        Deviance_Explained = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Process each sample in the results list
    for (sample_name in names(results_list)) {
        cat("Processing sample:", sample_name, "\n")
        sample_results <- results_list[[sample_name]]
        
        if (is.null(sample_results)) {
            cat("No results for this sample\n")
            next
        }
        
        # Process each gene set within the sample
        for (gene_set_name in names(sample_results)) {
            cat("Processing gene set:", gene_set_name, "\n")
            result <- sample_results[[gene_set_name]]
            
            if (is.null(result) || is.null(result$best_model)) {
                cat("No best model found\n")
                next
            }
            
            # Extract model and its summary
            model <- result$best_model
            model_summary <- summary(model)
            
            # Calculate smooth term EDF
            smooth_edf <- sum(model_summary$s.table[,"edf"])
            
            # Calculate parametric degrees of freedom
            parametric_df <- length(model$coefficients) - sum(grepl("s\\(TRPM4\\)", names(model$coefficients)))
            
            # Calculate total EDF
            total_edf <- sum(model_summary$edf)
            
            # Extract k value from best parameters or smooth term
            k_value <- NA
            if (!is.null(result$best_params) && "k" %in% names(result$best_params)) {
                k_value <- result$best_params$k
            } else if (length(model$smooth) > 0) {
                k_value <- model$smooth[[1]]$bs.dim
            }
            
            # Extract p-value for smooth term
            smooth_p_value <- NA
            if (nrow(model_summary$s.table) > 0) {
                smooth_p_value <- model_summary$s.table[,"p-value"]
            }
            
            # Extract deviance explained by model
            deviance_explained <- model_summary$dev.expl
            
            # Append results to summary dataframe
            edf_summary <- rbind(edf_summary, data.frame(
                Sample = sample_name,
                Gene_Set = gene_set_name,
                K_Value = k_value,
                Total_EDF = total_edf,
                Smooth_EDF = smooth_edf,
                Parametric_DF = parametric_df,
                Smooth_Term_P_Value = smooth_p_value,
                Deviance_Explained = deviance_explained,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    return(edf_summary)
}

# Extract EDF summary for PCa samples
pca_edf_summary <- extract_edf_summary(pca_results)
print(pca_edf_summary)

# Extract EDF summary for non-Ca samples
non_ca_edf_summary <- extract_edf_summary(non_ca_results)
print(non_ca_edf_summary)

# Apply FDR correction to p-values for PCa samples
if (nrow(pca_edf_summary) > 0) {
    pca_edf_summary$FDR_Adjusted_P_Value <- p.adjust(pca_edf_summary$Smooth_Term_P_Value, method = "BH")
}

# Apply FDR correction to p-values for non-Ca samples
if (nrow(non_ca_edf_summary) > 0) {
    non_ca_edf_summary$FDR_Adjusted_P_Value <- p.adjust(non_ca_edf_summary$Smooth_Term_P_Value, method = "BH")
}

# Identify significantly non-linear relationships in PCa samples
if (nrow(pca_edf_summary) > 0) {
    pca_edf_summary$Is_Significantly_Nonlinear <- pca_edf_summary$Smooth_EDF > 1 & 
                                                  pca_edf_summary$FDR_Adjusted_P_Value < 0.05
}

# Identify significantly non-linear relationships in non-Ca samples
if (nrow(non_ca_edf_summary) > 0) {
    non_ca_edf_summary$Is_Significantly_Nonlinear <- non_ca_edf_summary$Smooth_EDF > 1 & 
                                                     non_ca_edf_summary$FDR_Adjusted_P_Value < 0.05
}

# Combine PCa and non-Ca EDF summaries and export to Excel
all_edf_summary <- rbind(
    cbind(Group = "PCa", pca_edf_summary),
    cbind(Group = "Non-Ca", non_ca_edf_summary)
)

# Export EDF summary to Excel file
write_xlsx(list(EDF_Summary = all_edf_summary), path = "Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")

# Generate detailed EDF report with basis function counts
detailed_edf_report <- function(results_list) {
    # Initialize dataframe for detailed report
    detailed_report <- data.frame(
        Sample = character(),
        Gene_Set = character(),
        K_Value = numeric(),
        Total_EDF = numeric(),
        Smooth_EDF = numeric(),
        Parametric_DF = numeric(),
        Nonzero_Basis_Functions = numeric(), 
        Total_Basis_Functions = numeric(),
        Smooth_Coef_Count = numeric(),
        Smooth_Term_P_Value = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (sample_name in names(results_list)) {
        sample_results <- results_list[[sample_name]]
        
        if (is.null(sample_results)) next
        
        for (gene_set_name in names(sample_results)) {
            result <- sample_results[[gene_set_name]]
            
            if (is.null(result) || is.null(result$best_model)) next
            
            # Extract model and its summary
            model <- result$best_model
            model_summary <- summary(model)
            
            # Count basis functions and non-zero coefficients
            smooth_coefs <- coef(model)[grepl("s\\(TRPM4\\)", names(coef(model)))]
            nonzero_basis <- sum(abs(smooth_coefs) > 1e-8)
            total_basis <- length(smooth_coefs)
            
            # Extract k value from best parameters or smooth term
            k_value <- NA
            if (!is.null(result$best_params) && "k" %in% names(result$best_params)) {
                k_value <- result$best_params$k
            } else if (length(model$smooth) > 0) {
                k_value <- model$smooth[[1]]$bs.dim
            }
            
            # Extract p-value for smooth term
            smooth_p_value <- NA
            if (nrow(model_summary$s.table) > 0) {
                smooth_p_value <- model_summary$s.table[,"p-value"]
            }
            
            # Append results to detailed report
            detailed_report <- rbind(detailed_report, data.frame(
                Sample = sample_name,
                Gene_Set = gene_set_name,
                K_Value = k_value,
                Total_EDF = sum(model_summary$edf),
                Smooth_EDF = sum(model_summary$s.table[,"edf"]),
                Parametric_DF = length(model$coefficients) - length(smooth_coefs),
                Nonzero_Basis_Functions = nonzero_basis,
                Total_Basis_Functions = total_basis,
                Smooth_Coef_Count = length(smooth_coefs),
                Smooth_Term_P_Value = smooth_p_value,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    return(detailed_report)
}

# Generate detailed EDF reports for PCa samples
pca_detailed <- detailed_edf_report(pca_results)

# Generate detailed EDF reports for non-Ca samples
non_ca_detailed <- detailed_edf_report(non_ca_results)

# Apply FDR correction and non-linearity check to PCa detailed report
if (nrow(pca_detailed) > 0) {
    pca_detailed$FDR_Adjusted_P_Value <- p.adjust(pca_detailed$Smooth_Term_P_Value, method = "BH")
    pca_detailed$Is_Significantly_Nonlinear <- pca_detailed$Smooth_EDF > 1 & 
                                               pca_detailed$FDR_Adjusted_P_Value < 0.05
}

# Apply FDR correction and non-linearity check to non-Ca detailed report
if (nrow(non_ca_detailed) > 0) {
    non_ca_detailed$FDR_Adjusted_P_Value <- p.adjust(non_ca_detailed$Smooth_Term_P_Value, method = "BH")
    non_ca_detailed$Is_Significantly_Nonlinear <- non_ca_detailed$Smooth_EDF > 1 & 
                                                  non_ca_detailed$FDR_Adjusted_P_Value < 0.05
}

# Combine detailed reports for PCa and non-Ca samples
all_detailed <- rbind(
    cbind(Group = "PCa", pca_detailed),
    cbind(Group = "Non-Ca", non_ca_detailed)
)

# Add detailed report to existing Excel file or create new file
if (file.exists("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")) {
    # Read existing sheets from Excel file
    sheets <- readxl::excel_sheets("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
    sheet_data <- list()
    
    for (sheet in sheets) {
        sheet_data[[sheet]] <- readxl::read_excel("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx", sheet = sheet)
    }
    
    # Add detailed basis functions sheet
    sheet_data[["Detailed_Basis_Functions"]] <- all_detailed
    
    # Write updated sheets back to file
    write_xlsx(sheet_data, path = "Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
} else {
    write_xlsx(list(
        EDF_Summary = all_edf_summary,
        Detailed_Basis_Functions = all_detailed
    ), path = "Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
}

# Summarize non-linear relationships across samples
nonlinear_summary <- subset(all_edf_summary, Is_Significantly_Nonlinear == TRUE)
nonlinear_counts <- aggregate(Is_Significantly_Nonlinear ~ Group + Sample, data = all_edf_summary, FUN = sum)
colnames(nonlinear_counts)[3] <- "Nonlinear_Relationships_Count"
total_counts <- aggregate(Is_Significantly_Nonlinear ~ Group + Sample, data = all_edf_summary, FUN = length)
colnames(total_counts)[3] <- "Total_Relationships"
nonlinear_stats <- merge(nonlinear_counts, total_counts, by = c("Group", "Sample"))
nonlinear_stats$Nonlinear_Percentage <- (nonlinear_stats$Nonlinear_Relationships_Count / nonlinear_stats$Total_Relationships) * 100

# Add non-linear summary to Excel file
if (file.exists("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")) {
    # Read existing sheets from Excel file
    sheets <- readxl::excel_sheets("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
    sheet_data <- list()
    
    for (sheet in sheets) {
        sheet_data[[sheet]] <- readxl::read_excel("Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx", sheet = sheet)
    }
    
    # Add non-linear relationships and stats sheets
    sheet_data[["Nonlinear_Relationships"]] <- nonlinear_summary
    sheet_data[["Nonlinear_Stats"]] <- nonlinear_stats
    
    # Write updated sheets back to file
    write_xlsx(sheet_data, path = "Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
} else {
    write_xlsx(list(
        EDF_Summary = all_edf_summary,
        Detailed_Basis_Functions = all_detailed,
        Nonlinear_Relationships = nonlinear_summary,
        Nonlinear_Stats = nonlinear_stats
    ), path = "Ribo_AR_EDF_Summary_g1.5-all-lin.xlsx")
}

# Print summary statistics for EDF and basis functions
cat("\nSummary of EDF and basis functions:\n")
cat("-------------------------------------\n")

# Print summary for PCa samples
cat("PCa Samples:\n")
cat("Average K value:", mean(pca_detailed$K_Value, na.rm=TRUE), "\n")
cat("Average total EDF:", mean(pca_detailed$Total_EDF, na.rm=TRUE), "\n")
cat("Average smooth EDF:", mean(pca_detailed$Smooth_EDF, na.rm=TRUE), "\n")
cat("Average basis functions:", mean(pca_detailed$Total_Basis_Functions, na.rm=TRUE), "\n")
cat("Average non-zero basis functions:", mean(pca_detailed$Nonzero_Basis_Functions, na.rm=TRUE), "\n")
cat("Percentage of significantly non-linear relationships:", 
    sum(pca_detailed$Is_Significantly_Nonlinear, na.rm=TRUE) / nrow(pca_detailed) * 100, "%\n\n")

# Print summary for non-Ca samples
cat("Non-Ca Samples:\n")
cat("Average K value:", mean(non_ca_detailed$K_Value, na.rm=TRUE), "\n")
cat("Average total EDF:", mean(non_ca_detailed$Total_EDF, na.rm=TRUE), "\n")
cat("Average smooth EDF:", mean(non_ca_detailed$Smooth_EDF, na.rm=TRUE), "\n")
cat("Average basis functions:", mean(non_ca_detailed$Total_Basis_Functions, na.rm=TRUE), "\n")
cat("Average non-zero basis functions:", mean(non_ca_detailed$Nonzero_Basis_Functions, na.rm=TRUE), "\n")
cat("Percentage of significantly non-linear relationships:", 
    sum(non_ca_detailed$Is_Significantly_Nonlinear, na.rm=TRUE) / nrow(non_ca_detailed) * 100, "%\n")