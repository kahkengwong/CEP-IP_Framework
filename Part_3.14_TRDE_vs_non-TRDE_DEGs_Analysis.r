######################################################################################################
# Part 3.14: DEGs of TRDE vs non-TRDE in Pre-IP or Post-IP (for subsequent GO analysis using ToppGene)
######################################################################################################
library(future)
library(future.apply)
library(parallel)
library(dplyr)
library(writexl)
library(mgcv)

# Set up parallel processing with memory optimization for 23 CPU cores (Intel i9-14900KF)
n_cores <- max(1, min(detectCores() - 1, 23))
cat("Setting up parallel processing with", n_cores, "cores (memory optimized)\n")

# Configure future with memory limits (64GB DDR5 RAM)
options(future.globals.maxSize = 45000 * 1024^3)  # 45GB RAM
plan(multisession, workers = n_cores)

# Function to format DEG results for export with UP/DOWN separation
format_degs_for_export <- function(de_results, comparison_name) {
    if (is.null(de_results) || nrow(de_results) == 0) {
        # Create empty structure with UP/DOWN sections
        empty_df <- data.frame(
            Gene = c("No upregulated genes found", "", "No downregulated genes found"),
            log2FC = c(NA, NA, NA),
            Linear_FC = c(NA, NA, NA),
            p_value = c(NA, NA, NA),
            BH_corrected_p = c(NA, NA, NA),
            stringsAsFactors = FALSE
        )
        colnames(empty_df) <- c("UPREGULATED_Gene", "UPREGULATED_log2FC", "UPREGULATED_Linear_FC", 
                                "UPREGULATED_p_value", "UPREGULATED_BH_corrected_p")
        
        # Add empty column and downregulated columns
        empty_df$` ` <- ""  # Empty separator column
        empty_df$DOWNREGULATED_Gene <- c("", "", "No downregulated genes found")
        empty_df$DOWNREGULATED_log2FC <- NA
        empty_df$DOWNREGULATED_Linear_FC <- NA
        empty_df$DOWNREGULATED_p_value <- NA
        empty_df$DOWNREGULATED_BH_corrected_p <- NA
        
        return(empty_df)
    }
    
    # Add gene names and calculate linear FC
    de_results$Gene <- rownames(de_results)
    de_results$Linear_FC <- 2^de_results$avg_log2FC
    de_results$BH_corrected_p <- p.adjust(de_results$p_val, method = "BH")
    
    # Separate upregulated and downregulated genes
    up_genes <- de_results[de_results$avg_log2FC > 0.1 & de_results$p_val < 0.05, ] %>%
        dplyr::select(Gene, avg_log2FC, Linear_FC, p_val, BH_corrected_p) %>%
        rename(log2FC = avg_log2FC, p_value = p_val) %>%
        arrange(p_value)
    
    down_genes <- de_results[de_results$avg_log2FC < -0.1 & de_results$p_val < 0.05, ] %>%
        dplyr::select(Gene, avg_log2FC, Linear_FC, p_val, BH_corrected_p) %>%
        rename(log2FC = avg_log2FC, p_value = p_val) %>%
        arrange(p_value)
    
    # Determine maximum number of rows needed
    max_rows <- max(nrow(up_genes), nrow(down_genes))
    
    # Create final data frame with UP and DOWN sections
    final_df <- data.frame(
        UPREGULATED_Gene = character(max_rows),
        UPREGULATED_log2FC = numeric(max_rows),
        UPREGULATED_Linear_FC = numeric(max_rows),
        UPREGULATED_p_value = numeric(max_rows),
        UPREGULATED_BH_corrected_p = numeric(max_rows),
        ` ` = rep("", max_rows),  # Empty separator column
        DOWNREGULATED_Gene = character(max_rows),
        DOWNREGULATED_log2FC = numeric(max_rows),
        DOWNREGULATED_Linear_FC = numeric(max_rows),
        DOWNREGULATED_p_value = numeric(max_rows),
        DOWNREGULATED_BH_corrected_p = numeric(max_rows),
        stringsAsFactors = FALSE
    )
    
    # Fill upregulated genes
    if (nrow(up_genes) > 0) {
        final_df[1:nrow(up_genes), 1:5] <- up_genes
    }
    
    # Fill downregulated genes
    if (nrow(down_genes) > 0) {
        final_df[1:nrow(down_genes), 7:11] <- down_genes
    }
    
    return(final_df)
}

# Modified function to perform additional comparisons
perform_additional_deg_analysis <- function(current_sample) {
    cat("Processing additional DEG comparisons for:", current_sample, "\n")
    
    # Get PCa cells (reuse existing filtering logic)
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(current_sample, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    sample_subset <- subset(prostate_ca_seurat, cells = selected_cells)
    
    # Add metadata (reuse existing logic)
    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    sample_data <- pca_results[[current_sample]][["Ribo"]]$gam_data
    model_dev_explained <- summary(sample_gam_model)$dev.expl
    
    common_cells <- intersect(rownames(sample_data), colnames(sample_subset))
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    
    # Calculate explanatory power
    null_model <- gam(Expression ~ 1, data = sample_data)
    explanatory_power <- 1 - ((sample_data$Expression - fitted(sample_gam_model))^2 / 
                                  (sample_data$Expression - fitted(null_model))^2)
    
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    sorted_cell_names <- rownames(sample_data)[sorted_indices]
    target_cells <- round(nrow(sample_data) * model_dev_explained)
    
    deviance_cells <- sorted_cell_names[1:target_cells]
    non_deviance_cells <- sorted_cell_names[(target_cells+1):length(sorted_cell_names)]
    
    # Add metadata
    sample_subset@meta.data$dev_explained_group <- "unknown"
    sample_subset@meta.data$dev_explained_group[colnames(sample_subset) %in% deviance_cells] <- "dev_explained"
    sample_subset@meta.data$dev_explained_group[colnames(sample_subset) %in% non_deviance_cells] <- "non_dev_explained"
    
    # Add IP classification
    ip_value <- inflection_points[[current_sample]]
    trpm4_scores <- sample_data$TRPM4[match(colnames(sample_subset), rownames(sample_data))]
    sample_subset@meta.data$TRPM4_Score <- trpm4_scores
    sample_subset@meta.data$ip_group <- ifelse(trpm4_scores < ip_value, "pre_IP", "post_IP")
    sample_subset@meta.data$combined_group <- paste(sample_subset@meta.data$dev_explained_group, 
                                                    sample_subset@meta.data$ip_group, sep = "_")
    
    # Define new comparisons
    new_comparisons <- list(
        list(name = "postIP_dev_explained_vs_postIP_non_dev_explained",
             ident1 = "dev_explained_post_IP",
             ident2 = "non_dev_explained_post_IP"),
        list(name = "postIP_dev_explained_vs_ALL_non_dev_explained", 
             ident1 = "dev_explained_post_IP",
             ident2 = c("non_dev_explained_pre_IP", "non_dev_explained_post_IP")),
        list(name = "postIP_dev_explained_vs_preIP_dev_explained",
             ident1 = "dev_explained_post_IP",
             ident2 = "dev_explained_pre_IP"),
        list(name = "ALL_dev_explained_vs_ALL_non_dev_explained",
             ident1 = c("dev_explained_pre_IP", "dev_explained_post_IP"),
             ident2 = c("non_dev_explained_pre_IP", "non_dev_explained_post_IP")),
        list(name = "preIP_dev_explained_vs_preIP_non_dev_explained",
             ident1 = "dev_explained_pre_IP",
             ident2 = "non_dev_explained_pre_IP")
    )
    
    # Handle the special case for combining groups
    for (i in seq_along(new_comparisons)) {
        comp <- new_comparisons[[i]]
        if (length(comp$ident1) > 1 || length(comp$ident2) > 1) {
            # For comparisons involving combined groups
            new_comparisons[[i]]$special_case <- "combine_groups"
        }
    }
    
    results <- list()
    Idents(sample_subset) <- "combined_group"
    
    for (comp in new_comparisons) {
        if (!is.null(comp$special_case) && comp$special_case == "combine_groups") {
            # Handle combined group comparison
            group1_count <- if (length(comp$ident1) > 1) {
                sum(sample_subset@meta.data$combined_group %in% comp$ident1)
            } else {
                sum(sample_subset@meta.data$combined_group == comp$ident1)
            }
            
            group2_count <- if (length(comp$ident2) > 1) {
                sum(sample_subset@meta.data$combined_group %in% comp$ident2)
            } else {
                sum(sample_subset@meta.data$combined_group == comp$ident2)
            }
            
            cat("Comparison:", comp$name, "- Group counts:", group1_count, "vs", group2_count, "\n")
            
            if (group1_count < 3 || group2_count < 3) {
                cat("Insufficient cells. Skipping.\n")
                results[[comp$name]] <- NULL
                next
            }
            
            # Create temporary combined identities
            temp_ident <- sample_subset@meta.data$combined_group
            
            if (length(comp$ident1) > 1) {
                temp_ident[temp_ident %in% comp$ident1] <- "combined_group1"
                ident1_name <- "combined_group1"
            } else {
                ident1_name <- comp$ident1
            }
            
            if (length(comp$ident2) > 1) {
                temp_ident[temp_ident %in% comp$ident2] <- "combined_group2"
                ident2_name <- "combined_group2"
            } else {
                ident2_name <- comp$ident2
            }
            
            sample_subset@meta.data$temp_combined <- temp_ident
            
            de_results <- tryCatch({
                Idents(sample_subset) <- "temp_combined"
                FindMarkers(sample_subset, 
                            ident.1 = ident1_name, 
                            ident.2 = ident2_name,
                            min.pct = 0.1,
                            logfc.threshold = 0.1,
                            test.use = "wilcox",
                            verbose = FALSE)
            }, error = function(e) {
                cat("Error:", conditionMessage(e), "\n")
                return(NULL)
            })
            
            # Reset identity for next comparison
            Idents(sample_subset) <- "combined_group"
            
        } else {
            # Handle regular single-group comparisons
            group1_count <- sum(sample_subset@meta.data$combined_group == comp$ident1)
            group2_count <- sum(sample_subset@meta.data$combined_group == comp$ident2)
            
            cat("Comparison:", comp$name, "- Group counts:", group1_count, "vs", group2_count, "\n")
            
            if (group1_count < 3 || group2_count < 3) {
                cat("Insufficient cells. Skipping.\n")
                results[[comp$name]] <- NULL
                next
            }
            
            de_results <- tryCatch({
                FindMarkers(sample_subset, 
                            ident.1 = comp$ident1, 
                            ident.2 = comp$ident2,
                            min.pct = 0.1,
                            logfc.threshold = 0.1,
                            test.use = "wilcox",
                            verbose = FALSE)
            }, error = function(e) {
                cat("Error:", conditionMessage(e), "\n")
                return(NULL)
            })
        }
        
        results[[comp$name]] <- de_results
    }
    
    return(results)
}

# Process all samples for additional DEG analysis
cat("Running additional DEG analysis for all samples...\n")

all_deg_results <- future_lapply(tumor_samples, function(sample_name) {
    perform_additional_deg_analysis(sample_name)
}, future.seed = TRUE)

names(all_deg_results) <- tumor_samples

# Export DEG results to Excel
cat("Exporting DEG results to Excel...\n")

comparison_names <- c("postIP_dev_explained_vs_postIP_non_dev_explained", 
                      "postIP_dev_explained_vs_ALL_non_dev_explained",
                      "postIP_dev_explained_vs_preIP_dev_explained",
                      "ALL_dev_explained_vs_ALL_non_dev_explained",
                      "preIP_dev_explained_vs_preIP_non_dev_explained")

for (comp_name in comparison_names) {
    excel_list <- list()
    
    for (sample_name in tumor_samples) {
        if (!is.null(all_deg_results[[sample_name]][[comp_name]])) {
            excel_list[[sample_name]] <- format_degs_for_export(
                all_deg_results[[sample_name]][[comp_name]], 
                comp_name
            )
        } else {
            excel_list[[sample_name]] <- format_degs_for_export(NULL, comp_name)
        }
    }
    
    filename <- paste0("DEGs_", comp_name, ".xlsx")
    writexl::write_xlsx(excel_list, filename)
    cat("Exported:", filename, "\n")
}

# Print summary
cat("\nDEG Analysis Summary:\n")
for (sample_name in tumor_samples) {
    cat(sample_name, ":\n")
    for (comp_name in comparison_names) {
        if (!is.null(all_deg_results[[sample_name]][[comp_name]])) {
            n_degs <- nrow(all_deg_results[[sample_name]][[comp_name]])
            cat("  ", comp_name, ":", n_degs, "DEGs\n")
        } else {
            cat("  ", comp_name, ": 0 DEGs\n")
        }
    }
}






