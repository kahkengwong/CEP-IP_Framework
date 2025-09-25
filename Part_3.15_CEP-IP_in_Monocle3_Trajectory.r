##################################################################################################################
# Part 3.15: CEP-IP in Monocle3 Trajectory (Visualization for TREP and non-TREP Cells, and Quantitative Analysis)
##################################################################################################################
library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggridges)
library(openxlsx)
library(cluster)
library(FNN)
library(broom)

# Function to create Monocle3 trajectory visualization
create_monocle3_trajectory <- function(sample_id) {
    cat("=== CREATING MONOCLE3 TRAJECTORY FOR", sample_id, "===\n")
    
    # Get PCa cells for the sample
    pca_clusters <- c(6, 9, 11, 14, 19)
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = pca_clusters)
    sample_cells <- WhichCells(prostate_results$seurat_obj, 
                               cells = grep(sample_id, colnames(prostate_results$seurat_obj), value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    
    # Subset the Seurat object
    sample_seurat <- subset(prostate_ca_seurat, cells = selected_cells)
    
    cat("After subsetting:\n")
    cat("Selected cells:", length(selected_cells), "\n")
    cat("Seurat object cells:", ncol(sample_seurat), "\n")
    cat("First 5 selected cells:", head(selected_cells, 5), "\n")
    cat("First 5 seurat cells:", head(colnames(sample_seurat), 5), "\n")
    
    # Get GAM data for cell classifications
    original_gam_data <- pca_results[[sample_id]][["Ribo"]]$gam_data
    
    # Debug: Check if original_gam_data exists and matches the sample
    if (is.null(original_gam_data)) {
        stop("GAM data not found for sample ", sample_id)
    }
    
    cat("GAM data dimensions:", nrow(original_gam_data), "x", ncol(original_gam_data), "\n")
    
    # Check if GAM data is from the correct sample
    gam_sample_names <- rownames(original_gam_data)[1:5]
    cat("GAM data sample check - first 5 rownames:", gam_sample_names, "\n")
    
    # If GAM data doesn't match current sample, filter it
    if (!any(grepl(sample_id, gam_sample_names))) {
        cat("WARNING: GAM data doesn't match current sample. Filtering...\n")
        # Filter to only include cells from current sample
        matching_rows <- grepl(sample_id, rownames(original_gam_data))
        if (sum(matching_rows) == 0) {
            stop("No cells from sample ", sample_id, " found in GAM data")
        }
        original_gam_data <- original_gam_data[matching_rows, ]
        cat("Filtered GAM data to", nrow(original_gam_data), "cells from", sample_id, "\n")
    }
    
    # Get inflection point
    inflection_points <- c(3.800, 2.214, 3.179, 3.306, 2.636, 3.465, 3.476)
    patient_ids <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                     "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                     "HYW_5755_Tumor")
    names(inflection_points) <- patient_ids
    ip_value <- inflection_points[sample_id]
    
    # Calculate cell classifications (using original method)
    gam_model <- pca_results[[sample_id]][["Ribo"]]$best_model
    dev_explained <- summary(gam_model)$dev.expl
    
    # Calculate explanatory power using the original GAM data
    null_fitted <- mean(original_gam_data$Expression)
    gam_fitted <- fitted(gam_model)
    null_sq_diff <- (original_gam_data$Expression - null_fitted)^2
    gam_sq_diff <- (original_gam_data$Expression - gam_fitted)^2
    explanatory_power <- 1 - (gam_sq_diff / null_sq_diff)
    
    # Determine DE cells
    sorted_indices <- order(explanatory_power, decreasing = TRUE)
    target_cells <- round(nrow(original_gam_data) * dev_explained)
    
    # Ensure target_cells is at least 1 and not more than total cells
    target_cells <- max(1, min(target_cells, nrow(original_gam_data)))
    
    cat("Target DE cells:", target_cells, "out of", nrow(original_gam_data), "total cells\n")
    
    de_cells <- rownames(original_gam_data)[sorted_indices[1:target_cells]]
    
    # Create cell type classifications
    original_gam_data$is_de <- ifelse(rownames(original_gam_data) %in% de_cells, "DE", "Non-DE")
    original_gam_data$timing <- ifelse(original_gam_data$TRPM4 < ip_value, "Pre-IP", "Post-IP")
    
    # Create combined classification
    original_gam_data$cell_group <- paste0(original_gam_data$timing, "_", original_gam_data$is_de)
    
    # Match cells between Seurat and GAM data
    common_cells <- intersect(colnames(sample_seurat), rownames(original_gam_data))
    cat("Common cells found:", length(common_cells), "\n")
    
    if (length(common_cells) == 0) {
        stop("No common cells found between Seurat and GAM data")
    }
    
    # Filter both datasets to common cells
    sample_seurat <- subset(sample_seurat, cells = common_cells)
    sample_data_matched <- original_gam_data[common_cells, ]
    
    cat("Final processing:", nrow(sample_data_matched), "cells\n")
    cat("Cell group distribution:\n")
    print(table(sample_data_matched$cell_group))
    
    # Convert Seurat to Monocle3 cell_data_set
    cat("Converting to Monocle3 format...\n")
    
    # Extract count matrix
    count_matrix <- GetAssayData(sample_seurat, slot = "counts", assay = "RNA")
    
    # Create cell metadata
    cell_metadata <- data.frame(
        cell_id = colnames(sample_seurat),
        TRPM4 = sample_data_matched$TRPM4,
        Ribo = sample_data_matched$Expression,
        cell_group = sample_data_matched$cell_group,
        is_de = sample_data_matched$is_de,
        timing = sample_data_matched$timing,
        explanatory_power = explanatory_power[match(colnames(sample_seurat), rownames(original_gam_data))],
        stringsAsFactors = FALSE
    )
    rownames(cell_metadata) <- cell_metadata$cell_id
    
    # Create gene metadata
    gene_metadata <- data.frame(
        gene_id = rownames(count_matrix),
        gene_short_name = rownames(count_matrix),
        stringsAsFactors = FALSE
    )
    rownames(gene_metadata) <- gene_metadata$gene_id
    
    # Create Monocle3 cell_data_set
    cds <- new_cell_data_set(
        expression_data = count_matrix,
        cell_metadata = cell_metadata,
        gene_metadata = gene_metadata
    )
    
    # Preprocess the data
    cat("Preprocessing data...\n")
    cds <- preprocess_cds(cds, num_dim = 200, verbose = FALSE)
    
    # Reduce dimensions
    cat("Reducing dimensions...\n")
    cds <- reduce_dimension(cds, verbose = FALSE)
    
    # Cluster cells
    cat("Clustering cells...\n")
    cds <- cluster_cells(cds, verbose = FALSE, resolution = 0.0005, k = 250)
    
    # Learn trajectory with adjusted parameters for better extension
    cat("Learning trajectory...\n")
    cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE, close_loop = FALSE)
    
    # Find trajectory root automatically (cell with highest TRPM4 expression)
    trpm4_values <- colData(cds)$TRPM4
    root_cell <- colnames(cds)[which.max(trpm4_values)]
    
    # Order cells along trajectory
    cds <- order_cells(cds, root_cells = root_cell, verbose = FALSE)
    
    # Debug: Check cell group classification
    cat("Cell group breakdown:\n")
    print(table(sample_data_matched$cell_group))
    print(table(sample_data_matched$is_de, sample_data_matched$timing))
    
    # Define more distinct color palette for cell groups
    colors <- c(
        "Pre-IP_Non-DE" = "#CCCCCC",     # light gray
        "Post-IP_Non-DE" = "#666666",    # dark gray  
        "Pre-IP_DE" = "#DDA0DD",         # plum (light purple)
        "Post-IP_DE" = "#6666FF"         # changed to blue
    )
    
    # Create UMAP plot colored by cell groups
    cat("Creating trajectory plot...\n")
    
    trajectory_plot <- plot_cells(cds, 
                                  color_cells_by = "cell_group",
                                  label_cell_groups = FALSE,
                                  label_leaves = FALSE,
                                  label_branch_points = FALSE,
                                  label_roots = FALSE,
                                  show_trajectory_graph = TRUE,
                                  graph_label_size = 3,
                                  cell_size = 1.5,
                                  cell_stroke = 0.5,  # Add gray border to cells
                                  trajectory_graph_color = "#333333",  # Dark gray trajectory lines
                                  trajectory_graph_segment_size = 1) +  # Thicker trajectory lines
        scale_color_manual(values = colors) +
        labs(title = paste("Monocle3 Trajectory -", sample_id),
             subtitle = paste("DE cells (purple/blue) vs Non-DE cells (gray), Pre/Post IP =", round(ip_value, 2)),
             color = "Cell Group") +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10)
        )
    
    # Create pseudotime plot with plasma palette
    cat("Creating pseudotime plot...\n")
    
    pseudotime_plot <- plot_cells(cds, 
                                  color_cells_by = "pseudotime",
                                  label_cell_groups = FALSE,
                                  label_leaves = FALSE,
                                  label_branch_points = FALSE,
                                  label_roots = FALSE,
                                  show_trajectory_graph = TRUE,
                                  graph_label_size = 3,
                                  cell_size = 0.5,
                                  cell_stroke = 0.5,  # Add gray border to cells
                                  trajectory_graph_color = "#333333",  # Dark gray trajectory lines
                                  trajectory_graph_segment_size = 1) +  # Thicker trajectory lines
        scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
        labs(title = paste("Monocle3 Pseudotime -", sample_id),
             subtitle = paste("Lighter colors = higher pseudotime, Pre/Post IP =", round(ip_value, 2))) +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10)
        )
    
    # Print both plots
    print(trajectory_plot)
    print(pseudotime_plot)
    
    # Print summary statistics
    cat("\n=== TRAJECTORY SUMMARY ===\n")
    cat("Sample:", sample_id, "\n")
    cat("Total cells:", ncol(cds), "\n")
    cat("Inflection point:", ip_value, "\n")
    cat("Deviance explained:", round(dev_explained * 100, 1), "%\n")
    cat("Root cell:", root_cell, "\n")
    
    cat("\nCell group counts:\n")
    group_counts <- table(cell_metadata$cell_group)
    for (i in names(group_counts)) {
        cat("  ", i, ":", group_counts[i], "(", round(group_counts[i]/sum(group_counts)*100, 1), "%)\n")
    }
    
    # Return the cell_data_set for further analysis if needed
    return(cds)
}

# List of all tumor samples
tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", 
                   "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor", 
                   "HYW_5755_Tumor")

# Generate trajectory plots for all samples
cat("\n\nGenerating Monocle3 trajectory plots for all samples...\n")
cat(rep("=", 60), "\n")

# Create a list to store all trajectory results
all_trajectories <- list()

# Process each sample for trajectory plotting
for (sample_name in tumor_samples) {
    tryCatch({
        cds_result <- create_monocle3_trajectory(sample_name)
        if (!is.null(cds_result)) {
            all_trajectories[[sample_name]] <- cds_result
        }
    }, error = function(e) {
        cat("Error creating trajectory for", sample_name, ":", conditionMessage(e), "\n")
    })
}

cat("\nCompleted generating Monocle3 trajectory plots for all samples.\n")


#######################################################################################
# Quantitative Analysis of Pre-IP and Post-IP DE Cells in Monocle3 Trajectories 
#######################################################################################
# Load additional required libraries
library(ggridges)
library(openxlsx)
library(cluster)
library(FNN)
library(broom)

# Function to analyze cell type clustering and create ridgeline plots
analyze_cell_type_clustering <- function(cds, sample_id, colors) {
    cat("\n=== ANALYZING CELL TYPE CLUSTERING FOR", sample_id, "===\n")
    
    # Extract UMAP coordinates and cell metadata
    umap_coords <- reducedDims(cds)$UMAP
    cell_metadata <- colData(cds)
    
    # Create data frame for analysis
    analysis_data <- data.frame(
        cell_id = rownames(umap_coords),
        UMAP1 = umap_coords[,1],
        UMAP2 = umap_coords[,2],
        cell_group = cell_metadata$cell_group,
        is_de = cell_metadata$is_de,
        timing = cell_metadata$timing,
        TRPM4 = cell_metadata$TRPM4,
        Ribo = cell_metadata$Ribo,
        stringsAsFactors = FALSE
    )
    
    # Calculate distance matrix between all cells
    dist_matrix <- dist(umap_coords)
    
    # 1. WITHIN-GROUP vs BETWEEN-GROUP DISTANCE ANALYSIS
    cat("Calculating within-group and between-group distances...\n")
    
    distance_results <- list()
    cell_groups <- unique(analysis_data$cell_group)
    
    # Calculate within-group distances
    within_distances <- list()
    for (group in cell_groups) {
        group_cells <- which(analysis_data$cell_group == group)
        if (length(group_cells) > 1) {
            group_dist_matrix <- as.matrix(dist_matrix)[group_cells, group_cells]
            within_distances[[group]] <- group_dist_matrix[upper.tri(group_dist_matrix)]
        }
    }
    
    # Calculate between-group distances
    between_distances <- list()
    for (i in 1:(length(cell_groups)-1)) {
        for (j in (i+1):length(cell_groups)) {
            group1 <- cell_groups[i]
            group2 <- cell_groups[j]
            group1_cells <- which(analysis_data$cell_group == group1)
            group2_cells <- which(analysis_data$cell_group == group2)
            
            if (length(group1_cells) > 0 && length(group2_cells) > 0) {
                between_dist <- as.matrix(dist_matrix)[group1_cells, group2_cells]
                comparison_name <- paste(group1, "vs", group2, sep = "_")
                between_distances[[comparison_name]] <- as.vector(between_dist)
            }
        }
    }
    
    # 2. STATISTICAL TESTING
    cat("Performing statistical tests on UMAP1 distributions...\n")
    
    # Prepare data for statistical testing
    statistical_results <- data.frame()
    
    # Extract UMAP1 values for Post-IP_DE and Pre-IP_DE groups
    post_ip_de_umap1 <- analysis_data$UMAP1[analysis_data$cell_group == "Post-IP_DE"]
    pre_ip_de_umap1 <- analysis_data$UMAP1[analysis_data$cell_group == "Pre-IP_DE"]
    
    # Only proceed if both groups have cells
    if (length(post_ip_de_umap1) > 0 && length(pre_ip_de_umap1) > 0) {
        
        # Test for normality using Shapiro-Wilk test
        post_ip_de_normal <- ifelse(length(post_ip_de_umap1) >= 3 && length(post_ip_de_umap1) <= 5000,
                                    shapiro.test(post_ip_de_umap1)$p.value > 0.05, FALSE)
        pre_ip_de_normal <- ifelse(length(pre_ip_de_umap1) >= 3 && length(pre_ip_de_umap1) <= 5000,
                                   shapiro.test(pre_ip_de_umap1)$p.value > 0.05, FALSE)
        
        # Determine test type based on normality
        both_normal <- post_ip_de_normal && pre_ip_de_normal
        test_type <- ifelse(both_normal, "t-test", "Mann-Whitney")
        
        # Perform appropriate statistical test
        if (both_normal) {
            # Use t-test for normal distributions
            test_result <- t.test(post_ip_de_umap1, pre_ip_de_umap1)
            p_value <- test_result$p.value
            effect_size <- abs(mean(post_ip_de_umap1) - mean(pre_ip_de_umap1)) / 
                sqrt(((length(post_ip_de_umap1)-1)*var(post_ip_de_umap1) + 
                          (length(pre_ip_de_umap1)-1)*var(pre_ip_de_umap1)) / 
                         (length(post_ip_de_umap1) + length(pre_ip_de_umap1) - 2))
        } else {
            # Use Mann-Whitney U test for non-normal distributions
            test_result <- wilcox.test(post_ip_de_umap1, pre_ip_de_umap1)
            p_value <- test_result$p.value
            # Calculate Cliff's delta for effect size
            n1 <- length(post_ip_de_umap1)
            n2 <- length(pre_ip_de_umap1)
            u_stat <- test_result$statistic
            effect_size <- (2 * u_stat) / (n1 * n2) - 1
        }
        
        # Calculate descriptive statistics
        post_ip_de_stats <- list(
            median = median(post_ip_de_umap1),
            iqr = IQR(post_ip_de_umap1),
            mean = mean(post_ip_de_umap1),
            n = length(post_ip_de_umap1)
        )
        
        pre_ip_de_stats <- list(
            median = median(pre_ip_de_umap1),
            iqr = IQR(pre_ip_de_umap1),
            mean = mean(pre_ip_de_umap1),
            n = length(pre_ip_de_umap1)
        )
        
        # Store results
        result_row <- data.frame(
            Sample = sample_id,
            Cell_Group = "Post-IP_DE_vs_Pre-IP_DE",
            Test_Type = test_type,
            Post_IP_DE_Median = post_ip_de_stats$median,
            Post_IP_DE_IQR = post_ip_de_stats$iqr,
            Post_IP_DE_Mean = post_ip_de_stats$mean,
            Post_IP_DE_N = post_ip_de_stats$n,
            Post_IP_DE_Normal = post_ip_de_normal,
            Pre_IP_DE_Median = pre_ip_de_stats$median,
            Pre_IP_DE_IQR = pre_ip_de_stats$iqr,
            Pre_IP_DE_Mean = pre_ip_de_stats$mean,
            Pre_IP_DE_N = pre_ip_de_stats$n,
            Pre_IP_DE_Normal = pre_ip_de_normal,
            P_Value = p_value,
            Effect_Size = effect_size,
            stringsAsFactors = FALSE
        )
        
        statistical_results <- rbind(statistical_results, result_row)
        
        cat("Post-IP_DE: n =", post_ip_de_stats$n, ", Normal =", post_ip_de_normal, "\n")
        cat("Pre-IP_DE: n =", pre_ip_de_stats$n, ", Normal =", pre_ip_de_normal, "\n")
        cat("Test used:", test_type, ", p-value =", p_value, "\n")
    }
    
    # 3. SILHOUETTE ANALYSIS (Optional - for UMAP1 separation)
    cat("Calculating silhouette scores for UMAP1 separation...\n")
    
    # Create binary group labels: Post-IP_DE vs Pre-IP_DE only
    binary_group_labels <- ifelse(analysis_data$cell_group == "Post-IP_DE", 1, 
                                  ifelse(analysis_data$cell_group == "Pre-IP_DE", 2, NA))
    
    # Remove cells that are not Post-IP_DE or Pre-IP_DE
    valid_indices <- !is.na(binary_group_labels)
    binary_group_labels <- binary_group_labels[valid_indices]
    valid_coords <- umap_coords[valid_indices, ]
    
    # Calculate silhouette scores for binary classification (Post-IP_DE vs Pre-IP_DE only)
    if (length(unique(binary_group_labels)) == 2 && nrow(valid_coords) > 2) {
        valid_dist_matrix <- dist(valid_coords)
        silhouette_scores <- silhouette(binary_group_labels, valid_dist_matrix)
        avg_silhouette <- mean(silhouette_scores[,3])
        
        # Add silhouette results for overall separation
        silhouette_row <- data.frame(
            Sample = sample_id,
            Cell_Group = "Post-IP_DE_vs_Pre-IP_DE",
            Test_Type = "Silhouette_Score_UMAP_Separation",
            Post_IP_DE_Median = NA,
            Post_IP_DE_IQR = NA,
            Post_IP_DE_Mean = NA,
            Post_IP_DE_N = NA,
            Post_IP_DE_Normal = NA,
            Pre_IP_DE_Median = NA,
            Pre_IP_DE_IQR = NA,
            Pre_IP_DE_Mean = NA,
            Pre_IP_DE_N = NA,
            Pre_IP_DE_Normal = NA,
            P_Value = NA,
            Effect_Size = avg_silhouette,
            stringsAsFactors = FALSE
        )
        
        statistical_results <- rbind(statistical_results, silhouette_row)
    } else {
        avg_silhouette <- NA
    }
    
    # 4. APPLY BENJAMINI-HOCHBERG CORRECTION
    p_values <- statistical_results$P_Value[!is.na(statistical_results$P_Value)]
    if (length(p_values) > 0) {
        adjusted_p <- p.adjust(p_values, method = "BH")
        statistical_results$BH_Adjusted_P_Value <- NA
        statistical_results$BH_Adjusted_P_Value[!is.na(statistical_results$P_Value)] <- adjusted_p
    }
    
    # 5. CREATE RIDGELINE PLOT DATA - USING UMAP1 DISTRIBUTIONS
    cat("Preparing ridgeline plot data using UMAP1...\n")
    
    # Prepare data for ridgeline plots - use UMAP1 values for each cell group
    ridgeline_data <- data.frame(
        UMAP1 = analysis_data$UMAP1,
        Cell_Group = factor(analysis_data$cell_group, 
                            levels = c("Pre-IP_DE", "Pre-IP_Non-DE", "Post-IP_DE", "Post-IP_Non-DE")),
        stringsAsFactors = FALSE
    )
    
    # 6. CREATE RIDGELINE PLOT
    cat("Creating ridgeline plot...\n")
    
    ridgeline_plot <- ggplot(ridgeline_data, aes(x = UMAP1, y = Cell_Group, fill = Cell_Group)) +
        geom_density_ridges(alpha = 0.5, scale = 1.5, rel_min_height = 0.01) +
        scale_fill_manual(values = colors) +
        scale_y_discrete(limits = rev(c("Pre-IP_DE", "Pre-IP_Non-DE", "Post-IP_DE", "Post-IP_Non-DE")), 
                         expand = expansion(mult = c(0, 0.2))) +
        coord_cartesian(clip = "off") +
        labs(
            title = paste("Cell Type Clustering Analysis -", sample_id),
            subtitle = paste("UMAP1 distributions for each cell type"),
            x = "UMAP1 (Within Cell Type)",
            y = "Cell Type",
            caption = paste("Overall Silhouette Score:", round(avg_silhouette, 3), 
                            "| Lower spread = tighter clustering")
        ) +
        theme_ridges() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8, margin = margin(t = 2)),
            plot.caption = element_text(size = 9, hjust = 0),
            panel.grid.major.x = element_line(linetype = "dotted", color = "gray70"),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")
        ) +
        guides(fill = "none")  # Remove legend as colors match the trajectory plot
    
    print(ridgeline_plot)
    
    # 7. CREATE SUMMARY STATISTICS TABLE (FIXED)
    summary_stats <- statistical_results %>%
        filter(!grepl("Silhouette", Test_Type)) %>%
        select(Cell_Group, Post_IP_DE_Mean, Pre_IP_DE_Mean, 
               P_Value, BH_Adjusted_P_Value, Effect_Size, Test_Type) %>%
        mutate(
            Significant_BH = ifelse(BH_Adjusted_P_Value < 0.05, "Yes", "No"),
            Clustering_Quality = case_when(
                abs(Effect_Size) > 0.5 ~ "Strong",
                abs(Effect_Size) > 0.3 ~ "Moderate", 
                abs(Effect_Size) > 0.1 ~ "Weak",
                TRUE ~ "Very Weak"
            )
        )
    
    cat("\nClustering Summary:\n")
    print(summary_stats)
    
    # Return results
    return(list(
        statistical_results = statistical_results,
        ridgeline_data = ridgeline_data,
        ridgeline_plot = ridgeline_plot,
        summary_stats = summary_stats,
        overall_silhouette = avg_silhouette
    ))
}

# Function to process all samples and export results
process_all_samples_clustering <- function(all_trajectories, output_file = "monocle3_clustering_analysis_v6.xlsx") {
    cat("\n=== PROCESSING ALL SAMPLES FOR CLUSTERING ANALYSIS ===\n")
    
    # Define colors (same as in  main script)
    colors <- c(
        "Pre-IP_Non-DE" = "#CCCCCC",     # light gray
        "Post-IP_Non-DE" = "#666666",    # dark gray  
        "Pre-IP_DE" = "#DDA0DD",         # plum (light purple)
        "Post-IP_DE" = "#6666FF"         # blue
    )
    
    all_statistical_results <- data.frame()
    
    # Process each sample
    for (sample_name in names(all_trajectories)) {
        cat("\nProcessing", sample_name, "...\n")
        
        tryCatch({
            clustering_results <- analyze_cell_type_clustering(
                all_trajectories[[sample_name]], 
                sample_name, 
                colors
            )
            
            # Combine results
            all_statistical_results <- rbind(all_statistical_results, 
                                             clustering_results$statistical_results)
            
        }, error = function(e) {
            cat("Error analyzing clustering for", sample_name, ":", conditionMessage(e), "\n")
        })
    }
    
    # Apply Benjamini-Hochberg correction across ALL samples
    cat("Applying Benjamini-Hochberg correction across all samples...\n")
    
    if (nrow(all_statistical_results) > 0) {
        # Find rows with p-values (exclude silhouette rows)
        p_value_rows <- !is.na(all_statistical_results$P_Value)
        
        if (sum(p_value_rows) > 0) {
            # Extract p-values and apply BH correction
            p_values <- all_statistical_results$P_Value[p_value_rows]
            adjusted_p_values <- p.adjust(p_values, method = "BH")
            
            # Add BH_Adjusted_P_Value column and populate it
            all_statistical_results$BH_Adjusted_P_Value <- NA
            all_statistical_results$BH_Adjusted_P_Value[p_value_rows] <- adjusted_p_values
            
            cat("BH correction applied to", length(p_values), "p-values across all samples.\n")
        }
    }
    
    # Export results to Excel
    cat("\nExporting results to Excel file:", output_file, "\n")
    
    wb <- createWorkbook()
    
    # Add statistical results sheet
    addWorksheet(wb, "Statistical_Results")
    writeData(wb, "Statistical_Results", all_statistical_results)
    
    # Add metadata sheet
    metadata <- data.frame(
        Description = c(
            "Post_IP_DE_Median: Median UMAP1 value for Post-IP_DE cells",
            "Post_IP_DE_IQR: Interquartile range of UMAP1 values for Post-IP_DE cells",
            "Post_IP_DE_Mean: Mean UMAP1 value for Post-IP_DE cells", 
            "Post_IP_DE_N: Number of Post-IP_DE cells",
            "Post_IP_DE_Normal: Whether Post-IP_DE UMAP1 distribution is normal (Shapiro-Wilk p>0.05)",
            "Pre_IP_DE_Median: Median UMAP1 value for Pre-IP_DE cells",
            "Pre_IP_DE_IQR: Interquartile range of UMAP1 values for Pre-IP_DE cells",
            "Pre_IP_DE_Mean: Mean UMAP1 value for Pre-IP_DE cells",
            "Pre_IP_DE_N: Number of Pre-IP_DE cells", 
            "Pre_IP_DE_Normal: Whether Pre-IP_DE UMAP1 distribution is normal (Shapiro-Wilk p>0.05)",
            "P_Value: Raw p-value from t-test (if both normal) or Mann-Whitney test (if either non-normal)",
            "BH_Adjusted_P_Value: Benjamini-Hochberg corrected p-value (applied across all 7 samples)",
            "Effect_Size: Cohen's d (for t-test) or Cliff's delta (for Mann-Whitney test)",
            "Test_Type: Statistical test used (t-test or Mann-Whitney)"
        )
    )
    
    addWorksheet(wb, "Metadata")
    writeData(wb, "Metadata", metadata)
    
    # Save workbook
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    cat("Analysis complete! Results saved to:", output_file, "\n")
    
    # Print overall summary for Post-IP_DE vs Pre-IP_DE comparison
    cat("\n=== OVERALL CLUSTERING SUMMARY ===\n")
    
    if (nrow(all_statistical_results) > 0) {
        umap1_comparison_results <- all_statistical_results %>%
            filter(Cell_Group == "Post-IP_DE_vs_Pre-IP_DE" & !grepl("Silhouette", Test_Type))
        
        if (nrow(umap1_comparison_results) > 0) {
            overall_summary <- umap1_comparison_results %>%
                summarise(
                    Total_Samples = n(),
                    Mean_Effect_Size = mean(Effect_Size, na.rm = TRUE),
                    Mean_Post_IP_DE_UMAP1 = mean(Post_IP_DE_Mean, na.rm = TRUE),
                    Mean_Pre_IP_DE_UMAP1 = mean(Pre_IP_DE_Mean, na.rm = TRUE),
                    Significant_Samples_Raw = sum(P_Value < 0.05, na.rm = TRUE),
                    Significant_Samples_BH = sum(BH_Adjusted_P_Value < 0.05, na.rm = TRUE),
                    T_Test_Used = sum(Test_Type == "t-test", na.rm = TRUE),
                    Mann_Whitney_Used = sum(Test_Type == "Mann-Whitney", na.rm = TRUE),
                    .groups = 'drop'
                )
            
            print(overall_summary)
            
            cat("\nComparison of raw vs BH-adjusted significance:\n")
            sig_comparison <- umap1_comparison_results %>%
                select(Sample, P_Value, BH_Adjusted_P_Value) %>%
                mutate(
                    Raw_Significant = P_Value < 0.05,
                    BH_Significant = BH_Adjusted_P_Value < 0.05
                )
            print(sig_comparison)
            
        } else {
            cat("No UMAP1 comparison results found.\n")
        }
    } else {
        cat("No statistical results were generated.\n")
    }
    
    return(list(
        statistical_results = all_statistical_results
    ))
}

# Run the clustering analysis for all samples
cat("\n\nStarting clustering analysis for all samples...\n")
clustering_analysis_results <- process_all_samples_clustering(all_trajectories)


cat("\nClustering analysis completed!\n")
