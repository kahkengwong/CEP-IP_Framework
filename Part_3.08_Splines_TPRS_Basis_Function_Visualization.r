################################################
# Part 3.08: Building TPRS's Phi Visualization
################################################
# Define sample to be investigated
sample_data <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$gam_data

# Visualize accurate TPRS construction
visualize_tprs_construction <- function(sample_data, k = 10) {
    # Load required libraries
    library(ggplot2)
    library(mgcv)
    library(viridis)
    library(tidyr)
    library(gridExtra)
    
    # Set seed for reproducibility
    set.seed(123)
    
    # Create a sequence of x values for visualization
    x_seq <- seq(min(sample_data$TRPM4), max(sample_data$TRPM4), length.out = 500)
    
    # Fit a model to extract the basis functions
    model <- gam(Expression ~ s(TRPM4, bs = "tp", k = k), 
                 data = sample_data, method = "REML")
    
    # Get knot locations
    knots <- model$smooth[[1]]$xp
    if (is.null(knots)) {
        # Use mgcv's typical quantile-based knot placement
        knots <- quantile(sample_data$TRPM4, probs = seq(0, 1, length.out = k-2))
    }
    
    # Generate the X matrix with basis functions
    X_pred <- predict(model, newdata = data.frame(TRPM4 = x_seq), type = "lpmatrix")
    
    # Extract all basis functions
    smooth_cols <- grep("s\\(TRPM4", colnames(X_pred))
    phi_data <- data.frame(TRPM4 = x_seq)
    
    for(i in 1:length(smooth_cols)) {
        phi_data[paste0("phi", i)] <- X_pred[, smooth_cols[i]]
    }
    
    # Calculate data density for visualization
    density_data <- density(sample_data$TRPM4, n = 500)
    density_df <- data.frame(
        TRPM4 = density_data$x,
        density = density_data$y
    )
    # Scale density to fit nicely on the plot
    density_df$density_scaled <- scales::rescale(density_df$density, to = range(sample_data$Expression))
    
    # Plot 1: Data Distribution Analysis
    p1 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = Expression), 
                   alpha = 0.3, color = "gray50") +
        geom_line(data = density_df, aes(x = TRPM4, y = density_scaled), 
                  color = "#FFCC99", size = 1.2, alpha = 0.95) +
        labs(title = "1. Data Distribution Analysis",
             subtitle = "mgcv assesses data range and distribution patterns",
             x = "TRPM4 Expression (log2)", 
             y = "Expression\n(Red line: Data density)") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Plot 2: Knot Placement Strategy
    p2 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = Expression), 
                   alpha = 0.3, color = "gray50") +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        geom_point(data = data.frame(TRPM4 = knots, Expression = rep(mean(sample_data$Expression), length(knots))),
                   aes(x = TRPM4, y = Expression), color = "#F25B5B", size = 3, shape = 18, alpha = 0.75) +
        labs(title = "2. Knot Placement Strategy",
             subtitle = "mgcv places knots based on data distribution (quantiles/space-filling)",
             x = "TRPM4 Expression (log2)", 
             y = "Expression") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Plot 3: TPRS Radial Basis Construction
    # Create actual thin plate spline radial basis functions
    radial_basis_data <- data.frame()
    for (i in 1:min(8, length(knots))) {
        # TPRS radial basis function: |x - knot|^3 for univariate case
        radial_values <- abs(x_seq - knots[i])^3
        # Scale for visualization
        radial_values <- radial_values / max(radial_values) * 0.8
        
        temp_df <- data.frame(
            TRPM4 = x_seq,
            value = radial_values,
            basis_function = paste0("phi", i)
        )
        radial_basis_data <- rbind(radial_basis_data, temp_df)
    }
    
    # Define phi_colors using mako palette
    phi_colors <- viridis::mako(8, begin = 0.1, end = 0.9)
    names(phi_colors) <- paste0("phi", 1:8)
    
    p3 <- ggplot() +
        geom_line(data = radial_basis_data, 
                  aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "3. TPRS Radial Basis Construction",
             subtitle = "mgcv constructs |x-knot|³ radial basis functions",
             x = "TRPM4 Expression (log2)", 
             y = "Radial Basis Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Plot 4: Splines Formation (keep original accurate plot)
    transition_data_100_final <- data.frame()
    for (i in 1:min(8, length(knots))) {
        curve_x <- x_seq
        final_shape <- phi_data[, paste0("phi", i)]
        trans_df <- data.frame(
            TRPM4 = curve_x,
            value = final_shape,
            basis_function = paste0("phi", i)
        )
        transition_data_100_final <- rbind(transition_data_100_final, trans_df)
    }
    
    p4 <- ggplot() +
        geom_line(data = transition_data_100_final, 
                  aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "4. Splines Formation",
             subtitle = "Based on k, data distribution and data density at knots",
             x = "TRPM4 Expression (log2)", 
             y = "Basis Function Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "right",
              legend.title = element_blank())
    
    # Create a list of all plots for arrangeGrob
    all_plots <- list(p1, p2, p3, p4)
    
    # Filter out any NULL entries (just in case)
    all_plots <- all_plots[!sapply(all_plots, is.null)]
    
    # Arrange all plots
    final_plot <- do.call(arrangeGrob, c(all_plots, list(ncol = 1)))
    
    return(list(
        data_distribution = p1,
        knot_placement = p2,
        radial_basis = p3,
        splines_formation = p4,
        all_plots = final_plot
    ))
}

# Usage
results <- visualize_tprs_construction(sample_data, k = 10)
print(results$data_distribution)
print(results$knot_placement)
print(results$radial_basis)
print(results$splines_formation)

