##############################
# Building TPRS's Phi Visualization
##############################
# Define sample to be investigated
sample_data <- pca_results[["HYW_4881_Tumor"]][["Ribo"]]$gam_data

# Visualize phi emergence
visualize_phi1_emergence <- function(sample_data, k = 10) {
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
        # Fallback if knots are not directly available
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
    
    # Calculate residuals from linear fit
    linear_model <- lm(Expression ~ TRPM4, data = sample_data)
    sample_data$linear_pred <- predict(linear_model, newdata = sample_data)
    sample_data$residuals <- sample_data$Expression - sample_data$linear_pred
    linear_pred_seq <- predict(linear_model, newdata = data.frame(TRPM4 = x_seq))
    
    # Calculate a smoothed version of residuals for visualization
    smoother <- loess(residuals ~ TRPM4, data = sample_data, span = 0.5)
    residual_values <- predict(smoother, newdata = data.frame(TRPM4 = x_seq))
    flipped_residuals <- -residual_values
    
    # Scale to match φ1 range for comparison
    scale_to_range <- function(x, target_min, target_max) {
        x_min <- min(x, na.rm = TRUE)
        x_max <- max(x, na.rm = TRUE)
        scaled <- ((x - x_min) / (x_max - x_min)) * (target_max - target_min) + target_min
        return(scaled)
    }
    
    phi1_range <- range(phi_data$phi1)
    flipped_scaled <- scale_to_range(flipped_residuals, phi1_range[1], phi1_range[2])
    
    # New Plot 1: Linear Fit
    p1 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = Expression), alpha = 0.3, color = "gray40") +
        geom_line(data = data.frame(TRPM4 = x_seq, value = linear_pred_seq), 
                  aes(x = TRPM4, y = value), color = "#FFCC99", size = 1.2) +
        labs(title = "1. Linear Fit",
             subtitle = "Initial attempt to capture overall data shape",
             x = "TRPM4 Expression (log2)", 
             y = "Expression") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # New Plot 2: Actual Residual Plot
    p2 <- ggplot() +
        geom_point(data = sample_data, aes(x = TRPM4, y = residuals), alpha = 0.3, color = "gray40") +
        geom_line(data = data.frame(TRPM4 = x_seq, value = residual_values), 
                  aes(x = TRPM4, y = value), color = "#F25B5B", size = 1.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", alpha = 0.5) +
        labs(title = "2. Actual Residual Plot",
             subtitle = "Residuals from linear fit",
             x = "TRPM4 Expression (log2)", 
             y = "Residual Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Plot 3: Initial Flipped Residual Pattern
    p3 <- ggplot() +
        geom_line(data = data.frame(TRPM4 = x_seq, value = flipped_residuals), 
                  aes(x = TRPM4, y = value), color = "#F25B5B", size = 1.2) +
        geom_point(data = sample_data, aes(x = TRPM4, y = -residuals), 
                   alpha = 0.3, color = "gray40") +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        labs(title = "3. Initial Flipped Residual Pattern",
             subtitle = "Starting point for basis function development",
             x = "TRPM4 Expression (log2)", 
             y = "Flipped Residual Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
    
    # Define phi_colors using mako palette
    phi_colors <- viridis::mako(8, begin = 0.1, end = 0.9)
    names(phi_colors) <- paste0("phi", 1:8)
    
    # Plot 4: 100% Bell Shape
    transition_data_100_bell <- data.frame()
    for (i in 1:min(8, length(knots))) {
        width <- diff(range(sample_data$TRPM4)) / (1.5 * length(knots))
        curve_x <- x_seq
        bell_curve <- dnorm(curve_x, mean = knots[i], sd = width)
        # Scale bell curve to match the range of flipped residuals
        bell_curve <- scale_to_range(bell_curve, 0, max(abs(flipped_residuals)) * 0.8)
        trans_df <- data.frame(
            TRPM4 = curve_x,
            value = bell_curve,
            basis_function = paste0("phi", i)
        )
        transition_data_100_bell <- rbind(transition_data_100_bell, trans_df)
    }
    
    p4 <- ggplot() +
        geom_line(data = transition_data_100_bell, aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "4. 100% Bell Shape",
             subtitle = "Initial localized basis functions",
             x = "TRPM4 Expression (log2)", 
             y = "Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "none")
    
    # Plot 5: 100% Unpenalized Shape
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
    
    p5 <- ggplot() +
        geom_line(data = transition_data_100_final, aes(x = TRPM4, y = value, color = basis_function, group = basis_function), 
                  size = 1) +
        geom_vline(xintercept = knots, linetype = "dashed", color = "gray30", alpha = 0.5) +
        scale_color_manual(values = phi_colors) +
        labs(title = "5. Final Unpenalized Shapes",
             subtitle = "Computed basis functions used in the model",
             x = "TRPM4 Expression (log2)", 
             y = "Basis Function Value") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "right",
              legend.title = element_blank())
    
    # Create a list of all plots for arrangeGrob
    all_plots <- list(p1, p2, p3, p4, p5)
    
    # Filter out any NULL entries (just in case)
    all_plots <- all_plots[!sapply(all_plots, is.null)]
    
    # Arrange all plots
    final_plot <- do.call(arrangeGrob, c(all_plots, list(ncol = 1)))
    
    return(list(
        linear_fit = p1,
        actual_residuals = p2,
        flipped_residual = p3,
        transition_100_bell = p4,
        basis_final = p5,
        all_plots = final_plot
    ))
}

# Usage
results <- visualize_phi1_emergence(sample_data, k = 10)
print(results$linear_fit)
print(results$actual_residuals)
print(results$flipped_residual)
print(results$transition_100_bell)
print(results$basis_final)
