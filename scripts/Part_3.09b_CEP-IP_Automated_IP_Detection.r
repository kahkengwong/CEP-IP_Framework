##############################################################################
# Part 3.09b: CEP-IP IP Detection for TRPM4-Ribo (Prostate Cancer)
#            IP Reliability Score (IPRS) — Zero-crossing + Brent CI
#
# Dataset : Prostate cancer scRNA-seq (HYW cohort, 7 tumor samples)
# GOI     : TRPM4 (x-axis); Ribo composite (y-axis)
# Clusters: PCa clusters 6, 9, 11, 14, 19
#
# Method  : R(x) GAM zero-crossing with Brent 95% CI
#           + IP Reliability Score (IPRS) replacing categorical tiers
#
# Requires: pca_results (Part 3.04), prostate_results and prostate_ca_seurat (Part 1)
# Output  : ip_values                — named vector keyed by sample name
#           inflection_points        — same as ip_values (alias for Parts 3.10 to 3.15)
#           ip_results_prostate      — full IP result list (all IPRS fields)
#           ip_table_prostate        — data.frame summary
#           Excel: TRPM4_Ribo_IP_IPRS_Results_<timestamp>.xlsx (3 sheets)
#
# Algorithm (per sample):
#   1. Classify top (dev_expl × n) cells as "purple" (highest EP)
#      Null model uses same cell weights as fitted GAM (consistent EP)
#   2. Compute signed residuals for purple cells
#   3. Fit residual GAM R(x) via thin-plate REML spline on purple cells
#   4. Evaluate R(x) on 2000-pt grid; find first neg→pos zero-crossing
#      (linear interpolation for sub-grid precision)
#   5. Brent 95% CI: lower = zero of [R_hat + 1.96·SE]
#                    upper = zero of [R_hat − 1.96·SE]
#      Proximity fallback if one bound is NA (data-sparse edge)
#   6. Stat test (SW-gated Welch-t / Mann-Whitney U) pre vs post IP
#   7. IPRS — continuous 0–1 composite:
#        C1 (Precision)   = max(0, 1 − CI_Width / TRPM4_range)
#        C2 (Statistical) = tanh(atanh(0.5)/(-log10(0.05)) × (−log10(p_adj)))
#        C3 (Abruptness)  = tanh(|R'(IP)|)
#        C4 (GAM Ribbon)  = exp(−mean_ribbon_width / amplitude)
#        C5 (Smoothness)  = exp(−(TV/amplitude−1)) × (1−neg_post_frac)²
#        IPRS = (C1 + C2 + C3 + C4 + C5) / 5
#      Tiers: ≥0.80 Very Strong | 0.60–0.79 Strong |
#             0.40–0.59 Moderate (flag) | 0.20–0.39 Weak (flag) |
#             <0.20 Failed (no IP — do not proceed)
##############################################################################

library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)
library(Seurat)
library(openxlsx)

# ------------------------------------------------------------------
# SECTION 0: Pre-flight checks
# ------------------------------------------------------------------
if (!exists("pca_results"))        stop("pca_results not found. Run prior parts first.")
if (!exists("prostate_results"))   stop("prostate_results not found. Run prior parts first.")
if (!exists("prostate_ca_seurat")) stop("prostate_ca_seurat not found. Run prior parts first.")

tumor_samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor",
                   "HYW_4881_Tumor", "HYW_5386_Tumor", "HYW_5742_Tumor",
                   "HYW_5755_Tumor")

PCA_CLUSTERS <- c(6, 9, 11, 14, 19)

# ------------------------------------------------------------------
# SECTION 1: Tuning parameters
# ------------------------------------------------------------------
STAT_ALPHA        <- 0.05
SHAPIRO_ALPHA     <- 0.05
POST_IP_MIN_CELLS <- 5L
MIN_PURPLE_CELLS  <- 10L

IPRS_VERY_STRONG  <- 0.80
IPRS_STRONG       <- 0.60
IPRS_MODERATE     <- 0.40
IPRS_WEAK         <- 0.20

# ------------------------------------------------------------------
# SECTION 2: Helper — SW-gated two-group stat test
# ------------------------------------------------------------------
stat_test_zc <- function(pre_r, post_r) {
    out <- list(stat_test   = NA_character_,
                stat_p_raw  = NA_real_,
                stat_p_adj  = NA_real_,
                both_normal = NA)

    if (length(pre_r) < 3L || length(post_r) < POST_IP_MIN_CELLS)
        return(out)

    sw_pre  <- if (length(pre_r)  >= 3L)
        suppressWarnings(shapiro.test(pre_r))$p.value  else NA_real_
    sw_post <- if (length(post_r) >= 3L)
        suppressWarnings(shapiro.test(post_r))$p.value else NA_real_
    both_normal <- !is.na(sw_pre) && !is.na(sw_post) &&
        sw_pre > SHAPIRO_ALPHA && sw_post > SHAPIRO_ALPHA
    out$both_normal <- both_normal

    if (isTRUE(both_normal)) {
        tt <- suppressWarnings(
            t.test(post_r, pre_r, alternative = "greater", var.equal = FALSE))
        out$stat_test  <- "Welch_t"
        out$stat_p_raw <- tt$p.value
    } else {
        mw <- suppressWarnings(
            wilcox.test(post_r, pre_r, alternative = "greater", exact = FALSE))
        out$stat_test  <- "MannWhitneyU"
        out$stat_p_raw <- mw$p.value
    }
    out$stat_p_adj <- p.adjust(out$stat_p_raw, method = "BH")
    out
}

# ------------------------------------------------------------------
# SECTION 3: Helper — Brent CI on R(x) GAM confidence bands
# ------------------------------------------------------------------
find_ci_brent <- function(residual_gam, x_seq, r_hat, r_se,
                          multiplier, near_x) {
    band     <- r_hat + multiplier * r_se
    sign_chg <- which(diff(sign(band)) != 0L)
    if (length(sign_chg) == 0L) return(NA_real_)

    if (length(sign_chg) > 1L)
        sign_chg <- sign_chg[which.min(abs(x_seq[sign_chg] - near_x))]

    bracket <- c(x_seq[sign_chg], x_seq[sign_chg + 1L])

    tryCatch({
        uniroot(
            f = function(x) {
                p <- predict(residual_gam,
                             newdata = data.frame(TRPM4 = x), se.fit = TRUE)
                as.numeric(p$fit) + multiplier * as.numeric(p$se.fit)
            },
            interval = bracket,
            tol      = .Machine$double.eps^0.5
        )$root
    }, error = function(e) {
        x1 <- bracket[1L]; x2 <- bracket[2L]
        y1 <- band[sign_chg]; y2 <- band[sign_chg + 1L]
        x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
    })
}

# ------------------------------------------------------------------
# SECTION 4: Helper — Proximity fallback CI bound
# ------------------------------------------------------------------
find_ci_proximity_fallback <- function(x_seq, r_hat, r_se, multiplier,
                                       ip_value, side = c("lower", "upper")) {
    side <- match.arg(side)
    band <- r_hat + multiplier * r_se

    idx <- if (side == "lower") which(x_seq <= ip_value) else which(x_seq >= ip_value)
    if (length(idx) < 3L) return(NA_real_)

    band_region <- band[idx]
    if (any(diff(sign(band_region)) != 0L)) return(NA_real_)

    x_seq[idx[which.min(abs(band_region))]]
}

# ------------------------------------------------------------------
# SECTION 5: Helper — IPRS component calculation
# ------------------------------------------------------------------
compute_iprs <- function(ip_value, ip_ci_lower, ip_ci_upper,
                         stat_p_adj, x_seq, r_hat, r_upper, r_lower,
                         residual_data) {

    trpm4_range <- max(x_seq) - min(x_seq)

    both_ci  <- !is.na(ip_ci_lower) && !is.na(ip_ci_upper)
    lower_ok <- !is.na(ip_ci_lower)
    upper_ok <- !is.na(ip_ci_upper)
    ci_width <- if (both_ci) abs(ip_ci_upper - ip_ci_lower) else NA_real_

    # C1: Precision
    C1 <- if (both_ci && trpm4_range > 0) {
        max(0, 1 - ci_width / trpm4_range)
    } else if ((lower_ok || upper_ok) && !is.na(ip_value) && trpm4_range > 0) {
        half_width <- if (lower_ok) abs(ip_value - ip_ci_lower) else abs(ip_ci_upper - ip_value)
        0.5 * max(0, 1 - half_width / trpm4_range)
    } else {
        0
    }

    C1_note <- if (both_ci)
        sprintf("Both CI bounds: CI_width=%.4f / TRPM4_range=%.4f → C1=%.4f",
                ci_width, trpm4_range, C1)
    else if (lower_ok || upper_ok) {
        known_bound <- if (lower_ok) ip_ci_lower else ip_ci_upper
        bound_label <- if (lower_ok) "lower" else "upper"
        half_w      <- abs(ip_value - known_bound)
        sprintf("Partial CI (%s only): half_width=%.4f / TRPM4_range=%.4f → 0.5×max(0,1-ratio)=%.4f",
                bound_label, half_w, trpm4_range, C1)
    } else
        "Both CI bounds NA → C1=0"

    # C2: Statistical
    GW_ALPHA  <- 0.05
    scale_c2  <- atanh(0.5) / (-log10(GW_ALPHA))
    C2 <- if (!is.na(stat_p_adj))
        tanh(scale_c2 * (-log10(max(stat_p_adj, 1e-300))))
    else
        0

    C2_note <- if (is.na(stat_p_adj))
        "stat_p_adj not available → C2=0"
    else
        sprintf("tanh(%.4f * -log10(%.4e)) → C2=%.4f", scale_c2, stat_p_adj, C2)

    # C3: Abruptness
    deriv_at_ip <- if (!is.na(ip_value) && length(x_seq) >= 3L) {
        i_ip <- which.min(abs(x_seq - ip_value))
        n    <- length(x_seq)
        if (i_ip > 1L && i_ip < n) {
            abs((r_hat[i_ip + 1L] - r_hat[i_ip - 1L]) /
                    (x_seq[i_ip + 1L] - x_seq[i_ip - 1L]))
        } else if (i_ip == 1L) {
            abs((r_hat[2L] - r_hat[1L]) / (x_seq[2L] - x_seq[1L]))
        } else {
            abs((r_hat[n] - r_hat[n - 1L]) / (x_seq[n] - x_seq[n - 1L]))
        }
    } else NA_real_

    C3      <- if (!is.na(deriv_at_ip)) tanh(deriv_at_ip) else 0
    C3_note <- if (is.na(deriv_at_ip))
        "IP not available → C3=0"
    else
        sprintf("|R'(IP)|=%.4f → tanh(%.4f) → C3=%.4f", deriv_at_ip, deriv_at_ip, C3)

    # C4: GAM ribbon width penalty
    amplitude   <- diff(range(r_hat, na.rm = TRUE))
    mean_ribbon <- mean(r_upper - r_lower, na.rm = TRUE)
    C4          <- exp(-mean_ribbon / (amplitude + 1e-6))

    C4_note <- sprintf(
        "mean_ribbon=%.4f / amplitude=%.4f → exp(-ratio=%.4f) → C4=%.4f",
        mean_ribbon, amplitude, mean_ribbon / (amplitude + 1e-6), C4)

    # C5: R(x) smoothness + post-IP sign fidelity
    tv <- sum(abs(diff(r_hat)), na.rm = TRUE)

    if (!is.na(ip_value) && nrow(residual_data) > 0L) {
        post_cells    <- residual_data$residual[residual_data$TRPM4 >= ip_value]
        n_post        <- length(post_cells)
        neg_post_frac <- if (n_post > 0L) sum(post_cells < 0) / n_post else 0
    } else {
        n_post        <- 0L
        neg_post_frac <- 0
    }
    sign_fidelity <- (1 - neg_post_frac)^2

    if (amplitude < 1e-6) {
        C5      <- 0
        C5_note <- sprintf("amplitude=%.6f < 1e-6 (flat R(x)) → C5=0", amplitude)
    } else {
        excess  <- tv / (amplitude + 1e-6) - 1
        wave_sc <- exp(-max(0, excess))
        C5      <- wave_sc * sign_fidelity
        C5_note <- sprintf(
            "TV=%.4f/amp=%.4f → wave=%.4f | neg_post=%.3f → sign=(1-%.3f)^2=%.4f | C5=%.4f*%.4f=%.4f",
            tv, amplitude, wave_sc, neg_post_frac, neg_post_frac,
            sign_fidelity, wave_sc, sign_fidelity, C5)
    }

    IPRS      <- (C1 + C2 + C3 + C4 + C5) / 5
    any_zero  <- any(c(C1, C2, C3, C4, C5) == 0)
    IPRS_tier <- if (IPRS >= IPRS_VERY_STRONG) "Very Strong" else
        if (IPRS >= IPRS_STRONG)      "Strong"      else
            if (IPRS >= IPRS_MODERATE)    "Moderate"    else
                if (IPRS >= IPRS_WEAK)        "Weak"        else
                    "Failed"

    IPRS_flag <- if (any_zero) "Visual review needed [C=0]" else
        if (IPRS_tier %in% c("Very Strong", "Strong")) "None" else
            if (IPRS_tier %in% c("Moderate", "Weak"))  "Visual review needed" else
                "No IP - do not proceed"

    cat(sprintf("  IPRS components: C1=%.4f [%s]\n", C1, C1_note))
    cat(sprintf("                   C2=%.4f [%s]\n", C2, C2_note))
    cat(sprintf("                   C3=%.4f [%s]\n", C3, C3_note))
    cat(sprintf("                   C4=%.4f [%s]\n", C4, C4_note))
    cat(sprintf("                   C5=%.4f [%s]\n", C5, C5_note))
    cat(sprintf("  IPRS = (%.4f+%.4f+%.4f+%.4f+%.4f)/5 = %.4f  [%s]  %s\n",
                C1, C2, C3, C4, C5, IPRS, IPRS_tier,
                if (IPRS_flag == "None") "" else paste0("*** ", IPRS_flag, " ***")))

    list(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5,
         IPRS = IPRS, IPRS_tier = IPRS_tier, IPRS_flag = IPRS_flag,
         C1_note = C1_note, C2_note = C2_note, C3_note = C3_note,
         C4_note = C4_note, C5_note = C5_note,
         deriv_at_ip = deriv_at_ip, ci_width = ci_width,
         trpm4_range = trpm4_range, amplitude = amplitude,
         mean_ribbon = mean_ribbon, tv = tv)
}

# ------------------------------------------------------------------
# SECTION 6: Core IP detection function (per sample)
# ------------------------------------------------------------------
detect_ip_prostate <- function(current_sample) {
    cat("\n--- Processing:", current_sample, "---\n")

    # ---- 6a. Retrieve pre-fitted GAM model + data ----
    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    sample_data      <- pca_results[[current_sample]][["Ribo"]]$gam_data
    dev_explained    <- summary(sample_gam_model)$dev.expl

    cat("  Cells in gam_data:", nrow(sample_data),
        "| Dev explained:", sprintf("%.1f%%", dev_explained * 100), "\n")

    # ---- 6b. Restrict to cells present in Seurat object ----
    cluster_cells <- WhichCells(prostate_results$seurat_obj, idents = PCA_CLUSTERS)
    sample_cells  <- WhichCells(prostate_results$seurat_obj,
                                cells = grep(current_sample,
                                             colnames(prostate_results$seurat_obj),
                                             value = TRUE))
    selected_cells <- intersect(cluster_cells, sample_cells)
    common_cells   <- intersect(rownames(sample_data), selected_cells)

    if (length(common_cells) == 0) {
        cat("  *** NO COMMON CELLS — skipping ***\n")
        return(NULL)
    }

    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]
    cat("  Cells after Seurat intersection:", nrow(sample_data), "\n")

    # ---- 6c. FAILED_list helper ----
    FAILED_list <- function(status) {
        iprs_fail <- list(C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0,
                          IPRS = 0, IPRS_tier = "Failed",
                          IPRS_flag = "No IP — do not proceed",
                          C1_note = "IP detection failed", C2_note = "IP detection failed",
                          C3_note = "IP detection failed", C4_note = "IP detection failed",
                          C5_note = "IP detection failed",
                          deriv_at_ip = NA_real_, ci_width = NA_real_,
                          trpm4_range = NA_real_, amplitude = NA_real_,
                          mean_ribbon = NA_real_, tv = NA_real_)
        list(ip_value = NA_real_, manual_flag = TRUE, status = status,
             stat_test = NA_character_, stat_p_raw = NA_real_,
             stat_p_adj = NA_real_, both_normal = NA,
             ip_ci_lower = NA_real_, ip_ci_upper = NA_real_,
             ci_lower_method = NA_character_, ci_upper_method = NA_character_,
             neg_pos_count = 0L, iprs = iprs_fail,
             residual_data = data.frame(TRPM4 = numeric(0), residual = numeric(0)),
             residual_pred = data.frame(TRPM4 = numeric(0), r_hat = numeric(0),
                                        se = numeric(0), r_upper = numeric(0),
                                        r_lower = numeric(0)))
    }

    # ---- 6d. Purple-cell classification via EP (weighted null model) ----
    # Check if cell_weights are stored; fall back to uniform weights if not
    cell_weights <- if (!is.null(pca_results[[current_sample]][["Ribo"]]$cell_weights)) {
        wts <- pca_results[[current_sample]][["Ribo"]]$cell_weights
        wts[rownames(sample_data)]
    } else {
        rep(1, nrow(sample_data))
    }

    null_model   <- gam(Expression ~ 1, data = sample_data, weights = cell_weights)
    null_fitted  <- fitted(null_model)
    model_fitted <- fitted(sample_gam_model)

    expl_power <- 1 - ((sample_data$Expression - model_fitted)^2 /
                           (sample_data$Expression - null_fitted)^2)

    n_purple  <- round(nrow(sample_data) * dev_explained)
    is_purple <- rep(FALSE, nrow(sample_data))
    is_purple[order(expl_power, decreasing = TRUE)[seq_len(n_purple)]] <- TRUE

    # ---- 6e. Signed residuals for purple cells ----
    residual    <- sample_data$Expression -
        predict(sample_gam_model, newdata = sample_data)
    purple_idx  <- which(is_purple)
    residual_data <- data.frame(
        TRPM4    = sample_data$TRPM4[purple_idx],
        residual = residual[purple_idx]
    )

    cat("  Purple cells:", nrow(residual_data),
        "| Above:", sum(residual_data$residual > 0),
        "| Below:", sum(residual_data$residual < 0), "\n")

    if (nrow(residual_data) < MIN_PURPLE_CELLS ||
        length(unique(residual_data$TRPM4)) < 3L) {
        cat("  *** INSUFFICIENT PURPLE CELLS — skipping IP detection ***\n")
        return(FAILED_list("insufficient_purple_cells"))
    }

    # ---- 6f. Residual GAM R(x) — thin-plate REML, unweighted ----
    set.seed(123)
    n_unique     <- length(unique(residual_data$TRPM4))
    k_resid      <- min(20L, n_unique - 1L)
    residual_gam <- tryCatch(
        gam(residual ~ s(TRPM4, bs = "tp", k = k_resid),
            data = residual_data, method = "REML"),
        error = function(e) {
            cat("  Fallback k=8:", conditionMessage(e), "\n")
            gam(residual ~ s(TRPM4, bs = "tp",
                             k = min(8L, n_unique - 1L)),
                data = residual_data, method = "REML")
        }
    )
    gam_sum <- summary(residual_gam)
    cat("  R(x) GAM: edf =", sprintf("%.2f", gam_sum$edf),
        "| Dev =", sprintf("%.1f%%", gam_sum$dev.expl * 100),
        "| p =", sprintf("%.2e", gam_sum$s.table[4]), "\n")

    # ---- 6g. 2000-pt grid prediction ----
    x_seq <- seq(min(residual_data$TRPM4), max(residual_data$TRPM4),
                 length.out = 2000L)
    pred  <- predict(residual_gam,
                     newdata = data.frame(TRPM4 = x_seq), se.fit = TRUE)
    r_hat   <- as.numeric(pred$fit)
    r_se    <- as.numeric(pred$se.fit)
    r_upper <- r_hat + 1.96 * r_se
    r_lower <- r_hat - 1.96 * r_se

    # ---- 6h. First neg→pos zero-crossing (linear interpolation) ----
    sc_all <- which(diff(sign(r_hat)) != 0L)
    sc_np  <- sc_all[r_hat[sc_all] < 0 & r_hat[sc_all + 1L] > 0]

    ip_value    <- NA_real_
    manual_flag <- FALSE

    if (length(sc_np) > 0L) {
        i_cross  <- sc_np[1L]
        x0 <- x_seq[i_cross];      y0 <- r_hat[i_cross]
        x1 <- x_seq[i_cross + 1L]; y1 <- r_hat[i_cross + 1L]
        ip_value <- x0 - y0 * (x1 - x0) / (y1 - y0)
        cat("  IP (R(x) zero-crossing):", sprintf("%.6f", ip_value),
            "| neg\u2192pos crossings:", length(sc_np), "\n")
    } else {
        manual_flag <- TRUE
        cat("  *** R(x) has no neg\u2192pos zero-crossing",
            "— R(x) is entirely",
            ifelse(all(r_hat >= 0), "positive", "negative"), "***\n")
    }

    # ---- 6i. Brent 95% CI with proximity fallback ----
    ci_anchor <- if (!is.na(ip_value)) ip_value else median(x_seq)

    ip_ci_lower     <- tryCatch(
        find_ci_brent(residual_gam, x_seq, r_hat, r_se, +1.96, ci_anchor),
        error = function(e) NA_real_)
    ci_lower_method <- if (!is.na(ip_ci_lower)) "brent" else NA_character_

    ip_ci_upper     <- tryCatch(
        find_ci_brent(residual_gam, x_seq, r_hat, r_se, -1.96, ci_anchor),
        error = function(e) NA_real_)
    ci_upper_method <- if (!is.na(ip_ci_upper)) "brent" else NA_character_

    if (is.na(ip_ci_lower) && !is.na(ip_value)) {
        fb_lower <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                               +1.96, ip_value, side = "lower")
        if (!is.na(fb_lower)) {
            ip_ci_lower     <- fb_lower
            ci_lower_method <- "proximity_fallback"
            cat("  CI lower: Brent NA → proximity fallback =",
                sprintf("%.4f", ip_ci_lower), "\n")
        }
    }
    if (is.na(ip_ci_upper) && !is.na(ip_value)) {
        fb_upper <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                               -1.96, ip_value, side = "upper")
        if (!is.na(fb_upper)) {
            ip_ci_upper     <- fb_upper
            ci_upper_method <- "proximity_fallback"
            cat("  CI upper: Brent NA → proximity fallback =",
                sprintf("%.4f", ip_ci_upper), "\n")
        }
    }

    # Ensure lower ≤ upper
    if (!is.na(ip_ci_lower) && !is.na(ip_ci_upper) && ip_ci_lower > ip_ci_upper) {
        tmp <- ip_ci_lower; ip_ci_lower <- ip_ci_upper; ip_ci_upper <- tmp
        tmp_m <- ci_lower_method; ci_lower_method <- ci_upper_method; ci_upper_method <- tmp_m
    }

    cat("  95% CI: [",
        if (is.na(ip_ci_lower)) "NA" else sprintf("%.4f [%s]", ip_ci_lower, ci_lower_method),
        ",",
        if (is.na(ip_ci_upper)) "NA" else sprintf("%.4f [%s]", ip_ci_upper, ci_upper_method),
        "] | width:",
        if (is.na(ip_ci_lower) || is.na(ip_ci_upper)) "NA"
        else sprintf("%.4f", abs(ip_ci_upper - ip_ci_lower)), "\n")

    # ---- 6j. Stat test pre vs post IP ----
    if (!is.na(ip_value)) {
        pre_r  <- residual_data$residual[residual_data$TRPM4 <  ip_value]
        post_r <- residual_data$residual[residual_data$TRPM4 >= ip_value]
        st     <- stat_test_zc(pre_r, post_r)
        cat("  Stat:", st$stat_test,
            sprintf("| p_adj: %s",
                    if (is.na(st$stat_p_adj)) "NA"
                    else sprintf("%.3e", st$stat_p_adj)), "\n")
    } else {
        st <- list(stat_test = NA_character_, stat_p_raw = NA_real_,
                   stat_p_adj = NA_real_, both_normal = NA)
    }

    # ---- 6k. Compute IPRS ----
    iprs <- compute_iprs(ip_value       = ip_value,
                         ip_ci_lower    = ip_ci_lower,
                         ip_ci_upper    = ip_ci_upper,
                         stat_p_adj     = st$stat_p_adj,
                         x_seq          = x_seq,
                         r_hat          = r_hat,
                         r_upper        = r_upper,
                         r_lower        = r_lower,
                         residual_data  = residual_data)

    list(ip_value        = ip_value,
         manual_flag     = manual_flag,
         status          = "ok",
         stat_test       = st$stat_test,
         stat_p_raw      = st$stat_p_raw,
         stat_p_adj      = st$stat_p_adj,
         both_normal     = st$both_normal,
         ip_ci_lower     = ip_ci_lower,
         ip_ci_upper     = ip_ci_upper,
         ci_lower_method = ci_lower_method,
         ci_upper_method = ci_upper_method,
         neg_pos_count   = length(sc_np),
         iprs            = iprs,
         residual_data   = residual_data,
         residual_pred   = data.frame(
             TRPM4   = x_seq,
             r_hat   = r_hat,
             se      = r_se,
             r_upper = r_upper,
             r_lower = r_lower
         ))
}

# ------------------------------------------------------------------
# SECTION 7: Run IP detection across all samples
# ------------------------------------------------------------------
set.seed(123)
cat("\n=============================================\n")
cat(sprintf("IP DETECTION | Zero-crossing + Brent CI + IPRS | STAT_ALPHA=%.2f\n",
            STAT_ALPHA))
cat("Dataset  : Prostate Cancer scRNA-seq (HYW cohort)\n")
cat("GOI      : TRPM4 (x-axis) vs Ribo composite (y-axis)\n")
cat("Clusters : PCa clusters", paste(PCA_CLUSTERS, collapse = ", "), "\n")
cat("Samples  :", paste(tumor_samples, collapse = ", "), "\n")
cat("=============================================\n")

ip_results_prostate <- list()

for (sample_name in tumor_samples) {
    tryCatch({
        ip_results_prostate[[sample_name]] <- detect_ip_prostate(sample_name)
    }, error = function(e) {
        cat("  ERROR in", sample_name, ":", conditionMessage(e), "\n")
        ip_results_prostate[[sample_name]] <<- NULL
    })
}

# ------------------------------------------------------------------
# SECTION 8: Console summary
# ------------------------------------------------------------------
cat("\n\n=============================================\n")
cat("IP SUMMARY — IP Reliability Score (IPRS)\n")
cat("=============================================\n\n")
cat(sprintf("%-26s %8s %8s %8s %6s %6s %6s %6s %6s %6s  %-14s  %s\n",
            "Sample", "IP", "CI_Lo", "CI_Hi",
            "C1", "C2", "C3", "C4", "C5", "IPRS", "Tier", "Flag"))
cat(strrep("-", 130), "\n")

for (nm in tumor_samples) {
    r <- ip_results_prostate[[nm]]
    if (is.null(r)) {
        cat(sprintf("%-26s  %-80s\n", nm, "NULL result — check errors above"))
        next
    }
    iprs <- r$iprs
    cat(sprintf("%-26s %8s %8s %8s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  %-14s  %s\n",
                nm,
                if (is.na(r$ip_value))    "NA" else sprintf("%.4f", r$ip_value),
                if (is.na(r$ip_ci_lower)) "NA" else sprintf("%.4f", r$ip_ci_lower),
                if (is.na(r$ip_ci_upper)) "NA" else sprintf("%.4f", r$ip_ci_upper),
                iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5, iprs$IPRS,
                iprs$IPRS_tier, iprs$IPRS_flag))
}

cat("\n")
for (nm in tumor_samples) {
    r <- ip_results_prostate[[nm]]
    if (!is.null(r) && r$iprs$IPRS_flag != "None")
        cat(sprintf("  [FLAG] %s | IPRS=%.3f (%s) → %s\n",
                    nm, r$iprs$IPRS, r$iprs$IPRS_tier, r$iprs$IPRS_flag))
}

# ------------------------------------------------------------------
# SECTION 9: Export environment vectors
# (named to match what Parts 3.10 to 3.15 expect: ip_values + inflection_points)
# ------------------------------------------------------------------

# Build ip_table_prostate data frame
ip_table_prostate <- do.call(rbind, lapply(tumor_samples, function(nm) {
    r    <- ip_results_prostate[[nm]]
    iprs <- if (!is.null(r)) r$iprs else list(
        C1 = NA_real_, C2 = NA_real_, C3 = NA_real_, C4 = NA_real_, C5 = NA_real_,
        IPRS = NA_real_, IPRS_tier = NA_character_, IPRS_flag = NA_character_,
        ci_width = NA_real_, trpm4_range = NA_real_)
    data.frame(
        Sample          = nm,
        Patient_ID      = gsub("_Tumor$", "", nm),
        IP_Value        = if (!is.null(r)) r$ip_value      else NA_real_,
        IP_CI_Lower     = if (!is.null(r)) r$ip_ci_lower   else NA_real_,
        IP_CI_Upper     = if (!is.null(r)) r$ip_ci_upper   else NA_real_,
        CI_Width        = iprs$ci_width,
        TRPM4_Range     = iprs$trpm4_range,
        CI_Lower_Method = if (!is.null(r)) r$ci_lower_method else NA_character_,
        CI_Upper_Method = if (!is.null(r)) r$ci_upper_method else NA_character_,
        NegPos_Crossings = if (!is.null(r)) r$neg_pos_count  else NA_integer_,
        N_Purple_Cells   = if (!is.null(r)) nrow(r$residual_data) else NA_integer_,
        Stat_Test        = if (!is.null(r)) r$stat_test     else NA_character_,
        Stat_P_Raw       = if (!is.null(r)) r$stat_p_raw    else NA_real_,
        Stat_P_Adj       = if (!is.null(r)) r$stat_p_adj    else NA_real_,
        C1_Precision     = iprs$C1,
        C2_Statistical   = iprs$C2,
        C3_Abruptness    = iprs$C3,
        C4_GAM_Ribbon    = iprs$C4,
        C5_Smoothness    = iprs$C5,
        IPRS             = iprs$IPRS,
        IPRS_Tier        = iprs$IPRS_tier,
        IPRS_Flag        = iprs$IPRS_flag,
        stringsAsFactors = FALSE
    )
}))
rownames(ip_table_prostate) <- NULL

# ── Primary export vector: ip_values (keyed by sample name) ──────────────────
# This is what Parts 3.12, 3.13 and 3.15 check for with: if (!exists("ip_values")) stop(...)
ip_values <- setNames(ip_table_prostate$IP_Value, ip_table_prostate$Sample)

# ── Alias: inflection_points (same values, same names) ───────────────────────
# Part 3.13 also uses: inflection_points <- ip_values[patient_ids]
# We export inflection_points here directly so it is also available immediately
inflection_points <- ip_values

cat("\nEnvironment vectors exported:\n")
cat("  ip_values          :", length(ip_values), "entries (keyed by sample name)\n")
cat("  inflection_points  :", length(inflection_points), "entries (alias of ip_values)\n")
cat("  ip_results_prostate:", length(ip_results_prostate), "samples\n")
cat("  ip_table_prostate  : data.frame,", nrow(ip_table_prostate), "rows\n")

cat("\nip_values:\n")
print(ip_values)

# ------------------------------------------------------------------
# SECTION 10: Excel export (3 sheets)
# ------------------------------------------------------------------
cat("\nExporting Excel results...\n")

# Sheet 1: IP results
ip_sheet <- do.call(rbind, lapply(tumor_samples, function(nm) {
    r    <- ip_results_prostate[[nm]]
    iprs <- if (!is.null(r)) r$iprs else list(
        C1 = NA, C2 = NA, C3 = NA, C4 = NA, C5 = NA,
        IPRS = NA, IPRS_tier = NA, IPRS_flag = NA,
        C1_note = NA, C2_note = NA, C3_note = NA, C4_note = NA, C5_note = NA,
        ci_width = NA, trpm4_range = NA, amplitude = NA,
        mean_ribbon = NA, tv = NA, deriv_at_ip = NA)
    data.frame(
        Sample           = nm,
        Patient_ID       = gsub("_Tumor$", "", nm),
        Gene_Set         = "Ribo",
        Gene_Set_Label   = "Ribosomal composite",
        IP_Value         = if (!is.null(r)) r$ip_value      else NA_real_,
        CI_Lower_95      = if (!is.null(r)) r$ip_ci_lower   else NA_real_,
        CI_Upper_95      = if (!is.null(r)) r$ip_ci_upper   else NA_real_,
        CI_Width         = iprs$ci_width,
        TRPM4_Range      = iprs$trpm4_range,
        CI_Lower_Method  = if (!is.null(r)) r$ci_lower_method else NA_character_,
        CI_Upper_Method  = if (!is.null(r)) r$ci_upper_method else NA_character_,
        NegPos_Crossings = if (!is.null(r)) r$neg_pos_count  else NA_integer_,
        N_Purple_Cells   = if (!is.null(r)) nrow(r$residual_data) else NA_integer_,
        Stat_Test        = if (!is.null(r)) r$stat_test     else NA_character_,
        Stat_P_Raw       = if (!is.null(r)) r$stat_p_raw    else NA_real_,
        Stat_P_Adj       = if (!is.null(r)) r$stat_p_adj    else NA_real_,
        C1_Precision     = iprs$C1,
        C1_Note          = iprs$C1_note,
        C2_Statistical   = iprs$C2,
        C2_Note          = iprs$C2_note,
        C3_Abruptness    = iprs$C3,
        C3_Note          = iprs$C3_note,
        C4_GAM_Ribbon    = iprs$C4,
        C4_Note          = iprs$C4_note,
        C5_Smoothness    = iprs$C5,
        C5_Note          = iprs$C5_note,
        Rx_Amplitude     = iprs$amplitude,
        Rx_MeanRibbon    = iprs$mean_ribbon,
        Rx_TotalVar      = iprs$tv,
        Deriv_at_IP      = iprs$deriv_at_ip,
        IPRS             = iprs$IPRS,
        IPRS_Tier        = iprs$IPRS_tier,
        IPRS_Flag        = iprs$IPRS_flag,
        Method = "ZeroCrossing|BrentCI|IPRS=(C1+C2+C3+C4+C5)/5|C2=tanh(-log10p)|C5=waviness*sign_fidelity|SW_welch_t_OR_mannwhitney_u",
        stringsAsFactors = FALSE
    )
}))
rownames(ip_sheet) <- NULL

# Sheet 2: Residual predictions (downsampled 10×)
resid_pred_df <- do.call(rbind, lapply(tumor_samples, function(nm) {
    r <- ip_results_prostate[[nm]]
    if (is.null(r) || nrow(r$residual_pred) == 0L) return(NULL)
    p <- r$residual_pred[seq(1L, nrow(r$residual_pred), by = 10L), ]
    data.frame(Sample = nm, p, stringsAsFactors = FALSE)
}))
if (is.null(resid_pred_df)) resid_pred_df <- data.frame()
rownames(resid_pred_df) <- NULL

# Sheet 3: Method notes
method_notes <- data.frame(
    Step = paste0("Step ", 1:8),
    Description = c(
        "Cell classification: EP = 1 - model_residual² / null_residual²; top (dev_expl × n) cells are 'purple'. Null model uses same cell weights as fitted GAM (dynamic zero-weighting if available) for consistent EP.",
        "Signed residuals: residual = observed Expression − GAM_fitted, for each purple cell.",
        "Residual GAM R(x): thin-plate spline (REML) fitted to (TRPM4, residual) for purple cells; k = min(20, n_unique−1). Unweighted — all purple cells contribute equally to R(x).",
        "Fine grid: R(x) and SE evaluated at 2000 equally-spaced TRPM4 values.",
        "IP = first neg→pos zero-crossing of R(x) GAM (linear interpolation). If R(x) has no neg→pos crossing, IP=NA and all IPRS components default to 0.",
        paste0("95% CI primary (Brent): CI_Lower = zero of [R_hat + 1.96×SE]; CI_Upper = zero of [R_hat - 1.96×SE]; anchored near IP. ",
               "Fallback (proximity): if Brent returns NA for one bound, fallback finds x-value where band is closest to zero in the pre-IP (lower) or post-IP (upper) region. ",
               "CI_Lower_Method / CI_Upper_Method column records 'brent' or 'proximity_fallback'."),
        paste0("Stat test: Shapiro-Wilk (α=", SHAPIRO_ALPHA,
               ") per group; both normal → Welch t (one-sided post>pre); either non-normal → Mann-Whitney U; BH FDR applied."),
        paste0("IPRS = (C1 + C2 + C3 + C4 + C5) / 5  [range 0-1]. ",
               "C1 (Precision): both CI bounds → max(0, 1-CI_Width/TRPM4_range); one bound NA → 0.5×max(0, 1-half_width/TRPM4_range); both NA → 0. ",
               "C2 (Statistical) = tanh(scale_c2 * (-log10(p_adj))); scale_c2 = atanh(0.5)/(-log10(0.05)); C2=0.5 at standard significance threshold. ",
               "C3 (Abruptness) = tanh(|R'(IP)|) via central difference on 2000-pt grid. ",
               "C4 (GAM Ribbon) = exp(-mean_ribbon_width / amplitude). ",
               "C5 (Smoothness) = exp(-(TV/amplitude - 1)) × (1-neg_post_frac)²; C5=0 if amplitude < 1e-6. ",
               "Tiers: ≥0.80 Very Strong | 0.60-0.79 Strong | 0.40-0.59 Moderate [visual review] | 0.20-0.39 Weak [visual review] | <0.20 Failed [no IP].")
    ),
    stringsAsFactors = FALSE
)

# Write workbook
hdr_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#1F3864",
                         textDecoration = "bold", fontSize = 10,
                         halign = "center", valign = "center",
                         border = "TopBottomLeftRight",
                         borderStyle = "thin", borderColour = "#AAAAAA",
                         wrapText = TRUE)

wb <- createWorkbook()

addWorksheet(wb, "IP_Results")
writeData(wb, "IP_Results", ip_sheet)
addStyle(wb, "IP_Results", hdr_style,
         rows = 1, cols = seq_len(ncol(ip_sheet)), gridExpand = TRUE)
setColWidths(wb, "IP_Results", cols = 1:4,   widths = 18)
setColWidths(wb, "IP_Results", cols = 5:16,  widths = 16)
setColWidths(wb, "IP_Results", cols = 17:30, widths = 16)
setColWidths(wb, "IP_Results", cols = 31,    widths = 90)
setRowHeights(wb, "IP_Results", rows = 1, heights = 28)
freezePane(wb, "IP_Results", firstActiveRow = 2L, firstActiveCol = 3L)

addWorksheet(wb, "Residual_Predictions")
writeData(wb, "Residual_Predictions", resid_pred_df)
if (nrow(resid_pred_df) > 0) {
    addStyle(wb, "Residual_Predictions", hdr_style,
             rows = 1, cols = seq_len(ncol(resid_pred_df)), gridExpand = TRUE)
    setColWidths(wb, "Residual_Predictions",
                 cols = seq_len(ncol(resid_pred_df)), widths = 16)
    setRowHeights(wb, "Residual_Predictions", rows = 1, heights = 28)
}

addWorksheet(wb, "Method_Steps")
writeData(wb, "Method_Steps", method_notes)
addStyle(wb, "Method_Steps", hdr_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "Method_Steps", cols = 1, widths = 10)
setColWidths(wb, "Method_Steps", cols = 2, widths = 130)
setRowHeights(wb, "Method_Steps", rows = 1, heights = 28)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
filename  <- paste0("TRPM4_Ribo_IP_IPRS_Results_", timestamp, ".xlsx")
saveWorkbook(wb, filename, overwrite = TRUE)
cat("Exported:", filename, "\n")

# ------------------------------------------------------------------
# SECTION 11: Diagnostic plots
# ------------------------------------------------------------------
ip_line_color <- "#2E8B57"
manual_color  <- "#CC0000"

iprs_subtitle <- function(r) {
    iprs <- r$iprs
    sprintf("IPRS=%.3f (%s)  |  C1=%.3f  C2=%.3f  C3=%.3f  C4=%.3f  C5=%.3f  |  %s",
            iprs$IPRS, iprs$IPRS_tier,
            iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5,
            if (iprs$IPRS_flag == "None") "OK" else iprs$IPRS_flag)
}

# ---- Scatter plot (purple/gray + GAM line + IP vline) ----
create_scatter_prostate <- function(current_sample) {
    r <- ip_results_prostate[[current_sample]]
    if (is.null(r)) return(NULL)

    sample_gam_model <- pca_results[[current_sample]][["Ribo"]]$best_model
    sample_data      <- pca_results[[current_sample]][["Ribo"]]$gam_data
    dev_explained    <- summary(sample_gam_model)$dev.expl

    cluster_cells  <- WhichCells(prostate_results$seurat_obj, idents = PCA_CLUSTERS)
    sample_cells   <- WhichCells(prostate_results$seurat_obj,
                                 cells = grep(current_sample,
                                              colnames(prostate_results$seurat_obj),
                                              value = TRUE))
    common_cells   <- intersect(rownames(sample_data),
                                intersect(cluster_cells, sample_cells))
    if (length(common_cells) == 0) return(NULL)
    sample_data <- sample_data[rownames(sample_data) %in% common_cells, ]

    cell_weights <- if (!is.null(pca_results[[current_sample]][["Ribo"]]$cell_weights)) {
        wts <- pca_results[[current_sample]][["Ribo"]]$cell_weights
        wts[rownames(sample_data)]
    } else rep(1, nrow(sample_data))

    null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data, weights = cell_weights))
    model_fitted <- fitted(sample_gam_model)
    expl_power   <- 1 - ((sample_data$Expression - model_fitted)^2 /
                             (sample_data$Expression - null_fitted)^2)
    dev_cells    <- rownames(sample_data)[
        order(expl_power, decreasing = TRUE)[
            seq_len(round(nrow(sample_data) * dev_explained))]]

    plot_data <- data.frame(
        TRPM4      = sample_data$TRPM4,
        Expression = sample_data$Expression,
        Group      = ifelse(rownames(sample_data) %in% dev_cells,
                            "Dev explained", "Non-dev explained")
    )
    plot_data <- plot_data[order(plot_data$Group == "Dev explained"), ]

    pred_data        <- data.frame(TRPM4 = seq(min(plot_data$TRPM4),
                                               max(plot_data$TRPM4),
                                               length.out = 100L))
    pg               <- predict(sample_gam_model, newdata = pred_data, se.fit = TRUE)
    pred_data$fit    <- pg$fit
    pred_data$se.fit <- pg$se.fit

    title_col <- if (r$manual_flag) manual_color else "black"
    title <- if (r$manual_flag || is.na(r$ip_value)) {
        paste0(current_sample, "  *** NO IP DETECTED ***")
    } else {
        paste0(current_sample,
               "  IP: ", sprintf("%.4f", r$ip_value),
               "  95% CI [",
               ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
               ", ",
               ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
               "]")
    }

    lc <- ip_line_color  # "#2E8B57"

    p <- ggplot()

    # Green shaded CI band (drawn first, behind everything)
    if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
        p <- p + annotate("rect",
                          xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                          ymin = -Inf, ymax = Inf,
                          fill = lc, alpha = 0.08)

    p <- p +
        geom_point(data = plot_data,
                   aes(x = TRPM4, y = Expression, color = Group),
                   size = 1.8, alpha = 0.4) +
        geom_ribbon(data = pred_data,
                    aes(x = TRPM4, ymin = fit - 1.96 * se.fit,
                        ymax = fit + 1.96 * se.fit),
                    fill = "#FFCC99", alpha = 0.20) +
        geom_line(data = pred_data, aes(x = TRPM4, y = fit),
                  color = "#FFCC99", linewidth = 1.2)

    # Green dot-dash CI boundary lines
    if (!is.na(r$ip_ci_lower))
        p <- p + geom_vline(xintercept = r$ip_ci_lower,
                            linetype = "dotdash", color = lc, linewidth = 0.7)
    if (!is.na(r$ip_ci_upper))
        p <- p + geom_vline(xintercept = r$ip_ci_upper,
                            linetype = "dotdash", color = lc, linewidth = 0.7)

    # Green dashed IP line (on top)
    if (!is.na(r$ip_value))
        p <- p + geom_vline(xintercept = r$ip_value,
                            linetype = "dashed", color = lc, linewidth = 1.0)

    p + scale_color_manual(
        values = c("Dev explained"     = "#4B0082",
                   "Non-dev explained" = "#C0C0C0")) +
        theme_minimal() +
        theme(panel.grid.major = element_line(color = "#EEEEEE"),
              panel.grid.minor = element_line(color = "#F5F5F5"),
              legend.position  = "none",
              plot.title       = element_text(size = 10, face = "bold",
                                              hjust = 0.5, color = title_col),
              plot.subtitle    = element_text(size = 8, hjust = 0.5),
              panel.border     = element_rect(color = "black", fill = NA,
                                              linewidth = 0.5),
              axis.line        = element_blank()) +
        labs(title    = title,
             subtitle = paste0("Dev explained: ", round(dev_explained * 100, 2),
                               "%\n", iprs_subtitle(r)),
             x = "TRPM4 Expression",
             y = "Ribo Expression") +
        scale_x_continuous(
            limits = range(sample_data$TRPM4, na.rm = TRUE) +
                c(-1, 1) * 0.05 * diff(range(sample_data$TRPM4, na.rm = TRUE)))
}

# ---- R(x) diagnostic plot ----
create_rx_prostate <- function(current_sample) {
    r <- ip_results_prostate[[current_sample]]
    if (is.null(r)) return(NULL)

    pred      <- r$residual_pred
    rd        <- r$residual_data
    lc        <- ip_line_color
    title_col <- if (r$manual_flag) manual_color else "black"

    # Bin residuals for bar display
    n_bins_rx <- 100L
    p_trpm4   <- rd$TRPM4
    if (nrow(rd) > 0L) {
        if (length(unique(p_trpm4)) > n_bins_rx) {
            breaks_rx  <- seq(min(p_trpm4), max(p_trpm4), length.out = n_bins_rx + 1L)
            bin_ids_rx <- cut(p_trpm4, breaks = breaks_rx, include.lowest = TRUE, labels = FALSE)
            bin_ctr_rx <- (breaks_rx[-length(breaks_rx)] + breaks_rx[-1L]) / 2
            bin_w_rx   <- breaks_rx[2L] - breaks_rx[1L]
        } else {
            uq_rx      <- sort(unique(p_trpm4))
            bin_ids_rx <- match(p_trpm4, uq_rx)
            bin_ctr_rx <- uq_rx
            bin_w_rx   <- if (length(uq_rx) > 1L) min(diff(uq_rx)) else 0.1
        }
        bs_rx <- data.frame(
            foxm1_center  = bin_ctr_rx,
            mean_residual = sapply(seq_along(bin_ctr_rx), function(b) {
                idx <- which(bin_ids_rx == b)
                if (length(idx) > 0L) mean(rd$residual[idx]) else NA_real_
            })
        )
        bs_rx    <- bs_rx[!is.na(bs_rx$mean_residual), ]
        bar_w_rx <- if (nrow(bs_rx) > 1L)
            diff(range(bs_rx$foxm1_center)) / nrow(bs_rx) * 0.8 else bin_w_rx * 0.8
    } else {
        bs_rx    <- data.frame(foxm1_center = numeric(0), mean_residual = numeric(0))
        bar_w_rx <- 0.1
    }

    stat_txt <- if (!is.na(r$stat_test)) {
        sprintf("%s | BH FDR q = %.3e", r$stat_test,
                ifelse(is.na(r$stat_p_adj), NA, r$stat_p_adj))
    } else "Stat: NA"

    subtitle <- if (r$status == "insufficient_purple_cells") {
        paste0("*** FAILED — R(x) GAM not fitted (", nrow(rd),
               " purple cells, ", length(unique(rd$TRPM4)), " unique TRPM4 values",
               " — minimum ", MIN_PURPLE_CELLS, " cells / 3 unique values required) ***")
    } else if (r$manual_flag) {
        paste0("*** NO IP — R(x) has no neg\u2192pos crossing ***\n", iprs_subtitle(r))
    } else {
        paste0("IP = ", sprintf("%.4f", r$ip_value),
               "  95% CI [",
               ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
               ", ",
               ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
               "]  |  ", stat_txt, "\n", iprs_subtitle(r))
    }

    p <- ggplot()

    if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
        p <- p + annotate("rect",
                          xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                          ymin = -Inf, ymax = Inf, fill = lc, alpha = 0.08)

    if (nrow(bs_rx) > 0L)
        p <- p + geom_bar(data = bs_rx,
                          aes(x = foxm1_center, y = mean_residual),
                          stat = "identity", fill = "#4B0082", color = "#3D0065",
                          alpha = 0.40, width = bar_w_rx, linewidth = 0.4)

    if (nrow(pred) > 0L)
        p <- p +
            geom_ribbon(data = pred,
                        aes(x = TRPM4, ymin = r_lower, ymax = r_upper),
                        fill = "#FFCC99", alpha = 0.35) +
            geom_line(data = pred, aes(x = TRPM4, y = r_hat),
                      color = "#FFCC99", linewidth = 1.5)

    p <- p + geom_hline(yintercept = 0, linetype = "dotted",
                        color = "gray40", linewidth = 0.5)

    if (!is.na(r$ip_ci_lower))
        p <- p + geom_vline(xintercept = r$ip_ci_lower,
                            linetype = "dotdash", color = lc, linewidth = 0.7)
    if (!is.na(r$ip_ci_upper))
        p <- p + geom_vline(xintercept = r$ip_ci_upper,
                            linetype = "dotdash", color = lc, linewidth = 0.7)
    if (!is.na(r$ip_value))
        p <- p + geom_vline(xintercept = r$ip_value,
                            linetype = "dashed", color = lc, linewidth = 1.0)

    p + theme_minimal() +
        theme(panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
              axis.line     = element_blank(),
              plot.title    = element_text(size = 10, face = "bold",
                                           hjust = 0.5, color = title_col),
              plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "black"),
              axis.text     = element_text(size = 7)) +
        labs(title    = paste0(current_sample,
                               " \u2014 R(x) = E[residual | purple, TRPM4=x]",
                               "  [Zero-crossing | Brent CI | IPRS]"),
             subtitle = subtitle,
             x = "TRPM4 Expression",
             y = "Expected Signed Residual R(x)")
}

# ---- Generate all diagnostic plots ----
cat("\n\nGenerating scatter plots (TREP/non-TREP + IP)...\n")
cat(rep("=", 60), "\n")
for (nm in tumor_samples) {
    tryCatch({ p <- create_scatter_prostate(nm); if (!is.null(p)) print(p) },
             error = function(e) cat("Scatter error:", nm, ":", conditionMessage(e), "\n"))
}

cat("\nGenerating R(x) diagnostic plots...\n")
cat(rep("=", 60), "\n")
for (nm in tumor_samples) {
    tryCatch({ p <- create_rx_prostate(nm); if (!is.null(p)) print(p) },
             error = function(e) cat("R(x) error:", nm, ":", conditionMessage(e), "\n"))
}

# ---- IPRS tier summary ----
tier_counts <- table(sapply(ip_results_prostate, function(r) r$iprs$IPRS_tier))
cat(sprintf("\nIPRS tier summary — Ribo composite (%d samples):\n", length(ip_results_prostate)))
for (tier in c("Very Strong", "Strong", "Moderate", "Weak", "Failed"))
    cat(sprintf("  %-12s : %d\n", tier,
                if (tier %in% names(tier_counts)) tier_counts[[tier]] else 0L))

# ------------------------------------------------------------------
# SECTION 12: Final environment confirmation
# ------------------------------------------------------------------
cat("\n========== Part 3.09b complete ==========\n")
cat("\nObjects available for Parts 3.10 to 3.15:\n")
cat("  ip_values           — named numeric vector, length", length(ip_values), "\n")
cat("  inflection_points   — named numeric vector (alias), length", length(inflection_points), "\n")
cat("  ip_results_prostate — list of full IP results per sample\n")
cat("  ip_table_prostate   — data.frame with IPRS and all fields\n")
cat("\nParts 3.12, 3.13 and 3.15 check: if (!exists('ip_values')) → PASSES\n")
cat("Part 3.13 check: inflection_points <- ip_values[patient_ids] → READY\n")
