# ============================================================================
# Calibration helper for Scenario 4E confounding_strength
# ============================================================================
# Purpose:
#   Run a small data-generation grid over confounding_strength values for
#   S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_CONTINUOUS_X2_T100.
#
# Selection principle:
#   Choose the weakest confounding_strength that passes the S4E data checks.
#   This gives strong spatial/covariate confounding without making x2 nearly
#   deterministic from phi or accidentally creating another stress mechanism
#   such as sparsity or extreme counts.
#
# Usage:
#   source("s4e_spatial_covariate_confounding_continuous_x2.R")
#   source("calibrate_s4e_confounding_strength_continuous_x2.R")
#
#   s4e_cal <- calibrate_s4e_confounding_strength_continuous_x2(
#       root = ".",
#       strengths = c(0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00),
#       reps = 1:5,
#       overwrite_existing = TRUE,
#       verbose = TRUE
#   )
#
#   s4e_cal$summary
#   s4e_cal$recommended_strength
#
# After choosing the recommended strength, rerun official 10-rep generation
# using run_s4e_continuous_x2_data_generation(...).
# ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.s4e_strength_label <- function(x) {
    gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

.ensure_s4e_loaded <- function() {
    required <- c(
        "run_s4e_continuous_x2_data_generation",
        "check_s4e_continuous_x2_data_summary"
    )
    missing <- required[!vapply(required, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop(
            "Required S4E functions are not loaded: ",
            paste(missing, collapse = ", "),
            ". Source s4e_spatial_covariate_confounding_continuous_x2.R first.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.score_s4e_calibration <- function(df,
                                   target_abs_cell_cor = 0.65,
                                   target_abs_area_cor = 0.90,
                                   target_mean_count = 200,
                                   target_zero_prop = 0.02) {
    pass_penalty <- ifelse(isTRUE(df$passes_s4e_data_check), 0, 1000)

    count_penalty <- abs(df$mean_count_avg - target_mean_count) / target_mean_count
    zero_penalty <- abs(df$zero_prop_avg - target_zero_prop)
    cell_penalty <- abs(df$x2_phi_abs_cell_cor_avg - target_abs_cell_cor)
    area_penalty <- abs(df$x2_phi_abs_area_mean_cor_avg - target_abs_area_cor)

    # Additional penalty if the confounded x2 becomes too spatially deterministic.
    too_deterministic_penalty <- ifelse(
        is.finite(df$x2_phi_abs_cell_cor_avg) && df$x2_phi_abs_cell_cor_avg > 0.90,
        2 * (df$x2_phi_abs_cell_cor_avg - 0.90),
        0
    )

    pass_penalty +
        0.25 * count_penalty +
        1.00 * zero_penalty +
        1.00 * cell_penalty +
        0.50 * area_penalty +
        too_deterministic_penalty
}

calibrate_s4e_confounding_strength_continuous_x2 <- function(
        root = ".",
        strengths = c(0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00),
        reps = 1:5,
        TT_use = 100,
        beta0_reference_truth = -1.5,
        r_truth = 15,
        residual_weight = 1.00,
        preserve_reference_mean = TRUE,
        target_mean_count_range = c(80, 350),
        target_zero_prop_max = 0.10,
        minimum_abs_cell_cor = 0.50,
        minimum_abs_area_mean_cor = 0.80,
        target_abs_beta0_ident_max = 20,
        max_count_max_limit = Inf,
        data_dir = "calibration_s4e_confounding_strength_continuous_x2",
        overwrite_existing = TRUE,
        verbose = TRUE,
        ...) {

    .ensure_s4e_loaded()

    if (!is.numeric(strengths) || any(!is.finite(strengths)) || any(strengths <= 0)) {
        stop("strengths must be positive finite numeric values.", call. = FALSE)
    }
    strengths <- sort(unique(as.numeric(strengths)))

    rows <- list()
    manifests <- list()

    for (ss in strengths) {
        lab <- .s4e_strength_label(ss)
        scenario_id <- paste0(
            "S4E_CALIBRATION_CONFOUNDING_STRENGTH_",
            lab,
            "_CONTINUOUS_X2_T100"
        )

        if (isTRUE(verbose)) {
            message("\n--- Calibrating S4E confounding_strength = ", ss, " ---")
        }

        manifest <- run_s4e_continuous_x2_data_generation(
            root = root,
            reps = reps,
            TT_use = TT_use,
            beta0_reference_truth = beta0_reference_truth,
            r_truth = r_truth,
            confounding_strength = ss,
            residual_weight = residual_weight,
            preserve_reference_mean = preserve_reference_mean,
            data_dir = data_dir,
            scenario_id = scenario_id,
            overwrite_existing = overwrite_existing,
            verbose = verbose,
            x2_mode = "continuous_time",
            x2_ar = 0.50,
            x2_innov_sd = 0.80,
            ...
        )

        chk <- check_s4e_continuous_x2_data_summary(
            manifest,
            target_mean_count_range = target_mean_count_range,
            target_zero_prop_max = target_zero_prop_max,
            minimum_abs_cell_cor = minimum_abs_cell_cor,
            minimum_abs_area_mean_cor = minimum_abs_area_mean_cor,
            target_abs_beta0_ident_max = target_abs_beta0_ident_max,
            max_count_max_limit = max_count_max_limit
        )

        chk$confounding_strength <- ss
        chk$calibration_scenario_id <- scenario_id
        chk$calibration_score <- .score_s4e_calibration(chk)

        rows[[lab]] <- chk
        manifests[[lab]] <- manifest
    }

    summary <- do.call(rbind, rows)
    rownames(summary) <- NULL

    # The official rule is the weakest strength that passes all checks.
    pass_idx <- which(summary$passes_s4e_data_check)
    if (length(pass_idx) > 0L) {
        recommended_idx <- pass_idx[order(summary$confounding_strength[pass_idx])][1L]
        recommendation_rule <- "weakest strength passing all S4E data checks"
    } else {
        recommended_idx <- which.min(summary$calibration_score)
        recommendation_rule <- "no strength passed all checks; selected lowest calibration score"
    }

    recommended_strength <- summary$confounding_strength[recommended_idx]

    compact_cols <- c(
        "confounding_strength", "n_reps", "mean_count_avg", "zero_prop_avg",
        "x2_phi_abs_cell_cor_avg", "x2_phi_abs_area_mean_cor_avg",
        "reference_x2_phi_abs_cell_cor_avg", "reference_x2_phi_abs_area_mean_cor_avg",
        "x2_empirical_ar1_mean_avg", "x2_sd_avg", "x2_binary_like_prop_max",
        "beta0_shift_avg", "max_abs_beta0_ident", "passes_s4e_data_check",
        "calibration_score"
    )
    compact_cols <- compact_cols[compact_cols %in% names(summary)]
    compact_summary <- summary[, compact_cols, drop = FALSE]

    if (isTRUE(verbose)) {
        message("\n================ S4E confounding-strength calibration ================")
        print(compact_summary)
        message("\nRecommended confounding_strength: ", recommended_strength)
        message("Selection rule: ", recommendation_rule)
        message("=====================================================================")
    }

    list(
        summary = summary,
        compact_summary = compact_summary,
        manifests = manifests,
        recommended_strength = recommended_strength,
        recommended_row = summary[recommended_idx, , drop = FALSE],
        recommendation_rule = recommendation_rule,
        settings = list(
            strengths = strengths,
            reps = reps,
            TT_use = TT_use,
            beta0_reference_truth = beta0_reference_truth,
            r_truth = r_truth,
            residual_weight = residual_weight,
            preserve_reference_mean = preserve_reference_mean,
            target_mean_count_range = target_mean_count_range,
            target_zero_prop_max = target_zero_prop_max,
            minimum_abs_cell_cor = minimum_abs_cell_cor,
            minimum_abs_area_mean_cor = minimum_abs_area_mean_cor
        )
    )
}
