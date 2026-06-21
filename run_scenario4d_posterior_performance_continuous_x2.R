## ============================================================================
## run_scenario4d_posterior_performance_continuous_x2.R
##
## Posterior-performance summary for Scenario 4D:
## short-time-series / small-T stress test with gamma fixed.
##
## Design goal
##   Scenario 4D keeps the S3 mean, exposure, continuous-time covariate,
##   spatial, dispersion, and temporal latent-risk data-generating mechanism
##   unchanged, but shortens the temporal horizon from T = 100 to T = 30.
##   This isolates loss of temporal replication from the sparse-count,
##   low-exposure, and strong-overdispersion stresses studied in S4A--S4C.
##
## Typical use from project root:
##   source("run_scenario4d_posterior_performance_continuous_x2.R")
##
## Or programmatically:
##   out <- summarize_scenario4d_posterior_performance_continuous_x2(root = ".")
##
## Main input
##   output_s4d_short_time_fixed_gamma_continuous_x2/
##     S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2/
##       summary_S4D_short_time_continuous_x2_fixed_gamma_all_reps.csv
##
## Main output
##   analysis_s4d_short_time_continuous_x2/
##     S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2/tables/
##     S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2/figures/
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s4dctx_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s4dctx_ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

.s4dctx_read_required_csv <- function(path) {
    if (!file.exists(path)) {
        stop("Required CSV not found: ", path, call. = FALSE)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
}

.s4dctx_read_optional_csv <- function(path) {
    if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else NULL
}

.s4dctx_write_csv <- function(x, path) {
    .s4dctx_ensure_dir(dirname(path))
    utils::write.csv(x, path, row.names = FALSE)
    invisible(path)
}

.s4dctx_num <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    suppressWarnings(as.numeric(x[[1L]]))
}

.s4dctx_col <- function(d, nm, default = NA_real_) {
    if (!nm %in% names(d)) return(rep(default, nrow(d)))
    d[[nm]]
}

.s4dctx_first_existing_col <- function(d, candidates, default = NA_real_) {
    for (nm in candidates) {
        if (nm %in% names(d)) return(d[[nm]])
    }
    rep(default, nrow(d))
}

.s4dctx_safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.s4dctx_safe_sd <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (sum(is.finite(x)) <= 1L) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

.s4dctx_safe_sum <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(0)
    sum(x, na.rm = TRUE)
}

.s4dctx_safe_min <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    min(x, na.rm = TRUE)
}

.s4dctx_safe_max <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    max(x, na.rm = TRUE)
}

.s4dctx_safe_median <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    stats::median(x, na.rm = TRUE)
}

.s4dctx_bind_rows_aligned <- function(rows) {
    rows <- rows[!vapply(rows, is.null, logical(1))]
    rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
    if (length(rows) == 0L) return(data.frame())
    all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
    rows2 <- lapply(rows, function(x) {
        for (cc in setdiff(all_cols, names(x))) x[[cc]] <- NA
        x[, all_cols, drop = FALSE]
    })
    ans <- do.call(rbind, rows2)
    rownames(ans) <- NULL
    ans
}

.s4dctx_fix_summary_aliases <- function(d) {
    ## Harmonize source-data names and fit-summary names.
    if (!"mean_count" %in% names(d) && "observed_mean_count" %in% names(d)) d$mean_count <- d$observed_mean_count
    if (!"median_count" %in% names(d) && "observed_median_count" %in% names(d)) d$median_count <- d$observed_median_count
    if (!"zero_prop" %in% names(d) && "observed_zero_prop" %in% names(d)) d$zero_prop <- d$observed_zero_prop
    if (!"total_count" %in% names(d) && "observed_total_count" %in% names(d)) d$total_count <- d$observed_total_count
    if (!"max_count" %in% names(d) && "observed_max_count" %in% names(d)) d$max_count <- d$observed_max_count

    if (!"TT" %in% names(d) && "TT_short" %in% names(d)) d$TT <- d$TT_short
    if (!"TT_short" %in% names(d) && "TT" %in% names(d)) d$TT_short <- d$TT
    if (!"TT_reference" %in% names(d)) d$TT_reference <- 100L
    if (!"TT_ratio_to_reference" %in% names(d) && all(c("TT_short", "TT_reference") %in% names(d))) {
        d$TT_ratio_to_reference <- suppressWarnings(as.numeric(d$TT_short) / as.numeric(d$TT_reference))
    }
    if (!"n_cells" %in% names(d) && all(c("TT_short", "n1") %in% names(d))) {
        d$n_cells <- suppressWarnings(as.numeric(d$TT_short) * as.numeric(d$n1))
    }
    if (!"n_cells" %in% names(d) && "n_cells_avg" %in% names(d)) d$n_cells <- d$n_cells_avg

    ## Truth aliases.  The S4D beta0 truth on the fitted identified scale may be rep-specific.
    if (!"beta0_true" %in% names(d) && "beta0_star_ident" %in% names(d)) d$beta0_true <- d$beta0_star_ident
    if (!"beta1_true" %in% names(d)) d$beta1_true <- 0.5
    if (!"beta2_true" %in% names(d)) d$beta2_true <- -0.4
    if (!"r_true_mean" %in% names(d) && "r_truth_mean_avg" %in% names(d)) d$r_true_mean <- d$r_truth_mean_avg
    if (!"r_true_mean" %in% names(d) && "r_truth" %in% names(d)) d$r_true_mean <- d$r_truth
    if (!"r_true_mean" %in% names(d)) d$r_true_mean <- 15

    ## Guard columns.
    for (nm in c(
        "s4d_beta_guard_count", "s4d_kappa_guard_count",
        "s4d_lambda_input_guard_count", "s4d_lambda_output_guard_count"
    )) {
        if (!nm %in% names(d)) d[[nm]] <- 0
    }

    ## Continuous-time x2 diagnostic columns.
    if (!"x2_mode" %in% names(d)) d$x2_mode <- NA_character_
    if (!"mean_x2_within_sd" %in% names(d)) d$mean_x2_within_sd <- NA_real_
    if (!"min_x2_within_sd" %in% names(d)) d$min_x2_within_sd <- NA_real_
    if (!"mean_lag1_x2" %in% names(d)) d$mean_lag1_x2 <- NA_real_

    ## Derived kappa-CV bias.
    if (!"kappa_cv_bias" %in% names(d) && all(c("kappa_post_mean_cv", "kappa_truth_cv") %in% names(d))) {
        d$kappa_cv_bias <- d$kappa_post_mean_cv - d$kappa_truth_cv
    }

    d
}

.s4dctx_classify_one <- function(row) {
    beta_guard <- .s4dctx_num(row$s4d_beta_guard_count, 0)
    kappa_guard <- .s4dctx_num(row$s4d_kappa_guard_count, 0)
    lambda_in_guard <- .s4dctx_num(row$s4d_lambda_input_guard_count, 0)
    lambda_out_guard <- .s4dctx_num(row$s4d_lambda_output_guard_count, 0)

    beta0_mean <- .s4dctx_num(row$beta0_mean, NA_real_)
    beta1_mean <- .s4dctx_num(row$beta1_mean, NA_real_)
    beta2_mean <- .s4dctx_num(row$beta2_mean, NA_real_)
    r_mean <- .s4dctx_num(row$r_mean, NA_real_)
    lambda_rmse <- .s4dctx_num(row$lambda_rmse, NA_real_)
    log_lambda_rmse <- .s4dctx_num(row$log_lambda_rmse, NA_real_)
    kappa_rmse <- .s4dctx_num(row$kappa_rmse, NA_real_)
    log_kappa_rmse <- .s4dctx_num(row$log_kappa_rmse, NA_real_)

    severe <- (!is.finite(beta0_mean)) || (!is.finite(beta1_mean)) || (!is.finite(beta2_mean)) ||
        (!is.finite(r_mean)) || (!is.finite(lambda_rmse)) || (!is.finite(log_lambda_rmse)) ||
        abs(beta0_mean) > 30 || abs(beta1_mean) > 5 || abs(beta2_mean) > 5 ||
        r_mean <= 0 || r_mean > 80 ||
        lambda_rmse > 5 || log_lambda_rmse > 2 ||
        (is.finite(kappa_rmse) && kappa_rmse > 2) ||
        (is.finite(log_kappa_rmse) && log_kappa_rmse > 2) ||
        beta_guard > 1000 || kappa_guard > 1000 || lambda_in_guard > 0 || lambda_out_guard > 1000

    if (isTRUE(severe)) return("numerical_instability")

    warn <- beta_guard > 0 || kappa_guard > 0 || lambda_out_guard > 0 ||
        (is.finite(lambda_rmse) && lambda_rmse > 0.5) ||
        (is.finite(log_lambda_rmse) && log_lambda_rmse > 0.5)
    if (isTRUE(warn)) return("soft_warning")

    "stable"
}

.s4dctx_add_fit_status <- function(d) {
    d <- .s4dctx_fix_summary_aliases(d)
    d$fit_status <- vapply(seq_len(nrow(d)), function(i) .s4dctx_classify_one(d[i, , drop = FALSE]), character(1))
    d$nonfailed <- d$fit_status %in% c("stable", "soft_warning")
    d$stable_fit <- d$fit_status == "stable"
    d
}

.s4dctx_check_continuous_x2_input <- function(d, require_continuous_x2 = TRUE) {
    if (!isTRUE(require_continuous_x2)) return(invisible(TRUE))

    x2_mode <- unique(stats::na.omit(as.character(.s4dctx_col(d, "x2_mode", NA_character_))))
    if (length(x2_mode) == 0L || any(x2_mode != "continuous_time")) {
        stop(
            "S4D posterior-performance input is not uniformly continuous-time x2. Found x2_mode = ",
            paste(x2_mode, collapse = ", "),
            call. = FALSE
        )
    }

    within_sd <- suppressWarnings(as.numeric(.s4dctx_col(d, "mean_x2_within_sd", NA_real_)))
    min_within_sd <- suppressWarnings(as.numeric(.s4dctx_col(d, "min_x2_within_sd", NA_real_)))
    lag1 <- suppressWarnings(as.numeric(.s4dctx_col(d, "mean_lag1_x2", NA_real_)))

    if (!all(is.finite(within_sd)) || any(within_sd <= 0, na.rm = TRUE)) {
        stop("S4D posterior-performance input failed continuous-time x2 check: mean_x2_within_sd must be finite and positive.", call. = FALSE)
    }
    if (any(is.finite(min_within_sd)) && any(min_within_sd <= 0, na.rm = TRUE)) {
        stop("S4D posterior-performance input failed continuous-time x2 check: min_x2_within_sd includes non-positive values.", call. = FALSE)
    }
    if (!any(is.finite(lag1))) {
        stop("S4D posterior-performance input failed continuous-time x2 check: mean_lag1_x2 is entirely non-finite.", call. = FALSE)
    }

    invisible(TRUE)
}

.s4dctx_fit_status_counts <- function(d) {
    statuses <- c("stable", "soft_warning", "numerical_instability")
    tab <- table(factor(d$fit_status, levels = statuses))
    data.frame(
        fit_status = statuses,
        status_label = c("Stable", "Soft warning", "Numerical instability"),
        n_reps = as.integer(tab),
        prop = as.numeric(tab) / nrow(d),
        stringsAsFactors = FALSE
    )
}

.s4dctx_subset_indices <- function(d) {
    list(
        all_sampled_lambda = rep(TRUE, nrow(d)),
        stable_only = d$fit_status == "stable",
        stable_plus_soft_warning = d$fit_status %in% c("stable", "soft_warning"),
        numerical_instability = d$fit_status == "numerical_instability"
    )
}

.s4dctx_summarize_subset <- function(d, idx, subset_name) {
    dd <- d[idx, , drop = FALSE]
    if (nrow(dd) == 0L) {
        return(data.frame(subset = subset_name, n = 0, stringsAsFactors = FALSE))
    }
    data.frame(
        subset = subset_name,
        n = nrow(dd),
        TT_short_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "TT_short")),
        TT_reference_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "TT_reference")),
        TT_ratio_to_reference_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "TT_ratio_to_reference")),
        n_cells_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "n_cells")),
        mean_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "mean_count")),
        mean_count_min = .s4dctx_safe_min(.s4dctx_col(dd, "mean_count")),
        mean_count_max = .s4dctx_safe_max(.s4dctx_col(dd, "mean_count")),
        median_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "median_count")),
        zero_prop_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "zero_prop")),
        total_count_sum = .s4dctx_safe_sum(.s4dctx_col(dd, "total_count")),
        total_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "total_count")),
        count_cv_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "count_cv")),
        count_cv_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "count_cv")),
        variance_to_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "variance_to_mean")),
        q95_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "q95_count")),
        q99_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "q99_count")),
        max_count_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "max_count")),
        max_count_max = .s4dctx_safe_max(.s4dctx_col(dd, "max_count")),
        beta0_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "beta0_mean")),
        beta0_mean_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "beta0_mean")),
        beta1_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "beta1_mean")),
        beta1_mean_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "beta1_mean")),
        beta2_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "beta2_mean")),
        beta2_mean_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "beta2_mean")),
        r_true_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "r_true_mean")),
        r_mean_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "r_mean")),
        r_mean_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "r_mean")),
        r_region_coverage_95_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "r_region_coverage_95")),
        r_region_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "r_region_rmse")),
        r_region_mae_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "r_region_mae")),
        phi_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "phi_rmse")),
        phi_mae_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "phi_mae")),
        phi_cor_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "phi_cor")),
        lambda_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "lambda_rmse")),
        log_lambda_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "log_lambda_rmse")),
        cor_log_lambda_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "cor_log_lambda")),
        cor_delta_log_lambda_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "cor_delta_log_lambda")),
        lambda_coverage_95_avg = .s4dctx_safe_mean(.s4dctx_first_existing_col(dd, c("lambda_coverage_95", "lambda_region_coverage_95"))),
        log_lambda_sd_truth_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "log_lambda_sd_truth")),
        delta_log_lambda_sd_truth_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "delta_log_lambda_sd_truth")),
        delta_log_lambda_abs_mean_truth_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "delta_log_lambda_abs_mean_truth")),
        x2_mode_unique = paste(unique(stats::na.omit(as.character(.s4dctx_col(dd, "x2_mode", NA_character_)))), collapse = ";"),
        mean_x2_within_sd_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "mean_x2_within_sd")),
        min_x2_within_sd_min = .s4dctx_safe_min(.s4dctx_col(dd, "min_x2_within_sd")),
        mean_lag1_x2_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "mean_lag1_x2")),
        mean_lag1_x2_sd = .s4dctx_safe_sd(.s4dctx_col(dd, "mean_lag1_x2")),
        kappa_truth_cv_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "kappa_truth_cv")),
        kappa_post_mean_cv_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "kappa_post_mean_cv")),
        kappa_cv_bias_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "kappa_cv_bias")),
        kappa_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "kappa_rmse")),
        log_kappa_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "log_kappa_rmse")),
        cor_kappa_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "cor_kappa")),
        cor_log_kappa_avg = .s4dctx_safe_mean(.s4dctx_col(dd, "cor_log_kappa")),
        beta_guard_total = .s4dctx_safe_sum(.s4dctx_col(dd, "s4d_beta_guard_count", 0)),
        kappa_guard_total = .s4dctx_safe_sum(.s4dctx_col(dd, "s4d_kappa_guard_count", 0)),
        lambda_input_guard_total = .s4dctx_safe_sum(.s4dctx_col(dd, "s4d_lambda_input_guard_count", 0)),
        lambda_output_guard_total = .s4dctx_safe_sum(.s4dctx_col(dd, "s4d_lambda_output_guard_count", 0)),
        stringsAsFactors = FALSE
    )
}

.s4dctx_make_fit_status_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        TT_short = as.integer(.s4dctx_col(d, "TT_short")),
        n_cells = as.numeric(.s4dctx_col(d, "n_cells")),
        mean_count = as.numeric(.s4dctx_col(d, "mean_count")),
        zero_prop = as.numeric(.s4dctx_col(d, "zero_prop")),
        count_cv = as.numeric(.s4dctx_col(d, "count_cv")),
        x2_mode = as.character(.s4dctx_col(d, "x2_mode", NA_character_)),
        mean_x2_within_sd = as.numeric(.s4dctx_col(d, "mean_x2_within_sd")),
        min_x2_within_sd = as.numeric(.s4dctx_col(d, "min_x2_within_sd")),
        mean_lag1_x2 = as.numeric(.s4dctx_col(d, "mean_lag1_x2")),
        beta0_mean = as.numeric(.s4dctx_col(d, "beta0_mean")),
        beta1_mean = as.numeric(.s4dctx_col(d, "beta1_mean")),
        beta2_mean = as.numeric(.s4dctx_col(d, "beta2_mean")),
        r_mean = as.numeric(.s4dctx_col(d, "r_mean")),
        log_lambda_rmse = as.numeric(.s4dctx_col(d, "log_lambda_rmse")),
        cor_log_lambda = as.numeric(.s4dctx_col(d, "cor_log_lambda")),
        cor_delta_log_lambda = as.numeric(.s4dctx_col(d, "cor_delta_log_lambda")),
        phi_rmse = as.numeric(.s4dctx_col(d, "phi_rmse")),
        phi_cor = as.numeric(.s4dctx_col(d, "phi_cor")),
        log_kappa_rmse = as.numeric(.s4dctx_col(d, "log_kappa_rmse")),
        cor_log_kappa = as.numeric(.s4dctx_col(d, "cor_log_kappa")),
        beta_guard = as.numeric(.s4dctx_col(d, "s4d_beta_guard_count", 0)),
        kappa_guard = as.numeric(.s4dctx_col(d, "s4d_kappa_guard_count", 0)),
        lambda_input_guard = as.numeric(.s4dctx_col(d, "s4d_lambda_input_guard_count", 0)),
        lambda_output_guard = as.numeric(.s4dctx_col(d, "s4d_lambda_output_guard_count", 0)),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_beta_long <- function(d) {
    params <- c("beta0", "beta1", "beta2")
    rows <- lapply(params, function(p) {
        mean_col <- paste0(p, "_mean")
        sd_col <- paste0(p, "_sd")
        q025_col <- paste0(p, "_q025")
        q50_col <- paste0(p, "_q50")
        q975_col <- paste0(p, "_q975")
        true_col <- paste0(p, "_true")
        covered_col <- paste0(p, "_covered")
        truth <- .s4dctx_col(d, true_col)
        post_mean <- .s4dctx_col(d, mean_col)
        q025 <- .s4dctx_col(d, q025_col)
        q975 <- .s4dctx_col(d, q975_col)
        covered <- .s4dctx_col(d, covered_col)
        if (all(is.na(covered)) && any(is.finite(q025)) && any(is.finite(q975))) {
            covered <- as.integer(q025 <= truth & truth <= q975)
        }
        data.frame(
            scenario_id = d$scenario_id,
            rep_id = as.integer(d$rep_id),
            parameter = p,
            truth = as.numeric(truth),
            post_mean = as.numeric(post_mean),
            post_sd = as.numeric(.s4dctx_col(d, sd_col)),
            q025 = as.numeric(q025),
            q50 = as.numeric(.s4dctx_col(d, q50_col)),
            q975 = as.numeric(q975),
            bias = as.numeric(post_mean - truth),
            abs_bias = abs(as.numeric(post_mean - truth)),
            covered = as.numeric(covered),
            fit_status = d$fit_status,
            nonfailed = d$nonfailed,
            stable_fit = d$stable_fit,
            stringsAsFactors = FALSE
        )
    })
    .s4dctx_bind_rows_aligned(rows)
}

.s4dctx_beta_aggregate <- function(beta_df, subset_name) {
    if (is.null(beta_df) || nrow(beta_df) == 0L) return(data.frame())
    rows <- lapply(split(beta_df, beta_df$parameter), function(dd) {
        data.frame(
            subset = subset_name,
            parameter = unique(dd$parameter),
            n_reps = length(unique(dd$rep_id)),
            mean_truth = .s4dctx_safe_mean(dd$truth),
            mean_post_mean = .s4dctx_safe_mean(dd$post_mean),
            sd_post_mean_across_reps = .s4dctx_safe_sd(dd$post_mean),
            mean_bias = .s4dctx_safe_mean(dd$bias),
            median_bias = .s4dctx_safe_median(dd$bias),
            rmse_bias = sqrt(mean(dd$bias^2, na.rm = TRUE)),
            mean_abs_bias = .s4dctx_safe_mean(abs(dd$bias)),
            coverage = .s4dctx_safe_mean(dd$covered),
            mean_post_sd = .s4dctx_safe_mean(dd$post_sd),
            stringsAsFactors = FALSE
        )
    })
    .s4dctx_bind_rows_aligned(rows)
}

.s4dctx_lambda_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        TT_short = as.integer(.s4dctx_col(d, "TT_short")),
        TT_ratio_to_reference = as.numeric(.s4dctx_col(d, "TT_ratio_to_reference")),
        n_cells = as.numeric(.s4dctx_col(d, "n_cells")),
        lambda_rmse = as.numeric(.s4dctx_col(d, "lambda_rmse")),
        log_lambda_rmse = as.numeric(.s4dctx_col(d, "log_lambda_rmse")),
        cor_log_lambda = as.numeric(.s4dctx_col(d, "cor_log_lambda")),
        cor_delta_log_lambda = as.numeric(.s4dctx_col(d, "cor_delta_log_lambda")),
        lambda_coverage_95 = as.numeric(.s4dctx_first_existing_col(d, c("lambda_coverage_95", "lambda_region_coverage_95"))),
        log_lambda_sd_truth = as.numeric(.s4dctx_col(d, "log_lambda_sd_truth")),
        delta_log_lambda_sd_truth = as.numeric(.s4dctx_col(d, "delta_log_lambda_sd_truth")),
        lambda_raw_median = as.numeric(.s4dctx_col(d, "lambda_raw_median")),
        lambda_raw_max = as.numeric(.s4dctx_col(d, "lambda_raw_max")),
        x2_mode = as.character(.s4dctx_col(d, "x2_mode", NA_character_)),
        mean_x2_within_sd = as.numeric(.s4dctx_col(d, "mean_x2_within_sd")),
        min_x2_within_sd = as.numeric(.s4dctx_col(d, "min_x2_within_sd")),
        mean_lag1_x2 = as.numeric(.s4dctx_col(d, "mean_lag1_x2")),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_phi_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        TT_short = as.integer(.s4dctx_col(d, "TT_short")),
        n_cells = as.numeric(.s4dctx_col(d, "n_cells")),
        phi_rmse = as.numeric(.s4dctx_col(d, "phi_rmse")),
        phi_mae = as.numeric(.s4dctx_col(d, "phi_mae")),
        phi_cor = as.numeric(.s4dctx_col(d, "phi_cor")),
        phi_ident_sd = as.numeric(.s4dctx_col(d, "phi_ident_sd")),
        tau_phi_mean = as.numeric(.s4dctx_col(d, "tau_phi_mean")),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_r_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        r_truth = as.numeric(.s4dctx_col(d, "r_true_mean")),
        r_mean = as.numeric(.s4dctx_col(d, "r_mean")),
        r_q025 = as.numeric(.s4dctx_col(d, "r_q025")),
        r_q50 = as.numeric(.s4dctx_col(d, "r_q50")),
        r_q975 = as.numeric(.s4dctx_col(d, "r_q975")),
        r_bias = as.numeric(.s4dctx_col(d, "r_mean") - .s4dctx_col(d, "r_true_mean")),
        r_abs_bias = abs(as.numeric(.s4dctx_col(d, "r_mean") - .s4dctx_col(d, "r_true_mean"))),
        r_covered = as.numeric(.s4dctx_col(d, "r_covered")),
        r_region_coverage_95 = as.numeric(.s4dctx_col(d, "r_region_coverage_95")),
        r_region_rmse = as.numeric(.s4dctx_col(d, "r_region_rmse")),
        r_region_mae = as.numeric(.s4dctx_col(d, "r_region_mae")),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_kappa_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        kappa_truth_cv = as.numeric(.s4dctx_col(d, "kappa_truth_cv")),
        kappa_post_mean_cv = as.numeric(.s4dctx_col(d, "kappa_post_mean_cv")),
        kappa_cv_bias = as.numeric(.s4dctx_col(d, "kappa_cv_bias")),
        kappa_rmse = as.numeric(.s4dctx_col(d, "kappa_rmse")),
        log_kappa_rmse = as.numeric(.s4dctx_col(d, "log_kappa_rmse")),
        kappa_mae = as.numeric(.s4dctx_col(d, "kappa_mae")),
        log_kappa_mae = as.numeric(.s4dctx_col(d, "log_kappa_mae")),
        cor_kappa = as.numeric(.s4dctx_col(d, "cor_kappa")),
        cor_log_kappa = as.numeric(.s4dctx_col(d, "cor_log_kappa")),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_aggregate_numeric <- function(df, subset_name, prefix_keep = c("subset", "n")) {
    if (is.null(df) || nrow(df) == 0L) return(data.frame(subset = subset_name, n = 0, stringsAsFactors = FALSE))
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, c("rep_id"))
    out <- data.frame(subset = subset_name, n = nrow(df), stringsAsFactors = FALSE)
    for (cc in numeric_cols) out[[paste0(cc, "_avg")]] <- .s4dctx_safe_mean(df[[cc]])
    for (cc in numeric_cols) out[[paste0(cc, "_sd")]] <- .s4dctx_safe_sd(df[[cc]])
    out
}

.s4dctx_short_time_diagnostics_by_rep <- function(d) {
    data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        TT_short = as.integer(.s4dctx_col(d, "TT_short")),
        TT_reference = as.integer(.s4dctx_col(d, "TT_reference")),
        TT_ratio_to_reference = as.numeric(.s4dctx_col(d, "TT_ratio_to_reference")),
        n_cells = as.numeric(.s4dctx_col(d, "n_cells")),
        mean_count = as.numeric(.s4dctx_col(d, "mean_count")),
        median_count = as.numeric(.s4dctx_col(d, "median_count")),
        zero_prop = as.numeric(.s4dctx_col(d, "zero_prop")),
        total_count = as.numeric(.s4dctx_col(d, "total_count")),
        max_count = as.numeric(.s4dctx_col(d, "max_count")),
        count_cv = as.numeric(.s4dctx_col(d, "count_cv")),
        variance_to_mean = as.numeric(.s4dctx_col(d, "variance_to_mean")),
        q95_count = as.numeric(.s4dctx_col(d, "q95_count")),
        q99_count = as.numeric(.s4dctx_col(d, "q99_count")),
        beta0_star_ident = as.numeric(.s4dctx_col(d, "beta0_star_ident")),
        phi_ident_sd = as.numeric(.s4dctx_col(d, "phi_ident_sd")),
        log_lambda_sd_truth = as.numeric(.s4dctx_col(d, "log_lambda_sd_truth")),
        delta_log_lambda_sd_truth = as.numeric(.s4dctx_col(d, "delta_log_lambda_sd_truth")),
        delta_log_lambda_abs_mean_truth = as.numeric(.s4dctx_col(d, "delta_log_lambda_abs_mean_truth")),
        lambda_raw_median = as.numeric(.s4dctx_col(d, "lambda_raw_median")),
        lambda_raw_max = as.numeric(.s4dctx_col(d, "lambda_raw_max")),
        x2_mode = as.character(.s4dctx_col(d, "x2_mode", NA_character_)),
        mean_x2_within_sd = as.numeric(.s4dctx_col(d, "mean_x2_within_sd")),
        min_x2_within_sd = as.numeric(.s4dctx_col(d, "min_x2_within_sd")),
        mean_lag1_x2 = as.numeric(.s4dctx_col(d, "mean_lag1_x2")),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        stringsAsFactors = FALSE
    )
}

.s4dctx_try_make_comparison <- function(root, s4d_perf, s4d_status_counts, comparison_input_paths = NULL) {
    prev <- NULL
    if (!is.null(comparison_input_paths) && length(comparison_input_paths) > 0L) {
        for (p in comparison_input_paths) {
            if (file.exists(p)) {
                prev <- utils::read.csv(p, stringsAsFactors = FALSE)
                break
            }
        }
    }

    s4d_stable <- s4d_perf[s4d_perf$subset == "stable_only", , drop = FALSE]
    if (nrow(s4d_stable) == 0L) s4d_stable <- s4d_perf[s4d_perf$subset == "stable_plus_soft_warning", , drop = FALSE]
    if (nrow(s4d_stable) == 0L) s4d_stable <- s4d_perf[1L, , drop = FALSE]

    instab_rate <- s4d_status_counts$prop[s4d_status_counts$fit_status == "numerical_instability"]
    if (length(instab_rate) == 0L) instab_rate <- NA_real_

    s4d_row <- data.frame(
        scenario = "S4D short time",
        subset = as.character(s4d_stable$subset[1] %||% "stable_only"),
        n = as.integer(s4d_stable$n[1] %||% NA_integer_),
        mean = as.numeric(s4d_stable$mean_count_avg[1] %||% NA_real_),
        zero = as.numeric(s4d_stable$zero_prop_avg[1] %||% NA_real_),
        beta1_mean = as.numeric(s4d_stable$beta1_mean_avg[1] %||% NA_real_),
        beta2_mean = as.numeric(s4d_stable$beta2_mean_avg[1] %||% NA_real_),
        beta2_sd = as.numeric(s4d_stable$beta2_mean_sd[1] %||% NA_real_),
        r_mean = as.numeric(s4d_stable$r_mean_avg[1] %||% NA_real_),
        phi_rmse = as.numeric(s4d_stable$phi_rmse_avg[1] %||% NA_real_),
        phi_cor = as.numeric(s4d_stable$phi_cor_avg[1] %||% NA_real_),
        loglam_rmse = as.numeric(s4d_stable$log_lambda_rmse_avg[1] %||% NA_real_),
        cor_loglam = as.numeric(s4d_stable$cor_log_lambda_avg[1] %||% NA_real_),
        cor_delta_loglam = as.numeric(s4d_stable$cor_delta_log_lambda_avg[1] %||% NA_real_),
        kappa_log_cor = as.numeric(s4d_stable$cor_log_kappa_avg[1] %||% NA_real_),
        TT_short = as.numeric(s4d_stable$TT_short_avg[1] %||% NA_real_),
        TT_ratio = as.numeric(s4d_stable$TT_ratio_to_reference_avg[1] %||% NA_real_),
        n_cells = as.numeric(s4d_stable$n_cells_avg[1] %||% NA_real_),
        x2_mode = as.character(s4d_stable$x2_mode_unique[1] %||% NA_character_),
        mean_lag1_x2 = as.numeric(s4d_stable$mean_lag1_x2_avg[1] %||% NA_real_),
        mean_x2_within_sd = as.numeric(s4d_stable$mean_x2_within_sd_avg[1] %||% NA_real_),
        instab = as.numeric(instab_rate[1]),
        conclusion = "Short temporal horizon with continuous-time covariates; no numerical failure, regression recovery remains stable, and path-shape correlations are the main weakened component.",
        stringsAsFactors = FALSE
    )

    if (is.null(prev)) return(s4d_row)

    all_cols <- unique(c(names(prev), names(s4d_row)))
    for (cc in setdiff(all_cols, names(prev))) prev[[cc]] <- NA
    for (cc in setdiff(all_cols, names(s4d_row))) s4d_row[[cc]] <- NA
    ans <- rbind(prev[, all_cols, drop = FALSE], s4d_row[, all_cols, drop = FALSE])
    rownames(ans) <- NULL
    ans
}

.s4dctx_make_basic_plots <- function(d, out_fig_dir) {
    .s4dctx_ensure_dir(out_fig_dir)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 is not available; skipping S4D plots.", call. = FALSE)
        return(invisible(NULL))
    }
    ggplot2 <- asNamespace("ggplot2")

    p1 <- ggplot2$ggplot(d, ggplot2$aes(x = rep_id, y = beta1_mean)) +
        ggplot2$geom_hline(yintercept = 0.5, linetype = "dashed") +
        ggplot2$geom_point(ggplot2$aes(shape = fit_status), size = 2) +
        ggplot2$labs(title = "S4D beta1 recovery by replicate", x = "Replicate", y = "Posterior mean of beta1") +
        ggplot2$theme_bw()
    ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_beta1_mean_by_rep.png"), p1, width = 6, height = 4, dpi = 300)

    p2 <- ggplot2$ggplot(d, ggplot2$aes(x = rep_id, y = beta2_mean)) +
        ggplot2$geom_hline(yintercept = -0.4, linetype = "dashed") +
        ggplot2$geom_point(ggplot2$aes(shape = fit_status), size = 2) +
        ggplot2$labs(title = "S4D beta2 recovery by replicate", x = "Replicate", y = "Posterior mean of beta2") +
        ggplot2$theme_bw()
    ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_beta2_mean_by_rep.png"), p2, width = 6, height = 4, dpi = 300)

    p3 <- ggplot2$ggplot(d, ggplot2$aes(x = rep_id, y = log_lambda_rmse)) +
        ggplot2$geom_point(ggplot2$aes(shape = fit_status), size = 2) +
        ggplot2$labs(title = "S4D log-lambda RMSE by replicate", x = "Replicate", y = "Log-lambda RMSE") +
        ggplot2$theme_bw()
    ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_log_lambda_rmse_by_rep.png"), p3, width = 6, height = 4, dpi = 300)

    p4 <- ggplot2$ggplot(d, ggplot2$aes(x = cor_log_lambda, y = cor_delta_log_lambda)) +
        ggplot2$geom_hline(yintercept = 0, linetype = "dotted") +
        ggplot2$geom_vline(xintercept = 0, linetype = "dotted") +
        ggplot2$geom_point(ggplot2$aes(shape = fit_status), size = 2) +
        ggplot2$labs(title = "S4D path-shape recovery metrics", x = "cor(log lambda)", y = "cor(delta log lambda)") +
        ggplot2$theme_bw()
    ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_lambda_correlation_metrics.png"), p4, width = 6, height = 4, dpi = 300)

    p5 <- ggplot2$ggplot(d, ggplot2$aes(x = rep_id)) +
        ggplot2$geom_point(ggplot2$aes(y = kappa_truth_cv, shape = "Truth kappa CV"), size = 2) +
        ggplot2$geom_point(ggplot2$aes(y = kappa_post_mean_cv, shape = "Posterior mean kappa CV"), size = 2) +
        ggplot2$labs(title = "S4D kappa CV recovery", x = "Replicate", y = "Kappa CV", shape = "Quantity") +
        ggplot2$theme_bw()
    ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_kappa_cv_by_rep.png"), p5, width = 6, height = 4, dpi = 300)

    if ("mean_lag1_x2" %in% names(d) && any(is.finite(d$mean_lag1_x2))) {
        p6 <- ggplot2$ggplot(d, ggplot2$aes(x = mean_lag1_x2, y = beta2_mean)) +
            ggplot2$geom_hline(yintercept = -0.4, linetype = "dashed") +
            ggplot2$geom_point(ggplot2$aes(shape = fit_status), size = 2) +
            ggplot2$labs(title = "S4D beta2 recovery versus x2 temporal persistence",
                         x = "Mean lag-1 autocorrelation of x2",
                         y = "Posterior mean of beta2") +
            ggplot2$theme_bw()
        ggplot2$ggsave(file.path(out_fig_dir, "scenario4d_beta2_vs_x2_lag1.png"), p6, width = 6, height = 4, dpi = 300)
    }

    invisible(NULL)
}

summarize_scenario4d_posterior_performance_continuous_x2 <- function(
    root = ".",
    input_path = file.path(root, "output_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2", "summary_S4D_short_time_continuous_x2_fixed_gamma_all_reps.csv"),
    scenario_id = "S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2",
    output_dir = file.path(root, "analysis_s4d_short_time_continuous_x2", scenario_id),
    make_plots = TRUE,
    require_continuous_x2 = TRUE,
    comparison_input_paths = NULL
) {
    d <- .s4dctx_read_required_csv(input_path)
    d <- .s4dctx_add_fit_status(d)
    .s4dctx_check_continuous_x2_input(d, require_continuous_x2 = require_continuous_x2)

    table_dir <- file.path(output_dir, "tables")
    fig_dir <- file.path(output_dir, "figures")
    .s4dctx_ensure_dir(table_dir)
    .s4dctx_ensure_dir(fig_dir)

    fit_status_by_rep <- .s4dctx_make_fit_status_by_rep(d)
    fit_status_counts <- .s4dctx_fit_status_counts(d)

    subsets <- .s4dctx_subset_indices(d)
    perf_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_summarize_subset(d, subsets[[nm]], nm)
    }))

    beta_long <- .s4dctx_beta_long(d)
    beta_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        idx_reps <- d$rep_id[subsets[[nm]]]
        .s4dctx_beta_aggregate(beta_long[beta_long$rep_id %in% idx_reps, , drop = FALSE], nm)
    }))

    lambda_by_rep <- .s4dctx_lambda_by_rep(d)
    lambda_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(lambda_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    phi_by_rep <- .s4dctx_phi_by_rep(d)
    phi_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(phi_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    r_by_rep <- .s4dctx_r_by_rep(d)
    r_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(r_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    kappa_by_rep <- .s4dctx_kappa_by_rep(d)
    kappa_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(kappa_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    short_time_by_rep <- .s4dctx_short_time_diagnostics_by_rep(d)
    short_time_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(short_time_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    x2_covariate_by_rep <- data.frame(
        scenario_id = d$scenario_id,
        rep_id = as.integer(d$rep_id),
        x2_mode = as.character(.s4dctx_col(d, "x2_mode", NA_character_)),
        mean_x2_within_sd = as.numeric(.s4dctx_col(d, "mean_x2_within_sd")),
        min_x2_within_sd = as.numeric(.s4dctx_col(d, "min_x2_within_sd")),
        mean_lag1_x2 = as.numeric(.s4dctx_col(d, "mean_lag1_x2")),
        beta2_mean = as.numeric(.s4dctx_col(d, "beta2_mean")),
        log_lambda_rmse = as.numeric(.s4dctx_col(d, "log_lambda_rmse")),
        cor_log_lambda = as.numeric(.s4dctx_col(d, "cor_log_lambda")),
        fit_status = d$fit_status,
        stringsAsFactors = FALSE
    )

    x2_covariate_by_subset <- .s4dctx_bind_rows_aligned(lapply(names(subsets), function(nm) {
        .s4dctx_aggregate_numeric(x2_covariate_by_rep[subsets[[nm]], , drop = FALSE], nm)
    }))

    comparison <- .s4dctx_try_make_comparison(root, perf_by_subset, fit_status_counts, comparison_input_paths = comparison_input_paths)

    posterior_performance_diagnostics <- data.frame(
        scenario_id = scenario_id,
        n_reps = nrow(d),
        n_stable = sum(d$fit_status == "stable", na.rm = TRUE),
        n_soft_warning = sum(d$fit_status == "soft_warning", na.rm = TRUE),
        n_numerical_instability = sum(d$fit_status == "numerical_instability", na.rm = TRUE),
        stable_rate = mean(d$fit_status == "stable", na.rm = TRUE),
        numerical_instability_rate = mean(d$fit_status == "numerical_instability", na.rm = TRUE),
        TT_short = .s4dctx_safe_mean(.s4dctx_col(d, "TT_short")),
        TT_reference = .s4dctx_safe_mean(.s4dctx_col(d, "TT_reference")),
        TT_ratio_to_reference = .s4dctx_safe_mean(.s4dctx_col(d, "TT_ratio_to_reference")),
        n_cells_avg = .s4dctx_safe_mean(.s4dctx_col(d, "n_cells")),
        x2_mode_unique = paste(unique(stats::na.omit(as.character(.s4dctx_col(d, "x2_mode", NA_character_)))), collapse = ";"),
        mean_x2_within_sd_avg = .s4dctx_safe_mean(.s4dctx_col(d, "mean_x2_within_sd")),
        min_x2_within_sd_min = .s4dctx_safe_min(.s4dctx_col(d, "min_x2_within_sd")),
        mean_lag1_x2_avg = .s4dctx_safe_mean(.s4dctx_col(d, "mean_lag1_x2")),
        mean_count_avg = .s4dctx_safe_mean(.s4dctx_col(d, "mean_count")),
        zero_prop_avg = .s4dctx_safe_mean(.s4dctx_col(d, "zero_prop")),
        beta1_mean_avg = .s4dctx_safe_mean(.s4dctx_col(d, "beta1_mean")),
        beta2_mean_avg = .s4dctx_safe_mean(.s4dctx_col(d, "beta2_mean")),
        r_mean_avg = .s4dctx_safe_mean(.s4dctx_col(d, "r_mean")),
        log_lambda_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(d, "log_lambda_rmse")),
        cor_log_lambda_avg = .s4dctx_safe_mean(.s4dctx_col(d, "cor_log_lambda")),
        cor_delta_log_lambda_avg = .s4dctx_safe_mean(.s4dctx_col(d, "cor_delta_log_lambda")),
        log_kappa_rmse_avg = .s4dctx_safe_mean(.s4dctx_col(d, "log_kappa_rmse")),
        cor_log_kappa_avg = .s4dctx_safe_mean(.s4dctx_col(d, "cor_log_kappa")),
        total_beta_guards = .s4dctx_safe_sum(.s4dctx_col(d, "s4d_beta_guard_count", 0)),
        total_kappa_guards = .s4dctx_safe_sum(.s4dctx_col(d, "s4d_kappa_guard_count", 0)),
        total_lambda_input_guards = .s4dctx_safe_sum(.s4dctx_col(d, "s4d_lambda_input_guard_count", 0)),
        total_lambda_output_guards = .s4dctx_safe_sum(.s4dctx_col(d, "s4d_lambda_output_guard_count", 0)),
        main_interpretation = "T = 30 reduces temporal replication without creating sparse counts or numerical instability; lambda RMSE remains low but path-shape correlations are weaker.",
        stringsAsFactors = FALSE
    )

    ## Write tables.
    .s4dctx_write_csv(d, file.path(table_dir, "scenario4d_summary_with_fit_status.csv"))
    .s4dctx_write_csv(fit_status_by_rep, file.path(table_dir, "scenario4d_fit_status_by_rep.csv"))
    .s4dctx_write_csv(fit_status_counts, file.path(table_dir, "scenario4d_fit_status_counts.csv"))
    .s4dctx_write_csv(perf_by_subset, file.path(table_dir, "scenario4d_performance_by_subset.csv"))
    .s4dctx_write_csv(beta_long, file.path(table_dir, "posterior_beta_summary.csv"))
    .s4dctx_write_csv(beta_by_subset, file.path(table_dir, "scenario4d_beta_recovery_by_subset.csv"))
    .s4dctx_write_csv(lambda_by_rep, file.path(table_dir, "posterior_lambda_path_recovery.csv"))
    .s4dctx_write_csv(lambda_by_subset, file.path(table_dir, "scenario4d_lambda_recovery_by_subset.csv"))
    .s4dctx_write_csv(phi_by_rep, file.path(table_dir, "scenario4d_phi_recovery_by_rep.csv"))
    .s4dctx_write_csv(phi_by_subset, file.path(table_dir, "scenario4d_phi_recovery_by_subset.csv"))
    .s4dctx_write_csv(r_by_rep, file.path(table_dir, "scenario4d_r_recovery_by_rep.csv"))
    .s4dctx_write_csv(r_by_subset, file.path(table_dir, "scenario4d_r_recovery_by_subset.csv"))
    .s4dctx_write_csv(kappa_by_rep, file.path(table_dir, "scenario4d_kappa_recovery_by_rep.csv"))
    .s4dctx_write_csv(kappa_by_subset, file.path(table_dir, "scenario4d_kappa_recovery_by_subset.csv"))
    .s4dctx_write_csv(short_time_by_rep, file.path(table_dir, "scenario4d_short_time_diagnostics_by_rep.csv"))
    .s4dctx_write_csv(short_time_by_subset, file.path(table_dir, "scenario4d_short_time_diagnostics_by_subset.csv"))
    .s4dctx_write_csv(x2_covariate_by_rep, file.path(table_dir, "scenario4d_x2_covariate_diagnostics_by_rep.csv"))
    .s4dctx_write_csv(x2_covariate_by_subset, file.path(table_dir, "scenario4d_x2_covariate_diagnostics_by_subset.csv"))
    .s4dctx_write_csv(posterior_performance_diagnostics, file.path(table_dir, "posterior_performance_diagnostics.csv"))
    .s4dctx_write_csv(comparison, file.path(table_dir, "scenario4d_comparison_to_s3_s4a_s4b_s4c_summary.csv"))

    analysis_object <- list(
        scenario_id = scenario_id,
        input_path = .s4dctx_norm_path(input_path),
        table_dir = .s4dctx_norm_path(table_dir),
        fig_dir = .s4dctx_norm_path(fig_dir),
        summary_with_fit_status = d,
        fit_status_by_rep = fit_status_by_rep,
        fit_status_counts = fit_status_counts,
        performance_by_subset = perf_by_subset,
        beta_long = beta_long,
        beta_recovery_by_subset = beta_by_subset,
        lambda_by_rep = lambda_by_rep,
        lambda_recovery_by_subset = lambda_by_subset,
        phi_by_rep = phi_by_rep,
        phi_recovery_by_subset = phi_by_subset,
        r_by_rep = r_by_rep,
        r_recovery_by_subset = r_by_subset,
        kappa_by_rep = kappa_by_rep,
        kappa_recovery_by_subset = kappa_by_subset,
        short_time_by_rep = short_time_by_rep,
        short_time_by_subset = short_time_by_subset,
        x2_covariate_by_rep = x2_covariate_by_rep,
        x2_covariate_by_subset = x2_covariate_by_subset,
        posterior_performance_diagnostics = posterior_performance_diagnostics,
        comparison = comparison
    )
    saveRDS(analysis_object, file.path(output_dir, "scenario4d_posterior_performance_analysis.rds"))

    if (isTRUE(make_plots)) .s4dctx_make_basic_plots(d, fig_dir)

    cat("\n=== Scenario 4D posterior-performance summary ===\n")
    print(fit_status_counts)
    cat("\n--- Performance by subset ---\n")
    print(perf_by_subset[, intersect(c(
        "subset", "n", "TT_short_avg", "n_cells_avg", "mean_count_avg", "zero_prop_avg",
        "count_cv_avg", "x2_mode_unique", "mean_lag1_x2_avg", "mean_x2_within_sd_avg",
        "beta1_mean_avg", "beta2_mean_avg", "beta2_mean_sd", "r_mean_avg",
        "log_lambda_rmse_avg", "cor_log_lambda_avg", "cor_delta_log_lambda_avg",
        "kappa_post_mean_cv_avg", "log_kappa_rmse_avg", "cor_log_kappa_avg",
        "lambda_output_guard_total"
    ), names(perf_by_subset)), drop = FALSE])
    cat("\n--- Beta recovery by subset ---\n")
    print(beta_by_subset[, intersect(c(
        "subset", "parameter", "n_reps", "mean_truth", "mean_post_mean", "sd_post_mean_across_reps",
        "mean_bias", "rmse_bias", "coverage"
    ), names(beta_by_subset)), drop = FALSE])
    cat("\n--- Lambda recovery by subset ---\n")
    print(lambda_by_subset[, intersect(c(
        "subset", "n", "lambda_rmse_avg", "log_lambda_rmse_avg", "cor_log_lambda_avg", "cor_delta_log_lambda_avg"
    ), names(lambda_by_subset)), drop = FALSE])
    cat("\n--- R recovery by subset ---\n")
    print(r_by_subset[, intersect(c(
        "subset", "n", "r_truth_avg", "r_mean_avg", "r_bias_avg", "r_abs_bias_avg",
        "r_region_coverage_95_avg", "r_region_rmse_avg"
    ), names(r_by_subset)), drop = FALSE])
    cat("\n--- Kappa recovery by subset ---\n")
    print(kappa_by_subset[, intersect(c(
        "subset", "n", "kappa_truth_cv_avg", "kappa_post_mean_cv_avg", "kappa_cv_bias_avg",
        "log_kappa_rmse_avg", "cor_log_kappa_avg"
    ), names(kappa_by_subset)), drop = FALSE])
    cat("\nTables saved to: ", .s4dctx_norm_path(table_dir), "\n", sep = "")
    cat("Figures saved to: ", .s4dctx_norm_path(fig_dir), "\n", sep = "")

    invisible(analysis_object)
}

## Run by default when sourced interactively or from a driver script.
s4d_ctx2_posterior_performance <- summarize_scenario4d_posterior_performance_continuous_x2(root = ".")
