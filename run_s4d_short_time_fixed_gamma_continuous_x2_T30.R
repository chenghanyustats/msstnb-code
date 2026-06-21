## ============================================================================
## run_s4d_short_time_fixed_gamma_continuous_x2_T30.R
##
## Scenario 4D short-time-series stress test with continuous-time x2 and gamma fixed.
##
## Purpose
##   Fit the dynamic MSSTNB model used in Scenario 3 to Scenario 4D data, where
##   the data-generating mechanism is unchanged except that the temporal horizon
##   is shortened from T = 100 to T = 30.  This isolates short-time-series stress
##   from S4A sparse counts, S4B low exposure, and S4C strong overdispersion.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   s4d_short_time_series_continuous_x2.R
##   data_s4d_short_time_continuous_x2/S4D_SHORT_TIME_T30_CONTINUOUS_X2/data_repXX.rds
##
## Main outputs
##   output_s4d_short_time_fixed_gamma_continuous_x2/
##     S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2/
##
## Main user-facing function
##   run_s4d_short_time_fixed_gamma_continuous_x2_T30(
##       root = ".",
##       reps = 1,
##       run_profile = "short_test",
##       overwrite_existing = TRUE,
##       verbose = TRUE
##   )
##
## Notes
##   This script only runs fitting and writes one-row summaries/manifests.
##   It hard-checks that x2_mode is continuous_time and that x2 varies within regions.
##   Posterior-performance analysis should be handled later by a separate
##   run_scenario4d_posterior_performance.R script.
## ============================================================================

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir_s4d_fit <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s4d_fit <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true_s4d_fit <- function(x, message) {
    if (!isTRUE(x)) stop(message, call. = FALSE)
    invisible(TRUE)
}

cv_s4d_fit <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x, na.rm = TRUE)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) {
        return(NA_real_)
    }
    stats::sd(x, na.rm = TRUE) / mx
}

safe_ratio_s4d_fit <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}

bind_rows_aligned_s4d_fit <- function(x) {
    if (length(x) == 0L) return(data.frame())
    all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    out <- lapply(x, function(df) {
        missing <- setdiff(all_names, names(df))
        for (nm in missing) df[[nm]] <- NA
        df[, all_names, drop = FALSE]
    })
    do.call(rbind, out)
}

count_stats_s4d_fit <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(
        yy,
        probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
        names = FALSE,
        type = 7
    )
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        nonzero_prop = mean(yy > 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = cv_s4d_fit(yy),
        variance_to_mean = safe_ratio_s4d_fit(vv, mn),
        max_to_mean = safe_ratio_s4d_fit(max(yy), mn),
        q95_count = as.numeric(qs[6]),
        q99_count = as.numeric(qs[7]),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

kappa_stats_s4d_fit <- function(kappa, prefix = "") {
    kk <- as.numeric(kappa)
    qs <- stats::quantile(kk, probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        kappa_mean = mean(kk),
        kappa_sd = stats::sd(kk),
        kappa_cv = cv_s4d_fit(kk),
        kappa_min = min(kk),
        kappa_q01 = as.numeric(qs[1]),
        kappa_q05 = as.numeric(qs[2]),
        kappa_median = as.numeric(qs[3]),
        kappa_q95 = as.numeric(qs[4]),
        kappa_q99 = as.numeric(qs[5]),
        kappa_max = max(kk),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

lambda_stats_s4d_fit <- function(lambda_tilde, prefix = "") {
    lam <- as.matrix(lambda_tilde)
    log_lam <- log(pmax(lam, .Machine$double.xmin))
    dlog <- if (nrow(log_lam) >= 2L) diff(log_lam) else matrix(NA_real_, 0L, ncol(log_lam))
    qs <- stats::quantile(as.numeric(lam), probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    log_qs <- stats::quantile(as.numeric(log_lam), probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                              names = FALSE, type = 7)
    out <- data.frame(
        lambda_raw_min = min(lam),
        lambda_raw_q01 = as.numeric(qs[1]),
        lambda_raw_q05 = as.numeric(qs[2]),
        lambda_raw_median = as.numeric(qs[3]),
        lambda_raw_q95 = as.numeric(qs[4]),
        lambda_raw_q99 = as.numeric(qs[5]),
        lambda_raw_max = max(lam),
        log_lambda_mean = mean(log_lam),
        log_lambda_sd = stats::sd(as.numeric(log_lam)),
        log_lambda_q01 = as.numeric(log_qs[1]),
        log_lambda_q05 = as.numeric(log_qs[2]),
        log_lambda_median = as.numeric(log_qs[3]),
        log_lambda_q95 = as.numeric(log_qs[4]),
        log_lambda_q99 = as.numeric(log_qs[5]),
        delta_log_lambda_mean = if (length(dlog) > 0L) mean(dlog) else NA_real_,
        delta_log_lambda_sd = if (length(dlog) > 0L) stats::sd(as.numeric(dlog)) else NA_real_,
        delta_log_lambda_abs_mean = if (length(dlog) > 0L) mean(abs(dlog)) else NA_real_,
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

summarise_s4d_kappa_recovery <- function(kappa_draws, kappa_truth) {
    if (is.null(kappa_draws) || length(dim(kappa_draws)) != 3L) {
        return(data.frame(
            kappa_rmse = NA_real_,
            log_kappa_rmse = NA_real_,
            cor_kappa = NA_real_,
            cor_log_kappa = NA_real_,
            kappa_truth_cv = NA_real_,
            kappa_post_mean_cv = NA_real_,
            stringsAsFactors = FALSE
        ))
    }
    kappa_mean <- apply(kappa_draws, c(2, 3), mean, na.rm = TRUE)
    kappa_truth <- as.matrix(kappa_truth)
    if (!identical(dim(kappa_mean), dim(kappa_truth))) {
        stop("kappa posterior mean and kappa truth have incompatible dimensions.", call. = FALSE)
    }
    log_km <- log(pmax(kappa_mean, .Machine$double.xmin))
    log_kt <- log(pmax(kappa_truth, .Machine$double.xmin))
    data.frame(
        kappa_rmse = sqrt(mean((kappa_mean - kappa_truth)^2)),
        log_kappa_rmse = sqrt(mean((log_km - log_kt)^2)),
        cor_kappa = suppressWarnings(stats::cor(as.numeric(kappa_mean), as.numeric(kappa_truth))),
        cor_log_kappa = suppressWarnings(stats::cor(as.numeric(log_km), as.numeric(log_kt))),
        kappa_truth_cv = cv_s4d_fit(kappa_truth),
        kappa_post_mean_cv = cv_s4d_fit(kappa_mean),
        stringsAsFactors = FALSE
    )
}

r_recovery_s4d_fit <- function(r_draws, r_truth) {
    if (is.null(r_draws)) {
        return(data.frame(
            r_region_coverage_95 = NA_real_,
            r_region_rmse = NA_real_,
            r_region_mae = NA_real_,
            r_region_mean_sd = NA_real_,
            stringsAsFactors = FALSE
        ))
    }
    r_truth <- as.numeric(r_truth)
    if (length(r_truth) == 1L) r_truth <- rep(r_truth, ncol(r_draws))
    r_mean <- colMeans(r_draws, na.rm = TRUE)
    r_q <- apply(r_draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    r_sd <- apply(r_draws, 2, stats::sd, na.rm = TRUE)
    data.frame(
        r_region_coverage_95 = mean(r_q[1, ] <= r_truth & r_truth <= r_q[3, ]),
        r_region_rmse = sqrt(mean((r_mean - r_truth)^2)),
        r_region_mae = mean(abs(r_mean - r_truth)),
        r_region_mean_sd = mean(r_sd),
        stringsAsFactors = FALSE
    )
}

## ---- source-data validation -------------------------------------------------
s4d_source_data_path <- function(rep_id,
                                 root_arg = ".",
                                 data_dir_arg = "data_s4d_short_time_continuous_x2",
                                 data_scenario_id_arg = "S4D_SHORT_TIME_T30_CONTINUOUS_X2") {
    file.path(
        root_arg,
        data_dir_arg,
        data_scenario_id_arg,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

check_s4d_source_dataset <- function(data_file,
                                     TT_short_expected = 30L,
                                     TT_reference_expected = 100L,
                                     n1_use_arg = 9L,
                                     beta0_ident_abs_limit = 20) {
    assert_file_exists_s4d_fit(data_file, "S4D source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "stress_description", "y_coarse", "e",
        "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "gamma_star",
        "beta0_star", "beta0_star_ident", "beta_star_ident", "phi_star_ident",
        "kappa", "r_star", "r_truth", "TT", "TT_reference", "TT_short",
        "TT_ratio_to_reference", "n1", "n_cells"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s4d_fit(
        length(missing) == 0L,
        paste("Dataset is missing required S4D fields:", paste(missing, collapse = ", "))
    )

    assert_true_s4d_fit(
        identical(dat$stress_type, "short_time_series"),
        paste0("stress_type is not short_time_series. Got ", dat$stress_type, ".")
    )
    assert_true_s4d_fit(
        identical(as.integer(dat$TT), as.integer(TT_short_expected)) &&
            identical(as.integer(dat$TT_short), as.integer(TT_short_expected)),
        paste0("Dataset T is not ", TT_short_expected, ". Got TT=", dat$TT,
               ", TT_short=", dat$TT_short, ".")
    )
    assert_true_s4d_fit(
        identical(as.integer(dat$TT_reference), as.integer(TT_reference_expected)),
        paste0("Dataset TT_reference is not ", TT_reference_expected, ". Got ", dat$TT_reference, ".")
    )
    assert_true_s4d_fit(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_short_expected), as.integer(n1_use_arg))),
        paste0("y_coarse dimension is not TT_short by n1. Got ",
               paste(dim(dat$y_coarse), collapse = " x "), ".")
    )
    assert_true_s4d_fit(
        identical(as.integer(dat$n_cells), as.integer(dat$TT) * as.integer(dat$n1)),
        "n_cells is inconsistent with TT and n1."
    )

    for (nm in c("e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "kappa")) {
        assert_true_s4d_fit(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$e must be positive finite; dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }

    ## Corrected S4D requirement: x2 must be continuous-time and vary within regions.
    x2_mode_now <- as.character(dat$x2_mode %||% NA_character_)
    assert_true_s4d_fit(
        identical(x2_mode_now, "continuous_time"),
        paste0("Corrected S4D requires x2_mode = continuous_time. Got ", x2_mode_now, ".")
    )
    x2_within_sd <- apply(as.matrix(dat$x2), 2, stats::sd)
    assert_true_s4d_fit(
        all(is.finite(x2_within_sd)) && min(x2_within_sd) > 1e-8,
        paste0("Corrected S4D requires time-varying x2; min within-region sd = ",
               min(x2_within_sd, na.rm = TRUE), ".")
    )
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0) ||
        any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("lambda_tilde and lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$kappa)) || any(dat$kappa <= 0)) {
        stop("kappa must be positive and finite.", call. = FALSE)
    }
    assert_true_s4d_fit(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident, ".")
    )

    if (exists("validate_s4dctx_short_time_continuous_x2_data", mode = "function")) {
        validate_s4dctx_short_time_continuous_x2_data(dat)
    }

    y_stats <- count_stats_s4d_fit(dat$y_coarse, prefix = "")
    k_stats <- kappa_stats_s4d_fit(dat$kappa, prefix = "")
    lam_stats <- lambda_stats_s4d_fit(dat$lambda_tilde, prefix = "")
    lam_ident_stats <- lambda_stats_s4d_fit(dat$lambda_tilde_ident, prefix = "lambda_ident_")

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type,
        TT_short = as.integer(dat$TT_short),
        TT_reference = as.integer(dat$TT_reference),
        TT_ratio_to_reference = dat$TT_ratio_to_reference,
        n_cells = as.integer(dat$n_cells),
        r_truth_mean = mean(dat$r_star),
        gamma_truth_mean = mean(dat$gamma_star),
        y_stats = y_stats,
        k_stats = k_stats,
        lam_stats = lam_stats,
        lam_ident_stats = lam_ident_stats,
        beta0_star_ident = dat$beta0_star_ident,
        phi_ident_sd = stats::sd(as.numeric(dat$phi_star_ident)),
        x2_mode = as.character(dat$x2_mode %||% NA_character_),
        mean_x2_within_sd = mean(apply(as.matrix(dat$x2), 2, stats::sd), na.rm = TRUE),
        min_x2_within_sd = min(apply(as.matrix(dat$x2), 2, stats::sd), na.rm = TRUE),
        mean_lag1_x2 = mean(vapply(seq_len(ncol(as.matrix(dat$x2))), function(jj) {
            xx <- as.matrix(dat$x2)[, jj]
            if (length(xx) < 3L || stats::sd(xx[-length(xx)]) < 1e-12 || stats::sd(xx[-1L]) < 1e-12) {
                NA_real_
            } else {
                suppressWarnings(stats::cor(xx[-length(xx)], xx[-1L]))
            }
        }, numeric(1L)), na.rm = TRUE)
    )
}

make_s4d_source_data_manifest <- function(reps,
                                          root_arg = ".",
                                          data_dir_arg = "data_s4d_short_time_continuous_x2",
                                          data_scenario_id_arg = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                          TT_short_expected = 30L,
                                          TT_reference_expected = 100L) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s4d_source_data_path(
            rep_id = rep_id,
            root_arg = root_arg,
            data_dir_arg = data_dir_arg,
            data_scenario_id_arg = data_scenario_id_arg
        )
        chk <- check_s4d_source_dataset(
            data_file,
            TT_short_expected = TT_short_expected,
            TT_reference_expected = TT_reference_expected
        )
        data.frame(
            scenario_id = data_scenario_id_arg,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            TT_short = chk$TT_short,
            TT_reference = chk$TT_reference,
            TT_ratio_to_reference = chk$TT_ratio_to_reference,
            n_cells = chk$n_cells,
            r_truth_mean = chk$r_truth_mean,
            gamma_truth_mean = chk$gamma_truth_mean,
            mean_count = chk$y_stats$mean_count,
            median_count = chk$y_stats$median_count,
            zero_prop = chk$y_stats$zero_prop,
            total_count = chk$y_stats$total_count,
            max_count = chk$y_stats$max_count,
            count_sd = chk$y_stats$count_sd,
            count_cv = chk$y_stats$count_cv,
            variance_to_mean = chk$y_stats$variance_to_mean,
            q95_count = chk$y_stats$q95_count,
            q99_count = chk$y_stats$q99_count,
            kappa_cv = chk$k_stats$kappa_cv,
            x2_mode = chk$x2_mode,
            mean_x2_within_sd = chk$mean_x2_within_sd,
            min_x2_within_sd = chk$min_x2_within_sd,
            mean_lag1_x2 = chk$mean_lag1_x2,
            beta0_star_ident = chk$beta0_star_ident,
            phi_ident_sd = chk$phi_ident_sd,
            lambda_raw_min = chk$lam_stats$lambda_raw_min,
            lambda_raw_median = chk$lam_stats$lambda_raw_median,
            lambda_raw_max = chk$lam_stats$lambda_raw_max,
            log_lambda_sd = chk$lam_stats$log_lambda_sd,
            delta_log_lambda_sd = chk$lam_stats$delta_log_lambda_sd,
            delta_log_lambda_abs_mean = chk$lam_stats$delta_log_lambda_abs_mean,
            stringsAsFactors = FALSE
        )
    })
    bind_rows_aligned_s4d_fit(out)
}

## ---- stabilization overrides ------------------------------------------------
.s4d_fit_env <- new.env(parent = emptyenv())

reset_s4d_numeric_guards <- function() {
    .s4d_fit_env$n_beta_guard <- 0L
    .s4d_fit_env$n_kappa_guard <- 0L
    .s4d_fit_env$n_lambda_input_guard <- 0L
    .s4d_fit_env$n_lambda_output_guard <- 0L
    invisible(TRUE)
}

get_s4d_numeric_guards <- function() {
    list(
        n_beta_guard = .s4d_fit_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4d_fit_env$n_kappa_guard %||% 0L,
        n_lambda_input_guard = .s4d_fit_env$n_lambda_input_guard %||% 0L,
        n_lambda_output_guard = .s4d_fit_env$n_lambda_output_guard %||% 0L
    )
}

install_s4d_stabilization_overrides <- function(fixed_gamma_value = 0.8,
                                                beta0_bounds = c(-30, 10),
                                                beta_bounds = c(-5, 5),
                                                log_xi_bounds = c(-40, 40),
                                                kappa_bounds = c(1e-10, 1e10),
                                                lambda_bounds = c(1e-10, 1e10)) {
    reset_s4d_numeric_guards()

    .s4d_fit_env$beta0_bounds <- beta0_bounds
    .s4d_fit_env$beta_bounds <- beta_bounds
    .s4d_fit_env$log_xi_bounds <- log_xi_bounds
    .s4d_fit_env$kappa_bounds <- kappa_bounds
    .s4d_fit_env$lambda_bounds <- lambda_bounds

    .s4d_fit_env$update_beta_original <- get("update_beta", envir = .GlobalEnv)
    .s4d_fit_env$update_kappa_original <- get("update_kappa", envir = .GlobalEnv)
    .s4d_fit_env$ffbs_lambda_all_original <- get("ffbs_lambda_all", envir = .GlobalEnv)

    assign("compute_s3_xi", function(e, x1, x2, beta0, beta, phi) {
        if (!is.matrix(e) || !is.matrix(x1) || !is.matrix(x2)) {
            stop("e, x1, and x2 must be matrices.", call. = FALSE)
        }
        if (!identical(dim(e), dim(x1)) || !identical(dim(e), dim(x2))) {
            stop("e, x1, and x2 must have the same dimensions.", call. = FALSE)
        }
        if (length(beta) < 2L || length(phi) != ncol(e)) {
            stop("beta must have length at least 2 and phi must match ncol(e).", call. = FALSE)
        }
        if (any(!is.finite(e)) || any(e <= 0)) {
            stop("Exposure matrix e must be positive and finite.", call. = FALSE)
        }
        b <- .s4d_fit_env$log_xi_bounds
        TT_now <- nrow(e)
        n1_now <- ncol(e)
        xi <- matrix(NA_real_, TT_now, n1_now)
        for (j in seq_len(n1_now)) {
            eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
            log_xi_j <- log(e[, j]) + eta_j
            log_xi_j <- pmin(pmax(log_xi_j, b[1]), b[2])
            xi[, j] <- exp(log_xi_j)
        }
        if (any(!is.finite(xi)) || any(xi <= 0)) {
            stop("S4D safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
        }
        xi
    }, envir = .GlobalEnv)

    assign("update_beta", function(beta_current, ...) {
        res <- .s4d_fit_env$update_beta_original(beta_current = beta_current, ...)
        bad <- FALSE
        smp <- res$sample
        b0 <- .s4d_fit_env$beta0_bounds
        bb <- .s4d_fit_env$beta_bounds
        if (length(smp) < 3L || any(!is.finite(smp))) {
            bad <- TRUE
        } else {
            bad <- smp[1] < b0[1] || smp[1] > b0[2] ||
                any(smp[-1] < bb[1] | smp[-1] > bb[2])
        }
        if (isTRUE(bad)) {
            .s4d_fit_env$n_beta_guard <- (.s4d_fit_env$n_beta_guard %||% 0L) + 1L
            res$sample <- beta_current
            res$n_reject <- (res$n_reject %||% 0L) + 1L
            res$s4d_guard_rejected <- TRUE
        } else {
            res$s4d_guard_rejected <- FALSE
        }
        res
    }, envir = .GlobalEnv)

    assign("update_kappa", function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
        y <- as.matrix(y_coarse)
        L <- as.matrix(lambda_tilde)
        X <- as.matrix(xi)
        if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
            stop("S4D safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
        }
        r_vec <- as.numeric(r)
        if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
        if (length(r_vec) != ncol(y)) {
            stop("S4D safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
        }
        R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
        shape <- y + R
        rate <- X * L + R
        guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
        if (any(guard)) {
            .s4d_fit_env$n_kappa_guard <- (.s4d_fit_env$n_kappa_guard %||% 0L) + sum(guard)
        }
        kb <- .s4d_fit_env$kappa_bounds
        shape <- pmin(pmax(shape, 1e-10), 1e10)
        rate <- pmin(pmax(rate, 1e-10), 1e10)
        kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                        nrow = nrow(y), ncol = ncol(y))
        bad_k <- !is.finite(kappa) | kappa <= 0
        if (any(bad_k)) {
            .s4d_fit_env$n_kappa_guard <- (.s4d_fit_env$n_kappa_guard %||% 0L) + sum(bad_k)
            kappa[bad_k] <- 1
        }
        kappa <- pmin(pmax(kappa, kb[1]), kb[2])
        diag <- list(
            mean_kappa = mean(kappa),
            min_kappa = min(kappa),
            max_kappa = max(kappa),
            n_guarded = .s4d_fit_env$n_kappa_guard %||% 0L
        )
        if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
    }, envir = .GlobalEnv)

    assign("ffbs_lambda_all", function(gamma, y_coarse, xi, kappa, a0, b0, return_diag = TRUE) {
        X <- as.matrix(xi)
        K <- as.matrix(kappa)
        bad_in <- !is.finite(X) | X <= 0 | !is.finite(K) | K <= 0
        if (any(bad_in)) {
            .s4d_fit_env$n_lambda_input_guard <- (.s4d_fit_env$n_lambda_input_guard %||% 0L) + sum(bad_in)
        }
        kb <- .s4d_fit_env$kappa_bounds
        lb <- .s4d_fit_env$lambda_bounds
        X <- pmin(pmax(X, 1e-10), 1e10)
        K <- pmin(pmax(K, kb[1]), kb[2])
        out <- .s4d_fit_env$ffbs_lambda_all_original(
            gamma = gamma,
            y_coarse = y_coarse,
            xi = X,
            kappa = K,
            a0 = a0,
            b0 = b0,
            return_diag = return_diag
        )
        L <- out$lambda_tilde
        bad_out <- !is.finite(L) | L <= 0 | L < lb[1] | L > lb[2]
        if (any(bad_out)) {
            .s4d_fit_env$n_lambda_output_guard <- (.s4d_fit_env$n_lambda_output_guard %||% 0L) + sum(bad_out)
            L <- pmin(pmax(L, lb[1]), lb[2])
            out$lambda_tilde <- L
            if (!is.null(out$diag)) {
                out$diag$min_lambda <- min(L)
                out$diag$max_lambda <- max(L)
                out$diag$s4d_lambda_output_guarded <- TRUE
            }
        }
        out
    }, envir = .GlobalEnv)

    assign("update_gamma_common_s3", function(gamma_current,
                                              lambda_tilde,
                                              y_coarse,
                                              a0 = 10,
                                              gamma_prior = c(1, 1),
                                              proposal_sd = 0.15) {
        gamma_new <- rep(as.numeric(fixed_gamma_value), ncol(lambda_tilde))
        list(
            gamma = gamma_new,
            gamma_common = as.numeric(fixed_gamma_value),
            accept = FALSE,
            log_alpha = NA_real_,
            proposal_sd = proposal_sd,
            log_target_current = NA_real_,
            log_target_proposal = NA_real_
        )
    }, envir = .GlobalEnv)

    message(sprintf(
        "Using S4D stabilization guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].",
        beta0_bounds[1], beta0_bounds[2], beta_bounds[1], beta_bounds[2],
        log_xi_bounds[1], log_xi_bounds[2]
    ))
    message(sprintf("Using S4D fixed-gamma override: gamma fixed at %.3f.", fixed_gamma_value))
    invisible(TRUE)
}

## ---- fitting functions ------------------------------------------------------
fit_s4d_short_time_fixed_gamma_one_rep <- function(rep_id,
                                                   scenario_id = "S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2",
                                                   data_scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                                   data_dir = "data_s4d_short_time_continuous_x2",
                                                   output_dir = "output_s4d_short_time_fixed_gamma_continuous_x2",
                                                   root = ".",
                                                   settings_override = list(),
                                                   priors = MCMC_PRIORS,
                                                   spatial = build_s3_spatial(),
                                                   gamma_init = NULL,
                                                   fixed_gamma_value = 0.8,
                                                   gamma_prior = c(1, 1),
                                                   verbose = 1000L,
                                                   save_result = TRUE,
                                                   return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, data_scenario_id, paste0("data_rep", rr, ".rds"))
    chk <- check_s4d_source_dataset(dat_file)
    dat <- chk$dat

    validate_s3_data(dat)
    if (exists("validate_s4dctx_short_time_continuous_x2_data", mode = "function")) {
        validate_s4dctx_short_time_continuous_x2_data(dat)
    }

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) settings[[nm]] <- settings_override[[nm]]
    }

    gamma_start <- gamma_init %||% fixed_gamma_value
    if (length(gamma_start) == 1L) gamma_start <- rep(gamma_start, dat$n1)

    cat(sprintf("=== Scenario 4D short-time continuous-x2 fixed-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Stress    : %s; T = %d (reference %d; ratio %.2f)\n",
                dat$stress_type, dat$TT_short, dat$TT_reference, dat$TT_ratio_to_reference))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f, count CV = %.3f\n",
                chk$y_stats$mean_count, chk$y_stats$zero_prop, chk$y_stats$count_cv))
    cat(sprintf("Lambda    : log-lambda SD = %.3f; delta-log-lambda SD = %.3f\n",
                chk$lam_stats$log_lambda_sd, chk$lam_stats$delta_log_lambda_sd))
    cat(sprintf("Fixed     : gamma = %.3f\n", mean(gamma_start)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    reset_s4d_numeric_guards()

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = verbose
    )

    s4d_guard_counts <- get_s4d_numeric_guards()
    fit$diagnostics$s4d_numeric_guards <- s4d_guard_counts
    fit$diagnostics$s4d_beta_guard_count <- s4d_guard_counts$n_beta_guard
    fit$diagnostics$s4d_kappa_guard_count <- s4d_guard_counts$n_kappa_guard
    fit$diagnostics$s4d_lambda_input_guard_count <- s4d_guard_counts$n_lambda_input_guard
    fit$diagnostics$s4d_lambda_output_guard_count <- s4d_guard_counts$n_lambda_output_guard

    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    fit$diagnostics$gamma_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    fit$metadata$model <- "S4D short-time dynamic NB-ICAR with continuous-time x2 and fixed gamma"
    fit$metadata$fixed_gamma <- TRUE
    fit$metadata$learned_gamma <- FALSE
    fit$metadata$gamma_fixed_value <- fixed_gamma_value
    fit$metadata$updated_blocks <- setdiff(fit$metadata$updated_blocks, "gamma")
    fit$metadata$disabled_blocks <- unique(c(fit$metadata$disabled_blocks, "gamma", "delta", "omega"))

    summary <- summarise_s3_dynamic_learned_gamma_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id)
    )

    gamma_truth_mean <- mean(dat$gamma_star %||% fixed_gamma_value, na.rm = TRUE)
    gamma_fixed_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    summary$gamma_fixed_in_fit <- TRUE
    summary$gamma_learned_in_fit <- FALSE
    summary$gamma_truth_mean <- gamma_truth_mean
    summary$gamma_mean <- gamma_fixed_mean
    summary$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    summary$gamma_q025 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.025, na.rm = TRUE))
    summary$gamma_q50 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.500, na.rm = TRUE))
    summary$gamma_q975 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.975, na.rm = TRUE))
    summary$gamma_bias <- gamma_fixed_mean - gamma_truth_mean
    summary$gamma_covered <- as.integer(summary$gamma_q025 <= gamma_truth_mean && gamma_truth_mean <= summary$gamma_q975)
    summary$gamma_accept_rate <- NA_real_
    summary$gamma_proposal_sd_final <- NA_real_

    ## S4D-specific data context.
    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$TT_short <- as.integer(dat$TT_short)
    summary$TT_reference <- as.integer(dat$TT_reference)
    summary$TT_ratio_to_reference <- dat$TT_ratio_to_reference
    summary$n_cells <- as.integer(dat$n_cells)
    summary$observed_mean_count <- chk$y_stats$mean_count
    summary$observed_median_count <- chk$y_stats$median_count
    summary$observed_zero_prop <- chk$y_stats$zero_prop
    summary$observed_total_count <- chk$y_stats$total_count
    summary$observed_max_count <- chk$y_stats$max_count
    summary$count_sd <- chk$y_stats$count_sd
    summary$count_cv <- chk$y_stats$count_cv
    summary$variance_to_mean <- chk$y_stats$variance_to_mean
    summary$q95_count <- chk$y_stats$q95_count
    summary$q99_count <- chk$y_stats$q99_count
    summary$kappa_truth_cv <- chk$k_stats$kappa_cv
    summary$x2_mode <- chk$x2_mode
    summary$mean_x2_within_sd <- chk$mean_x2_within_sd
    summary$min_x2_within_sd <- chk$min_x2_within_sd
    summary$mean_lag1_x2 <- chk$mean_lag1_x2
    summary$beta0_star_ident <- chk$beta0_star_ident
    summary$phi_ident_sd <- chk$phi_ident_sd
    summary$lambda_raw_min <- chk$lam_stats$lambda_raw_min
    summary$lambda_raw_median <- chk$lam_stats$lambda_raw_median
    summary$lambda_raw_max <- chk$lam_stats$lambda_raw_max
    summary$log_lambda_sd_truth <- chk$lam_stats$log_lambda_sd
    summary$delta_log_lambda_sd_truth <- chk$lam_stats$delta_log_lambda_sd
    summary$delta_log_lambda_abs_mean_truth <- chk$lam_stats$delta_log_lambda_abs_mean

    r_extra <- r_recovery_s4d_fit(fit$samples$r, dat$r_star)
    for (nm in names(r_extra)) summary[[nm]] <- r_extra[[nm]]

    kappa_extra <- summarise_s4d_kappa_recovery(fit$samples$kappa, dat$kappa)
    for (nm in names(kappa_extra)) summary[[nm]] <- kappa_extra[[nm]]

    summary$s4d_beta_guard_count <- fit$diagnostics$s4d_beta_guard_count
    summary$s4d_kappa_guard_count <- fit$diagnostics$s4d_kappa_guard_count
    summary$s4d_lambda_input_guard_count <- fit$diagnostics$s4d_lambda_input_guard_count
    summary$s4d_lambda_output_guard_count <- fit$diagnostics$s4d_lambda_output_guard_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            stress_type = dat$stress_type,
            TT_short = as.integer(dat$TT_short),
            TT_reference = as.integer(dat$TT_reference),
            TT_ratio_to_reference = dat$TT_ratio_to_reference,
            n_cells = as.integer(dat$n_cells),
            model_source = "S3_DYNAMIC_FIXED_GAMMA",
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            fit_file_prefix = "fit_S4D_short_time_continuous_x2_fixed_gamma_rep",
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_fixed_value = fixed_gamma_value,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        ensure_dir_s4d_fit(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4D_short_time_continuous_x2_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4D_short_time_continuous_x2_fixed_gamma_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s4d_short_time_fixed_gamma_batch <- function(reps = 1:10,
                                                 scenario_id = "S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2",
                                                 data_scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                                 data_dir = "data_s4d_short_time_continuous_x2",
                                                 output_dir = "output_s4d_short_time_fixed_gamma_continuous_x2",
                                                 root = ".",
                                                 settings_override = list(),
                                                 gamma_init = NULL,
                                                 fixed_gamma_value = 0.8,
                                                 gamma_prior = c(1, 1),
                                                 verbose = 1000L,
                                                 overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    ensure_dir_s4d_fit(out_dir)

    summaries <- list()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4D_short_time_continuous_x2_fixed_gamma_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s4d_short_time_fixed_gamma_one_rep(
            rep_id = rep_id,
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            root = root,
            settings_override = settings_override,
            gamma_init = gamma_init,
            fixed_gamma_value = fixed_gamma_value,
            gamma_prior = gamma_prior,
            verbose = verbose,
            save_result = TRUE,
            return_result = TRUE
        )
        summaries[[rr]] <- fit$summary
    }

    summary_all <- bind_rows_aligned_s4d_fit(summaries)
    summary_file <- file.path(out_dir, "summary_S4D_short_time_continuous_x2_fixed_gamma_all_reps.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)
    invisible(summary_all)
}

print_s4d_fit_driver_summary <- function(source_manifest, fit_summary, fit_manifest) {
    cat("\n=== Scenario 4D short-time continuous-x2 fixed-gamma fitting summary ===\n")
    cat("Number of reps: ", nrow(source_manifest), "\n", sep = "")
    cat("T short / reference: ", unique(source_manifest$TT_short), " / ",
        unique(source_manifest$TT_reference), "\n", sep = "")
    cat("n_cells avg : ", round(mean(source_manifest$n_cells, na.rm = TRUE), 4), "\n", sep = "")
    cat("Mean count avg: ", round(mean(source_manifest$mean_count, na.rm = TRUE), 4), "\n", sep = "")
    cat("Zero prop avg : ", round(mean(source_manifest$zero_prop, na.rm = TRUE), 4), "\n", sep = "")
    cat("Count CV avg  : ", round(mean(source_manifest$count_cv, na.rm = TRUE), 4), "\n", sep = "")
    cat("log-lambda SD avg: ", round(mean(source_manifest$log_lambda_sd, na.rm = TRUE), 4), "\n", sep = "")
    cat("Fit files all present: ", all(fit_manifest$fit_exists), "\n", sep = "")

    if (!is.null(fit_summary) && nrow(fit_summary) > 0L) {
        selected <- intersect(
            c("beta0_mean", "beta1_mean", "beta2_mean", "r_mean",
              "r_region_coverage_95", "lambda_rmse", "log_lambda_rmse",
              "cor_log_lambda", "cor_delta_log_lambda", "phi_rmse", "phi_cor",
              "s4d_beta_guard_count", "s4d_kappa_guard_count", "s4d_lambda_output_guard_count"),
            names(fit_summary)
        )
        if (length(selected) > 0L) {
            cat("\nSelected fit-summary columns:\n")
            print(fit_summary[, c("scenario_id", "rep_id", selected), drop = FALSE])
        }
    }
    invisible(TRUE)
}

## ---- main driver ------------------------------------------------------------
run_s4d_short_time_fixed_gamma_continuous_x2_T30 <- function(root = ".",
                                               reps = NULL,
                                               run_profile = c("short_test", "formal"),
                                               overwrite_existing = NULL,
                                               verbose = TRUE,
                                               verbose_mcmc = 1000L,
                                               scenario_id = "S4D_SHORT_TIME_FIXED_GAMMA_T30_CONTINUOUS_X2",
                                               data_scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                               data_dir = "data_s4d_short_time_continuous_x2",
                                               output_dir = "output_s4d_short_time_fixed_gamma_continuous_x2",
                                               s3_core_file = "s3_dynamic_learned_gamma.R",
                                               s4d_data_file = "s4d_short_time_series_continuous_x2.R",
                                               fixed_gamma_value = 0.8,
                                               gamma_prior = c(1, 1)) {
    run_profile <- match.arg(run_profile)

    if (is.null(reps)) {
        reps <- if (identical(run_profile, "short_test")) 1:1 else 1:10
    }
    if (is.null(overwrite_existing)) {
        overwrite_existing <- identical(run_profile, "short_test")
    }

    if (identical(run_profile, "short_test")) {
        n_iter <- 6000L
        n_burnin <- 1000L
        n_thin <- 5L
    } else {
        n_iter <- 40000L
        n_burnin <- 20000L
        n_thin <- 5L
    }

    assert_file_exists_s4d_fit(file.path(root, s3_core_file), "Scenario 3 core script")
    assert_file_exists_s4d_fit(file.path(root, s4d_data_file), "Scenario 4D data-generation script")

    source(file.path(root, s4d_data_file))
    source_s4d_short_time_continuous_x2(root = root, verbose = isTRUE(verbose))

    install_s4d_stabilization_overrides(fixed_gamma_value = fixed_gamma_value)

    source_data_manifest <- make_s4d_source_data_manifest(
        reps = reps,
        root_arg = root,
        data_dir_arg = data_dir,
        data_scenario_id_arg = data_scenario_id,
        TT_short_expected = 30L,
        TT_reference_expected = 100L
    )

    if (isTRUE(verbose)) {
        cat("\n=== S4D source-data check ===\n")
        print(source_data_manifest[, c(
            "rep_id", "TT_short", "TT_reference", "n_cells", "mean_count",
            "median_count", "zero_prop", "total_count", "max_count", "count_cv",
            "variance_to_mean", "q99_count", "kappa_cv", "x2_mode",
            "mean_x2_within_sd", "min_x2_within_sd", "mean_lag1_x2",
            "beta0_star_ident", "lambda_raw_median", "lambda_raw_max",
            "log_lambda_sd", "delta_log_lambda_sd"
        ), drop = FALSE])
    }

    fit_summary <- fit_s4d_short_time_fixed_gamma_batch(
        reps = reps,
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root,
        settings_override = list(
            n_iter = n_iter,
            n_burnin = n_burnin,
            n_thin = n_thin
        ),
        fixed_gamma_value = fixed_gamma_value,
        gamma_prior = gamma_prior,
        verbose = verbose_mcmc,
        overwrite_existing = overwrite_existing
    )

    fit_dir_full <- file.path(root, output_dir, scenario_id)
    fit_files <- file.path(
        fit_dir_full,
        sprintf("fit_S4D_short_time_continuous_x2_fixed_gamma_rep%02d.rds", as.integer(reps))
    )
    fit_manifest <- data.frame(
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        rep_id = as.integer(reps),
        data_file = source_data_manifest$data_file,
        fit_file = fit_files,
        fit_exists = file.exists(fit_files),
        stringsAsFactors = FALSE
    )
    assert_true_s4d_fit(all(fit_manifest$fit_exists), "At least one S4D fit file was not created.")

    ensure_dir_s4d_fit(fit_dir_full)
    utils::write.csv(source_data_manifest, file.path(fit_dir_full, "s4d_source_data_manifest.csv"), row.names = FALSE)
    utils::write.csv(fit_manifest, file.path(fit_dir_full, "s4d_fit_manifest.csv"), row.names = FALSE)
    saveRDS(
        list(
            run_profile = run_profile,
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            reps = reps,
            mcmc = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin),
            fixed_gamma_value = fixed_gamma_value,
            gamma_prior = gamma_prior,
            source_data_manifest = source_data_manifest,
            fit_manifest = fit_manifest,
            fit_summary = fit_summary
        ),
        file.path(fit_dir_full, "run_info_S4D_short_time_continuous_x2_fixed_gamma_T30.rds")
    )

    if (isTRUE(verbose)) {
        print_s4d_fit_driver_summary(
            source_manifest = source_data_manifest,
            fit_summary = fit_summary,
            fit_manifest = fit_manifest
        )

        cat("\n=== Main output locations ===\n")
        cat("S4D data: ", file.path(data_dir, data_scenario_id), "\n", sep = "")
        cat("Fits    : ", file.path(output_dir, scenario_id), "\n", sep = "")
        cat("\nScenario 4D short-time continuous-x2 fixed-gamma T30 fitting finished successfully.\n")
    }

    invisible(list(
        source_manifest = source_data_manifest,
        fit_summary = fit_summary,
        fit_manifest = fit_manifest,
        run_profile = run_profile,
        scenario_id = scenario_id,
        output_dir = file.path(output_dir, scenario_id)
    ))
}

message("Loaded run_s4d_short_time_fixed_gamma_continuous_x2_T30.R. Use run_s4d_short_time_fixed_gamma_continuous_x2_T30(...) to fit corrected S4D.")
