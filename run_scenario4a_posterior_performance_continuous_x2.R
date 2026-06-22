## ============================================================================
## run_scenario4a_posterior_performance_continuous_x2.R
##
## Posterior performance summary for MSSTNB Scenario 4A, continuous-time x2.
##
## This script does NOT run MCMC.  It reads the already-fitted fixed-gamma S4A
## fit objects produced by
##
##   run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100_corrected.R
##
## and writes qmd-ready CSV tables for:
##   * replicate-level fit-status classification;
##   * beta0/beta1/beta2 posterior summaries;
##   * r posterior summaries;
##   * lambda/log-lambda path recovery;
##   * all/stable/stable+soft/unstable subset summaries;
##   * numerical guard summaries.
##
## Official S4A continuous-time x2 design:
##   scenario_id              = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100"
##   x2_mode                  = "continuous_time"
##   x2_ar                    = 0.50
##   x2_innov_sd              = 0.80
##   sparse_beta0_shift       = -4.25
##   beta0_reference_truth    = -1.50
##   beta0_sparse_truth       = -5.75
##   fixed gamma              = 0.80
## ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.s4a_pp_dir_create <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

.s4a_pp_assert_file <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop(label, " not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.s4a_pp_safe_mean <- function(x) {
    x <- as.numeric(x)
    if (!length(x) || all(!is.finite(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.s4a_pp_safe_sd <- function(x) {
    x <- as.numeric(x)
    if (sum(is.finite(x)) <= 1L) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

.s4a_pp_safe_quantile <- function(x, probs) {
    x <- as.numeric(x)
    if (!length(x) || all(!is.finite(x))) {
        return(rep(NA_real_, length(probs)))
    }
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

.s4a_pp_safe_cor <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) <= 2L) return(NA_real_)
    if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
    as.numeric(stats::cor(x[ok], y[ok]))
}

.s4a_pp_rmse <- function(est, truth) {
    est <- as.numeric(est)
    truth <- as.numeric(truth)
    ok <- is.finite(est) & is.finite(truth)
    if (!any(ok)) return(NA_real_)
    sqrt(mean((est[ok] - truth[ok])^2))
}

.s4a_pp_data_file <- function(root,
                              data_dir,
                              data_scenario_id,
                              rep_id) {
    file.path(root, data_dir, data_scenario_id, sprintf("data_rep%02d.rds", as.integer(rep_id)))
}

.s4a_pp_fit_file <- function(root,
                             output_dir,
                             scenario_id,
                             rep_id,
                             fit_subdir = "fits") {
    file.path(root, output_dir, scenario_id, fit_subdir, sprintf("fit_rep%02d.rds", as.integer(rep_id)))
}

.s4a_pp_validate_data <- function(dat,
                                  data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                  expected_TT = 100L,
                                  expected_n1 = 9L,
                                  expected_x2_mode = "continuous_time",
                                  expected_x2_ar = 0.50,
                                  expected_x2_innov_sd = 0.80,
                                  expected_sparse_beta0_shift = -4.25,
                                  expected_beta0_reference_truth = -1.5,
                                  expected_beta0_sparse_truth = -5.75,
                                  strict = TRUE) {
    required <- c(
        "scenario_id", "stress_type", "y_coarse", "x1", "x2", "e",
        "lambda_tilde", "lambda_tilde_ident", "beta0_star_ident",
        "beta_star_ident", "phi_star_ident", "r_star", "gamma_star",
        "x2_mode", "x2_ar", "x2_innov_sd", "sparse_beta0_shift",
        "beta0_reference_truth", "beta0_sparse_truth", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    if (length(missing) > 0L) {
        stop("S4A data is missing required fields: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    if (!identical(dat$scenario_id, data_scenario_id)) {
        stop("Unexpected data scenario_id: ", dat$scenario_id, call. = FALSE)
    }
    if (!identical(as.integer(dim(dat$y_coarse)), c(as.integer(expected_TT), as.integer(expected_n1)))) {
        stop("Unexpected y_coarse dimension: ", paste(dim(dat$y_coarse), collapse = " x "), call. = FALSE)
    }
    if (!identical(dat$x2_mode, expected_x2_mode)) {
        stop("Unexpected x2_mode: ", dat$x2_mode, call. = FALSE)
    }
    if (abs(as.numeric(dat$x2_ar) - expected_x2_ar) > 1e-10) {
        stop("Unexpected x2_ar: ", dat$x2_ar, call. = FALSE)
    }
    if (abs(as.numeric(dat$x2_innov_sd) - expected_x2_innov_sd) > 1e-10) {
        stop("Unexpected x2_innov_sd: ", dat$x2_innov_sd, call. = FALSE)
    }
    if (abs(as.numeric(dat$sparse_beta0_shift) - expected_sparse_beta0_shift) > 1e-10) {
        stop("Unexpected sparse_beta0_shift: ", dat$sparse_beta0_shift, call. = FALSE)
    }
    if (abs(as.numeric(dat$beta0_reference_truth) - expected_beta0_reference_truth) > 1e-10) {
        stop("Unexpected beta0_reference_truth: ", dat$beta0_reference_truth, call. = FALSE)
    }
    if (abs(as.numeric(dat$beta0_sparse_truth) - expected_beta0_sparse_truth) > 1e-10) {
        stop("Unexpected beta0_sparse_truth: ", dat$beta0_sparse_truth, call. = FALSE)
    }

    if (isTRUE(strict)) {
        for (nm in c("x1", "x2", "e", "lambda_tilde", "lambda_tilde_ident")) {
            if (!is.matrix(dat[[nm]]) || !identical(dim(dat[[nm]]), dim(dat$y_coarse))) {
                stop("dat$", nm, " must be a matrix with the same dimension as y_coarse.", call. = FALSE)
            }
        }
        if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0)) {
            stop("y_coarse contains invalid counts.", call. = FALSE)
        }
        if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
            stop("x1/x2 contains non-finite values.", call. = FALSE)
        }
        if (any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
            stop("lambda_tilde_ident must be positive and finite.", call. = FALSE)
        }
        x2_binary_like_prop <- mean(abs(as.numeric(dat$x2)) < 1e-10 |
                                       abs(as.numeric(dat$x2) - 1) < 1e-10)
        if (!is.finite(stats::sd(as.numeric(dat$x2))) || stats::sd(as.numeric(dat$x2)) < 0.10 ||
            x2_binary_like_prop > 0.25) {
            stop("x2 does not look like the required continuous-time covariate.", call. = FALSE)
        }
    }

    invisible(TRUE)
}

.s4a_pp_x2_ar_summary <- function(x2) {
    x2 <- as.matrix(x2)
    if (nrow(x2) < 3L) return(c(mean = NA_real_, median = NA_real_))
    out <- rep(NA_real_, ncol(x2))
    for (jj in seq_len(ncol(x2))) {
        a <- x2[-nrow(x2), jj]
        b <- x2[-1L, jj]
        out[[jj]] <- .s4a_pp_safe_cor(a, b)
    }
    c(mean = .s4a_pp_safe_mean(out), median = stats::median(out, na.rm = TRUE))
}

.s4a_pp_data_summary_one <- function(dat, rep_id) {
    ar <- .s4a_pp_x2_ar_summary(dat$x2)
    data.frame(
        rep = as.integer(rep_id),
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type %||% NA_character_,
        TT = as.integer(nrow(dat$y_coarse)),
        n1 = as.integer(ncol(dat$y_coarse)),
        mean_count = mean(dat$y_coarse),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        zero_prop = mean(dat$y_coarse == 0),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        x2_mode = dat$x2_mode %||% NA_character_,
        x2_ar = as.numeric(dat$x2_ar %||% NA_real_),
        x2_innov_sd = as.numeric(dat$x2_innov_sd %||% NA_real_),
        x2_mean = mean(dat$x2),
        x2_sd = stats::sd(as.numeric(dat$x2)),
        x2_empirical_ar1_mean = ar[["mean"]],
        x2_empirical_ar1_median = ar[["median"]],
        x2_binary_like_prop = mean(abs(as.numeric(dat$x2)) < 1e-10 |
                                       abs(as.numeric(dat$x2) - 1) < 1e-10),
        sparse_beta0_shift = as.numeric(dat$sparse_beta0_shift %||% NA_real_),
        expected_count_multiplier = as.numeric(dat$expected_count_multiplier %||% exp(dat$sparse_beta0_shift)),
        beta0_reference_truth = as.numeric(dat$beta0_reference_truth %||% NA_real_),
        beta0_sparse_truth = as.numeric(dat$beta0_sparse_truth %||% NA_real_),
        beta0_star_ident_truth = as.numeric(dat$beta0_star_ident %||% NA_real_),
        beta1_truth = as.numeric(dat$beta_star_ident[[1L]]),
        beta2_truth = as.numeric(dat$beta_star_ident[[2L]]),
        r_truth_mean = mean(as.numeric(dat$r_star), na.rm = TRUE),
        lambda_truth_min = min(dat$lambda_tilde_ident, na.rm = TRUE),
        lambda_truth_median = stats::median(as.numeric(dat$lambda_tilde_ident), na.rm = TRUE),
        lambda_truth_max = max(dat$lambda_tilde_ident, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

.s4a_pp_param_summary <- function(samples, truth, parameter, rep_id) {
    q <- .s4a_pp_safe_quantile(samples, c(0.025, 0.50, 0.975))
    mn <- .s4a_pp_safe_mean(samples)
    data.frame(
        rep = as.integer(rep_id),
        parameter = parameter,
        truth = as.numeric(truth),
        mean = mn,
        sd = .s4a_pp_safe_sd(samples),
        q025 = q[[1L]],
        q50 = q[[2L]],
        q975 = q[[3L]],
        bias = mn - as.numeric(truth),
        abs_bias = abs(mn - as.numeric(truth)),
        covered_95 = as.integer(is.finite(q[[1L]]) && is.finite(q[[3L]]) &&
                                    q[[1L]] <= as.numeric(truth) && as.numeric(truth) <= q[[3L]]),
        stringsAsFactors = FALSE
    )
}

.s4a_pp_beta_summary_one <- function(fit, dat, rep_id) {
    beta0 <- fit$samples$beta0
    beta <- fit$samples$beta
    if (is.null(beta) || !is.matrix(beta) || ncol(beta) < 2L) {
        stop("fit$samples$beta must be a matrix with at least two columns.", call. = FALSE)
    }
    do.call(rbind, list(
        .s4a_pp_param_summary(beta0, dat$beta0_star_ident, "beta0", rep_id),
        .s4a_pp_param_summary(beta[, 1L], dat$beta_star_ident[[1L]], "beta1", rep_id),
        .s4a_pp_param_summary(beta[, 2L], dat$beta_star_ident[[2L]], "beta2", rep_id)
    ))
}

.s4a_pp_r_summary_one <- function(fit, dat, rep_id) {
    r_samp <- fit$samples$r
    if (is.null(r_samp)) {
        return(.s4a_pp_param_summary(NA_real_, mean(as.numeric(dat$r_star), na.rm = TRUE), "r_mean", rep_id))
    }
    if (is.matrix(r_samp)) {
        r_mean_draw <- rowMeans(r_samp, na.rm = TRUE)
    } else {
        r_mean_draw <- as.numeric(r_samp)
    }
    .s4a_pp_param_summary(r_mean_draw, mean(as.numeric(dat$r_star), na.rm = TRUE), "r_mean", rep_id)
}

.s4a_pp_lambda_summary_one <- function(fit,
                                       dat,
                                       rep_id,
                                       eps = 1e-300,
                                       lower_boundary_tol = 1e-12,
                                       upper_boundary_tol = 1e12) {
    L <- fit$samples$lambda_tilde
    if (is.null(L) || length(dim(L)) != 3L) {
        stop("fit$samples$lambda_tilde must be a 3D array: draws x TT x n1.", call. = FALSE)
    }
    if (!identical(as.integer(dim(L)[2:3]), as.integer(dim(dat$lambda_tilde_ident)))) {
        stop("lambda_tilde sample dimensions do not match truth dimensions.", call. = FALSE)
    }

    truth <- as.matrix(dat$lambda_tilde_ident)
    log_truth <- log(pmax(truth, eps))

    L_safe <- pmax(L, eps)
    logL <- log(L_safe)

    lambda_mean <- apply(L_safe, c(2L, 3L), mean, na.rm = TRUE)
    lambda_q025 <- apply(L_safe, c(2L, 3L), stats::quantile, probs = 0.025, na.rm = TRUE)
    lambda_q975 <- apply(L_safe, c(2L, 3L), stats::quantile, probs = 0.975, na.rm = TRUE)

    log_lambda_mean <- apply(logL, c(2L, 3L), mean, na.rm = TRUE)
    log_lambda_q025 <- apply(logL, c(2L, 3L), stats::quantile, probs = 0.025, na.rm = TRUE)
    log_lambda_q975 <- apply(logL, c(2L, 3L), stats::quantile, probs = 0.975, na.rm = TRUE)

    d_est <- diff(log_lambda_mean)
    d_truth <- diff(log_truth)

    prop_lower_by_draw <- apply(L_safe <= lower_boundary_tol, 1L, mean, na.rm = TRUE)
    prop_upper_by_draw <- apply(L_safe >= upper_boundary_tol, 1L, mean, na.rm = TRUE)
    max_log_lambda_by_draw <- apply(logL, 1L, max, na.rm = TRUE)
    min_log_lambda_by_draw <- apply(logL, 1L, min, na.rm = TRUE)

    data.frame(
        rep = as.integer(rep_id),
        lambda_rmse = .s4a_pp_rmse(lambda_mean, truth),
        log_lambda_rmse = .s4a_pp_rmse(log_lambda_mean, log_truth),
        cor_log_lambda = .s4a_pp_safe_cor(log_lambda_mean, log_truth),
        cor_delta_log_lambda = .s4a_pp_safe_cor(d_est, d_truth),
        lambda_coverage_95 = mean(lambda_q025 <= truth & truth <= lambda_q975, na.rm = TRUE),
        log_lambda_coverage_95 = mean(log_lambda_q025 <= log_truth & log_truth <= log_lambda_q975, na.rm = TRUE),
        lambda_mean_min = min(lambda_mean, na.rm = TRUE),
        lambda_mean_median = stats::median(as.numeric(lambda_mean), na.rm = TRUE),
        lambda_mean_max = max(lambda_mean, na.rm = TRUE),
        log_lambda_mean_min = min(log_lambda_mean, na.rm = TRUE),
        log_lambda_mean_median = stats::median(as.numeric(log_lambda_mean), na.rm = TRUE),
        log_lambda_mean_max = max(log_lambda_mean, na.rm = TRUE),
        lambda_stored_min = min(L_safe, na.rm = TRUE),
        lambda_stored_max = max(L_safe, na.rm = TRUE),
        log_lambda_min_stored = min(logL, na.rm = TRUE),
        log_lambda_max_stored = max(logL, na.rm = TRUE),
        stored_prop_at_lower_mean = mean(prop_lower_by_draw, na.rm = TRUE),
        stored_prop_at_lower_max = max(prop_lower_by_draw, na.rm = TRUE),
        stored_prop_at_upper_mean = mean(prop_upper_by_draw, na.rm = TRUE),
        stored_prop_at_upper_max = max(prop_upper_by_draw, na.rm = TRUE),
        stored_log_lambda_max_mean = mean(max_log_lambda_by_draw, na.rm = TRUE),
        stored_log_lambda_max_q975 = as.numeric(stats::quantile(max_log_lambda_by_draw, 0.975, na.rm = TRUE)),
        stored_log_lambda_min_mean = mean(min_log_lambda_by_draw, na.rm = TRUE),
        stored_log_lambda_min_q025 = as.numeric(stats::quantile(min_log_lambda_by_draw, 0.025, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
}

.s4a_pp_guard_summary_one <- function(fit, rep_id, n_iter, TT, n1) {
    d <- fit$diagnostics %||% list()
    beta_guard <- as.numeric(d$s4a_beta_guard_count %||% 0)
    kappa_guard <- as.numeric(d$s4a_kappa_guard_count %||% 0)
    lambda_input_guard <- as.numeric(d$s4a_lambda_input_guard_count %||% 0)
    lambda_output_guard <- as.numeric(d$s4a_lambda_output_guard_count %||% 0)
    denom <- as.numeric(n_iter) * as.numeric(TT) * as.numeric(n1)
    if (!is.finite(denom) || denom <= 0) denom <- NA_real_
    data.frame(
        rep = as.integer(rep_id),
        beta_guard = beta_guard,
        kappa_guard = kappa_guard,
        lambda_input_guard = lambda_input_guard,
        lambda_output_guard = lambda_output_guard,
        beta_guard_rate_all_iter = beta_guard / as.numeric(n_iter),
        kappa_guard_rate_all_iter = kappa_guard / denom,
        lambda_input_guard_rate_all_iter = lambda_input_guard / denom,
        lambda_output_guard_rate_all_iter = lambda_output_guard / denom,
        stringsAsFactors = FALSE
    )
}

.s4a_pp_classify <- function(beta0_mean,
                             beta_guard,
                             kappa_guard,
                             lambda_output_guard,
                             log_lambda_max_stored,
                             lambda_rmse,
                             log_lambda_rmse,
                             beta0_abs_limit = 20,
                             log_lambda_max_limit = 12,
                             lambda_rmse_limit = Inf,
                             log_lambda_rmse_limit = 5) {
    if (!is.finite(beta0_mean) || !is.finite(log_lambda_max_stored) ||
        isTRUE(kappa_guard > 0) || isTRUE(lambda_output_guard > 0) ||
        isTRUE(abs(beta0_mean) > beta0_abs_limit) ||
        isTRUE(log_lambda_max_stored > log_lambda_max_limit) ||
        isTRUE(is.finite(lambda_rmse_limit) && lambda_rmse > lambda_rmse_limit) ||
        isTRUE(is.finite(log_lambda_rmse_limit) && log_lambda_rmse > log_lambda_rmse_limit)) {
        return("numerical_instability")
    }
    if (isTRUE(beta_guard > 0)) {
        return("soft_warning")
    }
    "stable"
}

.s4a_pp_summarise_one_rep <- function(rep_id,
                                      root,
                                      data_dir,
                                      output_dir,
                                      data_scenario_id,
                                      scenario_id,
                                      fit_subdir,
                                      strict_data,
                                      lower_boundary_tol,
                                      upper_boundary_tol,
                                      beta0_abs_limit,
                                      log_lambda_max_limit,
                                      lambda_rmse_limit,
                                      log_lambda_rmse_limit,
                                      verbose = TRUE) {
    data_file <- .s4a_pp_data_file(root, data_dir, data_scenario_id, rep_id)
    fit_file <- .s4a_pp_fit_file(root, output_dir, scenario_id, rep_id, fit_subdir = fit_subdir)
    .s4a_pp_assert_file(data_file, "S4A data file")
    .s4a_pp_assert_file(fit_file, "S4A fit file")

    if (isTRUE(verbose)) {
        message(sprintf("Reading rep %02d: %s", as.integer(rep_id), fit_file))
    }

    dat <- readRDS(data_file)
    fit <- readRDS(fit_file)
    .s4a_pp_validate_data(dat, data_scenario_id = data_scenario_id, strict = strict_data)

    n_iter <- as.integer((fit$settings %||% list())$n_iter %||%
                             (fit$metadata %||% list())$n_iter %||% NA_integer_)
    n_stored <- as.integer(fit$n_stored %||% length(fit$samples$beta0))

    data_summary <- .s4a_pp_data_summary_one(dat, rep_id)
    beta_summary <- .s4a_pp_beta_summary_one(fit, dat, rep_id)
    r_summary <- .s4a_pp_r_summary_one(fit, dat, rep_id)
    lambda_summary <- .s4a_pp_lambda_summary_one(
        fit = fit,
        dat = dat,
        rep_id = rep_id,
        lower_boundary_tol = lower_boundary_tol,
        upper_boundary_tol = upper_boundary_tol
    )
    guard_summary <- .s4a_pp_guard_summary_one(
        fit = fit,
        rep_id = rep_id,
        n_iter = n_iter,
        TT = data_summary$TT,
        n1 = data_summary$n1
    )

    beta0_row <- beta_summary[beta_summary$parameter == "beta0", , drop = FALSE]
    beta1_row <- beta_summary[beta_summary$parameter == "beta1", , drop = FALSE]
    beta2_row <- beta_summary[beta_summary$parameter == "beta2", , drop = FALSE]

    fit_status_label <- .s4a_pp_classify(
        beta0_mean = beta0_row$mean,
        beta_guard = guard_summary$beta_guard,
        kappa_guard = guard_summary$kappa_guard,
        lambda_output_guard = guard_summary$lambda_output_guard,
        log_lambda_max_stored = lambda_summary$log_lambda_max_stored,
        lambda_rmse = lambda_summary$lambda_rmse,
        log_lambda_rmse = lambda_summary$log_lambda_rmse,
        beta0_abs_limit = beta0_abs_limit,
        log_lambda_max_limit = log_lambda_max_limit,
        lambda_rmse_limit = lambda_rmse_limit,
        log_lambda_rmse_limit = log_lambda_rmse_limit
    )

    rep_summary <- cbind(
        data_summary,
        data.frame(
            fit_file = fit_file,
            data_file = data_file,
            n_iter = n_iter,
            n_stored = n_stored,
            fit_status_label = fit_status_label,
            gamma_fixed_in_fit = TRUE,
            gamma_mean = .s4a_pp_safe_mean(fit$samples$gamma_common %||% fit$samples$gamma),
            gamma_sd = .s4a_pp_safe_sd(fit$samples$gamma_common %||% fit$samples$gamma),
            beta0_mean = beta0_row$mean,
            beta0_sd = beta0_row$sd,
            beta0_q025 = beta0_row$q025,
            beta0_q50 = beta0_row$q50,
            beta0_q975 = beta0_row$q975,
            beta0_bias = beta0_row$bias,
            beta0_covered_95 = beta0_row$covered_95,
            beta1_mean = beta1_row$mean,
            beta1_sd = beta1_row$sd,
            beta1_q025 = beta1_row$q025,
            beta1_q50 = beta1_row$q50,
            beta1_q975 = beta1_row$q975,
            beta1_bias = beta1_row$bias,
            beta1_covered_95 = beta1_row$covered_95,
            beta2_mean = beta2_row$mean,
            beta2_sd = beta2_row$sd,
            beta2_q025 = beta2_row$q025,
            beta2_q50 = beta2_row$q50,
            beta2_q975 = beta2_row$q975,
            beta2_bias = beta2_row$bias,
            beta2_covered_95 = beta2_row$covered_95,
            r_mean = r_summary$mean,
            r_sd = r_summary$sd,
            r_q025 = r_summary$q025,
            r_q50 = r_summary$q50,
            r_q975 = r_summary$q975,
            r_bias = r_summary$bias,
            r_covered_95 = r_summary$covered_95,
            stringsAsFactors = FALSE
        ),
        guard_summary[, setdiff(names(guard_summary), "rep"), drop = FALSE],
        lambda_summary[, setdiff(names(lambda_summary), "rep"), drop = FALSE]
    )

    list(
        rep_summary = rep_summary,
        beta_summary = beta_summary,
        r_summary = r_summary,
        lambda_summary = lambda_summary,
        guard_summary = guard_summary
    )
}

.s4a_pp_status_counts <- function(rep_summary) {
    tab <- as.data.frame(table(rep_summary$fit_status_label), stringsAsFactors = FALSE)
    names(tab) <- c("fit_status_label", "n_reps")
    tab$prop <- tab$n_reps / sum(tab$n_reps)
    order_levels <- c("stable", "soft_warning", "numerical_instability")
    tab$order <- match(tab$fit_status_label, order_levels)
    tab <- tab[order(tab$order, tab$fit_status_label), c("fit_status_label", "n_reps", "prop")]
    rownames(tab) <- NULL
    tab
}

.s4a_pp_subset_index <- function(status, subset) {
    if (subset == "all") return(rep(TRUE, length(status)))
    if (subset == "stable") return(status == "stable")
    if (subset == "stable_plus_soft") return(status %in% c("stable", "soft_warning"))
    if (subset == "numerical_instability") return(status == "numerical_instability")
    rep(FALSE, length(status))
}

.s4a_pp_mean_col <- function(df, col) {
    if (!col %in% names(df) || nrow(df) == 0L) return(NA_real_)
    .s4a_pp_safe_mean(df[[col]])
}

.s4a_pp_sum_col <- function(df, col) {
    if (!col %in% names(df) || nrow(df) == 0L) return(NA_real_)
    sum(as.numeric(df[[col]]), na.rm = TRUE)
}

.s4a_pp_subset_summary <- function(rep_summary, beta_summary, r_summary, lambda_summary, guard_summary) {
    subsets <- c("all", "stable", "stable_plus_soft", "numerical_instability")
    out <- lapply(subsets, function(ss) {
        idx <- .s4a_pp_subset_index(rep_summary$fit_status_label, ss)
        reps_ss <- rep_summary$rep[idx]
        rs <- rep_summary[idx, , drop = FALSE]
        bs <- beta_summary[beta_summary$rep %in% reps_ss, , drop = FALSE]
        ls <- lambda_summary[lambda_summary$rep %in% reps_ss, , drop = FALSE]
        gs <- guard_summary[guard_summary$rep %in% reps_ss, , drop = FALSE]
        rr <- r_summary[r_summary$rep %in% reps_ss, , drop = FALSE]

        bget <- function(param, field) {
            tmp <- bs[bs$parameter == param, , drop = FALSE]
            if (nrow(tmp) == 0L || !field %in% names(tmp)) return(NA_real_)
            .s4a_pp_safe_mean(tmp[[field]])
        }
        bsd_between <- function(param) {
            tmp <- bs[bs$parameter == param, , drop = FALSE]
            if (nrow(tmp) <= 1L) return(NA_real_)
            stats::sd(tmp$mean, na.rm = TRUE)
        }

        data.frame(
            subset = ss,
            n_reps = length(reps_ss),
            reps = paste(reps_ss, collapse = ","),
            mean_count_avg = .s4a_pp_mean_col(rs, "mean_count"),
            median_count_avg = .s4a_pp_mean_col(rs, "median_count"),
            zero_prop_avg = .s4a_pp_mean_col(rs, "zero_prop"),
            total_count_sum = .s4a_pp_sum_col(rs, "total_count"),
            max_count_max = if (nrow(rs)) max(rs$max_count, na.rm = TRUE) else NA_real_,
            beta0_truth_avg = bget("beta0", "truth"),
            beta0_mean_avg = bget("beta0", "mean"),
            beta0_between_rep_sd = bsd_between("beta0"),
            beta0_bias_avg = bget("beta0", "bias"),
            beta0_abs_bias_avg = bget("beta0", "abs_bias"),
            beta0_coverage_rate = bget("beta0", "covered_95"),
            beta1_truth_avg = bget("beta1", "truth"),
            beta1_mean_avg = bget("beta1", "mean"),
            beta1_between_rep_sd = bsd_between("beta1"),
            beta1_bias_avg = bget("beta1", "bias"),
            beta1_abs_bias_avg = bget("beta1", "abs_bias"),
            beta1_coverage_rate = bget("beta1", "covered_95"),
            beta2_truth_avg = bget("beta2", "truth"),
            beta2_mean_avg = bget("beta2", "mean"),
            beta2_between_rep_sd = bsd_between("beta2"),
            beta2_bias_avg = bget("beta2", "bias"),
            beta2_abs_bias_avg = bget("beta2", "abs_bias"),
            beta2_coverage_rate = bget("beta2", "covered_95"),
            r_truth_avg = .s4a_pp_safe_mean(rr$truth),
            r_mean_avg = .s4a_pp_safe_mean(rr$mean),
            r_between_rep_sd = if (nrow(rr) > 1L) stats::sd(rr$mean, na.rm = TRUE) else NA_real_,
            r_bias_avg = .s4a_pp_safe_mean(rr$bias),
            r_abs_bias_avg = .s4a_pp_safe_mean(rr$abs_bias),
            r_coverage_rate = .s4a_pp_safe_mean(rr$covered_95),
            lambda_rmse_mean = .s4a_pp_mean_col(ls, "lambda_rmse"),
            log_lambda_rmse_mean = .s4a_pp_mean_col(ls, "log_lambda_rmse"),
            cor_log_lambda_mean = .s4a_pp_mean_col(ls, "cor_log_lambda"),
            cor_delta_log_lambda_mean = .s4a_pp_mean_col(ls, "cor_delta_log_lambda"),
            lambda_coverage_95_mean = .s4a_pp_mean_col(ls, "lambda_coverage_95"),
            log_lambda_coverage_95_mean = .s4a_pp_mean_col(ls, "log_lambda_coverage_95"),
            beta_guard_total = .s4a_pp_sum_col(gs, "beta_guard"),
            kappa_guard_total = .s4a_pp_sum_col(gs, "kappa_guard"),
            lambda_input_guard_total = .s4a_pp_sum_col(gs, "lambda_input_guard"),
            lambda_output_guard_total = .s4a_pp_sum_col(gs, "lambda_output_guard"),
            lambda_output_guard_rate_mean = .s4a_pp_mean_col(gs, "lambda_output_guard_rate_all_iter"),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

.s4a_pp_lambda_subset_summary <- function(lambda_summary, rep_summary) {
    subsets <- c("all", "stable", "stable_plus_soft", "numerical_instability")
    out <- lapply(subsets, function(ss) {
        idx <- .s4a_pp_subset_index(rep_summary$fit_status_label, ss)
        reps_ss <- rep_summary$rep[idx]
        ls <- lambda_summary[lambda_summary$rep %in% reps_ss, , drop = FALSE]
        data.frame(
            subset = ss,
            n_reps = length(reps_ss),
            reps = paste(reps_ss, collapse = ","),
            lambda_rmse_mean = .s4a_pp_mean_col(ls, "lambda_rmse"),
            lambda_rmse_median = if (nrow(ls)) stats::median(ls$lambda_rmse, na.rm = TRUE) else NA_real_,
            lambda_rmse_max = if (nrow(ls)) max(ls$lambda_rmse, na.rm = TRUE) else NA_real_,
            log_lambda_rmse_mean = .s4a_pp_mean_col(ls, "log_lambda_rmse"),
            log_lambda_rmse_median = if (nrow(ls)) stats::median(ls$log_lambda_rmse, na.rm = TRUE) else NA_real_,
            log_lambda_rmse_max = if (nrow(ls)) max(ls$log_lambda_rmse, na.rm = TRUE) else NA_real_,
            cor_log_lambda_mean = .s4a_pp_mean_col(ls, "cor_log_lambda"),
            cor_delta_log_lambda_mean = .s4a_pp_mean_col(ls, "cor_delta_log_lambda"),
            lambda_coverage_95_mean = .s4a_pp_mean_col(ls, "lambda_coverage_95"),
            log_lambda_coverage_95_mean = .s4a_pp_mean_col(ls, "log_lambda_coverage_95"),
            stored_prop_at_lower_mean = .s4a_pp_mean_col(ls, "stored_prop_at_lower_mean"),
            stored_prop_at_upper_mean = .s4a_pp_mean_col(ls, "stored_prop_at_upper_mean"),
            log_lambda_max_stored_max = if (nrow(ls)) max(ls$log_lambda_max_stored, na.rm = TRUE) else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

.s4a_pp_guard_subset_summary <- function(guard_summary, rep_summary) {
    subsets <- c("all", "stable", "stable_plus_soft", "numerical_instability")
    out <- lapply(subsets, function(ss) {
        idx <- .s4a_pp_subset_index(rep_summary$fit_status_label, ss)
        reps_ss <- rep_summary$rep[idx]
        gs <- guard_summary[guard_summary$rep %in% reps_ss, , drop = FALSE]
        data.frame(
            subset = ss,
            n_reps = length(reps_ss),
            reps = paste(reps_ss, collapse = ","),
            beta_guard_total = .s4a_pp_sum_col(gs, "beta_guard"),
            kappa_guard_total = .s4a_pp_sum_col(gs, "kappa_guard"),
            lambda_input_guard_total = .s4a_pp_sum_col(gs, "lambda_input_guard"),
            lambda_output_guard_total = .s4a_pp_sum_col(gs, "lambda_output_guard"),
            beta_guard_rate_mean = .s4a_pp_mean_col(gs, "beta_guard_rate_all_iter"),
            kappa_guard_rate_mean = .s4a_pp_mean_col(gs, "kappa_guard_rate_all_iter"),
            lambda_input_guard_rate_mean = .s4a_pp_mean_col(gs, "lambda_input_guard_rate_all_iter"),
            lambda_output_guard_rate_mean = .s4a_pp_mean_col(gs, "lambda_output_guard_rate_all_iter"),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

run_scenario4a_posterior_performance_continuous_x2 <- function(
    root = ".",
    reps = 1:10,
    data_dir = "data_s4a_sparse_counts_continuous_x2",
    output_dir = "output_s4a_sparse_counts_continuous_x2",
    data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
    scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
    fit_subdir = "fits",
    table_subdir = "tables",
    strict_data = TRUE,
    lower_boundary_tol = 1e-12,
    upper_boundary_tol = 1e12,
    beta0_abs_limit = 20,
    log_lambda_max_limit = 12,
    lambda_rmse_limit = Inf,
    log_lambda_rmse_limit = 5,
    write_csv = TRUE,
    verbose = TRUE
) {
    root <- normalizePath(root, mustWork = TRUE)
    out_root <- file.path(root, output_dir, scenario_id)
    fit_dir <- file.path(out_root, fit_subdir)
    table_dir <- file.path(out_root, table_subdir)
    .s4a_pp_dir_create(table_dir)

    if (isTRUE(verbose)) {
        cat("\n=== S4A posterior performance: continuous-time x2 ===\n")
        cat("scenario_id: ", scenario_id, "\n", sep = "")
        cat("root:        ", root, "\n", sep = "")
        cat("fit_dir:     ", fit_dir, "\n", sep = "")
        cat("table_dir:   ", table_dir, "\n", sep = "")
        cat("reps:        ", paste(reps, collapse = ","), "\n\n", sep = "")
    }

    pieces <- lapply(reps, function(rr) {
        .s4a_pp_summarise_one_rep(
            rep_id = rr,
            root = root,
            data_dir = data_dir,
            output_dir = output_dir,
            data_scenario_id = data_scenario_id,
            scenario_id = scenario_id,
            fit_subdir = fit_subdir,
            strict_data = strict_data,
            lower_boundary_tol = lower_boundary_tol,
            upper_boundary_tol = upper_boundary_tol,
            beta0_abs_limit = beta0_abs_limit,
            log_lambda_max_limit = log_lambda_max_limit,
            lambda_rmse_limit = lambda_rmse_limit,
            log_lambda_rmse_limit = log_lambda_rmse_limit,
            verbose = verbose
        )
    })

    rep_summary <- do.call(rbind, lapply(pieces, `[[`, "rep_summary"))
    beta_summary <- do.call(rbind, lapply(pieces, `[[`, "beta_summary"))
    r_summary <- do.call(rbind, lapply(pieces, `[[`, "r_summary"))
    lambda_summary <- do.call(rbind, lapply(pieces, `[[`, "lambda_summary"))
    guard_summary <- do.call(rbind, lapply(pieces, `[[`, "guard_summary"))

    status_counts <- .s4a_pp_status_counts(rep_summary)
    subset_summary <- .s4a_pp_subset_summary(
        rep_summary = rep_summary,
        beta_summary = beta_summary,
        r_summary = r_summary,
        lambda_summary = lambda_summary,
        guard_summary = guard_summary
    )
    lambda_subset_summary <- .s4a_pp_lambda_subset_summary(lambda_summary, rep_summary)
    guard_subset_summary <- .s4a_pp_guard_subset_summary(guard_summary, rep_summary)

    manifest <- data.frame(
        table_name = c(
            "s4a_fit_status_counts_continuous_x2",
            "s4a_replicate_level_summary_continuous_x2",
            "posterior_beta_summary_continuous_x2",
            "posterior_r_summary_continuous_x2",
            "lambda_recovery_by_rep_continuous_x2",
            "posterior_performance_by_subset_continuous_x2",
            "lambda_recovery_by_subset_continuous_x2",
            "numerical_guard_summary_continuous_x2",
            "numerical_guard_summary_by_subset_continuous_x2"
        ),
        file = file.path(table_dir, c(
            "s4a_fit_status_counts_continuous_x2.csv",
            "s4a_replicate_level_summary_continuous_x2.csv",
            "posterior_beta_summary_continuous_x2.csv",
            "posterior_r_summary_continuous_x2.csv",
            "lambda_recovery_by_rep_continuous_x2.csv",
            "posterior_performance_by_subset_continuous_x2.csv",
            "lambda_recovery_by_subset_continuous_x2.csv",
            "numerical_guard_summary_continuous_x2.csv",
            "numerical_guard_summary_by_subset_continuous_x2.csv"
        )),
        stringsAsFactors = FALSE
    )

    if (isTRUE(write_csv)) {
        utils::write.csv(status_counts, manifest$file[manifest$table_name == "s4a_fit_status_counts_continuous_x2"], row.names = FALSE)
        utils::write.csv(rep_summary, manifest$file[manifest$table_name == "s4a_replicate_level_summary_continuous_x2"], row.names = FALSE)
        utils::write.csv(beta_summary, manifest$file[manifest$table_name == "posterior_beta_summary_continuous_x2"], row.names = FALSE)
        utils::write.csv(r_summary, manifest$file[manifest$table_name == "posterior_r_summary_continuous_x2"], row.names = FALSE)
        utils::write.csv(lambda_summary, manifest$file[manifest$table_name == "lambda_recovery_by_rep_continuous_x2"], row.names = FALSE)
        utils::write.csv(subset_summary, manifest$file[manifest$table_name == "posterior_performance_by_subset_continuous_x2"], row.names = FALSE)
        utils::write.csv(lambda_subset_summary, manifest$file[manifest$table_name == "lambda_recovery_by_subset_continuous_x2"], row.names = FALSE)
        utils::write.csv(guard_summary, manifest$file[manifest$table_name == "numerical_guard_summary_continuous_x2"], row.names = FALSE)
        utils::write.csv(guard_subset_summary, manifest$file[manifest$table_name == "numerical_guard_summary_by_subset_continuous_x2"], row.names = FALSE)
        utils::write.csv(manifest, file.path(table_dir, "s4a_posterior_performance_manifest_continuous_x2.csv"), row.names = FALSE)
    }

    run_info <- list(
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        reps = reps,
        root = root,
        fit_dir = fit_dir,
        table_dir = table_dir,
        classification_rule = list(
            stable = "no beta/kappa/lambda-output guards and no extreme beta0/log-lambda metrics",
            soft_warning = "beta guard only",
            numerical_instability = paste(
                "kappa_guard > 0 OR lambda_output_guard > 0 OR",
                "abs(beta0_mean) > beta0_abs_limit OR",
                "log_lambda_max_stored > log_lambda_max_limit OR",
                "log_lambda_rmse > log_lambda_rmse_limit"
            ),
            beta0_abs_limit = beta0_abs_limit,
            log_lambda_max_limit = log_lambda_max_limit,
            lambda_rmse_limit = lambda_rmse_limit,
            log_lambda_rmse_limit = log_lambda_rmse_limit
        ),
        status_counts = status_counts,
        rep_summary = rep_summary,
        beta_summary = beta_summary,
        r_summary = r_summary,
        lambda_summary = lambda_summary,
        guard_summary = guard_summary,
        subset_summary = subset_summary,
        lambda_subset_summary = lambda_subset_summary,
        guard_subset_summary = guard_subset_summary,
        manifest = manifest
    )

    if (isTRUE(write_csv)) {
        saveRDS(run_info, file.path(table_dir, "s4a_posterior_performance_continuous_x2.rds"))
    }

    if (isTRUE(verbose)) {
        cat("\nFit-status counts:\n")
        print(status_counts)
        cat("\nSubset posterior-performance summary:\n")
        print(subset_summary)
        cat("\nWrote tables to:\n")
        print(manifest)
        cat("\nDone.\n")
    }

    invisible(run_info)
}

## Convenience alias with shorter name.
run_s4a_posterior_performance_continuous_x2 <- run_scenario4a_posterior_performance_continuous_x2

