## ============================================================================
## run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100.R
##
## Scenario 4A continuous-time x2 sparse-count stress test with gamma fixed.
##
## Purpose
##   Fit the dynamic-lambda MSSTNB model to the revised S4A observation-level
##   sparse-count data generated with the same continuous-time x2 covariate
##   design used for S1--S3 and revised S4D:
##
##       x2_mode     = "continuous_time"
##       x2_ar       = 0.50
##       x2_innov_sd = 0.80
##
##   Gamma is fixed at its truth value, gamma = 0.80.  This script is not a
##   thin wrapper around an external fitting engine.  It directly reuses the
##   Scenario 3 sampler after locally overriding the gamma update to be a no-op
##   and adding the sparse-count numerical guards used in the earlier S4A fixed
##   gamma scripts.
##
## Recommended use from project root
##
##   source("run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100.R")
##
##   out_test <- run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100(
##       root = ".",
##       reps = 1,
##       overwrite_existing = TRUE,
##       verbose = TRUE
##   )
##
##   out_s4a <- run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100(
##       root = ".",
##       reps = 1:10,
##       overwrite_existing = FALSE,
##       verbose = TRUE
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s4a_ct_dir_create <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

.s4a_ct_assert_file <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.s4a_ct_assert_true <- function(x, message) {
    if (!isTRUE(x)) stop(message, call. = FALSE)
    invisible(TRUE)
}

.s4a_ct_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
    invisible(NULL)
}

.s4a_ct_find_first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0L) return(NULL)
    hits[[1L]]
}

.s4a_ct_safe_source <- function(path, verbose = TRUE) {
    .s4a_ct_assert_file(path, "R script")
    .s4a_ct_msg("source: ", path, verbose = verbose)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.s4a_ct_source_s3_core <- function(root = ".",
                                   s3_script = NULL,
                                   verbose = TRUE) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    if (is.null(s3_script)) {
        s3_script <- .s4a_ct_find_first_existing(c(
            file.path(root, "s3_dynamic_learned_gamma.R"),
            file.path(root, "scripts", "s3_dynamic_learned_gamma.R"),
            file.path(root, "R", "s3_dynamic_learned_gamma.R"),
            file.path(root, "R", "scenarios", "s3_dynamic_learned_gamma.R")
        ))
    } else if (!grepl("^(/|[A-Za-z]:)", s3_script)) {
        s3_script <- file.path(root, s3_script)
    }
    .s4a_ct_safe_source(s3_script, verbose = verbose)

    if (!exists("source_s3_dynamic_learned_gamma", mode = "function", inherits = TRUE)) {
        stop("source_s3_dynamic_learned_gamma() was not found after sourcing Scenario 3 script.",
             call. = FALSE)
    }
    get("source_s3_dynamic_learned_gamma", mode = "function", inherits = TRUE)(
        root = root,
        verbose = verbose
    )

    needed <- c(
        "run_s3_dynamic_learned_gamma_mcmc",
        "summarise_s3_dynamic_learned_gamma_fit",
        "validate_s3_data",
        "build_s3_spatial",
        "MCMC_SETTINGS",
        "MCMC_PRIORS",
        "update_beta",
        "update_kappa",
        "ffbs_lambda_all",
        "compute_s3_xi"
    )
    missing <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
    if (length(missing) > 0L) {
        stop("Scenario 3 dependencies missing after source: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    invisible(TRUE)
}

.s4a_ct_data_file <- function(rep_id,
                              root = ".",
                              data_dir = "data_s4a_sparse_counts_continuous_x2",
                              data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100") {
    file.path(
        root,
        data_dir,
        data_scenario_id,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

.s4a_ct_check_source_dataset <- function(data_file,
                                         data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                         TT_use = 100L,
                                         n1_use = 9L,
                                         expected_stress_type = "observation_sparse_counts",
                                         expected_sparse_beta0_shift = -4.25,
                                         expected_beta0_reference_truth = -1.5,
                                         expected_beta0_sparse_truth = -5.75,
                                         expected_x2_mode = "continuous_time",
                                         expected_x2_ar = 0.50,
                                         expected_x2_innov_sd = 0.80,
                                         beta0_ident_abs_limit = 20) {
    .s4a_ct_assert_file(data_file, "S4A continuous-time x2 source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "y_coarse", "e", "x1", "x2",
        "lambda_tilde", "lambda_tilde_ident", "gamma_star", "beta0_star",
        "beta0_star_ident", "beta_star_ident", "phi_star_ident",
        "sparse_beta0_shift", "beta0_reference_truth", "beta0_sparse_truth",
        "expected_count_multiplier", "mean_count", "zero_prop", "TT", "n1",
        "x2_mode", "x2_ar", "x2_innov_sd"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    .s4a_ct_assert_true(
        length(missing) == 0L,
        paste("Dataset is missing required S4A continuous-time fields:",
              paste(missing, collapse = ", "))
    )

    .s4a_ct_assert_true(
        identical(dat$scenario_id, data_scenario_id),
        paste0("scenario_id is not ", data_scenario_id, ". Got ", dat$scenario_id, ".")
    )
    .s4a_ct_assert_true(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0("y_coarse dimension is not TT_use by n1_use. Got ",
               paste(dim(dat$y_coarse), collapse = " x "), ".")
    )
    for (nm in c("e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        .s4a_ct_assert_true(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    .s4a_ct_assert_true(
        identical(dat$stress_type, expected_stress_type),
        paste0("stress_type is not ", expected_stress_type, ". Got ", dat$stress_type, ".")
    )
    .s4a_ct_assert_true(
        abs(dat$sparse_beta0_shift - expected_sparse_beta0_shift) < 1e-12,
        paste0("sparse_beta0_shift is not ", expected_sparse_beta0_shift, ". Got ", dat$sparse_beta0_shift, ".")
    )
    .s4a_ct_assert_true(
        abs(dat$beta0_reference_truth - expected_beta0_reference_truth) < 1e-12,
        paste0("beta0_reference_truth is not ", expected_beta0_reference_truth, ". Got ",
               dat$beta0_reference_truth, ".")
    )
    .s4a_ct_assert_true(
        abs(dat$beta0_sparse_truth - expected_beta0_sparse_truth) < 1e-12,
        paste0("beta0_sparse_truth is not ", expected_beta0_sparse_truth, ". Got ",
               dat$beta0_sparse_truth, ".")
    )
    .s4a_ct_assert_true(
        abs(dat$expected_count_multiplier - exp(expected_sparse_beta0_shift)) < 1e-10,
        "expected_count_multiplier is inconsistent with sparse_beta0_shift."
    )

    .s4a_ct_assert_true(
        identical(dat$x2_mode, expected_x2_mode),
        paste0("x2_mode is not ", expected_x2_mode, ". Got ", dat$x2_mode, ".")
    )
    .s4a_ct_assert_true(
        abs(as.numeric(dat$x2_ar) - expected_x2_ar) < 1e-12,
        paste0("x2_ar is not ", expected_x2_ar, ". Got ", dat$x2_ar, ".")
    )
    .s4a_ct_assert_true(
        abs(as.numeric(dat$x2_innov_sd) - expected_x2_innov_sd) < 1e-12,
        paste0("x2_innov_sd is not ", expected_x2_innov_sd, ". Got ", dat$x2_innov_sd, ".")
    )

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0)) {
        stop("dat$e must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0)) {
        stop("dat$lambda_tilde must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("dat$lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    if (!is.finite(dat$beta0_star_ident) || abs(dat$beta0_star_ident) > beta0_ident_abs_limit) {
        stop("beta0_star_ident is outside the expected identified-scale guard range.", call. = FALSE)
    }

    x2_binary_like_prop <- mean(dat$x2 %in% c(0, 1))
    if (!is.finite(stats::sd(as.numeric(dat$x2))) || stats::sd(as.numeric(dat$x2)) < 0.10 ||
        x2_binary_like_prop > 0.25) {
        stop("x2 does not look like the required continuous-time covariate.", call. = FALSE)
    }

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat$lambda_tilde_ident, finite = TRUE)

    list(
        dat = dat,
        scenario_id = dat$scenario_id,
        stress_type = dat$stress_type,
        x2_mode = dat$x2_mode,
        x2_ar = dat$x2_ar,
        x2_innov_sd = dat$x2_innov_sd,
        x2_sd = stats::sd(as.numeric(dat$x2)),
        x2_binary_like_prop = x2_binary_like_prop,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_ident_min = lambda_ident_range[[1L]],
        lambda_ident_max = lambda_ident_range[[2L]]
    )
}

.s4a_ct_make_source_data_manifest <- function(reps,
                                              root = ".",
                                              data_dir = "data_s4a_sparse_counts_continuous_x2",
                                              data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100") {
    out <- lapply(reps, function(rep_id) {
        data_file <- .s4a_ct_data_file(rep_id, root = root, data_dir = data_dir,
                                       data_scenario_id = data_scenario_id)
        chk <- .s4a_ct_check_source_dataset(data_file, data_scenario_id = data_scenario_id)
        data.frame(
            scenario_id = data_scenario_id,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            x2_mode = chk$x2_mode,
            x2_ar = chk$x2_ar,
            x2_innov_sd = chk$x2_innov_sd,
            x2_sd = chk$x2_sd,
            x2_binary_like_prop = chk$x2_binary_like_prop,
            mean_count = chk$mean_count,
            median_count = chk$median_count,
            zero_prop = chk$zero_prop,
            total_count = chk$total_count,
            max_count = chk$max_count,
            beta0_star_ident = chk$beta0_star_ident,
            lambda_raw_min = chk$lambda_raw_min,
            lambda_raw_max = chk$lambda_raw_max,
            lambda_ident_min = chk$lambda_ident_min,
            lambda_ident_max = chk$lambda_ident_max,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

## ---- sparse-count numerical guards -----------------------------------------
.s4a_ct_guard_env <- new.env(parent = emptyenv())

reset_s4a_continuous_x2_numeric_guards <- function() {
    .s4a_ct_guard_env$n_beta_guard <- 0L
    .s4a_ct_guard_env$n_kappa_guard <- 0L
    .s4a_ct_guard_env$n_lambda_input_guard <- 0L
    .s4a_ct_guard_env$n_lambda_output_guard <- 0L
    invisible(TRUE)
}

get_s4a_continuous_x2_numeric_guards <- function() {
    list(
        n_beta_guard = .s4a_ct_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4a_ct_guard_env$n_kappa_guard %||% 0L,
        n_lambda_input_guard = .s4a_ct_guard_env$n_lambda_input_guard %||% 0L,
        n_lambda_output_guard = .s4a_ct_guard_env$n_lambda_output_guard %||% 0L
    )
}

.s4a_ct_install_fixed_gamma_and_guards <- function(fixed_gamma_value = 0.8,
                                                   verbose = TRUE) {
    reset_s4a_continuous_x2_numeric_guards()

    s4a_log_xi_lower <- -40
    s4a_log_xi_upper <-  40
    s4a_beta0_bounds <- c(-30, 10)
    s4a_beta_bounds <- c(-5, 5)
    s4a_kappa_bounds <- c(1e-10, 1e10)
    s4a_lambda_bounds <- c(1e-10, 1e10)

    compute_s3_xi_safe <- function(e, x1, x2, beta0, beta, phi) {
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

        TT_now <- nrow(e)
        n1_now <- ncol(e)
        xi <- matrix(NA_real_, TT_now, n1_now)
        for (j in seq_len(n1_now)) {
            eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
            log_xi_j <- log(e[, j]) + eta_j
            log_xi_j <- pmin(pmax(log_xi_j, s4a_log_xi_lower), s4a_log_xi_upper)
            xi[, j] <- exp(log_xi_j)
        }
        if (any(!is.finite(xi)) || any(xi <= 0)) {
            stop("S4A safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
        }
        xi
    }

    update_beta_original <- get("update_beta", mode = "function", inherits = TRUE)
    update_beta_safe <- function(beta_current, ...) {
        res <- update_beta_original(beta_current = beta_current, ...)
        bad <- FALSE
        smp <- res$sample
        if (length(smp) < 3L || any(!is.finite(smp))) {
            bad <- TRUE
        } else {
            bad <- smp[1] < s4a_beta0_bounds[1] || smp[1] > s4a_beta0_bounds[2] ||
                any(smp[-1] < s4a_beta_bounds[1] | smp[-1] > s4a_beta_bounds[2])
        }
        if (isTRUE(bad)) {
            .s4a_ct_guard_env$n_beta_guard <- (.s4a_ct_guard_env$n_beta_guard %||% 0L) + 1L
            res$sample <- beta_current
            res$n_reject <- (res$n_reject %||% 0L) + 1L
            res$s4a_guard_rejected <- TRUE
        } else {
            res$s4a_guard_rejected <- FALSE
        }
        res
    }

    update_kappa_safe <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
        y <- as.matrix(y_coarse)
        L <- as.matrix(lambda_tilde)
        X <- as.matrix(xi)
        if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
            stop("S4A safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.",
                 call. = FALSE)
        }
        r_vec <- as.numeric(r)
        if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
        if (length(r_vec) != ncol(y)) {
            stop("S4A safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
        }
        R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
        shape <- y + R
        rate <- X * L + R
        guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
        if (any(guard)) {
            .s4a_ct_guard_env$n_kappa_guard <- (.s4a_ct_guard_env$n_kappa_guard %||% 0L) + sum(guard)
        }
        shape <- pmin(pmax(shape, 1e-10), 1e10)
        rate <- pmin(pmax(rate, 1e-10), 1e10)
        kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                        nrow = nrow(y), ncol = ncol(y))
        bad_k <- !is.finite(kappa) | kappa <= 0
        if (any(bad_k)) {
            .s4a_ct_guard_env$n_kappa_guard <- (.s4a_ct_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
            kappa[bad_k] <- 1
        }
        kappa <- pmin(pmax(kappa, s4a_kappa_bounds[1]), s4a_kappa_bounds[2])
        diag <- list(
            mean_kappa = mean(kappa),
            min_kappa = min(kappa),
            max_kappa = max(kappa),
            n_guarded = .s4a_ct_guard_env$n_kappa_guard %||% 0L
        )
        if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
    }

    ffbs_lambda_original <- get("ffbs_lambda_all", mode = "function", inherits = TRUE)
    ffbs_lambda_safe <- function(gamma, y_coarse, xi, kappa, a0, b0, return_diag = TRUE) {
        X <- as.matrix(xi)
        K <- as.matrix(kappa)
        bad_in <- !is.finite(X) | X <= 0 | !is.finite(K) | K <= 0
        if (any(bad_in)) {
            .s4a_ct_guard_env$n_lambda_input_guard <- (.s4a_ct_guard_env$n_lambda_input_guard %||% 0L) + sum(bad_in)
        }
        X <- pmin(pmax(X, 1e-10), 1e10)
        K <- pmin(pmax(K, s4a_kappa_bounds[1]), s4a_kappa_bounds[2])
        out <- ffbs_lambda_original(
            gamma = gamma,
            y_coarse = y_coarse,
            xi = X,
            kappa = K,
            a0 = a0,
            b0 = b0,
            return_diag = return_diag
        )
        L <- out$lambda_tilde
        bad_out <- !is.finite(L) | L <= 0 | L < s4a_lambda_bounds[1] | L > s4a_lambda_bounds[2]
        if (any(bad_out)) {
            .s4a_ct_guard_env$n_lambda_output_guard <- (.s4a_ct_guard_env$n_lambda_output_guard %||% 0L) + sum(bad_out)
            L <- pmin(pmax(L, s4a_lambda_bounds[1]), s4a_lambda_bounds[2])
            out$lambda_tilde <- L
            if (!is.null(out$diag)) {
                out$diag$min_lambda <- min(L)
                out$diag$max_lambda <- max(L)
                out$diag$s4a_lambda_output_guarded <- TRUE
            }
        }
        out
    }

    update_gamma_fixed <- function(gamma_current,
                                   lambda_tilde,
                                   y_coarse,
                                   a0 = 10,
                                   gamma_prior = c(1, 1),
                                   proposal_sd = 0.15) {
        gamma_common_current <- as.numeric(fixed_gamma_value)
        gamma_common_current <- min(max(gamma_common_current, 1e-12), 1 - 1e-12)
        gamma_new <- rep(gamma_common_current, ncol(lambda_tilde))
        list(
            gamma = gamma_new,
            gamma_common = gamma_common_current,
            accept = FALSE,
            log_alpha = NA_real_,
            proposal_sd = proposal_sd,
            log_target_current = NA_real_,
            log_target_proposal = NA_real_
        )
    }

    assign("compute_s3_xi", compute_s3_xi_safe, envir = .GlobalEnv)
    assign("update_beta", update_beta_safe, envir = .GlobalEnv)
    assign("update_kappa", update_kappa_safe, envir = .GlobalEnv)
    assign("ffbs_lambda_all", ffbs_lambda_safe, envir = .GlobalEnv)
    assign("update_gamma_common_s3", update_gamma_fixed, envir = .GlobalEnv)

    .s4a_ct_msg(sprintf(
        "Using S4A continuous-time x2 guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].",
        s4a_beta0_bounds[1], s4a_beta0_bounds[2],
        s4a_beta_bounds[1], s4a_beta_bounds[2],
        s4a_log_xi_lower, s4a_log_xi_upper
    ), verbose = verbose)
    .s4a_ct_msg(sprintf("Using S4A fixed-gamma override: gamma fixed at %.3f.",
                        fixed_gamma_value), verbose = verbose)

    invisible(TRUE)
}

fit_s4a_sparse_counts_fixed_gamma_continuous_x2_one_rep <- function(rep_id,
                                                                    scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_FIXED_GAMMA_T100",
                                                                    data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                                                    data_dir = "data_s4a_sparse_counts_continuous_x2",
                                                                    output_dir = "output_s4a_sparse_counts_continuous_x2",
                                                                    root = ".",
                                                                    settings_override = list(),
                                                                    priors = NULL,
                                                                    spatial = NULL,
                                                                    gamma_init = NULL,
                                                                    fixed_gamma_value = 0.8,
                                                                    gamma_prior = c(1, 1),
                                                                    mcmc_verbose = 1000L,
                                                                    save_result = TRUE,
                                                                    return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    data_file <- .s4a_ct_data_file(rep_id, root = root, data_dir = data_dir,
                                   data_scenario_id = data_scenario_id)
    chk <- .s4a_ct_check_source_dataset(data_file, data_scenario_id = data_scenario_id)
    dat <- chk$dat

    validate_s3_data(dat)

    settings <- get("MCMC_SETTINGS", inherits = TRUE)
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) settings[[nm]] <- settings_override[[nm]]
    }
    if (is.null(priors)) priors <- get("MCMC_PRIORS", inherits = TRUE)
    if (is.null(spatial)) spatial <- build_s3_spatial()

    gamma_start <- gamma_init %||% fixed_gamma_value
    if (length(gamma_start) == 1L) gamma_start <- rep(gamma_start, dat$n1)

    cat(sprintf("=== Scenario 4A continuous-time x2 fixed-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", data_file))
    cat(sprintf("Stress    : %s; sparse_beta0_shift = %.2f\n", dat$stress_type, dat$sparse_beta0_shift))
    cat(sprintf("x2 design : %s; AR = %.2f; innov_sd = %.2f\n",
                dat$x2_mode, dat$x2_ar, dat$x2_innov_sd))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Fixed     : gamma = %.3f\n", fixed_gamma_value))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    reset_s4a_continuous_x2_numeric_guards()

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = mcmc_verbose
    )

    guard_counts <- get_s4a_continuous_x2_numeric_guards()
    fit$diagnostics$s4a_numeric_guards <- guard_counts
    fit$diagnostics$s4a_beta_guard_count <- guard_counts$n_beta_guard
    fit$diagnostics$s4a_kappa_guard_count <- guard_counts$n_kappa_guard
    fit$diagnostics$s4a_lambda_input_guard_count <- guard_counts$n_lambda_input_guard
    fit$diagnostics$s4a_lambda_output_guard_count <- guard_counts$n_lambda_output_guard

    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    fit$diagnostics$gamma_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)

    fit$metadata$model <- "S4A sparse-count continuous-time x2 dynamic NB-ICAR with fixed gamma"
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

    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$x2_mode <- dat$x2_mode %||% NA_character_
    summary$x2_ar <- dat$x2_ar %||% NA_real_
    summary$x2_innov_sd <- dat$x2_innov_sd %||% NA_real_
    summary$sparse_beta0_shift <- dat$sparse_beta0_shift %||% NA_real_
    summary$expected_count_multiplier <- dat$expected_count_multiplier %||% NA_real_
    summary$reference_mean_count <- dat$reference_mean_count %||%
        (dat$reference_count_summary$mean_count %||% NA_real_)
    summary$reference_zero_prop <- dat$reference_zero_prop %||%
        (dat$reference_count_summary$zero_prop %||% NA_real_)
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count
    summary$s4a_beta_guard_count <- fit$diagnostics$s4a_beta_guard_count
    summary$s4a_kappa_guard_count <- fit$diagnostics$s4a_kappa_guard_count
    summary$s4a_lambda_input_guard_count <- fit$diagnostics$s4a_lambda_input_guard_count
    summary$s4a_lambda_output_guard_count <- fit$diagnostics$s4a_lambda_output_guard_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            stress_type = dat$stress_type,
            x2_mode = dat$x2_mode,
            x2_ar = dat$x2_ar,
            x2_innov_sd = dat$x2_innov_sd,
            sparse_beta0_shift = dat$sparse_beta0_shift,
            beta0_reference_truth = dat$beta0_reference_truth,
            beta0_sparse_truth = dat$beta0_sparse_truth,
            expected_count_multiplier = dat$expected_count_multiplier,
            model_source = "S3_DYNAMIC_FIXED_GAMMA",
            rep_id = as.integer(rep_id),
            data_file = data_file,
            output_dir = output_dir,
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_fixed_value = fixed_gamma_value,
            run_time = Sys.time()
        )
    )

    fit_file <- NA_character_
    csv_file <- NA_character_
    if (isTRUE(save_result)) {
        fit_dir <- file.path(root, output_dir, data_scenario_id, "fits")
        table_dir <- file.path(root, output_dir, data_scenario_id, "tables")
        .s4a_ct_dir_create(fit_dir)
        .s4a_ct_dir_create(table_dir)
        fit_file <- file.path(fit_dir, paste0("fit_rep", rr, ".rds"))
        csv_file <- file.path(table_dir, paste0("summary_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) {
        fit$metadata$saved_fit_file <- fit_file
        fit$metadata$saved_summary_file <- csv_file
        return(fit)
    }
    invisible(NULL)
}

fit_s4a_sparse_counts_fixed_gamma_continuous_x2_batch <- function(reps = 1:10,
                                                                  scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_FIXED_GAMMA_T100",
                                                                  data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                                                  data_dir = "data_s4a_sparse_counts_continuous_x2",
                                                                  output_dir = "output_s4a_sparse_counts_continuous_x2",
                                                                  root = ".",
                                                                  settings_override = list(),
                                                                  gamma_init = NULL,
                                                                  fixed_gamma_value = 0.8,
                                                                  gamma_prior = c(1, 1),
                                                                  mcmc_verbose = 1000L,
                                                                  overwrite_existing = FALSE) {
    fit_dir <- file.path(root, output_dir, data_scenario_id, "fits")
    table_dir <- file.path(root, output_dir, data_scenario_id, "tables")
    .s4a_ct_dir_create(fit_dir)
    .s4a_ct_dir_create(table_dir)

    summaries <- list()
    fit_rows <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(fit_dir, paste0("fit_rep", rr, ".rds"))
        data_file <- .s4a_ct_data_file(rep_id, root = root, data_dir = data_dir,
                                       data_scenario_id = data_scenario_id)

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            status <- "skipped_existing"
        } else {
            fit <- fit_s4a_sparse_counts_fixed_gamma_continuous_x2_one_rep(
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
                mcmc_verbose = mcmc_verbose,
                save_result = TRUE,
                return_result = TRUE
            )
            summaries[[rr]] <- fit$summary
            status <- "fit_completed"
        }

        guards <- fit$diagnostics$s4a_numeric_guards %||% list()
        fit_rows[[rr]] <- data.frame(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            fit_file = fit_file,
            fit_exists = file.exists(fit_file),
            fit_status = status,
            n_beta_guard = guards$n_beta_guard %||% NA_integer_,
            n_kappa_guard = guards$n_kappa_guard %||% NA_integer_,
            n_lambda_input_guard = guards$n_lambda_input_guard %||% NA_integer_,
            n_lambda_output_guard = guards$n_lambda_output_guard %||% NA_integer_,
            stringsAsFactors = FALSE
        )
    }

    summary_all <- do.call(rbind, summaries)
    manifest <- do.call(rbind, fit_rows)

    summary_file <- file.path(table_dir, "summary_S4A_sparse_counts_fixed_gamma_continuous_x2_all_reps.csv")
    manifest_file <- file.path(table_dir, "manifest_S4A_sparse_counts_fixed_gamma_continuous_x2_fits.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    utils::write.csv(manifest, manifest_file, row.names = FALSE)

    message("Saved combined summary: ", summary_file)
    message("Saved fit manifest    : ", manifest_file)

    list(
        fit_summary = summary_all,
        fit_manifest = manifest,
        summary_file = summary_file,
        manifest_file = manifest_file
    )
}

.s4a_ct_make_status_counts <- function(fit_manifest) {
    tab <- as.data.frame(table(fit_manifest$fit_status), stringsAsFactors = FALSE)
    names(tab) <- c("fit_status", "n_reps")
    tab$prop <- tab$n_reps / sum(tab$n_reps)
    tab
}

.s4a_ct_make_compact_summary <- function(source_manifest, fit_summary, fit_manifest) {
    out <- data.frame(
        n_reps = length(unique(fit_manifest$rep_id)),
        n_fit_files = sum(fit_manifest$fit_exists),
        mean_count_avg = mean(source_manifest$mean_count, na.rm = TRUE),
        zero_prop_avg = mean(source_manifest$zero_prop, na.rm = TRUE),
        x2_mode_unique = paste(sort(unique(source_manifest$x2_mode)), collapse = ","),
        x2_ar_unique = paste(sort(unique(source_manifest$x2_ar)), collapse = ","),
        x2_innov_sd_unique = paste(sort(unique(source_manifest$x2_innov_sd)), collapse = ","),
        beta_guard_total = sum(fit_manifest$n_beta_guard, na.rm = TRUE),
        kappa_guard_total = sum(fit_manifest$n_kappa_guard, na.rm = TRUE),
        lambda_input_guard_total = sum(fit_manifest$n_lambda_input_guard, na.rm = TRUE),
        lambda_output_guard_total = sum(fit_manifest$n_lambda_output_guard, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    for (nm in c("beta0_mean", "beta1_mean", "beta2_mean", "r_mean",
                 "lambda_rmse_mean", "lambda_log_rmse_mean", "lambda_cor_mean",
                 "lambda_delta_cor_mean")) {
        if (nm %in% names(fit_summary)) {
            out[[paste0(nm, "_avg")]] <- mean(fit_summary[[nm]], na.rm = TRUE)
        }
    }
    out
}

run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100 <- function(root = ".",
                                                                reps = 1:10,
                                                                s3_script = NULL,
                                                                data_dir = "data_s4a_sparse_counts_continuous_x2",
                                                                data_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                                                output_dir = "output_s4a_sparse_counts_continuous_x2",
                                                                scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_FIXED_GAMMA_T100",
                                                                TT_use = 100L,
                                                                n1_use = 9L,
                                                                fixed_gamma_value = 0.8,
                                                                gamma_prior = c(1, 1),
                                                                n_iter = 12000L,
                                                                n_burnin = 2000L,
                                                                n_thin = 10L,
                                                                mcmc_verbose = 1000L,
                                                                overwrite_existing = FALSE,
                                                                verbose = TRUE) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    reps <- as.integer(reps)

    cat("\n=== MSSTNB S4A fixed-gamma fitting: continuous-time x2 ===\n")
    cat("data_scenario_id: ", data_scenario_id, "\n", sep = "")
    cat("fit_scenario_id : ", scenario_id, "\n", sep = "")
    cat("root:            ", root, "\n", sep = "")
    cat("data_dir:        ", file.path(root, data_dir, data_scenario_id), "\n", sep = "")
    cat("output_dir:      ", file.path(root, output_dir, data_scenario_id), "\n", sep = "")
    cat("fit_dir:         ", file.path(root, output_dir, data_scenario_id, "fits"), "\n", sep = "")
    cat("gamma_fixed:     ", fixed_gamma_value, "\n", sep = "")
    cat("n_iter/burn/thin:", paste(n_iter, n_burnin, n_thin, sep = "/"), "\n\n")

    .s4a_ct_source_s3_core(root = root, s3_script = s3_script, verbose = verbose)
    .s4a_ct_install_fixed_gamma_and_guards(fixed_gamma_value = fixed_gamma_value, verbose = verbose)

    source_manifest <- .s4a_ct_make_source_data_manifest(
        reps = reps,
        root = root,
        data_dir = data_dir,
        data_scenario_id = data_scenario_id
    )

    cat("\n=== S4A continuous-time x2 source-data check ===\n")
    print(source_manifest[, c(
        "rep_id", "x2_mode", "x2_ar", "x2_innov_sd", "mean_count",
        "median_count", "zero_prop", "total_count", "max_count",
        "beta0_star_ident", "lambda_raw_min", "lambda_raw_max"
    ), drop = FALSE])

    batch <- fit_s4a_sparse_counts_fixed_gamma_continuous_x2_batch(
        reps = reps,
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root,
        settings_override = list(
            n_iter = as.integer(n_iter),
            n_burnin = as.integer(n_burnin),
            n_thin = as.integer(n_thin)
        ),
        fixed_gamma_value = fixed_gamma_value,
        gamma_prior = gamma_prior,
        mcmc_verbose = mcmc_verbose,
        overwrite_existing = overwrite_existing
    )

    out_root <- file.path(root, output_dir, data_scenario_id)
    table_dir <- file.path(out_root, "tables")
    .s4a_ct_dir_create(table_dir)

    source_manifest_file <- file.path(table_dir, "s4a_continuous_x2_source_data_manifest.csv")
    utils::write.csv(source_manifest, source_manifest_file, row.names = FALSE)

    status_counts <- .s4a_ct_make_status_counts(batch$fit_manifest)
    compact_summary <- .s4a_ct_make_compact_summary(
        source_manifest = source_manifest,
        fit_summary = batch$fit_summary,
        fit_manifest = batch$fit_manifest
    )
    utils::write.csv(status_counts, file.path(table_dir, "s4a_fixed_gamma_fit_status_counts.csv"), row.names = FALSE)
    utils::write.csv(compact_summary, file.path(table_dir, "s4a_fixed_gamma_compact_summary.csv"), row.names = FALSE)

    run_info <- list(
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        reps = reps,
        mcmc = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin),
        fixed_gamma_value = fixed_gamma_value,
        gamma_prior = gamma_prior,
        source_manifest = source_manifest,
        fit_summary = batch$fit_summary,
        fit_manifest = batch$fit_manifest,
        status_counts = status_counts,
        compact_summary = compact_summary,
        files = list(
            source_manifest_file = source_manifest_file,
            fit_summary_file = batch$summary_file,
            fit_manifest_file = batch$manifest_file
        )
    )
    saveRDS(run_info, file.path(out_root, "run_info_S4A_sparse_counts_fixed_gamma_continuous_x2_T100.rds"))

    cat("\n=== S4A fixed-gamma continuous-time x2 fitting summary ===\n")
    print(status_counts)
    cat("\nCompact summary:\n")
    print(compact_summary)
    cat("\nMain output locations:\n")
    cat("Fits  : ", file.path(output_dir, data_scenario_id, "fits"), "\n", sep = "")
    cat("Tables: ", file.path(output_dir, data_scenario_id, "tables"), "\n", sep = "")
    cat("\nScenario 4A sparse-count continuous-time x2 fixed-gamma fitting finished.\n")

    invisible(run_info)
}
