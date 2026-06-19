## ============================================================================
## run_s4b_oracle_lambda_phi_fixed_gamma_diagnostic_T100.R
##
## Scenario 4B diagnostic: Oracle lambda + oracle phi + fixed gamma.
##
## Purpose
##   Diagnose whether the remaining S4B low/heterogeneous-exposure instability after
##   oracle-lambda fitting is caused by beta0--phi / spatial latent-effect
##   confounding or by the beta/r likelihood block itself.  This script fixes
##   lambda_tilde and phi at their true identified values and fixes gamma at truth
##   (0.8), while estimating beta0, beta1, beta2, r, and kappa.
##
## Interpretation
##   - If beta/r become stable after fixing both lambda and phi, then the original
##     failures are strong evidence of latent-effect confounding rather than a
##     generic beta-likelihood coding error.
##   - If beta/r remain unstable even with lambda and phi fixed, then inspect the
##     beta likelihood, covariate design, exposure scale, or kappa augmentation.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   data_s4b_low_exposure/S4B_LOW_HETEROGENEOUS_EXPOSURE_T100/data_rep04.rds
##   data_s4b_low_exposure/S4B_LOW_HETEROGENEOUS_EXPOSURE_T100/data_rep07.rds
##
## Main outputs
##   output_s4b_oracle_lambda_phi_fixed_gamma/
##     S4B_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4B_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100"
data_scenario_id <- "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100"

## This diagnostic targets the previously unstable S4B reps.
reps_formal <- c(4L, 7L)

## Use the same formal MCMC profile as S4B fixed-gamma runs.
n_iter <- 40000L
n_burnin <- 20000L
n_thin <- 5L

fixed_gamma_value <- 0.8
gamma_prior <- c(1, 1)
oracle_lambda_field <- "lambda_tilde_ident"
oracle_phi_field <- "phi_star_ident"

s3_core_file <- "s3_dynamic_learned_gamma.R"
s4b_data_generator_file <- "s4b_low_exposure.R"
data_dir <- "data_s4b_low_exposure"
output_dir <- "output_s4b_oracle_lambda_phi_fixed_gamma"
verbose <- 1000L
overwrite_fit <- TRUE

## Expected official S4B data settings.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "low_exposure"
expected_exposure_stress_type <- "low_heterogeneous_exposure"
expected_target_mean_multiplier <- 0.05
expected_area_log_sd <- 0.75
expected_time_log_sd <- 0.08
expected_lower_multiplier <- 0.005
expected_upper_multiplier <- 0.25

## Numerical guard settings: intentionally wide relative to the truth.
s4b_log_xi_lower <- -40
s4b_log_xi_upper <-  40
s4b_beta0_bounds <- c(-30, 10)
s4b_beta_bounds <- c(-5, 5)
s4b_kappa_bounds <- c(1e-10, 1e10)

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true <- function(x, message) {
    if (!isTRUE(x)) stop(message, call. = FALSE)
    invisible(TRUE)
}

s4b_data_file <- function(rep_id,
                          root_arg = root_dir,
                          data_dir_arg = data_dir,
                          data_scenario_id_arg = data_scenario_id) {
    file.path(
        root_arg,
        data_dir_arg,
        data_scenario_id_arg,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}


cv_s4b_oracle <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0L || !is.finite(mean(x, na.rm = TRUE)) || mean(x, na.rm = TRUE) == 0) {
        return(NA_real_)
    }
    stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

s4b_area_exposure_groups <- function(dat, n_groups = 3L) {
    n_groups <- as.integer(n_groups)
    if (n_groups < 2L) n_groups <- 2L
    area_mean_e <- colMeans(dat$e)
    ord <- order(area_mean_e, decreasing = FALSE)
    group_id <- rep(NA_integer_, length(area_mean_e))
    group_id[ord] <- cut(seq_along(ord), breaks = n_groups, labels = FALSE)
    if (n_groups == 3L) {
        labels <- c("low_exposure", "middle_exposure", "high_exposure")
    } else {
        labels <- paste0("group_", seq_len(n_groups))
    }
    list(group_id = group_id, labels = labels)
}

s4b_exposure_group_count_summary <- function(dat, n_groups = 3L) {
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        y_g <- dat$y_coarse[, jj, drop = FALSE]
        e_g <- dat$e[, jj, drop = FALSE]
        e_ref_g <- dat$e_reference[, jj, drop = FALSE]
        mult_g <- dat$exposure_multiplier[, jj, drop = FALSE]
        data.frame(
            exposure_group = g$labels[gg],
            n_areas = length(jj),
            area_ids = paste(jj, collapse = ","),
            mean_exposure = mean(e_g),
            mean_reference_exposure = mean(e_ref_g),
            mean_multiplier = mean(mult_g),
            mean_count = mean(y_g),
            zero_prop = mean(y_g == 0),
            total_count = sum(y_g),
            max_count = max(y_g),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

s4b_lambda_group_recovery <- function(lambda_draws, lambda_truth, dat, n_groups = 3L) {
    if (is.null(lambda_draws) || length(dim(lambda_draws)) != 3L) return(data.frame())
    lambda_mean <- apply(lambda_draws, c(2, 3), mean, na.rm = TRUE)
    lambda_truth <- as.matrix(lambda_truth)
    if (!identical(dim(lambda_mean), dim(lambda_truth))) {
        stop("lambda posterior mean and lambda truth have incompatible dimensions.", call. = FALSE)
    }
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        lm <- lambda_mean[, jj, drop = FALSE]
        lt <- lambda_truth[, jj, drop = FALSE]
        log_lm <- log(pmax(lm, .Machine$double.xmin))
        log_lt <- log(pmax(lt, .Machine$double.xmin))
        data.frame(
            exposure_group = g$labels[gg],
            lambda_rmse = sqrt(mean((lm - lt)^2)),
            log_lambda_rmse = sqrt(mean((log_lm - log_lt)^2)),
            cor_log_lambda = suppressWarnings(stats::cor(as.numeric(log_lm), as.numeric(log_lt))),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

s4b_phi_group_recovery <- function(phi_draws, phi_truth, dat, n_groups = 3L) {
    if (is.null(phi_draws) || length(dim(phi_draws)) != 2L) return(data.frame())
    phi_mean <- colMeans(phi_draws, na.rm = TRUE)
    phi_truth <- as.numeric(phi_truth)
    if (length(phi_mean) != length(phi_truth)) {
        stop("phi posterior mean and phi truth have incompatible lengths.", call. = FALSE)
    }
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        err <- phi_mean[jj] - phi_truth[jj]
        data.frame(
            exposure_group = g$labels[gg],
            phi_rmse = sqrt(mean(err^2)),
            phi_mae = mean(abs(err)),
            phi_cor = suppressWarnings(stats::cor(phi_mean[jj], phi_truth[jj])),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

prefix_group_metrics <- function(df, prefix) {
    if (is.null(df) || nrow(df) == 0L) return(list())
    out <- list()
    for (ii in seq_len(nrow(df))) {
        grp <- df$exposure_group[ii]
        for (nm in setdiff(names(df), "exposure_group")) {
            out[[paste0(prefix, "_", grp, "_", nm)]] <- df[[nm]][ii]
        }
    }
    out
}

check_s4b_oracle_source_dataset <- function(data_file,
                                            TT_expected = TT_use,
                                            n1_expected = n1_use,
                                            beta0_ident_abs_limit = 20) {
    assert_file_exists(data_file, "S4B source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "exposure_stress_type", "y_coarse", "e",
        "e_reference", "exposure_multiplier", "x1", "x2", "lambda_tilde",
        "lambda_tilde_ident", "gamma_star", "beta0_star", "beta0_star_ident",
        "beta_star_ident", "phi_star_ident", "target_mean_exposure_multiplier",
        "realized_mean_exposure_multiplier", "realized_min_exposure_multiplier",
        "realized_max_exposure_multiplier", "mean_count", "zero_prop", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true(
        length(missing) == 0L,
        paste("Dataset is missing required S4B fields:", paste(missing, collapse = ", "))
    )

    assert_true(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_expected), as.integer(n1_expected))),
        paste0("y_coarse dimension is not ", TT_expected, " x ", n1_expected, ".")
    )

    for (nm in c("e", "e_reference", "exposure_multiplier", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        assert_true(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    assert_true(identical(dat$stress_type, expected_stress_type),
                paste0("stress_type is not ", expected_stress_type, ". Got ", dat$stress_type, "."))
    assert_true(identical(dat$exposure_stress_type, expected_exposure_stress_type),
                paste0("exposure_stress_type is not ", expected_exposure_stress_type,
                       ". Got ", dat$exposure_stress_type, "."))
    assert_true(abs(dat$target_mean_exposure_multiplier - expected_target_mean_multiplier) < 1e-8,
                paste0("target_mean_exposure_multiplier is not ", expected_target_mean_multiplier,
                       ". Got ", dat$target_mean_exposure_multiplier, "."))

    if (!isTRUE(all.equal(dat$e, dat$e_reference * dat$exposure_multiplier,
                          tolerance = 1e-8, check.attributes = FALSE))) {
        stop("dat$e is not equal to dat$e_reference * dat$exposure_multiplier.", call. = FALSE)
    }
    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$e_reference)) || any(dat$e_reference <= 0) ||
        any(!is.finite(dat$exposure_multiplier)) || any(dat$exposure_multiplier <= 0) ||
        any(dat$exposure_multiplier >= 1)) {
        stop("S4B exposure values/multipliers must be positive, finite, and multipliers must be < 1.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("dat$lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    assert_true(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident)
    )

    if (exists("validate_s4b_low_exposure_data", mode = "function")) {
        validate_s4b_low_exposure_data(dat)
    }

    count_groups <- s4b_exposure_group_count_summary(dat, n_groups = 3L)
    low_g <- count_groups[count_groups$exposure_group == "low_exposure", , drop = FALSE]
    high_g <- count_groups[count_groups$exposure_group == "high_exposure", , drop = FALSE]

    list(
        dat = dat,
        mean_count = mean(dat$y_coarse),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        zero_prop = mean(dat$y_coarse == 0),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        reference_mean_count = if (!is.null(dat$reference_count_summary)) dat$reference_count_summary$mean_count else NA_real_,
        reference_zero_prop = if (!is.null(dat$reference_count_summary)) dat$reference_count_summary$zero_prop else NA_real_,
        target_mean_multiplier = dat$target_mean_exposure_multiplier,
        realized_mean_multiplier = mean(dat$exposure_multiplier),
        realized_min_multiplier = min(dat$exposure_multiplier),
        realized_max_multiplier = max(dat$exposure_multiplier),
        reference_mean_exposure = mean(dat$e_reference),
        mean_exposure = mean(dat$e),
        exposure_mean_ratio = mean(dat$e) / mean(dat$e_reference),
        exposure_cv = cv_s4b_oracle(dat$e),
        multiplier_cv = cv_s4b_oracle(dat$exposure_multiplier),
        area_multiplier_cv = cv_s4b_oracle(colMeans(dat$exposure_multiplier)),
        lowest_exposure_group_mean_count = low_g$mean_count %||% NA_real_,
        lowest_exposure_group_zero_prop = low_g$zero_prop %||% NA_real_,
        highest_exposure_group_mean_count = high_g$mean_count %||% NA_real_,
        highest_exposure_group_zero_prop = high_g$zero_prop %||% NA_real_,
        count_exposure_cor_area = suppressWarnings(stats::cor(colMeans(dat$y_coarse), colMeans(dat$e))),
        beta0_star_ident = dat$beta0_star_ident,
        lambda_oracle_min = min(dat[[oracle_lambda_field]]),
        lambda_oracle_max = max(dat[[oracle_lambda_field]])
    )
}

## ---- source Scenario 3/S4B model code --------------------------------------
assert_file_exists(file.path(root_dir, s3_core_file), "Scenario 3 core script")
assert_file_exists(file.path(root_dir, s4b_data_generator_file), "Scenario 4B data-generation script")
source(file.path(root_dir, s4b_data_generator_file))
source_s4b_low_exposure(root = root_dir, verbose = FALSE)

## ---- numerical guards and local overrides ----------------------------------
compute_s3_xi <- function(e, x1, x2, beta0, beta, phi) {
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

    xi <- matrix(NA_real_, nrow(e), ncol(e))
    for (j in seq_len(ncol(e))) {
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        log_xi_j <- log(e[, j]) + eta_j
        log_xi_j <- pmin(pmax(log_xi_j, s4b_log_xi_lower), s4b_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("Safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

.s4b_oracle_guard_env <- new.env(parent = emptyenv())
reset_oracle_guards <- function() {
    .s4b_oracle_guard_env$n_beta_guard <- 0L
    .s4b_oracle_guard_env$n_kappa_guard <- 0L
    invisible(TRUE)
}
get_oracle_guards <- function() {
    list(
        n_beta_guard = .s4b_oracle_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4b_oracle_guard_env$n_kappa_guard %||% 0L
    )
}
reset_oracle_guards()

.update_beta_s3_original <- update_beta
update_beta <- function(beta_current, ...) {
    res <- .update_beta_s3_original(beta_current = beta_current, ...)
    smp <- res$sample
    bad <- FALSE
    if (length(smp) < 3L || any(!is.finite(smp))) {
        bad <- TRUE
    } else {
        bad <- smp[1] < s4b_beta0_bounds[1] || smp[1] > s4b_beta0_bounds[2] ||
            any(smp[-1] < s4b_beta_bounds[1] | smp[-1] > s4b_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s4b_oracle_guard_env$n_beta_guard <- (.s4b_oracle_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s4b_oracle_guard_rejected <- TRUE
    } else {
        res$s4b_oracle_guard_rejected <- FALSE
    }
    res
}

.update_kappa_s3_original <- update_kappa
update_kappa <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
    y <- as.matrix(y_coarse)
    L <- as.matrix(lambda_tilde)
    X <- as.matrix(xi)
    if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
        stop("Oracle safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
    }
    r_vec <- as.numeric(r)
    if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
    if (length(r_vec) != ncol(y)) {
        stop("Oracle safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
    }
    R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
    shape <- y + R
    rate <- X * L + R
    guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
    if (any(guard)) {
        .s4b_oracle_guard_env$n_kappa_guard <- (.s4b_oracle_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s4b_oracle_guard_env$n_kappa_guard <- (.s4b_oracle_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s4b_kappa_bounds[1]), s4b_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s4b_oracle_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

cat(sprintf(
    "Using oracle-lambda-phi diagnostic guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s4b_beta0_bounds[1], s4b_beta0_bounds[2],
    s4b_beta_bounds[1], s4b_beta_bounds[2],
    s4b_log_xi_lower, s4b_log_xi_upper
))

## ---- oracle-lambda MCMC -----------------------------------------------------
run_s4b_oracle_lambda_phi_fixed_gamma_mcmc <- function(dat,
                                                   settings = MCMC_SETTINGS,
                                                   priors = MCMC_PRIORS,
                                                   spatial = build_s3_spatial(),
                                                   oracle_lambda = dat[[oracle_lambda_field]],
                                                   oracle_phi = dat[[oracle_phi_field]],
                                                   gamma_value = fixed_gamma_value,
                                                   verbose = 1000L) {
    validate_s3_fit_objects(dat, priors, spatial)

    n_iter_local <- as.integer(settings$n_iter %||% n_iter)
    n_burnin_local <- as.integer(settings$n_burnin %||% n_burnin)
    n_thin_local <- as.integer(settings$n_thin %||% n_thin)
    adapt_interval <- as.integer(settings$adapt_interval %||% 50L)

    if (n_iter_local <= 0L || n_burnin_local < 0L || n_burnin_local >= n_iter_local || n_thin_local <= 0L) {
        stop("Invalid MCMC iteration settings.", call. = FALSE)
    }

    TT_local <- as.integer(dat$TT)
    n1_local <- as.integer(dat$n1)
    p <- length(priors$beta_mean)
    n_stored <- as.integer(floor((n_iter_local - n_burnin_local) / n_thin_local))

    oracle_lambda <- as.matrix(oracle_lambda)
    if (!identical(dim(oracle_lambda), dim(dat$y_coarse)) || any(!is.finite(oracle_lambda)) || any(oracle_lambda <= 0)) {
        stop("oracle_lambda must be positive finite and have the same dimensions as y_coarse.", call. = FALSE)
    }
    oracle_phi <- as.numeric(oracle_phi)
    if (length(oracle_phi) != n1_local || any(!is.finite(oracle_phi))) {
        stop("oracle_phi must be finite and have length n1.", call. = FALSE)
    }

    gamma_init <- rep(gamma_value, n1_local)
    state <- initialise_s3_state(dat, settings, priors, spatial, gamma_init = gamma_init)
    state$lambda_tilde <- oracle_lambda
    state$phi <- oracle_phi
    state$gamma <- gamma_init
    state$gamma_common <- gamma_value
    state$xi <- compute_s3_xi(dat$e, dat$x1, dat$x2, state$beta0, state$beta, state$phi)

    samples <- list(
        beta0 = numeric(n_stored),
        beta = matrix(NA_real_, nrow = n_stored, ncol = p),
        phi = matrix(NA_real_, nrow = n_stored, ncol = n1_local),
        tau_phi = numeric(n_stored),
        r = matrix(NA_real_, nrow = n_stored, ncol = n1_local),
        gamma = matrix(NA_real_, nrow = n_stored, ncol = n1_local),
        gamma_common = numeric(n_stored),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT_local, n1_local)),
        kappa = array(NA_real_, dim = c(n_stored, TT_local, n1_local)),
        loglik = numeric(n_stored),
        loglik_nb = numeric(n_stored)
    )

    diagnostics <- list(
        loglik_trace = rep(NA_real_, n_iter_local),
        loglik_nb_trace = rep(NA_real_, n_iter_local),
        beta_n_reject = rep(NA_real_, n_iter_local),
        phi_accept_trace = rep(FALSE, n_iter_local),
        r_accept_trace = matrix(FALSE, nrow = n_iter_local, ncol = n1_local),
        gamma_accept_trace = rep(FALSE, n_iter_local),
        phi_log_alpha = rep(NA_real_, n_iter_local),
        r_log_alpha = matrix(NA_real_, nrow = n_iter_local, ncol = n1_local),
        gamma_log_alpha = rep(NA_real_, n_iter_local),
        gamma_common_trace = rep(gamma_value, n_iter_local),
        gamma_proposal_sd_trace = rep(NA_real_, n_iter_local),
        lambda_min_trace = rep(min(oracle_lambda), n_iter_local),
        lambda_max_trace = rep(max(oracle_lambda), n_iter_local),
        kappa_mean_trace = rep(NA_real_, n_iter_local),
        recenter_error = rep(0, n_iter_local),
        gamma_prior = gamma_prior,
        oracle_lambda_fixed = TRUE,
        oracle_phi_fixed = TRUE,
        oracle_lambda_source = oracle_lambda_field,
        oracle_phi_source = oracle_phi_field
    )

    phi_accept_window <- 0L
    r_accept_window <- rep(0L, n1_local)
    store_idx <- 0L
    start_time <- proc.time()

    for (iter in seq_len(n_iter_local)) {
        ## 1. beta update using marginal NB likelihood, conditional on oracle lambda.
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = oracle_lambda,
            phi = state$phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r,
            use_preconditioned = settings$use_preconditioned_beta %||% TRUE
        )
        state$beta0 <- beta_result$sample[1]
        state$beta <- beta_result$sample[2:3]
        state$xi <- compute_s3_xi(dat$e, dat$x1, dat$x2, state$beta0, state$beta, state$phi)

        ## 2. phi and tau_phi are fixed at the oracle identified truth in this diagnostic.
        state$phi <- oracle_phi
        state$xi <- compute_s3_xi(dat$e, dat$x1, dat$x2, state$beta0, state$beta, state$phi)

        ## 3. r update using marginal NB likelihood, conditional on oracle lambda and oracle phi.
        r_result <- update_r_s3(
            r_current = state$r,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = oracle_lambda,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            priors = priors,
            mh_sd = state$r_proposal_sd,
            method = "marginal_nb",
            return_diag = TRUE
        )
        state$r <- r_result$r

        ## 5. kappa update for diagnostic Poisson-augmentation likelihood only.
        kappa_result <- update_kappa(
            y_coarse = dat$y_coarse,
            lambda_tilde = oracle_lambda,
            xi = state$xi,
            r = state$r,
            return_diag = TRUE
        )
        state$kappa <- kappa_result$kappa

        ## 6. lambda and gamma are fixed in this diagnostic.
        state$lambda_tilde <- oracle_lambda
        state$gamma <- gamma_init
        state$gamma_common <- gamma_value

        loglik <- compute_s3_loglik_poisson_aug(
            y_coarse = dat$y_coarse,
            xi = state$xi,
            lambda_tilde = oracle_lambda,
            kappa = state$kappa
        )
        loglik_nb <- compute_s3_loglik_nb(
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            r = state$r,
            lambda_tilde = oracle_lambda
        )

        diagnostics$loglik_trace[iter] <- loglik
        diagnostics$loglik_nb_trace[iter] <- loglik_nb
        diagnostics$beta_n_reject[iter] <- beta_result$n_reject %||% NA_real_
        diagnostics$phi_accept_trace[iter] <- NA
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
        diagnostics$phi_log_alpha[iter] <- NA_real_
        diagnostics$r_log_alpha[iter, ] <- r_result$diag$log_alpha %||% rep(NA_real_, n1_local)
        diagnostics$kappa_mean_trace[iter] <- kappa_result$diag$mean_kappa %||% mean(state$kappa)

        r_accept_window <- r_accept_window + as.integer(r_result$accept)

        if (iter <= n_burnin_local && iter %% adapt_interval == 0L) {
            state$r_proposal_sd <- adapt_sd_s3(
                current_sd = state$r_proposal_sd,
                n_accept = r_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$r_target_accept %||% 0.30
            )
            r_accept_window <- rep(0L, n1_local)
        }

        if (iter > n_burnin_local && (iter - n_burnin_local) %% n_thin_local == 0L) {
            store_idx <- store_idx + 1L
            samples$beta0[store_idx] <- state$beta0
            samples$beta[store_idx, ] <- state$beta
            samples$phi[store_idx, ] <- state$phi
            samples$tau_phi[store_idx] <- state$tau_phi
            samples$r[store_idx, ] <- state$r
            samples$gamma[store_idx, ] <- state$gamma
            samples$gamma_common[store_idx] <- state$gamma_common
            samples$lambda_tilde[store_idx, , ] <- oracle_lambda
            samples$kappa[store_idx, , ] <- state$kappa
            samples$loglik[store_idx] <- loglik
            samples$loglik_nb[store_idx] <- loglik_nb
        }

        if (verbose > 0L && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1L, iter - 99L)
            phi_rate <- NA_real_
            r_rate <- mean(diagnostics$r_accept_trace[i0:iter, , drop = FALSE])
            beta_rej <- mean(diagnostics$beta_n_reject[i0:iter], na.rm = TRUE)
            cat(sprintf(
                paste0(
                    "iter %5d/%d [%.0fs] loglik_nb=%.1f beta0=%.3f ",
                    "beta=(%.3f, %.3f) r_mean=%.2f oracle_lambda=[%.3f, %.3f] ",
                    "gamma=%.3f | phi_acc=%.2f r_acc=%.2f beta_rej=%.1f\n"
                ),
                iter, n_iter_local, elapsed, loglik_nb, state$beta0, state$beta[1],
                state$beta[2], mean(state$r), min(oracle_lambda), max(oracle_lambda),
                state$gamma_common, phi_rate, r_rate, beta_rej
            ))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]
    diagnostics$elapsed_sec <- elapsed_total
    diagnostics$phi_accept_rate <- NA_real_
    diagnostics$r_accept_rate <- colMeans(diagnostics$r_accept_trace)
    diagnostics$gamma_accept_rate <- NA_real_
    diagnostics$beta_mean_n_reject <- mean(diagnostics$beta_n_reject, na.rm = TRUE)
    diagnostics$phi_proposal_sd_final <- state$phi_proposal_sd
    diagnostics$r_proposal_sd_final <- state$r_proposal_sd
    diagnostics$gamma_proposal_sd_final <- NA_real_
    diagnostics$gamma_mean <- mean(samples$gamma_common, na.rm = TRUE)
    diagnostics$gamma_sd <- stats::sd(samples$gamma_common, na.rm = TRUE)

    list(
        samples = samples,
        diagnostics = diagnostics,
        final_state = state,
        settings = settings,
        priors = priors,
        spatial = list(H = spatial$H, B_ICAR = spatial$B_ICAR, BHB = spatial$BHB),
        metadata = list(
            model = "S4B oracle-lambda-phi fixed-gamma diagnostic",
            fixed_lambda = TRUE,
            fixed_phi = TRUE,
            oracle_lambda = TRUE,
            oracle_phi = TRUE,
            oracle_lambda_source = oracle_lambda_field,
            oracle_phi_source = oracle_phi_field,
            dynamic_lambda = TRUE,
            lambda_sampled_in_fit = FALSE,
            fixed_gamma = TRUE,
            learned_gamma = FALSE,
            gamma_fixed_value = gamma_value,
            updated_blocks = c("beta", "r", "kappa"),
            disabled_blocks = c("gamma", "lambda", "phi", "tau_phi", "delta", "omega"),
            uses_marginal_nb_for_beta_phi_r = TRUE,
            uses_kappa_for_diagnostic_poisson_aug_loglik = TRUE,
            uses_recentered_identified_parameterization = TRUE
        ),
        n_stored = n_stored
    )
}

fit_s4b_oracle_lambda_phi_one_rep <- function(rep_id,
                                          scenario_id_arg = scenario_id,
                                          data_scenario_id_arg = data_scenario_id,
                                          data_dir_arg = data_dir,
                                          output_dir_arg = output_dir,
                                          root_arg = root_dir,
                                          settings_override = list(),
                                          priors = MCMC_PRIORS,
                                          spatial = build_s3_spatial(),
                                          fixed_gamma_value_arg = fixed_gamma_value,
                                          verbose_arg = verbose,
                                          save_result = TRUE,
                                          return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- s4b_data_file(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
    chk <- check_s4b_oracle_source_dataset(dat_file)
    dat <- chk$dat
    validate_s3_data(dat)

    settings <- MCMC_SETTINGS
    settings$n_iter <- n_iter
    settings$n_burnin <- n_burnin
    settings$n_thin <- n_thin
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) settings[[nm]] <- settings_override[[nm]]
    }

    oracle_lambda <- dat[[oracle_lambda_field]]
    oracle_phi <- dat[[oracle_phi_field]]

    cat(sprintf("=== S4B oracle-lambda-phi fixed-gamma diagnostic: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Fixed     : gamma = %.3f; lambda = dat$%s; phi = dat$%s\n", fixed_gamma_value_arg, oracle_lambda_field, oracle_phi_field))
    cat("Estimated : beta0, beta1, beta2, r, kappa\n")
    cat("Disabled  : lambda, phi, tau_phi, gamma, delta, omega updates\n\n")

    reset_oracle_guards()
    fit <- run_s4b_oracle_lambda_phi_fixed_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        oracle_lambda = oracle_lambda,
        oracle_phi = oracle_phi,
        gamma_value = fixed_gamma_value_arg,
        verbose = verbose_arg
    )

    guard_counts <- get_oracle_guards()
    fit$diagnostics$s4b_oracle_lamphi_numeric_guards <- guard_counts
    fit$diagnostics$s4b_oracle_lamphi_beta_guard_count <- guard_counts$n_beta_guard
    fit$diagnostics$s4b_oracle_lamphi_kappa_guard_count <- guard_counts$n_kappa_guard

    summary <- summarise_s3_dynamic_learned_gamma_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id_arg,
        rep_id = as.integer(rep_id)
    )

    gamma_truth_mean <- mean(dat$gamma_star %||% fixed_gamma_value_arg, na.rm = TRUE)
    gamma_fixed_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    summary$lambda_sampled_in_fit <- FALSE
    summary$oracle_lambda_fixed_in_fit <- TRUE
    summary$oracle_phi_fixed_in_fit <- TRUE
    summary$oracle_lambda_source <- oracle_lambda_field
    summary$oracle_phi_source <- oracle_phi_field
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
    summary$exposure_stress_type <- dat$exposure_stress_type %||% NA_character_
    summary$target_mean_exposure_multiplier <- dat$target_mean_exposure_multiplier %||% NA_real_
    summary$realized_mean_exposure_multiplier <- mean(dat$exposure_multiplier)
    summary$realized_min_exposure_multiplier <- min(dat$exposure_multiplier)
    summary$realized_max_exposure_multiplier <- max(dat$exposure_multiplier)
    summary$reference_mean_count <- chk$reference_mean_count
    summary$reference_zero_prop <- chk$reference_zero_prop
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count
    summary$reference_mean_exposure <- chk$reference_mean_exposure
    summary$mean_exposure <- chk$mean_exposure
    summary$exposure_mean_ratio <- chk$exposure_mean_ratio
    summary$exposure_cv <- chk$exposure_cv
    summary$multiplier_cv <- chk$multiplier_cv
    summary$area_multiplier_cv <- chk$area_multiplier_cv
    summary$lowest_exposure_group_mean_count <- chk$lowest_exposure_group_mean_count
    summary$lowest_exposure_group_zero_prop <- chk$lowest_exposure_group_zero_prop
    summary$highest_exposure_group_mean_count <- chk$highest_exposure_group_mean_count
    summary$highest_exposure_group_zero_prop <- chk$highest_exposure_group_zero_prop
    summary$count_exposure_cor_area <- chk$count_exposure_cor_area
    summary$lambda_oracle_min <- chk$lambda_oracle_min
    summary$lambda_oracle_max <- chk$lambda_oracle_max

    lambda_group <- s4b_lambda_group_recovery(
        lambda_draws = fit$samples$lambda_tilde,
        lambda_truth = dat$lambda_tilde_ident %||% dat$lambda_tilde,
        dat = dat,
        n_groups = 3L
    )
    phi_group <- s4b_phi_group_recovery(
        phi_draws = fit$samples$phi,
        phi_truth = dat$phi_star_ident %||% dat$phi_star,
        dat = dat,
        n_groups = 3L
    )
    count_group <- s4b_exposure_group_count_summary(dat, n_groups = 3L)
    group_metrics <- c(
        prefix_group_metrics(count_group[, c("exposure_group", "n_areas", "mean_exposure", "mean_multiplier", "mean_count", "zero_prop", "total_count"), drop = FALSE], "count"),
        prefix_group_metrics(lambda_group, "lambda"),
        prefix_group_metrics(phi_group, "phi")
    )
    if (length(group_metrics) > 0L) {
        for (nm in names(group_metrics)) summary[[nm]] <- group_metrics[[nm]]
    }

    summary$s4b_oracle_lamphi_beta_guard_count <- guard_counts$n_beta_guard
    summary$s4b_oracle_lamphi_kappa_guard_count <- guard_counts$n_kappa_guard

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id_arg,
            data_scenario_id = data_scenario_id_arg,
            stress_type = dat$stress_type,
            exposure_stress_type = dat$exposure_stress_type,
            target_mean_exposure_multiplier = dat$target_mean_exposure_multiplier,
            realized_mean_exposure_multiplier = mean(dat$exposure_multiplier),
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir_arg,
            fit_file_prefix = "fit_S4B_oracle_lambda_phi_fixed_gamma_rep",
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root_arg, output_dir_arg, scenario_id_arg)
        ensure_dir(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4B_oracle_lambda_phi_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4B_oracle_lambda_phi_fixed_gamma_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s4b_oracle_lambda_phi_batch <- function(reps = reps_formal,
                                        scenario_id_arg = scenario_id,
                                        data_scenario_id_arg = data_scenario_id,
                                        data_dir_arg = data_dir,
                                        output_dir_arg = output_dir,
                                        root_arg = root_dir,
                                        fixed_gamma_value_arg = fixed_gamma_value,
                                        verbose_arg = verbose,
                                        overwrite_existing = overwrite_fit) {
    out_dir <- file.path(root_arg, output_dir_arg, scenario_id_arg)
    ensure_dir(out_dir)

    summaries <- list()
    fit_manifest <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4B_oracle_lambda_phi_fixed_gamma_rep", rr, ".rds"))
        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
        } else {
            fit <- fit_s4b_oracle_lambda_phi_one_rep(
                rep_id = rep_id,
                scenario_id_arg = scenario_id_arg,
                data_scenario_id_arg = data_scenario_id_arg,
                data_dir_arg = data_dir_arg,
                output_dir_arg = output_dir_arg,
                root_arg = root_arg,
                fixed_gamma_value_arg = fixed_gamma_value_arg,
                verbose_arg = verbose_arg,
                save_result = TRUE,
                return_result = TRUE
            )
            summaries[[rr]] <- fit$summary
        }
        data_file <- s4b_data_file(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
        fit_manifest[[rr]] <- data.frame(
            scenario_id = scenario_id_arg,
            data_scenario_id = data_scenario_id_arg,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            fit_file = fit_file,
            fit_exists = file.exists(fit_file),
            stringsAsFactors = FALSE
        )
    }

    summary_all <- do.call(rbind, summaries)
    manifest_all <- do.call(rbind, fit_manifest)

    summary_file <- file.path(out_dir, "summary_S4B_oracle_lambda_phi_fixed_gamma_all_reps.csv")
    manifest_file <- file.path(out_dir, "s4b_oracle_lamphi_fit_manifest.csv")
    run_info_file <- file.path(out_dir, "run_info_S4B_oracle_lambda_phi_fixed_gamma_diagnostic_T100.rds")

    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    utils::write.csv(manifest_all, manifest_file, row.names = FALSE)
    saveRDS(
        list(
            scenario_id = scenario_id_arg,
            data_scenario_id = data_scenario_id_arg,
            reps = reps,
            n_iter = n_iter,
            n_burnin = n_burnin,
            n_thin = n_thin,
            fixed_gamma_value = fixed_gamma_value_arg,
            oracle_lambda_field = oracle_lambda_field,
            created_at = Sys.time(),
            output_dir = out_dir
        ),
        run_info_file
    )

    cat("\n=== S4B oracle-lambda-phi diagnostic summary ===\n")
    print(summary_all[, intersect(c(
        "rep_id", "mean_count", "zero_prop", "beta0_mean", "beta1_mean", "beta2_mean",
        "r_mean", "phi_rmse", "phi_cor", "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "s4b_oracle_lamphi_beta_guard_count", "s4b_oracle_lamphi_kappa_guard_count"
    ), names(summary_all)), drop = FALSE])
    cat("\nSaved combined summary: ", summary_file, "\n", sep = "")
    cat("Saved manifest        : ", manifest_file, "\n", sep = "")

    invisible(summary_all)
}

## ---- run diagnostic ---------------------------------------------------------
settings_override <- list(
    n_iter = n_iter,
    n_burnin = n_burnin,
    n_thin = n_thin
)

s4b_oracle_summary <- fit_s4b_oracle_lambda_phi_batch(
    reps = reps_formal,
    scenario_id_arg = scenario_id,
    data_scenario_id_arg = data_scenario_id,
    data_dir_arg = data_dir,
    output_dir_arg = output_dir,
    root_arg = root_dir,
    fixed_gamma_value_arg = fixed_gamma_value,
    verbose_arg = verbose,
    overwrite_existing = overwrite_fit
)

cat("\nS4B oracle-lambda-phi fixed-gamma diagnostic finished successfully.\n")
