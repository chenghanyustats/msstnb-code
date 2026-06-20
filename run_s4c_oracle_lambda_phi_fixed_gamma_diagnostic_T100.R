## ============================================================================
## run_s4c_oracle_lambda_phi_fixed_gamma_diagnostic_T100.R
##
## Scenario 4C diagnostic: Oracle lambda + oracle phi + fixed gamma.
##
## Purpose
##   Diagnose whether the remaining S4C strong-overdispersion instability after
##   oracle-lambda fitting is caused by beta0--phi / spatial-effect confounding
##   or by the beta/r/kappa likelihood block itself. This script fixes both
##   lambda_tilde and phi at their true identified values and fixes gamma at
##   truth (0.8), while estimating beta0, beta1, beta2, r, and kappa.
##
## Interpretation
##   - If beta/r/kappa become stable after fixing both lambda and phi, then the
##     rep 07 failure is strong evidence of latent-effect confounding rather
##     than a generic beta/r/kappa coding error.
##   - If beta/r/kappa remain unstable even with lambda and phi fixed, then the
##     S4C failure is not only a temporal/spatial latent-effect ridge and the
##     beta likelihood, r update, or kappa augmentation should be inspected.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   s4c_small_r_overdispersion_v2.R  (preferred) or s4c_small_r_overdispersion.R
##   data_s4c_overdispersion/S4C_STRONG_OVERDISPERSION_T100/data_rep07.rds
##
## Main outputs
##   output_s4c_oracle_lambda_phi_fixed_gamma/
##     S4C_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4C_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100"
data_scenario_id <- "S4C_STRONG_OVERDISPERSION_T100"

## This diagnostic targets the one numerical-instability S4C replicate after
## the sampled-lambda fit and the oracle-lambda diagnostic.
reps_formal <- c(7L)

## Use the same formal MCMC profile as S4C fixed-gamma/oracle-lambda runs.
n_iter <- 40000L
n_burnin <- 20000L
n_thin <- 5L

fixed_gamma_value <- 0.8
gamma_prior <- c(1, 1)
oracle_lambda_field <- "lambda_tilde_ident"
oracle_phi_field <- "phi_star_ident"

s3_core_file <- "s3_dynamic_learned_gamma.R"
s4c_data_generator_file <- if (file.exists(file.path(root_dir, "s4c_small_r_overdispersion_v2.R"))) {
    "s4c_small_r_overdispersion_v2.R"
} else {
    "s4c_small_r_overdispersion.R"
}

data_dir <- "data_s4c_overdispersion"
output_dir <- "output_s4c_oracle_lambda_phi_fixed_gamma"
verbose <- 1000L
overwrite_fit <- TRUE

## Expected official S4C data settings.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "small_r_overdispersion"
expected_r_reference_truth <- 15
expected_r_stress_truth <- 3

## Numerical guard settings: intentionally wide relative to truth.
s4c_oracle_log_xi_lower <- -40
s4c_oracle_log_xi_upper <-  40
s4c_oracle_beta0_bounds <- c(-30, 10)
s4c_oracle_beta_bounds <- c(-5, 5)
s4c_oracle_kappa_bounds <- c(1e-10, 1e10)

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_dir_s4c_oracle <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s4c_oracle <- function(path, label = "file") {
    if (!file.exists(path)) stop("Required ", label, " not found: ", path, call. = FALSE)
    invisible(TRUE)
}

assert_true_s4c_oracle <- function(x, message) {
    if (!isTRUE(x)) stop(message, call. = FALSE)
    invisible(TRUE)
}

cv_s4c_oracle <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x, na.rm = TRUE)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) return(NA_real_)
    stats::sd(x, na.rm = TRUE) / mx
}

safe_ratio_s4c_oracle <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}

count_stats_s4c_oracle <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(yy, probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = cv_s4c_oracle(yy),
        variance_to_mean = safe_ratio_s4c_oracle(vv, mn),
        max_to_mean = safe_ratio_s4c_oracle(max(yy), mn),
        q95_count = as.numeric(qs[6]),
        q99_count = as.numeric(qs[7]),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

kappa_stats_s4c_oracle <- function(kappa, prefix = "") {
    kk <- as.numeric(kappa)
    qs <- stats::quantile(kk, probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        kappa_mean = mean(kk),
        kappa_sd = stats::sd(kk),
        kappa_cv = cv_s4c_oracle(kk),
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

s4c_oracle_source_data_path <- function(rep_id,
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

summarise_s4c_oracle_kappa_recovery <- function(kappa_draws, kappa_truth) {
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
        kappa_truth_cv = cv_s4c_oracle(kappa_truth),
        kappa_post_mean_cv = cv_s4c_oracle(kappa_mean),
        stringsAsFactors = FALSE
    )
}

r_recovery_s4c_oracle <- function(r_draws, r_truth) {
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

check_s4c_oracle_source_dataset <- function(data_file,
                                            TT_expected = TT_use,
                                            n1_expected = n1_use,
                                            expected_stress_type_arg = expected_stress_type,
                                            expected_r_reference_truth_arg = expected_r_reference_truth,
                                            expected_r_stress_truth_arg = expected_r_stress_truth,
                                            beta0_ident_abs_limit = 20) {
    assert_file_exists_s4c_oracle(data_file, "S4C source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "y_coarse", "y_coarse_reference", "e",
        "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "gamma_star",
        "beta0_star", "beta0_star_ident", "beta_star_ident", "phi_star_ident",
        "kappa", "kappa_reference", "r_star", "r_reference_truth", "r_stress_truth",
        "reference_count_summary", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s4c_oracle(
        length(missing) == 0L,
        paste("Dataset is missing required S4C fields:", paste(missing, collapse = ", "))
    )

    assert_true_s4c_oracle(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_expected), as.integer(n1_expected))),
        paste0("y_coarse dimension is not TT by n1. Got ", paste(dim(dat$y_coarse), collapse = " x "), ".")
    )

    for (nm in c("y_coarse_reference", "e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "kappa", "kappa_reference")) {
        assert_true_s4c_oracle(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimensions as dat$y_coarse.")
        )
    }

    assert_true_s4c_oracle(
        identical(dat$stress_type, expected_stress_type_arg),
        paste0("stress_type is not ", expected_stress_type_arg, ". Got ", dat$stress_type, ".")
    )
    assert_true_s4c_oracle(
        abs(dat$r_reference_truth - expected_r_reference_truth_arg) < 1e-8,
        paste0("r_reference_truth is not ", expected_r_reference_truth_arg, ". Got ", dat$r_reference_truth, ".")
    )
    assert_true_s4c_oracle(
        abs(dat$r_stress_truth - expected_r_stress_truth_arg) < 1e-8,
        paste0("r_stress_truth is not ", expected_r_stress_truth_arg, ". Got ", dat$r_stress_truth, ".")
    )

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$e must be positive finite; dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("dat$lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$kappa)) || any(dat$kappa <= 0) ||
        any(!is.finite(dat$kappa_reference)) || any(dat$kappa_reference <= 0)) {
        stop("dat$kappa and dat$kappa_reference must be positive and finite.", call. = FALSE)
    }
    assert_true_s4c_oracle(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident, ".")
    )

    if (exists("validate_s4c_overdispersion_data", mode = "function")) {
        validate_s4c_overdispersion_data(dat)
    }

    y_stats <- count_stats_s4c_oracle(dat$y_coarse, prefix = "")
    y_ref_stats <- count_stats_s4c_oracle(dat$y_coarse_reference, prefix = "reference_")
    k_stats <- kappa_stats_s4c_oracle(dat$kappa, prefix = "")
    k_ref_stats <- kappa_stats_s4c_oracle(dat$kappa_reference, prefix = "reference_")

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat[[oracle_lambda_field]], finite = TRUE)

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type,
        r_reference_truth = dat$r_reference_truth,
        r_stress_truth = dat$r_stress_truth,
        r_ratio_to_reference = dat$r_stress_truth / dat$r_reference_truth,
        y_stats = y_stats,
        y_ref_stats = y_ref_stats,
        k_stats = k_stats,
        k_ref_stats = k_ref_stats,
        count_cv_increase = y_stats$count_cv - y_ref_stats$reference_count_cv,
        variance_to_mean_increase = y_stats$variance_to_mean - y_ref_stats$reference_variance_to_mean,
        kappa_cv_increase = k_stats$kappa_cv - k_ref_stats$reference_kappa_cv,
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_oracle_min = lambda_ident_range[[1L]],
        lambda_oracle_max = lambda_ident_range[[2L]]
    )
}

make_s4c_oracle_source_manifest <- function(reps,
                                            root_arg = root_dir,
                                            data_dir_arg = data_dir,
                                            data_scenario_id_arg = data_scenario_id) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s4c_oracle_source_data_path(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
        chk <- check_s4c_oracle_source_dataset(data_file)
        data.frame(
            scenario_id = data_scenario_id_arg,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            r_reference_truth = chk$r_reference_truth,
            r_stress_truth = chk$r_stress_truth,
            r_ratio_to_reference = chk$r_ratio_to_reference,
            mean_count = chk$y_stats$mean_count,
            reference_mean_count = chk$y_ref_stats$reference_mean_count,
            median_count = chk$y_stats$median_count,
            zero_prop = chk$y_stats$zero_prop,
            reference_zero_prop = chk$y_ref_stats$reference_zero_prop,
            total_count = chk$y_stats$total_count,
            max_count = chk$y_stats$max_count,
            reference_max_count = chk$y_ref_stats$reference_max_count,
            count_sd = chk$y_stats$count_sd,
            count_cv = chk$y_stats$count_cv,
            reference_count_cv = chk$y_ref_stats$reference_count_cv,
            count_cv_increase = chk$count_cv_increase,
            variance_to_mean = chk$y_stats$variance_to_mean,
            reference_variance_to_mean = chk$y_ref_stats$reference_variance_to_mean,
            variance_to_mean_increase = chk$variance_to_mean_increase,
            q95_count = chk$y_stats$q95_count,
            q99_count = chk$y_stats$q99_count,
            kappa_cv = chk$k_stats$kappa_cv,
            reference_kappa_cv = chk$k_ref_stats$reference_kappa_cv,
            kappa_cv_increase = chk$kappa_cv_increase,
            beta0_star_ident = chk$beta0_star_ident,
            lambda_raw_min = chk$lambda_raw_min,
            lambda_raw_max = chk$lambda_raw_max,
            lambda_oracle_min = chk$lambda_oracle_min,
            lambda_oracle_max = chk$lambda_oracle_max,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

## ---- source Scenario 3/S4C code -------------------------------------------
assert_file_exists_s4c_oracle(file.path(root_dir, s3_core_file), "Scenario 3 core script")
assert_file_exists_s4c_oracle(file.path(root_dir, s4c_data_generator_file), "Scenario 4C data-generation script")
source(file.path(root_dir, s4c_data_generator_file))
source_s4c_small_r_overdispersion(root = root_dir, verbose = FALSE)

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
        log_xi_j <- pmin(pmax(log_xi_j, s4c_oracle_log_xi_lower), s4c_oracle_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("S4C oracle safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

.s4c_oracle_guard_env <- new.env(parent = emptyenv())
reset_s4c_oracle_guards <- function() {
    .s4c_oracle_guard_env$n_beta_guard <- 0L
    .s4c_oracle_guard_env$n_kappa_guard <- 0L
    invisible(TRUE)
}
get_s4c_oracle_guards <- function() {
    list(
        n_beta_guard = .s4c_oracle_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4c_oracle_guard_env$n_kappa_guard %||% 0L
    )
}
reset_s4c_oracle_guards()

.update_beta_s3_original <- update_beta
update_beta <- function(beta_current, ...) {
    res <- .update_beta_s3_original(beta_current = beta_current, ...)
    smp <- res$sample
    bad <- FALSE
    if (length(smp) < 3L || any(!is.finite(smp))) {
        bad <- TRUE
    } else {
        bad <- smp[1] < s4c_oracle_beta0_bounds[1] || smp[1] > s4c_oracle_beta0_bounds[2] ||
            any(smp[-1] < s4c_oracle_beta_bounds[1] | smp[-1] > s4c_oracle_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s4c_oracle_guard_env$n_beta_guard <- (.s4c_oracle_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s4c_oracle_guard_rejected <- TRUE
    } else {
        res$s4c_oracle_guard_rejected <- FALSE
    }
    res
}

.update_kappa_s3_original <- update_kappa
update_kappa <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
    y <- as.matrix(y_coarse)
    L <- as.matrix(lambda_tilde)
    X <- as.matrix(xi)
    if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
        stop("S4C oracle safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
    }
    r_vec <- as.numeric(r)
    if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
    if (length(r_vec) != ncol(y)) {
        stop("S4C oracle safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
    }
    R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
    shape <- y + R
    rate <- X * L + R
    guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
    if (any(guard)) {
        .s4c_oracle_guard_env$n_kappa_guard <- (.s4c_oracle_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s4c_oracle_guard_env$n_kappa_guard <- (.s4c_oracle_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s4c_oracle_kappa_bounds[1]), s4c_oracle_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s4c_oracle_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

cat(sprintf(
    "Using S4C oracle-lambda+phi diagnostic guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s4c_oracle_beta0_bounds[1], s4c_oracle_beta0_bounds[2],
    s4c_oracle_beta_bounds[1], s4c_oracle_beta_bounds[2],
    s4c_oracle_log_xi_lower, s4c_oracle_log_xi_upper
))


## ---- oracle-lambda+phi MCMC -------------------------------------------------
run_s4c_oracle_lambda_phi_fixed_gamma_mcmc <- function(dat,
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
        phi_accept_trace = rep(NA, n_iter_local),
        r_accept_trace = matrix(FALSE, nrow = n_iter_local, ncol = n1_local),
        gamma_accept_trace = rep(NA, n_iter_local),
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

    r_accept_window <- rep(0L, n1_local)
    store_idx <- 0L
    start_time <- proc.time()

    for (iter in seq_len(n_iter_local)) {
        ## 1. beta update using marginal NB likelihood, conditional on oracle lambda and oracle phi.
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = oracle_lambda,
            phi = oracle_phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r,
            use_preconditioned = settings$use_preconditioned_beta %||% TRUE
        )
        state$beta0 <- beta_result$sample[1]
        state$beta <- beta_result$sample[2:3]

        ## 2. lambda, phi, tau_phi, and gamma are fixed in this diagnostic.
        state$lambda_tilde <- oracle_lambda
        state$phi <- oracle_phi
        state$gamma <- gamma_init
        state$gamma_common <- gamma_value
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
            phi = oracle_phi,
            priors = priors,
            mh_sd = state$r_proposal_sd,
            method = "marginal_nb",
            return_diag = TRUE
        )
        state$r <- r_result$r

        ## 4. kappa update for diagnostic Poisson-augmentation likelihood only.
        kappa_result <- update_kappa(
            y_coarse = dat$y_coarse,
            lambda_tilde = oracle_lambda,
            xi = state$xi,
            r = state$r,
            return_diag = TRUE
        )
        state$kappa <- kappa_result$kappa

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
            phi = oracle_phi,
            r = state$r,
            lambda_tilde = oracle_lambda
        )

        diagnostics$loglik_trace[iter] <- loglik
        diagnostics$loglik_nb_trace[iter] <- loglik_nb
        diagnostics$beta_n_reject[iter] <- beta_result$n_reject %||% NA_real_
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
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
            r_rate <- mean(diagnostics$r_accept_trace[i0:iter, , drop = FALSE])
            beta_rej <- mean(diagnostics$beta_n_reject[i0:iter], na.rm = TRUE)
            cat(sprintf(
                paste0(
                    "iter %5d/%d [%.0fs] loglik_nb=%.1f beta0=%.3f ",
                    "beta=(%.3f, %.3f) r_mean=%.2f oracle_lambda=[%.3f, %.3f] ",
                    "phi=fixed gamma=%.3f | r_acc=%.2f beta_rej=%.1f\n"
                ),
                iter, n_iter_local, elapsed, loglik_nb, state$beta0, state$beta[1],
                state$beta[2], mean(state$r), min(oracle_lambda), max(oracle_lambda),
                state$gamma_common, r_rate, beta_rej
            ))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]
    diagnostics$elapsed_sec <- elapsed_total
    diagnostics$phi_accept_rate <- NA_real_
    diagnostics$r_accept_rate <- colMeans(diagnostics$r_accept_trace)
    diagnostics$gamma_accept_rate <- NA_real_
    diagnostics$beta_mean_n_reject <- mean(diagnostics$beta_n_reject, na.rm = TRUE)
    diagnostics$phi_proposal_sd_final <- NA_real_
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
            model = "S4C oracle-lambda-phi fixed-gamma diagnostic",
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
            uses_marginal_nb_for_beta_r = TRUE,
            uses_kappa_for_diagnostic_poisson_aug_loglik = TRUE,
            uses_recentered_identified_parameterization = TRUE
        ),
        n_stored = n_stored
    )
}

## ---- fit wrapper ------------------------------------------------------------
fit_s4c_oracle_lambda_phi_one_rep <- function(rep_id,
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
    dat_file <- s4c_oracle_source_data_path(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
    chk <- check_s4c_oracle_source_dataset(dat_file)
    dat <- chk$dat

    validate_s3_data(dat)
    if (exists("validate_s4c_overdispersion_data", mode = "function")) {
        validate_s4c_overdispersion_data(dat)
    }

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) settings[[nm]] <- settings_override[[nm]]
    }

    oracle_lambda <- dat[[oracle_lambda_field]]
    oracle_phi <- as.numeric(dat[[oracle_phi_field]])
    if (length(oracle_phi) != dat$n1 || any(!is.finite(oracle_phi))) {
        stop("dat$", oracle_phi_field, " must be finite and have length n1.", call. = FALSE)
    }

    cat(sprintf("=== S4C oracle-lambda+phi fixed-gamma diagnostic: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Stress    : %s; r truth = %.3f (reference %.3f)\n",
                dat$stress_type, dat$r_stress_truth, dat$r_reference_truth))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f, count CV = %.3f\n",
                chk$y_stats$mean_count, chk$y_stats$zero_prop, chk$y_stats$count_cv))
    cat(sprintf("Oracle    : lambda = dat$%s in [%.3f, %.3f]; phi = dat$%s in [%.3f, %.3f]\n",
                oracle_lambda_field, min(oracle_lambda), max(oracle_lambda),
                oracle_phi_field, min(oracle_phi), max(oracle_phi)))
    cat(sprintf("Fixed     : gamma = %.3f\n", fixed_gamma_value_arg))
    cat("Estimated : beta0, beta1, beta2, r, kappa\n")
    cat("Disabled  : lambda, phi, tau_phi, gamma, delta, omega updates\n\n")

    reset_s4c_oracle_guards()
    fit <- run_s4c_oracle_lambda_phi_fixed_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        oracle_lambda = oracle_lambda,
        oracle_phi = oracle_phi,
        gamma_value = fixed_gamma_value_arg,
        verbose = verbose_arg
    )

    guard_counts <- get_s4c_oracle_guards()
    fit$diagnostics$s4c_oracle_numeric_guards <- guard_counts
    fit$diagnostics$s4c_oracle_beta_guard_count <- guard_counts$n_beta_guard
    fit$diagnostics$s4c_oracle_kappa_guard_count <- guard_counts$n_kappa_guard
    fit$diagnostics$phi_accept_rate <- NA_real_
    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$phi_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_

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
    summary$oracle_lambda_source <- oracle_lambda_field
    summary$oracle_phi_fixed_in_fit <- TRUE
    summary$oracle_phi_source <- oracle_phi_field
    summary$phi_sampled_in_fit <- FALSE
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
    summary$phi_accept_rate <- NA_real_
    summary$phi_proposal_sd_final <- NA_real_

    ## S4C-specific data context.
    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$r_reference_truth <- dat$r_reference_truth %||% NA_real_
    summary$r_stress_truth <- dat$r_stress_truth %||% mean(dat$r_star %||% NA_real_)
    summary$r_ratio_to_reference <- summary$r_stress_truth / summary$r_reference_truth
    summary$reference_mean_count <- chk$y_ref_stats$reference_mean_count
    summary$reference_zero_prop <- chk$y_ref_stats$reference_zero_prop
    summary$observed_mean_count <- chk$y_stats$mean_count
    summary$observed_zero_prop <- chk$y_stats$zero_prop
    summary$observed_total_count <- chk$y_stats$total_count
    summary$observed_max_count <- chk$y_stats$max_count
    summary$count_sd <- chk$y_stats$count_sd
    summary$count_cv <- chk$y_stats$count_cv
    summary$reference_count_cv <- chk$y_ref_stats$reference_count_cv
    summary$count_cv_increase <- chk$count_cv_increase
    summary$variance_to_mean <- chk$y_stats$variance_to_mean
    summary$reference_variance_to_mean <- chk$y_ref_stats$reference_variance_to_mean
    summary$variance_to_mean_increase <- chk$variance_to_mean_increase
    summary$q95_count <- chk$y_stats$q95_count
    summary$q99_count <- chk$y_stats$q99_count
    summary$reference_max_count <- chk$y_ref_stats$reference_max_count
    summary$kappa_truth_cv <- chk$k_stats$kappa_cv
    summary$reference_kappa_cv <- chk$k_ref_stats$reference_kappa_cv
    summary$kappa_cv_increase <- chk$kappa_cv_increase
    summary$lambda_oracle_min <- chk$lambda_oracle_min
    summary$lambda_oracle_max <- chk$lambda_oracle_max
    summary$phi_oracle_min <- min(oracle_phi)
    summary$phi_oracle_max <- max(oracle_phi)

    r_extra <- r_recovery_s4c_oracle(fit$samples$r, dat$r_star)
    for (nm in names(r_extra)) summary[[nm]] <- r_extra[[nm]]

    kappa_extra <- summarise_s4c_oracle_kappa_recovery(fit$samples$kappa, dat$kappa)
    for (nm in names(kappa_extra)) summary[[nm]] <- kappa_extra[[nm]]

    summary$s4c_oracle_beta_guard_count <- guard_counts$n_beta_guard
    summary$s4c_oracle_kappa_guard_count <- guard_counts$n_kappa_guard
    summary$s4c_oracle_lambda_input_guard_count <- 0L
    summary$s4c_oracle_lambda_output_guard_count <- 0L
    summary$s4c_oracle_phi_fixed <- TRUE

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id_arg,
            data_scenario_id = data_scenario_id_arg,
            stress_type = dat$stress_type,
            r_reference_truth = dat$r_reference_truth,
            r_stress_truth = dat$r_stress_truth,
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir_arg,
            fit_file_prefix = "fit_S4C_oracle_lambda_phi_fixed_gamma_rep",
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root_arg, output_dir_arg, scenario_id_arg)
        ensure_dir_s4c_oracle(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4C_oracle_lambda_phi_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4C_oracle_lambda_phi_fixed_gamma_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s4c_oracle_lambda_phi_batch <- function(reps = reps_formal,
                                            scenario_id_arg = scenario_id,
                                            data_scenario_id_arg = data_scenario_id,
                                            data_dir_arg = data_dir,
                                            output_dir_arg = output_dir,
                                            root_arg = root_dir,
                                            settings_override = list(),
                                            fixed_gamma_value_arg = fixed_gamma_value,
                                            verbose_arg = verbose,
                                            overwrite_existing = overwrite_fit) {
    out_dir <- file.path(root_arg, output_dir_arg, scenario_id_arg)
    ensure_dir_s4c_oracle(out_dir)

    summaries <- list()
    fit_manifest <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4C_oracle_lambda_phi_fixed_gamma_rep", rr, ".rds"))
        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
        } else {
            fit <- fit_s4c_oracle_lambda_phi_one_rep(
                rep_id = rep_id,
                scenario_id_arg = scenario_id_arg,
                data_scenario_id_arg = data_scenario_id_arg,
                data_dir_arg = data_dir_arg,
                output_dir_arg = output_dir_arg,
                root_arg = root_arg,
                settings_override = settings_override,
                fixed_gamma_value_arg = fixed_gamma_value_arg,
                verbose_arg = verbose_arg,
                save_result = TRUE,
                return_result = TRUE
            )
            summaries[[rr]] <- fit$summary
        }
        data_file <- s4c_oracle_source_data_path(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
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

    summary_file <- file.path(out_dir, "summary_S4C_oracle_lambda_phi_fixed_gamma_all_reps.csv")
    manifest_file <- file.path(out_dir, "s4c_oracle_lambda_phi_fit_manifest.csv")
    run_info_file <- file.path(out_dir, "run_info_S4C_oracle_lambda_phi_fixed_gamma_diagnostic_T100.rds")

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
            oracle_phi_field = oracle_phi_field,
            created_at = Sys.time(),
            output_dir = out_dir
        ),
        run_info_file
    )

    cat("\n=== S4C oracle-lambda+phi diagnostic summary ===\n")
    print(summary_all[, intersect(c(
        "rep_id", "mean_count", "zero_prop", "count_cv", "beta0_mean", "beta1_mean", "beta2_mean",
        "r_true_mean", "r_mean", "r_region_coverage_95", "phi_rmse", "phi_cor",
        "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "kappa_truth_cv", "kappa_post_mean_cv", "kappa_rmse", "log_kappa_rmse", "cor_log_kappa",
        "s4c_oracle_beta_guard_count", "s4c_oracle_kappa_guard_count"
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

s4c_oracle_lambda_phi_summary <- fit_s4c_oracle_lambda_phi_batch(
    reps = reps_formal,
    scenario_id_arg = scenario_id,
    data_scenario_id_arg = data_scenario_id,
    data_dir_arg = data_dir,
    output_dir_arg = output_dir,
    root_arg = root_dir,
    settings_override = settings_override,
    fixed_gamma_value_arg = fixed_gamma_value,
    verbose_arg = verbose,
    overwrite_existing = overwrite_fit
)

cat("\nS4C oracle-lambda+phi fixed-gamma diagnostic finished successfully.\n")
