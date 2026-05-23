## ============================================================================
## s1_baseline_only_clean.R
##
## Clean Scenario 1 script for the MSSTNB project.
##
## Scenario 1: baseline-only negative-binomial spatial regression.
##
## Purpose:
##   Isolate recovery of beta, phi, tau_phi, and r when the residual dynamic
##   risk is fixed at one.  This script intentionally does NOT update lambda,
##   gamma, delta, omega, or kappa during fitting.
##
## Fitted model:
##   Y_tj ~ NB(mu_tj, r_j)
##   log(mu_tj) = log(e_tj) + beta0 + beta1 x1_tj + beta2 x2_tj + phi_j
##   phi = B_ICAR u, sum_j phi_j = 0
##
## Important functions:
##   source_s1_baseline_only()
##   simulate_s1_baseline_only_one()
##   simulate_s1_baseline_only_batch()
##   run_s1_baseline_mcmc()
##   fit_s1_baseline_one_rep()
##   fit_s1_baseline_batch()
##   sanity_check_s1_baseline_only()
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.same_dim_s1 <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) && all(as.integer(d) == as.integer(target_dim))
}

.require_file_s1 <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked_s1 <- function(path, verbose = TRUE) {
    .require_file_s1(path)
    if (isTRUE(verbose)) {
        message("source: ", path)
    }
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source project dependencies
## -----------------------------------------------------------------------------
source_s1_baseline_only <- function(root = ".", verbose = TRUE) {
    r_dir <- file.path(root, "R")
    mcmc_dir <- file.path(r_dir, "mcmc")

    files <- c(
        file.path(r_dir, "00_setup.R"),
        file.path(r_dir, "01_helpers.R"),
        file.path(r_dir, "02_simulate_data.R"),
        file.path(mcmc_dir, "mcmc_config.R"),
        file.path(mcmc_dir, "mcmc_utils.R"),
        file.path(mcmc_dir, "update_regression_revised.R"),
        file.path(mcmc_dir, "update_icar_revised.R"),
        file.path(mcmc_dir, "update_dispersion_revised.R")
    )

    invisible(lapply(files, .source_checked_s1, verbose = verbose))

    needed <- c(
        "MCMC_SETTINGS", "MCMC_PRIORS",
        "TT", "N1", "N_CHILDREN", "C0", "H", "B_ICAR",
        "REP_SEEDS", "generate_inputs", "generate_icar", "rdirichlet",
        "fit_glm_for_ess", "update_beta", "update_phi",
        "update_tau_phi", "update_r"
    )

    missing <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop("After sourcing dependencies, missing objects: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Shared validators and utilities
## -----------------------------------------------------------------------------
validate_s1_data <- function(dat) {
    required <- c("y_coarse", "e", "x1", "x2", "TT", "n1")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (!is.matrix(dat$y_coarse)) {
        stop("dat$y_coarse must be a matrix.", call. = FALSE)
    }

    y_dim <- dim(dat$y_coarse)
    for (nm in c("e", "x1", "x2")) {
        if (!is.matrix(dat[[nm]]) || !identical(dim(dat[[nm]]), y_dim)) {
            stop("dat$", nm, " must be a matrix with the same dimension as dat$y_coarse.",
                 call. = FALSE)
        }
    }

    if (!identical(as.integer(dat$TT), as.integer(y_dim[1])) ||
        !identical(as.integer(dat$n1), as.integer(y_dim[2]))) {
        stop("dat$TT or dat$n1 is inconsistent with dat$y_coarse.",
             call. = FALSE)
    }

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.",
             call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0)) {
        stop("dat$e must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }

    invisible(TRUE)
}

build_s1_spatial <- function(H_obj = H, B_obj = B_ICAR) {
    if (!is.matrix(H_obj) || !is.matrix(B_obj)) {
        stop("H_obj and B_obj must be matrices.", call. = FALSE)
    }
    list(
        H = H_obj,
        B_ICAR = B_obj,
        BHB = crossprod(B_obj, H_obj %*% B_obj)
    )
}

validate_s1_fit_objects <- function(dat, priors, spatial) {
    validate_s1_data(dat)

    n1 <- dat$n1
    if (!is.matrix(spatial$H) || .same_dim_s1(spatial$H, c(n1, n1)) == FALSE) {
        stop("spatial$H must be an n1 by n1 matrix.", call. = FALSE)
    }
    if (!is.matrix(spatial$B_ICAR) ||
        .same_dim_s1(spatial$B_ICAR, c(n1, n1 - 1L)) == FALSE) {
        stop("spatial$B_ICAR must be an n1 by n1 - 1 matrix.", call. = FALSE)
    }
    if (!is.matrix(spatial$BHB) ||
        .same_dim_s1(spatial$BHB, c(n1 - 1L, n1 - 1L)) == FALSE) {
        stop("spatial$BHB must be an n1 - 1 by n1 - 1 matrix.", call. = FALSE)
    }

    required_priors <- c(
        "beta0_mean", "beta0_sd", "beta_mean", "beta_sd",
        "tau_phi_shape", "tau_phi_rate", "r_shape", "r_rate"
    )
    missing <- setdiff(required_priors, names(priors))
    if (length(missing) > 0L) {
        stop("priors is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    if (length(priors$beta_mean) != 2L || length(priors$beta_sd) != 2L) {
        stop("This Scenario 1 script expects exactly two covariates.",
             call. = FALSE)
    }

    invisible(TRUE)
}

make_s1_file_name <- function(prefix, rep_id, ext = "rds") {
    paste0(prefix, sprintf("%02d", as.integer(rep_id)), ".", ext)
}

## -----------------------------------------------------------------------------
## Scenario 1 simulation
## -----------------------------------------------------------------------------

## Internal input generator for formal Scenario 1 validation.
## This avoids accidental non-finite Poisson rates when T is increased.
generate_s1_baseline_inputs <- function(TT_use,
                                        n1_use,
                                        exposure_min = 300,
                                        exposure_max = 1200,
                                        exposure_jitter_sd = 0.02,
                                        x1_ar = 0.70,
                                        x1_innov_sd = 0.70,
                                        x2_mode = c("binary_region", "continuous_region", "continuous_time"),
                                        x2_ar = 0.50,
                                        x2_innov_sd = 0.80,
                                        standardize_covariates = TRUE) {
    x2_mode <- match.arg(x2_mode)

    e_base <- stats::runif(n1_use, min = exposure_min, max = exposure_max)
    e <- matrix(NA_real_, TT_use, n1_use)
    for (t in seq_len(TT_use)) {
        e[t, ] <- e_base * exp(stats::rnorm(n1_use, mean = 0, sd = exposure_jitter_sd))
    }

    x1_raw <- matrix(NA_real_, TT_use, n1_use)
    x1_raw[1, ] <- stats::rnorm(n1_use, mean = 0, sd = 1)
    if (TT_use >= 2L) {
        for (t in 2:TT_use) {
            x1_raw[t, ] <- x1_ar * x1_raw[t - 1, ] +
                stats::rnorm(n1_use, mean = 0, sd = x1_innov_sd)
        }
    }

    if (x2_mode == "binary_region") {
        x2_vec <- stats::rbinom(n1_use, size = 1, prob = 0.4)
        ## Make sure both groups are represented when possible.
        if (n1_use >= 2L && length(unique(x2_vec)) == 1L) {
            x2_vec[seq_len(floor(n1_use / 2))] <- 0L
            x2_vec[(floor(n1_use / 2) + 1L):n1_use] <- 1L
        }
        x2_raw <- matrix(rep(x2_vec, each = TT_use), nrow = TT_use, ncol = n1_use)
    } else if (x2_mode == "continuous_region") {
        x2_vec <- stats::rnorm(n1_use, mean = 0, sd = 1)
        x2_raw <- matrix(rep(x2_vec, each = TT_use), nrow = TT_use, ncol = n1_use)
    } else {
        x2_raw <- matrix(NA_real_, TT_use, n1_use)
        x2_raw[1, ] <- stats::rnorm(n1_use, mean = 0, sd = 1)
        if (TT_use >= 2L) {
            for (t in 2:TT_use) {
                x2_raw[t, ] <- x2_ar * x2_raw[t - 1, ] +
                    stats::rnorm(n1_use, mean = 0, sd = x2_innov_sd)
            }
        }
    }

    if (isTRUE(standardize_covariates)) {
        x1 <- standardize_matrix(x1_raw)
        x2 <- standardize_matrix(x2_raw)
    } else {
        x1 <- x1_raw
        x2 <- x2_raw
    }

    out <- list(
        e = e,
        x1 = x1,
        x2 = x2,
        x1_raw = x1_raw,
        x2_raw = x2_raw,
        x1_mean = mean(as.numeric(x1_raw)),
        x1_sd = stats::sd(as.numeric(x1_raw)),
        x2_mean = mean(as.numeric(x2_raw)),
        x2_sd = stats::sd(as.numeric(x2_raw)),
        x2_mode = x2_mode
    )

    for (nm in c("e", "x1", "x2")) {
        if (any(!is.finite(out[[nm]]))) {
            stop("Non-finite values generated in Scenario 1 input matrix: ", nm, call. = FALSE)
        }
    }

    out
}

## -----------------------------------------------------------------------------
## Scenario 1 simulation
## -----------------------------------------------------------------------------
simulate_s1_baseline_only_one <- function(seed = 1L,
                                          TT_use = TT,
                                          n1_use = N1,
                                          n_children_use = N_CHILDREN,
                                          c0 = C0,
                                          beta0_truth = -1.5,
                                          beta_truth = c(0.5, -0.4),
                                          r_truth = 15,
                                          tau_phi_truth = 2,
                                          gamma_truth = 0.8,
                                          delta_truth = 0.9,
                                          residual_risk_value = 1,
                                          omega_mode = c("fixed_prior_mean", "dirichlet_static"),
                                          H_obj = H,
                                          B_obj = B_ICAR,
                                          scenario_id = "S1_BASELINE_ONLY",
                                          rep_id = NA_integer_,
                                          use_internal_inputs = TRUE,
                                          exposure_min = 300,
                                          exposure_max = 1200,
                                          x1_ar = 0.70,
                                          x1_innov_sd = 0.70,
                                          x2_mode = c("binary_region", "continuous_region", "continuous_time"),
                                          x2_ar = 0.50,
                                          x2_innov_sd = 0.70,
                                          standardize_covariates = TRUE,
                                          max_poisson_rate = 1e7) {
    omega_mode <- match.arg(omega_mode)
    x2_mode <- match.arg(x2_mode)

    ## Force dimensions to be integer-valued. Several existing update functions
    ## use identical(dim(H), c(n1, n1)); this fails if n1 is numeric rather
    ## than integer even when the numerical values are the same.
    TT_use <- as.integer(TT_use)
    n1_use <- as.integer(n1_use)
    n_children_use <- as.integer(n_children_use)

    set.seed(seed)

    if (length(beta_truth) != 2L) {
        stop("beta_truth must have length 2.", call. = FALSE)
    }
    if (length(r_truth) == 1L) {
        r_truth <- rep(r_truth, n1_use)
    }
    if (length(gamma_truth) == 1L) {
        gamma_truth <- rep(gamma_truth, n1_use)
    }
    if (length(r_truth) != n1_use) {
        stop("r_truth must be scalar or length n1_use.", call. = FALSE)
    }

    if (isTRUE(use_internal_inputs)) {
        inputs <- generate_s1_baseline_inputs(
            TT_use = TT_use,
            n1_use = n1_use,
            exposure_min = exposure_min,
            exposure_max = exposure_max,
            x1_ar = x1_ar,
            x1_innov_sd = x1_innov_sd,
            x2_mode = x2_mode,
            x2_ar = x2_ar,
            x2_innov_sd = x2_innov_sd,
            standardize_covariates = standardize_covariates
        )
    } else {
        inputs <- generate_inputs(TT = TT_use, n1 = n1_use)
    }

    e <- inputs$e
    x1 <- inputs$x1
    x2 <- inputs$x2

    ## Defensive dimension checks. This is important when using T = 100.
    ## Use all(dim(.) == .) rather than identical(), because dim() is integer
    ## while c(TT_use, n1_use) may be numeric when passed by the user.
    target_dim <- c(as.integer(TT_use), as.integer(n1_use))
    if (!.same_dim_s1(e, target_dim) ||
        !.same_dim_s1(x1, target_dim) ||
        !.same_dim_s1(x2, target_dim)) {
        stop(sprintf(
            paste0(
                "Scenario 1 inputs do not have dimensions TT_use x n1_use. ",
                "Expected %d x %d; got e=%s, x1=%s, x2=%s."
            ),
            as.integer(TT_use), as.integer(n1_use),
            paste(dim(e), collapse = " x "),
            paste(dim(x1), collapse = " x "),
            paste(dim(x2), collapse = " x ")
        ), call. = FALSE)
    }

    ## The ICAR objects are created by R/00_setup.R. In the current project they
    ## may be fixed at N1 = 9. If n1_use is changed, H and B_ICAR must be changed
    ## consistently. Otherwise phi has the wrong dimension.
    if (!.same_dim_s1(H_obj, c(n1_use, n1_use)) ||
        !.same_dim_s1(B_obj, c(n1_use, n1_use - 1L))) {
        stop(sprintf(
            paste0(
                "The supplied ICAR objects are not compatible with n1_use. ",
                "n1_use=%d, dim(H)=%s, dim(B_ICAR)=%s. ",
                "Use n1_use = N1 from R/00_setup.R, or rebuild H and B_ICAR for the new n1_use."
            ),
            as.integer(n1_use),
            paste(dim(H_obj), collapse = " x "),
            paste(dim(B_obj), collapse = " x ")
        ), call. = FALSE)
    }
    if (any(!is.finite(e)) || any(!is.finite(x1)) || any(!is.finite(x2))) {
        stop("Scenario 1 inputs contain non-finite values before count generation.", call. = FALSE)
    }

    phi_truth <- generate_icar(tau_phi = tau_phi_truth, H = H_obj, B = B_obj)
    phi_truth <- phi_truth - mean(phi_truth)
    if (length(phi_truth) != n1_use || any(!is.finite(phi_truth))) {
        stop("Generated phi_truth is non-finite or has the wrong length.", call. = FALSE)
    }

    lambda_tilde <- matrix(residual_risk_value, nrow = TT_use, ncol = n1_use)
    lambda_tilde0 <- rep(residual_risk_value, n1_use)

    xi <- matrix(NA_real_, TT_use, n1_use)
    mu_nb <- matrix(NA_real_, TT_use, n1_use)
    kappa <- matrix(NA_real_, TT_use, n1_use)
    poisson_rate <- matrix(NA_real_, TT_use, n1_use)
    y_coarse <- matrix(NA_integer_, TT_use, n1_use)

    for (j in seq_len(n1_use)) {
        eta_j <- beta0_truth + beta_truth[1] * x1[, j] +
            beta_truth[2] * x2[, j] + phi_truth[j]
        xi[, j] <- e[, j] * exp(eta_j)
        mu_nb[, j] <- xi[, j] * residual_risk_value
        kappa[, j] <- stats::rgamma(TT_use, shape = r_truth[j], rate = r_truth[j])
        poisson_rate[, j] <- mu_nb[, j] * kappa[, j]
    }

    bad_rate <- !is.finite(poisson_rate) | poisson_rate < 0 | poisson_rate > max_poisson_rate
    if (any(bad_rate)) {
        bad_idx <- which(bad_rate, arr.ind = TRUE)[1, ]
        stop(sprintf(
            paste0(
                "Non-finite, negative, or excessively large Poisson rate in Scenario 1 DGP. ",
                "First bad cell: t=%d, j=%d, rate=%s, e=%s, x1=%s, x2=%s, phi=%s. ",
                "Try smaller exposure_max, less variable covariates, or a smaller beta0_truth."
            ),
            bad_idx[1], bad_idx[2],
            as.character(poisson_rate[bad_idx[1], bad_idx[2]]),
            as.character(e[bad_idx[1], bad_idx[2]]),
            as.character(x1[bad_idx[1], bad_idx[2]]),
            as.character(x2[bad_idx[1], bad_idx[2]]),
            as.character(phi_truth[bad_idx[2]])
        ), call. = FALSE)
    }

    y_coarse[,] <- stats::rpois(length(poisson_rate), lambda = as.numeric(poisson_rate))
    dim(y_coarse) <- c(TT_use, n1_use)

    if (anyNA(y_coarse)) {
        stop("rpois produced NA values in Scenario 1 DGP after rate validation.", call. = FALSE)
    }

    omega <- array(NA_real_, dim = c(TT_use, n1_use, n_children_use))
    y_fine <- array(0L, dim = c(TT_use, n1_use, n_children_use))
    omega_prior_mean <- c0 / sum(c0)

    for (t in seq_len(TT_use)) {
        for (j in seq_len(n1_use)) {
            if (omega_mode == "fixed_prior_mean") {
                omega[t, j, ] <- omega_prior_mean
            } else {
                omega[t, j, ] <- as.numeric(rdirichlet(1, c0))
            }

            if (y_coarse[t, j] > 0L) {
                y_fine[t, j, ] <- as.integer(stats::rmultinom(1, size = y_coarse[t, j],
                                                               prob = omega[t, j, ]))
            }
        }
    }

    coherent <- all(y_coarse == apply(y_fine, c(1, 2), sum))
    if (!coherent) {
        stop("Scenario 1 generated fine counts are not coherent with coarse counts.", call. = FALSE)
    }

    out <- list(
        scenario_id = scenario_id,
        rep_id = rep_id,
        y_coarse = y_coarse,
        y_fine = y_fine,
        y_levels = list(y_coarse, y_fine),
        e = e,
        x1 = x1,
        x2 = x2,
        X = array(c(x1, x2), dim = c(TT_use, n1_use, 2L)),
        x = array(c(x1, x2), dim = c(TT_use, n1_use, 2L)),
        tree = list(n1 = n1_use, n_children = n_children_use),
        H = H_obj,
        B_ICAR = B_obj,
        lambda_tilde = lambda_tilde,
        lambda_tilde0 = lambda_tilde0,
        omega = omega,
        xi = xi,
        mu_nb = mu_nb,
        kappa = kappa,
        poisson_rate = poisson_rate,
        beta0_star = beta0_truth,
        beta_star = beta_truth,
        phi_star = phi_truth,
        tau_phi_star = tau_phi_truth,
        r_star = r_truth,
        gamma_star = gamma_truth,
        delta_star = delta_truth,
        beta0_star_ident = beta0_truth,
        beta_star_ident = beta_truth,
        phi_star_ident = phi_truth,
        lambda_tilde_ident = lambda_tilde,
        residual_risk_value = residual_risk_value,
        omega_mode = omega_mode,
        data_type = "baseline_only_lambda_fixed",
        TT = TT_use,
        n1 = n1_use,
        N1 = n1_use,
        n_children = n_children_use,
        N_CHILDREN = n_children_use,
        mean_count = mean(y_coarse),
        median_count = stats::median(as.numeric(y_coarse)),
        zero_prop = mean(y_coarse == 0),
        mean_exposure = mean(e),
        min_exposure = min(e),
        max_exposure = max(e),
        coherent = coherent,
        input_generator = if (isTRUE(use_internal_inputs)) "generate_s1_baseline_inputs" else "generate_inputs"
    )

    if (!is.null(inputs$x1_raw)) out$x1_raw <- inputs$x1_raw
    if (!is.null(inputs$x2_raw)) out$x2_raw <- inputs$x2_raw
    if (!is.null(inputs$x2_mode)) out$x2_mode <- inputs$x2_mode

    validate_s1_data(out)
    out
}

simulate_s1_baseline_only_batch <- function(reps = 1:10,
                                             data_dir = "data_revised",
                                             scenario_id = "S1_BASELINE_ONLY",
                                             root = ".",
                                             seed_base = NULL,
                                             overwrite_existing = TRUE,
                                             verbose = TRUE,
                                             ...) {
    out_dir <- file.path(root, data_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    manifest <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        out_file <- file.path(out_dir, paste0("data_rep", rr, ".rds"))

        if (file.exists(out_file) && !isTRUE(overwrite_existing)) {
            if (isTRUE(verbose)) {
                message("Skipping existing file: ", out_file)
            }
            dat <- readRDS(out_file)
        } else {
            seed <- if (!is.null(seed_base)) {
                as.integer(seed_base + rep_id)
            } else if (exists("REP_SEEDS", envir = .GlobalEnv) && rep_id <= length(REP_SEEDS)) {
                as.integer(REP_SEEDS[rep_id])
            } else {
                as.integer(2026000L + rep_id)
            }

            dat <- simulate_s1_baseline_only_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                message(sprintf(
                    "Saved %s | mean_count=%.2f zero_prop=%.3f",
                    out_file, dat$mean_count, dat$zero_prop
                ))
            }
        }

        manifest[[rr]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            file = out_file,
            TT = dat$TT,
            n1 = dat$n1,
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            lambda_min = min(dat$lambda_tilde),
            lambda_max = max(dat$lambda_tilde),
            coherent = isTRUE(dat$coherent),
            stringsAsFactors = FALSE
        )
    }

    manifest_df <- do.call(rbind, manifest)
    manifest_file <- file.path(out_dir, "manifest_S1_baseline_only.csv")
    write.csv(manifest_df, manifest_file, row.names = FALSE)
    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
    }

    invisible(manifest_df)
}

## -----------------------------------------------------------------------------
## Scenario 1 fitting
## -----------------------------------------------------------------------------
adapt_sd_s1 <- function(current_sd, n_accept, n_trials,
                        target_rate = 0.30,
                        increase = 1.10,
                        decrease = 0.90,
                        min_sd = 1e-4,
                        max_sd = 5.0) {
    if (length(n_accept) == 1L && length(current_sd) > 1L) {
        n_accept <- rep(n_accept, length(current_sd))
    }
    rate <- n_accept / n_trials
    out <- current_sd
    out[rate > target_rate + 0.05] <- out[rate > target_rate + 0.05] * increase
    out[rate < target_rate - 0.05] <- out[rate < target_rate - 0.05] * decrease
    pmin(pmax(out, min_sd), max_sd)
}

compute_s1_loglik_nb <- function(y_coarse, e, x1, x2, beta0, beta, phi, r) {
    n1 <- ncol(y_coarse)
    out <- 0
    for (j in seq_len(n1)) {
        log_mu_j <- log(e[, j]) + beta0 + beta[1] * x1[, j] +
            beta[2] * x2[, j] + phi[j]
        out <- out + .nb_or_poisson_loglik(y = y_coarse[, j],
                                           log_mu = log_mu_j,
                                           r = r[j])
    }
    out
}

compute_s1_mu <- function(e, x1, x2, beta0, beta, phi) {
    TT_use <- nrow(e)
    n1_use <- ncol(e)
    mu <- matrix(NA_real_, TT_use, n1_use)
    for (j in seq_len(n1_use)) {
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        eta_j <- pmin(pmax(eta_j, -700), 700)
        mu[, j] <- e[, j] * exp(eta_j)
    }
    mu
}

initialise_s1_state <- function(dat, settings, priors, spatial) {
    n1 <- dat$n1
    p <- length(priors$beta_mean)

    ess_base <- fit_glm_for_ess(dat)
    beta0 <- ess_base$center[1]
    beta <- ess_base$center[2:(p + 1L)]

    list(
        beta0 = beta0,
        beta = beta,
        phi = rep(0, n1),
        u = rep(0, n1 - 1L),
        tau_phi = settings$tau_phi_init %||% 1,
        r = rep(settings$r_init %||% 10, n1),
        ess_base = ess_base,
        phi_proposal_sd = rep(settings$phi_proposal_sd_init %||% 0.01, n1 - 1L),
        r_proposal_sd = rep(settings$r_proposal_sd_init %||% 0.50, n1)
    )
}

run_s1_baseline_mcmc <- function(dat,
                                 settings = MCMC_SETTINGS,
                                 priors = MCMC_PRIORS,
                                 spatial = build_s1_spatial(),
                                 verbose = 1000L) {
    validate_s1_fit_objects(dat, priors, spatial)

    n_iter <- as.integer(settings$n_iter %||% 20000L)
    n_burnin <- as.integer(settings$n_burnin %||% 10000L)
    n_thin <- as.integer(settings$n_thin %||% 5L)
    adapt_interval <- as.integer(settings$adapt_interval %||% 50L)

    if (n_iter <= 0L || n_burnin < 0L || n_burnin >= n_iter || n_thin <= 0L) {
        stop("Invalid MCMC iteration settings.", call. = FALSE)
    }

    n_stored <- as.integer(floor((n_iter - n_burnin) / n_thin))
    TT_use <- as.integer(dat$TT)
    n1 <- as.integer(dat$n1)
    p <- length(priors$beta_mean)

    ## This is the key Scenario 1 restriction.
    lambda_fixed <- matrix(1, nrow = TT_use, ncol = n1)

    samples <- list(
        beta0 = numeric(n_stored),
        beta = matrix(NA_real_, nrow = n_stored, ncol = p),
        phi = matrix(NA_real_, nrow = n_stored, ncol = n1),
        tau_phi = numeric(n_stored),
        r = matrix(NA_real_, nrow = n_stored, ncol = n1),
        loglik = numeric(n_stored)
    )

    diagnostics <- list(
        loglik_trace = rep(NA_real_, n_iter),
        beta_n_reject = rep(NA_real_, n_iter),
        phi_accept_trace = rep(FALSE, n_iter),
        r_accept_trace = matrix(FALSE, nrow = n_iter, ncol = n1),
        phi_log_alpha = rep(NA_real_, n_iter),
        r_log_alpha = matrix(NA_real_, nrow = n_iter, ncol = n1)
    )

    state <- initialise_s1_state(dat, settings, priors, spatial)
    phi_accept_window <- 0L
    r_accept_window <- rep(0L, n1)
    store_idx <- 0L
    start_time <- proc.time()

    for (iter in seq_len(n_iter)) {
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = lambda_fixed,
            phi = state$phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r,
            use_preconditioned = settings$use_preconditioned_beta %||% TRUE
        )
        state$beta0 <- beta_result$sample[1]
        state$beta <- beta_result$sample[2:3]

        phi_result <- update_phi(
            u_current = state$u,
            B = spatial$B_ICAR,
            BHB = spatial$BHB,
            tau_phi = state$tau_phi,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            lambda_tilde = lambda_fixed,
            r = state$r,
            proposal_sd = state$phi_proposal_sd
        )
        state$u <- phi_result$u
        state$phi <- phi_result$phi

        state$tau_phi <- update_tau_phi(
            phi = state$phi,
            H = spatial$H,
            n1 = n1,
            priors = priors
        )

        r_result <- update_r(
            r_current = state$r,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = lambda_fixed,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            priors = priors,
            mh_sd = state$r_proposal_sd,
            method = "marginal_nb",
            return_diag = TRUE
        )
        state$r <- r_result$r

        loglik <- compute_s1_loglik_nb(
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            r = state$r
        )

        diagnostics$loglik_trace[iter] <- loglik
        diagnostics$beta_n_reject[iter] <- beta_result$n_reject %||% NA_real_
        diagnostics$phi_accept_trace[iter] <- isTRUE(phi_result$accept)
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
        diagnostics$phi_log_alpha[iter] <- phi_result$log_alpha %||% NA_real_
        diagnostics$r_log_alpha[iter, ] <- r_result$diag$log_alpha %||% rep(NA_real_, n1)

        phi_accept_window <- phi_accept_window + as.integer(isTRUE(phi_result$accept))
        r_accept_window <- r_accept_window + as.integer(r_result$accept)

        if (iter <= n_burnin && iter %% adapt_interval == 0L) {
            state$phi_proposal_sd <- adapt_sd_s1(
                current_sd = state$phi_proposal_sd,
                n_accept = phi_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$phi_target_accept %||% 0.25
            )
            state$r_proposal_sd <- adapt_sd_s1(
                current_sd = state$r_proposal_sd,
                n_accept = r_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$r_target_accept %||% 0.30
            )
            phi_accept_window <- 0L
            r_accept_window <- rep(0L, n1)
        }

        if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
            store_idx <- store_idx + 1L
            samples$beta0[store_idx] <- state$beta0
            samples$beta[store_idx, ] <- state$beta
            samples$phi[store_idx, ] <- state$phi
            samples$tau_phi[store_idx] <- state$tau_phi
            samples$r[store_idx, ] <- state$r
            samples$loglik[store_idx] <- loglik
        }

        if (verbose > 0L && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1L, iter - 99L)
            phi_rate <- mean(diagnostics$phi_accept_trace[i0:iter])
            r_rate <- mean(diagnostics$r_accept_trace[i0:iter, , drop = FALSE])
            beta_rej <- mean(diagnostics$beta_n_reject[i0:iter], na.rm = TRUE)
            cat(sprintf(
                "iter %5d/%d [%.0fs] loglik=%.1f beta0=%.3f beta=(%.3f, %.3f) r_mean=%.2f | phi_acc=%.2f r_acc=%.2f beta_rej=%.1f\n",
                iter, n_iter, elapsed, loglik, state$beta0, state$beta[1],
                state$beta[2], mean(state$r), phi_rate, r_rate, beta_rej
            ))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]
    diagnostics$elapsed_sec <- elapsed_total
    diagnostics$phi_accept_rate <- mean(diagnostics$phi_accept_trace)
    diagnostics$r_accept_rate <- colMeans(diagnostics$r_accept_trace)
    diagnostics$beta_mean_n_reject <- mean(diagnostics$beta_n_reject, na.rm = TRUE)
    diagnostics$phi_proposal_sd_final <- state$phi_proposal_sd
    diagnostics$r_proposal_sd_final <- state$r_proposal_sd

    list(
        samples = samples,
        diagnostics = diagnostics,
        final_state = state,
        settings = settings,
        priors = priors,
        spatial = list(H = spatial$H, B_ICAR = spatial$B_ICAR, BHB = spatial$BHB),
        metadata = list(
            model = "S1 baseline-only NB-ICAR",
            fixed_lambda = TRUE,
            lambda_fixed_value = 1,
            updated_blocks = c("beta", "phi", "tau_phi", "r"),
            disabled_blocks = c("lambda", "gamma", "delta", "omega", "kappa"),
            uses_marginal_nb_for_beta_phi_r = TRUE
        ),
        n_stored = n_stored
    )
}

summarise_s1_baseline_fit <- function(fit, dat, scenario_id = NULL, rep_id = NULL) {
    s <- fit$samples

    beta0_q <- as.numeric(quantile(s$beta0, c(0.025, 0.5, 0.975), na.rm = TRUE))
    beta_q <- apply(s$beta, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    phi_mean <- colMeans(s$phi)
    r_mean_by_region <- colMeans(s$r)
    r_q_by_region <- apply(s$r, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    beta0_true <- dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_
    beta_true <- dat$beta_star %||% rep(NA_real_, 2)
    phi_true <- dat$phi_star_ident %||% dat$phi_star %||% rep(NA_real_, dat$n1)
    r_true <- dat$r_star %||% rep(NA_real_, dat$n1)
    if (length(r_true) == 1L) r_true <- rep(r_true, dat$n1)

    phi_rmse <- sqrt(mean((phi_mean - phi_true)^2))
    phi_cor <- suppressWarnings(cor(phi_mean, phi_true))

    data.frame(
        scenario_id = scenario_id %||% dat$scenario_id %||% "S1_BASELINE_ONLY",
        rep_id = rep_id %||% dat$rep_id %||% NA_integer_,
        TT = dat$TT,
        n1 = dat$n1,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        lambda_fixed_in_fit = isTRUE(fit$metadata$fixed_lambda),
        lambda_fixed_value = fit$metadata$lambda_fixed_value %||% NA_real_,
        lambda_truth_min = if (!is.null(dat$lambda_tilde)) min(dat$lambda_tilde) else NA_real_,
        lambda_truth_max = if (!is.null(dat$lambda_tilde)) max(dat$lambda_tilde) else NA_real_,

        beta0_true = beta0_true,
        beta0_mean = mean(s$beta0),
        beta0_sd = sd(s$beta0),
        beta0_q025 = beta0_q[1],
        beta0_q50 = beta0_q[2],
        beta0_q975 = beta0_q[3],
        beta0_bias = mean(s$beta0) - beta0_true,
        beta0_covered = as.integer(beta0_q[1] <= beta0_true && beta0_true <= beta0_q[3]),

        beta1_true = beta_true[1],
        beta1_mean = mean(s$beta[, 1]),
        beta1_sd = sd(s$beta[, 1]),
        beta1_q025 = beta_q[1, 1],
        beta1_q50 = beta_q[2, 1],
        beta1_q975 = beta_q[3, 1],
        beta1_bias = mean(s$beta[, 1]) - beta_true[1],
        beta1_covered = as.integer(beta_q[1, 1] <= beta_true[1] && beta_true[1] <= beta_q[3, 1]),

        beta2_true = beta_true[2],
        beta2_mean = mean(s$beta[, 2]),
        beta2_sd = sd(s$beta[, 2]),
        beta2_q025 = beta_q[1, 2],
        beta2_q50 = beta_q[2, 2],
        beta2_q975 = beta_q[3, 2],
        beta2_bias = mean(s$beta[, 2]) - beta_true[2],
        beta2_covered = as.integer(beta_q[1, 2] <= beta_true[2] && beta_true[2] <= beta_q[3, 2]),

        phi_rmse = phi_rmse,
        phi_cor = phi_cor,
        phi_mean_bias = mean(phi_mean - phi_true),

        r_true_mean = mean(r_true),
        r_mean = mean(r_mean_by_region),
        r_bias = mean(r_mean_by_region) - mean(r_true),
        r_q025_mean = mean(r_q_by_region[1, ]),
        r_q50_mean = mean(r_q_by_region[2, ]),
        r_q975_mean = mean(r_q_by_region[3, ]),

        tau_phi_mean = mean(s$tau_phi),
        loglik_mean = mean(s$loglik),
        phi_accept_rate = fit$diagnostics$phi_accept_rate,
        r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate),
        beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject,
        elapsed_sec = fit$diagnostics$elapsed_sec,
        stringsAsFactors = FALSE
    )
}

fit_s1_baseline_one_rep <- function(rep_id,
                                    scenario_id = "S1_BASELINE_ONLY",
                                    data_dir = "data_revised",
                                    output_dir = "output_s1_baseline_only",
                                    root = ".",
                                    settings_override = list(),
                                    priors = MCMC_PRIORS,
                                    spatial = build_s1_spatial(),
                                    verbose = 1000L,
                                    save_result = TRUE,
                                    return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    .require_file_s1(dat_file)
    dat <- readRDS(dat_file)
    validate_s1_data(dat)

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    cat(sprintf("=== Scenario 1 baseline-only fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat("Fixed     : lambda_tilde = 1\n")
    cat("Disabled  : gamma, delta, lambda, omega, kappa updates\n")
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r\n\n")

    fit <- run_s1_baseline_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        verbose = verbose
    )

    summary <- summarise_s1_baseline_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id)
    )

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        fit_file <- file.path(out_dir, paste0("fit_S1_baseline_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S1_baseline_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) {
        return(fit)
    }
    invisible(NULL)
}

fit_s1_baseline_batch <- function(reps = 1:10,
                                  scenario_id = "S1_BASELINE_ONLY",
                                  data_dir = "data_revised",
                                  output_dir = "output_s1_baseline_only",
                                  root = ".",
                                  settings_override = list(),
                                  verbose = 1000L,
                                  overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    summaries <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S1_baseline_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s1_baseline_one_rep(
            rep_id = rep_id,
            scenario_id = scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            root = root,
            settings_override = settings_override,
            verbose = verbose,
            save_result = TRUE,
            return_result = TRUE
        )
        summaries[[rr]] <- fit$summary
    }

    summary_all <- do.call(rbind, summaries)
    summary_file <- file.path(out_dir, "summary_S1_baseline_all_reps.csv")
    write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)

    invisible(summary_all)
}

sanity_check_s1_baseline_only <- function(root = ".",
                                           rep_id = 1L,
                                           scenario_id = "S1_BASELINE_ONLY",
                                           data_dir = "data_revised",
                                           output_dir = "output_s1_baseline_only_sanity",
                                           simulate_if_missing = TRUE,
                                           n_iter = 2000L,
                                           n_burnin = 1000L,
                                           n_thin = 2L,
                                           verbose = 500L) {
    dat_dir <- file.path(root, data_dir, scenario_id)
    dat_file <- file.path(dat_dir, paste0("data_rep", sprintf("%02d", as.integer(rep_id)), ".rds"))

    if (!file.exists(dat_file)) {
        if (!isTRUE(simulate_if_missing)) {
            stop("Scenario 1 data file not found: ", dat_file, call. = FALSE)
        }
        simulate_s1_baseline_only_batch(
            reps = rep_id,
            data_dir = data_dir,
            scenario_id = scenario_id,
            root = root,
            overwrite_existing = TRUE,
            verbose = TRUE
        )
    }

    fit_s1_baseline_one_rep(
        rep_id = rep_id,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root,
        settings_override = list(
            n_iter = as.integer(n_iter),
            n_burnin = as.integer(n_burnin),
            n_thin = as.integer(n_thin)
        ),
        verbose = verbose,
        save_result = TRUE,
        return_result = TRUE
    )
}

check_s1_lambda_fixed <- function(fit = NULL, dat = NULL) {
    out <- list()
    if (!is.null(dat)) {
        out$truth_lambda_min <- if (!is.null(dat$lambda_tilde)) min(dat$lambda_tilde) else NA_real_
        out$truth_lambda_max <- if (!is.null(dat$lambda_tilde)) max(dat$lambda_tilde) else NA_real_
        out$truth_all_one <- isTRUE(all.equal(out$truth_lambda_min, 1)) &&
            isTRUE(all.equal(out$truth_lambda_max, 1))
    }
    if (!is.null(fit)) {
        out$fit_fixed_lambda <- isTRUE(fit$metadata$fixed_lambda)
        out$fit_lambda_fixed_value <- fit$metadata$lambda_fixed_value %||% NA_real_
        out$updated_blocks <- fit$metadata$updated_blocks
        out$disabled_blocks <- fit$metadata$disabled_blocks
    }
    out
}

## No automatic command-line execution is included on purpose.
## Source this file, then call the functions explicitly.
