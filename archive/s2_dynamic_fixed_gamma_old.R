## ============================================================================
## s2_dynamic_fixed_gamma.R
##
## Scenario 2 script for the MSSTNB project.
##
## Scenario 2: dynamic residual risk with fixed known gamma.
##
## Purpose:
##   Keep the Scenario 1 baseline validation structure as much as possible, but
##   replace the fixed residual risk lambda_tilde[t, j] = 1 with a genuinely
##   dynamic residual risk path.  During fitting, gamma is fixed at its true
##   known value so that this scenario isolates latent path recovery from
##   discount-factor learning.
##
## Fitted model:
##   Y_tj | kappa_tj, lambda_tilde_tj ~ Poisson(xi_tj * kappa_tj * lambda_tilde_tj)
##   kappa_tj | r_j ~ Gamma(r_j, r_j)
##   log(xi_tj) = log(e_tj) + beta0 + beta1 x1_tj + beta2 x2_tj + phi_j
##   lambda_tilde_{1:T,j} is sampled by Gamma FFBS with gamma_j fixed.
##   phi = B_ICAR u, sum_j phi_j = 0 after deterministic re-centering.
##
## Important functions:
##   source_s2_dynamic_fixed_gamma()
##   simulate_s2_dynamic_fixed_gamma_one()
##   simulate_s2_dynamic_fixed_gamma_batch()
##   run_s2_dynamic_fixed_gamma_mcmc()
##   fit_s2_dynamic_fixed_gamma_one_rep()
##   fit_s2_dynamic_fixed_gamma_batch()
##   sanity_check_s2_dynamic_fixed_gamma()
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.same_dim_s2 <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) && all(as.integer(d) == as.integer(target_dim))
}

.require_file_s2 <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked_s2 <- function(path, verbose = TRUE) {
    .require_file_s2(path)
    if (isTRUE(verbose)) {
        message("source: ", path)
    }
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source project dependencies
## -----------------------------------------------------------------------------
source_s2_dynamic_fixed_gamma <- function(root = ".", verbose = TRUE) {
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
        file.path(mcmc_dir, "update_dispersion_revised.R"),
        file.path(mcmc_dir, "update_kappa.R"),
        file.path(mcmc_dir, "ffbs_lambda_revised.R")
    )

    invisible(lapply(files, .source_checked_s2, verbose = verbose))

    needed <- c(
        "MCMC_SETTINGS", "MCMC_PRIORS",
        "TT", "N1", "N_CHILDREN", "C0", "H", "B_ICAR",
        "REP_SEEDS", "generate_inputs", "generate_icar", "rdirichlet",
        "fit_glm_for_ess", "update_beta", "update_phi",
        "update_tau_phi", "update_r", "update_kappa", "ffbs_lambda_all",
        "recenter"
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
validate_s2_data <- function(dat) {
    required <- c("y_coarse", "e", "x1", "x2", "lambda_tilde", "TT", "n1")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (!is.matrix(dat$y_coarse)) {
        stop("dat$y_coarse must be a matrix.", call. = FALSE)
    }

    y_dim <- dim(dat$y_coarse)
    for (nm in c("e", "x1", "x2", "lambda_tilde")) {
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
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0)) {
        stop("dat$lambda_tilde must be positive and finite.", call. = FALSE)
    }

    invisible(TRUE)
}

build_s2_spatial <- function(H_obj = H, B_obj = B_ICAR) {
    if (!is.matrix(H_obj) || !is.matrix(B_obj)) {
        stop("H_obj and B_obj must be matrices.", call. = FALSE)
    }
    list(
        H = H_obj,
        B_ICAR = B_obj,
        BHB = crossprod(B_obj, H_obj %*% B_obj)
    )
}

validate_s2_fit_objects <- function(dat, priors, spatial) {
    validate_s2_data(dat)

    n1 <- dat$n1
    if (!is.matrix(spatial$H) || .same_dim_s2(spatial$H, c(n1, n1)) == FALSE) {
        stop("spatial$H must be an n1 by n1 matrix.", call. = FALSE)
    }
    if (!is.matrix(spatial$B_ICAR) ||
        .same_dim_s2(spatial$B_ICAR, c(n1, n1 - 1L)) == FALSE) {
        stop("spatial$B_ICAR must be an n1 by n1 - 1 matrix.", call. = FALSE)
    }
    if (!is.matrix(spatial$BHB) ||
        .same_dim_s2(spatial$BHB, c(n1 - 1L, n1 - 1L)) == FALSE) {
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
        stop("This Scenario 2 script expects exactly two covariates.",
             call. = FALSE)
    }

    invisible(TRUE)
}

make_s2_file_name <- function(prefix, rep_id, ext = "rds") {
    paste0(prefix, sprintf("%02d", as.integer(rep_id)), ".", ext)
}

compute_s2_xi <- function(e, x1, x2, beta0, beta, phi) {
    TT_use <- nrow(e)
    n1_use <- ncol(e)
    xi <- matrix(NA_real_, TT_use, n1_use)
    for (j in seq_len(n1_use)) {
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        eta_j <- pmin(pmax(eta_j, -700), 700)
        xi[, j] <- e[, j] * exp(eta_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("compute_s2_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

compute_s2_loglik_poisson_aug <- function(y_coarse, xi, lambda_tilde, kappa) {
    mu <- pmax(xi * lambda_tilde * kappa, .Machine$double.xmin)
    sum(stats::dpois(y_coarse, lambda = mu, log = TRUE))
}

compute_s2_loglik_nb <- function(y_coarse, e, x1, x2, beta0, beta, phi, r, lambda_tilde) {
    n1 <- ncol(y_coarse)
    out <- 0
    for (j in seq_len(n1)) {
        log_mu_j <- log(e[, j]) + beta0 + beta[1] * x1[, j] +
            beta[2] * x2[, j] + phi[j] + log(pmax(lambda_tilde[, j], .Machine$double.xmin))
        out <- out + .nb_or_poisson_loglik(y = y_coarse[, j],
                                           log_mu = log_mu_j,
                                           r = r[j])
    }
    out
}

compute_u_from_phi_s2 <- function(phi, B) {
    BtB <- crossprod(B)
    rhs <- crossprod(B, phi)
    as.numeric(solve(BtB, rhs))
}

## Internal input generator. This intentionally mirrors Scenario 1 so that
## Scenario 2 changes only the dynamic residual risk and the fitting blocks
## required to recover it.
generate_s2_dynamic_inputs <- function(TT_use,
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
            stop("Non-finite values generated in Scenario 2 input matrix: ", nm, call. = FALSE)
        }
    }

    out
}

## -----------------------------------------------------------------------------
## Scenario 2 simulation
## -----------------------------------------------------------------------------
simulate_s2_dynamic_fixed_gamma_one <- function(seed = 1L,
                                                TT_use = TT,
                                                n1_use = N1,
                                                n_children_use = N_CHILDREN,
                                                c0 = C0,
                                                A0_use = if (exists("A0", envir = .GlobalEnv)) A0 else 10,
                                                B0_use = if (exists("B0", envir = .GlobalEnv)) B0 else 10,
                                                beta0_truth = -1.5,
                                                beta_truth = c(0.5, -0.4),
                                                r_truth = 15,
                                                tau_phi_truth = 2,
                                                gamma_truth = 0.8,
                                                delta_truth = 0.9,
                                                lambda_initial_value = 1,
                                                omega_mode = c("fixed_prior_mean", "dirichlet_static"),
                                                H_obj = H,
                                                B_obj = B_ICAR,
                                                scenario_id = "S2_DYNAMIC_FIXED_GAMMA",
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
                                                lambda_floor = .Machine$double.xmin,
                                                max_poisson_rate = 1e7) {
    omega_mode <- match.arg(omega_mode)
    x2_mode <- match.arg(x2_mode)

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
    if (length(gamma_truth) != n1_use || any(!is.finite(gamma_truth)) ||
        any(gamma_truth <= 0) || any(gamma_truth >= 1)) {
        stop("gamma_truth must be scalar or length n1_use, with entries in (0, 1).",
             call. = FALSE)
    }

    if (isTRUE(use_internal_inputs)) {
        inputs <- generate_s2_dynamic_inputs(
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

    target_dim <- c(as.integer(TT_use), as.integer(n1_use))
    if (!.same_dim_s2(e, target_dim) ||
        !.same_dim_s2(x1, target_dim) ||
        !.same_dim_s2(x2, target_dim)) {
        stop(sprintf(
            paste0(
                "Scenario 2 inputs do not have dimensions TT_use x n1_use. ",
                "Expected %d x %d; got e=%s, x1=%s, x2=%s."
            ),
            as.integer(TT_use), as.integer(n1_use),
            paste(dim(e), collapse = " x "),
            paste(dim(x1), collapse = " x "),
            paste(dim(x2), collapse = " x ")
        ), call. = FALSE)
    }

    if (!.same_dim_s2(H_obj, c(n1_use, n1_use)) ||
        !.same_dim_s2(B_obj, c(n1_use, n1_use - 1L))) {
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
        stop("Scenario 2 inputs contain non-finite values before count generation.", call. = FALSE)
    }

    phi_truth <- generate_icar(tau_phi = tau_phi_truth, H = H_obj, B = B_obj)
    phi_truth <- phi_truth - mean(phi_truth)
    if (length(phi_truth) != n1_use || any(!is.finite(phi_truth))) {
        stop("Generated phi_truth is non-finite or has the wrong length.", call. = FALSE)
    }

    xi <- compute_s2_xi(e = e, x1 = x1, x2 = x2,
                        beta0 = beta0_truth, beta = beta_truth, phi = phi_truth)

    lambda_tilde <- matrix(NA_real_, nrow = TT_use, ncol = n1_use)
    lambda_tilde0 <- rep(lambda_initial_value, n1_use)
    mu_nb <- matrix(NA_real_, TT_use, n1_use)
    kappa <- matrix(NA_real_, TT_use, n1_use)
    poisson_rate <- matrix(NA_real_, TT_use, n1_use)
    y_coarse <- matrix(NA_integer_, TT_use, n1_use)

    ## Sequential beta-gamma discount evolution.  The filtering shapes are
    ## updated after each generated observation.  This aligns the DGP with the
    ## Gamma FFBS block used during fitting while keeping gamma fixed and known.
    a_prev <- rep(A0_use, n1_use)
    b_prev <- rep(B0_use, n1_use)
    lambda_prev <- lambda_tilde0

    for (t in seq_len(TT_use)) {
        for (j in seq_len(n1_use)) {
            shape1 <- max(gamma_truth[j] * a_prev[j], 1e-12)
            shape2 <- max((1 - gamma_truth[j]) * a_prev[j], 1e-12)
            eta_tj <- stats::rbeta(1L, shape1 = shape1, shape2 = shape2)
            lambda_tilde[t, j] <- max(lambda_prev[j] * eta_tj / gamma_truth[j], lambda_floor)

            kappa[t, j] <- stats::rgamma(1L, shape = r_truth[j], rate = r_truth[j])
            mu_nb[t, j] <- xi[t, j] * lambda_tilde[t, j]
            poisson_rate[t, j] <- mu_nb[t, j] * kappa[t, j]

            bad_rate <- !is.finite(poisson_rate[t, j]) || poisson_rate[t, j] < 0 ||
                poisson_rate[t, j] > max_poisson_rate
            if (bad_rate) {
                stop(sprintf(
                    paste0(
                        "Non-finite, negative, or excessively large Poisson rate in Scenario 2 DGP. ",
                        "First bad cell: t=%d, j=%d, rate=%s, lambda=%s, xi=%s. ",
                        "Try smaller exposure_max, less variable covariates, or a smaller beta0_truth."
                    ),
                    t, j,
                    as.character(poisson_rate[t, j]),
                    as.character(lambda_tilde[t, j]),
                    as.character(xi[t, j])
                ), call. = FALSE)
            }

            y_coarse[t, j] <- stats::rpois(1L, lambda = poisson_rate[t, j])

            a_prev[j] <- gamma_truth[j] * a_prev[j] + y_coarse[t, j]
            b_prev[j] <- gamma_truth[j] * b_prev[j] + xi[t, j] * kappa[t, j]
            lambda_prev[j] <- lambda_tilde[t, j]
        }
    }

    if (anyNA(y_coarse)) {
        stop("rpois produced NA values in Scenario 2 DGP after rate validation.", call. = FALSE)
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
        stop("Scenario 2 generated fine counts are not coherent with coarse counts.", call. = FALSE)
    }

    ## Identified truth used for posterior comparison.  The model is recentered
    ## during MCMC, so comparing to raw beta0/phi/lambda is not meaningful.
    rc_truth <- recenter(
        beta0 = beta0_truth,
        phi = phi_truth,
        lambda_tilde = lambda_tilde,
        return_diag = TRUE
    )

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
        beta0_star_ident = rc_truth$beta0,
        beta_star_ident = beta_truth,
        phi_star_ident = rc_truth$phi,
        lambda_tilde_ident = rc_truth$lambda_tilde,
        lambda_recenter_diag = rc_truth$diag,
        lambda_dynamic = TRUE,
        lambda_initial_value = lambda_initial_value,
        A0 = A0_use,
        B0 = B0_use,
        omega_mode = omega_mode,
        data_type = "dynamic_lambda_fixed_gamma",
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
        input_generator = if (isTRUE(use_internal_inputs)) "generate_s2_dynamic_inputs" else "generate_inputs"
    )

    if (!is.null(inputs$x1_raw)) out$x1_raw <- inputs$x1_raw
    if (!is.null(inputs$x2_raw)) out$x2_raw <- inputs$x2_raw
    if (!is.null(inputs$x2_mode)) out$x2_mode <- inputs$x2_mode

    validate_s2_data(out)
    out
}

simulate_s2_dynamic_fixed_gamma_batch <- function(reps = 1:10,
                                                  data_dir = "data_revised",
                                                  scenario_id = "S2_DYNAMIC_FIXED_GAMMA",
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

            dat <- simulate_s2_dynamic_fixed_gamma_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                message(sprintf(
                    "Saved %s | mean_count=%.2f zero_prop=%.3f lambda_range=[%.3f, %.3f]",
                    out_file, dat$mean_count, dat$zero_prop,
                    min(dat$lambda_tilde_ident), max(dat$lambda_tilde_ident)
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
            gamma_truth_mean = mean(dat$gamma_star),
            lambda_raw_min = min(dat$lambda_tilde),
            lambda_raw_max = max(dat$lambda_tilde),
            lambda_ident_min = min(dat$lambda_tilde_ident),
            lambda_ident_max = max(dat$lambda_tilde_ident),
            lambda_ident_log_rm_mean = mean(abs(colMeans(log(dat$lambda_tilde_ident)))),
            coherent = isTRUE(dat$coherent),
            stringsAsFactors = FALSE
        )
    }

    manifest_df <- do.call(rbind, manifest)
    manifest_file <- file.path(out_dir, "manifest_S2_dynamic_fixed_gamma.csv")
    write.csv(manifest_df, manifest_file, row.names = FALSE)
    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
    }

    invisible(manifest_df)
}

## -----------------------------------------------------------------------------
## Scenario 2 fitting
## -----------------------------------------------------------------------------
adapt_sd_s2 <- function(current_sd, n_accept, n_trials,
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

initialise_s2_state <- function(dat, settings, priors, spatial, gamma_fixed) {
    n1 <- dat$n1
    TT_use <- dat$TT
    p <- length(priors$beta_mean)

    ess_base <- fit_glm_for_ess(dat)
    beta0 <- ess_base$center[1]
    beta <- ess_base$center[2:(p + 1L)]

    phi <- rep(0, n1)
    lambda_init <- matrix(1, TT_use, n1)
    xi <- compute_s2_xi(dat$e, dat$x1, dat$x2, beta0, beta, phi)

    list(
        beta0 = beta0,
        beta = beta,
        phi = phi,
        u = rep(0, n1 - 1L),
        tau_phi = settings$tau_phi_init %||% 1,
        r = rep(settings$r_init %||% 10, n1),
        kappa = matrix(1, TT_use, n1),
        lambda_tilde = lambda_init,
        gamma = gamma_fixed,
        xi = xi,
        ess_base = ess_base,
        phi_proposal_sd = rep(settings$phi_proposal_sd_init %||% 0.01, n1 - 1L),
        r_proposal_sd = rep(settings$r_proposal_sd_init %||% 0.50, n1)
    )
}

run_s2_dynamic_fixed_gamma_mcmc <- function(dat,
                                            settings = MCMC_SETTINGS,
                                            priors = MCMC_PRIORS,
                                            spatial = build_s2_spatial(),
                                            A0_use = dat$A0 %||% if (exists("A0", envir = .GlobalEnv)) A0 else 10,
                                            B0_use = dat$B0 %||% if (exists("B0", envir = .GlobalEnv)) B0 else 10,
                                            gamma_fixed = dat$gamma_star %||% 0.8,
                                            verbose = 1000L) {
    validate_s2_fit_objects(dat, priors, spatial)

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

    if (length(gamma_fixed) == 1L) {
        gamma_fixed <- rep(gamma_fixed, n1)
    }
    if (length(gamma_fixed) != n1 || any(!is.finite(gamma_fixed)) ||
        any(gamma_fixed <= 0) || any(gamma_fixed >= 1)) {
        stop("gamma_fixed must be scalar or length n1 with entries in (0, 1).",
             call. = FALSE)
    }

    samples <- list(
        beta0 = numeric(n_stored),
        beta = matrix(NA_real_, nrow = n_stored, ncol = p),
        phi = matrix(NA_real_, nrow = n_stored, ncol = n1),
        tau_phi = numeric(n_stored),
        r = matrix(NA_real_, nrow = n_stored, ncol = n1),
        gamma = matrix(rep(gamma_fixed, each = n_stored), nrow = n_stored, ncol = n1),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT_use, n1)),
        kappa = array(NA_real_, dim = c(n_stored, TT_use, n1)),
        loglik = numeric(n_stored),
        loglik_nb = numeric(n_stored)
    )

    diagnostics <- list(
        loglik_trace = rep(NA_real_, n_iter),
        loglik_nb_trace = rep(NA_real_, n_iter),
        beta_n_reject = rep(NA_real_, n_iter),
        phi_accept_trace = rep(FALSE, n_iter),
        r_accept_trace = matrix(FALSE, nrow = n_iter, ncol = n1),
        phi_log_alpha = rep(NA_real_, n_iter),
        r_log_alpha = matrix(NA_real_, nrow = n_iter, ncol = n1),
        lambda_min_trace = rep(NA_real_, n_iter),
        lambda_max_trace = rep(NA_real_, n_iter),
        kappa_mean_trace = rep(NA_real_, n_iter),
        recenter_error = rep(NA_real_, n_iter),
        fixed_gamma = gamma_fixed
    )

    state <- initialise_s2_state(dat, settings, priors, spatial, gamma_fixed)
    phi_accept_window <- 0L
    r_accept_window <- rep(0L, n1)
    store_idx <- 0L
    start_time <- proc.time()

    for (iter in seq_len(n_iter)) {
        ## 1. beta update using marginal NB likelihood, conditional on current lambda.
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = state$lambda_tilde,
            phi = state$phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r,
            use_preconditioned = settings$use_preconditioned_beta %||% TRUE
        )
        state$beta0 <- beta_result$sample[1]
        state$beta <- beta_result$sample[2:3]
        state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                                  state$beta0, state$beta, state$phi)

        ## 2. phi update using marginal NB likelihood, conditional on current lambda.
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
            lambda_tilde = state$lambda_tilde,
            r = state$r,
            proposal_sd = state$phi_proposal_sd
        )
        state$u <- phi_result$u
        state$phi <- phi_result$phi
        state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                                  state$beta0, state$beta, state$phi)

        ## 3. tau_phi update.
        state$tau_phi <- update_tau_phi(
            phi = state$phi,
            H = spatial$H,
            n1 = n1,
            priors = priors
        )

        ## 4. r update using marginal NB likelihood, conditional on current lambda.
        r_result <- update_r(
            r_current = state$r,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = state$lambda_tilde,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            priors = priors,
            mh_sd = state$r_proposal_sd,
            method = "marginal_nb",
            return_diag = TRUE
        )
        state$r <- r_result$r

        ## 5. kappa update, used only for the conditional lambda FFBS block.
        kappa_result <- update_kappa(
            y_coarse = dat$y_coarse,
            lambda_tilde = state$lambda_tilde,
            xi = state$xi,
            r = state$r,
            return_diag = TRUE
        )
        state$kappa <- kappa_result$kappa

        ## 6. lambda FFBS with gamma fixed at the known value.
        lambda_result <- ffbs_lambda_all(
            gamma = gamma_fixed,
            y_coarse = dat$y_coarse,
            xi = state$xi,
            kappa = state$kappa,
            a0 = A0_use,
            b0 = B0_use,
            return_diag = TRUE
        )
        state$lambda_tilde <- lambda_result$lambda_tilde

        ## 7. Deterministic re-centering for identified beta0, phi, and lambda.
        rc <- recenter(
            beta0 = state$beta0,
            phi = state$phi,
            lambda_tilde = state$lambda_tilde,
            return_diag = TRUE
        )
        state$beta0 <- rc$beta0
        state$phi <- rc$phi
        state$lambda_tilde <- rc$lambda_tilde
        state$u <- compute_u_from_phi_s2(state$phi, spatial$B_ICAR)
        state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                                  state$beta0, state$beta, state$phi)

        loglik <- compute_s2_loglik_poisson_aug(
            y_coarse = dat$y_coarse,
            xi = state$xi,
            lambda_tilde = state$lambda_tilde,
            kappa = state$kappa
        )
        loglik_nb <- compute_s2_loglik_nb(
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            r = state$r,
            lambda_tilde = state$lambda_tilde
        )

        diagnostics$loglik_trace[iter] <- loglik
        diagnostics$loglik_nb_trace[iter] <- loglik_nb
        diagnostics$beta_n_reject[iter] <- beta_result$n_reject %||% NA_real_
        diagnostics$phi_accept_trace[iter] <- isTRUE(phi_result$accept)
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
        diagnostics$phi_log_alpha[iter] <- phi_result$log_alpha %||% NA_real_
        diagnostics$r_log_alpha[iter, ] <- r_result$diag$log_alpha %||% rep(NA_real_, n1)
        diagnostics$lambda_min_trace[iter] <- lambda_result$diag$min_lambda %||% min(state$lambda_tilde)
        diagnostics$lambda_max_trace[iter] <- lambda_result$diag$max_lambda %||% max(state$lambda_tilde)
        diagnostics$kappa_mean_trace[iter] <- kappa_result$diag$mean_kappa %||% mean(state$kappa)
        diagnostics$recenter_error[iter] <- rc$diag$max_abs_log_core_difference %||% NA_real_

        phi_accept_window <- phi_accept_window + as.integer(isTRUE(phi_result$accept))
        r_accept_window <- r_accept_window + as.integer(r_result$accept)

        if (iter <= n_burnin && iter %% adapt_interval == 0L) {
            state$phi_proposal_sd <- adapt_sd_s2(
                current_sd = state$phi_proposal_sd,
                n_accept = phi_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$phi_target_accept %||% 0.25
            )
            state$r_proposal_sd <- adapt_sd_s2(
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
            samples$lambda_tilde[store_idx, , ] <- state$lambda_tilde
            samples$kappa[store_idx, , ] <- state$kappa
            samples$loglik[store_idx] <- loglik
            samples$loglik_nb[store_idx] <- loglik_nb
        }

        if (verbose > 0L && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1L, iter - 99L)
            phi_rate <- mean(diagnostics$phi_accept_trace[i0:iter])
            r_rate <- mean(diagnostics$r_accept_trace[i0:iter, , drop = FALSE])
            beta_rej <- mean(diagnostics$beta_n_reject[i0:iter], na.rm = TRUE)
            cat(sprintf(
                paste0(
                    "iter %5d/%d [%.0fs] loglik_nb=%.1f beta0=%.3f ",
                    "beta=(%.3f, %.3f) r_mean=%.2f lambda=[%.3f, %.3f] ",
                    "gamma_fixed=%.3f | phi_acc=%.2f r_acc=%.2f beta_rej=%.1f\n"
                ),
                iter, n_iter, elapsed, loglik_nb, state$beta0, state$beta[1],
                state$beta[2], mean(state$r), min(state$lambda_tilde),
                max(state$lambda_tilde), mean(gamma_fixed), phi_rate, r_rate, beta_rej
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
    diagnostics$gamma_fixed_mean <- mean(gamma_fixed)
    diagnostics$gamma_fixed_min <- min(gamma_fixed)
    diagnostics$gamma_fixed_max <- max(gamma_fixed)

    list(
        samples = samples,
        diagnostics = diagnostics,
        final_state = state,
        settings = settings,
        priors = priors,
        spatial = list(H = spatial$H, B_ICAR = spatial$B_ICAR, BHB = spatial$BHB),
        metadata = list(
            model = "S2 dynamic residual risk NB-ICAR with fixed gamma",
            fixed_lambda = FALSE,
            dynamic_lambda = TRUE,
            fixed_gamma = TRUE,
            gamma_fixed_value = gamma_fixed,
            A0 = A0_use,
            B0 = B0_use,
            updated_blocks = c("beta", "phi", "tau_phi", "r", "kappa", "lambda"),
            disabled_blocks = c("gamma", "delta", "omega"),
            uses_marginal_nb_for_beta_phi_r = TRUE,
            uses_kappa_for_lambda_ffbs = TRUE,
            uses_recentered_identified_parameterization = TRUE
        ),
        n_stored = n_stored
    )
}

summarise_s2_lambda_recovery <- function(lambda_draws, lambda_truth) {
    if (length(dim(lambda_draws)) != 3L) {
        stop("lambda_draws must be an array with dimensions draws x TT x n1.", call. = FALSE)
    }
    if (!identical(dim(lambda_draws)[2:3], dim(lambda_truth))) {
        stop("lambda_draws and lambda_truth dimensions are incompatible.", call. = FALSE)
    }

    lambda_floor <- .Machine$double.xmin
    lambda_draws <- pmax(lambda_draws, lambda_floor)
    lambda_truth <- pmax(lambda_truth, lambda_floor)

    lambda_mean <- apply(lambda_draws, c(2, 3), mean, na.rm = TRUE)
    lambda_q025 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
    lambda_q975 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)

    log_draws <- log(lambda_draws)
    log_truth <- log(lambda_truth)
    log_mean <- apply(log_draws, c(2, 3), mean, na.rm = TRUE)
    log_q025 <- apply(log_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
    log_q975 <- apply(log_draws, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)

    dlog_mean <- apply(log_mean, 2, diff)
    dlog_truth <- apply(log_truth, 2, diff)

    data.frame(
        lambda_rmse = sqrt(mean((lambda_mean - lambda_truth)^2)),
        lambda_mae = mean(abs(lambda_mean - lambda_truth)),
        lambda_coverage_95 = mean(lambda_q025 <= lambda_truth & lambda_truth <= lambda_q975),
        log_lambda_rmse = sqrt(mean((log_mean - log_truth)^2)),
        log_lambda_mae = mean(abs(log_mean - log_truth)),
        log_lambda_coverage_95 = mean(log_q025 <= log_truth & log_truth <= log_q975),
        cor_lambda = suppressWarnings(stats::cor(as.numeric(lambda_mean), as.numeric(lambda_truth))),
        cor_log_lambda = suppressWarnings(stats::cor(as.numeric(log_mean), as.numeric(log_truth))),
        rmse_delta_log_lambda = sqrt(mean((dlog_mean - dlog_truth)^2)),
        mae_delta_log_lambda = mean(abs(dlog_mean - dlog_truth)),
        cor_delta_log_lambda = suppressWarnings(stats::cor(as.numeric(dlog_mean), as.numeric(dlog_truth))),
        lambda_truth_min = min(lambda_truth),
        lambda_truth_max = max(lambda_truth),
        lambda_post_mean_min = min(lambda_mean),
        lambda_post_mean_max = max(lambda_mean),
        stringsAsFactors = FALSE
    )
}

summarise_s2_dynamic_fixed_gamma_fit <- function(fit, dat, scenario_id = NULL, rep_id = NULL) {
    s <- fit$samples

    beta0_q <- as.numeric(stats::quantile(s$beta0, c(0.025, 0.5, 0.975), na.rm = TRUE))
    beta_q <- apply(s$beta, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    phi_mean <- colMeans(s$phi)
    r_mean_by_region <- colMeans(s$r)
    r_q_by_region <- apply(s$r, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    beta0_true <- dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_
    beta_true <- dat$beta_star_ident %||% dat$beta_star %||% rep(NA_real_, 2)
    phi_true <- dat$phi_star_ident %||% dat$phi_star %||% rep(NA_real_, dat$n1)
    r_true <- dat$r_star %||% rep(NA_real_, dat$n1)
    if (length(r_true) == 1L) r_true <- rep(r_true, dat$n1)

    lambda_truth <- dat$lambda_tilde_ident %||% dat$lambda_tilde
    lambda_rec <- summarise_s2_lambda_recovery(s$lambda_tilde, lambda_truth)

    phi_rmse <- sqrt(mean((phi_mean - phi_true)^2))
    phi_cor <- suppressWarnings(stats::cor(phi_mean, phi_true))

    cbind(
        data.frame(
            scenario_id = scenario_id %||% dat$scenario_id %||% "S2_DYNAMIC_FIXED_GAMMA",
            rep_id = rep_id %||% dat$rep_id %||% NA_integer_,
            TT = dat$TT,
            n1 = dat$n1,
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            dynamic_lambda_in_truth = TRUE,
            lambda_sampled_in_fit = TRUE,
            gamma_fixed_in_fit = isTRUE(fit$metadata$fixed_gamma),
            gamma_truth_mean = mean(dat$gamma_star),
            gamma_fixed_mean = mean(fit$metadata$gamma_fixed_value),
            lambda_truth_min = if (!is.null(lambda_truth)) min(lambda_truth) else NA_real_,
            lambda_truth_max = if (!is.null(lambda_truth)) max(lambda_truth) else NA_real_,

            beta0_true = beta0_true,
            beta0_mean = mean(s$beta0),
            beta0_sd = stats::sd(s$beta0),
            beta0_q025 = beta0_q[1],
            beta0_q50 = beta0_q[2],
            beta0_q975 = beta0_q[3],
            beta0_bias = mean(s$beta0) - beta0_true,
            beta0_covered = as.integer(beta0_q[1] <= beta0_true && beta0_true <= beta0_q[3]),

            beta1_true = beta_true[1],
            beta1_mean = mean(s$beta[, 1]),
            beta1_sd = stats::sd(s$beta[, 1]),
            beta1_q025 = beta_q[1, 1],
            beta1_q50 = beta_q[2, 1],
            beta1_q975 = beta_q[3, 1],
            beta1_bias = mean(s$beta[, 1]) - beta_true[1],
            beta1_covered = as.integer(beta_q[1, 1] <= beta_true[1] && beta_true[1] <= beta_q[3, 1]),

            beta2_true = beta_true[2],
            beta2_mean = mean(s$beta[, 2]),
            beta2_sd = stats::sd(s$beta[, 2]),
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
            loglik_nb_mean = mean(s$loglik_nb),
            phi_accept_rate = fit$diagnostics$phi_accept_rate,
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate),
            beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject,
            elapsed_sec = fit$diagnostics$elapsed_sec,
            stringsAsFactors = FALSE
        ),
        lambda_rec
    )
}

fit_s2_dynamic_fixed_gamma_one_rep <- function(rep_id,
                                               scenario_id = "S2_DYNAMIC_FIXED_GAMMA",
                                               data_dir = "data_revised",
                                               output_dir = "output_s2_dynamic_fixed_gamma",
                                               root = ".",
                                               settings_override = list(),
                                               priors = MCMC_PRIORS,
                                               spatial = build_s2_spatial(),
                                               gamma_fixed = NULL,
                                               verbose = 1000L,
                                               save_result = TRUE,
                                               return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    .require_file_s2(dat_file)
    dat <- readRDS(dat_file)
    validate_s2_data(dat)

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    gamma_use <- gamma_fixed %||% dat$gamma_star

    cat(sprintf("=== Scenario 2 dynamic fixed-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Fixed     : gamma = %.3f on average\n", mean(gamma_use)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    fit <- run_s2_dynamic_fixed_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_fixed = gamma_use,
        verbose = verbose
    )

    summary <- summarise_s2_dynamic_fixed_gamma_fit(
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
        fit_file <- file.path(out_dir, paste0("fit_S2_dynamic_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S2_dynamic_fixed_gamma_rep", rr, ".csv"))
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

fit_s2_dynamic_fixed_gamma_batch <- function(reps = 1:10,
                                             scenario_id = "S2_DYNAMIC_FIXED_GAMMA",
                                             data_dir = "data_revised",
                                             output_dir = "output_s2_dynamic_fixed_gamma",
                                             root = ".",
                                             settings_override = list(),
                                             verbose = 1000L,
                                             overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    summaries <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S2_dynamic_fixed_gamma_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s2_dynamic_fixed_gamma_one_rep(
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
    summary_file <- file.path(out_dir, "summary_S2_dynamic_fixed_gamma_all_reps.csv")
    write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)

    invisible(summary_all)
}

sanity_check_s2_dynamic_fixed_gamma <- function(root = ".",
                                                rep_id = 1L,
                                                scenario_id = "S2_DYNAMIC_FIXED_GAMMA",
                                                data_dir = "data_revised",
                                                output_dir = "output_s2_dynamic_fixed_gamma_sanity",
                                                simulate_if_missing = TRUE,
                                                n_iter = 2000L,
                                                n_burnin = 1000L,
                                                n_thin = 2L,
                                                verbose = 500L) {
    dat_dir <- file.path(root, data_dir, scenario_id)
    dat_file <- file.path(dat_dir, paste0("data_rep", sprintf("%02d", as.integer(rep_id)), ".rds"))

    if (!file.exists(dat_file)) {
        if (!isTRUE(simulate_if_missing)) {
            stop("Scenario 2 data file not found: ", dat_file, call. = FALSE)
        }
        simulate_s2_dynamic_fixed_gamma_batch(
            reps = rep_id,
            data_dir = data_dir,
            scenario_id = scenario_id,
            root = root,
            overwrite_existing = TRUE,
            verbose = TRUE
        )
    }

    fit_s2_dynamic_fixed_gamma_one_rep(
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

check_s2_dynamic_fixed_gamma <- function(fit = NULL, dat = NULL) {
    out <- list()
    if (!is.null(dat)) {
        lambda_truth <- dat$lambda_tilde_ident %||% dat$lambda_tilde
        out$truth_lambda_min <- min(lambda_truth)
        out$truth_lambda_max <- max(lambda_truth)
        out$truth_lambda_dynamic <- isTRUE(sd(as.numeric(lambda_truth)) > 0)
        out$truth_gamma_mean <- mean(dat$gamma_star)
    }
    if (!is.null(fit)) {
        out$fit_fixed_gamma <- isTRUE(fit$metadata$fixed_gamma)
        out$fit_gamma_fixed_mean <- mean(fit$metadata$gamma_fixed_value)
        out$updated_blocks <- fit$metadata$updated_blocks
        out$disabled_blocks <- fit$metadata$disabled_blocks
        out$lambda_sampled <- "lambda" %in% fit$metadata$updated_blocks
    }
    out
}

## No automatic command-line execution is included on purpose.
## Source this file, then call the functions explicitly.
