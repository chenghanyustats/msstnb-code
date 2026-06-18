## ============================================================================
## run_s4a_oracle_lambda_fixed_gamma_diagnostic_T100.R
##
## Scenario 4A diagnostic: Oracle lambda + fixed gamma.
##
## Purpose
##   Diagnose whether the S4A sparse-count instability in reps 04 and 07 is
##   mainly caused by a beta--lambda ridge.  This script fixes lambda_tilde at
##   the true identified lambda path from the S4A data and fixes gamma at truth
##   (0.8), while estimating beta0, beta1, beta2, phi, tau_phi, r, and kappa.
##
## Interpretation
##   - If reps 04 and 07 become stable with oracle lambda, the original failures
##     are strong evidence of beta--lambda confounding / weak identification.
##   - If they remain unstable, the next diagnostic should fix both lambda and
##     phi and estimate only beta/r, to isolate beta-likelihood/covariate issues.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   data_s4a_sparse_counts/S4A_SPARSE_COUNTS_OBS_T100/data_rep04.rds
##   data_s4a_sparse_counts/S4A_SPARSE_COUNTS_OBS_T100/data_rep07.rds
##
## Main outputs
##   output_s4a_oracle_lambda_fixed_gamma/
##     S4A_ORACLE_LAMBDA_FIXED_GAMMA_DIAGNOSTIC_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4A_ORACLE_LAMBDA_FIXED_GAMMA_DIAGNOSTIC_T100"
data_scenario_id <- "S4A_SPARSE_COUNTS_OBS_T100"

## This diagnostic targets the previously unstable S4A reps.
reps_formal <- c(4L, 7L)

## Use the same formal MCMC profile as S4A fixed-gamma runs.
n_iter <- 40000L
n_burnin <- 20000L
n_thin <- 5L

fixed_gamma_value <- 0.8
gamma_prior <- c(1, 1)
oracle_lambda_field <- "lambda_tilde_ident"

s3_core_file <- "s3_dynamic_learned_gamma.R"
data_dir <- "data_s4a_sparse_counts"
output_dir <- "output_s4a_oracle_lambda_fixed_gamma"
verbose <- 1000L
overwrite_fit <- TRUE

## Expected official S4A data settings.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "observation_sparse_counts"
expected_sparse_beta0_shift <- -4.25
expected_beta0_reference_truth <- -1.5
expected_beta0_sparse_truth <- -5.75

## Numerical guard settings: intentionally wide relative to the truth.
s4a_log_xi_lower <- -40
s4a_log_xi_upper <-  40
s4a_beta0_bounds <- c(-30, 10)
s4a_beta_bounds <- c(-5, 5)
s4a_kappa_bounds <- c(1e-10, 1e10)

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

s4a_data_file <- function(rep_id,
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

check_s4a_oracle_source_dataset <- function(data_file,
                                            TT_expected = TT_use,
                                            n1_expected = n1_use,
                                            beta0_ident_abs_limit = 20) {
    assert_file_exists(data_file, "S4A source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "y_coarse", "e", "x1", "x2",
        "lambda_tilde", "lambda_tilde_ident", "gamma_star", "beta0_star_ident",
        "beta_star_ident", "phi_star_ident", "sparse_beta0_shift",
        "beta0_reference_truth", "beta0_sparse_truth", "expected_count_multiplier",
        "mean_count", "zero_prop", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true(
        length(missing) == 0L,
        paste("Dataset is missing required fields:", paste(missing, collapse = ", "))
    )

    assert_true(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_expected), as.integer(n1_expected))),
        paste0("y_coarse dimension is not ", TT_expected, " x ", n1_expected, ".")
    )

    for (nm in c("e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        assert_true(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    assert_true(identical(dat$stress_type, expected_stress_type),
                paste0("stress_type is not ", expected_stress_type, "."))
    assert_true(abs(dat$sparse_beta0_shift - expected_sparse_beta0_shift) < 1e-12,
                paste0("sparse_beta0_shift is not ", expected_sparse_beta0_shift, "."))
    assert_true(abs(dat$beta0_reference_truth - expected_beta0_reference_truth) < 1e-12,
                paste0("beta0_reference_truth is not ", expected_beta0_reference_truth, "."))
    assert_true(abs(dat$beta0_sparse_truth - expected_beta0_sparse_truth) < 1e-12,
                paste0("beta0_sparse_truth is not ", expected_beta0_sparse_truth, "."))

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
    if (any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("dat$lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    assert_true(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident)
    )

    list(
        dat = dat,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        beta0_star_ident = dat$beta0_star_ident,
        lambda_oracle_min = min(dat[[oracle_lambda_field]]),
        lambda_oracle_max = max(dat[[oracle_lambda_field]])
    )
}

## ---- source Scenario 3 model code ------------------------------------------
assert_file_exists(file.path(root_dir, s3_core_file), "Scenario 3 core script")
source(file.path(root_dir, s3_core_file))
source_s3_dynamic_learned_gamma(root = root_dir)

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
        log_xi_j <- pmin(pmax(log_xi_j, s4a_log_xi_lower), s4a_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("Safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

.s4a_oracle_guard_env <- new.env(parent = emptyenv())
reset_oracle_guards <- function() {
    .s4a_oracle_guard_env$n_beta_guard <- 0L
    .s4a_oracle_guard_env$n_kappa_guard <- 0L
    invisible(TRUE)
}
get_oracle_guards <- function() {
    list(
        n_beta_guard = .s4a_oracle_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4a_oracle_guard_env$n_kappa_guard %||% 0L
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
        bad <- smp[1] < s4a_beta0_bounds[1] || smp[1] > s4a_beta0_bounds[2] ||
            any(smp[-1] < s4a_beta_bounds[1] | smp[-1] > s4a_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s4a_oracle_guard_env$n_beta_guard <- (.s4a_oracle_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s4a_oracle_guard_rejected <- TRUE
    } else {
        res$s4a_oracle_guard_rejected <- FALSE
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
        .s4a_oracle_guard_env$n_kappa_guard <- (.s4a_oracle_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s4a_oracle_guard_env$n_kappa_guard <- (.s4a_oracle_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s4a_kappa_bounds[1]), s4a_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s4a_oracle_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

cat(sprintf(
    "Using oracle-lambda diagnostic guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s4a_beta0_bounds[1], s4a_beta0_bounds[2],
    s4a_beta_bounds[1], s4a_beta_bounds[2],
    s4a_log_xi_lower, s4a_log_xi_upper
))

## ---- oracle-lambda MCMC -----------------------------------------------------
run_s4a_oracle_lambda_fixed_gamma_mcmc <- function(dat,
                                                   settings = MCMC_SETTINGS,
                                                   priors = MCMC_PRIORS,
                                                   spatial = build_s3_spatial(),
                                                   oracle_lambda = dat[[oracle_lambda_field]],
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

    gamma_init <- rep(gamma_value, n1_local)
    state <- initialise_s3_state(dat, settings, priors, spatial, gamma_init = gamma_init)
    state$lambda_tilde <- oracle_lambda
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
        oracle_lambda_source = oracle_lambda_field
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

        ## 2. phi update using marginal NB likelihood, conditional on oracle lambda.
        phi_result <- update_phi_s3(
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
            lambda_tilde = oracle_lambda,
            r = state$r,
            proposal_sd = state$phi_proposal_sd
        )
        state$u <- phi_result$u
        state$phi <- phi_result$phi
        state$xi <- compute_s3_xi(dat$e, dat$x1, dat$x2, state$beta0, state$beta, state$phi)

        ## 3. tau_phi update.
        state$tau_phi <- update_tau_phi(
            phi = state$phi,
            H = spatial$H,
            n1 = n1_local,
            priors = priors
        )

        ## 4. r update using marginal NB likelihood, conditional on oracle lambda.
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
        diagnostics$phi_accept_trace[iter] <- isTRUE(phi_result$accept)
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
        diagnostics$phi_log_alpha[iter] <- phi_result$log_alpha %||% NA_real_
        diagnostics$r_log_alpha[iter, ] <- r_result$diag$log_alpha %||% rep(NA_real_, n1_local)
        diagnostics$kappa_mean_trace[iter] <- kappa_result$diag$mean_kappa %||% mean(state$kappa)

        phi_accept_window <- phi_accept_window + as.integer(isTRUE(phi_result$accept))
        r_accept_window <- r_accept_window + as.integer(r_result$accept)

        if (iter <= n_burnin_local && iter %% adapt_interval == 0L) {
            state$phi_proposal_sd <- adapt_sd_s3(
                current_sd = state$phi_proposal_sd,
                n_accept = phi_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$phi_target_accept %||% 0.25
            )
            state$r_proposal_sd <- adapt_sd_s3(
                current_sd = state$r_proposal_sd,
                n_accept = r_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$r_target_accept %||% 0.30
            )
            phi_accept_window <- 0L
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
            phi_rate <- mean(diagnostics$phi_accept_trace[i0:iter])
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
    diagnostics$phi_accept_rate <- mean(diagnostics$phi_accept_trace)
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
            model = "S4A oracle-lambda fixed-gamma diagnostic",
            fixed_lambda = TRUE,
            oracle_lambda = TRUE,
            oracle_lambda_source = oracle_lambda_field,
            dynamic_lambda = TRUE,
            lambda_sampled_in_fit = FALSE,
            fixed_gamma = TRUE,
            learned_gamma = FALSE,
            gamma_fixed_value = gamma_value,
            updated_blocks = c("beta", "phi", "tau_phi", "r", "kappa"),
            disabled_blocks = c("gamma", "lambda", "delta", "omega"),
            uses_marginal_nb_for_beta_phi_r = TRUE,
            uses_kappa_for_diagnostic_poisson_aug_loglik = TRUE,
            uses_recentered_identified_parameterization = TRUE
        ),
        n_stored = n_stored
    )
}

fit_s4a_oracle_lambda_one_rep <- function(rep_id,
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
    dat_file <- s4a_data_file(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
    chk <- check_s4a_oracle_source_dataset(dat_file)
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

    cat(sprintf("=== S4A oracle-lambda fixed-gamma diagnostic: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Fixed     : gamma = %.3f; lambda = dat$%s\n", fixed_gamma_value_arg, oracle_lambda_field))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa\n")
    cat("Disabled  : lambda, gamma, delta, omega updates\n\n")

    reset_oracle_guards()
    fit <- run_s4a_oracle_lambda_fixed_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        oracle_lambda = oracle_lambda,
        gamma_value = fixed_gamma_value_arg,
        verbose = verbose_arg
    )

    guard_counts <- get_oracle_guards()
    fit$diagnostics$s4a_oracle_numeric_guards <- guard_counts
    fit$diagnostics$s4a_oracle_beta_guard_count <- guard_counts$n_beta_guard
    fit$diagnostics$s4a_oracle_kappa_guard_count <- guard_counts$n_kappa_guard

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
    summary$sparse_beta0_shift <- dat$sparse_beta0_shift %||% NA_real_
    summary$expected_count_multiplier <- dat$expected_count_multiplier %||% NA_real_
    summary$reference_mean_count <- dat$reference_mean_count %||% NA_real_
    summary$reference_zero_prop <- dat$reference_zero_prop %||% NA_real_
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count
    summary$lambda_oracle_min <- chk$lambda_oracle_min
    summary$lambda_oracle_max <- chk$lambda_oracle_max
    summary$s4a_oracle_beta_guard_count <- guard_counts$n_beta_guard
    summary$s4a_oracle_kappa_guard_count <- guard_counts$n_kappa_guard

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id_arg,
            data_scenario_id = data_scenario_id_arg,
            stress_type = dat$stress_type,
            sparse_beta0_shift = dat$sparse_beta0_shift,
            beta0_reference_truth = dat$beta0_reference_truth,
            beta0_sparse_truth = dat$beta0_sparse_truth,
            expected_count_multiplier = dat$expected_count_multiplier,
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir_arg,
            fit_file_prefix = "fit_S4A_oracle_lambda_fixed_gamma_rep",
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root_arg, output_dir_arg, scenario_id_arg)
        ensure_dir(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4A_oracle_lambda_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4A_oracle_lambda_fixed_gamma_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s4a_oracle_lambda_batch <- function(reps = reps_formal,
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
        fit_file <- file.path(out_dir, paste0("fit_S4A_oracle_lambda_fixed_gamma_rep", rr, ".rds"))
        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
        } else {
            fit <- fit_s4a_oracle_lambda_one_rep(
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
        data_file <- s4a_data_file(rep_id, root_arg, data_dir_arg, data_scenario_id_arg)
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

    summary_file <- file.path(out_dir, "summary_S4A_oracle_lambda_fixed_gamma_all_reps.csv")
    manifest_file <- file.path(out_dir, "s4a_oracle_fit_manifest.csv")
    run_info_file <- file.path(out_dir, "run_info_S4A_oracle_lambda_fixed_gamma_diagnostic_T100.rds")

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

    cat("\n=== S4A oracle-lambda diagnostic summary ===\n")
    print(summary_all[, intersect(c(
        "rep_id", "mean_count", "zero_prop", "beta0_mean", "beta1_mean", "beta2_mean",
        "r_mean", "phi_rmse", "phi_cor", "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "s4a_oracle_beta_guard_count", "s4a_oracle_kappa_guard_count"
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

s4a_oracle_summary <- fit_s4a_oracle_lambda_batch(
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

cat("\nS4A oracle-lambda fixed-gamma diagnostic finished successfully.\n")
