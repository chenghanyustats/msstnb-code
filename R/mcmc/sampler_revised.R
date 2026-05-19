## ======================================================================
## sampler_revised.R
## Main MCMC driver for the MSSTNB model.
##
## This version is aligned with the revised update blocks:
##   * beta and phi are updated with kappa collapsed out.
##   * r is updated by the marginal NB likelihood by default.
##   * kappa is updated after r, beta, phi, and xi are refreshed.
##   * lambda is updated by Gamma FFBS conditional on kappa.
##   * beta0, phi, and lambda are deterministically recentered.
##   * delta is updated by the Dirichlet multinomial marginal likelihood.
##   * omega is updated by the Dirichlet FFBS smoother.
##
## Required update files should be sourced before this file:
##   update_regression_revised.R
##   update_icar_revised.R
##   update_kappa.R
##   update_dispersion_revised.R
##   update_gamma_revised.R
##   ffbs_lambda_revised.R
##   update_delta_revised.R
##   smooth_omega_revised.R
## ======================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

is_positive_finite <- function(x) {
    all(is.finite(x)) && all(x > 0)
}

check_same_dim <- function(reference, ..., reference_name = "reference") {
    ref_dim <- dim(reference)
    objects <- list(...)
    object_names <- names(objects)

    for (i in seq_along(objects)) {
        nm <- object_names[i]
        if (is.null(nm) || nm == "") nm <- paste0("object_", i)
        if (!identical(dim(objects[[i]]), ref_dim)) {
            stop(sprintf("%s must have the same dimension as %s.", nm, reference_name))
        }
    }

    invisible(TRUE)
}

validate_mcmc_inputs <- function(dat, settings, priors, spatial, constants) {
    required_dat <- c("y_coarse", "y_fine", "e", "x1", "x2", "TT", "n1", "n_children")
    missing_dat <- setdiff(required_dat, names(dat))
    if (length(missing_dat) > 0L) {
        stop("dat is missing required fields: ", paste(missing_dat, collapse = ", "))
    }

    if (is.null(dim(dat$y_coarse)) || length(dim(dat$y_coarse)) != 2L) {
        stop("dat$y_coarse must be a matrix.")
    }
    if (is.null(dim(dat$y_fine)) || length(dim(dat$y_fine)) != 3L) {
        stop("dat$y_fine must be a 3D array.")
    }

    TT <- nrow(dat$y_coarse)
    n1 <- ncol(dat$y_coarse)
    K <- dim(dat$y_fine)[3]

    if (dat$TT != TT || dat$n1 != n1 || dat$n_children != K) {
        stop("dat$TT, dat$n1, or dat$n_children are inconsistent with y_coarse/y_fine dimensions.")
    }
    if (!identical(dim(dat$y_fine), c(TT, n1, K))) {
        stop("dat$y_fine must have dimension c(TT, n1, n_children).")
    }

    check_same_dim(dat$y_coarse,
                   e = dat$e,
                   x1 = dat$x1,
                   x2 = dat$x2,
                   reference_name = "dat$y_coarse")

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0)) {
        stop("dat$y_coarse must be nonnegative and finite.")
    }
    if (any(!is.finite(dat$y_fine)) || any(dat$y_fine < 0)) {
        stop("dat$y_fine must be nonnegative and finite.")
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0)) {
        stop("dat$e must be positive and finite.")
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.")
    }

    required_settings <- c("n_iter", "n_burnin", "n_thin",
                           "include_covariates", "include_icar", "include_nb")
    missing_settings <- setdiff(required_settings, names(settings))
    if (length(missing_settings) > 0L) {
        stop("settings is missing required fields: ", paste(missing_settings, collapse = ", "))
    }

    if (settings$n_iter <= 0 || settings$n_burnin < 0 || settings$n_thin <= 0) {
        stop("settings$n_iter, settings$n_burnin, and settings$n_thin are invalid.")
    }
    if (settings$n_burnin >= settings$n_iter) {
        stop("settings$n_burnin must be smaller than settings$n_iter.")
    }

    required_constants <- c("A0", "B0", "C0")
    missing_constants <- setdiff(required_constants, names(constants))
    if (length(missing_constants) > 0L) {
        stop("constants is missing required fields: ", paste(missing_constants, collapse = ", "))
    }
    if (!is_positive_finite(constants$A0) || !is_positive_finite(constants$B0)) {
        stop("constants$A0 and constants$B0 must be positive and finite.")
    }
    if (length(constants$C0) != K || !is_positive_finite(constants$C0)) {
        stop("constants$C0 must be a positive finite vector of length n_children.")
    }

    if (settings$include_icar) {
        required_spatial <- c("B_ICAR", "BHB", "H")
        missing_spatial <- setdiff(required_spatial, names(spatial))
        if (length(missing_spatial) > 0L) {
            stop("spatial is missing required fields: ", paste(missing_spatial, collapse = ", "))
        }
        if (!identical(dim(spatial$B_ICAR), c(n1, n1 - 1L))) {
            stop("spatial$B_ICAR must have dimension c(n1, n1 - 1).")
        }
        if (!identical(dim(spatial$BHB), c(n1 - 1L, n1 - 1L))) {
            stop("spatial$BHB must have dimension c(n1 - 1, n1 - 1).")
        }
        if (!identical(dim(spatial$H), c(n1, n1))) {
            stop("spatial$H must have dimension c(n1, n1).")
        }
    }

    required_priors <- c("beta0_mean", "beta0_sd", "beta_mean", "beta_sd",
                         "gamma_a", "gamma_b", "delta_a", "delta_b")
    missing_priors <- setdiff(required_priors, names(priors))
    if (length(missing_priors) > 0L) {
        stop("priors is missing required fields: ", paste(missing_priors, collapse = ", "))
    }

    invisible(TRUE)
}

adapt_sd <- function(current_sd, n_accept, n_trials, target_rate = 0.30,
                     increase = 1.10, decrease = 0.90,
                     min_sd = 1e-4, max_sd = 5.0) {
    if (length(n_accept) == 1L && length(current_sd) > 1L) {
        n_accept <- rep(n_accept, length(current_sd))
    }
    rate <- n_accept / n_trials
    new_sd <- current_sd
    new_sd[rate > target_rate + 0.05] <- new_sd[rate > target_rate + 0.05] * increase
    new_sd[rate < target_rate - 0.05] <- new_sd[rate < target_rate - 0.05] * decrease
    pmin(pmax(new_sd, min_sd), max_sd)
}

method_label <- function(settings) {
    r_method <- settings$r_update_method %||% "marginal_nb"
    paste0("MSSTNB revised sampler; r_update=", r_method,
           "; beta_phi=kappa_collapsed")
}

compute_xi_mcmc <- function(e, x1, x2, beta0, beta, phi) {
    check_same_dim(e, x1 = x1, x2 = x2, reference_name = "e")

    if (length(beta) < 2L) {
        stop("beta must contain at least beta[1] and beta[2].")
    }
    if (length(phi) != ncol(e)) {
        stop("length(phi) must equal ncol(e).")
    }
    if (!is.finite(beta0) || any(!is.finite(beta)) || any(!is.finite(phi))) {
        stop("beta0, beta, and phi must be finite.")
    }

    TT <- nrow(e)
    n1 <- ncol(e)
    xi <- matrix(NA_real_, TT, n1)

    for (j in seq_len(n1)) {
        linpred_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        linpred_j <- pmin(pmax(linpred_j, -700), 700)
        xi[, j] <- e[, j] * exp(linpred_j)
    }

    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("compute_xi_mcmc produced nonpositive or nonfinite xi.")
    }

    xi
}

compute_loglik <- function(y_coarse, xi, lambda_tilde, kappa) {
    check_same_dim(y_coarse,
                   xi = xi,
                   lambda_tilde = lambda_tilde,
                   kappa = kappa,
                   reference_name = "y_coarse")
    mu <- pmax(xi * lambda_tilde * kappa, .Machine$double.xmin)
    sum(dpois(y_coarse, lambda = mu, log = TRUE))
}

compute_u_from_phi <- function(phi, B) {
    BtB <- crossprod(B)
    rhs <- crossprod(B, phi)
    as.numeric(solve(BtB, rhs))
}

initialise_state <- function(dat, settings, priors, spatial, constants) {
    TT <- dat$TT
    n1 <- dat$n1
    K <- dat$n_children
    p <- length(priors$beta_mean)

    if (p != 2L) {
        stop("This sampler currently expects exactly two covariates, so length(priors$beta_mean) must be 2.")
    }

    if (settings$include_covariates) {
        ess_base <- fit_glm_for_ess(dat)
        beta0 <- ess_base$center[1]
        beta <- ess_base$center[2:(p + 1L)]
    } else {
        ess_base <- NULL
        beta0 <- 0
        beta <- rep(0, p)
    }

    phi <- rep(0, n1)
    u <- rep(0, max(n1 - 1L, 0L))
    tau_phi <- 1

    if (settings$include_nb) {
        r <- rep(settings$r_init %||% 1, n1)
        kappa <- matrix(1, TT, n1)
    } else {
        r <- rep(Inf, n1)
        kappa <- matrix(1, TT, n1)
    }

    gamma_init <- settings$gamma_init %||% 0.9
    delta_init <- settings$delta_init %||% 0.9
    if (gamma_init <= 0 || gamma_init >= 1 || delta_init <= 0 || delta_init >= 1) {
        stop("gamma_init and delta_init must lie strictly in (0, 1).")
    }

    gamma_vec <- rep(gamma_init, n1)
    delta <- delta_init
    lambda_tilde <- matrix(1, TT, n1)

    omega <- array(NA_real_, dim = c(TT, n1, K))
    for (j in seq_len(n1)) {
        for (t in seq_len(TT)) {
            y_kids <- dat$y_fine[t, j, ]
            total <- sum(y_kids)
            omega[t, j, ] <- if (total > 0) {
                (y_kids + 0.1) / (total + 0.1 * K)
            } else {
                rep(1 / K, K)
            }
        }
    }

    xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2, beta0, beta, phi)

    list(
        beta0 = beta0,
        beta = beta,
        phi = phi,
        u = u,
        tau_phi = tau_phi,
        r = r,
        kappa = kappa,
        gamma = gamma_vec,
        delta = delta,
        lambda_tilde = lambda_tilde,
        omega = omega,
        xi = xi,
        ess_base = ess_base,
        phi_proposal_sd = rep(settings$phi_proposal_sd_init %||% 0.01, max(n1 - 1L, 0L)),
        gamma_proposal_sd = rep(settings$gamma_proposal_sd_init %||% 0.30, n1),
        delta_proposal_sd = settings$delta_proposal_sd_init %||% 0.30,
        r_proposal_sd = rep(settings$r_proposal_sd_init %||% 0.50, n1)
    )
}

run_one_iteration <- function(state, dat, settings, priors, spatial, constants) {
    diag <- list()
    r_update_method <- settings$r_update_method %||% "marginal_nb"

    ## 1. beta update using kappa collapsed marginal NB likelihood.
    if (settings$include_covariates) {
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
        diag$beta_n_reject <- beta_result$n_reject
        diag$beta_log_target <- beta_result$log_target %||% beta_result$log_lik %||% NA_real_
        diag$beta_log_marginal_lik <- beta_result$log_marginal_lik %||% NA_real_
        diag$beta_ess_mode <- beta_result$ess_mode %||% NA_character_
    } else {
        diag$beta_n_reject <- NA_real_
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## 2. phi update using kappa collapsed marginal NB likelihood.
    if (settings$include_icar) {
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
        diag$phi_accept <- phi_result$accept
        diag$phi_log_alpha <- phi_result$log_alpha %||% NA_real_
        diag$phi_log_posterior <- phi_result$log_posterior %||% NA_real_
    } else {
        diag$phi_accept <- NA
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## 3. tau_phi update.
    if (settings$include_icar) {
        state$tau_phi <- update_tau_phi(state$phi, spatial$H, dat$n1, priors)
    }

    ## 4. r update.  The default is marginal NB, so r is updated before kappa.
    if (settings$include_nb) {
        if (r_update_method == "marginal_nb") {
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
        } else if (r_update_method == "conditional_kappa") {
            r_result <- update_r(
                r_current = state$r,
                kappa = state$kappa,
                priors = priors,
                mh_sd = state$r_proposal_sd,
                method = "conditional_kappa",
                return_diag = TRUE
            )
        } else {
            stop("settings$r_update_method must be 'marginal_nb' or 'conditional_kappa'.")
        }
        state$r <- r_result$r
        diag$r_accept <- r_result$accept
        diag$r_log_alpha <- r_result$diag$log_alpha %||% rep(NA_real_, dat$n1)
        diag$r_method <- r_result$diag$method %||% r_update_method
    } else {
        diag$r_accept <- rep(FALSE, dat$n1)
    }

    ## 5. kappa exact conditional update, given refreshed beta, phi, r, xi, and lambda.
    if (settings$include_nb) {
        kappa_result <- update_kappa(
            y_coarse = dat$y_coarse,
            lambda_tilde = state$lambda_tilde,
            xi = state$xi,
            r = state$r,
            return_diag = TRUE
        )
        if (is.list(kappa_result) && !is.null(kappa_result$kappa)) {
            state$kappa <- kappa_result$kappa
            diag$kappa <- kappa_result$diag
        } else {
            state$kappa <- kappa_result
            diag$kappa <- NULL
        }
    } else {
        state$kappa <- matrix(1, dat$TT, dat$n1)
    }

    ## 6. gamma update.  This is lambda collapsed but kappa conditioned.
    gamma_result <- update_gamma(
        gamma_current = state$gamma,
        y_coarse = dat$y_coarse,
        xi = state$xi,
        kappa = state$kappa,
        a0 = constants$A0,
        b0 = constants$B0,
        priors = priors,
        mh_sd = state$gamma_proposal_sd,
        return_diag = TRUE
    )
    state$gamma <- gamma_result$gamma
    diag$gamma_accept <- gamma_result$accept
    diag$gamma_log_alpha <- gamma_result$diag$log_alpha %||% rep(NA_real_, dat$n1)

    ## 7. lambda FFBS conditional on kappa.
    lambda_result <- ffbs_lambda_all(
        gamma = state$gamma,
        y_coarse = dat$y_coarse,
        xi = state$xi,
        kappa = state$kappa,
        a0 = constants$A0,
        b0 = constants$B0,
        return_diag = TRUE
    )
    if (is.list(lambda_result) && !is.null(lambda_result$lambda_tilde)) {
        state$lambda_tilde <- lambda_result$lambda_tilde
        diag$lambda <- lambda_result$diag
    } else {
        state$lambda_tilde <- lambda_result
        diag$lambda <- NULL
    }

    ## 8. Deterministic recentering for beta0, phi, and lambda.
    if (settings$include_icar) {
        rc <- recenter(
            beta0 = state$beta0,
            phi = state$phi,
            lambda_tilde = state$lambda_tilde,
            return_diag = TRUE
        )
        state$beta0 <- rc$beta0
        state$phi <- rc$phi
        state$lambda_tilde <- rc$lambda_tilde
        state$u <- compute_u_from_phi(state$phi, spatial$B_ICAR)
        diag$recenter <- rc$diag
    } else if (settings$include_covariates) {
        lt_safe <- pmax(state$lambda_tilde, .Machine$double.xmin)
        s_global <- mean(log(lt_safe))
        state$lambda_tilde <- state$lambda_tilde * exp(-s_global)
        state$beta0 <- state$beta0 + s_global
        diag$recenter <- list(global_shift = s_global)
    } else {
        diag$recenter <- NULL
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## 9. delta update using omega collapsed Dirichlet multinomial likelihood.
    delta_result <- update_delta(
        delta_current = state$delta,
        y_fine = dat$y_fine,
        c0 = constants$C0,
        priors = priors,
        mh_sd = state$delta_proposal_sd,
        return_diag = TRUE
    )
    state$delta <- delta_result$delta
    diag$delta_accept <- delta_result$accept
    diag$delta_log_alpha <- delta_result$diag$log_alpha %||% NA_real_

    ## 10. omega FFBS smoother.
    omega_result <- smooth_omega_all(
        delta = state$delta,
        y_fine = dat$y_fine,
        c0 = constants$C0,
        return_diag = TRUE
    )
    if (is.list(omega_result) && !is.null(omega_result$omega)) {
        state$omega <- omega_result$omega
        diag$omega <- omega_result$diag
    } else {
        state$omega <- omega_result
        diag$omega <- NULL
    }

    ## 11. Complete data Poisson log likelihood for monitoring only.
    diag$loglik <- compute_loglik(dat$y_coarse, state$xi,
                                  state$lambda_tilde, state$kappa)

    list(state = state, diag = diag)
}

run_mcmc <- function(dat, settings, priors, spatial, constants,
                     verbose = 1000L) {
    validate_mcmc_inputs(dat, settings, priors, spatial, constants)

    n_iter <- settings$n_iter
    n_burnin <- settings$n_burnin
    n_thin <- settings$n_thin
    n_stored <- as.integer(floor((n_iter - n_burnin) / n_thin))

    TT <- dat$TT
    n1 <- dat$n1
    K <- dat$n_children
    p <- length(priors$beta_mean)

    adapt_interval <- settings$adapt_interval %||% 50L

    samples <- list(
        beta0 = numeric(n_stored),
        beta = matrix(NA_real_, n_stored, p),
        phi = matrix(NA_real_, n_stored, n1),
        tau_phi = numeric(n_stored),
        r = matrix(NA_real_, n_stored, n1),
        gamma = matrix(NA_real_, n_stored, n1),
        delta = numeric(n_stored),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT, n1)),
        kappa = array(NA_real_, dim = c(n_stored, TT, n1)),
        omega = array(NA_real_, dim = c(n_stored, TT, n1, K)),
        loglik = numeric(n_stored)
    )

    diag_all <- list(
        loglik_trace = numeric(n_iter),
        beta_n_reject = rep(NA_real_, n_iter),
        beta_log_target = rep(NA_real_, n_iter),
        beta_log_marginal_lik = rep(NA_real_, n_iter),
        phi_accept_trace = rep(NA, n_iter),
        phi_log_alpha = rep(NA_real_, n_iter),
        gamma_accept_trace = matrix(FALSE, n_iter, n1),
        gamma_log_alpha = matrix(NA_real_, n_iter, n1),
        delta_accept_trace = rep(NA, n_iter),
        delta_log_alpha = rep(NA_real_, n_iter),
        r_accept_trace = matrix(FALSE, n_iter, n1),
        r_log_alpha = matrix(NA_real_, n_iter, n1),
        lambda_min_trace = rep(NA_real_, n_iter),
        lambda_max_trace = rep(NA_real_, n_iter),
        omega_row_sum_error = rep(NA_real_, n_iter),
        recenter_error = rep(NA_real_, n_iter)
    )

    state <- initialise_state(dat, settings, priors, spatial, constants)

    phi_accept_window <- 0L
    gamma_accept_window <- rep(0L, n1)
    delta_accept_window <- 0L
    r_accept_window <- rep(0L, n1)

    start_time <- proc.time()
    store_idx <- 0L

    for (iter in seq_len(n_iter)) {
        result <- run_one_iteration(state, dat, settings, priors, spatial, constants)
        state <- result$state
        diag <- result$diag

        diag_all$loglik_trace[iter] <- diag$loglik
        diag_all$beta_n_reject[iter] <- diag$beta_n_reject %||% NA_real_
        diag_all$beta_log_target[iter] <- diag$beta_log_target %||% NA_real_
        diag_all$beta_log_marginal_lik[iter] <- diag$beta_log_marginal_lik %||% NA_real_

        if (!is.null(diag$phi_accept) && !is.na(diag$phi_accept)) {
            diag_all$phi_accept_trace[iter] <- diag$phi_accept
            if (isTRUE(diag$phi_accept)) phi_accept_window <- phi_accept_window + 1L
        }
        diag_all$phi_log_alpha[iter] <- diag$phi_log_alpha %||% NA_real_

        if (!is.null(diag$gamma_accept)) {
            diag_all$gamma_accept_trace[iter, ] <- diag$gamma_accept
            gamma_accept_window <- gamma_accept_window + as.integer(diag$gamma_accept)
        }
        if (!is.null(diag$gamma_log_alpha)) {
            diag_all$gamma_log_alpha[iter, ] <- diag$gamma_log_alpha
        }

        if (!is.null(diag$delta_accept) && !is.na(diag$delta_accept)) {
            diag_all$delta_accept_trace[iter] <- diag$delta_accept
            if (isTRUE(diag$delta_accept)) delta_accept_window <- delta_accept_window + 1L
        }
        diag_all$delta_log_alpha[iter] <- diag$delta_log_alpha %||% NA_real_

        if (!is.null(diag$r_accept)) {
            diag_all$r_accept_trace[iter, ] <- diag$r_accept
            r_accept_window <- r_accept_window + as.integer(diag$r_accept)
        }
        if (!is.null(diag$r_log_alpha)) {
            diag_all$r_log_alpha[iter, ] <- diag$r_log_alpha
        }

        if (!is.null(diag$lambda)) {
            diag_all$lambda_min_trace[iter] <- diag$lambda$min_lambda %||% NA_real_
            diag_all$lambda_max_trace[iter] <- diag$lambda$max_lambda %||% NA_real_
        }
        if (!is.null(diag$omega)) {
            diag_all$omega_row_sum_error[iter] <- diag$omega$max_row_sum_error %||% NA_real_
        }
        if (!is.null(diag$recenter)) {
            diag_all$recenter_error[iter] <- diag$recenter$max_abs_log_core_difference %||% NA_real_
        }

        if (iter <= n_burnin && iter %% adapt_interval == 0L) {
            if (settings$include_icar) {
                state$phi_proposal_sd <- adapt_sd(
                    current_sd = state$phi_proposal_sd,
                    n_accept = phi_accept_window,
                    n_trials = adapt_interval,
                    target_rate = settings$phi_target_accept %||% 0.25
                )
                phi_accept_window <- 0L
            }

            state$gamma_proposal_sd <- adapt_sd(
                current_sd = state$gamma_proposal_sd,
                n_accept = gamma_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$gamma_target_accept %||% 0.30
            )
            gamma_accept_window <- rep(0L, n1)

            state$delta_proposal_sd <- adapt_sd(
                current_sd = state$delta_proposal_sd,
                n_accept = delta_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$delta_target_accept %||% 0.30
            )
            delta_accept_window <- 0L

            if (settings$include_nb) {
                state$r_proposal_sd <- adapt_sd(
                    current_sd = state$r_proposal_sd,
                    n_accept = r_accept_window,
                    n_trials = adapt_interval,
                    target_rate = settings$r_target_accept %||% 0.30
                )
                r_accept_window <- rep(0L, n1)
            }
        }

        if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
            store_idx <- store_idx + 1L
            samples$beta0[store_idx] <- state$beta0
            samples$beta[store_idx, ] <- state$beta
            samples$phi[store_idx, ] <- state$phi
            samples$tau_phi[store_idx] <- state$tau_phi
            samples$r[store_idx, ] <- state$r
            samples$gamma[store_idx, ] <- state$gamma
            samples$delta[store_idx] <- state$delta
            samples$lambda_tilde[store_idx, , ] <- state$lambda_tilde
            samples$kappa[store_idx, , ] <- state$kappa
            samples$omega[store_idx, , , ] <- state$omega
            samples$loglik[store_idx] <- diag$loglik
        }

        if (verbose > 0 && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1L, iter - 99L)
            phi_rate <- if (settings$include_icar) {
                mean(diag_all$phi_accept_trace[i0:iter], na.rm = TRUE)
            } else {
                NA_real_
            }
            gamma_rate <- mean(diag_all$gamma_accept_trace[i0:iter, , drop = FALSE])
            r_rate <- if (settings$include_nb) {
                mean(diag_all$r_accept_trace[i0:iter, , drop = FALSE])
            } else {
                NA_real_
            }
            beta_rej <- mean(diag_all$beta_n_reject[i0:iter], na.rm = TRUE)

            cat(sprintf(
                "  iter %5d/%d [%.0fs] loglik=%.1f beta0=%.3f gamma_mean=%.3f delta=%.3f | phi_acc=%.2f gamma_acc=%.2f r_acc=%.2f beta_rej=%.1f\n",
                iter, n_iter, elapsed, diag$loglik,
                state$beta0, mean(state$gamma), state$delta,
                phi_rate, gamma_rate, r_rate, beta_rej
            ))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]

    diag_all$elapsed_sec <- elapsed_total
    diag_all$phi_accept_rate <- if (settings$include_icar) {
        mean(diag_all$phi_accept_trace, na.rm = TRUE)
    } else {
        NA_real_
    }
    diag_all$gamma_accept_rate <- colMeans(diag_all$gamma_accept_trace)
    diag_all$delta_accept_rate <- mean(diag_all$delta_accept_trace, na.rm = TRUE)
    diag_all$r_accept_rate <- if (settings$include_nb) {
        colMeans(diag_all$r_accept_trace)
    } else {
        rep(NA_real_, n1)
    }
    diag_all$beta_mean_n_reject <- mean(diag_all$beta_n_reject, na.rm = TRUE)
    diag_all$phi_proposal_sd_final <- state$phi_proposal_sd
    diag_all$gamma_proposal_sd_final <- state$gamma_proposal_sd
    diag_all$delta_proposal_sd_final <- state$delta_proposal_sd
    diag_all$r_proposal_sd_final <- state$r_proposal_sd
    diag_all$r_update_method <- settings$r_update_method %||% "marginal_nb"

    if (verbose > 0) {
        cat(sprintf("\nDone. %d iterations in %.1f sec (%.1f iter/sec)\n",
                    n_iter, elapsed_total, n_iter / elapsed_total))
        cat(sprintf("Stored %d post burn in samples (thin=%d)\n", n_stored, n_thin))
        cat(sprintf("Acceptance rates: phi=%.2f gamma=%.2f delta=%.2f r=%.2f\n",
                    diag_all$phi_accept_rate,
                    mean(diag_all$gamma_accept_rate),
                    diag_all$delta_accept_rate,
                    mean(diag_all$r_accept_rate, na.rm = TRUE)))
        cat(sprintf("Beta ESS: %.1f mean rejections per step\n",
                    diag_all$beta_mean_n_reject))
    }

    list(
        samples = samples,
        diagnostics = diag_all,
        final_state = state,
        settings = settings,
        priors = priors,
        method = method_label(settings),
        n_stored = n_stored
    )
}
