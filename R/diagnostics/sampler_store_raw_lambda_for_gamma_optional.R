## ==========================================================================
## sampler.R
## Kappa-collapsed ordinary sampler with storage of the raw lambda path used
## for gamma-focused diagnostics.
##
## New in this version:
##   - stores pre-recentering raw lambda_{0:T} at each saved draw
##   - saved as samples$lambda_raw_for_gamma
##   - fit$gamma_conditional_diag records how this raw path was defined
## ==========================================================================

compute_xi_mcmc <- function(e, x1, x2, beta0, beta, phi) {
    TT <- nrow(e); n1 <- ncol(e)
    xi <- matrix(NA_real_, TT, n1)
    for (j in seq_len(n1)) {
        linpred_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        xi[, j] <- e[, j] * exp(linpred_j)
    }
    xi
}

compute_loglik <- function(y_coarse, xi, lambda_tilde, kappa) {
    mu <- pmax(xi * lambda_tilde * kappa, .Machine$double.xmin)
    sum(dpois(y_coarse, lambda = mu, log = TRUE))
}

initialise_state <- function(dat, settings, priors, spatial, constants) {
    TT <- dat$TT; n1 <- dat$n1; K <- dat$n_children
    p <- length(priors$beta_mean)

    if (settings$include_covariates) {
        ess_base <- fit_glm_for_ess(dat)
        beta0 <- ess_base$center[1]
        beta  <- ess_base$center[2:(p + 1)]
    } else {
        ess_base <- NULL
        beta0 <- 0
        beta  <- rep(0, p)
    }

    phi <- rep(0, n1)
    u <- rep(0, n1 - 1)
    tau_phi <- 1

    if (settings$include_nb) {
        r <- rep(1, n1)
        kappa <- matrix(1, TT, n1)
    } else {
        r <- rep(Inf, n1)
        kappa <- matrix(1, TT, n1)
    }

    gamma_vec <- rep(0.9, n1)
    delta <- 0.9
    lambda_tilde <- matrix(1, TT, n1)
    lambda_raw_for_gamma <- rbind(rep(1, n1), lambda_tilde)

    omega <- array(NA_real_, dim = c(TT, n1, K))
    for (j in seq_len(n1)) {
        for (t in seq_len(TT)) {
            y_kids <- dat$y_fine[t, j, ]
            total <- sum(y_kids)
            omega[t, j, ] <- if (total > 0) (y_kids + 0.1) / (total + 0.1 * K) else rep(1 / K, K)
        }
    }

    xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2, beta0, beta, phi)

    phi_proposal_sd   <- rep(0.01, n1 - 1)
    gamma_proposal_sd <- rep(0.3, n1)
    delta_proposal_sd <- 0.3
    r_proposal_sd     <- rep(0.5, n1)

    list(beta0 = beta0, beta = beta, phi = phi, u = u,
         tau_phi = tau_phi, r = r, kappa = kappa,
         gamma = gamma_vec, delta = delta,
         lambda_tilde = lambda_tilde,
         lambda_raw_for_gamma = lambda_raw_for_gamma,
         omega = omega, xi = xi,
         ess_base = ess_base,
         phi_proposal_sd = phi_proposal_sd,
         gamma_proposal_sd = gamma_proposal_sd,
         delta_proposal_sd = delta_proposal_sd,
         r_proposal_sd = r_proposal_sd)
}

run_one_iteration <- function(state, dat, settings, priors, spatial, constants) {
    diag <- list()

    ## Step 3: beta, collapsed over kappa
    if (settings$include_covariates) {
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            lambda_tilde = state$lambda_tilde,
            phi = state$phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r
        )
        state$beta0 <- beta_result$sample[1]
        state$beta  <- beta_result$sample[2:3]
        diag$beta_n_reject <- beta_result$n_reject
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## Step 4: phi, collapsed over kappa
    if (settings$include_icar) {
        phi_result <- update_phi(
            u_current = state$u, B = spatial$B_ICAR,
            BHB = spatial$BHB, tau_phi = state$tau_phi,
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            beta0 = state$beta0, beta = state$beta,
            lambda_tilde = state$lambda_tilde,
            r = state$r,
            proposal_sd = state$phi_proposal_sd
        )
        state$u   <- phi_result$u
        state$phi <- phi_result$phi
        diag$phi_accept <- phi_result$accept
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## Step 5: tau_phi
    if (settings$include_icar) {
        state$tau_phi <- update_tau_phi(state$phi, spatial$H, dat$n1, priors)
    }

    ## Step 2 (moved here): kappa exact conditional
    if (settings$include_nb) {
        state$kappa <- update_kappa(dat$y_coarse, state$lambda_tilde,
                                    state$xi, state$r)
    }

    ## Step 6: r
    if (settings$include_nb) {
        r_result <- update_r(state$r, state$kappa, priors,
                             mh_sd = state$r_proposal_sd)
        state$r <- r_result$r
        diag$r_accept <- r_result$accept
    }

    ## Step 7: gamma
    gamma_result <- update_gamma(
        state$gamma, dat$y_coarse, state$xi, state$kappa,
        a0 = constants$A0, b0 = constants$B0,
        priors = priors, mh_sd = state$gamma_proposal_sd
    )
    state$gamma <- gamma_result$gamma
    diag$gamma_accept <- gamma_result$accept

    ## Step 8: lambda FFBS, capture raw 0:T path before re-centering
    ffbs_out <- ffbs_lambda_all_with_raw0(
        state$gamma, dat$y_coarse, state$xi, state$kappa,
        a0 = constants$A0, b0 = constants$B0
    )
    state$lambda_tilde <- ffbs_out$lambda_tilde
    state$lambda_raw_for_gamma <- ffbs_out$lambda_raw_0T
    diag$lambda_raw_for_gamma <- ffbs_out$lambda_raw_0T

    ## Step 9: re-centering
    if (settings$include_icar) {
        rc <- recenter(state$beta0, state$phi, state$lambda_tilde)
        state$beta0 <- rc$beta0
        state$phi   <- rc$phi
        state$lambda_tilde <- rc$lambda_tilde
        state$u <- as.numeric(t(spatial$B_ICAR) %*% state$phi)
    } else if (settings$include_covariates) {
        lt_safe <- pmax(state$lambda_tilde, .Machine$double.xmin)
        s_global <- mean(log(lt_safe))
        state$lambda_tilde <- state$lambda_tilde * exp(-s_global)
        state$beta0 <- state$beta0 + s_global
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## Step 10: delta
    delta_result <- update_delta(
        state$delta, dat$y_fine, c0 = constants$C0,
        priors = priors, mh_sd = state$delta_proposal_sd
    )
    state$delta <- delta_result$delta
    diag$delta_accept <- delta_result$accept

    ## Step 11: omega
    state$omega <- smooth_omega_all(state$delta, dat$y_fine, c0 = constants$C0)

    diag$loglik <- compute_loglik(dat$y_coarse, state$xi,
                                  state$lambda_tilde, state$kappa)

    list(state = state, diag = diag)
}

run_mcmc <- function(dat, settings, priors, spatial, constants,
                     verbose = 1000L) {
    n_iter <- settings$n_iter
    n_burnin <- settings$n_burnin
    n_thin <- settings$n_thin
    n_stored <- as.integer((n_iter - n_burnin) / n_thin)

    TT <- dat$TT
    n1 <- dat$n1
    K <- dat$n_children
    p <- length(priors$beta_mean)

    adapt_interval <- 50L

    samples <- list(
        beta0        = numeric(n_stored),
        beta         = matrix(NA_real_, n_stored, p),
        phi          = matrix(NA_real_, n_stored, n1),
        tau_phi      = numeric(n_stored),
        r            = matrix(NA_real_, n_stored, n1),
        gamma        = matrix(NA_real_, n_stored, n1),
        delta        = numeric(n_stored),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT, n1)),
        lambda_raw_for_gamma = array(NA_real_, dim = c(n_stored, TT + 1L, n1)),
        kappa        = array(NA_real_, dim = c(n_stored, TT, n1)),
        omega        = array(NA_real_, dim = c(n_stored, TT, n1, K)),
        loglik       = numeric(n_stored)
    )

    diag_all <- list(
        loglik_trace         = numeric(n_iter),
        beta_n_reject        = numeric(n_iter),
        phi_accept_trace     = logical(n_iter),
        gamma_accept_trace   = matrix(FALSE, n_iter, n1),
        delta_accept_trace   = logical(n_iter),
        r_accept_trace       = matrix(FALSE, n_iter, n1)
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
        if (!is.null(diag$beta_n_reject)) diag_all$beta_n_reject[iter] <- diag$beta_n_reject
        if (!is.null(diag$phi_accept)) {
            diag_all$phi_accept_trace[iter] <- diag$phi_accept
            if (diag$phi_accept) phi_accept_window <- phi_accept_window + 1L
        }
        if (!is.null(diag$gamma_accept)) {
            diag_all$gamma_accept_trace[iter, ] <- diag$gamma_accept
            gamma_accept_window <- gamma_accept_window + diag$gamma_accept
        }
        if (!is.null(diag$delta_accept)) {
            diag_all$delta_accept_trace[iter] <- diag$delta_accept
            if (diag$delta_accept) delta_accept_window <- delta_accept_window + 1L
        }
        if (!is.null(diag$r_accept)) {
            diag_all$r_accept_trace[iter, ] <- diag$r_accept
            r_accept_window <- r_accept_window + diag$r_accept
        }

        if (iter <= n_burnin && iter %% adapt_interval == 0L) {
            if (settings$include_icar) {
                state$phi_proposal_sd <- adapt_sd(state$phi_proposal_sd, phi_accept_window,
                                                  adapt_interval, target_rate = 0.25)
                phi_accept_window <- 0L
            }
            state$gamma_proposal_sd <- adapt_sd(state$gamma_proposal_sd, gamma_accept_window,
                                                adapt_interval, target_rate = 0.30)
            gamma_accept_window <- rep(0L, n1)
            state$delta_proposal_sd <- adapt_sd(state$delta_proposal_sd, delta_accept_window,
                                                adapt_interval, target_rate = 0.30)
            delta_accept_window <- 0L
            if (settings$include_nb) {
                state$r_proposal_sd <- adapt_sd(state$r_proposal_sd, r_accept_window,
                                                adapt_interval, target_rate = 0.30)
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
            samples$lambda_raw_for_gamma[store_idx, , ] <- state$lambda_raw_for_gamma
            samples$kappa[store_idx, , ] <- state$kappa
            samples$omega[store_idx, , , ] <- state$omega
            samples$loglik[store_idx] <- diag$loglik
        }

        if (verbose > 0 && iter %% verbose == 0) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1, iter - 99)
            phi_rate <- mean(diag_all$phi_accept_trace[i0:iter])
            gamma_rate <- mean(diag_all$gamma_accept_trace[i0:iter, 1])
            beta_rej <- mean(diag_all$beta_n_reject[i0:iter])
            cat(sprintf(
                "  iter %5d/%d [%.0fs] loglik=%.1f beta0=%.3f gamma1=%.3f | phi_acc=%.2f gamma_acc=%.2f beta_rej=%.1f\n",
                iter, n_iter, elapsed, diag$loglik,
                state$beta0, state$gamma[1],
                phi_rate, gamma_rate, beta_rej))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]

    diag_all$elapsed_sec         <- elapsed_total
    diag_all$phi_accept_rate     <- mean(diag_all$phi_accept_trace)
    diag_all$gamma_accept_rate   <- colMeans(diag_all$gamma_accept_trace)
    diag_all$delta_accept_rate   <- mean(diag_all$delta_accept_trace)
    diag_all$r_accept_rate       <- colMeans(diag_all$r_accept_trace)
    diag_all$beta_mean_n_reject  <- mean(diag_all$beta_n_reject)
    diag_all$phi_proposal_sd_final   <- state$phi_proposal_sd
    diag_all$gamma_proposal_sd_final <- state$gamma_proposal_sd
    diag_all$delta_proposal_sd_final <- state$delta_proposal_sd
    diag_all$r_proposal_sd_final     <- state$r_proposal_sd

    gamma_grid <- seq(0.01, 0.999, length.out = 200L)

    if (verbose > 0) {
        cat(sprintf("\nDone. %d iterations in %.1f sec (%.1f iter/sec)\n",
                    n_iter, elapsed_total, n_iter / elapsed_total))
        cat(sprintf("Stored %d post-burn-in samples (thin=%d)\n", n_stored, n_thin))
        cat(sprintf("Acceptance rates: phi=%.2f gamma=%.2f delta=%.2f r=%.2f\n",
                    diag_all$phi_accept_rate,
                    mean(diag_all$gamma_accept_rate),
                    diag_all$delta_accept_rate,
                    mean(diag_all$r_accept_rate)))
        cat(sprintf("Beta ESS: %.1f mean rejections/step\n",
                    diag_all$beta_mean_n_reject))
    }

    list(
        samples = samples,
        diagnostics = diag_all,
        gamma_conditional_diag = list(
            grid = gamma_grid,
            lambda_raw_source = "samples$lambda_raw_for_gamma",
            lambda_raw_definition = "FFBS draw before re-centering; dimension draws x (TT+1) x n1; includes sampled raw lambda_0 via one extra backward step from (a0,b0).",
            pre_recenter = TRUE
        ),
        settings = settings,
        priors = priors,
        method = method_label(settings),
        n_stored = n_stored
    )
}
