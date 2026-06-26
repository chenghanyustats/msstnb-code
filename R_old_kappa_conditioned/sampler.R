## ==========================================================================
## sampler.R
## Main MCMC sampler for the MSSTNB model
##
## Implements Algorithm 3 (one iteration of the partially collapsed blocked
## sampler) and the outer MCMC loop with storage, diagnostics, and progress.
## ==========================================================================


#' Compute the effective offset ξ_{t,j} = e·exp(β₀ + x'β + φ_j)
#'
#' @param e, x1, x2  T × n1 input matrices
#' @param beta0       Scalar intercept
#' @param beta        Length-p regression vector
#' @param phi         Length-n1 ICAR vector
#' @return            T × n1 matrix of ξ values
compute_xi_mcmc <- function(e, x1, x2, beta0, beta, phi) {
    TT <- nrow(e)
    n1 <- ncol(e)
    xi <- matrix(NA_real_, TT, n1)
    for (j in seq_len(n1)) {
        linpred_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        xi[, j] <- e[, j] * exp(linpred_j)
    }
    return(xi)
}


#' Compute the Poisson log-likelihood for the full model
#'
#' log L = Σ_{t,j} [ y log(μ) - μ - log(y!) ]
#' where μ = ξ · λ̃ · κ
#'
#' @param y_coarse     T × n1 count matrix
#' @param xi           T × n1 effective offset
#' @param lambda_tilde T × n1 residual risk
#' @param kappa        T × n1 NB random effect
#' @return             Scalar log-likelihood
compute_loglik <- function(y_coarse, xi, lambda_tilde, kappa) {
    mu <- xi * lambda_tilde * kappa
    sum(dpois(y_coarse, lambda = mu, log = TRUE))
}


#' Initialise MCMC state from data and priors
#'
#' @param dat       Simulated dataset (from simulate_one_dataset)
#' @param settings  MCMC_SETTINGS list
#' @param priors    MCMC_PRIORS list
#' @param spatial   List with H, B_ICAR, R_BHB (precomputed)
#' @param constants List with A0, B0, C0
#' @return          Named list of initial parameter values
initialise_state <- function(dat, settings, priors, spatial, constants) {

    TT <- dat$TT
    n1 <- dat$n1
    K  <- dat$n_children

    ## ---- regression ----
    if (settings$include_covariates) {
        beta0 <- 0
        beta  <- c(0, 0)
    } else {
        beta0 <- 0
        beta  <- c(0, 0)  # fixed at 0 throughout
    }

    ## ---- ICAR ----
    if (settings$include_icar) {
        phi     <- rep(0, n1)            # start at zero
        u       <- rep(0, n1 - 1)        # reduced coordinates
        tau_phi <- 1
    } else {
        phi     <- rep(0, n1)            # fixed at 0
        u       <- rep(0, n1 - 1)
        tau_phi <- 1                     # not updated
    }

    ## ---- NB dispersion ----
    if (settings$include_nb) {
        r     <- rep(1, n1)              # start moderate
        kappa <- matrix(1, TT, n1)       # start at prior mean
    } else {
        r     <- rep(Inf, n1)            # fixed: no overdispersion
        kappa <- matrix(1, TT, n1)       # fixed at 1
    }

    ## ---- discount factors ----
    gamma_vec <- rep(0.9, n1)            # start near 1 (smooth)
    delta     <- 0.9

    ## ---- residual risk: initialise near 1 ----
    lambda_tilde <- matrix(1, TT, n1)

    ## ---- split proportions: initialise from empirical fractions ----
    omega <- array(NA_real_, dim = c(TT, n1, K))
    for (j in seq_len(n1)) {
        for (t in seq_len(TT)) {
            y_kids <- dat$y_fine[t, j, ]
            total  <- sum(y_kids)
            if (total > 0) {
                omega[t, j, ] <- (y_kids + 0.1) / (total + 0.1 * K)
            } else {
                omega[t, j, ] <- rep(1 / K, K)
            }
        }
    }

    ## ---- effective offset ----
    xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2, beta0, beta, phi)

    list(
        beta0 = beta0, beta = beta, phi = phi, u = u,
        tau_phi = tau_phi, r = r, kappa = kappa,
        gamma = gamma_vec, delta = delta,
        lambda_tilde = lambda_tilde, omega = omega,
        xi = xi
    )
}


#' Run one iteration of Algorithm 3
#'
#' @param state     Current parameter state (named list)
#' @param dat       Data list (y_coarse, y_fine, e, x1, x2)
#' @param settings  MCMC_SETTINGS
#' @param priors    MCMC_PRIORS
#' @param spatial   List with H, B_ICAR, R_BHB
#' @param constants List with A0, B0, C0
#' @return          Updated state list, plus diagnostics
run_one_iteration <- function(state, dat, settings, priors, spatial, constants) {

    diag <- list()  # diagnostics for this iteration

    ## ====================================================================
    ## Step 2: Update κ (conjugate, Eq. 30)
    ## ====================================================================
    if (settings$include_nb) {
        state$kappa <- update_kappa(dat$y_coarse, state$lambda_tilde,
                                    state$xi, state$r)
    }
    ## (If !include_nb, kappa stays at 1)

    ## ====================================================================
    ## Step 3: Update (β₀, β) via ESS (Eq. 31)
    ## ====================================================================
    if (settings$include_covariates) {
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            kappa = state$kappa, lambda_tilde = state$lambda_tilde,
            phi = state$phi, priors = priors
        )
        state$beta0 <- beta_result$sample[1]
        state$beta  <- beta_result$sample[2:3]
        diag$beta_ess_reject <- beta_result$n_reject
    }

    ## Recompute ξ after updating β (and before updating φ)
    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## ====================================================================
    ## Step 4: Update φ via ESS (Eq. 32)
    ## ====================================================================
    if (settings$include_icar) {
        phi_result <- update_phi(
            u_current = state$u, B = spatial$B_ICAR,
            R_BHB = spatial$R_BHB, tau_phi = state$tau_phi,
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            beta0 = state$beta0, beta = state$beta,
            kappa = state$kappa, lambda_tilde = state$lambda_tilde
        )
        state$u   <- phi_result$u
        state$phi <- phi_result$phi
        diag$phi_ess_reject <- phi_result$n_reject
    }

    ## Recompute ξ after updating φ
    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## ====================================================================
    ## Step 5: Update τ_φ (conjugate, Eq. 34)
    ## ====================================================================
    if (settings$include_icar) {
        state$tau_phi <- update_tau_phi(state$phi, spatial$H,
                                       dat$n1, priors)
    }

    ## ====================================================================
    ## Step 6: Update r_{1j} via MH (Eq. 35)
    ## ====================================================================
    if (settings$include_nb) {
        r_result <- update_r(state$r, state$kappa, priors,
                             mh_sd = settings$mh_sd_log_r)
        state$r <- r_result$r
        diag$r_accept <- r_result$accept
    }

    ## ====================================================================
    ## Step 7: Update γ_j via MH (Algorithm 1, Eq. 43)
    ## ====================================================================
    gamma_result <- update_gamma(
        state$gamma, dat$y_coarse, state$xi, state$kappa,
        a0 = constants$A0, b0 = constants$B0,
        priors = priors, mh_sd = settings$mh_sd_gamma
    )
    state$gamma <- gamma_result$gamma
    diag$gamma_accept <- gamma_result$accept

    ## ====================================================================
    ## Step 8: FFBS for λ̃ (Algorithm 2)
    ## ====================================================================
    state$lambda_tilde <- ffbs_lambda_all(
        state$gamma, dat$y_coarse, state$xi, state$kappa,
        a0 = constants$A0, b0 = constants$B0
    )

    ## ====================================================================
    ## Step 9: Re-center (β₀, φ, λ̃) onto identified space C (Eqs. 24–27)
    ## ====================================================================
    rc <- recenter(state$beta0, state$phi, state$lambda_tilde)
    state$beta0        <- rc$beta0
    state$phi          <- rc$phi
    state$lambda_tilde <- rc$lambda_tilde

    ## Update u to match the re-centered φ (u = B' φ since B is orthonormal)
    if (settings$include_icar) {
        state$u <- as.numeric(t(spatial$B_ICAR) %*% state$phi)
    }

    ## Recompute ξ after re-centering (β₀ and φ changed)
    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    ## ====================================================================
    ## Step 10: Update δ via MH (Eq. 50)
    ## ====================================================================
    delta_result <- update_delta(
        state$delta, dat$y_fine, c0 = constants$C0,
        priors = priors, mh_sd = settings$mh_sd_delta
    )
    state$delta <- delta_result$delta
    diag$delta_accept <- delta_result$accept

    ## ====================================================================
    ## Step 11: Smooth ω via Dirichlet FFBS (Eqs. 45–48)
    ## ====================================================================
    state$omega <- smooth_omega_all(state$delta, dat$y_fine,
                                    c0 = constants$C0)

    ## ====================================================================
    ## Compute log-likelihood for monitoring
    ## ====================================================================
    diag$loglik <- compute_loglik(dat$y_coarse, state$xi,
                                 state$lambda_tilde, state$kappa)

    return(list(state = state, diag = diag))
}


#' Run the full MCMC sampler
#'
#' @param dat       Simulated dataset
#' @param settings  MCMC_SETTINGS
#' @param priors    MCMC_PRIORS
#' @param spatial   List with H, B_ICAR, R_BHB
#' @param constants List with A0, B0, C0
#' @param verbose   Print progress every this many iterations (0 = silent)
#' @return          List with: samples (thinned post-burn-in), diagnostics, timing
run_mcmc <- function(dat, settings, priors, spatial, constants,
                     verbose = 1000L) {

    n_iter   <- settings$n_iter
    n_burnin <- settings$n_burnin
    n_thin   <- settings$n_thin
    n_stored <- as.integer((n_iter - n_burnin) / n_thin)

    TT <- dat$TT
    n1 <- dat$n1
    K  <- dat$n_children

    ## ---- allocate storage for thinned post-burn-in samples ----
    samples <- list(
        beta0        = numeric(n_stored),
        beta         = matrix(NA_real_, n_stored, length(dat$beta_star)),
        phi          = matrix(NA_real_, n_stored, n1),
        tau_phi      = numeric(n_stored),
        r            = matrix(NA_real_, n_stored, n1),
        gamma        = matrix(NA_real_, n_stored, n1),
        delta        = numeric(n_stored),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT, n1)),
        kappa        = array(NA_real_, dim = c(n_stored, TT, n1)),
        omega        = array(NA_real_, dim = c(n_stored, TT, n1, K)),
        loglik       = numeric(n_stored)
    )

    ## ---- diagnostics: track acceptance rates and ESS rejections ----
    diag_all <- list(
        loglik_trace     = numeric(n_iter),  # every iteration (for trace plot)
        beta_ess_reject  = numeric(n_iter),
        phi_ess_reject   = numeric(n_iter),
        r_accept_rate    = numeric(n1),      # running acceptance count
        gamma_accept_rate = numeric(n1),
        delta_accept_count = 0L
    )

    ## ---- initialise ----
    state <- initialise_state(dat, settings, priors, spatial, constants)

    ## ---- MCMC loop ----
    start_time <- proc.time()
    store_idx  <- 0L

    for (iter in seq_len(n_iter)) {

        ## Run one iteration of Algorithm 3
        result <- run_one_iteration(state, dat, settings, priors,
                                    spatial, constants)
        state <- result$state
        diag  <- result$diag

        ## ---- accumulate diagnostics ----
        diag_all$loglik_trace[iter] <- diag$loglik

        if (!is.null(diag$beta_ess_reject)) {
            diag_all$beta_ess_reject[iter] <- diag$beta_ess_reject
        }
        if (!is.null(diag$phi_ess_reject)) {
            diag_all$phi_ess_reject[iter] <- diag$phi_ess_reject
        }
        if (!is.null(diag$r_accept)) {
            diag_all$r_accept_rate <- diag_all$r_accept_rate + diag$r_accept
        }
        if (!is.null(diag$gamma_accept)) {
            diag_all$gamma_accept_rate <- diag_all$gamma_accept_rate +
                diag$gamma_accept
        }
        if (!is.null(diag$delta_accept) && diag$delta_accept) {
            diag_all$delta_accept_count <- diag_all$delta_accept_count + 1L
        }

        ## ---- store thinned post-burn-in samples ----
        if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
            store_idx <- store_idx + 1L

            samples$beta0[store_idx]        <- state$beta0
            samples$beta[store_idx, ]       <- state$beta
            samples$phi[store_idx, ]        <- state$phi
            samples$tau_phi[store_idx]      <- state$tau_phi
            samples$r[store_idx, ]          <- state$r
            samples$gamma[store_idx, ]      <- state$gamma
            samples$delta[store_idx]        <- state$delta
            samples$lambda_tilde[store_idx, , ] <- state$lambda_tilde
            samples$kappa[store_idx, , ]    <- state$kappa
            samples$omega[store_idx, , , ]  <- state$omega
            samples$loglik[store_idx]       <- diag$loglik
        }

        ## ---- progress ----
        if (verbose > 0 && iter %% verbose == 0) {
            elapsed <- (proc.time() - start_time)[3]
            cat(sprintf("  iter %5d/%d  [%.0fs]  loglik=%.1f  beta0=%.3f  gamma[1]=%.3f\n",
                        iter, n_iter, elapsed, diag$loglik,
                        state$beta0, state$gamma[1]))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]

    ## ---- finalise diagnostics ----
    diag_all$r_accept_rate     <- diag_all$r_accept_rate / n_iter
    diag_all$gamma_accept_rate <- diag_all$gamma_accept_rate / n_iter
    diag_all$delta_accept_rate <- diag_all$delta_accept_count / n_iter
    diag_all$elapsed_sec       <- elapsed_total

    if (verbose > 0) {
        cat(sprintf("\nDone. %d iterations in %.1f sec (%.1f iter/sec)\n",
                    n_iter, elapsed_total, n_iter / elapsed_total))
        cat(sprintf("Stored %d post-burn-in samples (thin=%d)\n",
                    n_stored, n_thin))
        cat(sprintf("Acceptance rates: gamma=%.2f  delta=%.2f  r=%.2f\n",
                    mean(diag_all$gamma_accept_rate),
                    diag_all$delta_accept_rate,
                    mean(diag_all$r_accept_rate)))
    }

    list(
        samples     = samples,
        diagnostics = diag_all,
        settings    = settings,
        priors      = priors,
        method      = method_label(settings),
        n_stored    = n_stored
    )
}
