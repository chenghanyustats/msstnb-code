## ==========================================================================
## 02_simulate_data.R  -- Step 3 revision + numerical stability fix
## Core data-generating process for the MSSTNB simulation
##
## REVISION GOALS:
## 1. Use the matched Beta-innovation state law for lambda_tilde.
## 2. Store an identified-space version of the truth, consistent with the
##    deterministic re-centering used by the fitting algorithm.
## 3. Prevent machine-boundary zeros from rbeta()/underflow from breaking
##    identified-space re-centering.
## ==========================================================================

#' Re-center truth to the identified space used by the fitted model
#'
#' For each coarsest region j, define
#'   s_j = T^{-1} sum_{t=1}^T log lambda_tilde[t, j].
#' Then set
#'   lambda_ident[t, j]  = lambda[t, j] * exp(-s_j),
#'   lambda0_ident[j]    = lambda0[j] * exp(-s_j),
#'   phi_ident[j]        = phi[j] + s_j - sbar,
#'   beta0_ident         = beta0 + sbar,
#' where sbar is the spatial average of the region-specific shifts.
#'
#' This leaves the linear predictor unchanged while enforcing the temporal-mean
#' constraint on log lambda_tilde over t = 1, ..., T.
#'
#' @param beta0         Scalar intercept
#' @param phi           Length-n1 vector (baseline phi)
#' @param lambda_tilde  T x n1 matrix for observed times
#' @param lambda_tilde0 Length-n1 vector for t = 0 state
#' @param phi_after     Optional length-n1 vector for changepoint scenarios
#' @param eps_pos       Tiny positive floor used only for machine-boundary zeros
#' @return Named list with identified-space quantities and shifts
recenter_truth_to_identified <- function(beta0,
                                         phi,
                                         lambda_tilde,
                                         lambda_tilde0,
                                         phi_after = NULL,
                                         eps_pos = 1e-300) {
    if (any(lambda_tilde < 0) || any(lambda_tilde0 < 0)) {
        stop("Negative lambda_tilde values encountered before re-centering.")
    }

    if (any(lambda_tilde == 0) || any(lambda_tilde0 == 0)) {
        warning(
            paste(
                "Zero lambda_tilde values encountered before re-centering;",
                "replacing only machine-boundary zeros by a tiny positive epsilon."
            )
        )
        lambda_tilde  <- pmax(lambda_tilde,  eps_pos)
        lambda_tilde0 <- pmax(lambda_tilde0, eps_pos)
    }

    shift_region <- colMeans(log(lambda_tilde))
    shift_global <- mean(shift_region)

    lambda_tilde_ident  <- sweep(lambda_tilde,  2, exp(-shift_region), `*`)
    lambda_tilde0_ident <- lambda_tilde0 * exp(-shift_region)

    phi_star_ident   <- phi + shift_region - shift_global
    beta0_star_ident <- beta0 + shift_global

    phi_after_ident <- NULL
    if (!is.null(phi_after)) {
        phi_after_ident <- phi_after + shift_region - shift_global
    }

    list(
        beta0_star_ident = beta0_star_ident,
        phi_star_ident = phi_star_ident,
        phi_after_ident = phi_after_ident,
        lambda_tilde_ident = lambda_tilde_ident,
        lambda_tilde0_ident = lambda_tilde0_ident,
        recenter_shift_region = shift_region,
        recenter_shift_global = shift_global
    )
}


#' Numerically stable Beta draw on the open interval (0, 1)
#'
#' R's rbeta() can occasionally return exact 0 or 1 when shape parameters are
#' below 1 and the density is sharply concentrated near the boundaries. The
#' theoretical model has support on the open interval (0, 1), so here we apply
#' a tiny truncation only to machine-boundary values.
#'
#' @param shape1 First Beta shape parameter
#' @param shape2 Second Beta shape parameter
#' @param eps_lo Lower truncation level for boundary 0
#' @param eps_hi Upper truncation gap from 1 for boundary 1
#' @return       A scalar in (0, 1)
rbeta_open01 <- function(shape1, shape2,
                         eps_lo = 1e-300,
                         eps_hi = .Machine$double.eps) {
    x <- rbeta(1, shape1 = shape1, shape2 = shape2)
    if (!is.finite(x)) stop("Non-finite value returned by rbeta().")
    x <- max(x, eps_lo)
    x <- min(x, 1 - eps_hi)
    x
}


#' Simulate one complete dataset from the MSSTNB DGP
#'
#' Generates coarsest counts via Poisson–Gamma (NB) with Beta-innovation
#' discounted Gamma evolution for lambda_tilde, then allocates to children via
#' Dirichlet-discounted multinomial splitting.
#'
#' The coarse-level DGP follows the matched state construction:
#'   lambda_0 ~ Ga(a0, b0)
#'   eta_t | D_{t-1} ~ Beta(gamma * a_{t-1}, (1-gamma) * a_{t-1})
#'   lambda_t = lambda_{t-1} * eta_t / gamma
#'   y_t | lambda_t, kappa_t ~ Pois(xi_t * lambda_t * kappa_t)
#'   a_t = gamma * a_{t-1} + y_t
#'   b_t = gamma * b_{t-1} + zeta_t,   zeta_t = xi_t * kappa_t
#'
#' In addition, this function stores identified-space truth that matches the
#' sampler's deterministic re-centering convention.
#'
#' @param inputs    List from generate_inputs(): e, x1, x2 (and optionally raw versions)
#' @param params    List of true parameters
#' @param constants List with TT, N1, N_CHILDREN, A0, B0, C0
#' @return          A list with observed data, raw truth, and identified truth
simulate_one_dataset <- function(inputs, params, constants) {

    ## ---- unpack ----
    e  <- inputs$e
    x1 <- inputs$x1
    x2 <- inputs$x2

    TT         <- constants$TT
    n1         <- constants$N1
    n_children <- constants$N_CHILDREN
    a0         <- constants$A0
    b0         <- constants$B0
    c0         <- constants$C0

    beta0_star   <- params$beta0
    beta_star    <- params$beta
    phi_star     <- params$phi
    r_star       <- params$r
    gamma_star   <- params$gamma
    delta_star   <- params$delta
    tau_phi_star <- params$tau_phi

    is_changepoint <- !is.null(params$phi_after)
    phi_after      <- params$phi_after
    change_t       <- params$change_t

    ## ---- basic validation ----
    if (any(gamma_star <= 0 | gamma_star >= 1)) {
        stop("All gamma values must lie strictly in (0, 1) for Beta innovation DGP.")
    }
    if (a0 <= 0 || b0 <= 0) {
        stop("A0 and B0 must be strictly positive.")
    }
    if (any(r_star <= 0)) {
        stop("All r values must be strictly positive.")
    }
    if (length(c0) != n_children || any(c0 <= 0)) {
        stop("C0 must be a positive vector of length N_CHILDREN.")
    }

    ## ---- allocate storage ----
    y_coarse      <- matrix(NA_integer_, TT, n1)
    y_fine        <- array(NA_integer_, dim = c(TT, n1, n_children))
    lambda_tilde  <- matrix(NA_real_, TT, n1)
    lambda_tilde0 <- numeric(n1)
    eta_lambda    <- matrix(NA_real_, TT, n1)
    kappa         <- matrix(NA_real_, TT, n1)
    omega         <- array(NA_real_, dim = c(TT, n1, n_children))
    xi            <- matrix(NA_real_, TT, n1)
    mu            <- matrix(NA_real_, TT, n1)
    a_filt        <- matrix(NA_real_, TT, n1)
    b_filt        <- matrix(NA_real_, TT, n1)

    ## ---- initial state at t = 0 ----
    lambda_prev <- rgamma(n1, shape = a0, rate = b0)
    lambda_tilde0[] <- lambda_prev

    ## ---- filter states ----
    at <- rep(a0, n1)
    bt <- rep(b0, n1)
    ct <- matrix(rep(c0, each = n1), nrow = n1, ncol = n_children)

    ## ---- numerical-stability diagnostics ----
    n_eta_boundary <- 0L
    n_lambda_floor <- 0L

    ## ---- main DGP loop ----
    for (t in seq_len(TT)) {

        if (is_changepoint && t > change_t) {
            phi_t <- phi_after
        } else {
            phi_t <- phi_star
        }

        for (j in seq_len(n1)) {

            ## Step 1: matched Beta-innovation evolution for lambda_tilde
            alpha_eta <- gamma_star[j] * at[j]
            beta_eta  <- (1 - gamma_star[j]) * at[j]

            eta_raw <- rbeta(1, shape1 = alpha_eta, shape2 = beta_eta)
            eta_stable <- max(min(eta_raw, 1 - .Machine$double.eps), 1e-300)
            if (eta_stable != eta_raw) n_eta_boundary <- n_eta_boundary + 1L
            eta_lambda[t, j] <- eta_stable

            lambda_raw <- lambda_prev[j] * eta_lambda[t, j] / gamma_star[j]
            if (!is.finite(lambda_raw) || lambda_raw <= 0) {
                n_lambda_floor <- n_lambda_floor + 1L
                lambda_raw <- 1e-300
            }
            lambda_tilde[t, j] <- lambda_raw

            ## Step 2: draw NB random effect kappa
            kappa[t, j] <- rgamma(1, shape = r_star[j], rate = r_star[j])

            ## Step 3: compute Poisson rate and draw coarsest count
            xi[t, j] <- e[t, j] * exp(beta0_star +
                                       beta_star[1] * x1[t, j] +
                                       beta_star[2] * x2[t, j] +
                                       phi_t[j])
            mu[t, j] <- xi[t, j] * lambda_tilde[t, j] * kappa[t, j]
            y_coarse[t, j] <- rpois(1, lambda = mu[t, j])

            ## Step 4: update Gamma filter state
            zeta_tj <- xi[t, j] * kappa[t, j]
            at[j] <- gamma_star[j] * at[j] + y_coarse[t, j]
            bt[j] <- gamma_star[j] * bt[j] + zeta_tj

            a_filt[t, j] <- at[j]
            b_filt[t, j] <- bt[j]
            lambda_prev[j] <- lambda_tilde[t, j]
        }

        ## Step 5: split proportions and multinomial allocation
        for (j in seq_len(n1)) {
            c_prior <- delta_star * ct[j, ]
            omega[t, j, ] <- rdirichlet(1, c_prior)

            if (y_coarse[t, j] > 0L) {
                y_fine[t, j, ] <- rmultinom(1, size = y_coarse[t, j],
                                            prob = omega[t, j, ])
            } else {
                y_fine[t, j, ] <- rep(0L, n_children)
            }

            ct[j, ] <- c_prior + y_fine[t, j, ]
        }
    }

    ## ---- verify tree coherence ----
    coherent <- all(y_coarse == apply(y_fine, c(1, 2), sum))
    if (!coherent) {
        warning("Tree coherence violated! Check DGP logic.")
    }

    ## ---- identified-space truth storage ----
    ident <- recenter_truth_to_identified(
        beta0 = beta0_star,
        phi = phi_star,
        lambda_tilde = lambda_tilde,
        lambda_tilde0 = lambda_tilde0,
        phi_after = phi_after
    )

    ## ---- package output ----
    out <- list(
        ## ----- observed data -----
        y_coarse = y_coarse,
        y_fine   = y_fine,

        ## ----- known inputs -----
        e  = e,
        x1 = x1,
        x2 = x2,

        ## ----- raw true latent states -----
        lambda_tilde  = lambda_tilde,
        lambda_tilde0 = lambda_tilde0,
        eta_lambda    = eta_lambda,
        kappa         = kappa,
        omega         = omega,
        psi           = lambda_tilde * kappa,
        xi            = xi,
        mu            = mu,

        ## ----- identified-space true latent states -----
        lambda_tilde_ident  = ident$lambda_tilde_ident,
        lambda_tilde0_ident = ident$lambda_tilde0_ident,

        ## ----- filter diagnostics -----
        a_filt = a_filt,
        b_filt = b_filt,

        ## ----- raw true parameters -----
        beta0_star   = beta0_star,
        beta_star    = beta_star,
        phi_star     = phi_star,
        phi_after    = phi_after,
        change_t     = change_t,
        tau_phi_star = tau_phi_star,
        r_star       = r_star,
        gamma_star   = gamma_star,
        delta_star   = delta_star,

        ## ----- identified-space true parameters -----
        beta0_star_ident   = ident$beta0_star_ident,
        phi_star_ident     = ident$phi_star_ident,
        phi_after_ident    = ident$phi_after_ident,
        recenter_shift_region = ident$recenter_shift_region,
        recenter_shift_global = ident$recenter_shift_global,

        ## ----- structural info -----
        TT         = TT,
        n1         = n1,
        n_children = n_children,

        ## ----- diagnostics -----
        coherent = coherent,
        n_eta_boundary = n_eta_boundary,
        n_lambda_floor = n_lambda_floor
    )

    if (!is.null(inputs$x1_raw)) out$x1_raw <- inputs$x1_raw
    if (!is.null(inputs$x2_raw)) out$x2_raw <- inputs$x2_raw

    return(out)
}
