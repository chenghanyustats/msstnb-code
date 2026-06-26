## ==========================================================================
## update_icar.R
## Update ICAR field φ via ESS and precision τ_φ via conjugate Gamma
##
## Paper references:
##   Eq. 32 — log posterior kernel for φ (on constrained subspace 1'φ=0)
##   Eq. 33 — reduced representation φ = Bu, u ~ N(0, (τ_φ B'HB)^{-1})
##   Eq. 34 — conjugate full conditional for τ_φ
##
## The ESS operates on the reduced coordinates u ∈ R^{n1-1}, then
## converts back to φ = Bu.
## ==========================================================================


#' Compute the Poisson log-likelihood as a function of φ
#'
#' @param phi          Length-n1 candidate ICAR vector
#' @param y_coarse     T × n1 count matrix
#' @param e            T × n1 exposure matrix
#' @param x1           T × n1 covariate matrix
#' @param x2           T × n1 covariate matrix
#' @param beta0        Current intercept
#' @param beta         Current regression coefficients (length p)
#' @param kappa        T × n1 current κ
#' @param lambda_tilde T × n1 current λ̃
#' @return             Scalar log-likelihood
loglik_phi <- function(phi, y_coarse, e, x1, x2,
                       beta0, beta, kappa, lambda_tilde) {

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)
    ll <- 0

    for (j in seq_len(n1)) {
        # m_{t,j} = e * κ * λ̃
        m_j <- e[, j] * kappa[, j] * lambda_tilde[, j]

        # Poisson log-likelihood: y*φ_j - m*exp(β₀ + x'β + φ_j)
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        ll <- ll + sum(y_coarse[, j] * eta_j - m_j * exp(eta_j))
    }

    ## Note: the contribution y*(β₀+x'β) doesn't depend on φ, so it
    ## cancels in the ESS likelihood ratio. But including it is harmless
    ## and makes the function self-contained.

    return(ll)
}


#' Update φ (via reduced coordinates u) using Elliptical Slice Sampling
#'
#' @param u_current    Length-(n1-1) current reduced coordinates
#' @param B            n1 × (n1-1) orthonormal basis (1'Bv = 0 for all v)
#' @param R_BHB        Upper-triangular Cholesky of B'HB (precomputed)
#' @param tau_phi      Current ICAR precision
#' @param y_coarse     T × n1 count matrix
#' @param e, x1, x2   Inputs
#' @param beta0, beta  Current regression parameters
#' @param kappa        T × n1 current κ
#' @param lambda_tilde T × n1 current λ̃
#' @return             List with: u (updated), phi (= B %*% u), log_lik, n_reject
update_phi <- function(u_current, B, R_BHB, tau_phi,
                       y_coarse, e, x1, x2,
                       beta0, beta, kappa, lambda_tilde) {

    d <- length(u_current)   # n1 - 1

    ## ---- draw a centered prior sample ----
    ## Prior on u: N(0, (τ_φ B'HB)^{-1})
    ## Sample: u_nu = R_BHB^{-1} z / sqrt(τ_φ)
    z <- rnorm(d)
    u_prior_sample <- backsolve(R_BHB, z) / sqrt(tau_phi)

    ## ---- log-likelihood function (maps u → φ → loglik) ----
    log_lik_fn <- function(u_vec) {
        phi_candidate <- as.numeric(B %*% u_vec)
        loglik_phi(phi_candidate, y_coarse, e, x1, x2,
                   beta0, beta, kappa, lambda_tilde)
    }

    ## ---- run ESS (prior mean = 0) ----
    result <- ess_step(
        current      = u_current,
        prior_sample = u_prior_sample,
        log_lik_fn   = log_lik_fn,
        prior_mean   = rep(0, d)
    )

    ## ---- convert back to φ ----
    u_new   <- result$sample
    phi_new <- as.numeric(B %*% u_new)

    return(list(
        u       = u_new,
        phi     = phi_new,
        log_lik = result$log_lik,
        n_reject = result$n_reject
    ))
}


#' Update ICAR precision τ_φ from its conjugate full conditional
#'
#' Paper reference: Eq. 34
#'   τ_φ | φ ~ Ga(a_φ + (n1-1)/2,  b_φ + φ'Hφ/2)
#'
#' @param phi        Length-n1 current ICAR vector
#' @param H          n1 × n1 graph Laplacian
#' @param n1         Number of coarsest regions
#' @param priors     List with tau_phi_shape, tau_phi_rate
#' @return           Scalar: updated τ_φ
update_tau_phi <- function(phi, H, n1, priors) {

    shape_post <- priors$tau_phi_shape + (n1 - 1) / 2
    rate_post  <- priors$tau_phi_rate + as.numeric(t(phi) %*% H %*% phi) / 2

    tau_phi_new <- rgamma(1, shape = shape_post, rate = rate_post)

    return(tau_phi_new)
}
