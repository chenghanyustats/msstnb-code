## ==========================================================================
## update_regression.R
## Update regression coefficients (β₀, β) via Elliptical Slice Sampling
##
## Paper reference: Eq. 31
##   log p(β₀, β | ·) = -(β₀ - m₀)²/(2s₀²) - (β - m_β)'V_β^{-1}(β - m_β)/2
##                     + Σ_{t,j} [ y_{t,1j}(β₀ + x'_{t,1j}β) - m_{t,1j} exp(β₀ + x'β + φ_j) ]
##
## where m_{t,1j} = e_{t,1j} κ_{t,1j} λ̃_{t,1j}
##
## ESS uses the Gaussian prior as the ellipse generator and the Poisson
## log-linear log-likelihood as the slice target.
## ==========================================================================


#' Compute the Poisson log-linear log-likelihood for (β₀, β)
#'
#' @param beta_vec   Length-(1+p) vector: c(β₀, β₁, β₂)
#' @param y_coarse   T × n1 count matrix
#' @param e          T × n1 exposure matrix
#' @param x1         T × n1 covariate matrix
#' @param x2         T × n1 covariate matrix
#' @param kappa      T × n1 current κ matrix
#' @param lambda_tilde T × n1 current λ̃ matrix
#' @param phi        Length-n1 current ICAR vector
#' @return           Scalar log-likelihood (not including the prior)
loglik_beta <- function(beta_vec, y_coarse, e, x1, x2,
                        kappa, lambda_tilde, phi) {

    beta0 <- beta_vec[1]
    beta1 <- beta_vec[2]
    beta2 <- beta_vec[3]

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)
    ll <- 0

    for (j in seq_len(n1)) {
        # Linear predictor for region j across all times
        eta_j <- beta0 + beta1 * x1[, j] + beta2 * x2[, j] + phi[j]

        # m_{t,j} = e * κ * λ̃  (the multiplier in the Poisson rate)
        m_j <- e[, j] * kappa[, j] * lambda_tilde[, j]

        # Poisson log-likelihood contribution: y*η - m*exp(η)
        ll <- ll + sum(y_coarse[, j] * eta_j - m_j * exp(eta_j))
    }

    return(ll)
}


#' Update (β₀, β) using Elliptical Slice Sampling
#'
#' @param beta_current  Length-(1+p) current vector c(β₀, β₁, β₂)
#' @param y_coarse      T × n1 count matrix
#' @param e             T × n1 exposure matrix
#' @param x1            T × n1 covariate matrix
#' @param x2            T × n1 covariate matrix
#' @param kappa         T × n1 current κ
#' @param lambda_tilde  T × n1 current λ̃
#' @param phi           Length-n1 current φ
#' @param priors        List with beta0_mean, beta0_sd, beta_mean, beta_sd
#' @return              List with: sample, log_lik, n_reject
update_beta <- function(beta_current, y_coarse, e, x1, x2,
                        kappa, lambda_tilde, phi, priors) {

    ## ---- prior setup ----
    ## Prior: β₀ ~ N(m₀, s₀²),  β ~ N(m_β, diag(s_β²))
    ## Combined: (β₀, β) ~ N(prior_mean, diag(prior_sd²))
    prior_mean <- c(priors$beta0_mean, priors$beta_mean)
    prior_sd   <- c(priors$beta0_sd,   priors$beta_sd)

    ## Draw a centered prior sample: ν ~ N(0, diag(prior_sd²))
    prior_sample <- rnorm(length(prior_mean), mean = 0, sd = prior_sd)

    ## ---- log-likelihood function (fixed everything except beta) ----
    log_lik_fn <- function(beta_vec) {
        loglik_beta(beta_vec, y_coarse, e, x1, x2,
                    kappa, lambda_tilde, phi)
    }

    ## ---- run ESS ----
    result <- ess_step(
        current      = beta_current,
        prior_sample = prior_sample,
        log_lik_fn   = log_lik_fn,
        prior_mean   = prior_mean
    )

    return(result)
}
