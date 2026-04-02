## ==========================================================================
## 01_helpers.R
## Utility functions for the MSSTNB simulation DGP
## ==========================================================================

#' Sample from a Dirichlet distribution
#'
#' @param n     Number of draws
#' @param alpha Concentration vector (length K)
#' @return      n × K matrix; each row sums to 1
rdirichlet_local <- function(n, alpha) {
    K <- length(alpha)
    out <- matrix(rgamma(n * K, shape = rep(alpha, each = n),
                         rate = 1),
                  nrow = n, ncol = K)
    out / rowSums(out)
}

# Use MCMCpack::rdirichlet if available, otherwise our own
if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    rdirichlet <- rdirichlet_local
} else {
    rdirichlet <- MCMCpack::rdirichlet
}

#' Standardize a numeric matrix using all entries jointly
#'
#' @param x Numeric matrix
#' @return  Standardized matrix with global mean 0 and sd 1
standardize_matrix <- function(x) {
    mu <- mean(x)
    sd_x <- sd(as.numeric(x))
    if (!is.finite(sd_x) || sd_x <= 0) {
        stop("Cannot standardize matrix: standard deviation is not positive.")
    }
    (x - mu) / sd_x
}

#' Generate an ICAR realisation on the constrained subspace (1'φ = 0)
#'
#' Uses the reduced representation φ = B %*% u, where B is the (n1 × n1-1)
#' basis for the null-complement of 1, precomputed in 00_setup.R.
#'
#' @param tau_phi ICAR precision scalar
#' @param H       Graph Laplacian (n1 × n1)
#' @param B       Orthonormal basis for {φ : 1'φ = 0}
#' @return        Numeric vector of length n1, summing to zero
generate_icar <- function(tau_phi, H, B) {
    n1 <- nrow(H)
    d  <- ncol(B)   # = n1 - 1

    # Precision in reduced space: Q_u = τ_φ · B' H B
    Q_u <- tau_phi * crossprod(B, H %*% B)

    # Sample u ~ N(0, Q_u^{-1}) via Cholesky
    R <- chol(Q_u)                           # upper triangular
    u <- backsolve(R, rnorm(d))              # Q_u^{-1/2} z

    phi <- as.numeric(B %*% u)

    # Safety: enforce exact sum-to-zero (numerical)
    phi <- phi - mean(phi)
    return(phi)
}


#' Generate shared covariates and exposures for one replication
#'
#' These are shared across scenarios for the same replication seed, so that
#' the only difference between scenarios is the model parameters.
#'
#' @param TT  Number of time points
#' @param n1  Number of coarsest regions
#' @return    A list with components: e (T×n1), x1 (T×n1), x2 (T×n1)
generate_inputs <- function(TT, n1) {

    ## ---- exposures ----
    # Baseline population per region: Uniform(500, 2000)
    e_base <- runif(n1, min = 500, max = 2000)

    # Time-varying with ±2% jitter
    e <- matrix(NA_real_, TT, n1)
    for (t in seq_len(TT)) {
        e[t, ] <- e_base * exp(rnorm(n1, mean = 0, sd = 0.02))
    }

    ## ---- covariate x1: continuous, AR(1), region-specific ----
    x1 <- matrix(NA_real_, TT, n1)
    x1[1, ] <- rnorm(n1, mean = 0, sd = 1)
    for (t in 2:TT) {
        x1[t, ] <- 0.95 * x1[t - 1, ] + rnorm(n1, mean = 0, sd = 0.2)
    }

    # x1_raw <- matrix(NA_real_, TT, n1)
    # x1_raw[1, ] <- rnorm(n1, mean = 0, sd = 1)
    # for (t in 2:TT) {
    #     x1_raw[t, ] <- 0.95 * x1_raw[t - 1, ] + rnorm(n1, mean = 0, sd = 0.2)
    # }

    # ## ---- covariate x2: binary, region-specific, time-constant ----
    # x2_vec <- rbinom(n1, size = 1, prob = 0.4)
    # x2 <- matrix(rep(x2_vec, each = TT), nrow = TT, ncol = n1)

    # ## ---- covariate x2: continuous AR(1), region-specific ----
    # # Lower persistence and larger innovations than x1 so the two covariates
    # # are meaningfully distinct while remaining smooth over time.
    x2 <- matrix(NA_real_, TT, n1)
    x2[1, ] <- rnorm(n1, mean = 0, sd = 1)
    for (t in 2:TT) {
        x2[t, ] <- 0.90 * x2[t - 1, ] + rnorm(n1, mean = 0, sd = 1)
    }

    # may diverge if AR parameters are "wrong". This fails
    # x2 <- matrix(NA_real_, TT, n1)
    # x2[1, ] <- rnorm(n1, mean = 0, sd = 1)
    # for (t in 2:TT) {
    #     x2[t, ] <- 0.50 * x2[t - 1, ] + rnorm(n1, mean = 0, sd = .1)
    # }


    ## ---- standardize covariates globally within replication ----
    # x1 <- standardize_matrix(x1_raw)
    # x2 <- standardize_matrix(x2_raw)

    # list(
    #     e = e,
    #     x1 = x1,
    #     x1_raw = x1_raw,
    #     x1_mean = mean(x1_raw),
    #     x1_sd = sd(as.numeric(x1_raw)),
    #     x2 = x2
    # )

    list(e = e, x1 = x1, x2 = x2)
}


#' Compute the effective offset ξ_{t,j}
#'
#' ξ_{t,j} = e_{t,j} * exp(β0 + x1_{t,j}*β1 + x2_{t,j}*β2 + φ_j)
#'
#' @param e     T × n1 exposure matrix
#' @param x1    T × n1 covariate matrix
#' @param x2    T × n1 covariate matrix
#' @param beta0 Scalar intercept
#' @param beta  Length-2 vector of regression coefficients
#' @param phi   Length-n1 ICAR vector
#' @return      T × n1 matrix of effective offsets
compute_xi <- function(e, x1, x2, beta0, beta, phi) {
    TT <- nrow(e)
    n1 <- ncol(e)
    xi <- matrix(NA_real_, TT, n1)
    for (j in seq_len(n1)) {
        linpred_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        xi[, j] <- e[, j] * exp(linpred_j)
    }
    return(xi)
}
