## ==========================================================================
## mcmc_utils.R
## Utility functions for the MCMC sampler
## ==========================================================================


#' Elliptical Slice Sampling (Murray, Adams, MacKay, 2010)
#'
#' Samples from a posterior proportional to N(prior_mean, Sigma) × L(f),
#' where L(f) is an arbitrary likelihood function.
#'
#' The algorithm traces an ellipse defined by the current state and a
#' prior sample, and finds a point on the ellipse with higher likelihood.
#'
#' @param current      Current parameter vector
#' @param prior_sample A draw from N(0, Sigma) — the CENTERED prior
#' @param log_lik_fn   Function(param_vector) → scalar log-likelihood
#' @param prior_mean   Prior mean vector (default: zero vector)
#' @return List with: sample, log_lik, n_reject
ess_step <- function(current, prior_sample, log_lik_fn,
                     prior_mean = rep(0, length(current))) {

    # Center the current state
    f <- current - prior_mean

    # Step 1: Evaluate current log-likelihood
    cur_log_lik <- log_lik_fn(current)

    # Step 2: Set the log-likelihood threshold
    log_y <- cur_log_lik + log(runif(1))

    # Step 3: Draw initial angle and bracket
    theta <- runif(1, min = 0, max = 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta

    n_reject <- 0L

    # Step 4: Shrinking bracket loop
    repeat {
        # Propose on the ellipse
        f_proposal <- f * cos(theta) + prior_sample * sin(theta)
        proposal   <- f_proposal + prior_mean

        # Evaluate proposal
        prop_log_lik <- log_lik_fn(proposal)

        if (prop_log_lik > log_y) {
            # Accept
            return(list(
                sample   = proposal,
                log_lik  = prop_log_lik,
                n_reject = n_reject
            ))
        }

        # Shrink the bracket
        n_reject <- n_reject + 1L
        if (theta < 0) {
            theta_min <- theta
        } else {
            theta_max <- theta
        }
        theta <- runif(1, min = theta_min, max = theta_max)
    }
}


#' Draw a sample from N(0, Sigma) where Sigma = (tau * B'HB)^{-1}
#'
#' Uses the precomputed Cholesky R_base of B'HB:
#'   Q = tau * B'HB = tau * R_base' R_base
#'   Sigma = Q^{-1}
#'   sample = R_base^{-1} z / sqrt(tau),  z ~ N(0, I)
#'
#' @param R_base  Upper-triangular Cholesky of B'HB (precomputed)
#' @param tau     ICAR precision scalar
#' @return        Vector of length ncol(R_base)
draw_icar_prior_sample <- function(R_base, tau) {
    d <- ncol(R_base)
    z <- rnorm(d)
    u <- backsolve(R_base, z) / sqrt(tau)
    return(u)
}


#' Logit and inverse-logit transforms
logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))


#' Log-sum-exp (numerically stable)
log_sum_exp <- function(x) {
    mx <- max(x)
    mx + log(sum(exp(x - mx)))
}
