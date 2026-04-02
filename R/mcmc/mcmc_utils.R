## ==========================================================================
## mcmc_utils.R
## Utility functions for the MCMC sampler
## ==========================================================================


#' Elliptical Slice Sampling (Murray, Adams, MacKay, 2010)
#'
#' Samples from a target proportional to N(prior_mean, Sigma) × L(f).
#' The Gaussian N(prior_mean, Sigma) defines the ellipse; L(f) is evaluated
#' by log_lik_fn.  The prior is handled implicitly by the ellipse, so
#' log_lik_fn should NOT include log-prior terms for the Gaussian part.
#'
#' For preconditioned ESS, log_lik_fn should include a correction term
#' (true_prior / ess_base) — see update_regression.R for details.
#'
#' @param current      Current parameter vector
#' @param prior_sample A draw from N(0, Sigma) (centered, no mean)
#' @param log_lik_fn   Function(param_vector) → scalar log-likelihood
#' @param prior_mean   Mean of the Gaussian base distribution
#' @return             List with: sample, log_lik, n_reject
ess_step <- function(current, prior_sample, log_lik_fn,
                     prior_mean = rep(0, length(current))) {

    ## Center the current state relative to prior_mean
    f <- current - prior_mean

    ## Step 1: current log-likelihood
    cur_log_lik <- log_lik_fn(current)

    ## Step 2: threshold
    log_y <- cur_log_lik + log(runif(1))

    ## Step 3: initial angle and bracket
    theta <- runif(1, min = 0, max = 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta

    n_reject <- 0L

    ## Step 4: shrinking bracket
    repeat {
        f_proposal <- f * cos(theta) + prior_sample * sin(theta)
        proposal   <- f_proposal + prior_mean

        prop_log_lik <- log_lik_fn(proposal)

        if (is.finite(prop_log_lik) && prop_log_lik > log_y) {
            return(list(sample   = proposal,
                        log_lik  = prop_log_lik,
                        n_reject = n_reject))
        }

        ## Shrink bracket
        n_reject <- n_reject + 1L
        if (theta < 0) {
            theta_min <- theta
        } else {
            theta_max <- theta
        }
        theta <- runif(1, min = theta_min, max = theta_max)

        ## Safety: if bracket is tiny, return current (avoid infinite loop)
        if (theta_max - theta_min < 1e-10) {
            return(list(sample   = current,
                        log_lik  = cur_log_lik,
                        n_reject = n_reject))
        }
    }
}


#' Adapt a MH proposal SD using Robbins-Monro
#'
#' Called periodically during burn-in.
#' Adjusts log(sd): up if accepting too much, down if too little.
#'
#' @param current_sd     Current proposal SD (scalar or vector)
#' @param accept_count   Number of acceptances in recent window
#' @param window_size    Length of the adaptation window
#' @param target_rate    Target acceptance rate
#' @param adapt_factor   Robbins-Monro step size
#' @return               Updated proposal SD (same shape as input)
adapt_sd <- function(current_sd, accept_count, window_size,
                     target_rate = 0.30, adapt_factor = 0.5) {

    current_rate <- accept_count / window_size
    log_sd <- log(current_sd)
    log_sd <- log_sd + adapt_factor * (current_rate - target_rate)

    ## Safety bounds
    log_sd <- pmax(log_sd, log(1e-6))
    log_sd <- pmin(log_sd, log(100))

    return(exp(log_sd))
}


#' Logit and inverse-logit
logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))


#' Log-sum-exp (numerically stable)
log_sum_exp <- function(x) {
    mx <- max(x)
    mx + log(sum(exp(x - mx)))
}
