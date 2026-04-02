## ==========================================================================
## update_gamma_common.R
## Common-gamma Metropolis–Hastings update for a single coarse-level
## discount factor shared across all coarsest regions.
##
## This is a minimal modification of update_gamma.R:
##   - gamma is scalar instead of length-n1 vector
##   - the marginal log-likelihood sums the regionwise contributions
##     under the shared candidate gamma
## ==========================================================================

#' Compute the log marginal likelihood for a common gamma
#'
#' @param gamma    Proposed common discount factor
#' @param y_coarse T x n1 count matrix
#' @param xi       T x n1 effective offset matrix
#' @param kappa    T x n1 current kappa matrix
#' @param a0       Initial Gamma filter shape
#' @param b0       Initial Gamma filter rate
#' @return         Scalar log marginal likelihood summed over regions
log_marginal_gamma_common <- function(gamma, y_coarse, xi, kappa, a0, b0) {
    n1 <- ncol(y_coarse)
    log_ml <- 0

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]
        log_ml <- log_ml + log_marginal_gamma(gamma, y_coarse[, j], zeta_j, a0, b0)
    }

    log_ml
}


#' Update a common gamma via MH on the logit scale
#'
#' @param gamma_current Scalar current common gamma
#' @param y_coarse      T x n1 count matrix
#' @param xi            T x n1 effective offset matrix
#' @param kappa         T x n1 current kappa matrix
#' @param a0            Initial Gamma filter shape
#' @param b0            Initial Gamma filter rate
#' @param priors        List with gamma_a, gamma_b
#' @param mh_sd         Scalar proposal sd on the logit scale
#' @return              List with gamma (scalar) and accept (logical)
update_gamma_common <- function(gamma_current, y_coarse, xi, kappa,
                                a0, b0, priors, mh_sd) {
    logit_current  <- logit(gamma_current)
    logit_proposal <- logit_current + rnorm(1, mean = 0, sd = mh_sd)
    gamma_proposal <- expit(logit_proposal)

    log_ml_current  <- log_marginal_gamma_common(gamma_current, y_coarse, xi, kappa, a0, b0)
    log_ml_proposal <- log_marginal_gamma_common(gamma_proposal, y_coarse, xi, kappa, a0, b0)

    log_prior_current  <- dbeta(gamma_current, priors$gamma_a, priors$gamma_b, log = TRUE)
    log_prior_proposal <- dbeta(gamma_proposal, priors$gamma_a, priors$gamma_b, log = TRUE)

    log_jac_current  <- log(gamma_current) + log(1 - gamma_current)
    log_jac_proposal <- log(gamma_proposal) + log(1 - gamma_proposal)

    log_alpha <- (log_ml_proposal + log_prior_proposal + log_jac_proposal) -
                 (log_ml_current  + log_prior_current  + log_jac_current)

    gamma_new <- gamma_current
    accept <- FALSE
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        gamma_new <- gamma_proposal
        accept <- TRUE
    }

    list(gamma = gamma_new, accept = accept)
}
