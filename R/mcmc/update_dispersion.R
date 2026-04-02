## ==========================================================================
## update_dispersion.R
## Metropolis–Hastings update for the NB dispersion parameter r_{1j}
##
## Paper reference: Eq. 35
##   log π(r | κ) = (a_r - 1) log r - b_r r
##                + Σ_t [ r log r - log Γ(r) + (r-1) log κ_t - r κ_t ]
##
## We use a random-walk MH on log(r) for positivity.
## ==========================================================================


#' Log posterior kernel for r_{1j} given κ_{1:T,1j}
#'
#' @param r        Candidate dispersion value (positive scalar)
#' @param kappa_j  Length-T vector of κ values for region j
#' @param priors   List with r_shape (a_r) and r_rate (b_r)
#' @return         Scalar log posterior kernel
log_posterior_r <- function(r, kappa_j, priors) {

    if (r <= 0) return(-Inf)

    TT <- length(kappa_j)

    ## Prior: Ga(a_r, b_r) on r
    log_prior <- (priors$r_shape - 1) * log(r) - priors$r_rate * r

    ## Likelihood: κ_t | r ~ Ga(r, r) for t = 1,...,T
    ## log p(κ_t | r) = r log r - lgamma(r) + (r-1) log κ_t - r κ_t
    sum_log_kappa <- sum(log(kappa_j))
    sum_kappa     <- sum(kappa_j)

    log_lik <- TT * (r * log(r) - lgamma(r)) +
               (r - 1) * sum_log_kappa -
               r * sum_kappa

    return(log_prior + log_lik)
}


#' Update all r_{1j} via MH on the log scale
#'
#' @param r_current  Length-n1 current dispersion vector
#' @param kappa      T × n1 current κ matrix
#' @param priors     List with r_shape, r_rate
#' @param mh_sd      Proposal sd on the log scale
#' @return           List with: r (updated), accept (logical vector)
update_r <- function(r_current, kappa, priors, mh_sd) {

    n1 <- length(r_current)
    r_new  <- r_current
    accept <- logical(n1)

    for (j in seq_len(n1)) {

        ## ---- propose on log scale ----
        log_r_current  <- log(r_current[j])
        log_r_proposal <- log_r_current + rnorm(1, mean = 0,
                                                sd = if (length(mh_sd) > 1) mh_sd[j] else mh_sd)
        r_proposal     <- exp(log_r_proposal)

        ## ---- log posterior at current and proposal ----
        ## The MH target on log(r) is:  log π(r) + log(r)
        ## The extra log(r) is the Jacobian dr/d(log r) = r
        log_target_current  <- log_posterior_r(r_current[j], kappa[, j], priors) +
                               log_r_current
        log_target_proposal <- log_posterior_r(r_proposal, kappa[, j], priors) +
                               log_r_proposal

        ## ---- accept/reject ----
        log_alpha <- log_target_proposal - log_target_current
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            r_new[j]  <- r_proposal
            accept[j] <- TRUE
        }
    }

    return(list(r = r_new, accept = accept))
}
