## ==========================================================================
## update_delta.R
## Metropolis–Hastings update for split discount factors δ_{lj}
##
## Paper reference: Eqs. 49–50
##
## Same logic as the γ update: δ is updated from its marginal conditional
## posterior, obtained by integrating out ω_{1:T,lj} via the Dirichlet–
## multinomial predictive factors.
##
## In our simulation: L=2, each coarsest region has K=3 children, and
## δ is a scalar shared across all internal nodes.
## ==========================================================================


#' Compute the log marginal likelihood for δ (all internal nodes combined)
#'
#' Runs the Dirichlet filtering recursion under a given δ and accumulates
#' the log of the Dirichlet–multinomial predictive factors (Eq. 49).
#'
#' @param delta     Proposed split discount factor
#' @param y_fine    T × n1 × K array of child counts
#' @param c0        Length-K initial Dirichlet concentration
#' @return          Scalar log marginal likelihood (summed across all nodes)
log_marginal_delta <- function(delta, y_fine, c0) {

    TT <- dim(y_fine)[1]
    n1 <- dim(y_fine)[2]
    K  <- dim(y_fine)[3]

    log_ml <- 0

    for (j in seq_len(n1)) {

        ## Filter state (starts at t=0 values)
        ct <- c0   # length-K vector

        for (t in seq_len(TT)) {

            ## ---- discounted prior concentration ----
            c_prior   <- delta * ct          # δ c_{t-1}
            c_prior_p <- sum(c_prior)        # δ c_{t-1,+}
            y_parent  <- sum(y_fine[t, j, ]) # y_{t,lj}
            y_kids    <- y_fine[t, j, ]      # length-K vector

            ## ---- log Dirichlet–multinomial predictive (Eq. 49) ----
            ## log p = lgamma(δc_+) + lgamma(y_parent + 1) - lgamma(δc_+ + y_parent)
            ##       + Σ_k [ lgamma(δc_k + y_k) - lgamma(δc_k) - lgamma(y_k + 1) ]
            log_pred <- lgamma(c_prior_p) + lgamma(y_parent + 1) -
                        lgamma(c_prior_p + y_parent)

            for (k in seq_len(K)) {
                log_pred <- log_pred +
                    lgamma(c_prior[k] + y_kids[k]) -
                    lgamma(c_prior[k]) -
                    lgamma(y_kids[k] + 1)
            }

            log_ml <- log_ml + log_pred

            ## ---- update Dirichlet filter (Eq. 45) ----
            ct <- c_prior + y_kids    # c_t = δ c_{t-1} + y_children
        }
    }

    return(log_ml)
}


#' Update δ via MH on the logit scale
#'
#' @param delta_current  Current scalar discount factor
#' @param y_fine         T × n1 × K array of child counts
#' @param c0             Length-K initial Dirichlet concentration
#' @param priors         List with delta_a, delta_b (Beta prior)
#' @param mh_sd          Proposal sd on logit scale
#' @return               List with: delta (updated), accept (logical)
update_delta <- function(delta_current, y_fine, c0, priors, mh_sd) {

    ## ---- propose on logit scale ----
    logit_current  <- logit(delta_current)
    logit_proposal <- logit_current + rnorm(1, mean = 0, sd = mh_sd)
    delta_proposal <- expit(logit_proposal)

    ## ---- log marginal likelihoods ----
    log_ml_current  <- log_marginal_delta(delta_current, y_fine, c0)
    log_ml_proposal <- log_marginal_delta(delta_proposal, y_fine, c0)

    ## ---- log priors: Beta(a_δ, b_δ) ----
    log_prior_current  <- dbeta(delta_current, priors$delta_a,
                                priors$delta_b, log = TRUE)
    log_prior_proposal <- dbeta(delta_proposal, priors$delta_a,
                                priors$delta_b, log = TRUE)

    ## ---- log Jacobian for logit transform ----
    log_jac_current  <- log(delta_current) + log(1 - delta_current)
    log_jac_proposal <- log(delta_proposal) + log(1 - delta_proposal)

    ## ---- MH acceptance ----
    log_alpha <- (log_ml_proposal + log_prior_proposal + log_jac_proposal) -
                 (log_ml_current  + log_prior_current  + log_jac_current)

    accept <- FALSE
    delta_new <- delta_current

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        delta_new <- delta_proposal
        accept    <- TRUE
    }

    return(list(delta = delta_new, accept = accept))
}
