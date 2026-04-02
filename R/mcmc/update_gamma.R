## ==========================================================================
## update_gamma.R
## Metropolis–Hastings update for coarse-level discount factors γ_j
##
## Paper reference: Algorithm 1, Eqs. 40–43
##
## Key idea: γ_j is updated from its MARGINAL conditional posterior,
## obtained by integrating out the entire path λ̃_{1:T,1j} sequentially.
## This requires running the Gamma filtering recursion under a proposed
## γ and evaluating the product of one-step predictive factors (Eq. 41).
##
## The MH proposal is a random walk on logit(γ_j).
## ==========================================================================


#' Compute the log marginal likelihood for γ_j
#'
#' Runs the filtering recursion (Eq. 40) under a given γ and accumulates
#' the log of the one-step predictive factors (Eq. 41).
#'
#' Each predictive factor is a negative-binomial PMF:
#'   p(y_t | D_{t-1}, γ) with size α = γ a_{t-1} and
#'   "success prob" p = ζ_t / (ζ_t + γ b_{t-1})
#'
#' @param gamma_j   Proposed discount factor for region j
#' @param y_j       Length-T vector of coarsest counts for region j
#' @param zeta_j    Length-T vector of effective offsets ζ = ξ·κ for region j
#' @param a0        Initial Gamma filter shape
#' @param b0        Initial Gamma filter rate
#' @return          Scalar log marginal likelihood
log_marginal_gamma <- function(gamma_j, y_j, zeta_j, a0, b0) {

    TT <- length(y_j)
    log_ml <- 0

    ## Filter state (starts at t=0 values)
    at <- a0
    bt <- b0

    for (t in seq_len(TT)) {

        ## ---- discounted prior parameters ----
        alpha <- gamma_j * at     # γ a_{t-1}
        beta  <- gamma_j * bt     # γ b_{t-1}
        y_t   <- y_j[t]
        z_t   <- zeta_j[t]

        ## ---- log predictive factor (Eq. 41) ----
        ## This is log NB(y_t; size = α, prob = z_t / (z_t + β))
        ##   = lgamma(α + y) - lgamma(α) - lgamma(y+1)
        ##     + α log(β) + y log(z_t) - (α + y) log(β + z_t)
        log_pred <- lgamma(alpha + y_t) - lgamma(alpha) - lgamma(y_t + 1) +
                    alpha * log(beta) + y_t * log(z_t) -
                    (alpha + y_t) * log(beta + z_t)

        log_ml <- log_ml + log_pred

        ## ---- update filter (Eq. 38) ----
        at <- alpha + y_t       # = γ a_{t-1} + y_t
        bt <- beta  + z_t       # = γ b_{t-1} + ζ_t
    }

    return(log_ml)
}


#' Update all γ_j via MH on the logit scale (Algorithm 1)
#'
#' @param gamma_current Length-n1 current discount factors
#' @param y_coarse      T × n1 count matrix
#' @param xi            T × n1 effective offset matrix (Eq. 29)
#' @param kappa         T × n1 current κ matrix
#' @param a0            Initial Gamma filter shape
#' @param b0            Initial Gamma filter rate
#' @param priors        List with gamma_a, gamma_b (Beta prior parameters)
#' @param mh_sd         Proposal sd on logit scale
#' @return              List with: gamma (updated), accept (logical vector)
update_gamma <- function(gamma_current, y_coarse, xi, kappa,
                         a0, b0, priors, mh_sd) {

    n1 <- ncol(y_coarse)
    gamma_new <- gamma_current
    accept    <- logical(n1)

    for (j in seq_len(n1)) {

        ## Precompute ζ_j = ξ_j · κ_j for all t   (Eq. 36)
        zeta_j <- xi[, j] * kappa[, j]

        ## ---- propose on logit scale ----
        logit_current  <- logit(gamma_current[j])
        logit_proposal <- logit_current + rnorm(1, mean = 0,
                                                sd = if (length(mh_sd) > 1) mh_sd[j] else mh_sd)
        gamma_proposal <- expit(logit_proposal)

        ## ---- log marginal likelihoods ----
        log_ml_current  <- log_marginal_gamma(gamma_current[j], y_coarse[, j],
                                              zeta_j, a0, b0)
        log_ml_proposal <- log_marginal_gamma(gamma_proposal, y_coarse[, j],
                                              zeta_j, a0, b0)

        ## ---- log priors: Beta(a_γ, b_γ) ----
        log_prior_current  <- dbeta(gamma_current[j], priors$gamma_a,
                                    priors$gamma_b, log = TRUE)
        log_prior_proposal <- dbeta(gamma_proposal, priors$gamma_a,
                                    priors$gamma_b, log = TRUE)

        ## ---- log Jacobian for logit transform ----
        ## d(expit)/d(logit) = expit(x)(1-expit(x)) = γ(1-γ)
        log_jac_current  <- log(gamma_current[j]) + log(1 - gamma_current[j])
        log_jac_proposal <- log(gamma_proposal) + log(1 - gamma_proposal)

        ## ---- MH acceptance ----
        log_alpha <- (log_ml_proposal + log_prior_proposal + log_jac_proposal) -
                     (log_ml_current  + log_prior_current  + log_jac_current)

        ## NaN protection: if the target is undefined at the proposal
        ## (e.g. gamma ≈ 0 causing lgamma(0) = Inf), reject.
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            gamma_new[j] <- gamma_proposal
            accept[j]    <- TRUE
        }
    }

    return(list(gamma = gamma_new, accept = accept))
}
