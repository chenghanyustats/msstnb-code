## ==========================================================================
## update_gamma_revised.R
## Robust Metropolis-Hastings update for coarse-level discount factors gamma_j
##
## Main method:
##   lambda-collapsed update for gamma_j.
##   For each candidate gamma_j, the latent lambda_tilde path is integrated
##   out sequentially using the Gamma-Poisson filtering recursion.
##
## Observation model used in this update:
##   y_tj | lambda_tj, zeta_tj ~ Poisson(zeta_tj * lambda_tj),
##   lambda_tj | D_{t-1}, gamma_j ~ Gamma(gamma_j a_{t-1,j},
##                                        gamma_j b_{t-1,j}),
## where zeta_tj = xi_tj * kappa_tj.
##
## The resulting one-step predictive factor is Negative Binomial.
##
## Notes:
##   This update integrates out lambda_tilde, but still conditions on kappa.
##   This preserves the simple Gamma-Poisson filtering recursion.
## ==========================================================================

## --------------------------------------------------------------------------
## Safe transform helpers
## --------------------------------------------------------------------------

safe_logit <- function(p, eps = 1e-12) {
    p <- pmin(pmax(p, eps), 1 - eps)
    log(p) - log1p(-p)
}

safe_expit <- function(x) {
    out <- numeric(length(x))
    idx <- x >= 0
    out[idx] <- 1 / (1 + exp(-x[idx]))
    ex <- exp(x[!idx])
    out[!idx] <- ex / (1 + ex)
    out
}

check_gamma_inputs <- function(gamma_current, y_coarse, xi, kappa,
                               a0, b0, priors, mh_sd) {
    if (!is.matrix(y_coarse)) {
        stop("y_coarse must be a matrix.")
    }
    if (!is.matrix(xi)) {
        stop("xi must be a matrix.")
    }
    if (!is.matrix(kappa)) {
        stop("kappa must be a matrix.")
    }
    if (!identical(dim(y_coarse), dim(xi))) {
        stop("xi must have the same dimensions as y_coarse.")
    }
    if (!identical(dim(y_coarse), dim(kappa))) {
        stop("kappa must have the same dimensions as y_coarse.")
    }

    n1 <- ncol(y_coarse)
    if (length(gamma_current) != n1) {
        stop("length(gamma_current) must equal ncol(y_coarse).")
    }
    if (any(!is.finite(gamma_current)) || any(gamma_current <= 0) || any(gamma_current >= 1)) {
        stop("gamma_current must contain finite values strictly inside (0, 1).")
    }
    if (!is.finite(a0) || !is.finite(b0) || a0 <= 0 || b0 <= 0) {
        stop("a0 and b0 must be positive finite scalars.")
    }
    if (is.null(priors$gamma_a) || is.null(priors$gamma_b)) {
        stop("priors must contain gamma_a and gamma_b.")
    }
    if (!is.finite(priors$gamma_a) || !is.finite(priors$gamma_b) ||
        priors$gamma_a <= 0 || priors$gamma_b <= 0) {
        stop("priors$gamma_a and priors$gamma_b must be positive finite scalars.")
    }
    if (length(mh_sd) != 1 && length(mh_sd) != n1) {
        stop("mh_sd must be either a scalar or a length-n1 vector.")
    }
    if (any(!is.finite(mh_sd)) || any(mh_sd <= 0)) {
        stop("mh_sd must contain positive finite values.")
    }
    if (any(!is.finite(y_coarse)) || any(y_coarse < 0)) {
        stop("y_coarse must contain nonnegative finite counts.")
    }
    if (any(abs(y_coarse - round(y_coarse)) > 1e-8)) {
        stop("y_coarse must contain integer-valued counts.")
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("xi must contain positive finite values.")
    }
    if (any(!is.finite(kappa)) || any(kappa <= 0)) {
        stop("kappa must contain positive finite values.")
    }

    invisible(TRUE)
}

## --------------------------------------------------------------------------
## Log one-step predictive marginal likelihood for one region
## --------------------------------------------------------------------------

#' Compute the lambda-collapsed log marginal likelihood for one gamma_j.
#'
#' This function computes
#'
#'   log p(y_{1:T,j} | gamma_j, zeta_{1:T,j})
#'
#' by the sequential Gamma-Poisson filter.  It integrates out the full
#' lambda_tilde path, conditional on zeta_j = xi_j * kappa_j.
#'
#' @param gamma_j Proposed discount factor, strictly in (0, 1).
#' @param y_j Length-T nonnegative count vector.
#' @param zeta_j Length-T positive effective offset vector.
#' @param a0 Initial Gamma shape.
#' @param b0 Initial Gamma rate.
#' @param eps Numerical lower bound for positive quantities.
#' @param return_filter If TRUE, also return filtering states and diagnostics.
#' @return Scalar log marginal likelihood, or list if return_filter = TRUE.
log_marginal_gamma <- function(gamma_j, y_j, zeta_j, a0, b0,
                               eps = 1e-12,
                               return_filter = FALSE) {
    if (!is.finite(gamma_j) || gamma_j <= eps || gamma_j >= 1 - eps) {
        if (return_filter) {
            return(list(log_ml = -Inf, a_filt = NULL, b_filt = NULL,
                        ok = FALSE, reason = "gamma outside stable interior"))
        }
        return(-Inf)
    }
    if (length(y_j) != length(zeta_j)) {
        stop("y_j and zeta_j must have the same length.")
    }
    if (!is.finite(a0) || !is.finite(b0) || a0 <= 0 || b0 <= 0) {
        stop("a0 and b0 must be positive finite scalars.")
    }
    if (any(!is.finite(y_j)) || any(y_j < 0)) {
        if (return_filter) {
            return(list(log_ml = -Inf, a_filt = NULL, b_filt = NULL,
                        ok = FALSE, reason = "invalid y_j"))
        }
        return(-Inf)
    }
    if (any(!is.finite(zeta_j)) || any(zeta_j <= 0)) {
        if (return_filter) {
            return(list(log_ml = -Inf, a_filt = NULL, b_filt = NULL,
                        ok = FALSE, reason = "invalid zeta_j"))
        }
        return(-Inf)
    }

    TT <- length(y_j)
    log_ml <- 0

    at <- a0
    bt <- b0

    if (return_filter) {
        a_filt <- numeric(TT)
        b_filt <- numeric(TT)
        log_pred <- numeric(TT)
    }

    for (t in seq_len(TT)) {
        alpha <- gamma_j * at
        beta_rate <- gamma_j * bt
        y_t <- y_j[t]
        z_t <- zeta_j[t]

        if (!is.finite(alpha) || !is.finite(beta_rate) ||
            alpha <= eps || beta_rate <= eps) {
            if (return_filter) {
                return(list(log_ml = -Inf, a_filt = a_filt, b_filt = b_filt,
                            log_pred = log_pred, ok = FALSE,
                            reason = "nonpositive discounted filter parameter"))
            }
            return(-Inf)
        }

        ## Stable log predictive factor:
        ## lgamma(alpha + y) - lgamma(alpha) - lgamma(y + 1)
        ## + alpha log(beta_rate) + y log(z)
        ## - (alpha + y) log(beta_rate + z)
        lp <- lgamma(alpha + y_t) - lgamma(alpha) - lgamma(y_t + 1) +
              alpha * log(beta_rate) + y_t * log(z_t) -
              (alpha + y_t) * log(beta_rate + z_t)

        if (!is.finite(lp)) {
            if (return_filter) {
                return(list(log_ml = -Inf, a_filt = a_filt, b_filt = b_filt,
                            log_pred = log_pred, ok = FALSE,
                            reason = "nonfinite predictive factor"))
            }
            return(-Inf)
        }

        log_ml <- log_ml + lp

        at <- alpha + y_t
        bt <- beta_rate + z_t

        if (!is.finite(at) || !is.finite(bt) || at <= 0 || bt <= 0) {
            if (return_filter) {
                return(list(log_ml = -Inf, a_filt = a_filt, b_filt = b_filt,
                            log_pred = log_pred, ok = FALSE,
                            reason = "invalid updated filter state"))
            }
            return(-Inf)
        }

        if (return_filter) {
            a_filt[t] <- at
            b_filt[t] <- bt
            log_pred[t] <- lp
        }
    }

    if (return_filter) {
        return(list(log_ml = log_ml,
                    a_filt = a_filt,
                    b_filt = b_filt,
                    log_pred = log_pred,
                    ok = TRUE,
                    reason = NA_character_))
    }

    log_ml
}

## Backward-compatible alias, in case older files call this name.
log_marginal_gamma_lambda_collapsed <- log_marginal_gamma

## --------------------------------------------------------------------------
## MH update for gamma on the logit scale
## --------------------------------------------------------------------------

#' Update all gamma_j using logit-scale random-walk MH.
#'
#' @param gamma_current Length-n1 current discount factors.
#' @param y_coarse T x n1 count matrix.
#' @param xi T x n1 baseline offset matrix excluding kappa.
#' @param kappa T x n1 current NB random effect matrix.
#' @param a0 Initial Gamma shape.
#' @param b0 Initial Gamma rate.
#' @param priors List with gamma_a and gamma_b.
#' @param mh_sd Scalar or length-n1 proposal sd on logit scale.
#' @param eps Boundary tolerance for gamma and positive quantities.
#' @param return_diag If TRUE, return proposal diagnostics.
#' @return List with updated gamma and acceptance indicators.
update_gamma <- function(gamma_current, y_coarse, xi, kappa,
                         a0, b0, priors, mh_sd,
                         eps = 1e-12,
                         return_diag = TRUE) {
    check_gamma_inputs(gamma_current, y_coarse, xi, kappa,
                       a0, b0, priors, mh_sd)

    n1 <- ncol(y_coarse)
    gamma_new <- gamma_current
    accept <- logical(n1)

    gamma_proposal_vec <- rep(NA_real_, n1)
    log_alpha_vec <- rep(NA_real_, n1)
    log_ml_current_vec <- rep(NA_real_, n1)
    log_ml_proposal_vec <- rep(NA_real_, n1)
    log_target_current_vec <- rep(NA_real_, n1)
    log_target_proposal_vec <- rep(NA_real_, n1)

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]
        sd_j <- if (length(mh_sd) == 1) mh_sd else mh_sd[j]

        logit_current <- safe_logit(gamma_current[j], eps = eps)
        logit_proposal <- logit_current + rnorm(1, mean = 0, sd = sd_j)
        gamma_proposal <- as.numeric(safe_expit(logit_proposal))
        gamma_proposal <- pmin(pmax(gamma_proposal, eps), 1 - eps)
        gamma_proposal_vec[j] <- gamma_proposal

        log_ml_current <- log_marginal_gamma(gamma_current[j],
                                             y_coarse[, j], zeta_j,
                                             a0, b0, eps = eps)
        log_ml_proposal <- log_marginal_gamma(gamma_proposal,
                                              y_coarse[, j], zeta_j,
                                              a0, b0, eps = eps)

        log_ml_current_vec[j] <- log_ml_current
        log_ml_proposal_vec[j] <- log_ml_proposal

        log_prior_current <- dbeta(gamma_current[j], priors$gamma_a,
                                   priors$gamma_b, log = TRUE)
        log_prior_proposal <- dbeta(gamma_proposal, priors$gamma_a,
                                    priors$gamma_b, log = TRUE)

        ## Jacobian for transforming from gamma to logit(gamma):
        ## d gamma / d logit(gamma) = gamma * (1 - gamma).
        log_jac_current <- log(gamma_current[j]) + log1p(-gamma_current[j])
        log_jac_proposal <- log(gamma_proposal) + log1p(-gamma_proposal)

        log_target_current <- log_ml_current + log_prior_current + log_jac_current
        log_target_proposal <- log_ml_proposal + log_prior_proposal + log_jac_proposal

        log_target_current_vec[j] <- log_target_current
        log_target_proposal_vec[j] <- log_target_proposal

        log_alpha <- log_target_proposal - log_target_current
        log_alpha_vec[j] <- log_alpha

        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            gamma_new[j] <- gamma_proposal
            accept[j] <- TRUE
        }
    }

    out <- list(
        gamma = gamma_new,
        accept = accept
    )

    if (return_diag) {
        out$diag <- list(
            acceptance_rate = mean(accept),
            gamma_proposal = gamma_proposal_vec,
            log_alpha = log_alpha_vec,
            log_ml_current = log_ml_current_vec,
            log_ml_proposal = log_ml_proposal_vec,
            log_target_current = log_target_current_vec,
            log_target_proposal = log_target_proposal_vec,
            mean_gamma = mean(gamma_new),
            min_gamma = min(gamma_new),
            max_gamma = max(gamma_new)
        )
    }

    out
}

## Backward-compatible alias, if the sampler uses a generic discount name.
update_discount_gamma <- update_gamma
