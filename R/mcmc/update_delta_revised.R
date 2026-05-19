## ============================================================
## update_delta_revised.R
## ------------------------------------------------------------
## Update the multiscale splitting discount factor delta.
##
## Model block:
##   omega_{t,j} | D_{t-1}, delta ~ Dirichlet(delta * c_{t-1,j})
##   y_{t,j,1:K} | y_{t,j,+}, omega_{t,j} ~ Multinomial(y_{t,j,+}, omega_{t,j})
##
## The update below integrates out omega analytically and uses the
## Dirichlet-Multinomial one-step predictive likelihood.
##
## delta is updated on the logit scale by random-walk Metropolis-Hastings.
## The logit-scale Jacobian log(delta) + log(1 - delta) is included.
## ============================================================

## ---- numerical helpers ----
.safe_logit <- function(p, eps = 1e-12) {
    p <- pmin(pmax(p, eps), 1 - eps)
    log(p) - log1p(-p)
}

.safe_expit <- function(x) {
    out <- numeric(length(x))
    idx <- x >= 0
    out[idx] <- 1 / (1 + exp(-x[idx]))
    ex <- exp(x[!idx])
    out[!idx] <- ex / (1 + ex)
    out
}

.validate_delta_inputs <- function(delta_current, y_fine, c0, priors, mh_sd) {
    if (!is.numeric(delta_current) || length(delta_current) != 1L ||
        !is.finite(delta_current) || delta_current <= 0 || delta_current >= 1) {
        stop("delta_current must be a finite scalar strictly between 0 and 1.")
    }

    if (length(dim(y_fine)) != 3L) {
        stop("y_fine must be a 3-dimensional array with dimensions T x n1 x K.")
    }
    if (any(!is.finite(y_fine)) || any(y_fine < 0)) {
        stop("y_fine must contain finite nonnegative counts.")
    }
    if (any(abs(y_fine - round(y_fine)) > 1e-8)) {
        stop("y_fine must be integer-like count data.")
    }

    K <- dim(y_fine)[3]
    if (!is.numeric(c0) || length(c0) != K || any(!is.finite(c0)) || any(c0 <= 0)) {
        stop("c0 must be a finite positive numeric vector with length equal to dim(y_fine)[3].")
    }

    if (is.null(priors$delta_a) || is.null(priors$delta_b)) {
        stop("priors must contain delta_a and delta_b.")
    }
    if (!is.finite(priors$delta_a) || !is.finite(priors$delta_b) ||
        priors$delta_a <= 0 || priors$delta_b <= 0) {
        stop("priors$delta_a and priors$delta_b must be finite and positive.")
    }

    if (!is.numeric(mh_sd) || length(mh_sd) != 1L || !is.finite(mh_sd) || mh_sd <= 0) {
        stop("mh_sd must be a finite positive scalar.")
    }

    invisible(TRUE)
}

## ---- Dirichlet-Multinomial one-step log predictive ----
log_dirichlet_multinomial_predictive <- function(y_kids, c_prior) {
    if (any(!is.finite(y_kids)) || any(y_kids < 0)) return(-Inf)
    if (any(!is.finite(c_prior)) || any(c_prior <= 0)) return(-Inf)

    y_parent <- sum(y_kids)
    c_sum <- sum(c_prior)

    lgamma(c_sum) + lgamma(y_parent + 1) - lgamma(c_sum + y_parent) +
        sum(lgamma(c_prior + y_kids) - lgamma(c_prior) - lgamma(y_kids + 1))
}

## ---- marginal log likelihood for delta ----
log_marginal_delta <- function(delta, y_fine, c0,
                               eps = 1e-10,
                               return_diag = FALSE) {
    if (!is.numeric(delta) || length(delta) != 1L ||
        !is.finite(delta) || delta <= eps || delta >= 1 - eps) {
        if (return_diag) {
            return(list(
                log_marginal = -Inf,
                diag = list(reason = "delta outside stable open interval")
            ))
        }
        return(-Inf)
    }

    if (length(dim(y_fine)) != 3L) {
        stop("y_fine must be a 3-dimensional array with dimensions T x n1 x K.")
    }

    dims <- dim(y_fine)
    TT <- dims[1]
    n1 <- dims[2]
    K <- dims[3]

    if (length(c0) != K || any(!is.finite(c0)) || any(c0 <= 0)) {
        stop("c0 must be a finite positive vector with length equal to dim(y_fine)[3].")
    }

    ll <- 0
    min_c_prior <- Inf
    max_c_prior <- -Inf
    n_nonfinite_terms <- 0L

    for (j in seq_len(n1)) {
        ct <- as.numeric(c0)

        for (t in seq_len(TT)) {
            y_kids <- as.numeric(y_fine[t, j, ])
            c_prior <- delta * ct

            min_c_prior <- min(min_c_prior, c_prior)
            max_c_prior <- max(max_c_prior, c_prior)

            log_pred <- log_dirichlet_multinomial_predictive(
                y_kids = y_kids,
                c_prior = c_prior
            )

            if (!is.finite(log_pred)) {
                n_nonfinite_terms <- n_nonfinite_terms + 1L
                if (return_diag) {
                    return(list(
                        log_marginal = -Inf,
                        diag = list(
                            reason = "nonfinite Dirichlet-Multinomial predictive term",
                            t = t,
                            j = j,
                            min_c_prior = min_c_prior,
                            max_c_prior = max_c_prior,
                            n_nonfinite_terms = n_nonfinite_terms
                        )
                    ))
                }
                return(-Inf)
            }

            ll <- ll + log_pred
            ct <- c_prior + y_kids
        }
    }

    if (!return_diag) return(ll)

    list(
        log_marginal = ll,
        diag = list(
            min_c_prior = min_c_prior,
            max_c_prior = max_c_prior,
            n_nonfinite_terms = n_nonfinite_terms,
            total_parent_count = sum(y_fine),
            mean_parent_count = mean(apply(y_fine, c(1, 2), sum))
        )
    )
}

## ---- log posterior kernel on delta scale ----
log_posterior_delta <- function(delta, y_fine, c0, priors,
                                eps = 1e-10) {
    if (!is.numeric(delta) || length(delta) != 1L ||
        !is.finite(delta) || delta <= eps || delta >= 1 - eps) {
        return(-Inf)
    }

    log_ml <- log_marginal_delta(
        delta = delta,
        y_fine = y_fine,
        c0 = c0,
        eps = eps,
        return_diag = FALSE
    )

    if (!is.finite(log_ml)) return(-Inf)

    log_prior <- dbeta(
        delta,
        shape1 = priors$delta_a,
        shape2 = priors$delta_b,
        log = TRUE
    )

    log_ml + log_prior
}

## ---- update scalar delta by logit-scale MH ----
update_delta <- function(delta_current, y_fine, c0, priors, mh_sd,
                         eps = 1e-10,
                         return_diag = TRUE) {
    .validate_delta_inputs(
        delta_current = delta_current,
        y_fine = y_fine,
        c0 = c0,
        priors = priors,
        mh_sd = mh_sd
    )

    logit_current <- .safe_logit(delta_current, eps = eps)
    logit_proposal <- logit_current + rnorm(1, mean = 0, sd = mh_sd)
    delta_proposal <- as.numeric(.safe_expit(logit_proposal))

    ## Keep proposals away from exact numerical boundaries.
    if (!is.finite(delta_proposal) || delta_proposal <= eps || delta_proposal >= 1 - eps) {
        delta_new <- delta_current
        accept <- FALSE

        if (!return_diag) {
            return(list(delta = delta_new, accept = accept))
        }

        return(list(
            delta = delta_new,
            accept = accept,
            diag = list(
                acceptance_rate = as.numeric(accept),
                delta_proposal = delta_proposal,
                log_alpha = -Inf,
                log_ml_current = NA_real_,
                log_ml_proposal = -Inf,
                log_target_current = NA_real_,
                log_target_proposal = -Inf,
                reason = "proposal outside stable open interval"
            )
        ))
    }

    log_ml_current <- log_marginal_delta(
        delta = delta_current,
        y_fine = y_fine,
        c0 = c0,
        eps = eps,
        return_diag = FALSE
    )

    log_ml_proposal <- log_marginal_delta(
        delta = delta_proposal,
        y_fine = y_fine,
        c0 = c0,
        eps = eps,
        return_diag = FALSE
    )

    log_prior_current <- dbeta(
        delta_current,
        shape1 = priors$delta_a,
        shape2 = priors$delta_b,
        log = TRUE
    )
    log_prior_proposal <- dbeta(
        delta_proposal,
        shape1 = priors$delta_a,
        shape2 = priors$delta_b,
        log = TRUE
    )

    ## Jacobian for the logit transformation:
    ## d delta / d logit(delta) = delta * (1 - delta)
    log_jac_current <- log(delta_current) + log1p(-delta_current)
    log_jac_proposal <- log(delta_proposal) + log1p(-delta_proposal)

    log_target_current <- log_ml_current + log_prior_current + log_jac_current
    log_target_proposal <- log_ml_proposal + log_prior_proposal + log_jac_proposal
    log_alpha <- log_target_proposal - log_target_current

    accept <- is.finite(log_alpha) && (log(runif(1)) < log_alpha)
    delta_new <- if (accept) delta_proposal else delta_current

    if (!return_diag) {
        return(list(delta = delta_new, accept = accept))
    }

    list(
        delta = delta_new,
        accept = accept,
        diag = list(
            acceptance_rate = as.numeric(accept),
            delta_proposal = delta_proposal,
            log_alpha = log_alpha,
            log_ml_current = log_ml_current,
            log_ml_proposal = log_ml_proposal,
            log_prior_current = log_prior_current,
            log_prior_proposal = log_prior_proposal,
            log_target_current = log_target_current,
            log_target_proposal = log_target_proposal,
            mean_delta = delta_new,
            min_delta = delta_new,
            max_delta = delta_new
        )
    )
}

## Backward-compatible alias, in case older sampler files use this name.
update_split_discount <- update_delta
