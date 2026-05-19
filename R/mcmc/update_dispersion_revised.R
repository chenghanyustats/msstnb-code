## -----------------------------------------------------------------------------
## update_dispersion_revised.R
##
## Update region-specific negative-binomial dispersion parameters r_j.
##
## Main update: marginal negative-binomial update for r_j
##     y_{t,j} | mu_{t,j}, r_j ~ NB(mean = mu_{t,j}, size = r_j)
##
## where
##     mu_{t,j} = e_{t,j} * lambda_tilde_{t,j}
##               * exp(beta0 + beta1*x1_{t,j} + beta2*x2_{t,j} + phi_j).
##
## This update integrates out kappa_{t,j}, matching the kappa-collapsed
## beta and phi updates used in update_regression.R and update_icar.R.
##
## Optional legacy update: conditional-on-kappa update for r_j
##     kappa_{t,j} | r_j ~ Gamma(shape = r_j, rate = r_j).
##
## Both updates use a random-walk Metropolis step on log(r_j), with the
## required Jacobian adjustment + log(r_j) in the log target.
## -----------------------------------------------------------------------------

## Numerically stable log(exp(a) + exp(b)), vectorized.
logspace_add <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
}

check_dispersion_priors <- function(priors) {
    if (is.null(priors$r_shape) || is.null(priors$r_rate)) {
        stop("priors must contain r_shape and r_rate.")
    }
    if (!is.finite(priors$r_shape) || priors$r_shape <= 0) {
        stop("priors$r_shape must be positive and finite.")
    }
    if (!is.finite(priors$r_rate) || priors$r_rate <= 0) {
        stop("priors$r_rate must be positive and finite.")
    }
    invisible(TRUE)
}

check_mh_sd <- function(mh_sd, n1) {
    if (!(length(mh_sd) %in% c(1L, n1))) {
        stop("mh_sd must be either a scalar or a vector of length ncol(y_coarse).")
    }
    if (any(!is.finite(mh_sd)) || any(mh_sd <= 0)) {
        stop("All mh_sd values must be positive and finite.")
    }
    invisible(TRUE)
}

check_matrix_same_dim <- function(x, ref, name, ref_name = "y_coarse") {
    if (is.null(dim(x)) || length(dim(x)) != 2L) {
        stop(sprintf("%s must be a matrix.", name))
    }
    if (!identical(dim(x), dim(ref))) {
        stop(sprintf("%s must have the same dimension as %s.", name, ref_name))
    }
    invisible(TRUE)
}

check_marginal_nb_inputs <- function(y_coarse, e, x1, x2, lambda_tilde,
                                     beta0, beta, phi, r_current,
                                     priors, mh_sd) {
    if (is.null(dim(y_coarse)) || length(dim(y_coarse)) != 2L) {
        stop("y_coarse must be a matrix.")
    }

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    check_matrix_same_dim(e, y_coarse, "e")
    check_matrix_same_dim(x1, y_coarse, "x1")
    check_matrix_same_dim(x2, y_coarse, "x2")
    check_matrix_same_dim(lambda_tilde, y_coarse, "lambda_tilde")

    if (length(r_current) != n1) {
        stop("length(r_current) must equal ncol(y_coarse).")
    }
    if (length(phi) != n1) {
        stop("length(phi) must equal ncol(y_coarse).")
    }
    if (length(beta) < 2L) {
        stop("beta must contain at least beta[1] and beta[2].")
    }
    if (!is.finite(beta0) || any(!is.finite(beta[1:2])) || any(!is.finite(phi))) {
        stop("beta0, beta[1:2], and phi must be finite.")
    }
    if (any(!is.finite(y_coarse)) || any(y_coarse < 0) || any(abs(y_coarse - round(y_coarse)) > 1e-8)) {
        stop("y_coarse must contain nonnegative integer counts.")
    }
    if (any(!is.finite(e)) || any(e <= 0)) {
        stop("e must be positive and finite.")
    }
    if (any(!is.finite(lambda_tilde)) || any(lambda_tilde <= 0)) {
        stop("lambda_tilde must be positive and finite.")
    }
    if (any(!is.finite(x1)) || any(!is.finite(x2))) {
        stop("x1 and x2 must be finite.")
    }
    if (any(!is.finite(r_current)) || any(r_current <= 0)) {
        stop("r_current must be positive and finite for the marginal NB update.")
    }

    check_dispersion_priors(priors)
    check_mh_sd(mh_sd, n1)

    invisible(list(TT = TT, n1 = n1))
}

## -----------------------------------------------------------------------------
## Marginal NB log likelihood for one region j as a function of r_j.
##
## NB parameterization:
##     E(Y) = mu,
##     Var(Y) = mu + mu^2 / r.
##
## log p(y | mu, r)
##   = lgamma(y + r) - lgamma(r) - lgamma(y + 1)
##     + r * log(r)
##     + y * log(mu)
##     - (y + r) * log(r + mu).
##
## This function works on the log-mu scale to avoid overflow from exp(eta).
## -----------------------------------------------------------------------------
loglik_r_marginal_nb <- function(r, y_j, log_mu_j) {
    if (!is.finite(r) || r <= 0) {
        return(-Inf)
    }
    if (any(!is.finite(y_j)) || any(y_j < 0)) {
        return(-Inf)
    }
    if (any(!is.finite(log_mu_j))) {
        return(-Inf)
    }

    log_r <- log(r)
    log_r_plus_mu <- logspace_add(log_r, log_mu_j)

    sum(lgamma(y_j + r) - lgamma(r) - lgamma(y_j + 1) +
        r * log_r + y_j * log_mu_j - (y_j + r) * log_r_plus_mu)
}

log_posterior_r_marginal_nb <- function(r, y_j, log_mu_j, priors) {
    if (!is.finite(r) || r <= 0) {
        return(-Inf)
    }

    log_prior <- (priors$r_shape - 1) * log(r) - priors$r_rate * r
    log_lik <- loglik_r_marginal_nb(r, y_j, log_mu_j)

    log_prior + log_lik
}

## -----------------------------------------------------------------------------
## Legacy conditional-on-kappa log posterior for one region j.
##
## kappa_{t,j} | r_j ~ Gamma(shape = r_j, rate = r_j)
## r_j ~ Gamma(r_shape, r_rate)
## -----------------------------------------------------------------------------
log_posterior_r_conditional_kappa <- function(r, kappa_j, priors) {
    if (!is.finite(r) || r <= 0) {
        return(-Inf)
    }
    if (any(!is.finite(kappa_j)) || any(kappa_j <= 0)) {
        return(-Inf)
    }

    TT <- length(kappa_j)
    sum_log_kappa <- sum(log(kappa_j))
    sum_kappa <- sum(kappa_j)

    log_prior <- (priors$r_shape - 1) * log(r) - priors$r_rate * r
    log_lik <- TT * (r * log(r) - lgamma(r)) +
        (r - 1) * sum_log_kappa - r * sum_kappa

    log_prior + log_lik
}

## -----------------------------------------------------------------------------
## Main updater.
##
## method = "marginal_nb" is the recommended default.
## method = "conditional_kappa" reproduces the older conditional update and
## requires kappa.
##
## r_min and r_max impose a proper numerical support interval for r. The proposal
## is rejected if it falls outside this interval. A large r approximates Poisson.
## -----------------------------------------------------------------------------
update_r <- function(r_current,
                     y_coarse = NULL,
                     e = NULL,
                     x1 = NULL,
                     x2 = NULL,
                     lambda_tilde = NULL,
                     beta0 = NULL,
                     beta = NULL,
                     phi = NULL,
                     priors,
                     mh_sd,
                     kappa = NULL,
                     method = c("marginal_nb", "conditional_kappa"),
                     r_min = 1e-6,
                     r_max = 1e4,
                     return_diag = TRUE) {

    method <- match.arg(method)
    check_dispersion_priors(priors)

    if (any(!is.finite(r_current)) || any(r_current <= 0)) {
        stop("r_current must be positive and finite.")
    }
    if (!is.finite(r_min) || r_min <= 0) {
        stop("r_min must be positive and finite.")
    }
    if (!is.finite(r_max) || r_max <= r_min) {
        stop("r_max must be finite and larger than r_min.")
    }

    if (method == "marginal_nb") {
        dims <- check_marginal_nb_inputs(
            y_coarse = y_coarse,
            e = e,
            x1 = x1,
            x2 = x2,
            lambda_tilde = lambda_tilde,
            beta0 = beta0,
            beta = beta,
            phi = phi,
            r_current = r_current,
            priors = priors,
            mh_sd = mh_sd
        )
        TT <- dims$TT
        n1 <- dims$n1
    } else {
        if (is.null(kappa)) {
            stop("kappa must be provided when method = 'conditional_kappa'.")
        }
        if (is.null(dim(kappa)) || length(dim(kappa)) != 2L) {
            stop("kappa must be a matrix.")
        }
        TT <- nrow(kappa)
        n1 <- ncol(kappa)
        if (length(r_current) != n1) {
            stop("length(r_current) must equal ncol(kappa).")
        }
        if (any(!is.finite(kappa)) || any(kappa <= 0)) {
            stop("kappa must be positive and finite.")
        }
        check_mh_sd(mh_sd, n1)
    }

    r_new <- r_current
    accept <- rep(FALSE, n1)
    log_alpha_vec <- rep(NA_real_, n1)
    r_proposal_vec <- rep(NA_real_, n1)
    log_target_current_vec <- rep(NA_real_, n1)
    log_target_proposal_vec <- rep(NA_real_, n1)

    for (j in seq_len(n1)) {
        sd_j <- if (length(mh_sd) > 1L) mh_sd[j] else mh_sd

        log_r_current <- log(r_current[j])
        log_r_proposal <- log_r_current + rnorm(1, mean = 0, sd = sd_j)
        r_proposal <- exp(log_r_proposal)
        r_proposal_vec[j] <- r_proposal

        if (!is.finite(r_proposal) || r_proposal < r_min || r_proposal > r_max) {
            log_alpha_vec[j] <- -Inf
            next
        }

        if (method == "marginal_nb") {
            eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
            log_mu_j <- log(e[, j]) + log(lambda_tilde[, j]) + eta_j

            lp_current <- log_posterior_r_marginal_nb(
                r = r_current[j],
                y_j = y_coarse[, j],
                log_mu_j = log_mu_j,
                priors = priors
            )
            lp_proposal <- log_posterior_r_marginal_nb(
                r = r_proposal,
                y_j = y_coarse[, j],
                log_mu_j = log_mu_j,
                priors = priors
            )
        } else {
            lp_current <- log_posterior_r_conditional_kappa(
                r = r_current[j],
                kappa_j = kappa[, j],
                priors = priors
            )
            lp_proposal <- log_posterior_r_conditional_kappa(
                r = r_proposal,
                kappa_j = kappa[, j],
                priors = priors
            )
        }

        ## Jacobian adjustment for sampling s = log(r).
        log_target_current <- lp_current + log_r_current
        log_target_proposal <- lp_proposal + log_r_proposal

        log_target_current_vec[j] <- log_target_current
        log_target_proposal_vec[j] <- log_target_proposal

        log_alpha <- log_target_proposal - log_target_current
        log_alpha_vec[j] <- log_alpha

        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            r_new[j] <- r_proposal
            accept[j] <- TRUE
        }
    }

    if (!return_diag) {
        return(list(r = r_new, accept = accept))
    }

    list(
        r = r_new,
        accept = accept,
        diag = list(
            method = method,
            acceptance_rate = mean(accept),
            log_alpha = log_alpha_vec,
            r_proposal = r_proposal_vec,
            log_target_current = log_target_current_vec,
            log_target_proposal = log_target_proposal_vec,
            r_min = r_min,
            r_max = r_max,
            mean_r = mean(r_new),
            min_r = min(r_new),
            max_r = max(r_new)
        )
    )
}

## Backward-compatible wrapper name, if older sampler code calls update_dispersion().
update_dispersion <- update_r
