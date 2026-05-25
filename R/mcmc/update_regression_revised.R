## ============================================================================
## update_regression_revised.R
## Clean package-level regression update for MSSTNB.
##
## Main cleanup:
##   update_beta() now formally supports optional lambda_tilde.
##   If lambda_tilde = NULL, it defaults to a matrix of ones with the same
##   dimensions as y_coarse.  This makes the same updater compatible with:
##     - Scenario 1: fixed residual risk, lambda_tilde = 1
##     - Scenario 2: sampled dynamic residual risk, gamma fixed
##     - Full model: sampled dynamic residual risk, gamma learned
##
## Notes:
##   The beta update uses an elliptical slice sampler under the Gaussian prior.
##   The likelihood is the marginal NB likelihood with kappa integrated out,
##   conditional on the current lambda_tilde path.
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.ensure_matrix_same_dim <- function(obj, ref, obj_name = deparse(substitute(obj)),
                                    ref_name = deparse(substitute(ref))) {
    if (!is.matrix(obj) || !is.matrix(ref) || !identical(dim(obj), dim(ref))) {
        stop(obj_name, " must be a matrix with the same dimensions as ", ref_name, ".",
             call. = FALSE)
    }
    invisible(TRUE)
}

.prepare_lambda_tilde <- function(lambda_tilde, y_coarse) {
    if (is.null(lambda_tilde)) {
        return(matrix(1, nrow = nrow(y_coarse), ncol = ncol(y_coarse)))
    }
    .ensure_matrix_same_dim(lambda_tilde, y_coarse, "lambda_tilde", "y_coarse")
    if (any(!is.finite(lambda_tilde)) || any(lambda_tilde <= 0)) {
        stop("lambda_tilde must be positive and finite.", call. = FALSE)
    }
    lambda_tilde
}

.check_beta_inputs <- function(beta_current, y_coarse, e, x1, x2,
                               lambda_tilde = NULL, phi, r, priors) {
    if (!is.matrix(y_coarse)) {
        stop("y_coarse must be a matrix.", call. = FALSE)
    }
    .ensure_matrix_same_dim(e, y_coarse, "e", "y_coarse")
    .ensure_matrix_same_dim(x1, y_coarse, "x1", "y_coarse")
    .ensure_matrix_same_dim(x2, y_coarse, "x2", "y_coarse")
    lambda_tilde <- .prepare_lambda_tilde(lambda_tilde, y_coarse)

    n1 <- ncol(y_coarse)
    if (length(phi) != n1) {
        stop("phi must have length ncol(y_coarse).", call. = FALSE)
    }
    if (!(length(r) == 1L || length(r) == n1)) {
        stop("r must be scalar or have length ncol(y_coarse).", call. = FALSE)
    }
    if (length(r) == 1L) {
        r <- rep(r, n1)
    }
    if (any(!is.finite(r)) || any(r <= 0)) {
        stop("r must be positive and finite.", call. = FALSE)
    }

    required_priors <- c("beta0_mean", "beta0_sd", "beta_mean", "beta_sd")
    missing <- setdiff(required_priors, names(priors))
    if (length(missing) > 0L) {
        stop("priors is missing fields: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    if (length(priors$beta_mean) != length(priors$beta_sd)) {
        stop("priors$beta_mean and priors$beta_sd must have the same length.", call. = FALSE)
    }

    p <- length(priors$beta_mean)
    if (length(beta_current) != p + 1L) {
        stop("beta_current must have length 1 + length(priors$beta_mean).", call. = FALSE)
    }
    if (length(priors$beta0_sd) != 1L || priors$beta0_sd <= 0 ||
        any(priors$beta_sd <= 0)) {
        stop("beta prior standard deviations must be positive.", call. = FALSE)
    }

    list(lambda_tilde = lambda_tilde, r = r, p = p)
}

.nb_or_poisson_loglik <- function(y, log_mu, r) {
    log_mu <- pmin(pmax(log_mu, -700), 700)
    mu <- pmax(exp(log_mu), .Machine$double.xmin)
    if (is.null(r) || !is.finite(r) || r <= 0 || is.infinite(r)) {
        return(sum(stats::dpois(y, lambda = mu, log = TRUE)))
    }
    sum(stats::dnbinom(y, size = r, mu = mu, log = TRUE))
}

.compute_beta_nb_loglik <- function(beta_full, y_coarse, e, x1, x2,
                                    lambda_tilde, phi, r) {
    beta0 <- beta_full[1]
    beta <- beta_full[-1]
    n1 <- ncol(y_coarse)
    out <- 0
    for (j in seq_len(n1)) {
        log_mu_j <- log(e[, j]) + log(lambda_tilde[, j]) + beta0 +
            beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        out <- out + .nb_or_poisson_loglik(y = y_coarse[, j],
                                           log_mu = log_mu_j,
                                           r = r[j])
    }
    out
}

.logpost_beta <- function(beta_full, y_coarse, e, x1, x2, lambda_tilde,
                          phi, r, priors) {
    ll <- .compute_beta_nb_loglik(
        beta_full = beta_full,
        y_coarse = y_coarse,
        e = e,
        x1 = x1,
        x2 = x2,
        lambda_tilde = lambda_tilde,
        phi = phi,
        r = r
    )

    prior_mean <- c(priors$beta0_mean, priors$beta_mean)
    prior_sd <- c(priors$beta0_sd, priors$beta_sd)
    lp <- sum(stats::dnorm(beta_full, mean = prior_mean, sd = prior_sd, log = TRUE))
    ll + lp
}

fit_glm_for_ess <- function(dat, lambda_tilde = NULL) {
    if (is.null(dat$y_coarse) || is.null(dat$e) || is.null(dat$x1) || is.null(dat$x2)) {
        stop("dat must contain y_coarse, e, x1, and x2.", call. = FALSE)
    }
    y <- dat$y_coarse
    e <- dat$e
    x1 <- dat$x1
    x2 <- dat$x2
    lambda_tilde <- .prepare_lambda_tilde(lambda_tilde %||% dat$lambda_tilde, y)

    df <- data.frame(
        y = as.vector(y),
        x1 = as.vector(x1),
        x2 = as.vector(x2),
        offset_val = log(as.vector(e)) + log(as.vector(lambda_tilde))
    )

    fit <- try(
        stats::glm(y ~ x1 + x2, family = stats::poisson(),
                   offset = offset_val, data = df),
        silent = TRUE
    )

    if (inherits(fit, "try-error") || any(!is.finite(stats::coef(fit)))) {
        center <- c(0, 0, 0)
        vcov <- diag(100, 3)
    } else {
        center <- as.numeric(stats::coef(fit))
        vcov <- try(stats::vcov(fit), silent = TRUE)
        if (inherits(vcov, "try-error") || any(!is.finite(vcov))) {
            vcov <- diag(100, length(center))
        }
    }

    list(
        center = center,
        vcov = as.matrix(vcov),
        fit = fit
    )
}

update_beta <- function(beta_current, y_coarse, e, x1, x2,
                        lambda_tilde = NULL, phi, priors, ess_base = NULL, r,
                        use_preconditioned = TRUE, max_steps = 200L) {
    checked <- .check_beta_inputs(
        beta_current = beta_current,
        y_coarse = y_coarse,
        e = e,
        x1 = x1,
        x2 = x2,
        lambda_tilde = lambda_tilde,
        phi = phi,
        r = r,
        priors = priors
    )
    lambda_tilde <- checked$lambda_tilde
    r <- checked$r

    prior_mean <- c(priors$beta0_mean, priors$beta_mean)
    prior_sd <- c(priors$beta0_sd, priors$beta_sd)

    log_current <- .logpost_beta(
        beta_current, y_coarse, e, x1, x2, lambda_tilde, phi, r, priors
    )
    log_y <- log_current + log(stats::runif(1))

    ## Exact prior elliptical slice sampler.  The ess_base argument is kept for
    ## backward compatibility with existing Scenario 1/full-model state objects.
    nu <- stats::rnorm(length(beta_current), mean = 0, sd = prior_sd)
    theta <- stats::runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    centered <- beta_current - prior_mean
    n_reject <- 0L

    for (step in seq_len(max_steps)) {
        proposal <- prior_mean + centered * cos(theta) + nu * sin(theta)
        log_prop <- .logpost_beta(
            proposal, y_coarse, e, x1, x2, lambda_tilde, phi, r, priors
        )
        if (is.finite(log_prop) && log_prop >= log_y) {
            return(list(
                sample = proposal,
                n_reject = n_reject,
                accepted = TRUE,
                log_target = log_prop,
                log_lik = .compute_beta_nb_loglik(proposal, y_coarse, e, x1, x2,
                                                  lambda_tilde, phi, r),
                log_marginal_lik = .compute_beta_nb_loglik(proposal, y_coarse, e, x1, x2,
                                                           lambda_tilde, phi, r),
                ess_mode = "prior_ess_lambda_optional"
            ))
        }
        n_reject <- n_reject + 1L
        if (theta < 0) {
            theta_min <- theta
        } else {
            theta_max <- theta
        }
        theta <- stats::runif(1, theta_min, theta_max)
    }

    list(
        sample = beta_current,
        n_reject = n_reject,
        accepted = FALSE,
        log_target = log_current,
        log_lik = .compute_beta_nb_loglik(beta_current, y_coarse, e, x1, x2,
                                          lambda_tilde, phi, r),
        log_marginal_lik = .compute_beta_nb_loglik(beta_current, y_coarse, e, x1, x2,
                                                   lambda_tilde, phi, r),
        ess_mode = "prior_ess_lambda_optional"
    )
}
