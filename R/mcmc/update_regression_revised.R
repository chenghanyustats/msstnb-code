## ======================================================================
## update_regression_revised.R
## Beta update using a kappa-collapsed marginal negative-binomial likelihood.
##
## Main design:
##   beta is intentionally NOT updated conditional on sampled kappa.
##   Instead, kappa is integrated out analytically, so beta targets the
##   observed-data marginal likelihood conditional on lambda_tilde, phi, and r.
##
## Model component:
##   y_tj | beta, lambda_tilde_tj, phi_j, r_j ~ NB(mean = mu_tj, size = r_j)
##   mu_tj = e_tj * lambda_tilde_tj * exp(beta0 + beta1*x1_tj + beta2*x2_tj + phi_j)
## ======================================================================

flatten_region_major <- function(mat) {
    as.vector(mat)
}

.check_same_dim <- function(reference, ..., names = NULL) {
    objects <- list(...)
    if (is.null(names)) {
        names <- paste0("object", seq_along(objects))
    }
    ref_dim <- dim(reference)
    for (i in seq_along(objects)) {
        if (!identical(dim(objects[[i]]), ref_dim)) {
            stop(names[i], " must have the same dimensions as y_coarse.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

.check_positive_finite_matrix <- function(x, name, allow_zero = FALSE) {
    if (any(!is.finite(x))) {
        stop(name, " must contain only finite values.", call. = FALSE)
    }
    if (allow_zero) {
        if (any(x < 0)) stop(name, " must be nonnegative.", call. = FALSE)
    } else {
        if (any(x <= 0)) stop(name, " must be strictly positive.", call. = FALSE)
    }
    invisible(TRUE)
}

.check_beta_inputs <- function(beta_current, y_coarse, e, x1, x2,
                               lambda_tilde, phi, r, priors = NULL) {
    if (length(beta_current) != 3L) {
        stop("beta_current must have length 3: c(beta0, beta1, beta2).", call. = FALSE)
    }
    if (!is.matrix(y_coarse)) stop("y_coarse must be a matrix.", call. = FALSE)
    if (any(!is.finite(y_coarse)) || any(y_coarse < 0)) {
        stop("y_coarse must contain nonnegative finite counts.", call. = FALSE)
    }
    .check_same_dim(y_coarse, e, x1, x2, lambda_tilde,
                    names = c("e", "x1", "x2", "lambda_tilde"))
    .check_positive_finite_matrix(e, "e")
    .check_positive_finite_matrix(lambda_tilde, "lambda_tilde")
    if (any(!is.finite(x1)) || any(!is.finite(x2))) {
        stop("x1 and x2 must contain only finite values.", call. = FALSE)
    }
    n1 <- ncol(y_coarse)
    if (length(phi) != n1 || any(!is.finite(phi))) {
        stop("phi must be a finite vector with length ncol(y_coarse).", call. = FALSE)
    }
    if (is.null(r) || length(r) != n1) {
        stop("r must be supplied and have length ncol(y_coarse).", call. = FALSE)
    }
    if (any((r <= 0 | !is.finite(r)) & !is.infinite(r))) {
        stop("Each r[j] must be positive finite or Inf.", call. = FALSE)
    }
    if (!is.null(priors)) {
        prior_mean <- c(priors$beta0_mean, priors$beta_mean)
        prior_sd <- c(priors$beta0_sd, priors$beta_sd)
        if (length(prior_mean) != 3L || any(!is.finite(prior_mean))) {
            stop("Priors must define finite beta0_mean and beta_mean of total length 3.", call. = FALSE)
        }
        if (length(prior_sd) != 3L || any(!is.finite(prior_sd)) || any(prior_sd <= 0)) {
            stop("Priors must define positive finite beta0_sd and beta_sd of total length 3.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

.nb_or_poisson_loglik <- function(y, log_mu, r) {
    if (any(!is.finite(log_mu))) {
        return(-Inf)
    }

    if (is.infinite(r)) {
        mu <- exp(log_mu)
        if (any(!is.finite(mu))) return(-Inf)
        return(sum(y * log_mu - mu - lgamma(y + 1)))
    }

    if (!is.finite(r) || r <= 0) return(-Inf)

    log_r <- log(r)
    log_r_plus_mu <- pmax(log_r, log_mu) + log1p(exp(-abs(log_r - log_mu)))

    sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
            r * log_r + y * log_mu - (y + r) * log_r_plus_mu)
}

loglik_beta_marginal <- function(beta_vec, y_coarse, e, x1, x2,
                                 lambda_tilde, phi, r) {
    if (length(beta_vec) != 3L || any(!is.finite(beta_vec))) {
        return(-Inf)
    }

    beta0 <- beta_vec[1]
    beta1 <- beta_vec[2]
    beta2 <- beta_vec[3]
    n1 <- ncol(y_coarse)
    ll <- 0

    for (j in seq_len(n1)) {
        log_mu_j <- log(e[, j]) + log(lambda_tilde[, j]) +
            beta0 + beta1 * x1[, j] + beta2 * x2[, j] + phi[j]
        ll_j <- .nb_or_poisson_loglik(y = y_coarse[, j], log_mu = log_mu_j, r = r[j])
        if (!is.finite(ll_j)) return(-Inf)
        ll <- ll + ll_j
    }

    ll
}

ess_step <- function(current, prior_sample, log_lik_fn, prior_mean,
                     max_shrink = 1000L, min_width = 1e-12) {
    nu <- prior_mean + prior_sample
    cur_ll <- log_lik_fn(current)
    if (!is.finite(cur_ll)) {
        stop("Current value has non-finite ESS log target.", call. = FALSE)
    }

    log_y <- cur_ll + log(runif(1))
    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    n_reject <- 0L

    repeat {
        proposal <- prior_mean + cos(theta) * (current - prior_mean) +
            sin(theta) * (nu - prior_mean)
        prop_ll <- log_lik_fn(proposal)

        if (is.finite(prop_ll) && prop_ll >= log_y) {
            return(list(sample = proposal,
                        log_target = prop_ll,
                        n_reject = n_reject,
                        accepted = TRUE))
        }

        n_reject <- n_reject + 1L
        if (theta < 0) {
            theta_min <- theta
        } else {
            theta_max <- theta
        }

        if (n_reject >= max_shrink || theta_max - theta_min < min_width) {
            return(list(sample = current,
                        log_target = cur_ll,
                        n_reject = n_reject,
                        accepted = FALSE))
        }

        theta <- runif(1, theta_min, theta_max)
    }
}

fit_glm_for_ess <- function(dat,
                            min_sd = c(1.0, 0.25, 0.25),
                            inflate = 10) {
    y_vec  <- flatten_region_major(dat$y_coarse)
    x1_vec <- flatten_region_major(dat$x1)
    x2_vec <- flatten_region_major(dat$x2)
    e_vec <- pmax(flatten_region_major(dat$e), 1e-12)
    offset_vec <- log(e_vec)

    glm_fit <- try(suppressWarnings(glm(y_vec ~ x1_vec + x2_vec,
                                        family = poisson(link = "log"),
                                        offset = offset_vec,
                                        control = glm.control(maxit = 50))),
                   silent = TRUE)

    ok <- !inherits(glm_fit, "try-error")
    if (ok) {
        beta_hat <- as.numeric(coef(glm_fit))
        V <- try(vcov(glm_fit), silent = TRUE)
        ok <- length(beta_hat) == 3L && all(is.finite(beta_hat)) && !inherits(V, "try-error")
        if (ok) {
            beta_se <- sqrt(diag(V))
            ok <- length(beta_se) == 3L && all(is.finite(beta_se)) && all(beta_se > 0)
        }
    }

    if (!ok) {
        beta_hat <- c(log(mean(y_vec + 0.1) / mean(e_vec)), 0, 0)
        beta_se <- c(1, 1, 1)
    }

    list(center = beta_hat,
         scale = pmax(inflate * beta_se, min_sd))
}

build_beta_ess_base <- function(beta_current,
                                y_coarse, e, x1, x2,
                                lambda_tilde, phi,
                                fallback_base = NULL,
                                min_sd = c(1.0, 0.25, 0.25),
                                inflate = 10) {
    TT <- nrow(y_coarse)
    y_vec <- flatten_region_major(y_coarse)
    x1_vec <- flatten_region_major(x1)
    x2_vec <- flatten_region_major(x2)
    phi_vec <- rep(phi, each = TT)

    offset_core <- pmax(e * lambda_tilde, 1e-12)
    offset_vec <- log(flatten_region_major(offset_core)) + phi_vec

    glm_fit <- try(suppressWarnings(glm(y_vec ~ x1_vec + x2_vec,
                                        family = poisson(link = "log"),
                                        offset = offset_vec,
                                        control = glm.control(maxit = 50))),
                   silent = TRUE)

    ok <- !inherits(glm_fit, "try-error")
    if (ok) {
        beta_hat <- as.numeric(coef(glm_fit))
        V <- try(vcov(glm_fit), silent = TRUE)
        ok <- length(beta_hat) == 3L && all(is.finite(beta_hat)) && !inherits(V, "try-error")
        if (ok) {
            beta_se <- sqrt(diag(V))
            ok <- length(beta_se) == 3L && all(is.finite(beta_se)) && all(beta_se > 0)
        }
    }

    if (!ok) {
        if (!is.null(fallback_base) &&
            length(fallback_base$center) == 3L &&
            length(fallback_base$scale) == 3L &&
            all(is.finite(fallback_base$center)) &&
            all(is.finite(fallback_base$scale)) &&
            all(fallback_base$scale > 0)) {
            beta_hat <- as.numeric(fallback_base$center)
            beta_se <- as.numeric(fallback_base$scale)
            inflate <- 1
        } else {
            beta_hat <- as.numeric(beta_current)
            beta_se <- c(1, 0.5, 0.5)
            inflate <- 1
        }
    }

    list(center = beta_hat,
         scale = pmax(inflate * beta_se, min_sd))
}

update_beta <- function(beta_current, y_coarse, e, x1, x2,
                        kappa = NULL, lambda_tilde, phi, priors,
                        ess_base = NULL, r = NULL,
                        proposal_sd = NULL,
                        use_preconditioned = TRUE,
                        ess_min_sd = c(1.0, 0.25, 0.25),
                        ess_inflate = 10,
                        return_diagnostics = TRUE,
                        ...) {
    ## kappa is intentionally ignored.  The beta update uses the
    ## kappa-collapsed marginal NB likelihood.
    invisible(kappa)
    invisible(proposal_sd)

    .check_beta_inputs(beta_current, y_coarse, e, x1, x2,
                       lambda_tilde, phi, r, priors)

    d <- length(beta_current)
    prior_mean <- c(priors$beta0_mean, priors$beta_mean)
    prior_sd <- c(priors$beta0_sd, priors$beta_sd)

    log_true_prior <- function(beta_vec) {
        if (any(!is.finite(beta_vec))) return(-Inf)
        -0.5 * sum(((beta_vec - prior_mean) / prior_sd)^2)
    }

    if (use_preconditioned) {
        base <- build_beta_ess_base(beta_current = beta_current,
                                    y_coarse = y_coarse,
                                    e = e,
                                    x1 = x1,
                                    x2 = x2,
                                    lambda_tilde = lambda_tilde,
                                    phi = phi,
                                    fallback_base = ess_base,
                                    min_sd = ess_min_sd,
                                    inflate = ess_inflate)
        ess_center <- base$center
        ess_scale <- base$scale

        log_ess_base <- function(beta_vec) {
            if (any(!is.finite(beta_vec))) return(-Inf)
            -0.5 * sum(((beta_vec - ess_center) / ess_scale)^2)
        }

        log_target_fn <- function(beta_vec) {
            ll <- loglik_beta_marginal(beta_vec, y_coarse, e, x1, x2,
                                       lambda_tilde, phi, r)
            if (!is.finite(ll)) return(-Inf)
            ll + log_true_prior(beta_vec) - log_ess_base(beta_vec)
        }

        out <- ess_step(current = beta_current,
                        prior_sample = rnorm(d, mean = 0, sd = ess_scale),
                        log_lik_fn = log_target_fn,
                        prior_mean = ess_center)

        out$ess_center <- ess_center
        out$ess_scale <- ess_scale
        out$ess_mode <- "preconditioned_collapsed_kappa"
        out$log_marginal_lik <- loglik_beta_marginal(out$sample, y_coarse, e, x1, x2,
                                                     lambda_tilde, phi, r)
        out$log_true_prior <- log_true_prior(out$sample)
        out$log_ess_base <- log_ess_base(out$sample)
        ## Backward-compatible alias.  Prefer log_target in new diagnostics.
        out$log_lik <- out$log_target
        return(out)
    }

    log_target_fn <- function(beta_vec) {
        loglik_beta_marginal(beta_vec, y_coarse, e, x1, x2,
                             lambda_tilde, phi, r)
    }

    out <- ess_step(current = beta_current,
                    prior_sample = rnorm(d, mean = 0, sd = prior_sd),
                    log_lik_fn = log_target_fn,
                    prior_mean = prior_mean)

    out$ess_center <- prior_mean
    out$ess_scale <- prior_sd
    out$ess_mode <- "ordinary_collapsed_kappa"
    out$log_marginal_lik <- out$log_target
    out$log_true_prior <- log_true_prior(out$sample)
    out$log_ess_base <- NA_real_
    ## Backward-compatible alias.
    out$log_lik <- out$log_target
    out
}
