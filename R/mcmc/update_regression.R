## ======================================================================
## update_regression.R
## Beta update using a kappa-collapsed marginal NB likelihood.
##
## Key design change:
##   beta is NO LONGER updated conditional on sampled kappa.
##   Instead, kappa is integrated out analytically, so beta targets the
##   observed-data likelihood for the current lambda_tilde and r.
##
## This directly addresses the diagnostic finding that free kappa was
## absorbing baseline mean structure and pulling beta toward zero.
## ======================================================================

flatten_region_major <- function(mat) as.vector(mat)

if (!exists("ess_step", mode = "function")) {
    ess_step <- function(current, prior_sample, log_lik_fn, prior_mean) {
        nu <- prior_mean + prior_sample
        cur_ll <- log_lik_fn(current)
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
                return(list(sample = proposal, log_lik = prop_ll, n_reject = n_reject))
            }
            n_reject <- n_reject + 1L
            if (theta < 0) theta_min <- theta else theta_max <- theta
            if (theta_max - theta_min < 1e-12) {
                return(list(sample = current, log_lik = cur_ll, n_reject = n_reject))
            }
            theta <- runif(1, theta_min, theta_max)
        }
    }
}

fit_glm_for_ess <- function(dat,
                            min_sd = c(1.0, 0.25, 0.25),
                            inflate = 10) {
    y_vec  <- flatten_region_major(dat$y_coarse)
    x1_vec <- flatten_region_major(dat$x1)
    x2_vec <- flatten_region_major(dat$x2)
    offset_vec <- log(pmax(flatten_region_major(dat$e), 1e-12))
    glm_fit <- suppressWarnings(glm(y_vec ~ x1_vec + x2_vec,
                                    family = poisson(link = "log"),
                                    offset = offset_vec))
    beta_hat <- as.numeric(coef(glm_fit))
    beta_se  <- sqrt(diag(vcov(glm_fit)))
    if (length(beta_hat) != 3L || any(!is.finite(beta_hat)) || any(!is.finite(beta_se))) {
        beta_hat <- c(log(mean(y_vec + 0.1) / mean(flatten_region_major(dat$e))), 0, 0)
        beta_se  <- c(1, 1, 1)
    }
    list(center = beta_hat, scale = pmax(inflate * beta_se, min_sd))
}

build_beta_ess_base <- function(beta_current,
                                y_coarse, e, x1, x2,
                                lambda_tilde, phi,
                                fallback_base = NULL,
                                min_sd = c(1.0, 0.25, 0.25),
                                inflate = 10) {
    TT <- nrow(y_coarse)
    y_vec  <- flatten_region_major(y_coarse)
    x1_vec <- flatten_region_major(x1)
    x2_vec <- flatten_region_major(x2)
    phi_vec <- rep(phi, each = TT)
    offset_core <- pmax(e * lambda_tilde, 1e-12)
    offset_vec  <- log(flatten_region_major(offset_core)) + phi_vec
    glm_fit <- try(suppressWarnings(glm(y_vec ~ x1_vec + x2_vec,
                                        family = poisson(link = "log"),
                                        offset = offset_vec,
                                        control = glm.control(maxit = 50))), silent = TRUE)
    ok <- !inherits(glm_fit, "try-error")
    if (ok) {
        beta_hat <- as.numeric(coef(glm_fit))
        V <- try(vcov(glm_fit), silent = TRUE)
        ok <- length(beta_hat) == 3L && all(is.finite(beta_hat)) &&
            !inherits(V, "try-error")
        if (ok) {
            beta_se <- sqrt(diag(V))
            ok <- length(beta_se) == 3L && all(is.finite(beta_se))
        }
    }
    if (!ok) {
        if (!is.null(fallback_base) && length(fallback_base$center) == 3L && length(fallback_base$scale) == 3L) {
            beta_hat <- as.numeric(fallback_base$center)
            beta_se  <- as.numeric(fallback_base$scale)
        } else {
            beta_hat <- as.numeric(beta_current)
            beta_se  <- c(1, 0.5, 0.5)
        }
    }
    list(center = beta_hat, scale = pmax(inflate * beta_se, min_sd))
}

loglik_beta_marginal <- function(beta_vec, y_coarse, e, x1, x2,
                                 lambda_tilde, phi, r) {
    beta0 <- beta_vec[1]
    beta1 <- beta_vec[2]
    beta2 <- beta_vec[3]
    n1 <- ncol(y_coarse)
    ll <- 0
    for (j in seq_len(n1)) {
        eta_j <- beta0 + beta1 * x1[, j] + beta2 * x2[, j] + phi[j]
        mu_j  <- pmax(e[, j] * lambda_tilde[, j] * exp(eta_j), 1e-12)
        y_j   <- y_coarse[, j]
        if (is.infinite(r[j])) {
            ll <- ll + sum(dpois(y_j, lambda = mu_j, log = TRUE))
        } else {
            rj <- r[j]
            ll <- ll + sum(lgamma(y_j + rj) - lgamma(rj) - lgamma(y_j + 1) +
                           rj * log(rj) + y_j * log(mu_j) - (y_j + rj) * log(rj + mu_j))
        }
    }
    ll
}

update_beta <- function(beta_current, y_coarse, e, x1, x2,
                        kappa = NULL, lambda_tilde, phi, priors,
                        ess_base = NULL, r = NULL,
                        proposal_sd = NULL,
                        use_preconditioned = TRUE,
                        ess_min_sd = c(1.0, 0.25, 0.25),
                        ess_inflate = 10,
                        ...) {
    d <- length(beta_current)
    stopifnot(d == 3L)
    if (is.null(r)) stop("Collapsed beta update requires current r.")

    prior_mean <- c(priors$beta0_mean, priors$beta_mean)
    prior_sd   <- c(priors$beta0_sd,   priors$beta_sd)

    if (use_preconditioned) {
        base <- build_beta_ess_base(beta_current = beta_current,
                                    y_coarse = y_coarse, e = e,
                                    x1 = x1, x2 = x2,
                                    lambda_tilde = lambda_tilde,
                                    phi = phi,
                                    fallback_base = ess_base,
                                    min_sd = ess_min_sd,
                                    inflate = ess_inflate)
        ess_center <- base$center
        ess_scale  <- base$scale
        prior_sample <- rnorm(d, mean = 0, sd = ess_scale)
        log_lik_fn <- function(beta_vec) {
            ll <- loglik_beta_marginal(beta_vec, y_coarse, e, x1, x2,
                                       lambda_tilde, phi, r)
            log_true_prior <- -0.5 * sum(((beta_vec - prior_mean) / prior_sd)^2)
            log_ess_base   <- -0.5 * sum(((beta_vec - ess_center) / ess_scale)^2)
            ll + log_true_prior - log_ess_base
        }
        out <- ess_step(current = beta_current,
                        prior_sample = prior_sample,
                        log_lik_fn = log_lik_fn,
                        prior_mean = ess_center)
        out$ess_center <- ess_center
        out$ess_scale  <- ess_scale
        out$ess_mode   <- "preconditioned_collapsed_kappa"
        return(out)
    }

    prior_sample <- rnorm(d, mean = 0, sd = prior_sd)
    log_lik_fn <- function(beta_vec) {
        loglik_beta_marginal(beta_vec, y_coarse, e, x1, x2,
                             lambda_tilde, phi, r)
    }
    out <- ess_step(current = beta_current,
                    prior_sample = prior_sample,
                    log_lik_fn = log_lik_fn,
                    prior_mean = prior_mean)
    out$ess_center <- prior_mean
    out$ess_scale  <- prior_sd
    out$ess_mode   <- "ordinary_collapsed_kappa"
    out
}
