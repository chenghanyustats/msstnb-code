## ======================================================================
## update_icar_revised.R
## ICAR spatial random-effect update using a kappa-collapsed marginal
## negative-binomial likelihood.
##
## Main design:
##   phi is intentionally NOT updated conditional on sampled kappa.
##   Instead, kappa is integrated out analytically, so phi targets the
##   observed-data marginal likelihood conditional on beta, lambda_tilde, and r.
##
## Constraint handling:
##   phi = B u, where columns of B span the sum-to-zero ICAR subspace.
## ======================================================================

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

.check_phi_inputs <- function(u_vec, B, BHB, tau_phi,
                              y_coarse, e, x1, x2,
                              beta0, beta, lambda_tilde, r,
                              proposal_sd = NULL) {
    if (!is.matrix(y_coarse)) stop("y_coarse must be a matrix.", call. = FALSE)
    if (any(!is.finite(y_coarse)) || any(y_coarse < 0)) {
        stop("y_coarse must contain nonnegative finite counts.", call. = FALSE)
    }
    .check_same_dim(y_coarse, e, x1, x2, lambda_tilde,
                    names = c("e", "x1", "x2", "lambda_tilde"))
    if (any(!is.finite(e)) || any(e <= 0)) stop("e must be positive and finite.", call. = FALSE)
    if (any(!is.finite(lambda_tilde)) || any(lambda_tilde <= 0)) {
        stop("lambda_tilde must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(x1)) || any(!is.finite(x2))) {
        stop("x1 and x2 must contain only finite values.", call. = FALSE)
    }

    n1 <- ncol(y_coarse)
    d <- length(u_vec)
    if (!is.matrix(B) || nrow(B) != n1 || ncol(B) != d) {
        stop("B must be an n1 by length(u_vec) matrix.", call. = FALSE)
    }
    if (!is.matrix(BHB) || !identical(dim(BHB), c(d, d))) {
        stop("BHB must be a length(u_vec) by length(u_vec) matrix.", call. = FALSE)
    }
    if (any(!is.finite(B)) || any(!is.finite(BHB))) {
        stop("B and BHB must contain only finite values.", call. = FALSE)
    }
    if (!is.finite(tau_phi) || tau_phi <= 0) {
        stop("tau_phi must be positive and finite.", call. = FALSE)
    }
    if (!is.finite(beta0) || length(beta) != 2L || any(!is.finite(beta))) {
        stop("beta0 must be finite and beta must be a finite length-2 vector.", call. = FALSE)
    }
    if (length(r) != n1) {
        stop("r must have length ncol(y_coarse).", call. = FALSE)
    }
    if (any((r <= 0 | !is.finite(r)) & !is.infinite(r))) {
        stop("Each r[j] must be positive finite or Inf.", call. = FALSE)
    }
    if (!is.null(proposal_sd)) {
        if (length(proposal_sd) == 1L) {
            if (!is.finite(proposal_sd) || proposal_sd <= 0) {
                stop("proposal_sd must be positive and finite.", call. = FALSE)
            }
        } else if (length(proposal_sd) == d) {
            if (any(!is.finite(proposal_sd)) || any(proposal_sd <= 0)) {
                stop("proposal_sd must be positive and finite.", call. = FALSE)
            }
        } else {
            stop("proposal_sd must be length 1 or length(u_current).", call. = FALSE)
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

log_posterior_phi <- function(u_vec, B, BHB, tau_phi,
                              y_coarse, e, x1, x2,
                              beta0, beta, lambda_tilde, r,
                              kappa = NULL) {
    ## kappa is intentionally ignored.  The phi update uses the
    ## kappa-collapsed marginal NB likelihood.
    invisible(kappa)

    if (any(!is.finite(u_vec))) return(-Inf)

    phi <- as.numeric(B %*% u_vec)
    quad <- as.numeric(crossprod(u_vec, BHB %*% u_vec))
    if (!is.finite(quad) || quad < -1e-8) return(-Inf)
    quad <- max(quad, 0)
    log_prior <- -0.5 * tau_phi * quad

    n1 <- ncol(y_coarse)
    ll <- 0
    for (j in seq_len(n1)) {
        log_mu_j <- log(e[, j]) + log(lambda_tilde[, j]) +
            beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        ll_j <- .nb_or_poisson_loglik(y = y_coarse[, j], log_mu = log_mu_j, r = r[j])
        if (!is.finite(ll_j)) return(-Inf)
        ll <- ll + ll_j
    }

    log_prior + ll
}

update_phi <- function(u_current, B, BHB, tau_phi,
                       y_coarse, e, x1, x2,
                       beta0, beta, kappa = NULL, lambda_tilde,
                       r,
                       proposal_sd) {
    ## kappa is intentionally ignored.  The phi update uses the
    ## kappa-collapsed marginal NB likelihood.
    invisible(kappa)

    .check_phi_inputs(u_current, B, BHB, tau_phi,
                      y_coarse, e, x1, x2,
                      beta0, beta, lambda_tilde, r,
                      proposal_sd = proposal_sd)

    d <- length(u_current)
    u_proposal <- u_current + rnorm(d, mean = 0, sd = proposal_sd)

    lp_current <- log_posterior_phi(u_current, B, BHB, tau_phi,
                                    y_coarse, e, x1, x2,
                                    beta0, beta, lambda_tilde, r)
    lp_proposal <- log_posterior_phi(u_proposal, B, BHB, tau_phi,
                                     y_coarse, e, x1, x2,
                                     beta0, beta, lambda_tilde, r)

    log_alpha <- lp_proposal - lp_current

    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        phi_proposal <- as.numeric(B %*% u_proposal)
        return(list(u = u_proposal,
                    phi = phi_proposal,
                    accept = TRUE,
                    log_posterior = lp_proposal,
                    log_alpha = log_alpha))
    }

    phi_current <- as.numeric(B %*% u_current)
    list(u = u_current,
         phi = phi_current,
         accept = FALSE,
         log_posterior = lp_current,
         log_alpha = log_alpha)
}

update_tau_phi <- function(phi, H, n1, priors) {
    if (length(phi) != n1 || any(!is.finite(phi))) {
        stop("phi must be a finite vector of length n1.", call. = FALSE)
    }
    if (!is.matrix(H) || !identical(dim(H), c(n1, n1)) || any(!is.finite(H))) {
        stop("H must be a finite n1 by n1 matrix.", call. = FALSE)
    }
    if (!is.finite(priors$tau_phi_shape) || priors$tau_phi_shape <= 0 ||
        !is.finite(priors$tau_phi_rate) || priors$tau_phi_rate <= 0) {
        stop("tau_phi prior shape and rate must be positive and finite.", call. = FALSE)
    }

    roughness <- as.numeric(crossprod(phi, H %*% phi))
    if (!is.finite(roughness)) stop("phi' H phi is not finite.", call. = FALSE)
    roughness <- max(roughness, 0)

    shape_post <- priors$tau_phi_shape + (n1 - 1) / 2
    rate_post <- priors$tau_phi_rate + roughness / 2

    rgamma(1, shape = shape_post, rate = rate_post)
}
