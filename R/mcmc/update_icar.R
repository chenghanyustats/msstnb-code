## ======================================================================
## update_icar.R
## Phi update using a kappa-collapsed marginal NB likelihood.
##
## Key design change:
##   phi is NO LONGER updated conditional on sampled kappa.
##   Instead, kappa is integrated out analytically and phi targets the
##   observed-data likelihood for the current lambda_tilde and r.
## ======================================================================

log_posterior_phi <- function(u_vec, B, BHB, tau_phi,
                              y_coarse, e, x1, x2,
                              beta0, beta, lambda_tilde, r,
                              kappa = NULL) {
    phi <- as.numeric(B %*% u_vec)
    log_prior <- -0.5 * tau_phi * as.numeric(t(u_vec) %*% BHB %*% u_vec)
    n1 <- ncol(y_coarse)
    ll <- 0
    for (j in seq_len(n1)) {
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
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
    log_prior + ll
}

update_phi <- function(u_current, B, BHB, tau_phi,
                       y_coarse, e, x1, x2,
                       beta0, beta, kappa = NULL, lambda_tilde,
                       r,
                       proposal_sd) {
    d <- length(u_current)
    u_proposal <- u_current + rnorm(d, mean = 0, sd = proposal_sd)
    lp_current <- log_posterior_phi(u_current, B, BHB, tau_phi,
                                    y_coarse, e, x1, x2,
                                    beta0, beta, lambda_tilde, r)
    lp_prop <- log_posterior_phi(u_proposal, B, BHB, tau_phi,
                                 y_coarse, e, x1, x2,
                                 beta0, beta, lambda_tilde, r)
    log_alpha <- lp_prop - lp_current
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        return(list(u = u_proposal,
                    phi = as.numeric(B %*% u_proposal),
                    accept = TRUE))
    }
    list(u = u_current,
         phi = as.numeric(B %*% u_current),
         accept = FALSE)
}

update_tau_phi <- function(phi, H, n1, priors) {
    shape_post <- priors$tau_phi_shape + (n1 - 1) / 2
    rate_post  <- priors$tau_phi_rate + as.numeric(t(phi) %*% H %*% phi) / 2
    rgamma(1, shape = shape_post, rate = rate_post)
}
