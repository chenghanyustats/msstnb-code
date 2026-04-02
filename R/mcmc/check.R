## =========================================================
## Posterior mean / SD / CI / bias summary against truth
## For current M1 sampler
## =========================================================

posterior_summary_against_truth <- function(fit, data_obj) {
    s <- fit$samples

    out <- list()

    ## ----- beta0: compare to identified truth -----
    out$beta0 <- data.frame(
        parameter = "beta0",
        truth = data_obj$beta0_star_ident,
        mean = mean(s$beta0),
        sd = sd(s$beta0),
        q2.5 = quantile(s$beta0, 0.025),
        median = quantile(s$beta0, 0.5),
        q97.5 = quantile(s$beta0, 0.975)
    )
    out$beta0$bias <- out$beta0$mean - out$beta0$truth
    out$beta0$covered <- with(out$beta0, truth >= q2.5 & truth <= q97.5)

    ## ----- beta: compare to beta_star -----
    beta_tab <- data.frame(
        parameter = c("beta1", "beta2"),
        truth = as.numeric(data_obj$beta_star),
        mean = apply(s$beta, 2, mean),
        sd = apply(s$beta, 2, sd),
        q2.5 = apply(s$beta, 2, quantile, probs = 0.025),
        median = apply(s$beta, 2, quantile, probs = 0.5),
        q97.5 = apply(s$beta, 2, quantile, probs = 0.975)
    )
    beta_tab$bias <- beta_tab$mean - beta_tab$truth
    beta_tab$covered <- with(beta_tab, truth >= q2.5 & truth <= q97.5)
    out$beta <- beta_tab

    ## ----- gamma: one row per region -----
    gamma_truth <- as.numeric(data_obj$gamma_star)
    if (length(gamma_truth) == 1L) {
        gamma_truth <- rep(gamma_truth, ncol(s$gamma))
    }
    gamma_tab <- data.frame(
        parameter = paste0("gamma", seq_len(ncol(s$gamma))),
        truth = gamma_truth,
        mean = apply(s$gamma, 2, mean),
        sd = apply(s$gamma, 2, sd),
        q2.5 = apply(s$gamma, 2, quantile, probs = 0.025),
        median = apply(s$gamma, 2, quantile, probs = 0.5),
        q97.5 = apply(s$gamma, 2, quantile, probs = 0.975)
    )
    gamma_tab$bias <- gamma_tab$mean - gamma_tab$truth
    gamma_tab$covered <- with(gamma_tab, truth >= q2.5 & truth <= q97.5)
    out$gamma <- gamma_tab

    ## ----- r: one row per region -----
    r_tab <- data.frame(
        parameter = paste0("r", seq_len(ncol(s$r))),
        truth = as.numeric(data_obj$r_star),
        mean = apply(s$r, 2, mean),
        sd = apply(s$r, 2, sd),
        q2.5 = apply(s$r, 2, quantile, probs = 0.025),
        median = apply(s$r, 2, quantile, probs = 0.5),
        q97.5 = apply(s$r, 2, quantile, probs = 0.975)
    )
    r_tab$bias <- r_tab$mean - r_tab$truth
    r_tab$covered <- with(r_tab, truth >= q2.5 & truth <= q97.5)
    out$r <- r_tab

    ## ----- delta: scalar -----
    out$delta <- data.frame(
        parameter = "delta",
        truth = data_obj$delta_star,
        mean = mean(s$delta),
        sd = sd(s$delta),
        q2.5 = quantile(s$delta, 0.025),
        median = quantile(s$delta, 0.5),
        q97.5 = quantile(s$delta, 0.975)
    )
    out$delta$bias <- out$delta$mean - out$delta$truth
    out$delta$covered <- with(out$delta, truth >= q2.5 & truth <= q97.5)

    ## ----- phi: compare to identified truth -----
    phi_tab <- data.frame(
        parameter = paste0("phi", seq_len(ncol(s$phi))),
        truth = as.numeric(data_obj$phi_star_ident),
        mean = apply(s$phi, 2, mean),
        sd = apply(s$phi, 2, sd),
        q2.5 = apply(s$phi, 2, quantile, probs = 0.025),
        median = apply(s$phi, 2, quantile, probs = 0.5),
        q97.5 = apply(s$phi, 2, quantile, probs = 0.975)
    )
    phi_tab$bias <- phi_tab$mean - phi_tab$truth
    phi_tab$covered <- with(phi_tab, truth >= q2.5 & truth <= q97.5)
    out$phi <- phi_tab

    ## combine
    out$all <- do.call(rbind, out[c("beta0", "beta", "gamma", "r", "delta", "phi")])
    rownames(out$all) <- NULL
    out
}

lambda_region_summary <- function(fit, data_obj, j = 1) {
    post_mean <- apply(fit$samples$lambda_tilde[, , j, drop = FALSE], 2, mean)
    post_l <- apply(fit$samples$lambda_tilde[, , j, drop = FALSE], 2, quantile, probs = 0.025)
    post_u <- apply(fit$samples$lambda_tilde[, , j, drop = FALSE], 2, quantile, probs = 0.975)
    truth <- data_obj$lambda_tilde_ident[, j]

    data.frame(
        time = seq_along(truth),
        truth = as.numeric(truth),
        mean = as.numeric(post_mean),
        q2.5 = as.numeric(post_l),
        q97.5 = as.numeric(post_u),
        bias = as.numeric(post_mean) - as.numeric(truth),
        covered = truth >= post_l & truth <= post_u
    )
}

short_table <- function(res) {
    res$all[, c("parameter", "truth", "mean", "bias", "sd", "covered")]
}
