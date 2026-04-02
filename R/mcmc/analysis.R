## =========================================================
## simple_posterior_inference.R
## For current M1/M2 sampler:
##   beta0, phi, lambda_tilde  -> identified truth
##   beta, gamma, delta, r     -> raw truth
## =========================================================

## ---------- build truth map ----------
build_truth_map <- function(data_obj) {
    list(
        beta0 = data_obj$beta0_star_ident,
        beta = data_obj$beta_star,
        phi = data_obj$phi_star_ident,
        gamma = data_obj$gamma_star,
        delta = data_obj$delta_star,
        r = data_obj$r_star,
        lambda_tilde = data_obj$lambda_tilde_ident,
        kappa = data_obj$kappa,
        omega = data_obj$omega
    )
}

## ---------- helper: convert sample object to matrix ----------
## Assumes iteration is always the first dimension
as_sample_matrix <- function(x) {
    if (is.null(x)) stop("Input sample object is NULL.")

    if (is.atomic(x) && is.null(dim(x))) {
        out <- matrix(x, ncol = 1)
        colnames(out) <- "1"
        return(out)
    }

    if (is.matrix(x)) {
        out <- x
        if (is.null(colnames(out))) {
            colnames(out) <- paste0("V", seq_len(ncol(out)))
        }
        return(out)
    }

    if (length(dim(x)) >= 3) {
        d <- dim(x)
        n_iter <- d[1]
        out <- matrix(x, nrow = n_iter)

        idx <- arrayInd(seq_len(prod(d[-1])), .dim = d[-1])
        colnames(out) <- apply(idx, 1, function(ii) {
            paste0("(", paste(ii, collapse = ","), ")")
        })
        return(out)
    }

    stop("Unsupported sample structure.")
}

## ---------- helper: align truth with sample columns ----------
get_truth_vector <- function(sample_obj, truth_obj) {
    mat <- as_sample_matrix(sample_obj)
    n_comp <- ncol(mat)

    if (is.null(truth_obj)) {
        return(rep(NA_real_, n_comp))
    }

    if (is.atomic(truth_obj) && is.null(dim(truth_obj))) {
        return(rep(unname(as.numeric(truth_obj)), n_comp))
    }

    truth_vec <- unname(as.numeric(truth_obj))

    if (length(truth_vec) == 1L) {
        return(rep(truth_vec, n_comp))
    }

    if (length(truth_vec) != n_comp) {
        warning("Truth length does not match number of sample components. Returning NA.")
        return(rep(NA_real_, n_comp))
    }

    truth_vec
}
## ---------- summary for one parameter ----------
summarize_param <- function(fit, truth_map, param) {
    x <- fit$samples[[param]]
    if (is.null(x)) stop("Parameter not found in fit$samples: ", param)

    mat <- as_sample_matrix(x)
    truth_vec <- get_truth_vector(x, truth_map[[param]])

    n_comp <- ncol(mat)

    out <- data.frame(
        parameter = rep(param, n_comp),
        component = unname(colnames(mat)),
        truth = unname(truth_vec),
        mean = unname(apply(mat, 2, mean, na.rm = TRUE)),
        sd = unname(apply(mat, 2, sd, na.rm = TRUE)),
        q2.5 = unname(apply(mat, 2, quantile, probs = 0.025, na.rm = TRUE)),
        median = unname(apply(mat, 2, quantile, probs = 0.50, na.rm = TRUE)),
        q97.5 = unname(apply(mat, 2, quantile, probs = 0.975, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )

    out$bias <- out$mean - out$truth
    out$covered <- with(out, truth >= q2.5 & truth <= q97.5)
    out$covered[is.na(out$truth)] <- NA

    rownames(out) <- NULL
    out
}
## ---------- summary for several parameters ----------
posterior_summary <- function(fit,
                              data_obj,
                              params = c("beta0", "beta", "phi", "gamma", "delta", "r"),
                              max_components = 200L) {
    truth_map <- build_truth_map(data_obj)
    out_list <- list()

    for (param in params) {
        x <- fit$samples[[param]]
        if (is.null(x)) next

        dims <- dim(x)
        n_comp <- if (is.null(dims)) 1L else prod(dims[-1])

        if (n_comp > max_components) {
            message("Skipping ", param, " because it has ", n_comp, " components.")
            next
        }

        out_list[[param]] <- summarize_param(fit, truth_map, param)
    }

    out <- do.call(rbind, out_list)
    rownames(out) <- NULL
    out
}

## ---------- trace + density ----------
plot_trace_density <- function(fit, data_obj, param, component = 1L) {
    truth_map <- build_truth_map(data_obj)

    x <- fit$samples[[param]]
    if (is.null(x)) stop("Parameter not found in fit$samples: ", param)

    mat <- as_sample_matrix(x)
    if (component < 1L || component > ncol(mat)) stop("component out of range.")

    draw_vec <- mat[, component]
    truth_val <- get_truth_vector(x, truth_map[[param]])[component]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(1, 2))

    plot(draw_vec, type = "l",
         xlab = "Iteration", ylab = "Value",
         main = paste0("Trace: ", param, "[", colnames(mat)[component], "]"))
    if (!is.na(truth_val)) {
        abline(h = truth_val, lty = 2, lwd = 2)
        legend("topright", legend = "Truth", lty = 2, lwd = 2, bty = "n")
    }

    dens <- density(draw_vec, na.rm = TRUE)
    plot(dens,
         xlab = "Value", ylab = "Density",
         main = paste0("Density: ", param, "[", colnames(mat)[component], "]"))
    abline(v = mean(draw_vec, na.rm = TRUE), lty = 3, lwd = 2)
    if (!is.na(truth_val)) {
        abline(v = truth_val, lty = 2, lwd = 2)
        legend("topright",
               legend = c("Truth", "Posterior mean"),
               lty = c(2, 3), lwd = 2, bty = "n")
    } else {
        legend("topright",
               legend = "Posterior mean",
               lty = 3, lwd = 2, bty = "n")
    }
}

## ---------- posterior path vs truth ----------
## Supports:
##   sample = iter x T
##   sample = iter x T x J
plot_posterior_path <- function(fit, data_obj, param, index2 = NULL) {
    truth_map <- build_truth_map(data_obj)

    x <- fit$samples[[param]]
    if (is.null(x)) stop("Parameter not found in fit$samples: ", param)

    truth_x <- truth_map[[param]]
    d <- dim(x)
    if (is.null(d) || length(d) < 2) stop("Sample must be at least iter x T.")

    ## case: iter x T
    if (length(d) == 2) {
        post_mean <- apply(x, 2, mean, na.rm = TRUE)
        post_l <- apply(x, 2, quantile, probs = 0.025, na.rm = TRUE)
        post_u <- apply(x, 2, quantile, probs = 0.975, na.rm = TRUE)

        tt <- seq_along(post_mean)
        ylim_all <- range(c(post_l, post_u, truth_x), na.rm = TRUE)

        plot(tt, post_mean, type = "l", lwd = 2,
             xlab = "Time", ylab = param,
             ylim = ylim_all,
             main = paste0("Posterior path: ", param))
        lines(tt, post_l, lty = 3)
        lines(tt, post_u, lty = 3)

        if (!is.null(truth_x) && length(truth_x) == length(tt)) {
            lines(tt, truth_x, lty = 2, lwd = 2)
            legend("topright",
                   legend = c("Posterior mean", "95% CI", "Truth"),
                   lty = c(1, 3, 2), lwd = c(2, 1, 2), bty = "n")
        }
        return(invisible(NULL))
    }

    ## case: iter x T x J
    if (length(d) == 3) {
        if (is.null(index2)) index2 <- 1L
        if (index2 < 1L || index2 > d[3]) stop("index2 out of range.")

        post_mean <- apply(x[, , index2, drop = FALSE], 2, mean, na.rm = TRUE)
        post_l <- apply(x[, , index2, drop = FALSE], 2, quantile, probs = 0.025, na.rm = TRUE)
        post_u <- apply(x[, , index2, drop = FALSE], 2, quantile, probs = 0.975, na.rm = TRUE)

        post_mean <- as.numeric(post_mean)
        post_l <- as.numeric(post_l)
        post_u <- as.numeric(post_u)

        truth_line <- NULL
        if (!is.null(truth_x) && length(dim(truth_x)) == 2 && ncol(truth_x) >= index2) {
            truth_line <- truth_x[, index2]
        }

        tt <- seq_along(post_mean)
        ylim_all <- range(c(post_l, post_u, truth_line), na.rm = TRUE)

        plot(tt, post_mean, type = "l", lwd = 2,
             xlab = "Time", ylab = paste0(param, "[, ", index2, "]"),
             ylim = ylim_all,
             main = paste0("Posterior path: ", param, " region ", index2))
        lines(tt, post_l, lty = 3)
        lines(tt, post_u, lty = 3)

        if (!is.null(truth_line) && length(truth_line) == length(tt)) {
            lines(tt, truth_line, lty = 2, lwd = 2)
            legend("topright",
                   legend = c("Posterior mean", "95% CI", "Truth"),
                   lty = c(1, 3, 2), lwd = c(2, 1, 2), bty = "n")
        }
        return(invisible(NULL))
    }

    stop("Currently only supports iter x T or iter x T x J.")
}

## ---------- one simple driver ----------
simple_posterior_analysis <- function(fit,
                                      data_obj,
                                      out_dir = "output/posterior_analysis",
                                      lambda_regions = 1:3) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    ## 1. Main table
    tab_main <- posterior_summary(
        fit = fit,
        data_obj = data_obj,
        params = c("beta0", "beta", "phi", "gamma", "delta", "r")
    )
    write.csv(tab_main,
              file.path(out_dir, "posterior_summary_main.csv"),
              row.names = FALSE)

    ## 2. A few scalar/vector plots
    pdf(file.path(out_dir, "trace_density_main.pdf"), width = 10, height = 4)

    if ("beta0" %in% names(fit$samples)) plot_trace_density(fit, data_obj, "beta0", 1)
    if ("beta" %in% names(fit$samples)) {
        plot_trace_density(fit, data_obj, "beta", 1)
        if (ncol(fit$samples$beta) >= 2) plot_trace_density(fit, data_obj, "beta", 2)
    }
    if ("phi" %in% names(fit$samples)) plot_trace_density(fit, data_obj, "phi", 1)
    if ("gamma" %in% names(fit$samples)) plot_trace_density(fit, data_obj, "gamma", 1)
    if ("r" %in% names(fit$samples)) plot_trace_density(fit, data_obj, "r", 1)
    if ("delta" %in% names(fit$samples)) plot_trace_density(fit, data_obj, "delta", 1)

    dev.off()

    ## 3. Lambda path plots for selected regions
    if ("lambda_tilde" %in% names(fit$samples)) {
        pdf(file.path(out_dir, "lambda_paths.pdf"), width = 7, height = 5)
        for (j in lambda_regions) {
            plot_posterior_path(fit, data_obj, "lambda_tilde", index2 = j)
        }
        dev.off()
    }

    invisible(tab_main)
}

## =========================================================
## Simple MCMC plots only
## =========================================================

as_sample_matrix <- function(x) {
    if (is.null(x)) stop("Input sample object is NULL.")

    if (is.atomic(x) && is.null(dim(x))) {
        out <- matrix(x, ncol = 1)
        colnames(out) <- "1"
        return(out)
    }

    if (is.matrix(x)) {
        out <- x
        if (is.null(colnames(out))) {
            colnames(out) <- paste0("V", seq_len(ncol(out)))
        }
        return(out)
    }

    if (length(dim(x)) >= 3) {
        d <- dim(x)
        n_iter <- d[1]
        out <- matrix(x, nrow = n_iter)

        idx <- arrayInd(seq_len(prod(d[-1])), .dim = d[-1])
        colnames(out) <- apply(idx, 1, function(ii) {
            paste0("(", paste(ii, collapse = ","), ")")
        })
        return(out)
    }

    stop("Unsupported sample structure.")
}

plot_mcmc_one <- function(fit, param, component = 1L, breaks = 30) {
    x <- fit$samples[[param]]
    if (is.null(x)) stop("Parameter not found in fit$samples: ", param)

    mat <- as_sample_matrix(x)
    if (component < 1L || component > ncol(mat)) stop("component out of range.")

    draw_vec <- mat[, component]
    comp_name <- colnames(mat)[component]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mfrow = c(1, 3))

    plot(draw_vec, type = "l",
         xlab = "Iteration", ylab = "Value",
         main = paste0("Trace: ", param, "[", comp_name, "]"))

    hist(draw_vec, breaks = breaks,
         main = paste0("Histogram: ", param, "[", comp_name, "]"),
         xlab = "Value", border = "white")

    dens <- density(draw_vec, na.rm = TRUE)
    plot(dens,
         main = paste0("Density: ", param, "[", comp_name, "]"),
         xlab = "Value", ylab = "Density")
    abline(v = mean(draw_vec, na.rm = TRUE), lty = 2, lwd = 2)
}

plot_mcmc_batch <- function(fit, out_pdf = NULL) {
    if (!is.null(out_pdf)) {
        pdf(out_pdf, width = 12, height = 4)
        on.exit(dev.off(), add = TRUE)
    }

    if ("beta0" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "beta0", 1)
    }

    if ("beta" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "beta", 1)
        if (ncol(fit$samples$beta) >= 2) {
            plot_mcmc_one(fit, "beta", 2)
        }
    }

    if ("phi" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "phi", 1)
    }

    if ("gamma" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "gamma", 1)
    }

    if ("r" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "r", 1)
    }

    if ("delta" %in% names(fit$samples)) {
        plot_mcmc_one(fit, "delta", 1)
    }
}


## =========================================================
## Trace + density plots for gamma and r
## =========================================================

plot_trace_density_param <- function(samples_mat, truth_vec = NULL,
                                     idx = 1L, param_name = "param") {
    draw <- samples_mat[, idx]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(1, 2))

    ## trace
    plot(draw, type = "l",
         xlab = "Iteration", ylab = "Value",
         main = paste0("Trace: ", param_name, idx))
    if (!is.null(truth_vec)) {
        abline(h = truth_vec[idx], lty = 2, lwd = 2)
    }

    ## density
    dens <- density(draw, na.rm = TRUE)
    plot(dens,
         xlab = "Value", ylab = "Density",
         main = paste0("Density: ", param_name, idx))
    abline(v = mean(draw, na.rm = TRUE), lty = 3, lwd = 2)
    if (!is.null(truth_vec)) {
        abline(v = truth_vec[idx], lty = 2, lwd = 2)
        legend("topright",
               legend = c("Truth", "Posterior mean"),
               lty = c(2, 3), lwd = 2, bty = "n")
    } else {
        legend("topright",
               legend = "Posterior mean",
               lty = 3, lwd = 2, bty = "n")
    }
}

plot_gamma_r_key <- function(fit, data_obj,
                             gamma_idx = c(4, 6, 9),
                             r_idx = c(4, 3, 2),
                             out_pdf = NULL) {

    gamma_draws <- fit$samples$gamma
    r_draws <- fit$samples$r

    gamma_truth <- as.numeric(data_obj$gamma_star)
    if (length(gamma_truth) == 1L) {
        gamma_truth <- rep(gamma_truth, ncol(gamma_draws))
    }

    r_truth <- as.numeric(data_obj$r_star)

    if (!is.null(out_pdf)) {
        pdf(out_pdf, width = 10, height = 4)
        on.exit(dev.off(), add = TRUE)
    }

    ## gamma
    for (j in gamma_idx) {
        plot_trace_density_param(
            samples_mat = gamma_draws,
            truth_vec = gamma_truth,
            idx = j,
            param_name = "gamma"
        )
    }

    ## r
    for (j in r_idx) {
        plot_trace_density_param(
            samples_mat = r_draws,
            truth_vec = r_truth,
            idx = j,
            param_name = "r"
        )
    }
}