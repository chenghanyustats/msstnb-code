`%||%` <- function(x, y) if (is.null(x)) y else x

as_numeric_scalar <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) vals <- unlist(x, recursive = TRUE, use.names = FALSE) else vals <- x
    vals <- as.numeric(vals)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return(NULL)
    uniq <- unique(round(vals, 12))
    if (length(uniq) == 1L) return(as.numeric(uniq))
    as.numeric(vals[1L])
}

extract_truth_gamma_compat <- function(ds) {
    out <- as_numeric_scalar(ds$gamma_star %||% NULL)
    if (!is.null(out)) return(out)
    truth <- ds$truth %||% list()
    cand <- truth$gamma %||% truth$gamma_common %||% truth$gamma_lambda %||% NULL
    out <- as_numeric_scalar(cand)
    if (!is.null(out)) return(out)
    sr <- ds$scenario_row %||% ds$scenario %||% list()
    as_numeric_scalar(sr$gamma_lambda %||% sr$gamma %||% NULL)
}

extract_y_coarse_compat <- function(ds) {
    if (!is.null(ds$y_coarse)) return(as.matrix(ds$y_coarse))
    if (!is.null(ds$data) && !is.null(ds$data$y_levels) && length(ds$data$y_levels) >= 1L) {
        return(as.matrix(ds$data$y_levels[[1L]]))
    }
    stop("Could not find coarsest counts. Expected ds$y_coarse or ds$data$y_levels[[1]].")
}

extract_gamma_draws_compat <- function(fit) {
    if (!is.null(fit$samples$gamma)) {
        g <- fit$samples$gamma
        if (is.matrix(g)) return(as.numeric(g[, 1L]))
        return(as.numeric(g))
    }
    if (!is.null(fit$draws$gamma_common)) return(as.numeric(fit$draws$gamma_common))
    NULL
}

extract_gamma_grid_compat <- function(fit, M_default = 200L) {
    diag_obj <- fit$gamma_conditional_diag %||% list()
    grid <- diag_obj$grid %||% NULL
    if (!is.null(grid)) return(as.numeric(grid))
    seq(0.01, 0.999, length.out = M_default)
}

extract_prior_params_compat <- function(fit, a_default = 1, b_default = 1) {
    hyper <- fit$hyper %||% fit$meta$hyper %||% list()
    gamma_hyper <- hyper$gamma %||% list()
    a <- as_numeric_scalar(gamma_hyper$a_gamma %||% gamma_hyper$a %||% NULL) %||% a_default
    b <- as_numeric_scalar(gamma_hyper$b_gamma %||% gamma_hyper$b %||% NULL) %||% b_default
    c(a = a, b = b)
}

extract_lambda_raw_draws_compat <- function(fit, ds = NULL) {
    ## Highest priority: sampler-recorded raw path used for gamma diagnostics
    if (!is.null(fit$samples$lambda_raw_for_gamma)) {
        arr <- as.array(fit$samples$lambda_raw_for_gamma)
        if (length(dim(arr)) == 3L) return(arr)
    }

    draws <- fit$draws %||% list()
    cand_names <- c(
        "lambda_raw_for_gamma",
        "lambda_tilde_raw_for_gamma",
        "lambda_raw",
        "lambda_tilde_0T_raw",
        "lambda_tilde_0T"
    )
    for (nm in cand_names) {
        obj <- draws[[nm]]
        if (!is.null(obj)) {
            arr <- as.array(obj)
            if (length(dim(arr)) == 3L) return(arr)
        }
    }

    if (!is.null(fit$samples$lambda_tilde)) {
        arr <- as.array(fit$samples$lambda_tilde)
        if (length(dim(arr)) == 3L) {
            nd <- dim(arr)[1L]; TT <- dim(arr)[2L]; n1 <- dim(arr)[3L]
            out <- array(NA_real_, dim = c(nd, TT + 1L, n1))
            out[, 1L, ] <- 1
            out[, 2:(TT + 1L), ] <- arr
            warning("Using pseudo-raw lambda path with lambda_0 fixed at 1. ",
                    "For exact comparability, use a fit that stores samples$lambda_raw_for_gamma.")
            return(out)
        }
    }

    stop("Could not find lambda draws compatible with frozen-a diagnostic.")
}

discrete_weighted_quantile <- function(x, w, probs) {
    ord <- order(x); x <- x[ord]; w <- w[ord]
    cw <- cumsum(w) / sum(w)
    vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1L))
}

summarize_curve <- function(grid, weights) {
    weights <- weights / sum(weights)
    mode <- grid[which.max(weights)]
    meanv <- sum(grid * weights)
    qs <- discrete_weighted_quantile(grid, weights, c(0.025, 0.05, 0.5, 0.95, 0.975))
    list(mode = mode, mean = meanv, median = qs[3L], q05 = qs[2L], q95 = qs[4L], q025 = qs[1L], q975 = qs[5L])
}

compute_a_prev_fixed <- function(y, gamma_ref, a0j) {
    TT <- nrow(y); n1 <- ncol(y)
    a_prev <- matrix(NA_real_, nrow = TT, ncol = n1)
    a_cur <- rep_len(a0j, n1)
    for (t in seq_len(TT)) {
        a_prev[t, ] <- a_cur
        a_cur <- gamma_ref * a_cur + y[t, ]
    }
    a_prev
}

compute_gamma_curve_for_one_draw <- function(lambda_draw, y, grid, a0j, prior_ab, a_prev_fixed = NULL) {
    TT <- nrow(y); n1 <- ncol(y)
    if (nrow(lambda_draw) != TT + 1L || ncol(lambda_draw) != n1) {
        stop("lambda_draw must have dimension (TT+1) x n1")
    }
    logw <- rep(-Inf, length(grid))
    for (m in seq_along(grid)) {
        gm <- grid[m]
        lp <- dbeta(gm, prior_ab["a"], prior_ab["b"], log = TRUE)
        a_cur <- if (is.null(a_prev_fixed)) rep_len(a0j, n1) else NULL
        ok <- TRUE
        for (t in seq_len(TT)) {
            lam_prev <- lambda_draw[t, ]
            lam_cur  <- lambda_draw[t + 1L, ]
            u <- gm * lam_cur / lam_prev
            if (any(!is.finite(u)) || any(u <= 0 | u >= 1) || any(lam_prev <= 0)) {
                ok <- FALSE; break
            }
            a_prev_t <- if (is.null(a_prev_fixed)) a_cur else a_prev_fixed[t, ]
            alpha <- gm * a_prev_t
            beta  <- (1 - gm) * a_prev_t
            if (any(alpha <= 0) || any(beta <= 0)) {
                ok <- FALSE; break
            }
            lp <- lp + sum(dbeta(u, alpha, beta, log = TRUE) + log(gm) - log(lam_prev))
            if (is.null(a_prev_fixed)) a_cur <- gm * a_cur + y[t, ]
        }
        if (ok) logw[m] <- lp
    }
    if (all(!is.finite(logw))) stop("All grid points received -Inf log weight for at least one draw.")
    mx <- max(logw[is.finite(logw)])
    w <- exp(logw - mx); w[!is.finite(w)] <- 0; w <- w / sum(w)
    list(logw = logw, weights = w, curve = summarize_curve(grid, w))
}

run_frozen_a_gamma_diagnostic <- function(
    data_file,
    fit_file,
    out_prefix = "diag_frozen_a_gamma",
    gamma_ref = NULL,
    gamma_ref_mode = c("truth", "posterior_mean", "fixed"),
    lambda_a0 = 10,
    grid = NULL,
    prior_a = NULL,
    prior_b = NULL,
    save_png = TRUE,
    save_pdf = TRUE
) {
    gamma_ref_mode <- match.arg(gamma_ref_mode)
    ds  <- readRDS(data_file)
    fit <- readRDS(fit_file)

    y <- extract_y_coarse_compat(ds)
    lambda_arr <- extract_lambda_raw_draws_compat(fit, ds = ds)
    grid <- grid %||% extract_gamma_grid_compat(fit)
    prior_ab <- extract_prior_params_compat(fit)
    if (!is.null(prior_a)) prior_ab["a"] <- prior_a
    if (!is.null(prior_b)) prior_ab["b"] <- prior_b

    gamma_truth <- extract_truth_gamma_compat(ds)
    if (gamma_ref_mode == "truth") {
        gamma_ref <- gamma_truth
        if (is.null(gamma_ref)) stop("Truth gamma not found; supply gamma_ref manually.")
    } else if (gamma_ref_mode == "posterior_mean") {
        gamma_draws <- extract_gamma_draws_compat(fit)
        if (is.null(gamma_draws)) stop("Posterior gamma draws not found.")
        gamma_ref <- mean(gamma_draws)
    } else {
        if (is.null(gamma_ref)) stop("When gamma_ref_mode='fixed', you must provide gamma_ref.")
    }

    TT <- nrow(y); n1 <- ncol(y); nd <- dim(lambda_arr)[1L]
    if (dim(lambda_arr)[2L] != TT + 1L || dim(lambda_arr)[3L] != n1) stop("lambda draw array dimensions do not match y.")

    a_prev_fixed <- compute_a_prev_fixed(y, gamma_ref = gamma_ref, a0j = lambda_a0)

    exact_weights <- matrix(NA_real_, nrow = nd, ncol = length(grid))
    frozen_weights <- matrix(NA_real_, nrow = nd, ncol = length(grid))
    exact_mode <- frozen_mode <- exact_mean <- frozen_mean <- numeric(nd)

    for (d in seq_len(nd)) {
        lambda_draw <- lambda_arr[d, , ]
        exact_out <- compute_gamma_curve_for_one_draw(lambda_draw, y, grid, lambda_a0, prior_ab, a_prev_fixed = NULL)
        frozen_out <- compute_gamma_curve_for_one_draw(lambda_draw, y, grid, lambda_a0, prior_ab, a_prev_fixed = a_prev_fixed)
        exact_weights[d, ]  <- exact_out$weights
        frozen_weights[d, ] <- frozen_out$weights
        exact_mode[d]   <- exact_out$curve$mode
        exact_mean[d]   <- exact_out$curve$mean
        frozen_mode[d]  <- frozen_out$curve$mode
        frozen_mean[d]  <- frozen_out$curve$mean
    }

    exact_avg_w  <- colMeans(exact_weights)
    frozen_avg_w <- colMeans(frozen_weights)
    exact_avg_curve  <- summarize_curve(grid, exact_avg_w)
    frozen_avg_curve <- summarize_curve(grid, frozen_avg_w)

    exact_summary <- list(conditional_mode_mean = mean(exact_mode), conditional_mean_mean = mean(exact_mean))
    frozen_summary <- list(conditional_mode_mean = mean(frozen_mode), conditional_mean_mean = mean(frozen_mean))
    gamma_draws <- extract_gamma_draws_compat(fit)

    if (save_png) {
        png(sprintf("%s.png", out_prefix), width = 1400, height = 1100, res = 140)
        par(mfrow = c(2,2), mar = c(4,4,3,1))
        plot(grid, exact_avg_w, type = "l", lwd = 2, main = "Average conditional gamma curves",
             xlab = expression(gamma), ylab = "Posterior mass")
        lines(grid, frozen_avg_w, lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact a(gamma)", "frozen a", if (!is.null(gamma_truth)) sprintf("truth = %.3f", gamma_truth) else NULL),
               lty = c(1,2,if (!is.null(gamma_truth)) 1 else NULL), lwd = c(2,2,if (!is.null(gamma_truth)) 2 else NULL), bty = "n")
        plot(density(exact_mode), lwd = 2, main = "Drawwise conditional modes", xlab = expression(gamma), ylab = "Density")
        lines(density(frozen_mode), lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact", "frozen a"), lty = c(1,2), lwd = 2, bty = "n")
        plot(density(exact_mean), lwd = 2, main = "Drawwise conditional means", xlab = expression(gamma), ylab = "Density")
        lines(density(frozen_mean), lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact", "frozen a"), lty = c(1,2), lwd = 2, bty = "n")
        if (!is.null(gamma_draws)) {
            hist(gamma_draws, breaks = 25, freq = FALSE, main = "Full MCMC gamma draws", xlab = expression(gamma), border = "grey30", col = "grey85")
            lines(density(gamma_draws), lwd = 2)
            if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        } else { plot.new(); title("Full MCMC gamma draws\n(not available)") }
        dev.off()
    }

    if (save_pdf) {
        pdf(sprintf("%s.pdf", out_prefix), width = 10, height = 8)
        par(mfrow = c(2,2), mar = c(4,4,3,1))
        plot(grid, exact_avg_w, type = "l", lwd = 2, main = "Average conditional gamma curves",
             xlab = expression(gamma), ylab = "Posterior mass")
        lines(grid, frozen_avg_w, lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact a(gamma)", "frozen a", if (!is.null(gamma_truth)) sprintf("truth = %.3f", gamma_truth) else NULL),
               lty = c(1,2,if (!is.null(gamma_truth)) 1 else NULL), lwd = c(2,2,if (!is.null(gamma_truth)) 2 else NULL), bty = "n")
        plot(density(exact_mode), lwd = 2, main = "Drawwise conditional modes", xlab = expression(gamma), ylab = "Density")
        lines(density(frozen_mode), lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact", "frozen a"), lty = c(1,2), lwd = 2, bty = "n")
        plot(density(exact_mean), lwd = 2, main = "Drawwise conditional means", xlab = expression(gamma), ylab = "Density")
        lines(density(frozen_mean), lwd = 2, lty = 2)
        if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        legend("topright", legend = c("exact", "frozen a"), lty = c(1,2), lwd = 2, bty = "n")
        if (!is.null(gamma_draws)) {
            hist(gamma_draws, breaks = 25, freq = FALSE, main = "Full MCMC gamma draws", xlab = expression(gamma), border = "grey30", col = "grey85")
            lines(density(gamma_draws), lwd = 2)
            if (!is.null(gamma_truth)) abline(v = gamma_truth, lwd = 2)
        } else { plot.new(); title("Full MCMC gamma draws\n(not available)") }
        dev.off()
    }

    summary_df <- data.frame(
        diagnostic = c(
            "exact_avg_mode", "exact_avg_mean", "exact_avg_q025", "exact_avg_q975",
            "frozen_avg_mode", "frozen_avg_mean", "frozen_avg_q025", "frozen_avg_q975",
            "exact_draw_mean_of_modes", "exact_draw_mean_of_means",
            "frozen_draw_mean_of_modes", "frozen_draw_mean_of_means",
            "gamma_truth", "gamma_ref"
        ),
        value = c(
            exact_avg_curve$mode, exact_avg_curve$mean, exact_avg_curve$q025, exact_avg_curve$q975,
            frozen_avg_curve$mode, frozen_avg_curve$mean, frozen_avg_curve$q025, frozen_avg_curve$q975,
            exact_summary$conditional_mode_mean, exact_summary$conditional_mean_mean,
            frozen_summary$conditional_mode_mean, frozen_summary$conditional_mean_mean,
            gamma_truth %||% NA_real_, gamma_ref
        )
    )
    write.csv(summary_df, sprintf("%s_summary.csv", out_prefix), row.names = FALSE)

    out <- list(
        grid = grid,
        gamma_truth = gamma_truth,
        gamma_ref = gamma_ref,
        prior = prior_ab,
        lambda_raw_source = if (!is.null(fit$samples$lambda_raw_for_gamma)) "samples$lambda_raw_for_gamma" else "fallback",
        exact = list(weights_by_draw = exact_weights, average_weights = exact_avg_w, average_curve = exact_avg_curve, summary = exact_summary),
        frozen = list(weights_by_draw = frozen_weights, average_weights = frozen_avg_w, average_curve = frozen_avg_curve, summary = frozen_summary)
    )
    saveRDS(out, sprintf("%s_diagnostic.rds", out_prefix))

    cat("Frozen-a gamma diagnostic\n")
    cat("  gamma_ref used for frozen a:", sprintf("%.6f", gamma_ref), "\n")
    if (!is.null(gamma_truth)) cat("  gamma truth              :", sprintf("%.6f", gamma_truth), "\n")
    cat("  lambda raw source         :", out$lambda_raw_source, "\n")
    cat("\nExact complete-data target\n")
    cat("  average-curve mode       :", sprintf("%.6f", exact_avg_curve$mode), "\n")
    cat("  average-curve mean       :", sprintf("%.6f", exact_avg_curve$mean), "\n")
    cat("  average-curve 95% CI     : [", sprintf("%.6f", exact_avg_curve$q025), ", ", sprintf("%.6f", exact_avg_curve$q975), "]\n", sep = "")
    cat("  drawwise mode mean       :", sprintf("%.6f", exact_summary$conditional_mode_mean), "\n")
    cat("  drawwise mean mean       :", sprintf("%.6f", exact_summary$conditional_mean_mean), "\n")
    cat("\nFrozen-a target\n")
    cat("  average-curve mode       :", sprintf("%.6f", frozen_avg_curve$mode), "\n")
    cat("  average-curve mean       :", sprintf("%.6f", frozen_avg_curve$mean), "\n")
    cat("  average-curve 95% CI     : [", sprintf("%.6f", frozen_avg_curve$q025), ", ", sprintf("%.6f", frozen_avg_curve$q975), "]\n", sep = "")
    cat("  drawwise mode mean       :", sprintf("%.6f", frozen_summary$conditional_mode_mean), "\n")
    cat("  drawwise mean mean       :", sprintf("%.6f", frozen_summary$conditional_mean_mean), "\n")

    invisible(out)
}
