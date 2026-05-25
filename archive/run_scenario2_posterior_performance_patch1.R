## ============================================================================
## run_scenario2_posterior_performance.R
##
## Posterior performance tables and figures for Scenario 2:
## dynamic residual risk with fixed known gamma.
##
## This script mirrors the Scenario 1 posterior-performance workflow, but adds
## lambda path recovery diagnostics and plots.
##
## Typical use:
##   source("run_scenario2_posterior_performance.R")
##   out <- run_scenario2_posterior_performance(
##       root = ".",
##       reps = 1:20,
##       scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
##       data_dir = "data_revised",
##       output_dir = "output_s2_dynamic_fixed_gamma",
##       analysis_dir = "analysis_s2_dynamic_fixed_gamma",
##       plot_format = "png"
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s2run <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.safe_quantile_s2 <- function(x, probs = c(0.025, 0.5, 0.975)) {
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

.safe_cor_s2 <- function(x, y) {
    out <- suppressWarnings(stats::cor(as.vector(x), as.vector(y), use = "complete.obs"))
    if (!is.finite(out)) NA_real_ else out
}

.scalar_s2 <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) {
        return(default)
    }
    x <- x[1L]
    if (is.factor(x)) {
        x <- as.character(x)
    }
    x
}

.num_scalar_s2 <- function(x, default = NA_real_) {
    out <- .scalar_s2(x, default = default)
    suppressWarnings(as.numeric(out))
}

.char_scalar_s2 <- function(x, default = NA_character_) {
    out <- .scalar_s2(x, default = default)
    as.character(out)
}


.open_plot_device_s2 <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.close_plot_device_s2 <- function() {
    grDevices::dev.off()
}

.get_fit_file_s2 <- function(root, output_dir, scenario_id, rep_id) {
    file.path(
        root,
        output_dir,
        scenario_id,
        sprintf("fit_S2_dynamic_fixed_gamma_rep%02d.rds", as.integer(rep_id))
    )
}

.get_data_file_s2 <- function(root, data_dir, scenario_id, rep_id) {
    file.path(
        root,
        data_dir,
        scenario_id,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

.get_lambda_draws_s2 <- function(samples) {
    lam <- samples$lambda_tilde
    if (is.null(lam)) {
        stop("fit$samples$lambda_tilde is missing.", call. = FALSE)
    }
    if (length(dim(lam)) != 3L) {
        stop("fit$samples$lambda_tilde must be a 3D array: draws x time x region.", call. = FALSE)
    }
    lam
}

## ----------------------------------------------------------------------------
## Collect posterior summaries from saved fit/data files
## ----------------------------------------------------------------------------
collect_s2_posterior_performance <- function(reps = 1:20,
                                             scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
                                             data_dir = "data_revised",
                                             output_dir = "output_s2_dynamic_fixed_gamma",
                                             root = ".") {
    beta_list <- list()
    phi_list <- list()
    r_list <- list()
    lambda_list <- list()
    diag_list <- list()
    conf_list <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat_file <- .get_data_file_s2(root, data_dir, scenario_id, rep_id)
        fit_file <- .get_fit_file_s2(root, output_dir, scenario_id, rep_id)
        .require_file_s2run(dat_file)
        .require_file_s2run(fit_file)

        dat <- readRDS(dat_file)
        fit <- readRDS(fit_file)
        s <- fit$samples
        if (is.null(s)) {
            stop("Fit file has no $samples component: ", fit_file, call. = FALSE)
        }

        beta0_true <- fit$summary$beta0_true %||% dat$beta0_star_ident %||% dat$beta0_star
        beta_true <- c(
            beta1 = fit$summary$beta1_true %||% dat$beta_star[1],
            beta2 = fit$summary$beta2_true %||% dat$beta_star[2]
        )
        phi_true <- dat$phi_star_ident %||% dat$phi_star
        r_true <- dat$r_star
        if (length(r_true) == 1L) {
            r_true <- rep(r_true, dat$n1 %||% ncol(dat$y_coarse))
        }
        lambda_true <- dat$lambda_tilde_ident %||% dat$lambda_tilde

        beta_draws <- cbind(
            beta0 = s$beta0,
            beta1 = s$beta[, 1],
            beta2 = s$beta[, 2]
        )
        beta_truth <- c(beta0 = beta0_true, beta1 = beta_true[1], beta2 = beta_true[2])

        beta_df <- do.call(rbind, lapply(seq_along(beta_truth), function(k) {
            q <- .safe_quantile_s2(beta_draws[, k])
            data.frame(
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                parameter = names(beta_truth)[k],
                truth = as.numeric(beta_truth[k]),
                post_mean = mean(beta_draws[, k], na.rm = TRUE),
                post_sd = stats::sd(beta_draws[, k], na.rm = TRUE),
                q025 = q[1],
                q50 = q[2],
                q975 = q[3],
                bias = mean(beta_draws[, k], na.rm = TRUE) - as.numeric(beta_truth[k]),
                covered = as.integer(q[1] <= beta_truth[k] && beta_truth[k] <= q[3]),
                stringsAsFactors = FALSE
            )
        }))
        beta_list[[rr]] <- beta_df

        phi_mean <- colMeans(s$phi, na.rm = TRUE)
        phi_q <- apply(s$phi, 2, .safe_quantile_s2)
        phi_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            region = seq_along(phi_true),
            truth = as.numeric(phi_true),
            mean = as.numeric(phi_mean),
            q025 = as.numeric(phi_q[1, ]),
            q50 = as.numeric(phi_q[2, ]),
            q975 = as.numeric(phi_q[3, ]),
            bias = as.numeric(phi_mean - phi_true),
            covered = as.integer(phi_q[1, ] <= phi_true & phi_true <= phi_q[3, ]),
            stringsAsFactors = FALSE
        )
        phi_list[[rr]] <- phi_df

        r_mean <- colMeans(s$r, na.rm = TRUE)
        r_q <- apply(s$r, 2, .safe_quantile_s2)
        r_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            region = seq_along(r_true),
            truth = as.numeric(r_true),
            mean = as.numeric(r_mean),
            q025 = as.numeric(r_q[1, ]),
            q50 = as.numeric(r_q[2, ]),
            q975 = as.numeric(r_q[3, ]),
            bias = as.numeric(r_mean - r_true),
            covered = as.integer(r_q[1, ] <= r_true & r_true <= r_q[3, ]),
            stringsAsFactors = FALSE
        )
        r_list[[rr]] <- r_df

        lambda_draws <- .get_lambda_draws_s2(s)
        lambda_mean <- apply(lambda_draws, c(2, 3), mean, na.rm = TRUE)
        lambda_q025 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
        lambda_q50 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.50, na.rm = TRUE)
        lambda_q975 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)

        log_lambda_draws <- log(pmax(lambda_draws, .Machine$double.eps))
        log_lambda_mean <- apply(log_lambda_draws, c(2, 3), mean, na.rm = TRUE)
        log_lambda_q025 <- apply(log_lambda_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
        log_lambda_q975 <- apply(log_lambda_draws, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)

        lambda_true_vec <- as.vector(lambda_true)
        lambda_mean_vec <- as.vector(lambda_mean)
        log_lambda_true <- log(pmax(lambda_true, .Machine$double.eps))

        delta_true <- diff(log_lambda_true)
        delta_mean <- diff(log_lambda_mean)

        lambda_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            TT = nrow(lambda_true),
            n1 = ncol(lambda_true),
            gamma_truth_mean = .num_scalar_s2(mean(dat$gamma_star %||% NA_real_, na.rm = TRUE)),
            gamma_fixed_mean = .num_scalar_s2(mean(fit$samples$gamma %||% fit$summary$gamma_fixed_mean %||% NA_real_, na.rm = TRUE)),
            lambda_truth_min = min(lambda_true, na.rm = TRUE),
            lambda_truth_max = max(lambda_true, na.rm = TRUE),
            lambda_post_mean_min = min(lambda_mean, na.rm = TRUE),
            lambda_post_mean_max = max(lambda_mean, na.rm = TRUE),
            lambda_rmse = sqrt(mean((lambda_mean_vec - lambda_true_vec)^2, na.rm = TRUE)),
            lambda_mae = mean(abs(lambda_mean_vec - lambda_true_vec), na.rm = TRUE),
            lambda_coverage_95 = mean(as.vector(lambda_q025) <= lambda_true_vec & lambda_true_vec <= as.vector(lambda_q975), na.rm = TRUE),
            log_lambda_rmse = sqrt(mean((as.vector(log_lambda_mean) - as.vector(log_lambda_true))^2, na.rm = TRUE)),
            log_lambda_mae = mean(abs(as.vector(log_lambda_mean) - as.vector(log_lambda_true)), na.rm = TRUE),
            log_lambda_coverage_95 = mean(as.vector(log_lambda_q025) <= as.vector(log_lambda_true) &
                                              as.vector(log_lambda_true) <= as.vector(log_lambda_q975), na.rm = TRUE),
            cor_lambda = .safe_cor_s2(lambda_mean, lambda_true),
            cor_log_lambda = .safe_cor_s2(log_lambda_mean, log_lambda_true),
            rmse_delta_log_lambda = sqrt(mean((as.vector(delta_mean) - as.vector(delta_true))^2, na.rm = TRUE)),
            mae_delta_log_lambda = mean(abs(as.vector(delta_mean) - as.vector(delta_true)), na.rm = TRUE),
            cor_delta_log_lambda = .safe_cor_s2(delta_mean, delta_true),
            stringsAsFactors = FALSE
        )
        lambda_list[[rr]] <- lambda_df

        x1 <- as.vector(dat$x1)
        x2 <- as.vector(dat$x2)
        loglam <- as.vector(log_lambda_true)
        phi_mat <- matrix(phi_true, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)
        conf_list[[rr]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            x2_mode = .char_scalar_s2(dat$x2_mode %||% NA_character_),
            sd_x1 = stats::sd(x1),
            sd_x2 = stats::sd(x2),
            cor_x1_x2 = .safe_cor_s2(x1, x2),
            cor_x2_loglambda = .safe_cor_s2(x2, loglam),
            cor_x1_loglambda = .safe_cor_s2(x1, loglam),
            cor_x2_phi = .safe_cor_s2(x2, as.vector(phi_mat)),
            cor_x1_phi = .safe_cor_s2(x1, as.vector(phi_mat)),
            sd_beta2_x2 = stats::sd(dat$beta_star[2] * x2),
            stringsAsFactors = FALSE
        )

        diag_list[[rr]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            TT = .num_scalar_s2(dat$TT %||% nrow(dat$y_coarse)),
            n1 = .num_scalar_s2(dat$n1 %||% ncol(dat$y_coarse)),
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            dynamic_lambda_in_truth = TRUE,
            lambda_sampled_in_fit = TRUE,
            gamma_fixed_in_fit = TRUE,
            gamma_truth_mean = .num_scalar_s2(lambda_df$gamma_truth_mean),
            gamma_fixed_mean = .num_scalar_s2(lambda_df$gamma_fixed_mean),
            lambda_truth_min = .num_scalar_s2(lambda_df$lambda_truth_min),
            lambda_truth_max = .num_scalar_s2(lambda_df$lambda_truth_max),
            lambda_post_mean_min = .num_scalar_s2(lambda_df$lambda_post_mean_min),
            lambda_post_mean_max = .num_scalar_s2(lambda_df$lambda_post_mean_max),
            phi_rmse = sqrt(mean((phi_mean - phi_true)^2, na.rm = TRUE)),
            phi_cor = .safe_cor_s2(phi_mean, phi_true),
            r_rmse = sqrt(mean((r_mean - r_true)^2, na.rm = TRUE)),
            r_mean_bias = mean(r_mean - r_true, na.rm = TRUE),
            beta0_bias = .num_scalar_s2(beta_df$bias[beta_df$parameter == "beta0"]),
            beta1_bias = .num_scalar_s2(beta_df$bias[beta_df$parameter == "beta1"]),
            beta2_bias = .num_scalar_s2(beta_df$bias[beta_df$parameter == "beta2"]),
            beta0_covered = .num_scalar_s2(beta_df$covered[beta_df$parameter == "beta0"]),
            beta1_covered = .num_scalar_s2(beta_df$covered[beta_df$parameter == "beta1"]),
            beta2_covered = .num_scalar_s2(beta_df$covered[beta_df$parameter == "beta2"]),
            phi_coverage = mean(phi_df$covered, na.rm = TRUE),
            r_coverage = mean(r_df$covered, na.rm = TRUE),
            log_lambda_rmse = .num_scalar_s2(lambda_df$log_lambda_rmse),
            log_lambda_coverage_95 = .num_scalar_s2(lambda_df$log_lambda_coverage_95),
            cor_log_lambda = .num_scalar_s2(lambda_df$cor_log_lambda),
            cor_delta_log_lambda = .num_scalar_s2(lambda_df$cor_delta_log_lambda),
            phi_accept_rate = .num_scalar_s2(fit$diagnostics$phi_accept_rate %||% NA_real_),
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            beta_mean_n_reject = .num_scalar_s2(fit$diagnostics$beta_mean_n_reject %||% NA_real_),
            elapsed_sec = .num_scalar_s2(fit$diagnostics$elapsed_sec %||% NA_real_),
            stringsAsFactors = FALSE
        )
    }

    list(
        beta = do.call(rbind, beta_list),
        phi = do.call(rbind, phi_list),
        r = do.call(rbind, r_list),
        lambda = do.call(rbind, lambda_list),
        diagnostics = do.call(rbind, diag_list),
        confounding = do.call(rbind, conf_list)
    )
}

save_s2_performance_tables <- function(perf, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    files <- list(
        beta = file.path(out_dir, "posterior_beta_summary.csv"),
        phi = file.path(out_dir, "posterior_phi_summary.csv"),
        r = file.path(out_dir, "posterior_r_summary.csv"),
        lambda = file.path(out_dir, "posterior_lambda_path_recovery.csv"),
        diagnostics = file.path(out_dir, "posterior_performance_diagnostics.csv"),
        confounding = file.path(out_dir, "scenario2_confounding_diagnostics.csv")
    )
    utils::write.csv(perf$beta, files$beta, row.names = FALSE)
    utils::write.csv(perf$phi, files$phi, row.names = FALSE)
    utils::write.csv(perf$r, files$r, row.names = FALSE)
    utils::write.csv(perf$lambda, files$lambda, row.names = FALSE)
    utils::write.csv(perf$diagnostics, files$diagnostics, row.names = FALSE)
    utils::write.csv(perf$confounding, files$confounding, row.names = FALSE)
    invisible(files)
}

## ----------------------------------------------------------------------------
## Base R plots
## ----------------------------------------------------------------------------
plot_s2_beta_recovery <- function(beta_df, file) {
    .open_plot_device_s2(file, width = 10, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)

    params <- unique(beta_df$parameter)
    par(mfrow = c(1, length(params)), mar = c(4, 4, 3, 1))
    for (param in params) {
        d <- beta_df[beta_df$parameter == param, ]
        d <- d[order(d$rep_id), ]
        y_lim <- range(c(d$q025, d$q975, d$truth), finite = TRUE)
        pad <- diff(y_lim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        y_lim <- y_lim + c(-pad, pad)
        plot(d$rep_id, d$post_mean, ylim = y_lim,
             xlab = "Replicate", ylab = "Posterior mean and 95% interval",
             main = param, pch = 19)
        segments(d$rep_id, d$q025, d$rep_id, d$q975)
        if (length(unique(round(d$truth, 10))) == 1L) {
            abline(h = unique(d$truth)[1], lty = 2, lwd = 2)
        } else {
            lines(d$rep_id, d$truth, lty = 2, lwd = 2)
        }
        grid()
    }
    mtext("Scenario 2: beta recovery. Dashed line is truth.", outer = TRUE, line = -1.5, cex = 1.1)
}

plot_s2_beta_bias <- function(beta_df, file) {
    .open_plot_device_s2(file, width = 8, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    boxplot(bias ~ parameter, data = beta_df,
            xlab = "Parameter", ylab = "Posterior mean minus truth",
            main = "Scenario 2: beta bias across replicates")
    abline(h = 0, lty = 2, lwd = 2)
    grid()
}

plot_s2_phi_scatter <- function(phi_df, file, reps_to_show = NULL) {
    if (is.null(reps_to_show)) {
        reps_to_show <- sort(unique(phi_df$rep_id))[1:min(4, length(unique(phi_df$rep_id)))]
    }
    d_all <- phi_df[phi_df$rep_id %in% reps_to_show, ]
    n_panels <- length(reps_to_show)
    nr <- ceiling(sqrt(n_panels))
    nc <- ceiling(n_panels / nr)
    .open_plot_device_s2(file, width = 4.5 * nc, height = 4.2 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (rep_id in reps_to_show) {
        d <- d_all[d_all$rep_id == rep_id, ]
        lim <- range(c(d$truth, d$mean, d$q025, d$q975), finite = TRUE)
        pad <- diff(lim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        lim <- lim + c(-pad, pad)
        plot(d$truth, d$mean, xlim = lim, ylim = lim, pch = 19,
             xlab = "True phi", ylab = "Posterior mean phi", main = paste("rep", rep_id))
        segments(d$truth, d$q025, d$truth, d$q975)
        abline(0, 1, lty = 2, lwd = 2)
        grid()
    }
}

plot_s2_phi_metrics <- function(diag_df, file) {
    .open_plot_device_s2(file, width = 9, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$phi_rmse, type = "b", pch = 19,
         xlab = "Replicate", ylab = "RMSE", main = "phi RMSE")
    grid()
    plot(diag_df$rep_id, diag_df$phi_cor, type = "b", pch = 19,
         ylim = c(min(c(0, diag_df$phi_cor), na.rm = TRUE), 1),
         xlab = "Replicate", ylab = "Correlation", main = "cor(phi posterior mean, truth)")
    abline(h = 1, lty = 2)
    grid()
}

plot_s2_r_recovery <- function(r_df, file, reps_to_show = NULL) {
    if (is.null(reps_to_show)) {
        reps_to_show <- sort(unique(r_df$rep_id))[1:min(4, length(unique(r_df$rep_id)))]
    }
    d_all <- r_df[r_df$rep_id %in% reps_to_show, ]
    n_panels <- length(reps_to_show)
    nr <- ceiling(sqrt(n_panels))
    nc <- ceiling(n_panels / nr)
    .open_plot_device_s2(file, width = 4.8 * nc, height = 4.0 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (rep_id in reps_to_show) {
        d <- d_all[d_all$rep_id == rep_id, ]
        ylim <- range(c(d$q025, d$q975, d$truth), finite = TRUE)
        pad <- diff(ylim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        ylim <- ylim + c(-pad, pad)
        plot(d$region, d$mean, ylim = ylim, pch = 19,
             xlab = "Region", ylab = "r posterior mean and 95% interval", main = paste("rep", rep_id))
        segments(d$region, d$q025, d$region, d$q975)
        lines(d$region, d$truth, lty = 2, lwd = 2)
        grid()
    }
}

plot_s2_r_metrics <- function(diag_df, file) {
    .open_plot_device_s2(file, width = 9, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$r_rmse, type = "b", pch = 19,
         xlab = "Replicate", ylab = "RMSE", main = "r RMSE")
    grid()
    plot(diag_df$rep_id, diag_df$r_mean_bias, type = "b", pch = 19,
         xlab = "Replicate", ylab = "Mean r bias", main = "mean posterior r minus truth")
    abline(h = 0, lty = 2, lwd = 2)
    grid()
}

plot_s2_lambda_metrics <- function(lambda_df, file) {
    .open_plot_device_s2(file, width = 10, height = 8)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    plot(lambda_df$rep_id, lambda_df$log_lambda_rmse, type = "b", pch = 19,
         xlab = "Replicate", ylab = "RMSE", main = "log(lambda) RMSE")
    grid()
    plot(lambda_df$rep_id, lambda_df$cor_log_lambda, type = "b", pch = 19,
         ylim = c(min(c(0, lambda_df$cor_log_lambda), na.rm = TRUE), 1),
         xlab = "Replicate", ylab = "Correlation", main = "cor(log lambda mean, truth)")
    abline(h = 1, lty = 2)
    grid()
    plot(lambda_df$rep_id, lambda_df$log_lambda_coverage_95, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Coverage", main = "95% coverage for log(lambda)")
    abline(h = 0.95, lty = 2)
    grid()
    plot(lambda_df$rep_id, lambda_df$cor_delta_log_lambda, type = "b", pch = 19,
         ylim = c(min(c(0, lambda_df$cor_delta_log_lambda), na.rm = TRUE), 1),
         xlab = "Replicate", ylab = "Correlation", main = "cor(delta log lambda)")
    grid()
}

plot_s2_lambda_ranges <- function(lambda_df, file) {
    .open_plot_device_s2(file, width = 9, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    ylim <- range(c(lambda_df$lambda_truth_min, lambda_df$lambda_truth_max,
                    lambda_df$lambda_post_mean_min, lambda_df$lambda_post_mean_max), finite = TRUE)
    plot(lambda_df$rep_id, lambda_df$lambda_truth_min, type = "n",
         ylim = ylim, xlab = "Replicate", ylab = "lambda range",
         main = "Scenario 2: true range vs posterior mean range")
    segments(lambda_df$rep_id - 0.12, lambda_df$lambda_truth_min,
             lambda_df$rep_id - 0.12, lambda_df$lambda_truth_max, lwd = 2)
    points(lambda_df$rep_id - 0.12, lambda_df$lambda_truth_min, pch = 16)
    points(lambda_df$rep_id - 0.12, lambda_df$lambda_truth_max, pch = 16)
    segments(lambda_df$rep_id + 0.12, lambda_df$lambda_post_mean_min,
             lambda_df$rep_id + 0.12, lambda_df$lambda_post_mean_max, lwd = 2, lty = 2)
    points(lambda_df$rep_id + 0.12, lambda_df$lambda_post_mean_min, pch = 1)
    points(lambda_df$rep_id + 0.12, lambda_df$lambda_post_mean_max, pch = 1)
    legend("topright", legend = c("Truth", "Posterior mean"), lty = c(1, 2), pch = c(16, 1), bty = "n")
    grid()
}

plot_s2_lambda_representative_paths <- function(reps_to_show,
                                                selected_regions,
                                                root,
                                                data_dir,
                                                output_dir,
                                                scenario_id,
                                                file) {
    reps_to_show <- reps_to_show[seq_len(min(length(reps_to_show), 4L))]
    selected_regions <- selected_regions[seq_len(min(length(selected_regions), 3L))]
    n_panels <- length(reps_to_show) * length(selected_regions)
    nr <- length(reps_to_show)
    nc <- length(selected_regions)

    .open_plot_device_s2(file, width = 4.2 * nc, height = 3.2 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))

    for (rep_id in reps_to_show) {
        dat <- readRDS(.get_data_file_s2(root, data_dir, scenario_id, rep_id))
        fit <- readRDS(.get_fit_file_s2(root, output_dir, scenario_id, rep_id))
        lambda_true <- dat$lambda_tilde_ident %||% dat$lambda_tilde
        lambda_mean <- apply(.get_lambda_draws_s2(fit$samples), c(2, 3), mean, na.rm = TRUE)
        for (j in selected_regions) {
            ylim <- range(c(lambda_true[, j], lambda_mean[, j]), finite = TRUE)
            plot(seq_len(nrow(lambda_true)), lambda_true[, j], type = "l", lwd = 2,
                 ylim = ylim, xlab = "Time", ylab = "lambda",
                 main = paste0("rep ", rep_id, ", region ", j))
            lines(seq_len(nrow(lambda_true)), lambda_mean[, j], lty = 2, lwd = 2)
            legend("topright", legend = c("truth", "posterior mean"), lty = c(1, 2), bty = "n", cex = 0.8)
            grid()
        }
    }
}

plot_s2_mcmc_diagnostics <- function(diag_df, file) {
    .open_plot_device_s2(file, width = 10, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$phi_accept_rate, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Acceptance rate", main = "phi acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$r_accept_rate_mean, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Acceptance rate", main = "r acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$beta_mean_n_reject, type = "b", pch = 19,
         xlab = "Replicate", ylab = "Mean beta ESS rejections", main = "beta update diagnostic")
    grid()
}

plot_s2_confounding <- function(conf_df, file) {
    .open_plot_device_s2(file, width = 10, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2() }, add = TRUE)
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    plot(conf_df$rep_id, conf_df$cor_x1_x2, type = "b", pch = 19,
         ylim = range(c(-1, 1, conf_df$cor_x1_x2), finite = TRUE),
         xlab = "Replicate", ylab = "Correlation", main = "cor(x1, x2)")
    abline(h = 0, lty = 2); grid()
    plot(conf_df$rep_id, conf_df$cor_x2_loglambda, type = "b", pch = 19,
         ylim = range(c(-1, 1, conf_df$cor_x2_loglambda), finite = TRUE),
         xlab = "Replicate", ylab = "Correlation", main = "cor(x2, log lambda)")
    abline(h = 0, lty = 2); grid()
    plot(conf_df$rep_id, conf_df$cor_x2_phi, type = "b", pch = 19,
         ylim = range(c(-1, 1, conf_df$cor_x2_phi), finite = TRUE),
         xlab = "Replicate", ylab = "Correlation", main = "cor(x2, phi)")
    abline(h = 0, lty = 2); grid()
}

make_s2_posterior_performance_plots <- function(perf,
                                                fig_dir,
                                                plot_format = c("pdf", "png"),
                                                reps_to_show = NULL,
                                                selected_regions = c(1L, 5L, 9L),
                                                root = ".",
                                                data_dir = "data_revised",
                                                output_dir = "output_s2_dynamic_fixed_gamma",
                                                scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100") {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    plot_format <- match.arg(plot_format)
    ext <- plot_format
    if (is.null(reps_to_show)) {
        reps_to_show <- sort(unique(perf$diagnostics$rep_id))[1:min(4, length(unique(perf$diagnostics$rep_id)))]
    }

    files <- list(
        beta_recovery = file.path(fig_dir, paste0("s2_beta_recovery_intervals.", ext)),
        beta_bias = file.path(fig_dir, paste0("s2_beta_bias_boxplot.", ext)),
        phi_scatter = file.path(fig_dir, paste0("s2_phi_truth_vs_posterior_mean.", ext)),
        phi_metrics = file.path(fig_dir, paste0("s2_phi_metrics.", ext)),
        r_recovery = file.path(fig_dir, paste0("s2_r_recovery_intervals.", ext)),
        r_metrics = file.path(fig_dir, paste0("s2_r_metrics.", ext)),
        lambda_metrics = file.path(fig_dir, paste0("s2_lambda_path_recovery_metrics.", ext)),
        lambda_ranges = file.path(fig_dir, paste0("s2_lambda_range_shrinkage.", ext)),
        lambda_paths = file.path(fig_dir, paste0("s2_lambda_representative_paths.", ext)),
        mcmc_diagnostics = file.path(fig_dir, paste0("s2_mcmc_diagnostics.", ext)),
        confounding = file.path(fig_dir, paste0("s2_confounding_diagnostics.", ext))
    )

    plot_s2_beta_recovery(perf$beta, files$beta_recovery)
    plot_s2_beta_bias(perf$beta, files$beta_bias)
    plot_s2_phi_scatter(perf$phi, files$phi_scatter, reps_to_show = reps_to_show)
    plot_s2_phi_metrics(perf$diagnostics, files$phi_metrics)
    plot_s2_r_recovery(perf$r, files$r_recovery, reps_to_show = reps_to_show)
    plot_s2_r_metrics(perf$diagnostics, files$r_metrics)
    plot_s2_lambda_metrics(perf$lambda, files$lambda_metrics)
    plot_s2_lambda_ranges(perf$lambda, files$lambda_ranges)
    plot_s2_lambda_representative_paths(
        reps_to_show = reps_to_show,
        selected_regions = selected_regions,
        root = root,
        data_dir = data_dir,
        output_dir = output_dir,
        scenario_id = scenario_id,
        file = files$lambda_paths
    )
    plot_s2_mcmc_diagnostics(perf$diagnostics, files$mcmc_diagnostics)
    plot_s2_confounding(perf$confounding, files$confounding)

    invisible(files)
}

summarise_s2_performance_in_console <- function(perf) {
    beta <- perf$beta
    diag <- perf$diagnostics
    lambda <- perf$lambda

    cat("\n================ Scenario 2 posterior performance ================\n")
    cat("Dynamic lambda in truth: ", all(diag$dynamic_lambda_in_truth), "\n", sep = "")
    cat("Lambda sampled in fit  : ", all(diag$lambda_sampled_in_fit), "\n", sep = "")
    cat("Gamma fixed in fit     : ", all(diag$gamma_fixed_in_fit), "\n", sep = "")
    cat("Gamma truth range      : ", paste(range(diag$gamma_truth_mean), collapse = " to "), "\n", sep = "")

    cat("\nBeta recovery by parameter:\n")
    beta_summary <- do.call(rbind, lapply(split(beta, beta$parameter), function(d) {
        data.frame(
            parameter = unique(d$parameter),
            mean_bias = mean(d$bias, na.rm = TRUE),
            rmse_bias = sqrt(mean(d$bias^2, na.rm = TRUE)),
            coverage = mean(d$covered, na.rm = TRUE),
            mean_post_sd = mean(d$post_sd, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))
    print(beta_summary, row.names = FALSE)

    cat("\nLambda path recovery across replicates:\n")
    cat(sprintf("  mean log-lambda RMSE      : %.4f\n", mean(lambda$log_lambda_rmse, na.rm = TRUE)))
    cat(sprintf("  mean log-lambda coverage  : %.4f\n", mean(lambda$log_lambda_coverage_95, na.rm = TRUE)))
    cat(sprintf("  mean cor(log lambda)      : %.4f\n", mean(lambda$cor_log_lambda, na.rm = TRUE)))
    cat(sprintf("  mean cor(delta log lambda): %.4f\n", mean(lambda$cor_delta_log_lambda, na.rm = TRUE)))

    cat("\nPhi recovery across replicates:\n")
    cat(sprintf("  mean phi RMSE       : %.4f\n", mean(diag$phi_rmse, na.rm = TRUE)))
    cat(sprintf("  mean phi correlation: %.4f\n", mean(diag$phi_cor, na.rm = TRUE)))
    cat(sprintf("  mean phi coverage   : %.4f\n", mean(diag$phi_coverage, na.rm = TRUE)))

    cat("\nr recovery across replicates:\n")
    cat(sprintf("  mean r RMSE      : %.4f\n", mean(diag$r_rmse, na.rm = TRUE)))
    cat(sprintf("  mean r bias      : %.4f\n", mean(diag$r_mean_bias, na.rm = TRUE)))
    cat(sprintf("  mean r coverage  : %.4f\n", mean(diag$r_coverage, na.rm = TRUE)))

    cat("\nMCMC diagnostics:\n")
    cat(sprintf("  mean phi acceptance: %.4f\n", mean(diag$phi_accept_rate, na.rm = TRUE)))
    cat(sprintf("  mean r acceptance  : %.4f\n", mean(diag$r_accept_rate_mean, na.rm = TRUE)))
    cat("====================================================================\n\n")

    invisible(beta_summary)
}

## ----------------------------------------------------------------------------
## Main one-stop driver. It does not rerun MCMC by default. It collects saved fits.
## ----------------------------------------------------------------------------
run_scenario2_posterior_performance <- function(root = ".",
                                                reps = 1:20,
                                                scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
                                                data_dir = "data_revised",
                                                output_dir = "output_s2_dynamic_fixed_gamma",
                                                analysis_dir = "analysis_s2_dynamic_fixed_gamma",
                                                plot_format = c("pdf", "png"),
                                                reps_to_show = NULL,
                                                selected_regions = c(1L, 5L, 9L)) {
    plot_format <- match.arg(plot_format)

    analysis_root <- file.path(root, analysis_dir, scenario_id)
    tab_dir <- file.path(analysis_root, "tables")
    fig_dir <- file.path(analysis_root, "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    perf <- collect_s2_posterior_performance(
        reps = reps,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root
    )

    table_files <- save_s2_performance_tables(perf, tab_dir)
    plot_files <- make_s2_posterior_performance_plots(
        perf = perf,
        fig_dir = fig_dir,
        plot_format = plot_format,
        reps_to_show = reps_to_show,
        selected_regions = selected_regions,
        root = root,
        data_dir = data_dir,
        output_dir = output_dir,
        scenario_id = scenario_id
    )

    summary <- summarise_s2_performance_in_console(perf)

    invisible(list(
        perf = perf,
        tables = table_files,
        plots = plot_files,
        summary = summary,
        analysis_dir = analysis_root
    ))
}

