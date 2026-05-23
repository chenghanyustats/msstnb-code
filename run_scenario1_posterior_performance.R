## ============================================================================
## run_scenario1_posterior_performance.R
##
## One-stop driver for Scenario 1 baseline-only simulation, fitting, and
## posterior performance plots for beta, phi, and r.
##
## Required companion file:
##   s1_baseline_only_clean.R
##
## Recommended use from the MSSTNB project root:
##   source("run_scenario1_posterior_performance.R")
##   out <- run_scenario1_posterior_performance(
##       root = ".",
##       reps = 1:10,
##       n_iter = 20000,
##       n_burnin = 10000,
##       n_thin = 5,
##       overwrite_data = TRUE,
##       overwrite_fit = TRUE
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s1run <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.safe_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

.open_plot_device <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.close_plot_device <- function() {
    grDevices::dev.off()
}

## ----------------------------------------------------------------------------
## Collect posterior summaries from saved fit/data files
## ----------------------------------------------------------------------------
collect_s1_posterior_performance <- function(reps = 1:10,
                                             scenario_id = "S1_BASELINE_ONLY",
                                             data_dir = "data_revised",
                                             output_dir = "output_s1_baseline_only",
                                             root = ".") {
    beta_list <- list()
    phi_list <- list()
    r_list <- list()
    diag_list <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
        fit_file <- file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds"))
        .require_file_s1run(dat_file)
        .require_file_s1run(fit_file)

        dat <- readRDS(dat_file)
        fit <- readRDS(fit_file)
        s <- fit$samples

        beta0_true <- dat$beta0_star_ident %||% dat$beta0_star
        beta_true <- dat$beta_star
        phi_true <- dat$phi_star_ident %||% dat$phi_star
        r_true <- dat$r_star
        if (length(r_true) == 1L) {
            r_true <- rep(r_true, dat$n1)
        }

        beta_draws <- cbind(beta0 = s$beta0,
                            beta1 = s$beta[, 1],
                            beta2 = s$beta[, 2])
        beta_truth <- c(beta0 = beta0_true,
                        beta1 = beta_true[1],
                        beta2 = beta_true[2])

        beta_df <- do.call(rbind, lapply(seq_along(beta_truth), function(k) {
            q <- .safe_quantile(beta_draws[, k])
            data.frame(
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                parameter = names(beta_truth)[k],
                truth = as.numeric(beta_truth[k]),
                mean = mean(beta_draws[, k], na.rm = TRUE),
                sd = stats::sd(beta_draws[, k], na.rm = TRUE),
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
        phi_q <- apply(s$phi, 2, .safe_quantile)
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
        r_q <- apply(s$r, 2, .safe_quantile)
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

        diag_list[[rr]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            TT = dat$TT,
            n1 = dat$n1,
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            lambda_truth_min = if (!is.null(dat$lambda_tilde)) min(dat$lambda_tilde) else NA_real_,
            lambda_truth_max = if (!is.null(dat$lambda_tilde)) max(dat$lambda_tilde) else NA_real_,
            lambda_fixed_in_fit = isTRUE(fit$metadata$fixed_lambda),
            lambda_fixed_value = fit$metadata$lambda_fixed_value %||% NA_real_,
            phi_rmse = sqrt(mean((phi_mean - phi_true)^2)),
            phi_cor = suppressWarnings(stats::cor(phi_mean, phi_true)),
            r_rmse = sqrt(mean((r_mean - r_true)^2)),
            r_mean_bias = mean(r_mean - r_true),
            beta0_bias = beta_df$bias[beta_df$parameter == "beta0"],
            beta1_bias = beta_df$bias[beta_df$parameter == "beta1"],
            beta2_bias = beta_df$bias[beta_df$parameter == "beta2"],
            beta0_covered = beta_df$covered[beta_df$parameter == "beta0"],
            beta1_covered = beta_df$covered[beta_df$parameter == "beta1"],
            beta2_covered = beta_df$covered[beta_df$parameter == "beta2"],
            phi_coverage = mean(phi_df$covered),
            r_coverage = mean(r_df$covered),
            phi_accept_rate = fit$diagnostics$phi_accept_rate %||% NA_real_,
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject %||% NA_real_,
            elapsed_sec = fit$diagnostics$elapsed_sec %||% NA_real_,
            stringsAsFactors = FALSE
        )
    }

    list(
        beta = do.call(rbind, beta_list),
        phi = do.call(rbind, phi_list),
        r = do.call(rbind, r_list),
        diagnostics = do.call(rbind, diag_list)
    )
}

save_s1_performance_tables <- function(perf, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    beta_file <- file.path(out_dir, "posterior_beta_summary.csv")
    phi_file <- file.path(out_dir, "posterior_phi_summary.csv")
    r_file <- file.path(out_dir, "posterior_r_summary.csv")
    diag_file <- file.path(out_dir, "posterior_performance_diagnostics.csv")

    write.csv(perf$beta, beta_file, row.names = FALSE)
    write.csv(perf$phi, phi_file, row.names = FALSE)
    write.csv(perf$r, r_file, row.names = FALSE)
    write.csv(perf$diagnostics, diag_file, row.names = FALSE)

    invisible(list(beta = beta_file, phi = phi_file, r = r_file, diagnostics = diag_file))
}

## ----------------------------------------------------------------------------
## Base R plots.  These intentionally avoid extra package dependencies.
## ----------------------------------------------------------------------------
plot_s1_beta_recovery <- function(beta_df, file) {
    .open_plot_device(file, width = 10, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    params <- unique(beta_df$parameter)
    par(mfrow = c(1, length(params)), mar = c(4, 4, 3, 1))

    for (param in params) {
        d <- beta_df[beta_df$parameter == param, ]
        d <- d[order(d$rep_id), ]
        y_lim <- range(c(d$q025, d$q975, d$truth), finite = TRUE)
        pad <- diff(y_lim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        y_lim <- y_lim + c(-pad, pad)

        plot(d$rep_id, d$mean,
             ylim = y_lim,
             xlab = "Replicate",
             ylab = "Posterior mean and 95% interval",
             main = param,
             pch = 19)
        segments(d$rep_id, d$q025, d$rep_id, d$q975)
        abline(h = unique(d$truth)[1], lty = 2, lwd = 2)
        grid()
    }

    mtext("Scenario 1: beta recovery. Dashed line is truth.",
          outer = TRUE, line = -1.5, cex = 1.1)
}

plot_s1_beta_bias <- function(beta_df, file) {
    .open_plot_device(file, width = 8, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    boxplot(bias ~ parameter, data = beta_df,
            xlab = "Parameter",
            ylab = "Posterior mean minus truth",
            main = "Scenario 1: beta bias across replicates")
    abline(h = 0, lty = 2, lwd = 2)
    grid()
}

plot_s1_phi_scatter <- function(phi_df, file, reps_to_show = NULL) {
    if (is.null(reps_to_show)) {
        reps_to_show <- sort(unique(phi_df$rep_id))[1:min(4, length(unique(phi_df$rep_id)))]
    }
    d_all <- phi_df[phi_df$rep_id %in% reps_to_show, ]
    n_panels <- length(reps_to_show)
    nr <- ceiling(sqrt(n_panels))
    nc <- ceiling(n_panels / nr)

    .open_plot_device(file, width = 4.5 * nc, height = 4.2 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (rep_id in reps_to_show) {
        d <- d_all[d_all$rep_id == rep_id, ]
        lim <- range(c(d$truth, d$mean, d$q025, d$q975), finite = TRUE)
        pad <- diff(lim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        lim <- lim + c(-pad, pad)
        plot(d$truth, d$mean,
             xlim = lim, ylim = lim,
             pch = 19,
             xlab = "True phi",
             ylab = "Posterior mean phi",
             main = paste("rep", rep_id))
        segments(d$truth, d$q025, d$truth, d$q975)
        abline(0, 1, lty = 2, lwd = 2)
        grid()
    }
}

plot_s1_phi_metrics <- function(diag_df, file) {
    .open_plot_device(file, width = 9, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$phi_rmse,
         type = "b", pch = 19,
         xlab = "Replicate",
         ylab = "RMSE",
         main = "phi RMSE")
    grid()
    plot(diag_df$rep_id, diag_df$phi_cor,
         type = "b", pch = 19,
         ylim = c(min(c(0, diag_df$phi_cor), na.rm = TRUE), 1),
         xlab = "Replicate",
         ylab = "Correlation",
         main = "cor(phi posterior mean, truth)")
    abline(h = 1, lty = 2)
    grid()
}

plot_s1_r_recovery <- function(r_df, file, reps_to_show = NULL) {
    if (is.null(reps_to_show)) {
        reps_to_show <- sort(unique(r_df$rep_id))[1:min(4, length(unique(r_df$rep_id)))]
    }
    d_all <- r_df[r_df$rep_id %in% reps_to_show, ]
    n_panels <- length(reps_to_show)
    nr <- ceiling(sqrt(n_panels))
    nc <- ceiling(n_panels / nr)

    .open_plot_device(file, width = 4.8 * nc, height = 4.0 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (rep_id in reps_to_show) {
        d <- d_all[d_all$rep_id == rep_id, ]
        ylim <- range(c(d$q025, d$q975, d$truth), finite = TRUE)
        pad <- diff(ylim) * 0.08
        if (!is.finite(pad) || pad == 0) pad <- 0.1
        ylim <- ylim + c(-pad, pad)
        plot(d$region, d$mean,
             ylim = ylim,
             pch = 19,
             xlab = "Region",
             ylab = "r posterior mean and 95% interval",
             main = paste("rep", rep_id))
        segments(d$region, d$q025, d$region, d$q975)
        lines(d$region, d$truth, lty = 2, lwd = 2)
        grid()
    }
}

plot_s1_r_metrics <- function(diag_df, file) {
    .open_plot_device(file, width = 9, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$r_rmse,
         type = "b", pch = 19,
         xlab = "Replicate",
         ylab = "RMSE",
         main = "r RMSE")
    grid()
    plot(diag_df$rep_id, diag_df$r_mean_bias,
         type = "b", pch = 19,
         xlab = "Replicate",
         ylab = "Mean r bias",
         main = "mean posterior r minus truth")
    abline(h = 0, lty = 2, lwd = 2)
    grid()
}

plot_s1_mcmc_diagnostics <- function(diag_df, file) {
    .open_plot_device(file, width = 10, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device() }, add = TRUE)

    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$phi_accept_rate,
         type = "b", pch = 19,
         ylim = c(0, 1),
         xlab = "Replicate",
         ylab = "Acceptance rate",
         main = "phi acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$r_accept_rate_mean,
         type = "b", pch = 19,
         ylim = c(0, 1),
         xlab = "Replicate",
         ylab = "Acceptance rate",
         main = "r acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$beta_mean_n_reject,
         type = "b", pch = 19,
         xlab = "Replicate",
         ylab = "Mean beta ESS rejections",
         main = "beta update diagnostic")
    grid()
}

make_s1_posterior_performance_plots <- function(perf,
                                                fig_dir,
                                                plot_format = c("pdf", "png"),
                                                reps_to_show = NULL) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    plot_format <- match.arg(plot_format)
    ext <- plot_format

    files <- list(
        beta_recovery = file.path(fig_dir, paste0("s1_beta_recovery_intervals.", ext)),
        beta_bias = file.path(fig_dir, paste0("s1_beta_bias_boxplot.", ext)),
        phi_scatter = file.path(fig_dir, paste0("s1_phi_truth_vs_posterior_mean.", ext)),
        phi_metrics = file.path(fig_dir, paste0("s1_phi_metrics.", ext)),
        r_recovery = file.path(fig_dir, paste0("s1_r_recovery_intervals.", ext)),
        r_metrics = file.path(fig_dir, paste0("s1_r_metrics.", ext)),
        mcmc_diagnostics = file.path(fig_dir, paste0("s1_mcmc_diagnostics.", ext))
    )

    plot_s1_beta_recovery(perf$beta, files$beta_recovery)
    plot_s1_beta_bias(perf$beta, files$beta_bias)
    plot_s1_phi_scatter(perf$phi, files$phi_scatter, reps_to_show = reps_to_show)
    plot_s1_phi_metrics(perf$diagnostics, files$phi_metrics)
    plot_s1_r_recovery(perf$r, files$r_recovery, reps_to_show = reps_to_show)
    plot_s1_r_metrics(perf$diagnostics, files$r_metrics)
    plot_s1_mcmc_diagnostics(perf$diagnostics, files$mcmc_diagnostics)

    invisible(files)
}

summarise_s1_performance_in_console <- function(perf) {
    beta <- perf$beta
    diag <- perf$diagnostics

    cat("\n================ Scenario 1 posterior performance ================\n")
    cat("Lambda fixed in DGP range: ",
        paste(range(diag$lambda_truth_min, diag$lambda_truth_max), collapse = " to "), "\n", sep = "")
    cat("Lambda fixed in fit: ", all(diag$lambda_fixed_in_fit),
        "; fit value range: ", paste(range(diag$lambda_fixed_value), collapse = " to "), "\n", sep = "")

    cat("\nBeta recovery by parameter:\n")
    beta_summary <- aggregate(cbind(bias, covered) ~ parameter, data = beta, FUN = mean)
    beta_rmse <- aggregate(bias ~ parameter, data = beta, FUN = function(x) sqrt(mean(x^2)))
    names(beta_rmse)[2] <- "rmse"
    beta_summary <- merge(beta_summary, beta_rmse, by = "parameter")
    print(beta_summary, row.names = FALSE)

    cat("\nPhi recovery across replicates:\n")
    cat(sprintf("  mean phi RMSE      : %.4f\n", mean(diag$phi_rmse, na.rm = TRUE)))
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
## Main one-stop driver
## ----------------------------------------------------------------------------
run_scenario1_posterior_performance <- function(root = ".",
                                                reps = 1:10,
                                                scenario_id = "S1_BASELINE_ONLY",
                                                data_dir = "data_revised",
                                                output_dir = "output_s1_baseline_only",
                                                analysis_dir = "analysis_s1_baseline_only",
                                                s1_core_file = "s1_baseline_only_clean.R",
                                                n_iter = 20000L,
                                                n_burnin = 10000L,
                                                n_thin = 5L,
                                                overwrite_data = TRUE,
                                                overwrite_fit = TRUE,
                                                plot_format = c("pdf", "png"),
                                                reps_to_show = NULL,
                                                verbose = 1000L,
                                                ...) {
    plot_format <- match.arg(plot_format)

    core_path <- file.path(root, s1_core_file)
    .require_file_s1run(core_path)
    source(core_path, local = .GlobalEnv)

    source_s1_baseline_only(root = root, verbose = TRUE)

    cat("\n=== Simulating Scenario 1 baseline-only data ===\n")
    manifest <- simulate_s1_baseline_only_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        overwrite_existing = overwrite_data,
        verbose = TRUE,
        ...
    )

    cat("\n=== Fitting Scenario 1 baseline-only model ===\n")
    summary_all <- fit_s1_baseline_batch(
        reps = reps,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root,
        settings_override = list(
            n_iter = as.integer(n_iter),
            n_burnin = as.integer(n_burnin),
            n_thin = as.integer(n_thin)
        ),
        verbose = verbose,
        overwrite_existing = overwrite_fit
    )

    cat("\n=== Collecting posterior performance summaries ===\n")
    perf <- collect_s1_posterior_performance(
        reps = reps,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root
    )

    analysis_out_dir <- file.path(root, analysis_dir, scenario_id)
    fig_dir <- file.path(analysis_out_dir, "figures")
    tab_dir <- file.path(analysis_out_dir, "tables")
    dir.create(analysis_out_dir, recursive = TRUE, showWarnings = FALSE)

    table_files <- save_s1_performance_tables(perf, tab_dir)
    plot_files <- make_s1_posterior_performance_plots(
        perf = perf,
        fig_dir = fig_dir,
        plot_format = plot_format,
        reps_to_show = reps_to_show
    )

    saveRDS(
        list(
            manifest = manifest,
            summary_all = summary_all,
            posterior_performance = perf,
            table_files = table_files,
            plot_files = plot_files,
            settings = list(
                reps = reps,
                n_iter = n_iter,
                n_burnin = n_burnin,
                n_thin = n_thin,
                scenario_id = scenario_id,
                data_dir = data_dir,
                output_dir = output_dir,
                analysis_dir = analysis_dir
            )
        ),
        file = file.path(analysis_out_dir, "scenario1_posterior_performance_results.rds")
    )

    summarise_s1_performance_in_console(perf)

    cat("Saved tables to: ", tab_dir, "\n", sep = "")
    cat("Saved figures to: ", fig_dir, "\n", sep = "")

    invisible(list(
        manifest = manifest,
        summary_all = summary_all,
        posterior_performance = perf,
        table_files = table_files,
        plot_files = plot_files,
        analysis_out_dir = analysis_out_dir
    ))
}

## No automatic execution.  Source this file and call
## run_scenario1_posterior_performance() explicitly.
