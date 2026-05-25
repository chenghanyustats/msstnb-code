## ============================================================================
## run_scenario2_extra_mcmc_diagnostics.R
##
## Extra MCMC diagnostics for Scenario 2 dynamic fixed-gamma fits.
##
## Typical use:
##   source("run_scenario2_extra_mcmc_diagnostics.R")
##   diag_out <- run_scenario2_extra_mcmc_diagnostics(
##       root = ".",
##       reps = 1:20,
##       scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
##       output_dir = "output_s2_dynamic_fixed_gamma",
##       analysis_dir = "analysis_s2_dynamic_fixed_gamma",
##       plot_format = "png"
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s2mcmc <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.open_plot_device_s2mcmc <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.close_plot_device_s2mcmc <- function() {
    grDevices::dev.off()
}

.get_fit_file_s2mcmc <- function(root, output_dir, scenario_id, rep_id) {
    file.path(
        root,
        output_dir,
        scenario_id,
        sprintf("fit_S2_dynamic_fixed_gamma_rep%02d.rds", as.integer(rep_id))
    )
}

.load_s2_fits_for_mcmc <- function(root, output_dir, scenario_id, reps) {
    out <- lapply(reps, function(rep_id) {
        fit_file <- .get_fit_file_s2mcmc(root, output_dir, scenario_id, rep_id)
        .require_file_s2mcmc(fit_file)
        list(rep_id = as.integer(rep_id), fit = readRDS(fit_file), fit_file = fit_file)
    })
    names(out) <- sprintf("rep%02d", as.integer(reps))
    out
}

.plot_trace_matrix <- function(series_list, file, main, ylab) {
    n <- length(series_list)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    .open_plot_device_s2mcmc(file, width = 4.2 * nc, height = 3.2 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (nm in names(series_list)) {
        z <- as.numeric(series_list[[nm]])
        plot(z, type = "l", xlab = "Saved draw", ylab = ylab, main = paste(main, nm))
        abline(h = mean(z, na.rm = TRUE), lty = 2)
        grid()
    }
}

plot_s2_mcmc_beta_traces <- function(fits, file) {
    reps <- names(fits)
    n <- length(fits)
    .open_plot_device_s2mcmc(file, width = 12, height = max(7, 2.4 * n))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(n, 3), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        beta_series <- list(
            beta0 = s$beta0,
            beta1 = s$beta[, 1],
            beta2 = s$beta[, 2]
        )
        truths <- c(
            beta0 = fits[[nm]]$fit$summary$beta0_true %||% NA_real_,
            beta1 = fits[[nm]]$fit$summary$beta1_true %||% NA_real_,
            beta2 = fits[[nm]]$fit$summary$beta2_true %||% NA_real_
        )
        for (parname in names(beta_series)) {
            z <- as.numeric(beta_series[[parname]])
            plot(z, type = "l", xlab = "Saved draw", ylab = parname,
                 main = paste(nm, parname))
            abline(h = mean(z, na.rm = TRUE), lty = 2)
            if (is.finite(truths[parname])) {
                abline(h = truths[parname], lty = 3, lwd = 2)
            }
            grid()
        }
    }
}

plot_s2_mcmc_loglik_traces <- function(fits, file) {
    series <- lapply(fits, function(x) x$fit$samples$loglik_nb %||% x$fit$samples$loglik)
    .plot_trace_matrix(series, file = file, main = "log-likelihood", ylab = "log-likelihood")
}

plot_s2_mcmc_phi_traces <- function(fits, file, selected_regions = c(1L, 5L, 9L)) {
    reps <- names(fits)
    nr <- length(reps)
    nc <- length(selected_regions)
    .open_plot_device_s2mcmc(file, width = 4.2 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        for (j in selected_regions) {
            if (j <= ncol(s$phi)) {
                z <- s$phi[, j]
                plot(z, type = "l", xlab = "Saved draw", ylab = paste0("phi[", j, "]"),
                     main = paste(nm, "region", j))
                abline(h = mean(z, na.rm = TRUE), lty = 2)
                grid()
            } else {
                plot.new(); title(main = paste(nm, "region", j, "missing"))
            }
        }
    }
}

plot_s2_mcmc_r_traces <- function(fits, file, selected_regions = c(1L, 5L, 9L)) {
    reps <- names(fits)
    nr <- length(reps)
    nc <- length(selected_regions)
    .open_plot_device_s2mcmc(file, width = 4.2 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        for (j in selected_regions) {
            if (j <= ncol(s$r)) {
                z <- s$r[, j]
                plot(z, type = "l", xlab = "Saved draw", ylab = paste0("r[", j, "]"),
                     main = paste(nm, "region", j))
                abline(h = mean(z, na.rm = TRUE), lty = 2)
                grid()
            } else {
                plot.new(); title(main = paste(nm, "region", j, "missing"))
            }
        }
    }
}

plot_s2_mcmc_lambda_traces <- function(fits,
                                       file,
                                       selected_regions = c(1L, 5L, 9L),
                                       selected_times = c(25L, 50L, 75L)) {
    combos <- expand.grid(time = selected_times, region = selected_regions)
    reps <- names(fits)
    nr <- length(reps)
    nc <- nrow(combos)
    .open_plot_device_s2mcmc(file, width = 3.4 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        lam <- fits[[nm]]$fit$samples$lambda_tilde
        for (kk in seq_len(nrow(combos))) {
            tt <- combos$time[kk]
            jj <- combos$region[kk]
            if (length(dim(lam)) == 3L && tt <= dim(lam)[2] && jj <= dim(lam)[3]) {
                z <- lam[, tt, jj]
                plot(z, type = "l", xlab = "Saved draw", ylab = "lambda",
                     main = paste0(nm, " t", tt, " r", jj))
                abline(h = mean(z, na.rm = TRUE), lty = 2)
                grid()
            } else {
                plot.new(); title(main = paste0(nm, " t", tt, " r", jj, " missing"))
            }
        }
    }
}

plot_s2_mcmc_acf <- function(fits, file, max_reps = 4L) {
    fits_show <- fits[seq_len(min(length(fits), max_reps))]
    nr <- length(fits_show)
    nc <- 3L
    .open_plot_device_s2mcmc(file, width = 4.2 * nc, height = max(4, 3 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in names(fits_show)) {
        s <- fits_show[[nm]]$fit$samples
        stats::acf(s$beta0, main = paste(nm, "ACF beta0"), na.action = na.pass)
        stats::acf(s$beta[, 1], main = paste(nm, "ACF beta1"), na.action = na.pass)
        stats::acf(s$beta[, 2], main = paste(nm, "ACF beta2"), na.action = na.pass)
    }
}

plot_s2_mcmc_acceptance_summary <- function(diag_df, file) {
    .open_plot_device_s2mcmc(file, width = 10, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s2mcmc() }, add = TRUE)
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

collect_s2_mcmc_diagnostics <- function(fits) {
    out <- lapply(fits, function(x) {
        fit <- x$fit
        data.frame(
            scenario_id = fit$metadata$scenario_id %||% NA_character_,
            rep_id = x$rep_id,
            n_stored = fit$n_stored %||% nrow(fit$samples$beta),
            phi_accept_rate = fit$diagnostics$phi_accept_rate %||% NA_real_,
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject %||% NA_real_,
            elapsed_sec = fit$diagnostics$elapsed_sec %||% NA_real_,
            loglik_nb_mean = mean(fit$samples$loglik_nb %||% NA_real_, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

run_scenario2_extra_mcmc_diagnostics <- function(root = ".",
                                                 reps = 1:20,
                                                 scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
                                                 output_dir = "output_s2_dynamic_fixed_gamma",
                                                 analysis_dir = "analysis_s2_dynamic_fixed_gamma",
                                                 plot_format = c("pdf", "png"),
                                                 reps_to_show = NULL,
                                                 selected_phi_regions = c(1L, 5L, 9L),
                                                 selected_r_regions = c(1L, 5L, 9L),
                                                 selected_lambda_regions = c(1L, 5L, 9L),
                                                 selected_lambda_times = c(25L, 50L, 75L)) {
    plot_format <- match.arg(plot_format)
    if (is.null(reps_to_show)) {
        reps_to_show <- reps[seq_len(min(length(reps), 4L))]
    }

    all_fits <- .load_s2_fits_for_mcmc(root, output_dir, scenario_id, reps)
    fits_show <- all_fits[sprintf("rep%02d", as.integer(reps_to_show))]

    analysis_root <- file.path(root, analysis_dir, scenario_id)
    fig_dir <- file.path(analysis_root, "figures_mcmc")
    tab_dir <- file.path(analysis_root, "tables")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

    diag_df <- collect_s2_mcmc_diagnostics(all_fits)
    diag_file <- file.path(tab_dir, "scenario2_extra_mcmc_diagnostics.csv")
    utils::write.csv(diag_df, diag_file, row.names = FALSE)

    ext <- plot_format
    plot_files <- list(
        beta_traces = file.path(fig_dir, paste0("s2_mcmc_beta_traces.", ext)),
        loglik_traces = file.path(fig_dir, paste0("s2_mcmc_loglik_traces.", ext)),
        phi_traces = file.path(fig_dir, paste0("s2_mcmc_phi_traces.", ext)),
        r_traces = file.path(fig_dir, paste0("s2_mcmc_r_traces.", ext)),
        lambda_traces = file.path(fig_dir, paste0("s2_mcmc_lambda_traces.", ext)),
        beta_acf = file.path(fig_dir, paste0("s2_mcmc_beta_acf.", ext)),
        acceptance_summary = file.path(fig_dir, paste0("s2_mcmc_acceptance_summary.", ext))
    )

    plot_s2_mcmc_beta_traces(fits_show, plot_files$beta_traces)
    plot_s2_mcmc_loglik_traces(fits_show, plot_files$loglik_traces)
    plot_s2_mcmc_phi_traces(fits_show, plot_files$phi_traces, selected_regions = selected_phi_regions)
    plot_s2_mcmc_r_traces(fits_show, plot_files$r_traces, selected_regions = selected_r_regions)
    plot_s2_mcmc_lambda_traces(
        fits_show,
        plot_files$lambda_traces,
        selected_regions = selected_lambda_regions,
        selected_times = selected_lambda_times
    )
    plot_s2_mcmc_acf(fits_show, plot_files$beta_acf)
    plot_s2_mcmc_acceptance_summary(diag_df, plot_files$acceptance_summary)

    cat("\nScenario 2 extra MCMC diagnostics written to:\n")
    cat("  Figures: ", fig_dir, "\n", sep = "")
    cat("  Table  : ", diag_file, "\n", sep = "")

    invisible(list(
        diagnostics = diag_df,
        table = diag_file,
        plots = plot_files,
        figures_mcmc_dir = fig_dir
    ))
}

