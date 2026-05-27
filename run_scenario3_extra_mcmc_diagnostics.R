## ============================================================================
## VERSION: S3_EXTRA_MCMC_DIAGNOSTICS_V3_TABLES_FIGURES_2026_05_26
## run_scenario3_extra_mcmc_diagnostics.R
##
## Extra MCMC diagnostics for Scenario 3 dynamic learned-gamma fits.
## This version is deliberately parallel to run_scenario2_extra_mcmc_diagnostics.R,
## with the additional common-gamma diagnostics needed for Scenario 3.
##
## Important implementation detail:
##   samples$gamma may be an n_draws by n1 matrix because the common gamma is
##   replicated across regions. For diagnostics we use samples$gamma_common when
##   available, and otherwise fall back to the first gamma column.
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s3mcmc <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.open_plot_device_s3mcmc <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.close_plot_device_s3mcmc <- function() {
    grDevices::dev.off()
}

.parse_rep_s3mcmc <- function(path) {
    base <- basename(path)
    m <- regmatches(base, regexpr("rep[0-9]+", base))
    if (length(m) == 0L || is.na(m)) return(NA_integer_)
    as.integer(sub("rep", "", m))
}

.get_fit_file_s3mcmc <- function(root, output_dir, scenario_id, rep_id) {
    file.path(
        root,
        output_dir,
        scenario_id,
        sprintf("fit_S3_dynamic_learned_gamma_rep%02d.rds", as.integer(rep_id))
    )
}

.load_s3_fits_for_mcmc <- function(root,
                                   output_dir = "output_s3_dynamic_learned_gamma",
                                   scenario_id = "S3_DYNAMIC_LEARNED_GAMMA_T100",
                                   reps = 1:20,
                                   fit_dir = NULL,
                                   fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$") {
    if (!is.null(fit_dir)) {
        fit_dir_full <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
        fit_files <- list.files(fit_dir_full, pattern = fit_pattern, full.names = TRUE)
        if (length(fit_files) == 0L) {
            stop("No Scenario 3 fit files found in: ", fit_dir_full, call. = FALSE)
        }
        rep_ids <- vapply(fit_files, .parse_rep_s3mcmc, integer(1L))
        if (!is.null(reps)) {
            keep <- rep_ids %in% as.integer(reps)
            fit_files <- fit_files[keep]
            rep_ids <- rep_ids[keep]
        }
        ord <- order(rep_ids)
        fit_files <- fit_files[ord]
        rep_ids <- rep_ids[ord]
        out <- lapply(seq_along(fit_files), function(i) {
            .require_file_s3mcmc(fit_files[[i]])
            list(rep_id = as.integer(rep_ids[[i]]), fit = readRDS(fit_files[[i]]), fit_file = fit_files[[i]])
        })
        names(out) <- sprintf("rep%02d", as.integer(rep_ids))
        return(out)
    }

    out <- lapply(reps, function(rep_id) {
        fit_file <- .get_fit_file_s3mcmc(root, output_dir, scenario_id, rep_id)
        .require_file_s3mcmc(fit_file)
        list(rep_id = as.integer(rep_id), fit = readRDS(fit_file), fit_file = fit_file)
    })
    names(out) <- sprintf("rep%02d", as.integer(reps))
    out
}

.get_analysis_root_s3mcmc <- function(root, analysis_dir, scenario_id = NULL) {
    if (grepl("^/", analysis_dir)) {
        base <- analysis_dir
    } else {
        base <- file.path(root, analysis_dir)
    }
    if (!is.null(scenario_id) && nzchar(scenario_id) && basename(base) != scenario_id) {
        base <- file.path(base, scenario_id)
    }
    base
}

.get_gamma_common_s3mcmc <- function(fit) {
    s <- fit$samples %||% fit$draws %||% list()
    if (!is.null(s$gamma_common)) {
        return(as.numeric(s$gamma_common))
    }
    if (!is.null(s$gamma)) {
        g <- s$gamma
        if (is.matrix(g) || is.data.frame(g)) {
            return(as.numeric(g[, 1L]))
        }
        if (is.array(g)) {
            return(as.numeric(g[, 1L]))
        }
        return(as.numeric(g))
    }
    numeric(0L)
}

.acf1_s3mcmc <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 3L) return(NA_real_)
    suppressWarnings(stats::cor(x[-length(x)], x[-1L]))
}

.ess1_s3mcmc <- function(x, max_lag = NULL) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 10L) return(NA_real_)
    x <- x - mean(x)
    v <- sum(x^2)
    if (!is.finite(v) || v <= 0) return(NA_real_)
    if (is.null(max_lag)) max_lag <- min(1000L, n - 1L)

    rho_sum <- 0
    for (lag in seq_len(max_lag)) {
        rho <- sum(x[seq_len(n - lag)] * x[(lag + 1L):n]) / v
        if (!is.finite(rho) || rho <= 0) break
        rho_sum <- rho_sum + rho
    }
    n / (1 + 2 * rho_sum)
}

.parameter_diag_s3mcmc <- function(x, parameter, rep_id, scenario_id, n_stored = NULL) {
    x <- as.numeric(x)
    data.frame(
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id),
        parameter = parameter,
        n_draws = sum(is.finite(x)),
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        acf1 = .acf1_s3mcmc(x),
        ess = .ess1_s3mcmc(x),
        stringsAsFactors = FALSE
    )
}

.plot_trace_matrix_s3mcmc <- function(series_list, file, main, ylab) {
    n <- length(series_list)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n / nr)
    .open_plot_device_s3mcmc(file, width = 4.2 * nc, height = 3.2 * nr)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    for (nm in names(series_list)) {
        z <- as.numeric(series_list[[nm]])
        plot(z, type = "l", xlab = "Saved draw", ylab = ylab, main = paste(main, nm))
        abline(h = mean(z, na.rm = TRUE), lty = 2)
        grid()
    }
}

plot_s3_mcmc_beta_traces <- function(fits, file) {
    reps <- names(fits)
    n <- length(fits)
    .open_plot_device_s3mcmc(file, width = 12, height = max(7, 2.4 * n))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(n, 3), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        beta_series <- list(
            beta0 = s$beta0,
            beta1 = if (!is.null(s$beta) && ncol(s$beta) >= 1L) s$beta[, 1L] else NA_real_,
            beta2 = if (!is.null(s$beta) && ncol(s$beta) >= 2L) s$beta[, 2L] else NA_real_
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

plot_s3_mcmc_gamma_traces <- function(fits, file) {
    series <- lapply(fits, function(x) .get_gamma_common_s3mcmc(x$fit))
    .plot_trace_matrix_s3mcmc(series, file = file, main = "common gamma", ylab = "gamma")
}

plot_s3_mcmc_loglik_traces <- function(fits, file) {
    series <- lapply(fits, function(x) x$fit$samples$loglik_nb %||% x$fit$samples$loglik)
    .plot_trace_matrix_s3mcmc(series, file = file, main = "log-likelihood", ylab = "log-likelihood")
}

plot_s3_mcmc_phi_traces <- function(fits, file, selected_regions = c(1L, 5L, 9L)) {
    reps <- names(fits)
    nr <- length(reps)
    nc <- length(selected_regions)
    .open_plot_device_s3mcmc(file, width = 4.2 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        for (j in selected_regions) {
            if (!is.null(s$phi) && j <= ncol(s$phi)) {
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

plot_s3_mcmc_r_traces <- function(fits, file, selected_regions = c(1L, 5L, 9L)) {
    reps <- names(fits)
    nr <- length(reps)
    nc <- length(selected_regions)
    .open_plot_device_s3mcmc(file, width = 4.2 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in reps) {
        s <- fits[[nm]]$fit$samples
        for (j in selected_regions) {
            if (!is.null(s$r) && j <= ncol(s$r)) {
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

plot_s3_mcmc_lambda_traces <- function(fits,
                                       file,
                                       selected_regions = c(1L, 5L, 9L),
                                       selected_times = c(25L, 50L, 75L)) {
    combos <- expand.grid(time = selected_times, region = selected_regions)
    reps <- names(fits)
    nr <- length(reps)
    nc <- nrow(combos)
    .open_plot_device_s3mcmc(file, width = 3.4 * nc, height = max(4, 2.8 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
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

plot_s3_mcmc_acf <- function(fits, file, max_reps = 4L) {
    fits_show <- fits[seq_len(min(length(fits), max_reps))]
    nr <- length(fits_show)
    nc <- 4L
    .open_plot_device_s3mcmc(file, width = 4.2 * nc, height = max(4, 3 * nr))
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(nr, nc), mar = c(3.5, 3.5, 2, 1))
    for (nm in names(fits_show)) {
        s <- fits_show[[nm]]$fit$samples
        stats::acf(s$beta0, main = paste(nm, "ACF beta0"), na.action = na.pass)
        stats::acf(s$beta[, 1L], main = paste(nm, "ACF beta1"), na.action = na.pass)
        stats::acf(s$beta[, 2L], main = paste(nm, "ACF beta2"), na.action = na.pass)
        stats::acf(.get_gamma_common_s3mcmc(fits_show[[nm]]$fit), main = paste(nm, "ACF gamma"), na.action = na.pass)
    }
}

plot_s3_mcmc_acceptance_summary <- function(diag_df, file) {
    .open_plot_device_s3mcmc(file, width = 12, height = 5)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); .close_plot_device_s3mcmc() }, add = TRUE)
    par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
    plot(diag_df$rep_id, diag_df$phi_accept_rate, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Acceptance rate", main = "phi acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$r_accept_rate_mean, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Acceptance rate", main = "r acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$gamma_accept_rate, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "Replicate", ylab = "Acceptance rate", main = "gamma acceptance")
    grid()
    plot(diag_df$rep_id, diag_df$beta_mean_n_reject, type = "b", pch = 19,
         xlab = "Replicate", ylab = "Mean beta ESS rejections", main = "beta update diagnostic")
    grid()
}

collect_s3_mcmc_diagnostics <- function(fits, scenario_id = "S3_DYNAMIC_LEARNED_GAMMA_T100") {
    out <- lapply(fits, function(x) {
        fit <- x$fit
        gamma_common <- .get_gamma_common_s3mcmc(fit)
        data.frame(
            scenario_id = fit$metadata$scenario_id %||% scenario_id,
            rep_id = x$rep_id,
            n_stored = fit$n_stored %||% length(gamma_common),
            phi_accept_rate = fit$diagnostics$phi_accept_rate %||% NA_real_,
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            gamma_accept_rate = fit$diagnostics$gamma_accept_rate %||% NA_real_,
            beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject %||% NA_real_,
            elapsed_sec = fit$diagnostics$elapsed_sec %||% NA_real_,
            loglik_nb_mean = mean(fit$samples$loglik_nb %||% NA_real_, na.rm = TRUE),
            gamma_mean = mean(gamma_common, na.rm = TRUE),
            gamma_sd = stats::sd(gamma_common, na.rm = TRUE),
            gamma_acf1 = .acf1_s3mcmc(gamma_common),
            gamma_ess = .ess1_s3mcmc(gamma_common),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

collect_s3_parameter_mcmc_diagnostics <- function(fits, scenario_id = "S3_DYNAMIC_LEARNED_GAMMA_T100") {
    rows <- list()
    for (nm in names(fits)) {
        fit <- fits[[nm]]$fit
        rep_id <- fits[[nm]]$rep_id
        s <- fit$samples
        rows[[length(rows) + 1L]] <- .parameter_diag_s3mcmc(s$beta0, "beta0", rep_id, scenario_id)
        if (!is.null(s$beta) && ncol(s$beta) >= 1L) {
            rows[[length(rows) + 1L]] <- .parameter_diag_s3mcmc(s$beta[, 1L], "beta1", rep_id, scenario_id)
        }
        if (!is.null(s$beta) && ncol(s$beta) >= 2L) {
            rows[[length(rows) + 1L]] <- .parameter_diag_s3mcmc(s$beta[, 2L], "beta2", rep_id, scenario_id)
        }
        rows[[length(rows) + 1L]] <- .parameter_diag_s3mcmc(.get_gamma_common_s3mcmc(fit), "gamma_common", rep_id, scenario_id)
    }
    do.call(rbind, rows)
}

run_scenario3_extra_mcmc_diagnostics <- function(root = ".",
                                                 reps = NULL,
                                                 scenario_id = NULL,
                                                 output_dir = "output_s3_dynamic_learned_gamma",
                                                 analysis_dir = "analysis_s3_dynamic_learned_gamma",
                                                 fit_dir = NULL,
                                                 fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$",
                                                 make_plots = TRUE,
                                                 plot_format = c("pdf", "png"),
                                                 verbose = TRUE,
                                                 reps_to_show = NULL,
                                                 selected_phi_regions = c(1L, 5L, 9L),
                                                 selected_r_regions = c(1L, 5L, 9L),
                                                 selected_lambda_regions = c(1L, 5L, 9L),
                                                 selected_lambda_times = c(25L, 50L, 75L)) {
    plot_format <- match.arg(plot_format)
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)

    all_fits <- .load_s3_fits_for_mcmc(
        root = root,
        output_dir = output_dir,
        scenario_id = scenario_id %||% "S3_DYNAMIC_LEARNED_GAMMA_T100",
        reps = reps,
        fit_dir = fit_dir,
        fit_pattern = fit_pattern
    )

    if (is.null(reps_to_show)) {
        reps_to_show <- vapply(all_fits, function(x) x$rep_id, integer(1L))
        reps_to_show <- reps_to_show[seq_len(min(length(reps_to_show), 4L))]
    }
    fits_show <- all_fits[sprintf("rep%02d", as.integer(reps_to_show))]
    fits_show <- fits_show[!vapply(fits_show, is.null, logical(1L))]

    analysis_root <- .get_analysis_root_s3mcmc(root, analysis_dir, scenario_id)
    fig_dir <- file.path(analysis_root, "figures_mcmc")
    tab_dir <- file.path(analysis_root, "tables")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

    diag_df <- collect_s3_mcmc_diagnostics(all_fits, scenario_id = scenario_id %||% "S3_DYNAMIC_LEARNED_GAMMA_T100")
    param_diag_df <- collect_s3_parameter_mcmc_diagnostics(all_fits, scenario_id = scenario_id %||% "S3_DYNAMIC_LEARNED_GAMMA_T100")

    diag_file <- file.path(tab_dir, "scenario3_extra_mcmc_diagnostics.csv")
    param_diag_file <- file.path(tab_dir, "scenario3_parameter_mcmc_diagnostics.csv")
    ## S2-style aliases for easier cross-scenario navigation.
    accept_file <- file.path(tab_dir, "s3_mcmc_acceptance_summary.csv")
    param_file_alias <- file.path(tab_dir, "s3_mcmc_extra_diagnostics.csv")
    utils::write.csv(diag_df, diag_file, row.names = FALSE)
    utils::write.csv(param_diag_df, param_diag_file, row.names = FALSE)
    utils::write.csv(diag_df, accept_file, row.names = FALSE)
    utils::write.csv(param_diag_df, param_file_alias, row.names = FALSE)

    plot_files <- list()
    if (isTRUE(make_plots) && length(fits_show) > 0L) {
        ext <- plot_format
        plot_files <- list(
            beta_traces = file.path(fig_dir, paste0("s3_mcmc_beta_traces.", ext)),
            gamma_traces = file.path(fig_dir, paste0("s3_mcmc_gamma_traces.", ext)),
            loglik_traces = file.path(fig_dir, paste0("s3_mcmc_loglik_traces.", ext)),
            phi_traces = file.path(fig_dir, paste0("s3_mcmc_phi_traces.", ext)),
            r_traces = file.path(fig_dir, paste0("s3_mcmc_r_traces.", ext)),
            lambda_traces = file.path(fig_dir, paste0("s3_mcmc_lambda_traces.", ext)),
            beta_gamma_acf = file.path(fig_dir, paste0("s3_mcmc_beta_gamma_acf.", ext)),
            acceptance_summary = file.path(fig_dir, paste0("s3_mcmc_acceptance_summary.", ext))
        )

        plot_s3_mcmc_beta_traces(fits_show, plot_files$beta_traces)
        plot_s3_mcmc_gamma_traces(fits_show, plot_files$gamma_traces)
        plot_s3_mcmc_loglik_traces(fits_show, plot_files$loglik_traces)
        plot_s3_mcmc_phi_traces(fits_show, plot_files$phi_traces, selected_regions = selected_phi_regions)
        plot_s3_mcmc_r_traces(fits_show, plot_files$r_traces, selected_regions = selected_r_regions)
        plot_s3_mcmc_lambda_traces(
            fits_show,
            plot_files$lambda_traces,
            selected_regions = selected_lambda_regions,
            selected_times = selected_lambda_times
        )
        plot_s3_mcmc_acf(fits_show, plot_files$beta_gamma_acf)
        plot_s3_mcmc_acceptance_summary(diag_df, plot_files$acceptance_summary)
    }

    if (isTRUE(verbose)) {
        cat("\nScenario 3 extra MCMC diagnostics written to:\n")
        cat("  Table  : ", diag_file, "\n", sep = "")
        cat("  Params : ", param_diag_file, "\n", sep = "")
        if (isTRUE(make_plots)) cat("  Figures: ", fig_dir, "\n", sep = "")
    }

    invisible(list(
        mcmc_diagnostics = param_diag_df,
        diagnostics = diag_df,
        table = diag_file,
        parameter_table = param_diag_file,
        acceptance_table = accept_file,
        extra_diagnostics_table = param_file_alias,
        plots = plot_files,
        figures_mcmc_dir = fig_dir
    ))
}

