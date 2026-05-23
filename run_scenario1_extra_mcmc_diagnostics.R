## ============================================================================
## run_scenario1_extra_mcmc_diagnostics.R
##
## Extra MCMC convergence diagnostics and posterior distribution plots for
## Scenario 1 baseline-only MSSTNB fits.
##
## This script does NOT rerun simulation or MCMC. It reads saved files from:
##   output_s1_baseline_only/S1_BASELINE_ONLY/fit_S1_baseline_repXX.rds
##   data_revised/S1_BASELINE_ONLY/data_repXX.rds
##
## It creates:
##   1. trace plots for beta, tau_phi, log-likelihood, phi, and r
##   2. running-mean plots for beta and r_mean
##   3. ACF plots
##   4. posterior histograms/densities
##   5. scalar MCMC diagnostics table: approximate ESS, MCSE, split z-score
##
## Recommended use:
##   source("run_scenario1_extra_mcmc_diagnostics.R")
##   out <- run_scenario1_extra_mcmc_diagnostics(
##       root = ".",
##       reps = 1:10,
##       plot_format = "pdf"
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s1diag <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.open_plot_device_s1diag <- function(file, width = 8.5, height = 6.5, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height, onefile = TRUE)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.close_plot_device_s1diag <- function() {
    grDevices::dev.off()
}

.safe_true_s1diag <- function(dat, name, fallback = NA_real_) {
    if (!is.null(dat[[name]])) dat[[name]] else fallback
}

.get_s1_truths_diag <- function(dat) {
    beta0_true <- dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_
    beta_true <- dat$beta_star %||% c(NA_real_, NA_real_)
    phi_true <- dat$phi_star_ident %||% dat$phi_star %||% rep(NA_real_, dat$n1)
    r_true <- dat$r_star %||% rep(NA_real_, dat$n1)
    if (length(r_true) == 1L) {
        r_true <- rep(r_true, dat$n1)
    }
    list(
        beta0 = beta0_true,
        beta = beta_true,
        phi = phi_true,
        r = r_true
    )
}

.as_scalar_matrix_s1diag <- function(fit, dat,
                                     selected_phi = NULL,
                                     selected_r = NULL) {
    s <- fit$samples
    n_draw <- length(s$beta0)
    n1 <- ncol(s$phi)

    if (is.null(selected_phi)) {
        selected_phi <- unique(round(seq(1, n1, length.out = min(4L, n1))))
    }
    if (is.null(selected_r)) {
        selected_r <- unique(round(seq(1, n1, length.out = min(4L, n1))))
    }
    selected_phi <- selected_phi[selected_phi >= 1L & selected_phi <= n1]
    selected_r <- selected_r[selected_r >= 1L & selected_r <= n1]

    out <- data.frame(
        draw = seq_len(n_draw),
        beta0 = as.numeric(s$beta0),
        beta1 = as.numeric(s$beta[, 1]),
        beta2 = as.numeric(s$beta[, 2]),
        tau_phi = as.numeric(s$tau_phi),
        loglik = as.numeric(s$loglik),
        r_mean = rowMeans(s$r)
    )

    for (j in selected_phi) {
        out[[paste0("phi_", j)]] <- s$phi[, j]
    }
    for (j in selected_r) {
        out[[paste0("r_", j)]] <- s$r[, j]
    }

    out
}

.approx_ess_s1diag <- function(x, max_lag = NULL) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 20L || stats::sd(x) == 0) {
        return(NA_real_)
    }
    if (is.null(max_lag)) {
        max_lag <- min(1000L, floor(n / 2))
    }
    ac <- stats::acf(x, lag.max = max_lag, plot = FALSE, na.action = na.pass)$acf[-1]
    if (length(ac) == 0L) {
        return(NA_real_)
    }

    ## Initial positive sequence. This is intentionally simple and conservative.
    keep <- which(ac > 0)
    if (length(keep) == 0L || keep[1] != 1L) {
        tau <- 1
    } else {
        last <- 0L
        for (k in seq_along(ac)) {
            if (!is.finite(ac[k]) || ac[k] <= 0) break
            last <- k
        }
        tau <- 1 + 2 * sum(ac[seq_len(last)], na.rm = TRUE)
    }

    ess <- n / max(tau, 1)
    max(1, min(n, ess))
}

.split_z_s1diag <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 40L) {
        return(NA_real_)
    }

    n_first <- max(10L, floor(0.10 * n))
    n_last <- max(10L, floor(0.50 * n))
    x1 <- x[seq_len(n_first)]
    x2 <- x[(n - n_last + 1L):n]

    v1 <- stats::var(x1) / length(x1)
    v2 <- stats::var(x2) / length(x2)
    if (!is.finite(v1 + v2) || (v1 + v2) <= 0) {
        return(NA_real_)
    }

    (mean(x1) - mean(x2)) / sqrt(v1 + v2)
}

.scalar_diag_table_s1diag <- function(draw_df, rep_id, max_lag = 1000L) {
    vars <- setdiff(names(draw_df), "draw")
    pieces <- lapply(vars, function(v) {
        x <- draw_df[[v]]
        ess <- .approx_ess_s1diag(x, max_lag = max_lag)
        data.frame(
            rep_id = rep_id,
            parameter = v,
            n_draw = sum(is.finite(x)),
            mean = mean(x, na.rm = TRUE),
            sd = stats::sd(x, na.rm = TRUE),
            q025 = as.numeric(stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE)),
            q500 = as.numeric(stats::quantile(x, 0.500, na.rm = TRUE, names = FALSE)),
            q975 = as.numeric(stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)),
            ess_approx = ess,
            mcse_approx = stats::sd(x, na.rm = TRUE) / sqrt(ess),
            split_z_approx = .split_z_s1diag(x)
        )
    })
    do.call(rbind, pieces)
}

.plot_trace_panel_s1diag <- function(draw_df, truth = NULL, main_prefix = "") {
    vars <- c("beta0", "beta1", "beta2", "tau_phi", "r_mean", "loglik")
    vars <- vars[vars %in% names(draw_df)]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (v in vars) {
        graphics::plot(draw_df$draw, draw_df[[v]], type = "l",
                       xlab = "Stored draw", ylab = v,
                       main = paste("Trace:", v))
        if (!is.null(truth) && !is.null(truth[[v]]) && is.finite(truth[[v]])) {
            graphics::abline(h = truth[[v]], lty = 2, lwd = 1.5)
        }
    }
    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
}

.plot_running_mean_s1diag <- function(draw_df, truth = NULL, main_prefix = "") {
    vars <- c("beta0", "beta1", "beta2", "r_mean")
    vars <- vars[vars %in% names(draw_df)]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (v in vars) {
        x <- draw_df[[v]]
        rm <- cumsum(x) / seq_along(x)
        graphics::plot(draw_df$draw, rm, type = "l",
                       xlab = "Stored draw", ylab = paste("Running mean of", v),
                       main = paste("Running mean:", v))
        if (!is.null(truth) && !is.null(truth[[v]]) && is.finite(truth[[v]])) {
            graphics::abline(h = truth[[v]], lty = 2, lwd = 1.5)
        }
    }
    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
}

.plot_acf_panel_s1diag <- function(draw_df, main_prefix = "", lag_max = 60L) {
    vars <- c("beta0", "beta1", "beta2", "tau_phi", "r_mean", "loglik")
    vars <- vars[vars %in% names(draw_df)]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (v in vars) {
        stats::acf(draw_df[[v]], lag.max = lag_max, main = paste("ACF:", v), na.action = na.pass)
    }
    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
}

.plot_density_panel_s1diag <- function(draw_df, truth = NULL, main_prefix = "") {
    vars <- c("beta0", "beta1", "beta2", "tau_phi", "r_mean")
    vars <- vars[vars %in% names(draw_df)]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (v in vars) {
        x <- draw_df[[v]]
        graphics::hist(x, breaks = 35, freq = FALSE,
                       xlab = v, main = paste("Posterior:", v),
                       border = "gray70")
        if (length(unique(x[is.finite(x)])) > 5L) {
            graphics::lines(stats::density(x, na.rm = TRUE), lwd = 2)
        }
        if (!is.null(truth) && !is.null(truth[[v]]) && is.finite(truth[[v]])) {
            graphics::abline(v = truth[[v]], lty = 2, lwd = 1.5)
        }
    }
    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
}

.plot_selected_phi_r_trace_s1diag <- function(draw_df, dat, main_prefix = "") {
    phi_vars <- grep("^phi_[0-9]+$", names(draw_df), value = TRUE)
    r_vars <- grep("^r_[0-9]+$", names(draw_df), value = TRUE)
    vars <- c(phi_vars, r_vars)
    if (length(vars) == 0L) {
        return(invisible(FALSE))
    }

    truth <- .get_s1_truths_diag(dat)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    n_panel <- min(8L, length(vars))
    nr <- ceiling(n_panel / 2)
    graphics::par(mfrow = c(nr, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

    for (v in vars[seq_len(n_panel)]) {
        graphics::plot(draw_df$draw, draw_df[[v]], type = "l",
                       xlab = "Stored draw", ylab = v,
                       main = paste("Trace:", v))
        idx <- as.integer(sub("^[a-z]+_", "", v))
        if (startsWith(v, "phi_")) {
            tv <- truth$phi[idx]
        } else {
            tv <- truth$r[idx]
        }
        if (is.finite(tv)) {
            graphics::abline(h = tv, lty = 2, lwd = 1.5)
        }
    }

    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
    invisible(TRUE)
}

.plot_selected_phi_r_density_s1diag <- function(draw_df, dat, main_prefix = "") {
    phi_vars <- grep("^phi_[0-9]+$", names(draw_df), value = TRUE)
    r_vars <- grep("^r_[0-9]+$", names(draw_df), value = TRUE)
    vars <- c(phi_vars, r_vars)
    if (length(vars) == 0L) {
        return(invisible(FALSE))
    }

    truth <- .get_s1_truths_diag(dat)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    n_panel <- min(8L, length(vars))
    nr <- ceiling(n_panel / 2)
    graphics::par(mfrow = c(nr, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

    for (v in vars[seq_len(n_panel)]) {
        x <- draw_df[[v]]
        graphics::hist(x, breaks = 35, freq = FALSE,
                       xlab = v, main = paste("Posterior:", v),
                       border = "gray70")
        if (length(unique(x[is.finite(x)])) > 5L) {
            graphics::lines(stats::density(x, na.rm = TRUE), lwd = 2)
        }
        idx <- as.integer(sub("^[a-z]+_", "", v))
        if (startsWith(v, "phi_")) {
            tv <- truth$phi[idx]
        } else {
            tv <- truth$r[idx]
        }
        if (is.finite(tv)) {
            graphics::abline(v = tv, lty = 2, lwd = 1.5)
        }
    }

    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
    invisible(TRUE)
}

.plot_pooled_beta_density_s1diag <- function(all_draws, beta_truth = NULL, main_prefix = "") {
    vars <- c("beta0", "beta1", "beta2")
    vars <- vars[vars %in% names(all_draws)]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)

    graphics::par(mfrow = c(1, 3), mar = c(3.5, 3.6, 2.2, 1), oma = c(0, 0, 2, 0))
    for (v in vars) {
        graphics::plot(stats::density(all_draws[[v]], na.rm = TRUE),
                       main = paste("Pooled posterior:", v), xlab = v, lwd = 2)
        if (!is.null(beta_truth) && !is.null(beta_truth[[v]]) && is.finite(beta_truth[[v]])) {
            graphics::abline(v = beta_truth[[v]], lty = 2, lwd = 1.5)
        }
    }
    graphics::mtext(main_prefix, outer = TRUE, cex = 1.1, font = 2)
}

## ----------------------------------------------------------------------------
## Main user-facing function
## ----------------------------------------------------------------------------
run_scenario1_extra_mcmc_diagnostics <- function(root = ".",
                                                 reps = 1:10,
                                                 scenario_id = "S1_BASELINE_ONLY",
                                                 data_dir = "data_revised",
                                                 output_dir = "output_s1_baseline_only",
                                                 analysis_dir = "analysis_s1_baseline_only",
                                                 selected_phi_regions = NULL,
                                                 selected_r_regions = NULL,
                                                 acf_lag_max = 60L,
                                                 ess_lag_max = 1000L,
                                                 plot_format = c("pdf", "png"),
                                                 verbose = TRUE) {
    plot_format <- match.arg(plot_format)

    fig_dir <- file.path(root, analysis_dir, scenario_id, "figures_mcmc")
    tab_dir <- file.path(root, analysis_dir, scenario_id, "tables")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

    diag_tables <- list()
    all_beta_draws <- list()
    first_truth <- NULL

    trace_file <- file.path(fig_dir, paste0("s1_mcmc_traceplots.", plot_format))
    runmean_file <- file.path(fig_dir, paste0("s1_mcmc_running_means.", plot_format))
    acf_file <- file.path(fig_dir, paste0("s1_mcmc_acf_plots.", plot_format))
    density_file <- file.path(fig_dir, paste0("s1_posterior_histograms_main.", plot_format))
    selected_trace_file <- file.path(fig_dir, paste0("s1_selected_phi_r_traceplots.", plot_format))
    selected_density_file <- file.path(fig_dir, paste0("s1_selected_phi_r_posteriors.", plot_format))
    pooled_beta_file <- file.path(fig_dir, paste0("s1_pooled_beta_posteriors.", plot_format))

    ## Multi-page PDF works naturally. For PNG, each plot file is overwritten by
    ## repeated pages, so we create one PNG per replicate instead.
    use_multipage <- identical(plot_format, "pdf")

    if (use_multipage) {
        .open_plot_device_s1diag(trace_file, width = 8.5, height = 7.0)
    }

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
        fit_file <- file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds"))
        .require_file_s1diag(dat_file)
        .require_file_s1diag(fit_file)

        dat <- readRDS(dat_file)
        fit <- readRDS(fit_file)
        truths <- .get_s1_truths_diag(dat)

        if (is.null(first_truth)) {
            first_truth <- list(
                beta0 = truths$beta0,
                beta1 = truths$beta[1],
                beta2 = truths$beta[2]
            )
        }

        draw_df <- .as_scalar_matrix_s1diag(
            fit = fit,
            dat = dat,
            selected_phi = selected_phi_regions,
            selected_r = selected_r_regions
        )

        draw_df$rep_id <- rep_id
        all_beta_draws[[length(all_beta_draws) + 1L]] <- draw_df[, c("rep_id", "draw", "beta0", "beta1", "beta2")]

        truth_scalar <- list(
            beta0 = truths$beta0,
            beta1 = truths$beta[1],
            beta2 = truths$beta[2],
            r_mean = mean(truths$r)
        )

        diag_tables[[length(diag_tables) + 1L]] <-
            .scalar_diag_table_s1diag(draw_df[, setdiff(names(draw_df), "rep_id"), drop = FALSE],
                                      rep_id = rep_id,
                                      max_lag = ess_lag_max)

        main_prefix <- paste0("Scenario 1 baseline-only, replicate ", rr)

        if (use_multipage) {
            .plot_trace_panel_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
        } else {
            f <- file.path(fig_dir, paste0("s1_mcmc_traceplots_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 7.0)
            .plot_trace_panel_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
            .close_plot_device_s1diag()
        }

        if (isTRUE(verbose)) {
            message("Processed MCMC trace diagnostics for replicate ", rr)
        }
    }

    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## Running means
    if (use_multipage) {
        .open_plot_device_s1diag(runmean_file, width = 8.5, height = 6.5)
    }
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat <- readRDS(file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds")))
        fit <- readRDS(file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds")))
        truths <- .get_s1_truths_diag(dat)
        draw_df <- .as_scalar_matrix_s1diag(fit, dat, selected_phi_regions, selected_r_regions)
        truth_scalar <- list(beta0 = truths$beta0, beta1 = truths$beta[1],
                             beta2 = truths$beta[2], r_mean = mean(truths$r))
        main_prefix <- paste0("Scenario 1 running means, replicate ", rr)
        if (use_multipage) {
            .plot_running_mean_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
        } else {
            f <- file.path(fig_dir, paste0("s1_mcmc_running_means_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 6.5)
            .plot_running_mean_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
            .close_plot_device_s1diag()
        }
    }
    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## ACF plots
    if (use_multipage) {
        .open_plot_device_s1diag(acf_file, width = 8.5, height = 7.0)
    }
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat <- readRDS(file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds")))
        fit <- readRDS(file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds")))
        draw_df <- .as_scalar_matrix_s1diag(fit, dat, selected_phi_regions, selected_r_regions)
        main_prefix <- paste0("Scenario 1 ACF, replicate ", rr)
        if (use_multipage) {
            .plot_acf_panel_s1diag(draw_df, main_prefix = main_prefix, lag_max = acf_lag_max)
        } else {
            f <- file.path(fig_dir, paste0("s1_mcmc_acf_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 7.0)
            .plot_acf_panel_s1diag(draw_df, main_prefix = main_prefix, lag_max = acf_lag_max)
            .close_plot_device_s1diag()
        }
    }
    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## Main posterior histograms
    if (use_multipage) {
        .open_plot_device_s1diag(density_file, width = 8.5, height = 7.0)
    }
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat <- readRDS(file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds")))
        fit <- readRDS(file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds")))
        truths <- .get_s1_truths_diag(dat)
        draw_df <- .as_scalar_matrix_s1diag(fit, dat, selected_phi_regions, selected_r_regions)
        truth_scalar <- list(beta0 = truths$beta0, beta1 = truths$beta[1],
                             beta2 = truths$beta[2], r_mean = mean(truths$r))
        main_prefix <- paste0("Scenario 1 posterior histograms, replicate ", rr)
        if (use_multipage) {
            .plot_density_panel_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
        } else {
            f <- file.path(fig_dir, paste0("s1_posterior_histograms_main_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 7.0)
            .plot_density_panel_s1diag(draw_df, truth = truth_scalar, main_prefix = main_prefix)
            .close_plot_device_s1diag()
        }
    }
    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## Selected phi/r trace plots
    if (use_multipage) {
        .open_plot_device_s1diag(selected_trace_file, width = 8.5, height = 8.0)
    }
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat <- readRDS(file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds")))
        fit <- readRDS(file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds")))
        draw_df <- .as_scalar_matrix_s1diag(fit, dat, selected_phi_regions, selected_r_regions)
        main_prefix <- paste0("Selected phi/r traceplots, replicate ", rr)
        if (use_multipage) {
            .plot_selected_phi_r_trace_s1diag(draw_df, dat, main_prefix = main_prefix)
        } else {
            f <- file.path(fig_dir, paste0("s1_selected_phi_r_traceplots_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 8.0)
            .plot_selected_phi_r_trace_s1diag(draw_df, dat, main_prefix = main_prefix)
            .close_plot_device_s1diag()
        }
    }
    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## Selected phi/r posterior histograms
    if (use_multipage) {
        .open_plot_device_s1diag(selected_density_file, width = 8.5, height = 8.0)
    }
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat <- readRDS(file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds")))
        fit <- readRDS(file.path(root, output_dir, scenario_id, paste0("fit_S1_baseline_rep", rr, ".rds")))
        draw_df <- .as_scalar_matrix_s1diag(fit, dat, selected_phi_regions, selected_r_regions)
        main_prefix <- paste0("Selected phi/r posteriors, replicate ", rr)
        if (use_multipage) {
            .plot_selected_phi_r_density_s1diag(draw_df, dat, main_prefix = main_prefix)
        } else {
            f <- file.path(fig_dir, paste0("s1_selected_phi_r_posteriors_rep", rr, ".", plot_format))
            .open_plot_device_s1diag(f, width = 8.5, height = 8.0)
            .plot_selected_phi_r_density_s1diag(draw_df, dat, main_prefix = main_prefix)
            .close_plot_device_s1diag()
        }
    }
    if (use_multipage) {
        .close_plot_device_s1diag()
    }

    ## Pooled beta posterior density. This is diagnostic only because each
    ## replicate has a different simulated data set.
    all_beta <- do.call(rbind, all_beta_draws)
    .open_plot_device_s1diag(pooled_beta_file, width = 8.5, height = 3.5)
    .plot_pooled_beta_density_s1diag(all_beta, beta_truth = first_truth,
                                     main_prefix = "Pooled beta posterior draws across replicates")
    .close_plot_device_s1diag()

    diag_table <- do.call(rbind, diag_tables)
    diag_file <- file.path(tab_dir, "s1_mcmc_scalar_diagnostics.csv")
    utils::write.csv(diag_table, diag_file, row.names = FALSE)

    ## A compact flag table for quick screening.
    flag_table <- within(diag_table, {
        low_ess <- is.finite(ess_approx) & ess_approx < 100
        large_split_z <- is.finite(split_z_approx) & abs(split_z_approx) > 2
        high_mcse_ratio <- is.finite(mcse_approx) & is.finite(sd) & sd > 0 & (mcse_approx / sd) > 0.10
    })
    flag_file <- file.path(tab_dir, "s1_mcmc_scalar_diagnostic_flags.csv")
    utils::write.csv(flag_table, flag_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("")
        message("Saved MCMC diagnostic tables:")
        message("  ", diag_file)
        message("  ", flag_file)
        message("")
        message("Saved MCMC diagnostic figures under:")
        message("  ", fig_dir)
        message("")
        message("Rule-of-thumb screening:")
        message("  ESS < 100: inspect trace and ACF.")
        message("  |split z| > 2: possible nonstationarity or slow mixing.")
        message("  MCSE / posterior SD > 0.10: posterior mean may be too noisy.")
    }

    invisible(list(
        diagnostics = diag_table,
        flags = flag_table,
        files = list(
            traceplots = trace_file,
            running_means = runmean_file,
            acf = acf_file,
            posterior_histograms = density_file,
            selected_phi_r_traceplots = selected_trace_file,
            selected_phi_r_posteriors = selected_density_file,
            pooled_beta_posteriors = pooled_beta_file,
            scalar_diagnostics = diag_file,
            scalar_flags = flag_file
        )
    ))
}
