## ============================================================================
## VERSION: S3_POSTERIOR_PERFORMANCE_V4_TABLES_FIGURES_2026_05_26
## run_scenario3_posterior_performance.R
##
## Posterior performance summaries for Scenario 3:
## dynamic residual risk with learned common gamma.
##
## Design goal
##   Mirror the Scenario 2 posterior-performance workflow as closely as possible,
##   but add gamma posterior summaries.  Unlike the earlier thin version, this
##   script reads the original data file recorded in fit$metadata$data_file, so
##   truth values do not need to be stored as fit$truth.
##
## Typical use from run_s3_dynamic_learned_gamma_T100.R:
##   perf_out <- summarize_scenario3_posterior_performance(
##       root = root_dir,
##       fit_dir = file.path(output_dir, scenario_id),
##       analysis_dir = file.path(analysis_dir, scenario_id),
##       fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\\\.rds$"
##   )
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s3run_require_file <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.s3run_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}


.s3run_open_plot_device <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.s3run_close_plot_device <- function() {
    grDevices::dev.off()
}

.s3run_safe_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(rep(NA_real_, length(probs)))
    }
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

.s3run_safe_cor <- function(x, y) {
    out <- suppressWarnings(stats::cor(as.vector(x), as.vector(y), use = "complete.obs"))
    if (!is.finite(out)) NA_real_ else out
}


.s3run_safe_range <- function(x, default = c(-1, 1)) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(default)
    }
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || r[1L] == r[2L]) {
        eps <- ifelse(is.finite(r[1L]) && abs(r[1L]) > 0, 0.05 * abs(r[1L]), 0.05)
        return(c(r[1L] - eps, r[2L] + eps))
    }
    r
}

.s3run_num_scalar <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) {
        return(default)
    }
    suppressWarnings(as.numeric(x[[1L]]))
}

.s3run_char_scalar <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0L) {
        return(default)
    }
    as.character(x[[1L]])
}

.s3run_parse_rep <- function(path) {
    base <- basename(path)
    m <- regmatches(base, regexpr("rep[0-9]+", base))
    if (length(m) == 0L || is.na(m)) {
        return(NA_integer_)
    }
    as.integer(sub("rep", "", m))
}

.s3run_get_fit_files <- function(fit_dir,
                                 pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$") {
    sort(list.files(fit_dir, pattern = pattern, full.names = TRUE))
}

.s3run_get_data_file <- function(fit) {
    out <- fit$metadata$data_file %||% fit$data_file %||% NULL
    if (is.null(out)) {
        return(NULL)
    }
    as.character(out[[1L]])
}

.s3run_get_lambda_draws <- function(samples) {
    lam <- samples$lambda_tilde
    if (is.null(lam)) {
        stop("fit$samples$lambda_tilde is missing.", call. = FALSE)
    }
    if (length(dim(lam)) != 3L) {
        stop("fit$samples$lambda_tilde must be a 3D array: draws x time x region.", call. = FALSE)
    }
    lam
}

.s3run_get_gamma_draws <- function(samples) {
    if (!is.null(samples$gamma_common)) {
        return(as.numeric(samples$gamma_common))
    }
    if (!is.null(samples$gamma)) {
        g <- samples$gamma
        if (is.matrix(g) || is.data.frame(g)) {
            return(rowMeans(as.matrix(g), na.rm = TRUE))
        }
        return(as.numeric(g))
    }
    rep(NA_real_, length(samples$beta0 %||% numeric(0L)))
}

.s3run_summarize_parameter <- function(draws, truth, scenario_id, rep_id, parameter) {
    q <- .s3run_safe_quantile(draws)
    post_mean <- mean(draws, na.rm = TRUE)
    data.frame(
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id),
        parameter = parameter,
        truth = as.numeric(truth),
        post_mean = post_mean,
        post_sd = stats::sd(draws, na.rm = TRUE),
        q025 = q[1],
        q50 = q[2],
        q975 = q[3],
        bias = post_mean - as.numeric(truth),
        covered = as.integer(q[1] <= truth && truth <= q[3]),
        covered_95 = as.integer(q[1] <= truth && truth <= q[3]),
        stringsAsFactors = FALSE
    )
}

collect_s3_posterior_performance <- function(fit_files,
                                             scenario_id = "S3_DYNAMIC_LEARNED_GAMMA_T100",
                                             root = ".") {
    beta_list <- list()
    phi_list <- list()
    r_list <- list()
    gamma_list <- list()
    lambda_list <- list()
    diag_list <- list()
    conf_list <- list()

    for (ff in fit_files) {
        message("Reading: ", ff)
        fit <- readRDS(ff)
        s <- fit$samples
        if (is.null(s)) {
            stop("Fit file has no $samples component: ", ff, call. = FALSE)
        }

        rep_id <- fit$metadata$rep_id %||% fit$rep_id %||% .s3run_parse_rep(ff)
        rep_id <- as.integer(rep_id[[1L]])

        dat_file <- .s3run_get_data_file(fit)
        if (is.null(dat_file)) {
            stop(
                "Cannot locate the data file for fit: ", ff,
                ". Expected fit$metadata$data_file or fit$data_file.",
                call. = FALSE
            )
        }
        if (!file.exists(dat_file)) {
            dat_file_alt <- file.path(root, dat_file)
            if (file.exists(dat_file_alt)) {
                dat_file <- dat_file_alt
            }
        }
        .s3run_require_file(dat_file)
        dat <- readRDS(dat_file)

        TT <- dat$TT %||% nrow(dat$y_coarse)
        n1 <- dat$n1 %||% ncol(dat$y_coarse)

        beta0_true <- fit$summary$beta0_true %||% dat$beta0_star_ident %||% dat$beta0_star
        beta_true <- c(
            beta1 = fit$summary$beta1_true %||% dat$beta_star_ident[1] %||% dat$beta_star[1],
            beta2 = fit$summary$beta2_true %||% dat$beta_star_ident[2] %||% dat$beta_star[2]
        )
        phi_true <- dat$phi_star_ident %||% dat$phi_star
        r_true <- dat$r_star %||% rep(NA_real_, n1)
        if (length(r_true) == 1L) {
            r_true <- rep(r_true, n1)
        }
        lambda_true <- dat$lambda_tilde_ident %||% dat$lambda_tilde
        gamma_true <- mean(dat$gamma_star %||% fit$summary$gamma_truth_mean %||% NA_real_, na.rm = TRUE)

        beta_draws <- cbind(
            beta0 = as.numeric(s$beta0),
            beta1 = as.numeric(s$beta[, 1]),
            beta2 = as.numeric(s$beta[, 2])
        )
        beta_true <- as.numeric(beta_true)
        beta_truth <- c(beta0 = beta0_true, beta1 = beta_true[1], beta2 = beta_true[2])
        beta_df <- do.call(rbind, lapply(seq_along(beta_truth), function(k) {
            .s3run_summarize_parameter(
                draws = beta_draws[, k],
                truth = beta_truth[[k]],
                scenario_id = scenario_id,
                rep_id = rep_id,
                parameter = names(beta_truth)[[k]]
            )
        }))
        beta_list[[sprintf("%02d", rep_id)]] <- beta_df

        gamma_draws <- .s3run_get_gamma_draws(s)
        gamma_df <- .s3run_summarize_parameter(
            draws = gamma_draws,
            truth = gamma_true,
            scenario_id = scenario_id,
            rep_id = rep_id,
            parameter = "gamma_common"
        )
        names(gamma_df)[names(gamma_df) == "post_mean"] <- "mean"
        names(gamma_df)[names(gamma_df) == "post_sd"] <- "sd"
        gamma_df$median <- gamma_df$q50
        gamma_df$gamma_accept_rate <- .s3run_num_scalar(fit$diagnostics$gamma_accept_rate)
        gamma_df$gamma_proposal_sd_final <- .s3run_num_scalar(fit$diagnostics$gamma_proposal_sd_final)
        gamma_list[[sprintf("%02d", rep_id)]] <- gamma_df

        phi_mean <- colMeans(s$phi, na.rm = TRUE)
        phi_q <- apply(s$phi, 2, .s3run_safe_quantile)
        phi_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = rep_id,
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
        phi_list[[sprintf("%02d", rep_id)]] <- phi_df

        r_mean <- colMeans(s$r, na.rm = TRUE)
        r_q <- apply(s$r, 2, .s3run_safe_quantile)
        r_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = rep_id,
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
        r_list[[sprintf("%02d", rep_id)]] <- r_df

        lambda_draws <- .s3run_get_lambda_draws(s)
        lambda_mean <- apply(lambda_draws, c(2, 3), mean, na.rm = TRUE)
        lambda_q025 <- apply(lambda_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
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
            rep_id = rep_id,
            TT = as.integer(TT),
            n1 = as.integer(n1),
            gamma_truth_mean = gamma_true,
            gamma_learned_mean = mean(gamma_draws, na.rm = TRUE),
            gamma_learned_sd = stats::sd(gamma_draws, na.rm = TRUE),
            lambda_truth_min = min(lambda_true, na.rm = TRUE),
            lambda_truth_max = max(lambda_true, na.rm = TRUE),
            lambda_post_mean_min = min(lambda_mean, na.rm = TRUE),
            lambda_post_mean_max = max(lambda_mean, na.rm = TRUE),
            lambda_rmse = sqrt(mean((lambda_mean_vec - lambda_true_vec)^2, na.rm = TRUE)),
            lambda_mae = mean(abs(lambda_mean_vec - lambda_true_vec), na.rm = TRUE),
            lambda_coverage_95 = mean(as.vector(lambda_q025) <= lambda_true_vec &
                                          lambda_true_vec <= as.vector(lambda_q975), na.rm = TRUE),
            log_lambda_rmse = sqrt(mean((as.vector(log_lambda_mean) - as.vector(log_lambda_true))^2, na.rm = TRUE)),
            log_lambda_mae = mean(abs(as.vector(log_lambda_mean) - as.vector(log_lambda_true)), na.rm = TRUE),
            log_lambda_coverage_95 = mean(as.vector(log_lambda_q025) <= as.vector(log_lambda_true) &
                                              as.vector(log_lambda_true) <= as.vector(log_lambda_q975), na.rm = TRUE),
            cor_lambda = .s3run_safe_cor(lambda_mean, lambda_true),
            cor_log_lambda = .s3run_safe_cor(log_lambda_mean, log_lambda_true),
            rmse_delta_log_lambda = sqrt(mean((as.vector(delta_mean) - as.vector(delta_true))^2, na.rm = TRUE)),
            mae_delta_log_lambda = mean(abs(as.vector(delta_mean) - as.vector(delta_true)), na.rm = TRUE),
            cor_delta_log_lambda = .s3run_safe_cor(delta_mean, delta_true),
            stringsAsFactors = FALSE
        )
        lambda_list[[sprintf("%02d", rep_id)]] <- lambda_df

        x1 <- as.vector(dat$x1)
        x2 <- as.vector(dat$x2)
        loglam <- as.vector(log_lambda_true)
        phi_mat <- matrix(phi_true, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)
        conf_list[[sprintf("%02d", rep_id)]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = rep_id,
            x2_mode = .s3run_char_scalar(dat$x2_mode),
            sd_x1 = stats::sd(x1),
            sd_x2 = stats::sd(x2),
            cor_x1_x2 = .s3run_safe_cor(x1, x2),
            cor_x2_loglambda = .s3run_safe_cor(x2, loglam),
            cor_x1_loglambda = .s3run_safe_cor(x1, loglam),
            cor_x2_phi = .s3run_safe_cor(x2, as.vector(phi_mat)),
            cor_x1_phi = .s3run_safe_cor(x1, as.vector(phi_mat)),
            sd_beta2_x2 = stats::sd(dat$beta_star[2] * x2),
            stringsAsFactors = FALSE
        )

        diag_list[[sprintf("%02d", rep_id)]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = rep_id,
            TT = as.integer(TT),
            n1 = as.integer(n1),
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            dynamic_lambda_in_truth = TRUE,
            lambda_sampled_in_fit = TRUE,
            gamma_fixed_in_fit = FALSE,
            gamma_learned_in_fit = TRUE,
            gamma_truth_mean = gamma_true,
            gamma_learned_mean = mean(gamma_draws, na.rm = TRUE),
            gamma_learned_sd = stats::sd(gamma_draws, na.rm = TRUE),
            gamma_bias = mean(gamma_draws, na.rm = TRUE) - gamma_true,
            gamma_covered = gamma_df$covered,
            gamma_accept_rate = .s3run_num_scalar(fit$diagnostics$gamma_accept_rate),
            lambda_truth_min = .s3run_num_scalar(lambda_df$lambda_truth_min),
            lambda_truth_max = .s3run_num_scalar(lambda_df$lambda_truth_max),
            lambda_post_mean_min = .s3run_num_scalar(lambda_df$lambda_post_mean_min),
            lambda_post_mean_max = .s3run_num_scalar(lambda_df$lambda_post_mean_max),
            phi_rmse = sqrt(mean((phi_mean - phi_true)^2, na.rm = TRUE)),
            phi_cor = .s3run_safe_cor(phi_mean, phi_true),
            r_rmse = sqrt(mean((r_mean - r_true)^2, na.rm = TRUE)),
            r_mean_bias = mean(r_mean - r_true, na.rm = TRUE),
            beta0_bias = .s3run_num_scalar(beta_df$bias[beta_df$parameter == "beta0"]),
            beta1_bias = .s3run_num_scalar(beta_df$bias[beta_df$parameter == "beta1"]),
            beta2_bias = .s3run_num_scalar(beta_df$bias[beta_df$parameter == "beta2"]),
            beta0_covered = .s3run_num_scalar(beta_df$covered[beta_df$parameter == "beta0"]),
            beta1_covered = .s3run_num_scalar(beta_df$covered[beta_df$parameter == "beta1"]),
            beta2_covered = .s3run_num_scalar(beta_df$covered[beta_df$parameter == "beta2"]),
            phi_coverage = mean(phi_df$covered, na.rm = TRUE),
            r_coverage = mean(r_df$covered, na.rm = TRUE),
            log_lambda_rmse = .s3run_num_scalar(lambda_df$log_lambda_rmse),
            log_lambda_coverage_95 = .s3run_num_scalar(lambda_df$log_lambda_coverage_95),
            cor_log_lambda = .s3run_num_scalar(lambda_df$cor_log_lambda),
            cor_delta_log_lambda = .s3run_num_scalar(lambda_df$cor_delta_log_lambda),
            phi_accept_rate = .s3run_num_scalar(fit$diagnostics$phi_accept_rate),
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            beta_mean_n_reject = .s3run_num_scalar(fit$diagnostics$beta_mean_n_reject),
            elapsed_sec = .s3run_num_scalar(fit$diagnostics$elapsed_sec),
            stringsAsFactors = FALSE
        )
    }

    list(
        beta = do.call(rbind, beta_list),
        phi = do.call(rbind, phi_list),
        r = do.call(rbind, r_list),
        gamma = do.call(rbind, gamma_list),
        lambda = do.call(rbind, lambda_list),
        diagnostics = do.call(rbind, diag_list),
        confounding = do.call(rbind, conf_list)
    )
}

save_s3_performance_tables <- function(perf, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    files <- list(
        beta = file.path(out_dir, "posterior_beta_summary.csv"),
        phi = file.path(out_dir, "posterior_phi_summary.csv"),
        r = file.path(out_dir, "posterior_r_summary.csv"),
        gamma = file.path(out_dir, "posterior_gamma_summary.csv"),
        lambda = file.path(out_dir, "posterior_lambda_path_recovery.csv"),
        latent = file.path(out_dir, "posterior_latent_risk_performance.csv"),
        diagnostics = file.path(out_dir, "posterior_performance_diagnostics.csv"),
        confounding = file.path(out_dir, "scenario3_confounding_diagnostics_from_performance.csv")
    )
    utils::write.csv(perf$beta, files$beta, row.names = FALSE)
    utils::write.csv(perf$phi, files$phi, row.names = FALSE)
    utils::write.csv(perf$r, files$r, row.names = FALSE)
    utils::write.csv(perf$gamma, files$gamma, row.names = FALSE)
    utils::write.csv(perf$lambda, files$lambda, row.names = FALSE)
    utils::write.csv(perf$lambda, files$latent, row.names = FALSE)
    utils::write.csv(perf$diagnostics, files$diagnostics, row.names = FALSE)
    utils::write.csv(perf$confounding, files$confounding, row.names = FALSE)
    invisible(files)
}

make_s3_posterior_performance_plots <- function(perf,
                                                fig_dir,
                                                plot_format = c("pdf", "png"),
                                                reps_to_show = NULL) {
    plot_format <- match.arg(plot_format)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    ext <- plot_format

    files <- list(
        beta_bias = file.path(fig_dir, paste0("s3_beta_bias_by_rep.", ext)),
        gamma_posterior = file.path(fig_dir, paste0("s3_gamma_posterior_by_rep.", ext)),
        lambda_recovery = file.path(fig_dir, paste0("s3_lambda_recovery_by_rep.", ext))
    )

    beta <- perf$beta
    gamma <- perf$gamma
    lambda <- perf$lambda

    if (!is.null(reps_to_show)) {
        reps_to_show <- as.integer(reps_to_show)
        beta <- beta[beta$rep_id %in% reps_to_show, , drop = FALSE]
        gamma <- gamma[gamma$rep_id %in% reps_to_show, , drop = FALSE]
        lambda <- lambda[lambda$rep_id %in% reps_to_show, , drop = FALSE]
    }

    .s3run_open_plot_device(files$beta_bias, width = 8, height = 5)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    if (nrow(beta) > 0L) {
        beta$parameter <- as.character(beta$parameter)
        y_lim <- .s3run_safe_range(beta$bias)
        graphics::plot(
            NA, NA,
            xlim = range(beta$rep_id, na.rm = TRUE),
            ylim = y_lim,
            xlab = "Replicate", ylab = "Posterior mean bias",
            main = "Scenario 3 beta bias by replicate"
        )
        graphics::abline(h = 0, lty = 2)
        p_levels <- unique(beta$parameter)
        pch_vals <- seq_along(p_levels)
        for (i in seq_along(p_levels)) {
            d <- beta[beta$parameter == p_levels[[i]], , drop = FALSE]
            graphics::points(d$rep_id, d$bias, type = "b", pch = pch_vals[[i]])
        }
        graphics::legend("topright", legend = p_levels, pch = pch_vals, bty = "n")
    }
    .s3run_close_plot_device()

    .s3run_open_plot_device(files$gamma_posterior, width = 8, height = 5)
    if (nrow(gamma) > 0L) {
        y_lim <- .s3run_safe_range(c(gamma$q025, gamma$q975, gamma$truth), default = c(0, 1))
        graphics::plot(
            gamma$rep_id, gamma$mean,
            ylim = y_lim, pch = 19, type = "b",
            xlab = "Replicate", ylab = "Gamma",
            main = "Scenario 3 learned common gamma"
        )
        graphics::arrows(gamma$rep_id, gamma$q025, gamma$rep_id, gamma$q975,
                         angle = 90, code = 3, length = 0.04)
        graphics::abline(h = unique(gamma$truth)[1L], lty = 2)
    }
    .s3run_close_plot_device()

    .s3run_open_plot_device(files$lambda_recovery, width = 8, height = 5)
    if (nrow(lambda) > 0L) {
        y <- lambda$log_lambda_rmse
        graphics::plot(
            lambda$rep_id, y, type = "b", pch = 19,
            xlab = "Replicate", ylab = "log-lambda RMSE",
            main = "Scenario 3 latent path recovery"
        )
    }
    .s3run_close_plot_device()

    invisible(files)
}

summarise_s3_performance_in_console <- function(perf) {
    beta <- perf$beta
    diag <- perf$diagnostics
    lambda <- perf$lambda
    gamma <- perf$gamma

    cat("\n================ Scenario 3 posterior performance ================\n")
    cat("Dynamic lambda in truth: ", all(diag$dynamic_lambda_in_truth), "\n", sep = "")
    cat("Lambda sampled in fit  : ", all(diag$lambda_sampled_in_fit), "\n", sep = "")
    cat("Gamma learned in fit   : ", all(diag$gamma_learned_in_fit), "\n", sep = "")
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

    cat("\nGamma learning across replicates:\n")
    cat(sprintf("  mean posterior gamma      : %.4f\n", mean(gamma$mean, na.rm = TRUE)))
    cat(sprintf("  mean gamma bias           : %.4f\n", mean(gamma$bias, na.rm = TRUE)))
    cat(sprintf("  gamma 95%% coverage        : %.4f\n", mean(gamma$covered, na.rm = TRUE)))
    cat(sprintf("  mean gamma acceptance     : %.4f\n", mean(gamma$gamma_accept_rate, na.rm = TRUE)))

    cat("\nLambda path recovery across replicates:\n")
    cat(sprintf("  mean log-lambda RMSE      : %.4f\n", mean(lambda$log_lambda_rmse, na.rm = TRUE)))
    cat(sprintf("  mean log-lambda coverage  : %.4f\n", mean(lambda$log_lambda_coverage_95, na.rm = TRUE)))
    cat(sprintf("  mean cor(log lambda)      : %.4f\n", mean(lambda$cor_log_lambda, na.rm = TRUE)))
    cat(sprintf("  mean cor(delta log lambda): %.4f\n", mean(lambda$cor_delta_log_lambda, na.rm = TRUE)))

    cat("\nMCMC diagnostics:\n")
    cat(sprintf("  mean phi acceptance  : %.4f\n", mean(diag$phi_accept_rate, na.rm = TRUE)))
    cat(sprintf("  mean r acceptance    : %.4f\n", mean(diag$r_accept_rate_mean, na.rm = TRUE)))
    cat(sprintf("  mean gamma acceptance: %.4f\n", mean(diag$gamma_accept_rate, na.rm = TRUE)))
    cat("====================================================================\n\n")

    invisible(beta_summary)
}

summarize_scenario3_posterior_performance <- function(root = ".",
                                                      fit_dir = "output_s3_dynamic_learned_gamma/S3_DYNAMIC_LEARNED_GAMMA_T100",
                                                      analysis_dir = "analysis_s3_dynamic_learned_gamma/S3_DYNAMIC_LEARNED_GAMMA_T100",
                                                      fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$",
                                                      make_plots = TRUE,
                                                      plot_format = c("pdf", "png"),
                                                      reps_to_show = NULL,
                                                      verbose = TRUE) {
    plot_format <- match.arg(plot_format)
    root <- .s3run_norm_path(root)
    fit_dir <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)

    tab_dir <- file.path(analysis_dir, "tables")
    fig_dir <- file.path(analysis_dir, "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    fit_files <- .s3run_get_fit_files(fit_dir, pattern = fit_pattern)
    if (length(fit_files) == 0L) {
        stop("No Scenario 3 fit files found in: ", fit_dir, call. = FALSE)
    }

    first_fit <- readRDS(fit_files[[1L]])
    scenario_id <- first_fit$metadata$scenario_id %||% "S3_DYNAMIC_LEARNED_GAMMA_T100"
    scenario_id <- as.character(scenario_id[[1L]])

    perf <- collect_s3_posterior_performance(
        fit_files = fit_files,
        scenario_id = scenario_id,
        root = root
    )
    table_files <- save_s3_performance_tables(perf, out_dir = tab_dir)
    plot_files <- list()
    if (isTRUE(make_plots)) {
        plot_files <- make_s3_posterior_performance_plots(
            perf = perf,
            fig_dir = fig_dir,
            plot_format = plot_format,
            reps_to_show = reps_to_show
        )
    }

    rds_file <- file.path(analysis_dir, "scenario3_posterior_performance_results.rds")
    saveRDS(perf, rds_file)

    if (isTRUE(verbose)) {
        summarise_s3_performance_in_console(perf)
        cat("Saved Scenario 3 posterior performance RDS: ", rds_file, "\n", sep = "")
        cat("Tables: ", tab_dir, "\n", sep = "")
        if (isTRUE(make_plots)) cat("Figures: ", fig_dir, "\n", sep = "")
    }

    invisible(list(
        perf = perf,
        beta_summary = perf$beta,
        phi_summary = perf$phi,
        r_summary = perf$r,
        gamma_summary = perf$gamma,
        lambda_summary = perf$lambda,
        latent_summary = perf$lambda,
        performance = perf$diagnostics,
        diagnostics = perf$diagnostics,
        confounding = perf$confounding,
        tables = table_files,
        plots = plot_files,
        rds_file = rds_file,
        analysis_dir = analysis_dir,
        tables_dir = tab_dir,
        figures_dir = fig_dir
    ))
}

## Backward-compatible alias used by some older notes.
run_scenario3_posterior_performance <- function(root = ".",
                                                reps = NULL,
                                                scenario_id = "S3_DYNAMIC_LEARNED_GAMMA_T100",
                                                data_scenario_id = NULL,
                                                data_dir = "data_revised",
                                                output_dir = "output_s3_dynamic_learned_gamma",
                                                analysis_dir = "analysis_s3_dynamic_learned_gamma",
                                                plot_format = c("pdf", "png"),
                                                make_plots = TRUE,
                                                reps_to_show = NULL,
                                                verbose = TRUE,
                                                ...) {
    plot_format <- match.arg(plot_format)
    fit_dir <- file.path(output_dir, scenario_id)
    out_dir <- file.path(analysis_dir, scenario_id)
    summarize_scenario3_posterior_performance(
        root = root,
        fit_dir = fit_dir,
        analysis_dir = out_dir,
        fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$",
        make_plots = make_plots,
        plot_format = plot_format,
        reps_to_show = reps_to_show,
        verbose = verbose,
        ...
    )
}

if (sys.nframe() == 0L) {
    summarize_scenario3_posterior_performance(root = ".")
}
