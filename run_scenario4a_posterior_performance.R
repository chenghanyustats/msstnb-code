## ============================================================================
## VERSION: S4A_POSTERIOR_PERFORMANCE_V1_STATUS_ORACLE_2026_06_18
## run_scenario4a_posterior_performance.R
##
## Posterior-performance and failure-mode summaries for Scenario 4A:
## observation-level sparse-count stress test with fixed gamma.
##
## Design goal
##   Mirror the Scenario 3 posterior-performance workflow, but adapt the output
##   to the S4A stress-test logic:
##     1. Do NOT average failed numerical-instability fits as ordinary fits.
##     2. Classify sampled-lambda fits into stable / soft_warning /
##        numerical_instability using guard counts and recovery diagnostics.
##     3. Save stable-only and nonfailed summaries.
##     4. Save oracle-lambda and oracle-lambda+phi diagnostic summaries when
##        those diagnostic outputs are available.
##
## Typical use:
##   source("run_scenario4a_posterior_performance.R")
##
## Or programmatically:
##   out <- summarize_scenario4a_posterior_performance(root = ".")
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s4a_require_file <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.s4a_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s4a_num_scalar <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    suppressWarnings(as.numeric(x[[1L]]))
}

.s4a_char_scalar <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0L) return(default)
    as.character(x[[1L]])
}

.s4a_bool_scalar <- function(x, default = NA) {
    if (is.null(x) || length(x) == 0L) return(default)
    as.logical(x[[1L]])
}

.s4a_safe_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(rep(NA_real_, length(probs)))
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

.s4a_safe_cor <- function(x, y) {
    out <- suppressWarnings(stats::cor(as.vector(x), as.vector(y), use = "complete.obs"))
    if (!is.finite(out)) NA_real_ else out
}

.s4a_parse_rep <- function(path) {
    base <- basename(path)
    m <- regmatches(base, regexpr("rep[0-9]+", base))
    if (length(m) == 0L || is.na(m)) return(NA_integer_)
    as.integer(sub("rep", "", m))
}

.s4a_get_fit_files <- function(fit_dir,
                               pattern = "fit_S4A_sparse_counts_obs_fixed_gamma_rep[0-9]+\\.rds$") {
    sort(list.files(fit_dir, pattern = pattern, full.names = TRUE))
}

.s4a_get_data_file <- function(fit, root = ".") {
    out <- fit$metadata$data_file %||% fit$data_file %||% NULL
    if (is.null(out)) return(NULL)
    out <- as.character(out[[1L]])
    if (file.exists(out)) return(out)
    alt <- file.path(root, out)
    if (file.exists(alt)) return(alt)
    out
}

.s4a_get_lambda_draws <- function(samples) {
    lam <- samples$lambda_tilde
    if (is.null(lam)) stop("fit$samples$lambda_tilde is missing.", call. = FALSE)
    if (length(dim(lam)) != 3L) {
        stop("fit$samples$lambda_tilde must be a 3D array: draws x time x region.", call. = FALSE)
    }
    lam
}

.s4a_get_gamma_draws <- function(samples) {
    if (!is.null(samples$gamma_common)) return(as.numeric(samples$gamma_common))
    if (!is.null(samples$gamma)) {
        g <- samples$gamma
        if (is.matrix(g) || is.data.frame(g)) return(rowMeans(as.matrix(g), na.rm = TRUE))
        return(as.numeric(g))
    }
    rep(NA_real_, length(samples$beta0 %||% numeric(0L)))
}

.s4a_summarize_parameter <- function(draws, truth, scenario_id, rep_id, parameter) {
    draws <- as.numeric(draws)
    q <- .s4a_safe_quantile(draws)
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

.s4a_classify_fit_status <- function(d) {
    beta_guard <- d$s4a_beta_guard_count %||% 0
    kappa_guard <- d$s4a_kappa_guard_count %||% 0
    lambda_in_guard <- d$s4a_lambda_input_guard_count %||% 0
    lambda_out_guard <- d$s4a_lambda_output_guard_count %||% 0
    beta0_mean <- d$beta0_mean %||% NA_real_
    lambda_rmse <- d$lambda_rmse %||% NA_real_
    log_lambda_rmse <- d$log_lambda_rmse %||% NA_real_

    bad <- (!is.finite(beta0_mean)) ||
        (!is.finite(lambda_rmse)) ||
        (!is.finite(log_lambda_rmse)) ||
        abs(beta0_mean) > 30 ||
        lambda_rmse > 10 ||
        log_lambda_rmse > 5 ||
        (kappa_guard > 0) ||
        (lambda_in_guard > 0) ||
        (lambda_out_guard > 0)

    if (isTRUE(bad)) return("numerical_instability")

    warn <- (beta_guard > 0)
    if (isTRUE(warn)) return("soft_warning")

    "stable"
}

.s4a_add_fit_status <- function(diag_df) {
    if (is.null(diag_df) || nrow(diag_df) == 0L) return(diag_df)
    status <- character(nrow(diag_df))
    for (i in seq_len(nrow(diag_df))) {
        status[[i]] <- .s4a_classify_fit_status(diag_df[i, , drop = FALSE])
    }
    diag_df$fit_status <- status
    diag_df$nonfailed <- diag_df$fit_status %in% c("stable", "soft_warning")
    diag_df$stable_fit <- diag_df$fit_status == "stable"
    diag_df
}

.collect_s4a_sampled_lambda_performance <- function(fit_files,
                                                     scenario_id = "S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100",
                                                     root = ".") {
    beta_list <- list()
    phi_list <- list()
    r_list <- list()
    gamma_list <- list()
    lambda_list <- list()
    diag_list <- list()
    conf_list <- list()

    for (ff in fit_files) {
        message("Reading sampled-lambda fit: ", ff)
        fit <- readRDS(ff)
        s <- fit$samples
        if (is.null(s)) stop("Fit file has no $samples component: ", ff, call. = FALSE)

        rep_id <- fit$metadata$rep_id %||% fit$rep_id %||% .s4a_parse_rep(ff)
        rep_id <- as.integer(rep_id[[1L]])

        dat_file <- .s4a_get_data_file(fit, root = root)
        if (is.null(dat_file)) {
            stop("Cannot locate the data file for fit: ", ff, call. = FALSE)
        }
        .s4a_require_file(dat_file)
        dat <- readRDS(dat_file)

        TT <- dat$TT %||% nrow(dat$y_coarse)
        n1 <- dat$n1 %||% ncol(dat$y_coarse)

        beta0_true <- fit$summary$beta0_true %||% dat$beta0_star_ident %||% dat$beta0_star
        beta_truth <- c(
            beta0 = beta0_true,
            beta1 = fit$summary$beta1_true %||% dat$beta_star_ident[1] %||% dat$beta_star[1],
            beta2 = fit$summary$beta2_true %||% dat$beta_star_ident[2] %||% dat$beta_star[2]
        )

        phi_true <- dat$phi_star_ident %||% dat$phi_star
        r_true <- dat$r_star %||% rep(NA_real_, n1)
        if (length(r_true) == 1L) r_true <- rep(r_true, n1)
        lambda_true <- dat$lambda_tilde_ident %||% dat$lambda_tilde
        gamma_true <- mean(dat$gamma_star %||% fit$summary$gamma_truth_mean %||% fit$metadata$gamma_fixed_value %||% NA_real_, na.rm = TRUE)

        beta_draws <- cbind(
            beta0 = as.numeric(s$beta0),
            beta1 = as.numeric(s$beta[, 1]),
            beta2 = as.numeric(s$beta[, 2])
        )
        beta_df <- do.call(rbind, lapply(seq_along(beta_truth), function(k) {
            .s4a_summarize_parameter(
                draws = beta_draws[, k],
                truth = beta_truth[[k]],
                scenario_id = scenario_id,
                rep_id = rep_id,
                parameter = names(beta_truth)[[k]]
            )
        }))
        beta_list[[sprintf("%02d", rep_id)]] <- beta_df

        gamma_draws <- .s4a_get_gamma_draws(s)
        gamma_df <- .s4a_summarize_parameter(
            draws = gamma_draws,
            truth = gamma_true,
            scenario_id = scenario_id,
            rep_id = rep_id,
            parameter = "gamma_common"
        )
        names(gamma_df)[names(gamma_df) == "post_mean"] <- "mean"
        names(gamma_df)[names(gamma_df) == "post_sd"] <- "sd"
        gamma_df$median <- gamma_df$q50
        gamma_df$gamma_fixed_in_fit <- TRUE
        gamma_df$gamma_learned_in_fit <- FALSE
        gamma_df$gamma_accept_rate <- NA_real_
        gamma_df$gamma_proposal_sd_final <- NA_real_
        gamma_list[[sprintf("%02d", rep_id)]] <- gamma_df

        phi_mean <- colMeans(s$phi, na.rm = TRUE)
        phi_q <- apply(s$phi, 2, .s4a_safe_quantile)
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
        r_q <- apply(s$r, 2, .s4a_safe_quantile)
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

        lambda_draws <- .s4a_get_lambda_draws(s)
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
            gamma_fixed_mean = mean(gamma_draws, na.rm = TRUE),
            gamma_fixed_sd = stats::sd(gamma_draws, na.rm = TRUE),
            lambda_truth_min = min(lambda_true, na.rm = TRUE),
            lambda_truth_max = max(lambda_true, na.rm = TRUE),
            lambda_post_mean_min = min(lambda_mean, na.rm = TRUE),
            lambda_post_mean_max = max(lambda_mean, na.rm = TRUE),
            lambda_rmse = sqrt(mean((lambda_mean_vec - lambda_true_vec)^2, na.rm = TRUE)),
            lambda_mae = mean(abs(lambda_mean_vec - lambda_true_vec), na.rm = TRUE),
            lambda_coverage_95 = mean(as.vector(lambda_q025) <= lambda_true_vec & lambda_true_vec <= as.vector(lambda_q975), na.rm = TRUE),
            log_lambda_rmse = sqrt(mean((as.vector(log_lambda_mean) - as.vector(log_lambda_true))^2, na.rm = TRUE)),
            log_lambda_mae = mean(abs(as.vector(log_lambda_mean) - as.vector(log_lambda_true)), na.rm = TRUE),
            log_lambda_coverage_95 = mean(as.vector(log_lambda_q025) <= as.vector(log_lambda_true) & as.vector(log_lambda_true) <= as.vector(log_lambda_q975), na.rm = TRUE),
            cor_lambda = .s4a_safe_cor(lambda_mean, lambda_true),
            cor_log_lambda = .s4a_safe_cor(log_lambda_mean, log_lambda_true),
            rmse_delta_log_lambda = sqrt(mean((as.vector(delta_mean) - as.vector(delta_true))^2, na.rm = TRUE)),
            mae_delta_log_lambda = mean(abs(as.vector(delta_mean) - as.vector(delta_true)), na.rm = TRUE),
            cor_delta_log_lambda = .s4a_safe_cor(delta_mean, delta_true),
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
            x2_mode = .s4a_char_scalar(dat$x2_mode),
            sd_x1 = stats::sd(x1),
            sd_x2 = stats::sd(x2),
            cor_x1_x2 = .s4a_safe_cor(x1, x2),
            cor_x2_loglambda = .s4a_safe_cor(x2, loglam),
            cor_x1_loglambda = .s4a_safe_cor(x1, loglam),
            cor_x2_phi = .s4a_safe_cor(x2, as.vector(phi_mat)),
            cor_x1_phi = .s4a_safe_cor(x1, as.vector(phi_mat)),
            sd_beta2_x2 = stats::sd(dat$beta_star[2] * x2),
            stringsAsFactors = FALSE
        )

        beta0_row <- beta_df[beta_df$parameter == "beta0", , drop = FALSE]
        beta1_row <- beta_df[beta_df$parameter == "beta1", , drop = FALSE]
        beta2_row <- beta_df[beta_df$parameter == "beta2", , drop = FALSE]

        diag_df <- data.frame(
            scenario_id = scenario_id,
            rep_id = rep_id,
            TT = as.integer(TT),
            n1 = as.integer(n1),
            mean_count = mean(dat$y_coarse),
            zero_prop = mean(dat$y_coarse == 0),
            stress_type = dat$stress_type %||% NA_character_,
            sparse_beta0_shift = dat$sparse_beta0_shift %||% NA_real_,
            expected_count_multiplier = dat$expected_count_multiplier %||% NA_real_,
            reference_mean_count = dat$reference_mean_count %||% NA_real_,
            reference_zero_prop = dat$reference_zero_prop %||% NA_real_,
            observed_total_count = sum(dat$y_coarse),
            observed_max_count = max(dat$y_coarse),
            dynamic_lambda_in_truth = TRUE,
            lambda_sampled_in_fit = TRUE,
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_truth_mean = gamma_true,
            gamma_fixed_mean = mean(gamma_draws, na.rm = TRUE),
            gamma_fixed_sd = stats::sd(gamma_draws, na.rm = TRUE),
            lambda_truth_min = .s4a_num_scalar(lambda_df$lambda_truth_min),
            lambda_truth_max = .s4a_num_scalar(lambda_df$lambda_truth_max),
            lambda_post_mean_min = .s4a_num_scalar(lambda_df$lambda_post_mean_min),
            lambda_post_mean_max = .s4a_num_scalar(lambda_df$lambda_post_mean_max),
            beta0_true = .s4a_num_scalar(beta0_row$truth),
            beta0_mean = .s4a_num_scalar(beta0_row$post_mean),
            beta0_sd = .s4a_num_scalar(beta0_row$post_sd),
            beta0_q025 = .s4a_num_scalar(beta0_row$q025),
            beta0_q975 = .s4a_num_scalar(beta0_row$q975),
            beta0_bias = .s4a_num_scalar(beta0_row$bias),
            beta0_covered = .s4a_num_scalar(beta0_row$covered),
            beta1_true = .s4a_num_scalar(beta1_row$truth),
            beta1_mean = .s4a_num_scalar(beta1_row$post_mean),
            beta1_sd = .s4a_num_scalar(beta1_row$post_sd),
            beta1_q025 = .s4a_num_scalar(beta1_row$q025),
            beta1_q975 = .s4a_num_scalar(beta1_row$q975),
            beta1_bias = .s4a_num_scalar(beta1_row$bias),
            beta1_covered = .s4a_num_scalar(beta1_row$covered),
            beta2_true = .s4a_num_scalar(beta2_row$truth),
            beta2_mean = .s4a_num_scalar(beta2_row$post_mean),
            beta2_sd = .s4a_num_scalar(beta2_row$post_sd),
            beta2_q025 = .s4a_num_scalar(beta2_row$q025),
            beta2_q975 = .s4a_num_scalar(beta2_row$q975),
            beta2_bias = .s4a_num_scalar(beta2_row$bias),
            beta2_covered = .s4a_num_scalar(beta2_row$covered),
            phi_rmse = sqrt(mean((phi_mean - phi_true)^2, na.rm = TRUE)),
            phi_cor = .s4a_safe_cor(phi_mean, phi_true),
            phi_coverage = mean(phi_df$covered, na.rm = TRUE),
            r_true_mean = mean(r_true, na.rm = TRUE),
            r_mean = mean(r_mean, na.rm = TRUE),
            r_rmse = sqrt(mean((r_mean - r_true)^2, na.rm = TRUE)),
            r_mean_bias = mean(r_mean - r_true, na.rm = TRUE),
            r_coverage = mean(r_df$covered, na.rm = TRUE),
            lambda_rmse = .s4a_num_scalar(lambda_df$lambda_rmse),
            lambda_coverage_95 = .s4a_num_scalar(lambda_df$lambda_coverage_95),
            log_lambda_rmse = .s4a_num_scalar(lambda_df$log_lambda_rmse),
            log_lambda_coverage_95 = .s4a_num_scalar(lambda_df$log_lambda_coverage_95),
            cor_log_lambda = .s4a_num_scalar(lambda_df$cor_log_lambda),
            cor_delta_log_lambda = .s4a_num_scalar(lambda_df$cor_delta_log_lambda),
            phi_accept_rate = .s4a_num_scalar(fit$diagnostics$phi_accept_rate),
            r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate %||% NA_real_, na.rm = TRUE),
            beta_mean_n_reject = .s4a_num_scalar(fit$diagnostics$beta_mean_n_reject),
            elapsed_sec = .s4a_num_scalar(fit$diagnostics$elapsed_sec),
            s4a_beta_guard_count = .s4a_num_scalar(fit$diagnostics$s4a_beta_guard_count, 0),
            s4a_kappa_guard_count = .s4a_num_scalar(fit$diagnostics$s4a_kappa_guard_count, 0),
            s4a_lambda_input_guard_count = .s4a_num_scalar(fit$diagnostics$s4a_lambda_input_guard_count, 0),
            s4a_lambda_output_guard_count = .s4a_num_scalar(fit$diagnostics$s4a_lambda_output_guard_count, 0),
            stringsAsFactors = FALSE
        )
        diag_list[[sprintf("%02d", rep_id)]] <- diag_df
    }

    beta <- do.call(rbind, beta_list)
    phi <- do.call(rbind, phi_list)
    r <- do.call(rbind, r_list)
    gamma <- do.call(rbind, gamma_list)
    lambda <- do.call(rbind, lambda_list)
    diagnostics <- .s4a_add_fit_status(do.call(rbind, diag_list))

    ## Attach fit status to long tables.
    status_key <- diagnostics[, c("rep_id", "fit_status", "nonfailed", "stable_fit"), drop = FALSE]
    beta <- merge(beta, status_key, by = "rep_id", all.x = TRUE, sort = FALSE)
    phi <- merge(phi, status_key, by = "rep_id", all.x = TRUE, sort = FALSE)
    r <- merge(r, status_key, by = "rep_id", all.x = TRUE, sort = FALSE)
    gamma <- merge(gamma, status_key, by = "rep_id", all.x = TRUE, sort = FALSE)
    lambda <- merge(lambda, status_key, by = "rep_id", all.x = TRUE, sort = FALSE)

    list(
        beta = beta,
        phi = phi,
        r = r,
        gamma = gamma,
        lambda = lambda,
        diagnostics = diagnostics,
        confounding = do.call(rbind, conf_list)
    )
}

.s4a_parameter_aggregate <- function(beta_df, subset_name = "all") {
    if (is.null(beta_df) || nrow(beta_df) == 0L) return(data.frame())
    out <- do.call(rbind, lapply(split(beta_df, beta_df$parameter), function(d) {
        data.frame(
            subset = subset_name,
            parameter = unique(d$parameter),
            n_reps = length(unique(d$rep_id)),
            mean_bias = mean(d$bias, na.rm = TRUE),
            median_bias = stats::median(d$bias, na.rm = TRUE),
            rmse_bias = sqrt(mean(d$bias^2, na.rm = TRUE)),
            mean_abs_bias = mean(abs(d$bias), na.rm = TRUE),
            coverage = mean(d$covered, na.rm = TRUE),
            mean_post_sd = mean(d$post_sd, na.rm = TRUE),
            mean_post_mean = mean(d$post_mean, na.rm = TRUE),
            sd_post_mean_across_reps = stats::sd(d$post_mean, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

.s4a_diagnostics_aggregate <- function(diag_df, subset_name = "all") {
    if (is.null(diag_df) || nrow(diag_df) == 0L) return(data.frame())
    data.frame(
        subset = subset_name,
        n_reps = nrow(diag_df),
        mean_count_avg = mean(diag_df$mean_count, na.rm = TRUE),
        zero_prop_avg = mean(diag_df$zero_prop, na.rm = TRUE),
        beta1_mean_avg = mean(diag_df$beta1_mean, na.rm = TRUE),
        beta1_mean_sd = stats::sd(diag_df$beta1_mean, na.rm = TRUE),
        beta2_mean_avg = mean(diag_df$beta2_mean, na.rm = TRUE),
        beta2_mean_sd = stats::sd(diag_df$beta2_mean, na.rm = TRUE),
        beta0_coverage = mean(diag_df$beta0_covered, na.rm = TRUE),
        beta1_coverage = mean(diag_df$beta1_covered, na.rm = TRUE),
        beta2_coverage = mean(diag_df$beta2_covered, na.rm = TRUE),
        r_mean_avg = mean(diag_df$r_mean, na.rm = TRUE),
        phi_rmse_avg = mean(diag_df$phi_rmse, na.rm = TRUE),
        phi_cor_avg = mean(diag_df$phi_cor, na.rm = TRUE),
        lambda_rmse_avg = mean(diag_df$lambda_rmse, na.rm = TRUE),
        log_lambda_rmse_avg = mean(diag_df$log_lambda_rmse, na.rm = TRUE),
        cor_log_lambda_avg = mean(diag_df$cor_log_lambda, na.rm = TRUE),
        beta_guard_total = sum(diag_df$s4a_beta_guard_count, na.rm = TRUE),
        kappa_guard_total = sum(diag_df$s4a_kappa_guard_count, na.rm = TRUE),
        lambda_output_guard_total = sum(diag_df$s4a_lambda_output_guard_count, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

.s4a_save_performance_tables <- function(perf,
                                         out_dir,
                                         oracle_lambda_summary = NULL,
                                         oracle_lamphi_summary = NULL,
                                         s3_control_summary = NULL,
                                         easy_s4a_summary = NULL) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    files <- list(
        beta = file.path(out_dir, "posterior_beta_summary.csv"),
        phi = file.path(out_dir, "posterior_phi_summary.csv"),
        r = file.path(out_dir, "posterior_r_summary.csv"),
        gamma = file.path(out_dir, "posterior_gamma_summary.csv"),
        lambda = file.path(out_dir, "posterior_lambda_path_recovery.csv"),
        latent = file.path(out_dir, "posterior_latent_risk_performance.csv"),
        diagnostics = file.path(out_dir, "posterior_performance_diagnostics.csv"),
        confounding = file.path(out_dir, "scenario4a_confounding_diagnostics_from_performance.csv"),
        fit_status = file.path(out_dir, "scenario4a_fit_status_by_rep.csv"),
        fit_status_counts = file.path(out_dir, "scenario4a_fit_status_counts.csv"),
        beta_aggregate = file.path(out_dir, "scenario4a_beta_recovery_by_subset.csv"),
        performance_aggregate = file.path(out_dir, "scenario4a_performance_by_subset.csv"),
        diagnostic_ladder = file.path(out_dir, "scenario4a_diagnostic_ladder_summary.csv")
    )

    utils::write.csv(perf$beta, files$beta, row.names = FALSE)
    utils::write.csv(perf$phi, files$phi, row.names = FALSE)
    utils::write.csv(perf$r, files$r, row.names = FALSE)
    utils::write.csv(perf$gamma, files$gamma, row.names = FALSE)
    utils::write.csv(perf$lambda, files$lambda, row.names = FALSE)
    utils::write.csv(perf$lambda, files$latent, row.names = FALSE)
    utils::write.csv(perf$diagnostics, files$diagnostics, row.names = FALSE)
    utils::write.csv(perf$confounding, files$confounding, row.names = FALSE)

    status_tbl <- perf$diagnostics[, c(
        "scenario_id", "rep_id", "fit_status", "nonfailed", "stable_fit",
        "mean_count", "zero_prop", "beta0_mean", "beta1_mean", "beta2_mean",
        "r_mean", "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "s4a_beta_guard_count", "s4a_kappa_guard_count",
        "s4a_lambda_input_guard_count", "s4a_lambda_output_guard_count"
    ), drop = FALSE]
    utils::write.csv(status_tbl, files$fit_status, row.names = FALSE)

    status_counts <- as.data.frame(table(perf$diagnostics$fit_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("fit_status", "n_reps")
    status_counts$prop_reps <- status_counts$n_reps / sum(status_counts$n_reps)
    utils::write.csv(status_counts, files$fit_status_counts, row.names = FALSE)

    beta_all <- .s4a_parameter_aggregate(perf$beta, "all_sampled_lambda")
    beta_stable <- .s4a_parameter_aggregate(subset(perf$beta, fit_status == "stable"), "stable_only")
    beta_nonfailed <- .s4a_parameter_aggregate(subset(perf$beta, fit_status %in% c("stable", "soft_warning")), "stable_plus_soft_warning")
    beta_failed <- .s4a_parameter_aggregate(subset(perf$beta, fit_status == "numerical_instability"), "numerical_instability")
    beta_agg <- do.call(rbind, list(beta_all, beta_stable, beta_nonfailed, beta_failed))
    utils::write.csv(beta_agg, files$beta_aggregate, row.names = FALSE)

    diag_all <- .s4a_diagnostics_aggregate(perf$diagnostics, "all_sampled_lambda")
    diag_stable <- .s4a_diagnostics_aggregate(subset(perf$diagnostics, fit_status == "stable"), "stable_only")
    diag_nonfailed <- .s4a_diagnostics_aggregate(subset(perf$diagnostics, fit_status %in% c("stable", "soft_warning")), "stable_plus_soft_warning")
    diag_failed <- .s4a_diagnostics_aggregate(subset(perf$diagnostics, fit_status == "numerical_instability"), "numerical_instability")
    diag_agg <- do.call(rbind, list(diag_all, diag_stable, diag_nonfailed, diag_failed))
    utils::write.csv(diag_agg, files$performance_aggregate, row.names = FALSE)

    ladder <- list()
    ladder[["sampled_lambda"]] <- data.frame(
        diagnostic_step = "S4A sampled lambda, fixed gamma",
        n_reps = nrow(perf$diagnostics),
        n_stable = sum(perf$diagnostics$fit_status == "stable", na.rm = TRUE),
        n_soft_warning = sum(perf$diagnostics$fit_status == "soft_warning", na.rm = TRUE),
        n_numerical_instability = sum(perf$diagnostics$fit_status == "numerical_instability", na.rm = TRUE),
        beta1_mean_avg = mean(perf$diagnostics$beta1_mean[perf$diagnostics$nonfailed], na.rm = TRUE),
        beta2_mean_avg = mean(perf$diagnostics$beta2_mean[perf$diagnostics$nonfailed], na.rm = TRUE),
        beta2_mean_sd = stats::sd(perf$diagnostics$beta2_mean[perf$diagnostics$nonfailed], na.rm = TRUE),
        lambda_rmse_avg = mean(perf$diagnostics$lambda_rmse[perf$diagnostics$nonfailed], na.rm = TRUE),
        cor_log_lambda_avg = mean(perf$diagnostics$cor_log_lambda[perf$diagnostics$nonfailed], na.rm = TRUE),
        conclusion = "Main stress-test fit; classify failures rather than averaging all reps.",
        stringsAsFactors = FALSE
    )

    if (!is.null(s3_control_summary) && nrow(s3_control_summary) > 0L) {
        ladder[["s3_control"]] <- data.frame(
            diagnostic_step = "S3 control, fixed gamma stabilized",
            n_reps = nrow(s3_control_summary),
            n_stable = nrow(s3_control_summary),
            n_soft_warning = 0L,
            n_numerical_instability = 0L,
            beta1_mean_avg = mean(s3_control_summary$beta1_mean, na.rm = TRUE),
            beta2_mean_avg = mean(s3_control_summary$beta2_mean, na.rm = TRUE),
            beta2_mean_sd = stats::sd(s3_control_summary$beta2_mean, na.rm = TRUE),
            lambda_rmse_avg = mean(s3_control_summary$lambda_rmse, na.rm = TRUE),
            cor_log_lambda_avg = mean(s3_control_summary$cor_log_lambda, na.rm = TRUE),
            conclusion = "Same machinery stable in non-sparse reference setting.",
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(easy_s4a_summary) && nrow(easy_s4a_summary) > 0L) {
        ladder[["easy_s4a"]] <- data.frame(
            diagnostic_step = "Easier S4A shift -4.0, fixed gamma stabilized",
            n_reps = nrow(easy_s4a_summary),
            n_stable = NA_integer_,
            n_soft_warning = NA_integer_,
            n_numerical_instability = NA_integer_,
            beta1_mean_avg = mean(easy_s4a_summary$beta1_mean, na.rm = TRUE),
            beta2_mean_avg = mean(easy_s4a_summary$beta2_mean, na.rm = TRUE),
            beta2_mean_sd = stats::sd(easy_s4a_summary$beta2_mean, na.rm = TRUE),
            lambda_rmse_avg = mean(easy_s4a_summary$lambda_rmse[is.finite(easy_s4a_summary$lambda_rmse) & easy_s4a_summary$lambda_rmse < 10], na.rm = TRUE),
            cor_log_lambda_avg = mean(easy_s4a_summary$cor_log_lambda, na.rm = TRUE),
            conclusion = "Same problem reps can remain unstable even with slightly easier sparsity.",
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(oracle_lambda_summary) && nrow(oracle_lambda_summary) > 0L) {
        files$oracle_lambda <- file.path(out_dir, "scenario4a_oracle_lambda_diagnostic_summary.csv")
        utils::write.csv(oracle_lambda_summary, files$oracle_lambda, row.names = FALSE)
        ladder[["oracle_lambda"]] <- data.frame(
            diagnostic_step = "Oracle lambda, estimate beta/phi/r",
            n_reps = nrow(oracle_lambda_summary),
            n_stable = NA_integer_,
            n_soft_warning = NA_integer_,
            n_numerical_instability = NA_integer_,
            beta1_mean_avg = mean(oracle_lambda_summary$beta1_mean, na.rm = TRUE),
            beta2_mean_avg = mean(oracle_lambda_summary$beta2_mean, na.rm = TRUE),
            beta2_mean_sd = stats::sd(oracle_lambda_summary$beta2_mean, na.rm = TRUE),
            lambda_rmse_avg = mean(oracle_lambda_summary$lambda_rmse, na.rm = TRUE),
            cor_log_lambda_avg = mean(oracle_lambda_summary$cor_log_lambda, na.rm = TRUE),
            conclusion = "Fixing lambda restores beta2, but beta0/phi/kappa can remain unstable.",
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(oracle_lamphi_summary) && nrow(oracle_lamphi_summary) > 0L) {
        files$oracle_lamphi <- file.path(out_dir, "scenario4a_oracle_lambda_phi_diagnostic_summary.csv")
        utils::write.csv(oracle_lamphi_summary, files$oracle_lamphi, row.names = FALSE)
        ladder[["oracle_lamphi"]] <- data.frame(
            diagnostic_step = "Oracle lambda + oracle phi, estimate beta/r",
            n_reps = nrow(oracle_lamphi_summary),
            n_stable = NA_integer_,
            n_soft_warning = NA_integer_,
            n_numerical_instability = NA_integer_,
            beta1_mean_avg = mean(oracle_lamphi_summary$beta1_mean, na.rm = TRUE),
            beta2_mean_avg = mean(oracle_lamphi_summary$beta2_mean, na.rm = TRUE),
            beta2_mean_sd = stats::sd(oracle_lamphi_summary$beta2_mean, na.rm = TRUE),
            lambda_rmse_avg = mean(oracle_lamphi_summary$lambda_rmse, na.rm = TRUE),
            cor_log_lambda_avg = mean(oracle_lamphi_summary$cor_log_lambda, na.rm = TRUE),
            conclusion = "Fixing both lambda and phi restores beta/r and nearly eliminates guards.",
            stringsAsFactors = FALSE
        )
    }

    ladder_df <- do.call(rbind, ladder)
    rownames(ladder_df) <- NULL
    utils::write.csv(ladder_df, files$diagnostic_ladder, row.names = FALSE)

    invisible(files)
}

.s4a_read_optional_csv <- function(path) {
    if (file.exists(path)) return(utils::read.csv(path))
    NULL
}

.s4a_open_plot_device <- function(file, width = 8, height = 6, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.s4a_close_plot_device <- function() {
    grDevices::dev.off()
}

.s4a_make_plots <- function(perf, fig_dir, plot_format = c("pdf", "png")) {
    plot_format <- match.arg(plot_format)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    ext <- plot_format
    files <- list(
        beta2_by_rep = file.path(fig_dir, paste0("s4a_beta2_mean_by_fit_status.", ext)),
        lambda_rmse_by_rep = file.path(fig_dir, paste0("s4a_log_lambda_rmse_by_fit_status.", ext)),
        guard_counts_by_rep = file.path(fig_dir, paste0("s4a_guard_counts_by_rep.", ext))
    )

    d <- perf$diagnostics
    pch_map <- c(stable = 19, soft_warning = 17, numerical_instability = 4)
    pch_vals <- pch_map[d$fit_status]
    pch_vals[is.na(pch_vals)] <- 1

    .s4a_open_plot_device(files$beta2_by_rep, width = 8, height = 5)
    graphics::plot(
        d$rep_id, d$beta2_mean,
        pch = pch_vals,
        xlab = "Replicate", ylab = "Posterior mean of beta2",
        main = "S4A beta2 recovery by fit status"
    )
    graphics::abline(h = unique(d$beta2_true)[1L], lty = 2)
    graphics::legend("topright", legend = names(pch_map), pch = pch_map, bty = "n")
    .s4a_close_plot_device()

    .s4a_open_plot_device(files$lambda_rmse_by_rep, width = 8, height = 5)
    y <- d$log_lambda_rmse
    y_plot <- pmin(y, 10)
    graphics::plot(
        d$rep_id, y_plot,
        pch = pch_vals,
        xlab = "Replicate", ylab = "log-lambda RMSE (capped at 10)",
        main = "S4A latent path recovery by fit status"
    )
    graphics::legend("topright", legend = names(pch_map), pch = pch_map, bty = "n")
    .s4a_close_plot_device()

    .s4a_open_plot_device(files$guard_counts_by_rep, width = 8, height = 5)
    total_guards <- d$s4a_beta_guard_count + d$s4a_kappa_guard_count + d$s4a_lambda_output_guard_count
    graphics::plot(
        d$rep_id, log10(total_guards + 1),
        type = "h", lwd = 3,
        xlab = "Replicate", ylab = "log10(total guards + 1)",
        main = "S4A numerical guard activity"
    )
    graphics::points(d$rep_id, log10(total_guards + 1), pch = pch_vals)
    .s4a_close_plot_device()

    invisible(files)
}

.s4a_console_summary <- function(perf,
                                 oracle_lambda_summary = NULL,
                                 oracle_lamphi_summary = NULL) {
    d <- perf$diagnostics
    beta <- perf$beta
    cat("\n================ Scenario 4A posterior performance ================\n")
    cat("Main fit: observation-level sparse counts, fixed gamma, sampled lambda\n")
    cat("Fit status counts:\n")
    print(as.data.frame(table(d$fit_status)), row.names = FALSE)

    cat("\nData sparsity:\n")
    cat(sprintf("  mean count average: %.4f\n", mean(d$mean_count, na.rm = TRUE)))
    cat(sprintf("  zero proportion avg: %.4f\n", mean(d$zero_prop, na.rm = TRUE)))

    cat("\nBeta recovery by subset:\n")
    beta_agg <- do.call(rbind, list(
        .s4a_parameter_aggregate(beta, "all_sampled_lambda"),
        .s4a_parameter_aggregate(subset(beta, fit_status == "stable"), "stable_only"),
        .s4a_parameter_aggregate(subset(beta, fit_status %in% c("stable", "soft_warning")), "stable_plus_soft_warning")
    ))
    print(beta_agg[, c("subset", "parameter", "n_reps", "mean_bias", "rmse_bias", "coverage", "mean_post_sd")], row.names = FALSE)

    cat("\nLatent path recovery by subset:\n")
    print(do.call(rbind, list(
        .s4a_diagnostics_aggregate(d, "all_sampled_lambda"),
        .s4a_diagnostics_aggregate(subset(d, fit_status == "stable"), "stable_only"),
        .s4a_diagnostics_aggregate(subset(d, fit_status %in% c("stable", "soft_warning")), "stable_plus_soft_warning")
    ))[, c("subset", "n_reps", "lambda_rmse_avg", "log_lambda_rmse_avg", "cor_log_lambda_avg", "beta_guard_total", "kappa_guard_total", "lambda_output_guard_total")], row.names = FALSE)

    if (!is.null(oracle_lambda_summary)) {
        cat("\nOracle lambda diagnostic available for reps: ", paste(oracle_lambda_summary$rep_id, collapse = ", "), "\n", sep = "")
    }
    if (!is.null(oracle_lamphi_summary)) {
        cat("Oracle lambda+phi diagnostic available for reps: ", paste(oracle_lamphi_summary$rep_id, collapse = ", "), "\n", sep = "")
    }
    cat("====================================================================\n\n")
}

summarize_scenario4a_posterior_performance <- function(root = ".",
                                                       fit_dir = "output_s4a_sparse_counts_fixed_gamma/S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100",
                                                       analysis_dir = "analysis_s4a_sparse_counts/S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100",
                                                       fit_pattern = "fit_S4A_sparse_counts_obs_fixed_gamma_rep[0-9]+\\.rds$",
                                                       oracle_lambda_summary_file = "output_s4a_oracle_lambda_fixed_gamma/S4A_ORACLE_LAMBDA_FIXED_GAMMA_DIAGNOSTIC_T100/summary_S4A_oracle_lambda_fixed_gamma_all_reps.csv",
                                                       oracle_lamphi_summary_file = "output_s4a_oracle_lambda_phi_fixed_gamma/S4A_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100/summary_S4A_oracle_lambda_phi_fixed_gamma_all_reps.csv",
                                                       s3_control_summary_file = "output_s3_control_fixed_gamma_stabilized/S3_CONTROL_FIXED_GAMMA_STABILIZED_T100/summary_S3_control_fixed_gamma_stabilized_all_reps.csv",
                                                       easy_s4a_summary_file = "output_s4a_easy_sparse_counts_fixed_gamma/S4A_EASY_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100/summary_S4A_easy_sparse_counts_obs_fixed_gamma_all_reps.csv",
                                                       make_plots = TRUE,
                                                       plot_format = c("pdf", "png"),
                                                       verbose = TRUE) {
    plot_format <- match.arg(plot_format)
    root <- .s4a_norm_path(root)
    fit_dir <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)

    tab_dir <- file.path(analysis_dir, "tables")
    fig_dir <- file.path(analysis_dir, "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    fit_files <- .s4a_get_fit_files(fit_dir, pattern = fit_pattern)
    if (length(fit_files) == 0L) {
        stop("No Scenario 4A sampled-lambda fit files found in: ", fit_dir, call. = FALSE)
    }

    first_fit <- readRDS(fit_files[[1L]])
    scenario_id <- first_fit$metadata$scenario_id %||% "S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100"
    scenario_id <- as.character(scenario_id[[1L]])

    perf <- .collect_s4a_sampled_lambda_performance(
        fit_files = fit_files,
        scenario_id = scenario_id,
        root = root
    )

    oracle_lambda_summary <- .s4a_read_optional_csv(if (grepl("^/", oracle_lambda_summary_file)) oracle_lambda_summary_file else file.path(root, oracle_lambda_summary_file))
    oracle_lamphi_summary <- .s4a_read_optional_csv(if (grepl("^/", oracle_lamphi_summary_file)) oracle_lamphi_summary_file else file.path(root, oracle_lamphi_summary_file))
    s3_control_summary <- .s4a_read_optional_csv(if (grepl("^/", s3_control_summary_file)) s3_control_summary_file else file.path(root, s3_control_summary_file))
    easy_s4a_summary <- .s4a_read_optional_csv(if (grepl("^/", easy_s4a_summary_file)) easy_s4a_summary_file else file.path(root, easy_s4a_summary_file))

    table_files <- .s4a_save_performance_tables(
        perf = perf,
        out_dir = tab_dir,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        s3_control_summary = s3_control_summary,
        easy_s4a_summary = easy_s4a_summary
    )

    plot_files <- list()
    if (isTRUE(make_plots)) {
        plot_files <- .s4a_make_plots(perf, fig_dir = fig_dir, plot_format = plot_format)
    }

    rds_file <- file.path(analysis_dir, "scenario4a_posterior_performance_results.rds")
    saveRDS(list(
        perf = perf,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        s3_control_summary = s3_control_summary,
        easy_s4a_summary = easy_s4a_summary
    ), rds_file)

    if (isTRUE(verbose)) {
        .s4a_console_summary(
            perf = perf,
            oracle_lambda_summary = oracle_lambda_summary,
            oracle_lamphi_summary = oracle_lamphi_summary
        )
        cat("Saved Scenario 4A posterior performance RDS: ", rds_file, "\n", sep = "")
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
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        s3_control_summary = s3_control_summary,
        easy_s4a_summary = easy_s4a_summary,
        tables = table_files,
        plots = plot_files,
        rds_file = rds_file,
        analysis_dir = analysis_dir,
        tables_dir = tab_dir,
        figures_dir = fig_dir
    ))
}

run_scenario4a_posterior_performance <- function(root = ".",
                                                 reps = NULL,
                                                 scenario_id = "S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100",
                                                 output_dir = "output_s4a_sparse_counts_fixed_gamma",
                                                 analysis_dir = "analysis_s4a_sparse_counts",
                                                 plot_format = c("pdf", "png"),
                                                 make_plots = TRUE,
                                                 verbose = TRUE,
                                                 ...) {
    plot_format <- match.arg(plot_format)
    fit_dir <- file.path(output_dir, scenario_id)
    out_dir <- file.path(analysis_dir, scenario_id)
    summarize_scenario4a_posterior_performance(
        root = root,
        fit_dir = fit_dir,
        analysis_dir = out_dir,
        fit_pattern = "fit_S4A_sparse_counts_obs_fixed_gamma_rep[0-9]+\\.rds$",
        make_plots = make_plots,
        plot_format = plot_format,
        verbose = verbose,
        ...
    )
}

if (sys.nframe() == 0L) {
    summarize_scenario4a_posterior_performance(root = ".")
}
