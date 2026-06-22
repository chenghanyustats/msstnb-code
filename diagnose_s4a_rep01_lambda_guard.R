## ============================================================================
## diagnose_s4a_rep01_lambda_guard.R
##
## Diagnostic helper for S4A continuous-time x2 fixed-gamma fit.
## It checks whether lambda output guards are only numerical/burn-in artifacts
## or whether the stored posterior draws still sit on the beta0--lambda ridge.
## ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.s4a_diag_as_num <- function(x) as.numeric(x)

.s4a_diag_split_summary <- function(v, name) {
    v <- as.numeric(v)
    n <- length(v)
    i1 <- seq_len(floor(n / 2))
    i2 <- seq.int(floor(n / 2) + 1L, n)
    data.frame(
        quantity = name,
        part = c("all", "first_half", "second_half"),
        mean = c(mean(v, na.rm = TRUE), mean(v[i1], na.rm = TRUE), mean(v[i2], na.rm = TRUE)),
        sd = c(stats::sd(v, na.rm = TRUE), stats::sd(v[i1], na.rm = TRUE), stats::sd(v[i2], na.rm = TRUE)),
        q025 = c(stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE),
                 stats::quantile(v[i1], 0.025, na.rm = TRUE, names = FALSE),
                 stats::quantile(v[i2], 0.025, na.rm = TRUE, names = FALSE)),
        q50 = c(stats::quantile(v, 0.50, na.rm = TRUE, names = FALSE),
                stats::quantile(v[i1], 0.50, na.rm = TRUE, names = FALSE),
                stats::quantile(v[i2], 0.50, na.rm = TRUE, names = FALSE)),
        q975 = c(stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE),
                 stats::quantile(v[i1], 0.975, na.rm = TRUE, names = FALSE),
                 stats::quantile(v[i2], 0.975, na.rm = TRUE, names = FALSE)),
        stringsAsFactors = FALSE
    )
}

run_s4a_lambda_guard_diagnostic <- function(
    fit_file = "output_s4a_sparse_counts_continuous_x2/S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100/fits/fit_rep01.rds",
    data_file = "data_s4a_sparse_counts_continuous_x2/S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100/data_rep01.rds",
    out_dir = "output_s4a_sparse_counts_continuous_x2/S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100/diagnostics",
    out_prefix = "s4a_rep01_lambda_guard",
    lower = 1e-10,
    upper = 1e10,
    make_plots = TRUE
) {
    fit <- readRDS(fit_file)
    dat <- readRDS(data_file)

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    beta0 <- as.numeric(fit$samples$beta0)
    beta <- as.matrix(fit$samples$beta)
    phi <- as.matrix(fit$samples$phi)
    lambda <- fit$samples$lambda_tilde

    if (length(dim(lambda)) != 3L) {
        stop("fit$samples$lambda_tilde must be a draw x time x area array.", call. = FALSE)
    }
    n_draw <- dim(lambda)[1L]
    TT <- dim(lambda)[2L]
    n1 <- dim(lambda)[3L]

    lambda_mat <- matrix(as.numeric(lambda), nrow = n_draw)
    log_lambda_mat <- log(lambda_mat)

    lambda_min <- apply(lambda_mat, 1L, min, na.rm = TRUE)
    lambda_median <- apply(lambda_mat, 1L, stats::median, na.rm = TRUE)
    lambda_max <- apply(lambda_mat, 1L, max, na.rm = TRUE)
    log_lambda_mean <- rowMeans(log_lambda_mat, na.rm = TRUE)
    log_lambda_median <- apply(log_lambda_mat, 1L, stats::median, na.rm = TRUE)
    prop_at_lower <- rowMeans(lambda_mat <= lower * (1 + 1e-12), na.rm = TRUE)
    prop_at_upper <- rowMeans(lambda_mat >= upper * (1 - 1e-12), na.rm = TRUE)
    prop_nonfinite <- rowMeans(!is.finite(lambda_mat), na.rm = TRUE)

    ## Recenter each stored draw if the project recenter() function is available.
    beta0_ident <- rep(NA_real_, n_draw)
    lambda_ident_max <- rep(NA_real_, n_draw)
    lambda_ident_median <- rep(NA_real_, n_draw)
    lambda_ident_log_colmean_abs_mean <- rep(NA_real_, n_draw)

    if (exists("recenter", mode = "function", inherits = TRUE)) {
        recenter_fun <- get("recenter", mode = "function", inherits = TRUE)
        for (d in seq_len(n_draw)) {
            rc <- recenter_fun(
                beta0 = beta0[d],
                phi = phi[d, ],
                lambda_tilde = lambda[d, , ],
                return_diag = TRUE
            )
            beta0_ident[d] <- rc$beta0
            lambda_ident_max[d] <- max(rc$lambda_tilde, na.rm = TRUE)
            lambda_ident_median[d] <- stats::median(as.numeric(rc$lambda_tilde), na.rm = TRUE)
            lambda_ident_log_colmean_abs_mean[d] <- mean(abs(colMeans(log(rc$lambda_tilde))), na.rm = TRUE)
        }
    } else {
        warning("Project recenter() function was not found. Identified-scale summaries are omitted.")
    }

    draw_df <- data.frame(
        draw = seq_len(n_draw),
        beta0_raw = beta0,
        beta1 = beta[, 1L],
        beta2 = beta[, 2L],
        lambda_min = lambda_min,
        lambda_median = lambda_median,
        lambda_max = lambda_max,
        log_lambda_mean = log_lambda_mean,
        log_lambda_median = log_lambda_median,
        prop_at_lower = prop_at_lower,
        prop_at_upper = prop_at_upper,
        prop_nonfinite = prop_nonfinite,
        beta0_plus_mean_log_lambda = beta0 + log_lambda_mean,
        beta0_plus_median_log_lambda = beta0 + log_lambda_median,
        beta0_ident = beta0_ident,
        lambda_ident_median = lambda_ident_median,
        lambda_ident_max = lambda_ident_max,
        lambda_ident_log_colmean_abs_mean = lambda_ident_log_colmean_abs_mean
    )

    guard_counts <- fit$diagnostics$s4a_numeric_guards %||% list()
    n_iter <- fit$settings$n_iter %||% fit$settings$iter %||% NA_integer_
    total_cells_all_iter <- if (is.finite(n_iter)) n_iter * TT * n1 else NA_real_
    lambda_output_guard_rate_all_iter <- if (is.finite(total_cells_all_iter)) {
        (guard_counts$n_lambda_output_guard %||% NA_real_) / total_cells_all_iter
    } else NA_real_

    scalar_df <- data.frame(
        fit_file = fit_file,
        data_file = data_file,
        n_iter = n_iter,
        n_stored = n_draw,
        TT = TT,
        n1 = n1,
        lambda_output_guard_count = guard_counts$n_lambda_output_guard %||% NA_real_,
        lambda_output_guard_rate_all_iter = lambda_output_guard_rate_all_iter,
        stored_prop_at_lower_mean = mean(prop_at_lower, na.rm = TRUE),
        stored_prop_at_upper_mean = mean(prop_at_upper, na.rm = TRUE),
        stored_prop_at_lower_max = max(prop_at_lower, na.rm = TRUE),
        stored_prop_at_upper_max = max(prop_at_upper, na.rm = TRUE),
        beta0_raw_mean = mean(beta0, na.rm = TRUE),
        beta0_ident_mean = mean(beta0_ident, na.rm = TRUE),
        beta0_sparse_truth = dat$beta0_sparse_truth %||% NA_real_,
        beta0_star_ident_truth = dat$beta0_star_ident %||% NA_real_,
        beta1_mean = mean(beta[, 1L], na.rm = TRUE),
        beta1_truth = if (!is.null(dat$beta_star)) dat$beta_star[1L] else NA_real_,
        beta2_mean = mean(beta[, 2L], na.rm = TRUE),
        beta2_truth = if (!is.null(dat$beta_star)) dat$beta_star[2L] else NA_real_,
        r_mean = mean(as.numeric(fit$samples$r), na.rm = TRUE),
        r_truth_mean = mean(dat$r_star %||% NA_real_, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    split_df <- do.call(rbind, list(
        .s4a_diag_split_summary(draw_df$beta0_raw, "beta0_raw"),
        .s4a_diag_split_summary(draw_df$beta0_ident, "beta0_ident"),
        .s4a_diag_split_summary(draw_df$beta0_plus_mean_log_lambda, "beta0_plus_mean_log_lambda"),
        .s4a_diag_split_summary(log(draw_df$lambda_max), "log_lambda_max"),
        .s4a_diag_split_summary(draw_df$prop_at_upper, "prop_at_upper"),
        .s4a_diag_split_summary(draw_df$beta1, "beta1"),
        .s4a_diag_split_summary(draw_df$beta2, "beta2")
    ))

    utils::write.csv(draw_df, file.path(out_dir, paste0(out_prefix, "_draw_trace.csv")), row.names = FALSE)
    utils::write.csv(scalar_df, file.path(out_dir, paste0(out_prefix, "_scalar_summary.csv")), row.names = FALSE)
    utils::write.csv(split_df, file.path(out_dir, paste0(out_prefix, "_split_summary.csv")), row.names = FALSE)

    if (isTRUE(make_plots)) {
        pdf_file <- file.path(out_dir, paste0(out_prefix, "_trace_plots.pdf"))
        grDevices::pdf(pdf_file, width = 8, height = 5)
        on.exit(grDevices::dev.off(), add = TRUE)

        graphics::plot(draw_df$draw, draw_df$beta0_raw, type = "l",
                       xlab = "Stored draw", ylab = "Raw beta0",
                       main = "S4A rep 01: raw beta0 trace")
        if (is.finite(dat$beta0_sparse_truth %||% NA_real_)) graphics::abline(h = dat$beta0_sparse_truth, lty = 2)

        graphics::plot(draw_df$draw, draw_df$beta0_ident, type = "l",
                       xlab = "Stored draw", ylab = "Identified beta0",
                       main = "S4A rep 01: identified beta0 trace")
        if (is.finite(dat$beta0_star_ident %||% NA_real_)) graphics::abline(h = dat$beta0_star_ident, lty = 2)

        graphics::plot(draw_df$draw, draw_df$beta0_plus_mean_log_lambda, type = "l",
                       xlab = "Stored draw", ylab = "beta0 + mean(log lambda)",
                       main = "S4A rep 01: beta0-lambda combined scale")
        if (is.finite(dat$beta0_star_ident %||% NA_real_)) graphics::abline(h = dat$beta0_star_ident, lty = 2)

        graphics::plot(draw_df$draw, log(draw_df$lambda_max), type = "l",
                       xlab = "Stored draw", ylab = "log(max lambda_tilde)",
                       main = "S4A rep 01: lambda upper-tail trace")

        graphics::plot(draw_df$draw, draw_df$prop_at_upper, type = "l",
                       xlab = "Stored draw", ylab = "Proportion at upper clamp",
                       main = "S4A rep 01: stored lambda upper-clamp mass")

        graphics::plot(draw_df$draw, draw_df$beta1, type = "l",
                       xlab = "Stored draw", ylab = "beta1",
                       main = "S4A rep 01: beta1 trace")
        if (!is.null(dat$beta_star)) graphics::abline(h = dat$beta_star[1L], lty = 2)

        graphics::plot(draw_df$draw, draw_df$beta2, type = "l",
                       xlab = "Stored draw", ylab = "beta2",
                       main = "S4A rep 01: beta2 trace")
        if (!is.null(dat$beta_star)) graphics::abline(h = dat$beta_star[2L], lty = 2)
    }

    message("Saved diagnostic files to: ", out_dir)
    list(
        scalar_summary = scalar_df,
        split_summary = split_df,
        draw_trace = draw_df,
        out_dir = out_dir
    )
}
