# ============================================================
# Scenario 3 posterior performance summaries
# ============================================================
# This script summarizes fitted Scenario 3 MCMC objects. It is intentionally
# robust to small naming differences in fit$draws and fit$truth.
# ============================================================

.s3p_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) cat(..., "\n")
}

.s3p_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s3p_q <- function(x, probs) {
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

.s3p_summarize_vector <- function(x, truth = NA_real_) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(data.frame(
            mean = NA_real_, sd = NA_real_, median = NA_real_,
            q025 = NA_real_, q25 = NA_real_, q75 = NA_real_, q975 = NA_real_,
            truth = truth, bias = NA_real_, covered_95 = NA,
            stringsAsFactors = FALSE
        ))
    }
    qs <- .s3p_q(x, c(0.025, 0.25, 0.50, 0.75, 0.975))
    mn <- mean(x)
    data.frame(
        mean = mn,
        sd = stats::sd(x),
        median = qs[[3L]],
        q025 = qs[[1L]],
        q25 = qs[[2L]],
        q75 = qs[[4L]],
        q975 = qs[[5L]],
        truth = truth,
        bias = ifelse(is.finite(truth), mn - truth, NA_real_),
        covered_95 = ifelse(is.finite(truth), truth >= qs[[1L]] && truth <= qs[[5L]], NA),
        stringsAsFactors = FALSE
    )
}

.s3p_get_nested <- function(x, names_vec) {
    for (nm in names_vec) {
        if (!is.null(x[[nm]])) return(x[[nm]])
    }
    NULL
}

.s3p_get_draws <- function(fit) {
    if (!is.null(fit$draws)) return(fit$draws)
    if (!is.null(fit$samples)) return(fit$samples)
    if (!is.null(fit$mcmc)) return(fit$mcmc)
    list()
}

.s3p_as_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.vector(x) && !is.list(x)) return(matrix(as.numeric(x), ncol = 1L))
    if (is.data.frame(x)) return(as.matrix(x))
    if (is.matrix(x)) return(x)
    if (is.array(x)) return(matrix(as.numeric(x), nrow = dim(x)[[1L]]))
    NULL
}

.s3p_truth_beta <- function(truth, p) {
    if (is.null(truth)) return(rep(NA_real_, p))

    beta_candidates <- c("beta", "beta_star", "beta_ident", "beta_star_ident")
    for (nm in beta_candidates) {
        if (!is.null(truth[[nm]])) {
            b <- as.numeric(truth[[nm]])
            if (length(b) >= p) return(b[seq_len(p)])
        }
    }

    out <- rep(NA_real_, p)
    names_to_try <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5")
    for (j in seq_len(min(p, length(names_to_try)))) {
        if (!is.null(truth[[names_to_try[[j]]]])) out[[j]] <- as.numeric(truth[[names_to_try[[j]]]])
    }

    if (p >= 1L && !is.null(truth$beta0_star_ident)) out[[1L]] <- as.numeric(truth$beta0_star_ident)
    if (p >= 2L && !is.null(truth$beta1_truth)) out[[2L]] <- as.numeric(truth$beta1_truth)
    if (p >= 3L && !is.null(truth$beta2_truth)) out[[3L]] <- as.numeric(truth$beta2_truth)
    out
}

.s3p_truth_gamma <- function(truth) {
    if (is.null(truth)) return(NA_real_)
    candidates <- c("gamma", "gamma_truth", "discount_gamma", "gamma_true")
    for (nm in candidates) {
        if (!is.null(truth[[nm]])) return(as.numeric(truth[[nm]])[[1L]])
    }
    NA_real_
}

.s3p_truth_matrix <- function(truth, names_vec) {
    if (is.null(truth)) return(NULL)
    x <- .s3p_get_nested(truth, names_vec)
    if (is.null(x)) return(NULL)
    as.matrix(x)
}

.s3p_get_fit_files <- function(fit_dir, pattern = "s3_dynamic_learned_gamma.*\\.rds$") {
    list.files(fit_dir, pattern = pattern, full.names = TRUE)
}

.s3p_parse_rep <- function(path) {
    base <- basename(path)
    m <- regmatches(base, regexpr("rep[0-9]+", base))
    if (length(m) == 0L || is.na(m)) return(NA_integer_)
    as.integer(sub("rep", "", m))
}

.s3p_matrix_rmse_summary <- function(draw_obj, truth_mat, parameter_name) {
    draw_mat <- .s3p_as_matrix(draw_obj)
    if (is.null(draw_mat) || is.null(truth_mat)) return(NULL)

    truth_vec <- as.numeric(truth_mat)
    if (ncol(draw_mat) != length(truth_vec)) {
        return(NULL)
    }

    post_mean <- colMeans(draw_mat, na.rm = TRUE)
    err <- post_mean - truth_vec
    data.frame(
        parameter = parameter_name,
        rmse = sqrt(mean(err^2, na.rm = TRUE)),
        mae = mean(abs(err), na.rm = TRUE),
        bias = mean(err, na.rm = TRUE),
        n_elements = length(truth_vec),
        stringsAsFactors = FALSE
    )
}

summarize_scenario3_posterior_performance <- function(
    root = ".",
    fit_dir = "fits_s3_dynamic_learned_gamma",
    analysis_dir = "analysis_s3_dynamic_learned_gamma",
    fit_pattern = "s3_dynamic_learned_gamma.*\\.rds$",
    verbose = TRUE
) {
    root <- .s3p_norm_path(root)
    fit_dir <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

    fit_files <- .s3p_get_fit_files(fit_dir, pattern = fit_pattern)
    if (length(fit_files) == 0L) {
        stop("No Scenario 3 fit files found in: ", fit_dir, call. = FALSE)
    }

    beta_rows <- list()
    gamma_rows <- list()
    latent_rows <- list()

    for (ff in fit_files) {
        .s3p_msg("Reading: ", ff, verbose = verbose)
        fit <- readRDS(ff)
        rep_id <- if (!is.null(fit$rep_id)) as.integer(fit$rep_id) else .s3p_parse_rep(ff)
        TT <- if (!is.null(fit$TT)) as.integer(fit$TT) else NA_integer_
        truth <- fit$truth
        draws <- .s3p_get_draws(fit)

        beta_draws <- .s3p_as_matrix(.s3p_get_nested(draws, c("beta", "beta_draws", "beta_save")))
        if (!is.null(beta_draws)) {
            beta_truth <- .s3p_truth_beta(truth, ncol(beta_draws))
            for (j in seq_len(ncol(beta_draws))) {
                tmp <- .s3p_summarize_vector(beta_draws[, j], truth = beta_truth[[j]])
                tmp$scenario <- "S3_DYNAMIC_LEARNED_GAMMA"
                tmp$rep_id <- rep_id
                tmp$TT <- TT
                tmp$parameter <- if (!is.null(colnames(beta_draws))) colnames(beta_draws)[[j]] else paste0("beta", j - 1L)
                beta_rows[[length(beta_rows) + 1L]] <- tmp
            }
        }

        gamma_draws <- .s3p_get_nested(draws, c("gamma", "gamma_draws", "gamma_save"))
        if (!is.null(gamma_draws)) {
            gamma_truth <- .s3p_truth_gamma(truth)
            tmp <- .s3p_summarize_vector(gamma_draws, truth = gamma_truth)
            tmp$scenario <- "S3_DYNAMIC_LEARNED_GAMMA"
            tmp$rep_id <- rep_id
            tmp$TT <- TT
            tmp$parameter <- "gamma"
            gamma_rows[[length(gamma_rows) + 1L]] <- tmp
        }

        lambda_truth <- .s3p_truth_matrix(truth, c("lambda_tilde_ident", "lambda_tilde", "lambda", "lambda_truth"))
        lambda_draws <- .s3p_get_nested(draws, c("lambda_tilde", "lambda_tilde_draws", "lambda", "lambda_draws"))
        tmp <- .s3p_matrix_rmse_summary(lambda_draws, lambda_truth, "lambda_tilde")
        if (!is.null(tmp)) {
            tmp$scenario <- "S3_DYNAMIC_LEARNED_GAMMA"
            tmp$rep_id <- rep_id
            tmp$TT <- TT
            latent_rows[[length(latent_rows) + 1L]] <- tmp
        }

        phi_truth <- .s3p_truth_matrix(truth, c("phi_star_ident", "phi_ident", "phi", "phi_truth"))
        phi_draws <- .s3p_get_nested(draws, c("phi", "phi_draws", "phi_save"))
        tmp <- .s3p_matrix_rmse_summary(phi_draws, phi_truth, "phi")
        if (!is.null(tmp)) {
            tmp$scenario <- "S3_DYNAMIC_LEARNED_GAMMA"
            tmp$rep_id <- rep_id
            tmp$TT <- TT
            latent_rows[[length(latent_rows) + 1L]] <- tmp
        }
    }

    beta_summary <- if (length(beta_rows)) do.call(rbind, beta_rows) else data.frame()
    gamma_summary <- if (length(gamma_rows)) do.call(rbind, gamma_rows) else data.frame()
    latent_summary <- if (length(latent_rows)) do.call(rbind, latent_rows) else data.frame()

    beta_file <- file.path(analysis_dir, "posterior_beta_summary.csv")
    gamma_file <- file.path(analysis_dir, "posterior_gamma_summary.csv")
    latent_file <- file.path(analysis_dir, "posterior_latent_risk_performance.csv")

    write.csv(beta_summary, beta_file, row.names = FALSE)
    write.csv(gamma_summary, gamma_file, row.names = FALSE)
    write.csv(latent_summary, latent_file, row.names = FALSE)

    # Replicate level performance table. This is intentionally conservative:
    # beta and gamma summaries are direct posterior summaries; latent metrics
    # are added only when the corresponding truth arrays are available.
    perf <- data.frame()
    if (nrow(beta_summary) > 0L) {
        perf <- aggregate(
            cbind(abs_bias = abs(beta_summary$bias), covered_95 = as.numeric(beta_summary$covered_95)) ~ rep_id + TT,
            data = beta_summary,
            FUN = function(x) mean(x, na.rm = TRUE)
        )
        names(perf)[names(perf) == "abs_bias"] <- "mean_abs_beta_bias"
        names(perf)[names(perf) == "covered_95"] <- "mean_beta_coverage_95"
    }
    perf_file <- file.path(analysis_dir, "posterior_performance_diagnostics.csv")
    write.csv(perf, perf_file, row.names = FALSE)

    .s3p_msg("Saved: ", beta_file, verbose = verbose)
    .s3p_msg("Saved: ", gamma_file, verbose = verbose)
    .s3p_msg("Saved: ", latent_file, verbose = verbose)
    .s3p_msg("Saved: ", perf_file, verbose = verbose)

    invisible(list(
        beta_summary = beta_summary,
        gamma_summary = gamma_summary,
        latent_summary = latent_summary,
        performance = perf,
        files = list(
            beta = beta_file,
            gamma = gamma_file,
            latent = latent_file,
            performance = perf_file
        )
    ))
}

# Allow command line execution by Rscript.
if (sys.nframe() == 0L) {
    summarize_scenario3_posterior_performance(root = ".")
}
