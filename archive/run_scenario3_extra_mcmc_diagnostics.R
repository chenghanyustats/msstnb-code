# ============================================================
# Scenario 3 extra MCMC diagnostics
# ============================================================
# Produces simple base R diagnostics that do not require extra packages.
# The key purpose is to check whether learning gamma harms the chain or the
# main inferential targets relative to Scenario 2.
# ============================================================

.s3d_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) cat(..., "\n")
}

.s3d_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s3d_get_nested <- function(x, names_vec) {
    for (nm in names_vec) {
        if (!is.null(x[[nm]])) return(x[[nm]])
    }
    NULL
}

.s3d_get_draws <- function(fit) {
    if (!is.null(fit$draws)) return(fit$draws)
    if (!is.null(fit$samples)) return(fit$samples)
    if (!is.null(fit$mcmc)) return(fit$mcmc)
    list()
}

.s3d_as_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.vector(x) && !is.list(x)) return(matrix(as.numeric(x), ncol = 1L))
    if (is.data.frame(x)) return(as.matrix(x))
    if (is.matrix(x)) return(x)
    if (is.array(x)) return(matrix(as.numeric(x), nrow = dim(x)[[1L]]))
    NULL
}

.s3d_parse_rep <- function(path) {
    base <- basename(path)
    m <- regmatches(base, regexpr("rep[0-9]+", base))
    if (length(m) == 0L || is.na(m)) return(NA_integer_)
    as.integer(sub("rep", "", m))
}

.s3d_acf1 <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 3L) return(NA_real_)
    suppressWarnings(stats::cor(x[-length(x)], x[-1L]))
}

.s3d_ess1 <- function(x, max_lag = NULL) {
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

.s3d_parameter_rows <- function(mat, prefix, rep_id, TT) {
    if (is.null(mat)) return(NULL)
    out <- vector("list", ncol(mat))
    for (j in seq_len(ncol(mat))) {
        nm <- if (!is.null(colnames(mat))) colnames(mat)[[j]] else paste0(prefix, j)
        x <- mat[, j]
        out[[j]] <- data.frame(
            scenario = "S3_DYNAMIC_LEARNED_GAMMA",
            rep_id = rep_id,
            TT = TT,
            parameter = nm,
            n_draws = sum(is.finite(x)),
            mean = mean(x, na.rm = TRUE),
            sd = stats::sd(x, na.rm = TRUE),
            acf1 = .s3d_acf1(x),
            ess = .s3d_ess1(x),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, out)
}

.s3d_extract_acceptance <- function(fit, rep_id, TT) {
    candidates <- list(fit$acceptance, fit$accept, fit$diagnostics, fit$mcmc_diagnostics)
    rows <- list()
    for (obj in candidates) {
        if (is.null(obj)) next
        if (is.data.frame(obj)) {
            obj$scenario <- "S3_DYNAMIC_LEARNED_GAMMA"
            obj$rep_id <- rep_id
            obj$TT <- TT
            return(obj)
        }
        if (is.list(obj)) {
            nums <- unlist(obj)
            nums <- nums[is.finite(nums)]
            if (length(nums) > 0L) {
                rows[[length(rows) + 1L]] <- data.frame(
                    scenario = "S3_DYNAMIC_LEARNED_GAMMA",
                    rep_id = rep_id,
                    TT = TT,
                    diagnostic = names(nums),
                    value = as.numeric(nums),
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    if (length(rows)) do.call(rbind, rows) else data.frame()
}

.s3d_plot_trace <- function(draw_vec, title, file) {
    draw_vec <- as.numeric(draw_vec)
    draw_vec <- draw_vec[is.finite(draw_vec)]
    if (length(draw_vec) < 2L) return(FALSE)
    grDevices::pdf(file, width = 8, height = 4.5)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(draw_vec, type = "l", xlab = "Saved draw", ylab = "Value", main = title)
    invisible(TRUE)
}

run_scenario3_extra_mcmc_diagnostics <- function(
    root = ".",
    fit_dir = "fits_s3_dynamic_learned_gamma",
    analysis_dir = "analysis_s3_dynamic_learned_gamma",
    fit_pattern = "s3_dynamic_learned_gamma.*\\.rds$",
    make_plots = TRUE,
    verbose = TRUE
) {
    root <- .s3d_norm_path(root)
    fit_dir <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)
    plot_dir <- file.path(analysis_dir, "mcmc_plots")
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    fit_files <- list.files(fit_dir, pattern = fit_pattern, full.names = TRUE)
    if (length(fit_files) == 0L) {
        stop("No Scenario 3 fit files found in: ", fit_dir, call. = FALSE)
    }

    diag_rows <- list()
    acc_rows <- list()

    for (ff in fit_files) {
        .s3d_msg("Reading: ", ff, verbose = verbose)
        fit <- readRDS(ff)
        rep_id <- if (!is.null(fit$rep_id)) as.integer(fit$rep_id) else .s3d_parse_rep(ff)
        TT <- if (!is.null(fit$TT)) as.integer(fit$TT) else NA_integer_
        draws <- .s3d_get_draws(fit)

        gamma_mat <- .s3d_as_matrix(.s3d_get_nested(draws, c("gamma", "gamma_draws", "gamma_save")))
        beta_mat <- .s3d_as_matrix(.s3d_get_nested(draws, c("beta", "beta_draws", "beta_save")))
        delta_mat <- .s3d_as_matrix(.s3d_get_nested(draws, c("delta", "delta_draws", "delta_save")))
        r_mat <- .s3d_as_matrix(.s3d_get_nested(draws, c("r", "r_draws", "dispersion", "dispersion_draws")))

        if (!is.null(gamma_mat)) {
            colnames(gamma_mat) <- "gamma"
            diag_rows[[length(diag_rows) + 1L]] <- .s3d_parameter_rows(gamma_mat, "gamma", rep_id, TT)
            if (isTRUE(make_plots)) {
                .s3d_plot_trace(
                    gamma_mat[, 1L],
                    sprintf("Scenario 3 gamma trace, rep %s", rep_id),
                    file.path(plot_dir, sprintf("trace_gamma_rep%02d.pdf", rep_id))
                )
            }
        }
        if (!is.null(beta_mat)) {
            diag_rows[[length(diag_rows) + 1L]] <- .s3d_parameter_rows(beta_mat, "beta", rep_id, TT)
            if (isTRUE(make_plots)) {
                for (j in seq_len(min(3L, ncol(beta_mat)))) {
                    nm <- if (!is.null(colnames(beta_mat))) colnames(beta_mat)[[j]] else paste0("beta", j - 1L)
                    .s3d_plot_trace(
                        beta_mat[, j],
                        sprintf("Scenario 3 %s trace, rep %s", nm, rep_id),
                        file.path(plot_dir, sprintf("trace_%s_rep%02d.pdf", nm, rep_id))
                    )
                }
            }
        }
        if (!is.null(delta_mat)) {
            diag_rows[[length(diag_rows) + 1L]] <- .s3d_parameter_rows(delta_mat, "delta", rep_id, TT)
        }
        if (!is.null(r_mat)) {
            diag_rows[[length(diag_rows) + 1L]] <- .s3d_parameter_rows(r_mat, "r", rep_id, TT)
        }

        acc <- .s3d_extract_acceptance(fit, rep_id, TT)
        if (nrow(acc) > 0L) acc_rows[[length(acc_rows) + 1L]] <- acc
    }

    mcmc_diag <- if (length(diag_rows)) do.call(rbind, diag_rows) else data.frame()
    acc_diag <- if (length(acc_rows)) do.call(rbind, acc_rows) else data.frame()

    diag_file <- file.path(analysis_dir, "mcmc_extra_diagnostics.csv")
    acc_file <- file.path(analysis_dir, "mcmc_acceptance_diagnostics.csv")
    write.csv(mcmc_diag, diag_file, row.names = FALSE)
    write.csv(acc_diag, acc_file, row.names = FALSE)

    .s3d_msg("Saved: ", diag_file, verbose = verbose)
    .s3d_msg("Saved: ", acc_file, verbose = verbose)
    if (isTRUE(make_plots)) .s3d_msg("Saved trace plots in: ", plot_dir, verbose = verbose)

    invisible(list(
        mcmc_diagnostics = mcmc_diag,
        acceptance_diagnostics = acc_diag,
        files = list(
            mcmc = diag_file,
            acceptance = acc_file,
            plot_dir = plot_dir
        )
    ))
}

if (sys.nframe() == 0L) {
    run_scenario3_extra_mcmc_diagnostics(root = ".")
}
