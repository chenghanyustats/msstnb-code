# ============================================================
# Long-chain rerun script for S1 rep01:03, 2 chains each,
# with gamma prior changed to Beta(1,1)
#
# Purpose:
#   1. Re-run S1 replications 01:03
#   2. Run 2 chains per replication
#   3. Use the current ordinary sampler
#   4. Override ONLY gamma prior to Beta(1,1)
#   5. Save fit objects and a compact summary
#
# Typical use:
#
# source("s1_longchains_beta11.R")
#
# res <- run_s1_longchains_beta11(
#     project_dir = "/Users/chenghanyu/Dropbox/academia/research/msstnb_new/simulation",
#     reps = 1:3,
#     chain_ids = 1:2,
#     data_dir = "data_revised",
#     fit_dir = "longchain_fits_S1_beta11",
#     out_prefix = "s1_longchain_beta11",
#     n_iter = 10000L,
#     n_burnin = 4000L,
#     n_thin = 5L,
#     overwrite_existing = TRUE,
#     verbose = TRUE
# )
# ============================================================

run_s1_longchains_beta11 <- function(
    project_dir,
    reps = 1:3,
    chain_ids = 1:2,
    data_dir = "data_revised",
    fit_dir = "longchain_fits_S1_beta11",
    out_prefix = "s1_longchain_beta11",
    n_iter = 10000L,
    n_burnin = 4000L,
    n_thin = 5L,
    overwrite_existing = TRUE,
    verbose = TRUE
) {
    stopifnot(length(project_dir) == 1L, dir.exists(project_dir))

    find_file <- function(filename, roots) {
        for (root in roots) {
            cand <- file.path(root, filename)
            if (file.exists(cand)) {
                return(normalizePath(cand, winslash = "/", mustWork = TRUE))
            }
        }
        stop("Could not find required file: ", filename,
             "\nSearched in:\n", paste(roots, collapse = "\n"))
    }

    safe_mean <- function(x) {
        if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_real_)
        mean(x, na.rm = TRUE)
    }

    safe_sd <- function(x) {
        if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_real_)
        stats::sd(x, na.rm = TRUE)
    }

    roots <- c(
        project_dir,
        file.path(project_dir, "R"),
        file.path(project_dir, "R", "mcmc"),
        file.path(project_dir, "mcmc"),
        file.path(project_dir, "scripts"),
        file.path(project_dir, "code"),
        file.path(project_dir, "R_gpt"),
        file.path(project_dir, "R_gpt", "mcmc")
    )

    f_setup   <- find_file("00_setup.R", roots)
    f_cfg     <- find_file("mcmc_config.R", roots)
    f_utils   <- find_file("mcmc_utils.R", roots)
    f_kappa   <- find_file("update_kappa.R", roots)
    f_reg     <- find_file("update_regression.R", roots)
    f_icar    <- find_file("update_icar.R", roots)
    f_disp    <- find_file("update_dispersion.R", roots)
    f_gamma   <- find_file("update_gamma.R", roots)
    f_ffbs    <- find_file("ffbs_lambda.R", roots)
    f_delta   <- find_file("update_delta.R", roots)
    f_omega   <- find_file("smooth_omega.R", roots)
    f_sampler <- find_file("sampler.R", roots)

    source(f_setup, local = FALSE)
    source(f_cfg, local = FALSE)
    source(f_utils, local = FALSE)
    source(f_kappa, local = FALSE)
    source(f_reg, local = FALSE)
    source(f_icar, local = FALSE)
    source(f_disp, local = FALSE)
    source(f_gamma, local = FALSE)
    source(f_ffbs, local = FALSE)
    source(f_delta, local = FALSE)
    source(f_omega, local = FALSE)
    source(f_sampler, local = FALSE)

    settings_long <- MCMC_SETTINGS
    settings_long$n_iter   <- as.integer(n_iter)
    settings_long$n_burnin <- as.integer(n_burnin)
    settings_long$n_thin   <- as.integer(n_thin)

    priors <- MCMC_PRIORS
    priors$gamma_a <- 1
    priors$gamma_b <- 1

    spatial <- list(
        H      = H,
        B_ICAR = B_ICAR,
        BHB    = crossprod(B_ICAR, H %*% B_ICAR)
    )

    constants <- list(
        A0 = A0,
        B0 = B0,
        C0 = C0
    )

    fit_dir_abs <- file.path(project_dir, fit_dir)
    dir.create(fit_dir_abs, recursive = TRUE, showWarnings = FALSE)

    rows <- list()
    idx <- 0L

    for (rep_id in reps) {
        rr <- sprintf("%02d", rep_id)
        dat_file <- file.path(project_dir, data_dir, "S1", paste0("data_rep", rr, ".rds"))

        if (verbose) {
            cat("\n============================================================\n")
            cat(sprintf("S1 rep %s\n", rr))
            cat("============================================================\n")
        }

        if (!file.exists(dat_file)) {
            if (verbose) cat("Missing data file:", dat_file, "\n")
            for (chain_id in chain_ids) {
                idx <- idx + 1L
                rows[[idx]] <- data.frame(
                    rep_id = rep_id,
                    chain_id = chain_id,
                    fit_file = NA_character_,
                    status = "missing_data",
                    gamma1_mean = NA_real_,
                    gamma1_sd = NA_real_,
                    beta1_mean = NA_real_,
                    beta2_mean = NA_real_,
                    phi_cor = NA_real_,
                    lambda_log_rmse = NA_real_,
                    lambda_log_cor = NA_real_,
                    lambda_cover = NA_real_,
                    gamma_accept = NA_real_,
                    beta_mean_n_reject = NA_real_,
                    stringsAsFactors = FALSE
                )
            }
            next
        }

        dat <- readRDS(dat_file)

        for (chain_id in chain_ids) {
            fit_file <- file.path(
                fit_dir_abs,
                sprintf("fit_S1_rep%s_chain%02d_long_beta11.rds", rr, chain_id)
            )

            if (file.exists(fit_file) && !overwrite_existing) {
                if (verbose) cat(sprintf("Chain %02d: using existing fit %s\n", chain_id, basename(fit_file)))
                fit <- readRDS(fit_file)
                status <- "ok_existing"
            } else {
                if ("seed" %in% names(settings_long)) {
                    settings_long$seed <- as.integer(200000L + rep_id * 100L + chain_id)
                }

                if (verbose) {
                    cat(sprintf("Chain %02d: running long chain with gamma prior Beta(1,1)...\n", chain_id))
                }

                fit <- tryCatch(
                    run_mcmc(
                        dat,
                        settings_long,
                        priors,
                        spatial,
                        constants,
                        verbose = 1000L
                    ),
                    error = function(e) e
                )

                if (inherits(fit, "error")) {
                    if (verbose) {
                        cat(sprintf("Chain %02d FAILED: %s\n", chain_id, conditionMessage(fit)))
                    }
                    idx <- idx + 1L
                    rows[[idx]] <- data.frame(
                        rep_id = rep_id,
                        chain_id = chain_id,
                        fit_file = fit_file,
                        status = paste0("fit_error: ", conditionMessage(fit)),
                        gamma1_mean = NA_real_,
                        gamma1_sd = NA_real_,
                        beta1_mean = NA_real_,
                        beta2_mean = NA_real_,
                        phi_cor = NA_real_,
                        lambda_log_rmse = NA_real_,
                        lambda_log_cor = NA_real_,
                        lambda_cover = NA_real_,
                        gamma_accept = NA_real_,
                        beta_mean_n_reject = NA_real_,
                        stringsAsFactors = FALSE
                    )
                    next
                }

                saveRDS(fit, fit_file)
                status <- "ok"
                if (verbose) cat(sprintf("Chain %02d: saved %s\n", chain_id, basename(fit_file)))
            }

            s <- fit$samples

            gamma1 <- if (!is.null(s$gamma)) s$gamma[, 1] else numeric(0)
            beta1  <- if (!is.null(s$beta))  s$beta[, 1] else numeric(0)
            beta2  <- if (!is.null(s$beta))  s$beta[, 2] else numeric(0)

            phi_cor <- NA_real_
            if (!is.null(s$phi) && !is.null(dat$phi_star_ident)) {
                phi_pm <- colMeans(s$phi)
                phi_cor <- suppressWarnings(stats::cor(phi_pm, dat$phi_star_ident))
            }

            lambda_log_rmse <- NA_real_
            lambda_log_cor  <- NA_real_
            lambda_cover    <- NA_real_
            if (!is.null(s$lambda_tilde) && !is.null(dat$lambda_tilde_ident)) {
                lt_draws <- s$lambda_tilde
                truth_lt <- dat$lambda_tilde_ident

                if (length(dim(lt_draws)) == 3L) {
                    post_mean_log <- apply(log(lt_draws), c(2, 3), mean, na.rm = TRUE)
                    post_q025 <- apply(lt_draws, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
                    post_q975 <- apply(lt_draws, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)
                    truth_log <- log(truth_lt)

                    lambda_log_cor <- suppressWarnings(stats::cor(c(post_mean_log), c(truth_log)))
                    lambda_log_rmse <- sqrt(mean((c(post_mean_log) - c(truth_log))^2, na.rm = TRUE))
                    lambda_cover <- mean(truth_lt >= post_q025 & truth_lt <= post_q975, na.rm = TRUE)
                }
            }

            gamma_accept <- if (!is.null(fit$diagnostics$gamma_accept_rate)) {
                safe_mean(fit$diagnostics$gamma_accept_rate)
            } else {
                NA_real_
            }

            beta_mean_n_reject <- if (!is.null(fit$diagnostics$beta_mean_n_reject)) {
                as.numeric(fit$diagnostics$beta_mean_n_reject)
            } else {
                NA_real_
            }

            idx <- idx + 1L
            rows[[idx]] <- data.frame(
                rep_id = rep_id,
                chain_id = chain_id,
                fit_file = fit_file,
                status = status,
                gamma1_mean = safe_mean(gamma1),
                gamma1_sd = safe_sd(gamma1),
                beta1_mean = safe_mean(beta1),
                beta2_mean = safe_mean(beta2),
                phi_cor = phi_cor,
                lambda_log_rmse = lambda_log_rmse,
                lambda_log_cor = lambda_log_cor,
                lambda_cover = lambda_cover,
                gamma_accept = gamma_accept,
                beta_mean_n_reject = beta_mean_n_reject,
                stringsAsFactors = FALSE
            )
        }
    }

    out_df <- do.call(rbind, rows)

    out_csv <- file.path(fit_dir_abs, paste0(out_prefix, "_summary.csv"))
    out_txt <- file.path(fit_dir_abs, paste0(out_prefix, "_summary.txt"))
    utils::write.csv(out_df, out_csv, row.names = FALSE)

    ok <- out_df$status %in% c("ok", "ok_existing")

    lines <- c(
        "============================================================",
        "Long-chain rerun summary: S1 ordinary fit with gamma prior Beta(1,1)",
        "============================================================",
        paste("Project dir                  :", project_dir),
        paste("Replications requested       :", length(reps)),
        paste("Chains per replication       :", length(chain_ids)),
        paste("Rows in summary              :", nrow(out_df)),
        paste("Successful fits              :", sum(ok)),
        paste("Mean gamma1 posterior mean   :", round(safe_mean(out_df$gamma1_mean[ok]), 4)),
        paste("Mean gamma1 posterior SD     :", round(safe_mean(out_df$gamma1_sd[ok]), 4)),
        paste("Mean beta1 posterior mean    :", round(safe_mean(out_df$beta1_mean[ok]), 4)),
        paste("Mean beta2 posterior mean    :", round(safe_mean(out_df$beta2_mean[ok]), 4)),
        paste("Mean lambda log-RMSE         :", round(safe_mean(out_df$lambda_log_rmse[ok]), 4)),
        "",
        "Status counts:",
        capture.output(print(table(out_df$status, useNA = "ifany")))
    )
    writeLines(lines, con = out_txt)

    if (verbose) {
        cat("\n")
        cat(paste(lines, collapse = "\n"))
        cat("\n\nWrote:\n", out_csv, "\n", out_txt, "\n", sep = "")
    }

    invisible(list(
        summary = out_df,
        csv = out_csv,
        txt = out_txt
    ))
}
