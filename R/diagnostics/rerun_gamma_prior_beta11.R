# ============================================================
# Minimal rerun script with gamma prior changed to Beta(1,1)
#
# Purpose:
#   Re-run one ordinary fit after replacing only the gamma prior by
#   Beta(1,1), while keeping the current sampler and other priors
#   unchanged.
#
# Typical use:
#
# source("rerun_gamma_prior_beta11.R")
#
# out <- rerun_gamma_prior_beta11(
#     project_dir = "/Users/chenghanyu/Dropbox/academia/research/msstnb_new/simulation",
#     data_file   = "data_revised/S1/data_rep01.rds",
#     fit_file    = "fit_S1_rep01_gamma_prior_beta11_short.rds",
#     n_iter      = 2000L,
#     n_burnin    = 1000L,
#     n_thin      = 2L,
#     verbose     = TRUE
# )
# ============================================================

rerun_gamma_prior_beta11 <- function(
    project_dir,
    data_file = "data_revised/S1/data_rep01.rds",
    fit_file = "fit_S1_rep01_gamma_prior_beta11_short.rds",
    n_iter = 2000L,
    n_burnin = 1000L,
    n_thin = 2L,
    overwrite_fit = TRUE,
    verbose = TRUE
) {
    stopifnot(length(project_dir) == 1L, dir.exists(project_dir))

    find_file <- function(filename, roots) {
        for (root in roots) {
            cand <- file.path(root, filename)
            if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
        }
        stop("Could not find required file: ", filename,
             "\nSearched in:\n", paste(roots, collapse = "\n"))
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

    data_path <- if (file.exists(data_file)) {
        normalizePath(data_file, winslash = "/", mustWork = TRUE)
    } else {
        normalizePath(file.path(project_dir, data_file), winslash = "/", mustWork = TRUE)
    }
    dat <- readRDS(data_path)

    settings_run <- MCMC_SETTINGS
    settings_run$n_iter   <- as.integer(n_iter)
    settings_run$n_burnin <- as.integer(n_burnin)
    settings_run$n_thin   <- as.integer(n_thin)

    priors <- MCMC_PRIORS

    # Only change gamma prior
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

    fit_path <- if (dirname(fit_file) == ".") {
        file.path(project_dir, fit_file)
    } else {
        fit_file
    }

    dir.create(dirname(fit_path), recursive = TRUE, showWarnings = FALSE)

    if (file.exists(fit_path) && overwrite_fit) {
        file.remove(fit_path)
    }

    if (verbose) {
        cat("============================================================\n")
        cat("Running ordinary fit with gamma prior Beta(1,1)\n")
        cat("============================================================\n")
        cat("Project dir   :", project_dir, "\n")
        cat("Data file     :", data_path, "\n")
        cat("Fit file      :", fit_path, "\n")
        cat("gamma prior   : Beta(1,1)\n")
        cat("n_iter        :", n_iter, "\n")
        cat("n_burnin      :", n_burnin, "\n")
        cat("n_thin        :", n_thin, "\n\n")
    }

    fit <- run_mcmc(
        dat,
        settings_run,
        priors,
        spatial,
        constants,
        verbose = 500L
    )

    saveRDS(fit, fit_path)

    s <- fit$samples
    gamma_mean <- if (!is.null(s$gamma)) mean(s$gamma[, 1], na.rm = TRUE) else NA_real_
    gamma_sd   <- if (!is.null(s$gamma)) stats::sd(s$gamma[, 1], na.rm = TRUE) else NA_real_
    beta1_mean <- if (!is.null(s$beta)) mean(s$beta[, 1], na.rm = TRUE) else NA_real_
    beta2_mean <- if (!is.null(s$beta)) mean(s$beta[, 2], na.rm = TRUE) else NA_real_

    summary_txt <- sub("\\.rds$", "_summary.txt", fit_path)
    lines <- c(
        "============================================================",
        "Ordinary rerun with gamma prior Beta(1,1)",
        "============================================================",
        paste("Fit file                :", fit_path),
        paste("gamma mean              :", round(gamma_mean, 6)),
        paste("gamma sd                :", round(gamma_sd, 6)),
        paste("beta1 mean              :", round(beta1_mean, 6)),
        paste("beta2 mean              :", round(beta2_mean, 6))
    )
    writeLines(lines, con = summary_txt)

    if (verbose) {
        cat(paste(lines, collapse = "\n"), "\n")
        cat("\nSaved:\n")
        cat("  ", fit_path, "\n", sep = "")
        cat("  ", summary_txt, "\n", sep = "")
    }

    invisible(list(
        fit_file = fit_path,
        summary_txt = summary_txt,
        fit = fit
    ))
}
