# ============================================================
# Minimal gamma audit rerun script -- v2
# Fix:
#   resets internal gamma audit counter before running
# ============================================================

run_minimal_gamma_audit <- function(
    project_dir,
    data_file = "data_revised/S1/data_rep01.rds",
    fit_file = "gamma_audit_fit_S1_rep01.rds",
    iter_audit = 2000L,
    region_audit = 1L,
    grid = seq(0.50, 0.99, length.out = 200L),
    n_iter = 2500L,
    n_burnin = 1000L,
    n_thin = 5L,
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

    settings_audit <- MCMC_SETTINGS
    settings_audit$n_iter   <- as.integer(n_iter)
    settings_audit$n_burnin <- as.integer(n_burnin)
    settings_audit$n_thin   <- as.integer(n_thin)

    priors <- MCMC_PRIORS

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

    gamma_audit_env <- new.env(parent = emptyenv())

    options(
        MSSTNB_GAMMA_AUDIT = list(
            enabled = TRUE,
            iter = as.integer(iter_audit),
            region = as.integer(region_audit),
            grid = as.numeric(grid),
            env = gamma_audit_env,
            label = "minimal_gamma_audit"
        ),
        MSSTNB_GAMMA_AUDIT_COUNTER = 0L
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
        cat("Running minimal gamma audit rerun\n")
        cat("============================================================\n")
        cat("Project dir   :", project_dir, "\n")
        cat("Data file     :", data_path, "\n")
        cat("Fit file      :", fit_path, "\n")
        cat("Audit iter    :", iter_audit, "\n")
        cat("Audit region  :", region_audit, "\n")
        cat("n_iter        :", n_iter, "\n")
        cat("n_burnin      :", n_burnin, "\n")
        cat("n_thin        :", n_thin, "\n\n")
    }

    fit <- run_mcmc(
        dat,
        settings_audit,
        priors,
        spatial,
        constants,
        verbose = 500L
    )

    saveRDS(fit, fit_path)

    recs <- gamma_audit_env$.gamma_audit_records
    nrec <- if (is.null(recs)) 0L else length(recs)

    if (nrec == 0L) {
        warning("No gamma audit record was captured. ",
                "Check that update_gamma.R is the audit version and that the ",
                "requested iter/region were actually visited.")
        return(invisible(list(
            fit_file = fit_path,
            audit_records = NULL,
            audit_env = gamma_audit_env
        )))
    }

    audit_rds <- sub("\\.rds$", "_audit_record.rds", fit_path)
    saveRDS(recs, audit_rds)

    rec1 <- recs[[1L]]
    audit_txt <- sub("\\.rds$", "_audit_summary.txt", fit_path)
    lines <- c(
        "============================================================",
        "Minimal gamma audit summary",
        "============================================================",
        paste("Fit file                :", fit_path),
        paste("Audit record count      :", nrec),
        paste("Recorded iter           :", rec1$iter),
        paste("Recorded region         :", rec1$region),
        paste("gamma_current          :", round(rec1$gamma_current, 6)),
        paste("gamma_proposal         :", round(rec1$gamma_proposal, 6)),
        paste("log_target_current     :", round(rec1$log_target_current, 6)),
        paste("log_target_proposal    :", round(rec1$log_target_proposal, 6)),
        paste("internal_curve_mode    :", round(rec1$internal_curve_mode, 6)),
        paste("internal_curve_mean    :", round(rec1$internal_curve_mean, 6)),
        paste("accepted               :", rec1$accepted),
        paste("mh_sd_used             :", round(rec1$mh_sd_used, 6))
    )
    writeLines(lines, con = audit_txt)

    if (verbose) {
        cat("\n")
        cat(paste(lines, collapse = "\n"))
        cat("\n\nSaved audit files:\n")
        cat("  ", audit_rds, "\n", sep = "")
        cat("  ", audit_txt, "\n", sep = "")
    }

    invisible(list(
        fit_file = fit_path,
        audit_records = recs,
        audit_rds = audit_rds,
        audit_txt = audit_txt,
        audit_env = gamma_audit_env
    ))
}
