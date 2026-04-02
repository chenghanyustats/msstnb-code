## ==========================================================================
## run_single_fit_common_gamma.R
## Source common-gamma dependencies, load data, run MCMC, save results
## ==========================================================================

source_mcmc_common_gamma <- function() {
    source("R/00_setup.R")
    source("R/01_helpers.R")
    source("R/mcmc/mcmc_config.R")
    source("R/mcmc/mcmc_utils.R")
    source("R/mcmc/update_kappa.R")
    source("R/mcmc/update_regression.R")
    source("R/mcmc/update_icar.R")
    source("R/mcmc/update_dispersion.R")
    source("R/mcmc/update_gamma.R")
    source("R/mcmc/update_gamma_common.R")
    source("R/mcmc/ffbs_lambda.R")
    source("R/mcmc/update_delta.R")
    source("R/mcmc/smooth_omega.R")
    source("R/mcmc/sampler_common_gamma.R")
}

method_settings <- function(method) {
    s <- MCMC_SETTINGS
    if (method == "M1") {
        s$include_nb <- TRUE;  s$include_icar <- TRUE;  s$include_covariates <- TRUE
    } else if (method == "M2") {
        s$include_nb <- FALSE; s$include_icar <- TRUE;  s$include_covariates <- TRUE
    } else if (method == "M3") {
        s$include_nb <- TRUE;  s$include_icar <- FALSE; s$include_covariates <- TRUE
    } else if (method == "M4") {
        s$include_nb <- FALSE; s$include_icar <- FALSE; s$include_covariates <- FALSE
    } else {
        stop("Unknown method: ", method)
    }
    s
}

fit_one_common_gamma <- function(scenario_id, rep_id, method = "M1",
                                 data_dir = "data", output_dir = "output_common_gamma",
                                 verbose = 1000L) {
    rr <- sprintf("%02d", rep_id)
    dat_file <- file.path(data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    if (!file.exists(dat_file)) stop("Data file not found: ", dat_file)
    dat <- readRDS(dat_file)

    cat(sprintf("=== Fitting %s (common gamma) on %s rep %s ===\n", method, scenario_id, rr))

    settings  <- method_settings(method)
    priors    <- MCMC_PRIORS
    constants <- list(A0 = A0, B0 = B0, C0 = C0)

    BHB <- crossprod(B_ICAR, H %*% B_ICAR)
    spatial <- list(H = H, B_ICAR = B_ICAR, BHB = BHB)

    result <- run_mcmc_common_gamma(dat, settings, priors, spatial, constants,
                                    verbose = verbose)

    out_dir <- file.path(output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0("fit_", method, "_common_gamma_rep", rr, ".rds"))
    saveRDS(result, out_file)
    cat(sprintf("Saved: %s\n\n", out_file))

    result
}

fit_all_methods_common_gamma <- function(scenario_id, rep_id,
                                         methods = c("M1", "M2", "M3", "M4"), ...) {
    fits <- list()
    for (m in methods) {
        fits[[m]] <- fit_one_common_gamma(scenario_id, rep_id, method = m, ...)
    }
    fits
}
