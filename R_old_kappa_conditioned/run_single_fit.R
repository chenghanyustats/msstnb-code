## ==========================================================================
## run_single_fit.R
## Load data, configure, run one MCMC fit, and save results
##
## Usage (from the simulation/ directory):
##   source("R/mcmc/run_single_fit.R")
##
##   # Fit M1 (full model) on S1, replication 1:
##   fit <- fit_one("S1", rep_id = 1, method = "M1")
##
##   # Fit all 4 methods on S1, replication 1:
##   fits <- fit_all_methods("S1", rep_id = 1)
## ==========================================================================


## ---- source all dependencies ----
source_mcmc <- function() {
    ## DGP helpers (for rdirichlet, etc.)
    source("R/01_helpers.R")

    ## Setup (for constants, spatial structure)
    source("R/00_setup.R")

    ## MCMC components
    source("R/mcmc/mcmc_config.R")
    source("R/mcmc/mcmc_utils.R")
    source("R/mcmc/update_kappa.R")
    source("R/mcmc/update_regression.R")
    source("R/mcmc/update_icar.R")
    source("R/mcmc/update_dispersion.R")
    source("R/mcmc/update_gamma.R")
    source("R/mcmc/ffbs_lambda.R")
    source("R/mcmc/update_delta.R")
    source("R/mcmc/smooth_omega.R")
    source("R/mcmc/sampler.R")
}


#' Build method-specific MCMC settings
#'
#' @param method  One of "M1", "M2", "M3", "M4"
#' @return        Modified MCMC_SETTINGS list
method_settings <- function(method) {
    s <- MCMC_SETTINGS
    if (method == "M1") {
        s$include_nb         <- TRUE
        s$include_icar       <- TRUE
        s$include_covariates <- TRUE
    } else if (method == "M2") {
        s$include_nb         <- FALSE
        s$include_icar       <- TRUE
        s$include_covariates <- TRUE
    } else if (method == "M3") {
        s$include_nb         <- TRUE
        s$include_icar       <- FALSE
        s$include_covariates <- TRUE
    } else if (method == "M4") {
        s$include_nb         <- FALSE
        s$include_icar       <- FALSE
        s$include_covariates <- FALSE
    } else {
        stop("Unknown method: ", method)
    }
    return(s)
}


#' Fit one method on one dataset
#'
#' @param scenario_id  "S1"..."S5"
#' @param rep_id       Replication number (1–50)
#' @param method       "M1", "M2", "M3", or "M4"
#' @param data_dir     Base directory for data files
#' @param output_dir   Base directory for output files
#' @param verbose      Progress print interval (0 = silent)
#' @return             MCMC result list (also saved to disk)
fit_one <- function(scenario_id, rep_id, method = "M1",
                    data_dir = "data", output_dir = "output",
                    verbose = 1000L) {

    rr <- sprintf("%02d", rep_id)

    ## ---- load data ----
    dat_file <- file.path(data_dir, scenario_id,
                          paste0("data_rep", rr, ".rds"))
    if (!file.exists(dat_file)) {
        stop("Data file not found: ", dat_file)
    }
    dat <- readRDS(dat_file)

    cat(sprintf("=== Fitting %s on %s rep %s ===\n",
                method, scenario_id, rr))

    ## ---- configure ----
    settings  <- method_settings(method)
    priors    <- MCMC_PRIORS
    constants <- list(A0 = A0, B0 = B0, C0 = C0)

    ## Precompute Cholesky of B'HB (for ICAR ESS)
    BHB   <- crossprod(B_ICAR, H %*% B_ICAR)
    R_BHB <- chol(BHB)

    spatial <- list(
        H      = H,
        B_ICAR = B_ICAR,
        R_BHB  = R_BHB
    )

    ## ---- run MCMC ----
    result <- run_mcmc(dat, settings, priors, spatial, constants,
                       verbose = verbose)

    ## ---- save ----
    out_dir <- file.path(output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    out_file <- file.path(out_dir,
                          paste0("fit_", method, "_rep", rr, ".rds"))
    saveRDS(result, out_file)
    cat(sprintf("Saved: %s\n\n", out_file))

    return(result)
}


#' Fit all 4 methods on one dataset
#'
#' @param scenario_id  "S1"..."S5"
#' @param rep_id       Replication number (1–50)
#' @param methods      Character vector of methods to fit
#' @param ...          Passed to fit_one()
#' @return             Named list of MCMC result lists
fit_all_methods <- function(scenario_id, rep_id,
                            methods = c("M1", "M2", "M3", "M4"),
                            ...) {
    fits <- list()
    for (m in methods) {
        fits[[m]] <- fit_one(scenario_id, rep_id, method = m, ...)
    }
    return(fits)
}


#' Quick test: fit M1 on S1 rep 01 with reduced iterations
#'
#' @return MCMC result list
sanity_check_mcmc <- function() {
    cat("=== MCMC Sanity Check ===\n")
    cat("Fitting M1 on S1 rep01, 2000 iterations (1000 burn-in)...\n\n")

    ## Temporarily reduce iterations
    orig_settings <- MCMC_SETTINGS
    MCMC_SETTINGS$n_iter   <<- 2000L
    MCMC_SETTINGS$n_burnin <<- 1000L
    MCMC_SETTINGS$n_thin   <<- 2L

    result <- fit_one("S1", rep_id = 1, method = "M1", verbose = 500L)

    ## Restore original settings
    MCMC_SETTINGS <<- orig_settings

    ## Print quick diagnostics
    s <- result$samples
    cat("\n--- Quick posterior check ---\n")
    cat(sprintf("  beta0: mean=%.3f, sd=%.3f  (true=check data file)\n",
                mean(s$beta0), sd(s$beta0)))
    cat(sprintf("  beta1: mean=%.3f, sd=%.3f\n",
                mean(s$beta[, 1]), sd(s$beta[, 1])))
    cat(sprintf("  beta2: mean=%.3f, sd=%.3f\n",
                mean(s$beta[, 2]), sd(s$beta[, 2])))
    cat(sprintf("  gamma[1]: mean=%.3f, sd=%.3f\n",
                mean(s$gamma[, 1]), sd(s$gamma[, 1])))
    cat(sprintf("  r[1]: mean=%.3f, sd=%.3f\n",
                mean(s$r[, 1]), sd(s$r[, 1])))
    cat(sprintf("  delta: mean=%.3f, sd=%.3f\n",
                mean(s$delta), sd(s$delta)))

    return(result)
}
