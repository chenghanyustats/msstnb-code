
## ======================================================================
## run_single_fit_fixed_common_gamma.R
## Source dependencies, fit one dataset with a fixed common gamma,
## or compare a small grid such as gamma = 0.7, 0.8, 0.9.
## ======================================================================

source_mcmc_fixed_common_gamma <- function() {
    source("R/00_setup.R")
    source("R/01_helpers.R")
    source("R/mcmc/mcmc_config.R")
    source("R/mcmc/mcmc_utils.R")
    source("R/mcmc/update_kappa.R")
    source("R/mcmc/update_regression.R")
    source("R/mcmc/update_icar.R")
    source("R/mcmc/update_dispersion.R")
    source("R/mcmc/ffbs_lambda.R")
    source("R/mcmc/update_delta.R")
    source("R/mcmc/smooth_omega.R")
    source("R/mcmc/sampler_fixed_common_gamma.R")
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

fit_one_fixed_common_gamma <- function(scenario_id, rep_id, fixed_gamma,
                                       method = "M1",
                                       data_dir = "data",
                                       output_dir = "output_fixed_common_gamma",
                                       verbose = 1000L) {
    rr <- sprintf("%02d", rep_id)
    dat_file <- file.path(data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    if (!file.exists(dat_file)) stop("Data file not found: ", dat_file)
    dat <- readRDS(dat_file)

    cat(sprintf("=== Fitting %s (fixed common gamma = %.3f) on %s rep %s ===\n",
                method, fixed_gamma, scenario_id, rr))

    settings  <- method_settings(method)
    priors    <- MCMC_PRIORS
    constants <- list(A0 = A0, B0 = B0, C0 = C0)

    BHB <- crossprod(B_ICAR, H %*% B_ICAR)
    spatial <- list(H = H, B_ICAR = B_ICAR, BHB = BHB)

    result <- run_mcmc_fixed_common_gamma(dat, settings, priors, spatial,
                                          constants, fixed_gamma = fixed_gamma,
                                          verbose = verbose)

    out_dir <- file.path(output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    gtag <- gsub("\\.", "p", sprintf("%.2f", fixed_gamma))
    out_file <- file.path(out_dir,
                          paste0("fit_", method, "_fixed_common_gamma_", gtag,
                                 "_rep", rr, ".rds"))
    saveRDS(result, out_file)
    cat(sprintf("Saved: %s\n\n", out_file))

    result
}

summarise_fixed_common_gamma_fit <- function(fit, data_obj) {
    s <- fit$samples

    data.frame(
        fixed_gamma = unique(s$gamma)[1],
        mean_loglik = mean(s$loglik),
        beta0_mean = mean(s$beta0),
        beta1_mean = mean(s$beta[, 1]),
        beta2_mean = mean(s$beta[, 2]),
        delta_mean = mean(s$delta),
        gamma_truth = if (length(data_obj$gamma_star) == 1L) data_obj$gamma_star else data_obj$gamma_star[1],
        beta1_truth = data_obj$beta_star[1],
        beta2_truth = data_obj$beta_star[2],
        delta_truth = data_obj$delta_star,
        phi_accept_rate = fit$diagnostics$phi_accept_rate,
        delta_accept_rate = fit$diagnostics$delta_accept_rate,
        r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate),
        beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject,
        stringsAsFactors = FALSE
    )
}

compare_fixed_common_gamma_grid <- function(scenario_id, rep_id,
                                            gamma_grid = c(0.7, 0.8, 0.9),
                                            method = "M1",
                                            data_dir = "data",
                                            output_dir = "output_fixed_common_gamma",
                                            verbose = 1000L) {
    rr <- sprintf("%02d", rep_id)
    dat_file <- file.path(data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    if (!file.exists(dat_file)) stop("Data file not found: ", dat_file)
    data_obj <- readRDS(dat_file)

    fits <- vector("list", length(gamma_grid))
    names(fits) <- paste0("g", gamma_grid)
    tabs <- vector("list", length(gamma_grid))

    for (i in seq_along(gamma_grid)) {
        g <- gamma_grid[i]
        fits[[i]] <- fit_one_fixed_common_gamma(
            scenario_id = scenario_id,
            rep_id = rep_id,
            fixed_gamma = g,
            method = method,
            data_dir = data_dir,
            output_dir = output_dir,
            verbose = verbose
        )
        tabs[[i]] <- summarise_fixed_common_gamma_fit(fits[[i]], data_obj)
    }

    summary_tab <- do.call(rbind, tabs)
    rownames(summary_tab) <- NULL

    out_dir <- file.path(output_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    csv_file <- file.path(out_dir,
                          paste0("compare_fixed_common_gamma_", method, "_rep", rr, ".csv"))
    write.csv(summary_tab, csv_file, row.names = FALSE)
    cat(sprintf("Comparison summary saved: %s\n", csv_file))

    list(fits = fits, summary = summary_tab)
}
