## ============================================================================
## run_single_fit_revised.R
## Driver script for the revised MSSTNB MCMC implementation
##
## Expected working directory:
##   the project/simulation directory containing R/ and data/
##
## Main changes relative to the original run_single_fit.R:
##   1. Sources the revised MCMC block files.
##   2. Uses marginal-NB r update by default.
##   3. Adds input/path checks and cleaner settings overrides.
##   4. Keeps fit_one(), fit_all_methods(), and sanity_check_mcmc().
## ============================================================================


## ---- small utilities --------------------------------------------------------
.require_file <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked <- function(path, verbose = TRUE) {
    .require_file(path)
    if (isTRUE(verbose)) {
        message("source: ", path)
    }
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}


## ---- source all revised dependencies ---------------------------------------
source_mcmc_revised <- function(root = ".", verbose = TRUE) {
    r_dir <- file.path(root, "R")
    mcmc_dir <- file.path(r_dir, "mcmc")

    files <- c(
        file.path(r_dir, "00_setup.R"),
        file.path(r_dir, "01_helpers.R"),
        file.path(mcmc_dir, "mcmc_config.R"),
        file.path(mcmc_dir, "mcmc_utils.R"),
        file.path(mcmc_dir, "update_kappa.R"),
        file.path(mcmc_dir, "update_regression_revised.R"),
        file.path(mcmc_dir, "update_icar_revised.R"),
        file.path(mcmc_dir, "update_dispersion_revised.R"),
        file.path(mcmc_dir, "update_gamma_revised.R"),
        file.path(mcmc_dir, "ffbs_lambda_revised.R"),
        file.path(mcmc_dir, "update_delta_revised.R"),
        file.path(mcmc_dir, "smooth_omega_revised.R"),
        file.path(mcmc_dir, "sampler_revised.R")
    )

    invisible(lapply(files, .source_checked, verbose = verbose))

    needed_objects <- c(
        "MCMC_SETTINGS", "MCMC_PRIORS", "A0", "B0", "C0",
        "H", "B_ICAR", "run_mcmc"
    )
    missing_objects <- needed_objects[!vapply(needed_objects, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing_objects) > 0L) {
        stop(
            "After sourcing dependencies, these required objects are missing: ",
            paste(missing_objects, collapse = ", "),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

## Backward-compatible alias if you are used to calling source_mcmc().
source_mcmc <- source_mcmc_revised


## ---- model settings ---------------------------------------------------------
method_settings <- function(method,
                            base_settings = MCMC_SETTINGS,
                            r_update_method = "marginal_nb",
                            overrides = list()) {
    method <- toupper(method)
    s <- base_settings

    ## Revised default: use marginal NB update for r when include_nb = TRUE.
    s$r_update_method <- r_update_method

    if (method == "M1") {
        s$include_nb <- TRUE
        s$include_icar <- TRUE
        s$include_covariates <- TRUE
    } else if (method == "M2") {
        s$include_nb <- FALSE
        s$include_icar <- TRUE
        s$include_covariates <- TRUE
    } else if (method == "M3") {
        s$include_nb <- TRUE
        s$include_icar <- FALSE
        s$include_covariates <- TRUE
    } else if (method == "M4") {
        s$include_nb <- FALSE
        s$include_icar <- FALSE
        s$include_covariates <- FALSE
    } else {
        stop("Unknown method: ", method, call. = FALSE)
    }

    ## If NB is disabled, r/kappa updates should be skipped inside sampler.
    ## The sampler should set r = Inf and kappa = 1 in this case.
    if (!isTRUE(s$include_nb)) {
        s$r_update_method <- "none"
    }

    if (length(overrides) > 0L) {
        for (nm in names(overrides)) {
            s[[nm]] <- overrides[[nm]]
        }
    }

    s$method <- method
    return(s)
}


## ---- data and object validation --------------------------------------------
validate_fit_inputs <- function(dat, settings, spatial, constants) {
    required_dat <- c("y_coarse", "y_fine", "e", "x1", "x2")
    missing_dat <- required_dat[!vapply(required_dat, function(z) !is.null(dat[[z]]), logical(1))]
    if (length(missing_dat) > 0L) {
        stop("dat is missing required fields: ", paste(missing_dat, collapse = ", "), call. = FALSE)
    }

    y_dim <- dim(dat$y_coarse)
    if (length(y_dim) != 2L) {
        stop("dat$y_coarse must be a matrix.", call. = FALSE)
    }
    if (!all(dim(dat$e) == y_dim) || !all(dim(dat$x1) == y_dim) || !all(dim(dat$x2) == y_dim)) {
        stop("dat$e, dat$x1, and dat$x2 must have the same dimensions as dat$y_coarse.", call. = FALSE)
    }
    if (length(dim(dat$y_fine)) != 3L) {
        stop("dat$y_fine must be a T x n1 x K array.", call. = FALSE)
    }
    if (!all(dim(dat$y_fine)[1:2] == y_dim)) {
        stop("The first two dimensions of dat$y_fine must match dim(dat$y_coarse).", call. = FALSE)
    }
    if (!is.null(dat$n1) && dat$n1 != y_dim[2]) {
        stop("dat$n1 does not match ncol(dat$y_coarse).", call. = FALSE)
    }

    if (isTRUE(settings$include_icar)) {
        if (is.null(spatial$H) || is.null(spatial$B_ICAR) || is.null(spatial$BHB)) {
            stop("settings$include_icar is TRUE, but spatial$H, spatial$B_ICAR, or spatial$BHB is missing.", call. = FALSE)
        }
        if (nrow(spatial$H) != y_dim[2] || ncol(spatial$H) != y_dim[2]) {
            stop("spatial$H must be n1 x n1.", call. = FALSE)
        }
        if (nrow(spatial$B_ICAR) != y_dim[2]) {
            stop("spatial$B_ICAR must have n1 rows.", call. = FALSE)
        }
    }

    if (is.null(constants$A0) || is.null(constants$B0) || is.null(constants$C0)) {
        stop("constants must contain A0, B0, and C0.", call. = FALSE)
    }
    if (constants$A0 <= 0 || constants$B0 <= 0) {
        stop("constants$A0 and constants$B0 must be positive.", call. = FALSE)
    }
    if (length(constants$C0) != dim(dat$y_fine)[3] || any(constants$C0 <= 0)) {
        stop("constants$C0 must be a positive vector of length K = dim(dat$y_fine)[3].", call. = FALSE)
    }

    invisible(TRUE)
}


## ---- spatial object builder -------------------------------------------------
build_spatial_objects <- function(H_obj = H, B_obj = B_ICAR) {
    if (is.null(H_obj) || is.null(B_obj)) {
        return(list(H = NULL, B_ICAR = NULL, BHB = NULL))
    }
    BHB <- crossprod(B_obj, H_obj %*% B_obj)
    list(H = H_obj, B_ICAR = B_obj, BHB = BHB)
}


## ---- fit one method on one dataset -----------------------------------------
fit_one <- function(scenario_id,
                    rep_id,
                    method = "M1",
                    data_dir = "data",
                    output_dir = "output_final_revised",
                    verbose = 1000L,
                    root = ".",
                    settings_override = list(),
                    priors = MCMC_PRIORS,
                    constants = list(A0 = A0, B0 = B0, C0 = C0),
                    spatial = NULL,
                    save_result = TRUE,
                    return_result = TRUE) {

    method <- toupper(method)
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    .require_file(dat_file)

    dat <- readRDS(dat_file)

    settings <- method_settings(
        method = method,
        base_settings = MCMC_SETTINGS,
        r_update_method = "marginal_nb",
        overrides = settings_override
    )

    if (is.null(spatial)) {
        spatial <- build_spatial_objects()
    }

    validate_fit_inputs(dat, settings, spatial, constants)

    cat(sprintf("=== Revised MCMC: fitting %s on %s rep %s ===\n", method, scenario_id, rr))
    cat(sprintf("Data file       : %s\n", dat_file))
    cat(sprintf("NB              : %s\n", isTRUE(settings$include_nb)))
    cat(sprintf("ICAR            : %s\n", isTRUE(settings$include_icar)))
    cat(sprintf("Covariates      : %s\n", isTRUE(settings$include_covariates)))
    cat(sprintf("r update method : %s\n\n", settings$r_update_method %||% "NULL"))

    result <- run_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        constants = constants,
        verbose = verbose
    )

    ## Add reproducibility metadata without disturbing sampler output.
    result$metadata <- c(result$metadata %||% list(), list(
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id),
        method = method,
        data_file = dat_file,
        output_dir = output_dir,
        r_update_method = settings$r_update_method,
        run_time = Sys.time()
    ))

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(out_dir, paste0("fit_", method, "_rep", rr, ".rds"))
        saveRDS(result, out_file)
        cat(sprintf("Saved: %s\n\n", out_file))
        result$metadata$output_file <- out_file
    }

    if (isTRUE(return_result)) {
        return(result)
    }

    invisible(NULL)
}


## ---- fit multiple methods ---------------------------------------------------
fit_all_methods <- function(scenario_id,
                            rep_id,
                            methods = c("M1", "M2", "M3", "M4"),
                            ...) {
    fits <- list()
    for (m in methods) {
        fits[[toupper(m)]] <- fit_one(
            scenario_id = scenario_id,
            rep_id = rep_id,
            method = m,
            ...
        )
    }
    return(fits)
}


## ---- concise posterior summary ---------------------------------------------
print_quick_summary <- function(result) {
    if (is.null(result$samples)) {
        warning("result$samples is missing; cannot print posterior summary.")
        return(invisible(NULL))
    }

    s <- result$samples
    cat("\n--- Quick posterior check ---\n")

    if (!is.null(s$beta0)) {
        cat(sprintf("  beta0   : mean=% .3f  sd=%.3f\n", mean(s$beta0), sd(s$beta0)))
    }
    if (!is.null(s$beta) && ncol(s$beta) >= 1L) {
        cat(sprintf("  beta1   : mean=% .3f  sd=%.3f\n", mean(s$beta[, 1]), sd(s$beta[, 1])))
    }
    if (!is.null(s$beta) && ncol(s$beta) >= 2L) {
        cat(sprintf("  beta2   : mean=% .3f  sd=%.3f\n", mean(s$beta[, 2]), sd(s$beta[, 2])))
    }
    if (!is.null(s$gamma)) {
        cat(sprintf("  gamma[1]: mean=% .3f  sd=%.3f\n", mean(s$gamma[, 1]), sd(s$gamma[, 1])))
    }
    if (!is.null(s$r) && all(is.finite(s$r[, 1]))) {
        cat(sprintf("  r[1]    : mean=% .3f  sd=%.3f\n", mean(s$r[, 1]), sd(s$r[, 1])))
    }
    if (!is.null(s$delta)) {
        cat(sprintf("  delta   : mean=% .3f  sd=%.3f\n", mean(s$delta), sd(s$delta)))
    }
    if (!is.null(s$loglik)) {
        cat(sprintf("  loglik  : final=% .1f  mean=% .1f\n", tail(s$loglik, 1), mean(s$loglik)))
    }

    invisible(NULL)
}


## ---- quick sanity check ------------------------------------------------------
sanity_check_mcmc <- function(scenario_id = "S1",
                              rep_id = 1L,
                              method = "M1",
                              n_iter = 2000L,
                              n_burnin = 1000L,
                              n_thin = 2L,
                              verbose = 500L,
                              ...) {
    cat("=== Revised MCMC Sanity Check ===\n")
    cat(sprintf(
        "Fitting %s on %s rep%02d, %d iter, %d burn-in, thin=%d.\n\n",
        toupper(method), scenario_id, as.integer(rep_id),
        as.integer(n_iter), as.integer(n_burnin), as.integer(n_thin)
    ))

    overrides <- list(n_iter = as.integer(n_iter),
                      n_burnin = as.integer(n_burnin),
                      n_thin = as.integer(n_thin))

    result <- fit_one(
        scenario_id = scenario_id,
        rep_id = rep_id,
        method = method,
        verbose = verbose,
        settings_override = overrides,
        ...
    )

    print_quick_summary(result)
    return(result)
}


## ---- optional command-line interface ---------------------------------------
## Example:
##   Rscript run_single_fit_revised.R S1 1 M1 data output_final_revised
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L && identical(environment(), globalenv())) {
    scenario_id <- args[[1]]
    rep_id <- as.integer(args[[2]] %||% 1L)
    method <- args[[3]] %||% "M1"
    data_dir <- args[[4]] %||% "data"
    output_dir <- args[[5]] %||% "output_final_revised"

    source_mcmc_revised(verbose = TRUE)
    fit_one(
        scenario_id = scenario_id,
        rep_id = rep_id,
        method = method,
        data_dir = data_dir,
        output_dir = output_dir,
        verbose = 1000L
    )
}
