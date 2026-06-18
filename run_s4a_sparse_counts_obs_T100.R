## ============================================================================
## run_s4a_sparse_counts_obs_T100_v3.R
##
## Scenario 4A sparse-count observation-level stress test.
## v3 adds a local boundary-safe gamma MH update for sparse-count fits.
##
## Purpose
##   Fit the same dynamic learned-gamma MSSTNB model used in Scenario 3 to
##   S4A sparse-count observation-level data.
##
## Design
##   S4A keeps the Scenario 3 latent spatiotemporal structure fixed and lowers
##   only the observation-level intensity.  This creates sparse counts without
##   forcing the latent lambda evolution to collapse through the y-driven
##   filtering recursion.
##
## Main comparison
##   Scenario 3 : dynamic DGP, dynamic fit, gamma learned by MCMC.
##   Scenario 4A: same model fit, but sparse observation-level counts.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   data_s4a_sparse_counts/S4A_SPARSE_COUNTS_OBS_T100/data_repXX.rds
##
## Main outputs
##   output_s4a_sparse_counts/S4A_SPARSE_COUNTS_OBS_T100/
##
## Notes
##   This script only runs fitting and writes fit manifests.  Posterior
##   performance and extra MCMC diagnostics should be handled by dedicated S4A
##   analysis scripts so that output table and figure names are scenario-specific.
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4A_SPARSE_COUNTS_OBS_T100"
data_scenario_id <- "S4A_SPARSE_COUNTS_OBS_T100"

## Use "short_test" first.  After rep 1 succeeds, change to "formal" for all
## 10 S4A replicates using the same MCMC profile as Scenario 3 formal runs.
run_profile <- "short_test"

if (identical(run_profile, "short_test")) {
    reps_formal <- 1:1
} else if (identical(run_profile, "formal")) {
    reps_formal <- 1:10
} else {
    stop("Unknown run_profile: ", run_profile, call. = FALSE)
}

## Core Scenario 3 model code.  S4A intentionally reuses the same learned-gamma
## sampler so the only substantive change is the data-generating regime.
s3_core_file <- "s3_dynamic_learned_gamma.R"

## Data-generation constants expected in the official S4A data.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "observation_sparse_counts"
expected_sparse_beta0_shift <- -4.25
expected_beta0_reference_truth <- -1.5
expected_beta0_sparse_truth <- -5.75
expected_count_multiplier <- exp(expected_sparse_beta0_shift)

## MCMC settings for learned-gamma fitting.  These mirror Scenario 3.
if (identical(run_profile, "short_test")) {
    n_iter <- 6000L
    n_burnin <- 1000L
    n_thin <- 5L
} else if (identical(run_profile, "formal")) {
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
}

gamma_prior <- c(1, 1)

## Output settings.
data_dir <- "data_s4a_sparse_counts"
output_dir <- "output_s4a_sparse_counts"
verbose <- 1000L

## In short_test mode, overwrite by default.  In formal mode, preserve existing
## fits unless you explicitly set overwrite_fit <- TRUE.
overwrite_fit <- identical(run_profile, "short_test")

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir_s4a <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s4a <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true_s4a <- function(x, message) {
    if (!isTRUE(x)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

s4a_data_file <- function(rep_id,
                          root = root_dir,
                          data_dir = data_dir,
                          data_scenario_id = data_scenario_id) {
    file.path(
        root,
        data_dir,
        data_scenario_id,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

s4a_fit_file <- function(rep_id,
                         root = root_dir,
                         output_dir = output_dir,
                         scenario_id = scenario_id) {
    file.path(
        root,
        output_dir,
        scenario_id,
        sprintf("fit_S4A_sparse_counts_obs_rep%02d.rds", as.integer(rep_id))
    )
}

s4a_summary_file <- function(rep_id,
                             root = root_dir,
                             output_dir = output_dir,
                             scenario_id = scenario_id) {
    file.path(
        root,
        output_dir,
        scenario_id,
        sprintf("summary_S4A_sparse_counts_obs_rep%02d.csv", as.integer(rep_id))
    )
}

check_s4a_source_dataset <- function(data_file,
                                     TT_use = 100L,
                                     n1_use = 9L,
                                     expected_stress_type = "observation_sparse_counts",
                                     expected_sparse_beta0_shift = -4.25,
                                     expected_beta0_reference_truth = -1.5,
                                     expected_beta0_sparse_truth = -5.75,
                                     beta0_ident_abs_limit = 20) {
    assert_file_exists_s4a(data_file, "S4A source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "y_coarse", "e", "x1", "x2",
        "lambda_tilde", "lambda_tilde_ident", "gamma_star", "beta0_star",
        "beta0_star_ident", "beta_star_ident", "phi_star_ident",
        "sparse_beta0_shift", "beta0_reference_truth", "beta0_sparse_truth",
        "expected_count_multiplier", "mean_count", "zero_prop", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s4a(
        length(missing) == 0L,
        paste("Dataset is missing required S4A fields:", paste(missing, collapse = ", "))
    )

    assert_true_s4a(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0(
            "y_coarse dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$y_coarse), collapse = " x "), "."
        )
    )

    for (nm in c("e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        assert_true_s4a(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    assert_true_s4a(
        identical(dat$stress_type, expected_stress_type),
        paste0("stress_type is not ", expected_stress_type, ". Got ", dat$stress_type, ".")
    )
    assert_true_s4a(
        abs(dat$sparse_beta0_shift - expected_sparse_beta0_shift) < 1e-12,
        paste0("sparse_beta0_shift is not ", expected_sparse_beta0_shift, ". Got ", dat$sparse_beta0_shift, ".")
    )
    assert_true_s4a(
        abs(dat$beta0_reference_truth - expected_beta0_reference_truth) < 1e-12,
        paste0("beta0_reference_truth is not ", expected_beta0_reference_truth, ". Got ", dat$beta0_reference_truth, ".")
    )
    assert_true_s4a(
        abs(dat$beta0_sparse_truth - expected_beta0_sparse_truth) < 1e-12,
        paste0("beta0_sparse_truth is not ", expected_beta0_sparse_truth, ". Got ", dat$beta0_sparse_truth, ".")
    )
    assert_true_s4a(
        abs(dat$expected_count_multiplier - exp(expected_sparse_beta0_shift)) < 1e-10,
        paste0("expected_count_multiplier is inconsistent with sparse_beta0_shift.")
    )

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0)) {
        stop("dat$e must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0) ||
        any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("lambda_tilde and lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }

    assert_true_s4a(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0(
            "beta0_star_ident appears pathological: ", dat$beta0_star_ident,
            ". Check whether this is the clean observation-level S4A data, not the older dynamic sparse data."
        )
    )

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat$lambda_tilde_ident, finite = TRUE)

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_ident_min = lambda_ident_range[[1L]],
        lambda_ident_max = lambda_ident_range[[2L]]
    )
}

make_s4a_source_data_manifest <- function(reps,
                                          root = root_dir,
                                          data_dir = data_dir,
                                          data_scenario_id = data_scenario_id) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s4a_data_file(
            rep_id = rep_id,
            root = root,
            data_dir = data_dir,
            data_scenario_id = data_scenario_id
        )
        chk <- check_s4a_source_dataset(data_file)
        data.frame(
            scenario_id = data_scenario_id,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            mean_count = chk$mean_count,
            median_count = chk$median_count,
            zero_prop = chk$zero_prop,
            total_count = chk$total_count,
            max_count = chk$max_count,
            beta0_star_ident = chk$beta0_star_ident,
            lambda_raw_min = chk$lambda_raw_min,
            lambda_raw_max = chk$lambda_raw_max,
            lambda_ident_min = chk$lambda_ident_min,
            lambda_ident_max = chk$lambda_ident_max,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

## ---- S4A fitting functions --------------------------------------------------
fit_s4a_sparse_counts_obs_one_rep <- function(rep_id,
                                              scenario_id = "S4A_SPARSE_COUNTS_OBS_T100",
                                              data_scenario_id = scenario_id,
                                              data_dir = "data_s4a_sparse_counts",
                                              output_dir = "output_s4a_sparse_counts",
                                              root = ".",
                                              settings_override = list(),
                                              priors = MCMC_PRIORS,
                                              spatial = build_s3_spatial(),
                                              gamma_init = NULL,
                                              gamma_prior = c(1, 1),
                                              verbose = 1000L,
                                              save_result = TRUE,
                                              return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, data_scenario_id, paste0("data_rep", rr, ".rds"))
    chk <- check_s4a_source_dataset(dat_file)
    dat <- chk$dat

    ## Use Scenario 3 validators and sampler.  This is intentional: S4A changes
    ## the data regime, not the fitted model.
    validate_s3_data(dat)

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    gamma_start <- gamma_init %||% dat$gamma_star %||% 0.8

    cat(sprintf("=== Scenario 4A sparse-count learned-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Stress    : %s; sparse_beta0_shift = %.2f\n", dat$stress_type, dat$sparse_beta0_shift))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Initial   : gamma = %.3f on average\n", mean(gamma_start)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, gamma, lambda\n")
    cat("Disabled  : delta, omega updates\n\n")

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = verbose
    )

    summary <- summarise_s3_dynamic_learned_gamma_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id)
    )

    ## Add S4A-specific data-level context to the one-row fit summary.
    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$sparse_beta0_shift <- dat$sparse_beta0_shift %||% NA_real_
    summary$expected_count_multiplier <- dat$expected_count_multiplier %||% NA_real_
    summary$reference_mean_count <- dat$reference_mean_count %||% NA_real_
    summary$reference_zero_prop <- dat$reference_zero_prop %||% NA_real_
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            stress_type = dat$stress_type,
            sparse_beta0_shift = dat$sparse_beta0_shift,
            beta0_reference_truth = dat$beta0_reference_truth,
            beta0_sparse_truth = dat$beta0_sparse_truth,
            expected_count_multiplier = dat$expected_count_multiplier,
            model_source = "S3_DYNAMIC_LEARNED_GAMMA",
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            fit_file_prefix = "fit_S4A_sparse_counts_obs_rep",
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        ensure_dir_s4a(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4A_sparse_counts_obs_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4A_sparse_counts_obs_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) {
        return(fit)
    }
    invisible(NULL)
}

fit_s4a_sparse_counts_obs_batch <- function(reps = 1:10,
                                            scenario_id = "S4A_SPARSE_COUNTS_OBS_T100",
                                            data_scenario_id = scenario_id,
                                            data_dir = "data_s4a_sparse_counts",
                                            output_dir = "output_s4a_sparse_counts",
                                            root = ".",
                                            settings_override = list(),
                                            gamma_init = NULL,
                                            gamma_prior = c(1, 1),
                                            verbose = 1000L,
                                            overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    ensure_dir_s4a(out_dir)

    summaries <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4A_sparse_counts_obs_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s4a_sparse_counts_obs_one_rep(
            rep_id = rep_id,
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            root = root,
            settings_override = settings_override,
            gamma_init = gamma_init,
            gamma_prior = gamma_prior,
            verbose = verbose,
            save_result = TRUE,
            return_result = TRUE
        )
        summaries[[rr]] <- fit$summary
    }

    summary_all <- do.call(rbind, summaries)
    summary_file <- file.path(out_dir, "summary_S4A_sparse_counts_obs_all_reps.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)

    invisible(summary_all)
}

print_s4a_fit_driver_summary <- function(source_manifest, fit_summary, fit_manifest) {
    cat("\n=== Scenario 4A sparse-count observation-level fitting summary ===\n")
    cat("Number of reps: ", nrow(source_manifest), "\n", sep = "")
    cat("Mean count avg: ", round(mean(source_manifest$mean_count, na.rm = TRUE), 4), "\n", sep = "")
    cat("Zero prop avg : ", round(mean(source_manifest$zero_prop, na.rm = TRUE), 4), "\n", sep = "")
    cat("Fit files all present: ", all(fit_manifest$fit_exists), "\n", sep = "")

    if (!is.null(fit_summary) && nrow(fit_summary) > 0L) {
        beta_cols <- intersect(c("beta0_mean", "beta1_mean", "beta2_mean", "gamma_common_mean", "gamma_mean"), names(fit_summary))
        if (length(beta_cols) > 0L) {
            cat("\nSelected fit-summary columns:\n")
            print(fit_summary[, c("scenario_id", "rep_id", beta_cols), drop = FALSE])
        }
    }
    invisible(TRUE)
}

## ---- pre-flight checks ------------------------------------------------------
assert_file_exists_s4a(file.path(root_dir, s3_core_file), "Scenario 3 core script")

source(file.path(root_dir, s3_core_file))
source_s3_dynamic_learned_gamma(root = root_dir)


## ---- S4A robust gamma update override --------------------------------------
## In sparse-count S4A fits, the current gamma can lie numerically on the
## transition-support boundary implied by the current lambda path.  The original
## S3 update used `is.finite(log_alpha)` before accepting, which rejects the
## important case `log_alpha = Inf` that occurs when the current state has
## target -Inf but a downward proposal is valid.  The override below keeps the
## same complete-data conditional target as S3, but uses boundary-safe MH logic:
##   - finite log_alpha: ordinary MH decision;
##   - +Inf log_alpha : accept;
##   - -Inf/NaN        : reject.
## This is intentionally local to the S4A driver; it does not modify the S3
## source file.
use_s4a_robust_gamma_update <- TRUE

if (isTRUE(use_s4a_robust_gamma_update)) {
    update_gamma_common_s3 <- function(gamma_current,
                                       lambda_tilde,
                                       y_coarse,
                                       a0 = 10,
                                       gamma_prior = c(1, 1),
                                       proposal_sd = 0.15) {
        gamma_current <- as.numeric(gamma_current)
        gamma_common_current <- if (length(gamma_current) == 1L) {
            gamma_current
        } else {
            mean(gamma_current)
        }

        z_current <- .logit_s3(gamma_common_current)
        z_proposal <- stats::rnorm(1L, mean = z_current, sd = proposal_sd)
        gamma_common_proposal <- .inv_logit_s3(z_proposal)

        log_target_current <- .log_gamma_transition_target_s3(
            gamma_common = gamma_common_current,
            lambda_tilde = lambda_tilde,
            y_coarse = y_coarse,
            a0 = a0,
            gamma_prior = gamma_prior
        )
        log_target_proposal <- .log_gamma_transition_target_s3(
            gamma_common = gamma_common_proposal,
            lambda_tilde = lambda_tilde,
            y_coarse = y_coarse,
            a0 = a0,
            gamma_prior = gamma_prior
        )

        log_jac_current <- log(gamma_common_current) + log1p(-gamma_common_current)
        log_jac_proposal <- log(gamma_common_proposal) + log1p(-gamma_common_proposal)

        target_current_z <- log_target_current + log_jac_current
        target_proposal_z <- log_target_proposal + log_jac_proposal
        log_alpha <- target_proposal_z - target_current_z

        ## Boundary-safe MH decision.
        ## Important: log_alpha = +Inf means proposal has positive target while
        ## current is outside/at numerical boundary; it should be accepted.
        if (is.nan(log_alpha) || is.na(log_alpha)) {
            accept <- FALSE
        } else if (is.infinite(log_alpha) && log_alpha > 0) {
            accept <- TRUE
        } else if (is.infinite(log_alpha) && log_alpha < 0) {
            accept <- FALSE
        } else {
            accept <- log(stats::runif(1L)) < min(0, log_alpha)
        }

        gamma_common_new <- if (accept) gamma_common_proposal else gamma_common_current
        gamma_new <- rep(gamma_common_new, ncol(lambda_tilde))

        list(
            gamma = gamma_new,
            gamma_common = gamma_common_new,
            accept = accept,
            log_alpha = log_alpha,
            proposal_sd = proposal_sd,
            log_target_current = log_target_current,
            log_target_proposal = log_target_proposal
        )
    }

    cat("Using S4A robust gamma update override with boundary-safe MH acceptance.\n")
}

source_data_manifest <- make_s4a_source_data_manifest(
    reps = reps_formal,
    root = root_dir,
    data_dir = data_dir,
    data_scenario_id = data_scenario_id
)

cat("\n=== S4A source-data check ===\n")
print(source_data_manifest[, c(
    "rep_id", "mean_count", "median_count", "zero_prop", "total_count",
    "max_count", "beta0_star_ident", "lambda_raw_min", "lambda_raw_max"
), drop = FALSE])

## ---- run Scenario 4A learned-gamma fit -------------------------------------
fit_summary <- fit_s4a_sparse_counts_obs_batch(
    reps = reps_formal,
    scenario_id = scenario_id,
    data_scenario_id = data_scenario_id,
    data_dir = data_dir,
    output_dir = output_dir,
    root = root_dir,
    settings_override = list(
        n_iter = n_iter,
        n_burnin = n_burnin,
        n_thin = n_thin
    ),
    gamma_prior = gamma_prior,
    verbose = verbose,
    overwrite_existing = overwrite_fit
)

fit_dir_full <- file.path(root_dir, output_dir, scenario_id)
fit_files <- file.path(
    fit_dir_full,
    sprintf("fit_S4A_sparse_counts_obs_rep%02d.rds", as.integer(reps_formal))
)
fit_manifest <- data.frame(
    scenario_id = scenario_id,
    data_scenario_id = data_scenario_id,
    rep_id = as.integer(reps_formal),
    data_file = source_data_manifest$data_file,
    fit_file = fit_files,
    fit_exists = file.exists(fit_files),
    stringsAsFactors = FALSE
)
assert_true_s4a(all(fit_manifest$fit_exists), "At least one S4A fit file was not created.")

ensure_dir_s4a(fit_dir_full)
utils::write.csv(
    source_data_manifest,
    file.path(fit_dir_full, "s4a_source_data_manifest.csv"),
    row.names = FALSE
)
utils::write.csv(
    fit_manifest,
    file.path(fit_dir_full, "s4a_fit_manifest.csv"),
    row.names = FALSE
)
saveRDS(
    list(
        run_profile = run_profile,
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        reps = reps_formal,
        mcmc = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin),
        gamma_prior = gamma_prior,
        source_data_manifest = source_data_manifest,
        fit_manifest = fit_manifest,
        fit_summary = fit_summary
    ),
    file.path(fit_dir_full, "run_info_S4A_sparse_counts_obs_T100.rds")
)

## ---- final summary ----------------------------------------------------------
print_s4a_fit_driver_summary(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
)

cat("\n=== Main output locations ===\n")
cat("S4A data: ", file.path(data_dir, data_scenario_id), "\n", sep = "")
cat("Fits    : ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("\nScenario 4A sparse-count observation-level T100 fitting finished successfully.\n")

invisible(list(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
))
