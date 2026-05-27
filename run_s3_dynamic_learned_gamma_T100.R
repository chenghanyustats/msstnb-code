## ============================================================================
## run_s3_dynamic_learned_gamma_T100.R
##
## Formal Scenario 3 dynamic residual-risk validation with learned gamma.
##
## Purpose
##   Fit the full dynamic MSSTNB model to a shared T100 dynamic DGP dataset,
##   but learn the common gamma discount parameter inside the MCMC instead of
##   fixing it at the truth.
##
## Main comparison
##   Scenario 2: dynamic DGP, dynamic fit, gamma fixed at truth.
##   Scenario 3: dynamic DGP, dynamic fit, gamma learned by MCMC.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   run_scenario3_posterior_performance.R
##   run_scenario3_extra_mcmc_diagnostics.R
##
## Data policy
##   This driver does not modify or rerun Scenario 2. It reads shared dynamic
##   DGP data from data_revised/DGP_DYNAMIC_T100/. If that folder is missing,
##   it copies existing files from data_revised/S2_DYNAMIC_FIXED_GAMMA_T100/.
##
## Main outputs
##   output_s3_dynamic_learned_gamma/S3_DYNAMIC_LEARNED_GAMMA_T100/
##   analysis_s3_dynamic_learned_gamma/S3_DYNAMIC_LEARNED_GAMMA_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA_T100"

## Shared data folder used by Scenario 3.
## This avoids making S3 appear to depend on the Scenario 2 model.
source_data_scenario_id <- "DGP_DYNAMIC_T100"

## Existing dynamic data folder produced earlier by the Scenario 2 workflow.
## It is used only as a source to copy data into DGP_DYNAMIC_T100 if needed.
legacy_source_data_scenario_id <- "S2_DYNAMIC_FIXED_GAMMA_T100"

## Run profile. Use "short_test" first. After the whole workflow succeeds,
## change this to "formal" for the 20-replicate production run.
run_profile <- "formal"

if (identical(run_profile, "short_test")) {
    reps_formal <- 1:1
} else if (identical(run_profile, "formal")) {
    reps_formal <- 1:20
} else {
    stop("Unknown run_profile: ", run_profile, call. = FALSE)
}

## Core Scenario 3 scripts
s3_core_file <- "s3_dynamic_learned_gamma.R"
s3_perf_file <- "run_scenario3_posterior_performance.R"
s3_diag_file <- "run_scenario3_extra_mcmc_diagnostics.R"

## Do not rerun or modify Scenario 2 here. If shared DGP data are missing,
## copy the already generated Scenario 2 dynamic data into DGP_DYNAMIC_T100.
copy_legacy_s2_data_if_missing <- TRUE

## Data-generating settings. These must match Scenario 2 T100 exactly.
TT_use <- 100L
n1_use <- 9L
beta0_truth <- -1.5
beta_truth <- c(0.5, -0.4)
r_truth <- 15
tau_phi_truth <- 2
omega_mode <- "fixed_prior_mean"

gamma_truth <- 0.80
A0_use <- 10
B0_use <- 10

## Clean covariate design. This must match Scenario 2 T100.
x2_mode <- "continuous_time"
x2_ar <- 0.50
x2_innov_sd <- 0.80

## MCMC settings for learned-gamma Scenario 3.
## These are controlled only here through run_profile. There should be no
## second n_iter/n_burnin/n_thin assignment later in this file.
if (identical(run_profile, "short_test")) {
    n_iter <- 6000L
    n_burnin <- 1000L
    n_thin <- 5L
} else if (identical(run_profile, "formal")) {
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
}

## Gamma prior for the learned common gamma.
## Beta(1, 1) is intentionally weak. If posterior mass piles up near 1,
## diagnose before imposing a stronger prior.
gamma_prior <- c(1, 1)

## Output settings
data_dir <- "data_revised"
output_dir <- "output_s3_dynamic_learned_gamma"
analysis_dir <- "analysis_s3_dynamic_learned_gamma"
plot_format <- "png"
verbose <- 1000L

## Set overwrite_fit = FALSE to reuse existing Scenario 3 fits and regenerate
## summaries only. Set overwrite_data = FALSE unless you intentionally want to
## regenerate the matched dynamic datasets. In short_test mode, overwrite is
## TRUE by default so the test run is not accidentally skipped.
overwrite_data <- FALSE
overwrite_fit <- identical(run_profile, "short_test")

## ---- helper functions -------------------------------------------------------
assert_file_exists <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true <- function(x, message) {
    if (!isTRUE(x)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

scenario_data_file <- function(rep_id) {
    file.path(
        root_dir,
        data_dir,
        source_data_scenario_id,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

check_s3_source_dataset <- function(data_file, TT_use, n1_use,
                                    expected_x2_mode = "continuous_time") {
    assert_file_exists(data_file, "Scenario 3 source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "y_coarse", "x1", "x2", "lambda_tilde", "lambda_tilde_ident",
        "gamma_star", "beta_star", "phi_star"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true(
        length(missing) == 0L,
        paste("Dataset is missing required fields:", paste(missing, collapse = ", "))
    )

    assert_true(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0(
            "y_coarse dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$y_coarse), collapse = " x "), "."
        )
    )
    assert_true(
        identical(as.integer(dim(dat$x1)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0("x1 dimension is not TT_use by n1_use. Got ", paste(dim(dat$x1), collapse = " x "), ".")
    )
    assert_true(
        identical(as.integer(dim(dat$x2)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0("x2 dimension is not TT_use by n1_use. Got ", paste(dim(dat$x2), collapse = " x "), ".")
    )
    assert_true(
        identical(as.integer(dim(dat$lambda_tilde_ident)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0(
            "lambda_tilde_ident dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$lambda_tilde_ident), collapse = " x "), "."
        )
    )

    if (!is.null(dat$x2_mode)) {
        assert_true(
            identical(dat$x2_mode, expected_x2_mode),
            paste0("x2_mode is not ", expected_x2_mode, ". Got ", dat$x2_mode, ".")
        )
    }

    gamma_range <- range(dat$gamma_star, finite = TRUE)
    assert_true(
        all(abs(gamma_range - gamma_truth) < 1e-12),
        paste0("gamma_star is not fixed at gamma_truth. Range is ", paste(gamma_range, collapse = ", "), ".")
    )

    x2_unique_counts <- apply(dat$x2, 2, function(z) length(unique(round(z, 8))))
    assert_true(
        all(x2_unique_counts >= max(10L, floor(TT_use / 2L))),
        paste0(
            "x2 does not appear to be time-varying in every region. Unique counts: ",
            paste(x2_unique_counts, collapse = ", "), "."
        )
    )

    lambda_range <- range(dat$lambda_tilde_ident, finite = TRUE)
    assert_true(
        length(lambda_range) == 2L && is.finite(lambda_range[1]) &&
            is.finite(lambda_range[2]) && lambda_range[1] > 0 && diff(lambda_range) > 1e-8,
        paste0(
            "lambda_tilde_ident does not appear to be positive and dynamic. Range is ",
            paste(lambda_range, collapse = ", "), "."
        )
    )

    x1 <- as.vector(dat$x1)
    x2 <- as.vector(dat$x2)
    loglam <- as.vector(log(dat$lambda_tilde_ident))
    phi_mat <- matrix(dat$phi_star, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)

    list(
        dat = dat,
        x2_unique_counts = x2_unique_counts,
        x_cor = stats::cor(x1, x2, use = "complete.obs"),
        x1_sd = stats::sd(x1, na.rm = TRUE),
        x2_sd = stats::sd(x2, na.rm = TRUE),
        cor_x2_loglambda = stats::cor(x2, loglam, use = "complete.obs"),
        cor_x1_loglambda = stats::cor(x1, loglam, use = "complete.obs"),
        cor_x2_phi = stats::cor(x2, as.vector(phi_mat), use = "complete.obs"),
        cor_x1_phi = stats::cor(x1, as.vector(phi_mat), use = "complete.obs"),
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        lambda_truth_min = lambda_range[1],
        lambda_truth_max = lambda_range[2]
    )
}

make_s3_confounding_table <- function(reps, source_data_scenario_id, data_dir, root = ".") {
    out <- lapply(reps, function(rep_id) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat_file <- file.path(root, data_dir, source_data_scenario_id, paste0("data_rep", rr, ".rds"))
        dat <- readRDS(dat_file)

        x1 <- as.vector(dat$x1)
        x2 <- as.vector(dat$x2)
        loglam <- as.vector(log(dat$lambda_tilde_ident %||% dat$lambda_tilde))
        phi_mat <- matrix(dat$phi_star, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)

        data.frame(
            scenario_id = scenario_id,
            source_data_scenario_id = source_data_scenario_id,
            rep_id = as.integer(rep_id),
            x2_mode = dat$x2_mode %||% NA_character_,
            sd_x1 = stats::sd(x1),
            sd_x2 = stats::sd(x2),
            cor_x1_x2 = stats::cor(x1, x2, use = "complete.obs"),
            cor_x2_loglambda = stats::cor(x2, loglam, use = "complete.obs"),
            cor_x1_loglambda = stats::cor(x1, loglam, use = "complete.obs"),
            cor_x2_phi = stats::cor(x2, as.vector(phi_mat), use = "complete.obs"),
            cor_x1_phi = stats::cor(x1, as.vector(phi_mat), use = "complete.obs"),
            sd_beta2_x2 = stats::sd(dat$beta_star[2] * x2),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

print_s3_summary <- function(dataset_check, perf_out, diag_out = NULL) {
    cat("\n=== Scenario 3 dynamic learned-gamma validation checks ===\n")
    cat("Scenario: ", scenario_id, "\n", sep = "")
    cat("Source data scenario: ", source_data_scenario_id, "\n", sep = "")
    cat("Replicates: ", min(reps_formal), " to ", max(reps_formal), " (n = ", length(reps_formal), ")\n", sep = "")
    cat("TT_use: ", TT_use, "\n", sep = "")
    cat("n1_use: ", n1_use, "\n", sep = "")
    cat("lambda_tilde: dynamic in truth and sampled in fit\n")
    cat("gamma: learned inside MCMC; truth gamma = ", gamma_truth, "\n", sep = "")
    cat("gamma prior: Beta(", gamma_prior[1], ", ", gamma_prior[2], ")\n", sep = "")
    cat("x2 mode: ", x2_mode, ", AR = ", x2_ar, ", innovation SD = ", x2_innov_sd, "\n", sep = "")

    cat("\nData check from rep 1:\n")
    cat("  mean count: ", round(dataset_check$mean_count, 3), "\n", sep = "")
    cat("  zero proportion: ", round(dataset_check$zero_prop, 3), "\n", sep = "")
    cat("  lambda truth range: [", round(dataset_check$lambda_truth_min, 3), ", ",
        round(dataset_check$lambda_truth_max, 3), "]\n", sep = "")
    cat("  x2 unique counts by region: ", paste(dataset_check$x2_unique_counts, collapse = ", "), "\n", sep = "")
    cat("  cor(vec(x1), vec(x2)): ", round(dataset_check$x_cor, 4), "\n", sep = "")
    cat("  cor(vec(x2), vec(log lambda)): ", round(dataset_check$cor_x2_loglambda, 4), "\n", sep = "")
    cat("  cor(vec(x2), vec(phi)): ", round(dataset_check$cor_x2_phi, 4), "\n", sep = "")
    cat("  sd(vec(x1)): ", round(dataset_check$x1_sd, 4), "\n", sep = "")
    cat("  sd(vec(x2)): ", round(dataset_check$x2_sd, 4), "\n", sep = "")

    if (!is.null(perf_out$gamma_summary) && nrow(perf_out$gamma_summary) > 0L) {
        cat("\nGamma posterior summary across replicates:\n")
        gamma_sum <- perf_out$gamma_summary
        keep <- intersect(c("mean", "sd", "median", "q025", "q975", "bias", "covered_95"), names(gamma_sum))
        print(summary(gamma_sum[, keep, drop = FALSE]))
    }

    if (!is.null(perf_out$beta_summary) && nrow(perf_out$beta_summary) > 0L) {
        cat("\nRegression recovery across replicates:\n")
        beta_sum <- perf_out$beta_summary
        keep <- intersect(c("post_mean", "bias", "covered_95"), names(beta_sum))
        beta_agg <- aggregate(
            beta_sum[, keep, drop = FALSE],
            by = list(parameter = beta_sum$parameter),
            FUN = function(z) c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
        )
        print(beta_agg)
    }

    if (!is.null(perf_out$latent_summary) && nrow(perf_out$latent_summary) > 0L) {
        cat("\nLatent path recovery metrics:\n")
        print(perf_out$latent_summary)
    }

    if (!is.null(diag_out) && !is.null(diag_out$mcmc_diagnostics) && nrow(diag_out$mcmc_diagnostics) > 0L) {
        cat("\nExtra MCMC diagnostic summary:\n")
        diag <- diag_out$mcmc_diagnostics
        print(summary(diag[, intersect(c("acf1", "ess"), names(diag)), drop = FALSE]))
    }
}

## ---- required source files --------------------------------------------------
assert_file_exists(file.path(root_dir, s3_core_file), "Scenario 3 core script")
assert_file_exists(file.path(root_dir, s3_perf_file), "Scenario 3 posterior performance script")
assert_file_exists(file.path(root_dir, s3_diag_file), "Scenario 3 extra MCMC diagnostics script")

message("Scenario 3 T100 run profile: ", run_profile)
message("Replicates: ", paste(reps_formal, collapse = ", "))
message("MCMC: n_iter = ", n_iter, ", n_burnin = ", n_burnin, ", n_thin = ", n_thin)
message("overwrite_fit = ", overwrite_fit)

## ---- clean old output if requested -----------------------------------------
if (isTRUE(overwrite_fit)) {
    unlink(file.path(output_dir, scenario_id), recursive = TRUE)
    unlink(file.path(analysis_dir, scenario_id), recursive = TRUE)
}

## ---- create or verify shared DGP source data -------------------------------
source_data_files <- vapply(reps_formal, scenario_data_file, character(1L))
missing_data_files <- source_data_files[!file.exists(source_data_files)]

if (length(missing_data_files) > 0L) {
    if (!isTRUE(copy_legacy_s2_data_if_missing)) {
        stop(
            "The shared DGP T100 data files are missing. First missing file: ",
            missing_data_files[[1L]],
            call. = FALSE
        )
    }

    legacy_dir <- file.path(root_dir, data_dir, legacy_source_data_scenario_id)
    shared_dir <- file.path(root_dir, data_dir, source_data_scenario_id)

    assert_true(
        dir.exists(legacy_dir),
        paste0(
            "Shared DGP data are missing and legacy Scenario 2 data folder does not exist: ",
            legacy_dir,
            ". Run Scenario 2 data generation first, or manually create ",
            shared_dir,
            "."
        )
    )

    ensure_dir(shared_dir)

    legacy_files <- file.path(
        legacy_dir,
        sprintf("data_rep%02d.rds", as.integer(reps_formal))
    )
    missing_legacy <- legacy_files[!file.exists(legacy_files)]
    assert_true(
        length(missing_legacy) == 0L,
        paste0(
            "Some legacy Scenario 2 dynamic data files are missing. First missing file: ",
            missing_legacy[[1L]]
        )
    )

    message("Copying existing dynamic DGP data from ", legacy_dir, " to ", shared_dir)
    copied <- file.copy(
        from = legacy_files,
        to = shared_dir,
        overwrite = FALSE
    )
    assert_true(
        all(copied | file.exists(source_data_files)),
        "Failed to copy at least one shared DGP data file."
    )

    source_data_files <- vapply(reps_formal, scenario_data_file, character(1L))
    missing_data_files <- source_data_files[!file.exists(source_data_files)]
    assert_true(
        length(missing_data_files) == 0L,
        paste0("Shared DGP data are still missing after copy. First missing file: ", missing_data_files[[1L]])
    )
}

source_manifest <- data.frame(
    scenario_id = source_data_scenario_id,
    legacy_source_data_scenario_id = legacy_source_data_scenario_id,
    rep_id = as.integer(reps_formal),
    data_file = source_data_files,
    status = "existing_or_copied",
    stringsAsFactors = FALSE
)

## ---- run Scenario 3 learned gamma ------------------------------------------
source(file.path(root_dir, s3_core_file))
source_s3_dynamic_learned_gamma(root = root_dir)

fit_summary <- fit_s3_dynamic_learned_gamma_batch(
    reps = reps_formal,
    scenario_id = scenario_id,
    data_scenario_id = source_data_scenario_id,
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
    sprintf("fit_S3_dynamic_learned_gamma_rep%02d.rds", as.integer(reps_formal))
)
fit_manifest <- data.frame(
    scenario_id = scenario_id,
    source_data_scenario_id = source_data_scenario_id,
    rep_id = as.integer(reps_formal),
    data_file = source_data_files,
    fit_file = fit_files,
    fit_exists = file.exists(fit_files),
    stringsAsFactors = FALSE
)
assert_true(all(fit_manifest$fit_exists), "At least one Scenario 3 fit file was not created.")

## ---- strict source-data sanity checks --------------------------------------
data_rep1_file <- scenario_data_file(reps_formal[1])
dataset_check <- check_s3_source_dataset(
    data_file = data_rep1_file,
    TT_use = TT_use,
    n1_use = n1_use,
    expected_x2_mode = x2_mode
)

## ---- run posterior performance analysis ------------------------------------
source(file.path(root_dir, s3_perf_file))

perf_out <- summarize_scenario3_posterior_performance(
    root = root_dir,
    fit_dir = file.path(output_dir, scenario_id),
    analysis_dir = file.path(analysis_dir, scenario_id),
    fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$",
    verbose = TRUE
)

s3_perf_rds_file <- file.path(
    analysis_dir,
    scenario_id,
    "scenario3_posterior_performance_results.rds"
)
saveRDS(perf_out, s3_perf_rds_file)
message("Saved Scenario 3 posterior performance RDS: ", s3_perf_rds_file)

## ---- add source-data confounding diagnostics --------------------------------
analysis_tables_dir <- file.path(root_dir, analysis_dir, scenario_id)
ensure_dir(analysis_tables_dir)

confounding_sum <- make_s3_confounding_table(
    reps = reps_formal,
    source_data_scenario_id = source_data_scenario_id,
    data_dir = data_dir,
    root = root_dir
)
utils::write.csv(
    confounding_sum,
    file.path(analysis_tables_dir, "scenario3_confounding_diagnostics.csv"),
    row.names = FALSE
)
utils::write.csv(
    source_manifest,
    file.path(analysis_tables_dir, "scenario3_source_data_manifest.csv"),
    row.names = FALSE
)
utils::write.csv(
    fit_manifest,
    file.path(analysis_tables_dir, "scenario3_fit_manifest.csv"),
    row.names = FALSE
)

## ---- run extra MCMC diagnostics --------------------------------------------
source(file.path(root_dir, s3_diag_file))

diag_out <- run_scenario3_extra_mcmc_diagnostics(
    root = root_dir,
    fit_dir = file.path(output_dir, scenario_id),
    analysis_dir = file.path(analysis_dir, scenario_id),
    fit_pattern = "fit_S3_dynamic_learned_gamma_rep[0-9]+\\.rds$",
    make_plots = TRUE,
    verbose = TRUE
)

s3_diag_rds_file <- file.path(
    analysis_dir,
    scenario_id,
    "scenario3_extra_mcmc_diagnostics_results.rds"
)
saveRDS(diag_out, s3_diag_rds_file)
message("Saved Scenario 3 extra MCMC diagnostics RDS: ", s3_diag_rds_file)

## ---- final summary ----------------------------------------------------------
print_s3_summary(
    dataset_check = dataset_check,
    perf_out = perf_out,
    diag_out = diag_out
)

cat("\n=== Main output locations ===\n")
cat("Shared DGP data: ", file.path(data_dir, source_data_scenario_id), "\n", sep = "")
cat("Fits:        ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("Analysis:    ", file.path(analysis_dir, scenario_id), "\n", sep = "")
cat("\nScenario 3 dynamic learned-gamma T100 validation finished successfully.\n")

invisible(list(
    source_manifest = source_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest,
    posterior_performance = perf_out,
    extra_mcmc_diagnostics = diag_out,
    confounding_sum = confounding_sum,
    dataset_check = dataset_check
))
