## ============================================================================
## run_s1_baseline_validation_T100_v2.R
##
## Formal Scenario 1 clean baseline validation.
##
## Purpose
##   Validate the baseline negative-binomial spatial regression component under
##   a clean data-generating mechanism with
##       lambda_tilde[t, j] = 1.
##   The fitted model fixes lambda at 1 and focuses on recovery of
##       beta0, beta1, beta2, phi, tau_phi, and r.
##
## Required files in project root
##   s1_baseline.R
##   run_scenario1_posterior_performance.R
##   run_scenario1_extra_mcmc_diagnostics.R   optional but recommended
##
## Main outputs
##   data_revised/S1_BASELINE_VALIDATION_T100/
##   output_s1_baseline_only/S1_BASELINE_VALIDATION_T100/
##   analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S1_BASELINE_VALIDATION_T100"
reps_formal <- 1:20

## Data-generating settings
TT_use <- 100L
n1_use <- 9L
beta0_truth <- -1.5
beta_truth <- c(0.5, -0.4)
r_truth <- 15
tau_phi_truth <- 2
residual_risk_value <- 1
omega_mode <- "fixed_prior_mean"

## Clean covariate design
## x1 is generated inside s1_baseline.R using its defaults.
## x2 is forced to be time-varying to avoid region-level x2 vs ICAR phi
## spatial confounding in the clean baseline validation.
x2_mode <- "continuous_time"
x2_ar <- 0.50
x2_innov_sd <- 0.80

## MCMC settings
n_iter <- 40000L
n_burnin <- 20000L
n_thin <- 5L

## Output settings
data_dir <- "data_revised"
output_dir <- "output_s1_baseline_only"
analysis_dir <- "analysis_s1_baseline_only"
plot_format <- "png"
verbose <- 1000L

## Set both to FALSE to reuse existing data and fits and only regenerate summaries.
overwrite_data <- TRUE
overwrite_fit <- TRUE

## Extra diagnostics
run_extra_mcmc_diagnostics <- TRUE
selected_phi_regions <- c(1L, 5L, 9L)
selected_r_regions <- c(1L, 5L, 9L)

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

read_required_csv <- function(path) {
    assert_file_exists(path, "CSV")
    utils::read.csv(path, stringsAsFactors = FALSE)
}

check_s1_dataset <- function(data_file, TT_use, n1_use) {
    assert_file_exists(data_file, "Scenario 1 dataset")
    dat <- readRDS(data_file)

    assert_true(!is.null(dat$y_coarse), "Dataset is missing y_coarse.")
    assert_true(!is.null(dat$x), "Dataset is missing x.")
    assert_true(!is.null(dat$lambda_tilde), "Dataset is missing lambda_tilde.")

    assert_true(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0(
            "y_coarse dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$y_coarse), collapse = " x "), "."
        )
    )

    assert_true(
        identical(as.integer(dim(dat$x)), c(as.integer(TT_use), as.integer(n1_use), 2L)),
        paste0(
            "x dimension is not TT_use by n1_use by 2. Got ",
            paste(dim(dat$x), collapse = " x "), "."
        )
    )

    assert_true(
        identical(as.integer(dim(dat$lambda_tilde)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0(
            "lambda_tilde dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$lambda_tilde), collapse = " x "), "."
        )
    )

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    assert_true(
        length(lambda_range) == 2L && all(abs(lambda_range - c(1, 1)) < 1e-12),
        paste0(
            "lambda_tilde is not fixed at 1. Range is ",
            paste(lambda_range, collapse = ", "), "."
        )
    )

    x1 <- dat$x[, , 1]
    x2 <- dat$x[, , 2]
    x2_unique_counts <- apply(x2, 2, function(z) length(unique(round(z, 8))))

    assert_true(
        all(x2_unique_counts >= max(10L, floor(TT_use / 2L))),
        paste0(
            "x2 does not appear to be time-varying in every region. Unique counts: ",
            paste(x2_unique_counts, collapse = ", "), "."
        )
    )

    x_cor <- stats::cor(as.vector(x1), as.vector(x2), use = "complete.obs")
    x1_sd <- stats::sd(as.vector(x1), na.rm = TRUE)
    x2_sd <- stats::sd(as.vector(x2), na.rm = TRUE)

    assert_true(is.finite(x1_sd) && abs(x1_sd - 1) < 1e-8, "x1 is not standardized to sd 1.")
    assert_true(is.finite(x2_sd) && abs(x2_sd - 1) < 1e-8, "x2 is not standardized to sd 1.")

    list(
        dat = dat,
        x2_unique_counts = x2_unique_counts,
        x_cor = x_cor,
        x1_sd = x1_sd,
        x2_sd = x2_sd,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0)
    )
}

check_performance_table <- function(perf, scenario_id, reps_formal, TT_use, n1_use) {
    required_cols <- c(
        "scenario_id", "rep_id", "TT", "n1", "mean_count", "zero_prop",
        "lambda_truth_min", "lambda_truth_max", "lambda_fixed_in_fit",
        "lambda_fixed_value"
    )
    missing_cols <- setdiff(required_cols, names(perf))
    assert_true(
        length(missing_cols) == 0L,
        paste("posterior_performance_diagnostics.csv is missing columns:", paste(missing_cols, collapse = ", "))
    )

    assert_true(nrow(perf) == length(reps_formal), "Number of completed replicates does not match reps_formal.")
    assert_true(all(sort(unique(perf$rep_id)) == sort(as.integer(reps_formal))), "rep_id values do not match reps_formal.")
    assert_true(all(perf$scenario_id == scenario_id), "scenario_id in diagnostics does not match requested scenario_id.")
    assert_true(all(as.integer(perf$TT) == as.integer(TT_use)), "TT in diagnostics does not match TT_use.")
    assert_true(all(as.integer(perf$n1) == as.integer(n1_use)), "n1 in diagnostics does not match n1_use.")

    lambda_ok <- all(abs(perf$lambda_truth_min - 1) < 1e-12) &&
        all(abs(perf$lambda_truth_max - 1) < 1e-12) &&
        all(as.logical(perf$lambda_fixed_in_fit)) &&
        all(abs(perf$lambda_fixed_value - 1) < 1e-12)

    assert_true(lambda_ok, "Lambda fixed-at-one check failed. Inspect posterior_performance_diagnostics.csv.")

    invisible(TRUE)
}

check_beta_summary <- function(beta_sum, reps_formal) {
    required_cols <- c("rep_id", "parameter")
    missing_cols <- setdiff(required_cols, names(beta_sum))
    assert_true(
        length(missing_cols) == 0L,
        paste("posterior_beta_summary.csv is missing columns:", paste(missing_cols, collapse = ", "))
    )

    beta_tab <- table(beta_sum$parameter)
    expected_names <- c("beta0", "beta1", "beta2")
    assert_true(
        all(expected_names %in% names(beta_tab)),
        "posterior_beta_summary.csv does not contain beta0, beta1, and beta2."
    )
    assert_true(
        all(beta_tab[expected_names] == length(reps_formal)),
        "posterior_beta_summary.csv does not contain one row per beta parameter per replicate."
    )

    invisible(TRUE)
}

print_s1_summary <- function(perf, dataset_check) {
    cat("\n=== Scenario 1 clean baseline validation checks ===\n")
    cat("Scenario: ", scenario_id, "\n", sep = "")
    cat("Replicates: ", min(reps_formal), " to ", max(reps_formal), " (n = ", length(reps_formal), ")\n", sep = "")
    cat("TT_use: ", TT_use, "\n", sep = "")
    cat("n1_use: ", n1_use, "\n", sep = "")
    cat("lambda_tilde: fixed at 1 in truth and fit\n")
    cat("x2 mode: ", x2_mode, ", AR = ", x2_ar, ", innovation SD = ", x2_innov_sd, "\n", sep = "")

    cat("\nData check from rep 1:\n")
    cat("  mean count: ", round(dataset_check$mean_count, 3), "\n", sep = "")
    cat("  zero proportion: ", round(dataset_check$zero_prop, 3), "\n", sep = "")
    cat("  x2 unique counts by region: ", paste(dataset_check$x2_unique_counts, collapse = ", "), "\n", sep = "")
    cat("  cor(vec(x1), vec(x2)): ", round(dataset_check$x_cor, 4), "\n", sep = "")
    cat("  sd(vec(x1)): ", round(dataset_check$x1_sd, 4), "\n", sep = "")
    cat("  sd(vec(x2)): ", round(dataset_check$x2_sd, 4), "\n", sep = "")

    cat("\nPerformance diagnostics across replicates:\n")
    cat("  completed replicates: ", nrow(perf), "\n", sep = "")
    cat("  mean count summary:\n")
    print(summary(perf$mean_count))
    cat("  zero proportion summary:\n")
    print(summary(perf$zero_prop))
}

## ---- required source files --------------------------------------------------
perf_script <- file.path(root_dir, "run_scenario1_posterior_performance.R")
core_script <- file.path(root_dir, "s1_baseline.R")
extra_diag_script <- file.path(root_dir, "run_scenario1_extra_mcmc_diagnostics.R")

assert_file_exists(perf_script, "posterior performance script")
assert_file_exists(core_script, "Scenario 1 core script")

## ---- clean old output if requested -----------------------------------------
if (isTRUE(overwrite_data)) {
    unlink(file.path(data_dir, scenario_id), recursive = TRUE)
}
if (isTRUE(overwrite_fit)) {
    unlink(file.path(output_dir, scenario_id), recursive = TRUE)
    unlink(file.path(analysis_dir, scenario_id), recursive = TRUE)
}

## ---- run Scenario 1 posterior performance analysis -------------------------
source(perf_script)

out_s1 <- run_scenario1_posterior_performance(
    root = root_dir,
    s1_core_file = "s1_baseline.R",

    scenario_id = scenario_id,
    reps = reps_formal,

    TT_use = TT_use,
    n1_use = n1_use,

    beta0_truth = beta0_truth,
    beta_truth = beta_truth,
    r_truth = r_truth,
    tau_phi_truth = tau_phi_truth,
    residual_risk_value = residual_risk_value,
    omega_mode = omega_mode,

    x2_mode = x2_mode,
    x2_ar = x2_ar,
    x2_innov_sd = x2_innov_sd,

    n_iter = n_iter,
    n_burnin = n_burnin,
    n_thin = n_thin,

    overwrite_data = overwrite_data,
    overwrite_fit = overwrite_fit,
    plot_format = plot_format,
    verbose = verbose
)

## ---- strict sanity checks ---------------------------------------------------
perf_file <- file.path(
    analysis_dir,
    scenario_id,
    "tables",
    "posterior_performance_diagnostics.csv"
)

beta_file <- file.path(
    analysis_dir,
    scenario_id,
    "tables",
    "posterior_beta_summary.csv"
)

data_rep1_file <- file.path(
    data_dir,
    scenario_id,
    sprintf("data_rep%02d.rds", reps_formal[1])
)

perf <- read_required_csv(perf_file)
beta_sum <- read_required_csv(beta_file)
dataset_check <- check_s1_dataset(data_rep1_file, TT_use = TT_use, n1_use = n1_use)

check_performance_table(
    perf = perf,
    scenario_id = scenario_id,
    reps_formal = reps_formal,
    TT_use = TT_use,
    n1_use = n1_use
)

check_beta_summary(beta_sum, reps_formal = reps_formal)
print_s1_summary(perf, dataset_check)

## ---- optional extra MCMC diagnostics ---------------------------------------
if (isTRUE(run_extra_mcmc_diagnostics)) {
    if (file.exists(extra_diag_script)) {
        source(extra_diag_script)

        diag_out <- run_scenario1_extra_mcmc_diagnostics(
            root = root_dir,
            reps = reps_formal,
            scenario_id = scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            analysis_dir = analysis_dir,
            selected_phi_regions = selected_phi_regions,
            selected_r_regions = selected_r_regions,
            plot_format = plot_format
        )

        mcmc_diag_file <- file.path(
            analysis_dir,
            scenario_id,
            "tables",
            "s1_mcmc_scalar_diagnostics.csv"
        )
        mcmc_flags_file <- file.path(
            analysis_dir,
            scenario_id,
            "tables",
            "s1_mcmc_scalar_diagnostic_flags.csv"
        )

        mcmc_diag <- read_required_csv(mcmc_diag_file)
        mcmc_flags <- read_required_csv(mcmc_flags_file)

        assert_true(
            length(unique(mcmc_diag$rep_id)) == length(reps_formal),
            "MCMC diagnostics do not contain all requested replicates."
        )
        assert_true(
            all(sort(unique(mcmc_diag$rep_id)) == sort(as.integer(reps_formal))),
            "MCMC diagnostics rep_id values do not match reps_formal."
        )

        cat("\nMCMC diagnostics check:\n")
        cat("  scalar diagnostic rows: ", nrow(mcmc_diag), "\n", sep = "")
        cat("  flag rows: ", nrow(mcmc_flags), "\n", sep = "")
        if ("low_ess" %in% names(mcmc_flags)) {
            cat("  low ESS flags: ", sum(as.logical(mcmc_flags$low_ess), na.rm = TRUE), "\n", sep = "")
        }
        if ("high_mcse_ratio" %in% names(mcmc_flags)) {
            cat("  high MCSE-ratio flags: ", sum(as.logical(mcmc_flags$high_mcse_ratio), na.rm = TRUE), "\n", sep = "")
        }
        if ("large_split_z" %in% names(mcmc_flags)) {
            cat("  large split-z flags: ", sum(as.logical(mcmc_flags$large_split_z), na.rm = TRUE), "\n", sep = "")
        }
    } else {
        diag_out <- NULL
        message("Optional MCMC diagnostics script not found. Skipping extra diagnostics: ", extra_diag_script)
    }
} else {
    diag_out <- NULL
}

## ---- final output summary ---------------------------------------------------
cat("\n=== Main output locations ===\n")
cat("Data:     ", file.path(data_dir, scenario_id), "\n", sep = "")
cat("Fits:     ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("Analysis: ", file.path(analysis_dir, scenario_id), "\n", sep = "")
cat("\nScenario 1 clean baseline validation finished successfully.\n")

invisible(list(
    out_s1 = out_s1,
    perf = perf,
    beta_sum = beta_sum,
    dataset_check = dataset_check,
    diag_out = diag_out
))
