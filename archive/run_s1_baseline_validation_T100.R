## ============================================================================
## run_s1_baseline_validation_T100.R
##
## Formal Scenario 1 baseline validation.
##
## Purpose:
##   Validate the baseline NB spatial regression component under a clean DGP:
##       lambda_tilde[t, j] = 1
##   so that recovery focuses on beta, phi, tau_phi, and r.
##
## Required files in project root:
##   s1_baseline.R
##   run_scenario1_posterior_performance.R
##   run_scenario1_extra_mcmc_diagnostics.R   optional, for extra diagnostics
## ============================================================================

# unlink("data_revised/S1_BASELINE_VALIDATION_T100", recursive = TRUE)
# unlink("output_s1_baseline_only/S1_BASELINE_VALIDATION_T100", recursive = TRUE)
# unlink("analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100", recursive = TRUE)
#
# dir.exists("data_revised/S1_BASELINE_VALIDATION_T100")
# dir.exists("output_s1_baseline_only/S1_BASELINE_VALIDATION_T100")
# dir.exists("analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100")

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S1_BASELINE_VALIDATION_T100"
reps_formal <- 1:20

## Data generating settings
TT_use <- 100L
n1_use <- 9L
beta0_truth <- -1.5
beta_truth <- c(0.5, -0.4)
r_truth <- 15
tau_phi_truth <- 2
residual_risk_value <- 1
omega_mode <- "fixed_prior_mean"

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

## Set these to FALSE if you want to resume without overwriting previous work.
overwrite_data <- TRUE
overwrite_fit <- TRUE

## ---- clean old output if overwriting ---------------------------------------
if (isTRUE(overwrite_data)) {
    unlink(file.path(data_dir, scenario_id), recursive = TRUE)
}
if (isTRUE(overwrite_fit)) {
    unlink(file.path(output_dir, scenario_id), recursive = TRUE)
    unlink(file.path(analysis_dir, scenario_id), recursive = TRUE)
}

## ---- run posterior performance analysis ------------------------------------
source(file.path(root_dir, "run_scenario1_posterior_performance.R"))

out_s1 <- run_scenario1_posterior_performance(
    root = ".",
    s1_core_file = "s1_baseline.R",

    scenario_id = "S1_BASELINE_VALIDATION_T100",
    reps = 1:20,

    TT_use = 100,
    n1_use = 9,

    beta0_truth = -1.5,
    beta_truth = c(0.5, -0.4),
    r_truth = 15,
    tau_phi_truth = 2,
    residual_risk_value = 1,
    omega_mode = "fixed_prior_mean",

    x2_mode = "continuous_time",
    x2_ar = 0.50,
    x2_innov_sd = 0.80,

    n_iter = 40000,
    n_burnin = 20000,
    n_thin = 5,

    overwrite_data = TRUE,
    overwrite_fit = TRUE,
    plot_format = "png",
    verbose = 500
)

## ---- sanity checks ----------------------------------------------------------
perf_file <- file.path(
    analysis_dir,
    scenario_id,
    "tables",
    "posterior_performance_diagnostics.csv"
)

if (!file.exists(perf_file)) {
    stop("Expected performance diagnostics CSV was not created: ", perf_file)
}

perf <- read.csv(perf_file)

cat("\n=== Scenario 1 formal T100 validation checks ===\n")
cat("Scenario:", scenario_id, "\n")
cat("Number of completed replicates:", nrow(perf), "\n")
cat("Unique TT values:", paste(unique(perf$TT), collapse = ", "), "\n")
cat("Unique n1 values:", paste(unique(perf$n1), collapse = ", "), "\n")
cat("Mean count summary:\n")
print(summary(perf$mean_count))
cat("Zero proportion summary:\n")
print(summary(perf$zero_prop))

lambda_ok <- all(perf$lambda_truth_min == 1) &&
    all(perf$lambda_truth_max == 1) &&
    all(perf$lambda_fixed_in_fit) &&
    all(perf$lambda_fixed_value == 1)

cat("Lambda fixed at 1 in truth and fit:", lambda_ok, "\n")

if (!lambda_ok) {
    warning("Lambda fixed-at-one check failed. Inspect posterior_performance_diagnostics.csv.")
}

## ---- optional extra MCMC diagnostics ---------------------------------------
if (file.exists(file.path(root_dir, "run_scenario1_extra_mcmc_diagnostics.R"))) {
    source(file.path(root_dir, "run_scenario1_extra_mcmc_diagnostics.R"))

    diag_out <- run_scenario1_extra_mcmc_diagnostics(
        # root = ".",
        # reps = 1:20,
        # scenario_id = "S1_BASELINE_VALIDATION_T100",
        # data_dir = "data_revised",
        # output_dir = "output_s1_baseline_only",
        # analysis_dir = "analysis_s1_baseline_only",
        # plot_format = "png"
        #
        # selected_phi_regions = c(1, 5, 9),
        # selected_r_regions = c(1, 5, 9)
        root = root_dir,
        reps = reps_formal,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        analysis_dir = analysis_dir,
        plot_format = plot_format

    )
} else {
    diag_out <- NULL
    message("Optional MCMC diagnostics script not found. Skipping extra diagnostics.")
}

cat("\n=== Main output locations ===\n")
cat("Data:     ", file.path(data_dir, scenario_id), "\n")
cat("Fits:     ", file.path(output_dir, scenario_id), "\n")
cat("Analysis: ", file.path(analysis_dir, scenario_id), "\n")

invisible(list(out_s1 = out_s1, perf = perf, diag_out = diag_out))




# diag <- read.csv(
#     "analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100/tables/posterior_performance_diagnostics.csv"
# )
#
# unique(diag$scenario_id)
# unique(diag$TT)
# unique(diag$n1)
# nrow(diag)
#
# dat <- readRDS("data_revised/S1_BASELINE_VALIDATION_T100/data_rep01.rds")
# x2 <- dat$x[, , 2]
#
# apply(x2, 2, function(z) length(unique(round(z, 8))))
#
#
#
# beta_sum <- read.csv(
#     "analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100/tables/posterior_beta_summary.csv"
# )
#
# diag <- read.csv(
#     "analysis_s1_baseline_only/S1_BASELINE_VALIDATION_T100/tables/posterior_performance_diagnostics.csv"
# )
#
# unique(diag$TT)
# unique(diag$n1)
# nrow(diag)
# table(beta_sum$parameter)
# nrow(beta_sum)