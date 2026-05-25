## ============================================================================
## run_s2_dynamic_fixed_gamma_T100_v1.R
##
## Formal Scenario 2 dynamic residual-risk validation with fixed gamma.
##
## Purpose
##   Validate recovery of dynamic residual-risk paths under a controlled
##   data-generating mechanism with time-varying lambda_tilde[t, j].
##   The fitted model samples lambda_tilde but fixes gamma at the true value.
##   This separates latent path recovery from discount-factor learning.
##
## Required files in project root
##   s2_dynamic_fixed_gamma_v2_with_regression_dispersion_patch1.R
##
## Main outputs
##   data_revised/S2_DYNAMIC_FIXED_GAMMA_T100/
##   output_s2_dynamic_fixed_gamma/S2_DYNAMIC_FIXED_GAMMA_T100/
##   analysis_s2_dynamic_fixed_gamma/S2_DYNAMIC_FIXED_GAMMA_T100/
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S2_DYNAMIC_FIXED_GAMMA_T100"
reps_formal <- 1:20

## Core Scenario 2 script
s2_core_file <- "s2_dynamic_fixed_gamma_v2_with_regression_dispersion_patch1.R"

## Data-generating settings
TT_use <- 100L
n1_use <- 9L
beta0_truth <- -1.5
beta_truth <- c(0.5, -0.4)
r_truth <- 15
tau_phi_truth <- 2
omega_mode <- "fixed_prior_mean"

## Dynamic residual-risk settings
gamma_truth <- 0.80
A0_use <- 10
B0_use <- 10

## Clean covariate design
## This must remain continuous_time to match Scenario 1's clean covariate design
## and avoid region-level x2 versus ICAR phi spatial confounding.
x2_mode <- "continuous_time"
x2_ar <- 0.50
x2_innov_sd <- 0.80

## MCMC settings
n_iter <- 40000L
n_burnin <- 20000L
n_thin <- 5L

## Output settings
data_dir <- "data_revised"
output_dir <- "output_s2_dynamic_fixed_gamma"
analysis_dir <- "analysis_s2_dynamic_fixed_gamma"
plot_format <- "png"
verbose <- 1000L

## Set both to FALSE to reuse existing data and fits and only regenerate summaries.
overwrite_data <- TRUE
overwrite_fit <- TRUE

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

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

check_s2_dataset <- function(data_file, TT_use, n1_use,
                             expected_x2_mode = "continuous_time") {
    assert_file_exists(data_file, "Scenario 2 dataset")
    dat <- readRDS(data_file)

    assert_true(!is.null(dat$y_coarse), "Dataset is missing y_coarse.")
    assert_true(!is.null(dat$x1), "Dataset is missing x1.")
    assert_true(!is.null(dat$x2), "Dataset is missing x2.")
    assert_true(!is.null(dat$lambda_tilde), "Dataset is missing lambda_tilde.")
    assert_true(!is.null(dat$lambda_tilde_ident), "Dataset is missing lambda_tilde_ident.")
    assert_true(!is.null(dat$gamma_star), "Dataset is missing gamma_star.")

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
        length(lambda_range) == 2L && is.finite(lambda_range[1]) && is.finite(lambda_range[2]) &&
            lambda_range[1] > 0 && diff(lambda_range) > 1e-8,
        paste0("lambda_tilde_ident does not appear to be positive and dynamic. Range is ",
               paste(lambda_range, collapse = ", "), ".")
    )

    gamma_range <- range(dat$gamma_star, finite = TRUE)
    assert_true(
        all(abs(gamma_range - gamma_truth) < 1e-12),
        paste0("gamma_star is not fixed at gamma_truth. Range is ", paste(gamma_range, collapse = ", "), ".")
    )

    x1 <- as.vector(dat$x1)
    x2 <- as.vector(dat$x2)
    loglam <- as.vector(log(dat$lambda_tilde_ident))
    phi_mat <- matrix(dat$phi_star, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)

    x_cor <- stats::cor(x1, x2, use = "complete.obs")
    x1_sd <- stats::sd(x1, na.rm = TRUE)
    x2_sd <- stats::sd(x2, na.rm = TRUE)

    assert_true(is.finite(x1_sd) && abs(x1_sd - 1) < 1e-8, "x1 is not standardized to sd 1.")
    assert_true(is.finite(x2_sd) && abs(x2_sd - 1) < 1e-8, "x2 is not standardized to sd 1.")

    list(
        dat = dat,
        x2_unique_counts = x2_unique_counts,
        x_cor = x_cor,
        x1_sd = x1_sd,
        x2_sd = x2_sd,
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

check_s2_summary <- function(summary_all, scenario_id, reps_formal, TT_use, n1_use) {
    required_cols <- c(
        "scenario_id", "rep_id", "TT", "n1", "mean_count", "zero_prop",
        "dynamic_lambda_in_truth", "lambda_sampled_in_fit", "gamma_fixed_in_fit",
        "gamma_truth_mean", "gamma_fixed_mean",
        "lambda_truth_min", "lambda_truth_max",
        "beta0_mean", "beta1_mean", "beta2_mean",
        "lambda_rmse", "log_lambda_rmse", "log_lambda_coverage_95",
        "cor_log_lambda", "cor_delta_log_lambda"
    )
    missing_cols <- setdiff(required_cols, names(summary_all))
    assert_true(
        length(missing_cols) == 0L,
        paste("S2 combined summary is missing columns:", paste(missing_cols, collapse = ", "))
    )

    assert_true(nrow(summary_all) == length(reps_formal), "Number of completed replicates does not match reps_formal.")
    assert_true(all(sort(unique(summary_all$rep_id)) == sort(as.integer(reps_formal))), "rep_id values do not match reps_formal.")
    assert_true(all(summary_all$scenario_id == scenario_id), "scenario_id in summary does not match requested scenario_id.")
    assert_true(all(as.integer(summary_all$TT) == as.integer(TT_use)), "TT in summary does not match TT_use.")
    assert_true(all(as.integer(summary_all$n1) == as.integer(n1_use)), "n1 in summary does not match n1_use.")
    assert_true(all(as.logical(summary_all$dynamic_lambda_in_truth)), "dynamic_lambda_in_truth is not TRUE for all rows.")
    assert_true(all(as.logical(summary_all$lambda_sampled_in_fit)), "lambda_sampled_in_fit is not TRUE for all rows.")
    assert_true(all(as.logical(summary_all$gamma_fixed_in_fit)), "gamma_fixed_in_fit is not TRUE for all rows.")
    assert_true(all(abs(summary_all$gamma_truth_mean - gamma_truth) < 1e-12), "gamma_truth_mean does not equal gamma_truth.")
    assert_true(all(abs(summary_all$gamma_fixed_mean - gamma_truth) < 1e-12), "gamma_fixed_mean does not equal gamma_truth.")
    assert_true(all(is.finite(summary_all$log_lambda_rmse)), "log_lambda_rmse contains non-finite values.")
    assert_true(all(is.finite(summary_all$cor_log_lambda)), "cor_log_lambda contains non-finite values.")

    invisible(TRUE)
}

make_s2_beta_summary_table <- function(summary_all) {
    params <- list(
        beta0 = c(true = "beta0_true", mean = "beta0_mean", sd = "beta0_sd", q025 = "beta0_q025",
                  q50 = "beta0_q50", q975 = "beta0_q975", bias = "beta0_bias", covered = "beta0_covered"),
        beta1 = c(true = "beta1_true", mean = "beta1_mean", sd = "beta1_sd", q025 = "beta1_q025",
                  q50 = "beta1_q50", q975 = "beta1_q975", bias = "beta1_bias", covered = "beta1_covered"),
        beta2 = c(true = "beta2_true", mean = "beta2_mean", sd = "beta2_sd", q025 = "beta2_q025",
                  q50 = "beta2_q50", q975 = "beta2_q975", bias = "beta2_bias", covered = "beta2_covered")
    )

    out <- lapply(names(params), function(par) {
        cols <- params[[par]]
        data.frame(
            scenario_id = summary_all$scenario_id,
            rep_id = summary_all$rep_id,
            parameter = par,
            truth = summary_all[[cols["true"]]],
            post_mean = summary_all[[cols["mean"]]],
            post_sd = summary_all[[cols["sd"]]],
            q025 = summary_all[[cols["q025"]]],
            q50 = summary_all[[cols["q50"]]],
            q975 = summary_all[[cols["q975"]]],
            bias = summary_all[[cols["bias"]]],
            covered = summary_all[[cols["covered"]]],
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, out)
}

make_s2_lambda_summary_table <- function(summary_all) {
    keep <- c(
        "scenario_id", "rep_id", "TT", "n1", "gamma_truth_mean", "gamma_fixed_mean",
        "lambda_truth_min", "lambda_truth_max", "lambda_post_mean_min", "lambda_post_mean_max",
        "lambda_rmse", "lambda_mae", "lambda_coverage_95",
        "log_lambda_rmse", "log_lambda_mae", "log_lambda_coverage_95",
        "cor_lambda", "cor_log_lambda",
        "rmse_delta_log_lambda", "mae_delta_log_lambda", "cor_delta_log_lambda"
    )
    keep <- intersect(keep, names(summary_all))
    summary_all[, keep, drop = FALSE]
}

make_s2_confounding_table <- function(reps, scenario_id, data_dir, root = ".") {
    out <- lapply(reps, function(rep_id) {
        rr <- sprintf("%02d", as.integer(rep_id))
        dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
        dat <- readRDS(dat_file)

        x1 <- as.vector(dat$x1)
        x2 <- as.vector(dat$x2)
        loglam <- as.vector(log(dat$lambda_tilde_ident %||% dat$lambda_tilde))
        phi_mat <- matrix(dat$phi_star, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)

        data.frame(
            scenario_id = scenario_id,
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

print_s2_summary <- function(summary_all, dataset_check, beta_sum, lambda_sum) {
    cat("\n=== Scenario 2 dynamic fixed-gamma validation checks ===\n")
    cat("Scenario: ", scenario_id, "\n", sep = "")
    cat("Replicates: ", min(reps_formal), " to ", max(reps_formal), " (n = ", length(reps_formal), ")\n", sep = "")
    cat("TT_use: ", TT_use, "\n", sep = "")
    cat("n1_use: ", n1_use, "\n", sep = "")
    cat("lambda_tilde: dynamic in truth and sampled in fit\n")
    cat("gamma: fixed at truth, gamma = ", gamma_truth, "\n", sep = "")
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

    cat("\nRegression recovery across replicates:\n")
    beta_agg <- aggregate(
        cbind(post_mean, bias, covered) ~ parameter,
        data = beta_sum,
        FUN = function(z) c(mean = mean(z), sd = stats::sd(z))
    )
    print(beta_agg)

    cat("\nLambda path recovery across replicates:\n")
    print(summary(lambda_sum[, c("log_lambda_rmse", "log_lambda_coverage_95", "cor_log_lambda", "cor_delta_log_lambda")]))

    cat("\nMCMC diagnostic summaries:\n")
    diag_cols <- intersect(c("phi_accept_rate", "r_accept_rate_mean", "beta_mean_n_reject"), names(summary_all))
    print(summary(summary_all[, diag_cols, drop = FALSE]))
}

## ---- required source files --------------------------------------------------
core_script <- file.path(root_dir, s2_core_file)
assert_file_exists(core_script, "Scenario 2 core script")

## ---- clean old output if requested -----------------------------------------
if (isTRUE(overwrite_data)) {
    unlink(file.path(data_dir, scenario_id), recursive = TRUE)
}
if (isTRUE(overwrite_fit)) {
    unlink(file.path(output_dir, scenario_id), recursive = TRUE)
    unlink(file.path(analysis_dir, scenario_id), recursive = TRUE)
}

## ---- run Scenario 2 ---------------------------------------------------------
source(core_script)
source_s2_dynamic_fixed_gamma(root = root_dir)

## Generate data explicitly so x2_mode and T100 are guaranteed.
s2_manifest <- simulate_s2_dynamic_fixed_gamma_batch(
    reps = reps_formal,
    data_dir = data_dir,
    scenario_id = scenario_id,
    root = root_dir,
    overwrite_existing = overwrite_data,
    verbose = TRUE,

    TT_use = TT_use,
    n1_use = n1_use,
    beta0_truth = beta0_truth,
    beta_truth = beta_truth,
    r_truth = r_truth,
    tau_phi_truth = tau_phi_truth,
    gamma_truth = gamma_truth,
    A0_use = A0_use,
    B0_use = B0_use,
    omega_mode = omega_mode,
    x2_mode = x2_mode,
    x2_ar = x2_ar,
    x2_innov_sd = x2_innov_sd
)

settings_override <- list(
    n_iter = n_iter,
    n_burnin = n_burnin,
    n_thin = n_thin
)

summary_all <- fit_s2_dynamic_fixed_gamma_batch(
    reps = reps_formal,
    scenario_id = scenario_id,
    data_dir = data_dir,
    output_dir = output_dir,
    root = root_dir,
    settings_override = settings_override,
    verbose = verbose,
    overwrite_existing = overwrite_fit
)

## ---- strict sanity checks ---------------------------------------------------
data_rep1_file <- file.path(
    root_dir,
    data_dir,
    scenario_id,
    sprintf("data_rep%02d.rds", reps_formal[1])
)

dataset_check <- check_s2_dataset(
    data_file = data_rep1_file,
    TT_use = TT_use,
    n1_use = n1_use,
    expected_x2_mode = x2_mode
)

check_s2_summary(
    summary_all = summary_all,
    scenario_id = scenario_id,
    reps_formal = reps_formal,
    TT_use = TT_use,
    n1_use = n1_use
)

## ---- organize analysis tables ---------------------------------------------
analysis_tables_dir <- file.path(root_dir, analysis_dir, scenario_id, "tables")
ensure_dir(analysis_tables_dir)

posterior_perf <- summary_all
beta_sum <- make_s2_beta_summary_table(summary_all)
lambda_sum <- make_s2_lambda_summary_table(summary_all)
confounding_sum <- make_s2_confounding_table(
    reps = reps_formal,
    scenario_id = scenario_id,
    data_dir = data_dir,
    root = root_dir
)

utils::write.csv(
    posterior_perf,
    file.path(analysis_tables_dir, "posterior_performance_diagnostics.csv"),
    row.names = FALSE
)
utils::write.csv(
    beta_sum,
    file.path(analysis_tables_dir, "posterior_beta_summary.csv"),
    row.names = FALSE
)
utils::write.csv(
    lambda_sum,
    file.path(analysis_tables_dir, "posterior_lambda_path_recovery.csv"),
    row.names = FALSE
)
utils::write.csv(
    confounding_sum,
    file.path(analysis_tables_dir, "scenario2_confounding_diagnostics.csv"),
    row.names = FALSE
)
utils::write.csv(
    s2_manifest,
    file.path(analysis_tables_dir, "scenario2_data_manifest.csv"),
    row.names = FALSE
)

## ---- final checks and summary ---------------------------------------------
check_s2_summary(
    summary_all = posterior_perf,
    scenario_id = scenario_id,
    reps_formal = reps_formal,
    TT_use = TT_use,
    n1_use = n1_use
)

assert_true(
    all(c("beta0", "beta1", "beta2") %in% unique(beta_sum$parameter)),
    "posterior_beta_summary.csv does not contain beta0, beta1, and beta2."
)

print_s2_summary(
    summary_all = posterior_perf,
    dataset_check = dataset_check,
    beta_sum = beta_sum,
    lambda_sum = lambda_sum
)

cat("\n=== Main output locations ===\n")
cat("Data:     ", file.path(data_dir, scenario_id), "\n", sep = "")
cat("Fits:     ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("Analysis: ", file.path(analysis_dir, scenario_id), "\n", sep = "")
cat("\nScenario 2 dynamic fixed-gamma validation finished successfully.\n")

invisible(list(
    manifest = s2_manifest,
    posterior_perf = posterior_perf,
    beta_sum = beta_sum,
    lambda_sum = lambda_sum,
    confounding_sum = confounding_sum,
    dataset_check = dataset_check
))
