# Scenario 4C posterior performance summary
# Integrates sampled-lambda formal fit and S4C oracle diagnostics.
# Source from project root after S4C sampled and oracle diagnostics have been run.
# This version intentionally avoids row-binding helpers.

root <- getwd()

sampled_path <- file.path(
  root,
  "output_s4c_small_r_fixed_gamma",
  "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100",
  "summary_S4C_small_r_fixed_gamma_all_reps.csv"
)

oracle_lambda_path <- file.path(
  root,
  "output_s4c_oracle_lambda_fixed_gamma",
  "S4C_ORACLE_LAMBDA_FIXED_GAMMA_DIAGNOSTIC_T100",
  "summary_S4C_oracle_lambda_fixed_gamma_all_reps.csv"
)

oracle_lam_phi_path <- file.path(
  root,
  "output_s4c_oracle_lambda_phi_fixed_gamma",
  "S4C_ORACLE_LAMBDA_PHI_FIXED_GAMMA_DIAGNOSTIC_T100",
  "summary_S4C_oracle_lambda_phi_fixed_gamma_all_reps.csv"
)

out_dir <- file.path(
  root,
  "analysis_s4c_overdispersion",
  "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100",
  "tables"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  message("Found ", label, ": ", normalizePath(path, winslash = "/", mustWork = FALSE))
}

stop_if_missing(sampled_path, "sampled summary")
stop_if_missing(oracle_lambda_path, "oracle-lambda summary")
stop_if_missing(oracle_lam_phi_path, "oracle-lambda+phi summary")
message("Writing repaired S4C tables to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

sampled <- read.csv(sampled_path, stringsAsFactors = FALSE)
oracle_lambda <- read.csv(oracle_lambda_path, stringsAsFactors = FALSE)
oracle_lam_phi <- read.csv(oracle_lam_phi_path, stringsAsFactors = FALSE)

get_num <- function(df, col) {
  if (col %in% names(df)) {
    return(suppressWarnings(as.numeric(df[[col]])))
  }
  rep(NA_real_, nrow(df))
}

get_chr <- function(df, col) {
  if (col %in% names(df)) {
    return(as.character(df[[col]]))
  }
  rep(NA_character_, nrow(df))
}

mean_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

sd_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) <= 1 || all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

sum_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(0)
  sum(x, na.rm = TRUE)
}

first_safe <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

reps_label <- function(df) {
  if (!("rep_id" %in% names(df)) || nrow(df) == 0) return(NA_character_)
  paste(df$rep_id, collapse = ";")
}

# Classify sampled-lambda fits for repair purposes.
beta_guard <- get_num(sampled, "s4c_beta_guard_count")
kappa_guard <- get_num(sampled, "s4c_kappa_guard_count")
lambda_input_guard <- get_num(sampled, "s4c_lambda_input_guard_count")
lambda_output_guard <- get_num(sampled, "s4c_lambda_output_guard_count")
log_lambda_rmse <- get_num(sampled, "log_lambda_rmse")
lambda_rmse <- get_num(sampled, "lambda_rmse")
beta0_mean <- get_num(sampled, "beta0_mean")

sampled_failure <-
  (!is.finite(log_lambda_rmse)) |
  (!is.finite(lambda_rmse)) |
  log_lambda_rmse > 5 |
  lambda_rmse > 1e6 |
  abs(beta0_mean) > 20 |
  beta_guard > 100 |
  kappa_guard > 100 |
  lambda_input_guard > 100 |
  lambda_output_guard > 100

sampled$repair_fit_status <- ifelse(sampled_failure, "numerical_instability", "stable")
stable_df <- sampled[sampled$repair_fit_status == "stable", , drop = FALSE]
failure_df <- sampled[sampled$repair_fit_status == "numerical_instability", , drop = FALSE]

summarize_sampled <- function(df, label) {
  data.frame(
    diagnostic_label = label,
    n = nrow(df),
    reps = reps_label(df),
    mean_count_avg = mean_safe(get_num(df, "mean_count")),
    zero_prop_avg = mean_safe(get_num(df, "zero_prop")),
    count_cv_avg = mean_safe(get_num(df, "count_cv")),
    beta0_mean_avg = mean_safe(get_num(df, "beta0_mean")),
    beta1_mean_avg = mean_safe(get_num(df, "beta1_mean")),
    beta2_mean_avg = mean_safe(get_num(df, "beta2_mean")),
    beta2_mean_sd = sd_safe(get_num(df, "beta2_mean")),
    r_mean_avg = mean_safe(get_num(df, "r_mean")),
    r_region_coverage_95_avg = mean_safe(get_num(df, "r_region_coverage_95")),
    r_region_rmse_avg = mean_safe(get_num(df, "r_region_rmse")),
    r_region_mae_avg = mean_safe(get_num(df, "r_region_mae")),
    phi_rmse_avg = mean_safe(get_num(df, "phi_rmse")),
    phi_cor_avg = mean_safe(get_num(df, "phi_cor")),
    lambda_rmse_avg = mean_safe(get_num(df, "lambda_rmse")),
    log_lambda_rmse_avg = mean_safe(get_num(df, "log_lambda_rmse")),
    cor_log_lambda_avg = mean_safe(get_num(df, "cor_log_lambda")),
    kappa_truth_cv_avg = mean_safe(get_num(df, "kappa_truth_cv")),
    kappa_post_mean_cv_avg = mean_safe(get_num(df, "kappa_post_mean_cv")),
    kappa_rmse_avg = mean_safe(get_num(df, "kappa_rmse")),
    log_kappa_rmse_avg = mean_safe(get_num(df, "log_kappa_rmse")),
    cor_log_kappa_avg = mean_safe(get_num(df, "cor_log_kappa")),
    beta_guard_total = sum_safe(get_num(df, "s4c_beta_guard_count")),
    kappa_guard_total = sum_safe(get_num(df, "s4c_kappa_guard_count")),
    lambda_guard_total = sum_safe(get_num(df, "s4c_lambda_input_guard_count")) +
      sum_safe(get_num(df, "s4c_lambda_output_guard_count")),
    stringsAsFactors = FALSE
  )
}

summarize_oracle <- function(df, label) {
  data.frame(
    diagnostic_label = label,
    n = nrow(df),
    reps = reps_label(df),
    mean_count_avg = mean_safe(get_num(df, "mean_count")),
    zero_prop_avg = mean_safe(get_num(df, "zero_prop")),
    count_cv_avg = mean_safe(get_num(df, "count_cv")),
    beta0_mean_avg = mean_safe(get_num(df, "beta0_mean")),
    beta1_mean_avg = mean_safe(get_num(df, "beta1_mean")),
    beta2_mean_avg = mean_safe(get_num(df, "beta2_mean")),
    beta2_mean_sd = sd_safe(get_num(df, "beta2_mean")),
    r_mean_avg = mean_safe(get_num(df, "r_mean")),
    r_region_coverage_95_avg = mean_safe(get_num(df, "r_region_coverage_95")),
    r_region_rmse_avg = mean_safe(get_num(df, "r_region_rmse")),
    r_region_mae_avg = mean_safe(get_num(df, "r_region_mae")),
    phi_rmse_avg = mean_safe(get_num(df, "phi_rmse")),
    phi_cor_avg = mean_safe(get_num(df, "phi_cor")),
    lambda_rmse_avg = mean_safe(get_num(df, "lambda_rmse")),
    log_lambda_rmse_avg = mean_safe(get_num(df, "log_lambda_rmse")),
    cor_log_lambda_avg = mean_safe(get_num(df, "cor_log_lambda")),
    kappa_truth_cv_avg = mean_safe(get_num(df, "kappa_truth_cv")),
    kappa_post_mean_cv_avg = mean_safe(get_num(df, "kappa_post_mean_cv")),
    kappa_rmse_avg = mean_safe(get_num(df, "kappa_rmse")),
    log_kappa_rmse_avg = mean_safe(get_num(df, "log_kappa_rmse")),
    cor_log_kappa_avg = mean_safe(get_num(df, "cor_log_kappa")),
    beta_guard_total = sum_safe(get_num(df, "s4c_oracle_beta_guard_count")),
    kappa_guard_total = sum_safe(get_num(df, "s4c_oracle_kappa_guard_count")),
    lambda_guard_total = 0,
    stringsAsFactors = FALSE
  )
}

# Build detail table with a fixed schema.
oracle_detail <- data.frame(
  diagnostic_label = c("oracle_lambda", "oracle_lambda_phi"),
  rep_id = c(first_safe(get_num(oracle_lambda, "rep_id")), first_safe(get_num(oracle_lam_phi, "rep_id"))),
  mean_count = c(first_safe(get_num(oracle_lambda, "mean_count")), first_safe(get_num(oracle_lam_phi, "mean_count"))),
  zero_prop = c(first_safe(get_num(oracle_lambda, "zero_prop")), first_safe(get_num(oracle_lam_phi, "zero_prop"))),
  count_cv = c(first_safe(get_num(oracle_lambda, "count_cv")), first_safe(get_num(oracle_lam_phi, "count_cv"))),
  beta0_mean = c(first_safe(get_num(oracle_lambda, "beta0_mean")), first_safe(get_num(oracle_lam_phi, "beta0_mean"))),
  beta1_mean = c(first_safe(get_num(oracle_lambda, "beta1_mean")), first_safe(get_num(oracle_lam_phi, "beta1_mean"))),
  beta2_mean = c(first_safe(get_num(oracle_lambda, "beta2_mean")), first_safe(get_num(oracle_lam_phi, "beta2_mean"))),
  r_mean = c(first_safe(get_num(oracle_lambda, "r_mean")), first_safe(get_num(oracle_lam_phi, "r_mean"))),
  r_region_coverage_95 = c(first_safe(get_num(oracle_lambda, "r_region_coverage_95")), first_safe(get_num(oracle_lam_phi, "r_region_coverage_95"))),
  r_region_rmse = c(first_safe(get_num(oracle_lambda, "r_region_rmse")), first_safe(get_num(oracle_lam_phi, "r_region_rmse"))),
  r_region_mae = c(first_safe(get_num(oracle_lambda, "r_region_mae")), first_safe(get_num(oracle_lam_phi, "r_region_mae"))),
  phi_rmse = c(first_safe(get_num(oracle_lambda, "phi_rmse")), first_safe(get_num(oracle_lam_phi, "phi_rmse"))),
  phi_cor = c(first_safe(get_num(oracle_lambda, "phi_cor")), first_safe(get_num(oracle_lam_phi, "phi_cor"))),
  lambda_rmse = c(first_safe(get_num(oracle_lambda, "lambda_rmse")), first_safe(get_num(oracle_lam_phi, "lambda_rmse"))),
  log_lambda_rmse = c(first_safe(get_num(oracle_lambda, "log_lambda_rmse")), first_safe(get_num(oracle_lam_phi, "log_lambda_rmse"))),
  cor_log_lambda = c(first_safe(get_num(oracle_lambda, "cor_log_lambda")), first_safe(get_num(oracle_lam_phi, "cor_log_lambda"))),
  kappa_truth_cv = c(first_safe(get_num(oracle_lambda, "kappa_truth_cv")), first_safe(get_num(oracle_lam_phi, "kappa_truth_cv"))),
  kappa_post_mean_cv = c(first_safe(get_num(oracle_lambda, "kappa_post_mean_cv")), first_safe(get_num(oracle_lam_phi, "kappa_post_mean_cv"))),
  kappa_rmse = c(first_safe(get_num(oracle_lambda, "kappa_rmse")), first_safe(get_num(oracle_lam_phi, "kappa_rmse"))),
  log_kappa_rmse = c(first_safe(get_num(oracle_lambda, "log_kappa_rmse")), first_safe(get_num(oracle_lam_phi, "log_kappa_rmse"))),
  cor_log_kappa = c(first_safe(get_num(oracle_lambda, "cor_log_kappa")), first_safe(get_num(oracle_lam_phi, "cor_log_kappa"))),
  beta_guard = c(sum_safe(get_num(oracle_lambda, "s4c_oracle_beta_guard_count")), sum_safe(get_num(oracle_lam_phi, "s4c_oracle_beta_guard_count"))),
  kappa_guard = c(sum_safe(get_num(oracle_lambda, "s4c_oracle_kappa_guard_count")), sum_safe(get_num(oracle_lam_phi, "s4c_oracle_kappa_guard_count"))),
  stringsAsFactors = FALSE
)

# Build summary table with a fixed schema.
oracle_summary <- data.frame(
  diagnostic_label = c("oracle_lambda", "oracle_lambda_phi"),
  n = c(nrow(oracle_lambda), nrow(oracle_lam_phi)),
  reps = c(reps_label(oracle_lambda), reps_label(oracle_lam_phi)),
  mean_count_avg = c(mean_safe(get_num(oracle_lambda, "mean_count")), mean_safe(get_num(oracle_lam_phi, "mean_count"))),
  zero_prop_avg = c(mean_safe(get_num(oracle_lambda, "zero_prop")), mean_safe(get_num(oracle_lam_phi, "zero_prop"))),
  count_cv_avg = c(mean_safe(get_num(oracle_lambda, "count_cv")), mean_safe(get_num(oracle_lam_phi, "count_cv"))),
  beta0_mean_avg = c(mean_safe(get_num(oracle_lambda, "beta0_mean")), mean_safe(get_num(oracle_lam_phi, "beta0_mean"))),
  beta1_mean_avg = c(mean_safe(get_num(oracle_lambda, "beta1_mean")), mean_safe(get_num(oracle_lam_phi, "beta1_mean"))),
  beta2_mean_avg = c(mean_safe(get_num(oracle_lambda, "beta2_mean")), mean_safe(get_num(oracle_lam_phi, "beta2_mean"))),
  beta2_mean_sd = c(sd_safe(get_num(oracle_lambda, "beta2_mean")), sd_safe(get_num(oracle_lam_phi, "beta2_mean"))),
  r_mean_avg = c(mean_safe(get_num(oracle_lambda, "r_mean")), mean_safe(get_num(oracle_lam_phi, "r_mean"))),
  r_region_coverage_95_avg = c(mean_safe(get_num(oracle_lambda, "r_region_coverage_95")), mean_safe(get_num(oracle_lam_phi, "r_region_coverage_95"))),
  r_region_rmse_avg = c(mean_safe(get_num(oracle_lambda, "r_region_rmse")), mean_safe(get_num(oracle_lam_phi, "r_region_rmse"))),
  r_region_mae_avg = c(mean_safe(get_num(oracle_lambda, "r_region_mae")), mean_safe(get_num(oracle_lam_phi, "r_region_mae"))),
  phi_rmse_avg = c(mean_safe(get_num(oracle_lambda, "phi_rmse")), mean_safe(get_num(oracle_lam_phi, "phi_rmse"))),
  phi_cor_avg = c(mean_safe(get_num(oracle_lambda, "phi_cor")), mean_safe(get_num(oracle_lam_phi, "phi_cor"))),
  lambda_rmse_avg = c(mean_safe(get_num(oracle_lambda, "lambda_rmse")), mean_safe(get_num(oracle_lam_phi, "lambda_rmse"))),
  log_lambda_rmse_avg = c(mean_safe(get_num(oracle_lambda, "log_lambda_rmse")), mean_safe(get_num(oracle_lam_phi, "log_lambda_rmse"))),
  cor_log_lambda_avg = c(mean_safe(get_num(oracle_lambda, "cor_log_lambda")), mean_safe(get_num(oracle_lam_phi, "cor_log_lambda"))),
  kappa_truth_cv_avg = c(mean_safe(get_num(oracle_lambda, "kappa_truth_cv")), mean_safe(get_num(oracle_lam_phi, "kappa_truth_cv"))),
  kappa_post_mean_cv_avg = c(mean_safe(get_num(oracle_lambda, "kappa_post_mean_cv")), mean_safe(get_num(oracle_lam_phi, "kappa_post_mean_cv"))),
  kappa_rmse_avg = c(mean_safe(get_num(oracle_lambda, "kappa_rmse")), mean_safe(get_num(oracle_lam_phi, "kappa_rmse"))),
  log_kappa_rmse_avg = c(mean_safe(get_num(oracle_lambda, "log_kappa_rmse")), mean_safe(get_num(oracle_lam_phi, "log_kappa_rmse"))),
  cor_log_kappa_avg = c(mean_safe(get_num(oracle_lambda, "cor_log_kappa")), mean_safe(get_num(oracle_lam_phi, "cor_log_kappa"))),
  beta_guard_total = c(sum_safe(get_num(oracle_lambda, "s4c_oracle_beta_guard_count")), sum_safe(get_num(oracle_lam_phi, "s4c_oracle_beta_guard_count"))),
  kappa_guard_total = c(sum_safe(get_num(oracle_lambda, "s4c_oracle_kappa_guard_count")), sum_safe(get_num(oracle_lam_phi, "s4c_oracle_kappa_guard_count"))),
  lambda_guard_total = c(0, 0),
  stringsAsFactors = FALSE
)

sampled_stable_summary <- summarize_sampled(stable_df, "sampled_lambda_stable")
sampled_failure_summary <- summarize_sampled(failure_df, "sampled_lambda_failure")

# Build ladder directly from fixed vectors.
ladder <- data.frame(
  step_order = 1:4,
  step = c(
    "S4C sampled lambda, stable reps",
    "S4C sampled lambda, failure rep",
    "Oracle lambda, estimate beta/phi/r/kappa",
    "Oracle lambda + phi, estimate beta/r/kappa"
  ),
  n = c(
    sampled_stable_summary$n,
    sampled_failure_summary$n,
    oracle_summary$n[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$n[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  reps = c(
    sampled_stable_summary$reps,
    sampled_failure_summary$reps,
    oracle_summary$reps[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$reps[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  beta2 = c(
    sampled_stable_summary$beta2_mean_avg,
    sampled_failure_summary$beta2_mean_avg,
    oracle_summary$beta2_mean_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$beta2_mean_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  r_mean = c(
    sampled_stable_summary$r_mean_avg,
    sampled_failure_summary$r_mean_avg,
    oracle_summary$r_mean_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$r_mean_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  loglam_rmse = c(
    sampled_stable_summary$log_lambda_rmse_avg,
    sampled_failure_summary$log_lambda_rmse_avg,
    oracle_summary$log_lambda_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$log_lambda_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  phi_rmse = c(
    sampled_stable_summary$phi_rmse_avg,
    sampled_failure_summary$phi_rmse_avg,
    oracle_summary$phi_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$phi_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  kappa_post_cv = c(
    sampled_stable_summary$kappa_post_mean_cv_avg,
    sampled_failure_summary$kappa_post_mean_cv_avg,
    oracle_summary$kappa_post_mean_cv_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$kappa_post_mean_cv_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  log_kappa_rmse = c(
    sampled_stable_summary$log_kappa_rmse_avg,
    sampled_failure_summary$log_kappa_rmse_avg,
    oracle_summary$log_kappa_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$log_kappa_rmse_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  cor_log_kappa = c(
    sampled_stable_summary$cor_log_kappa_avg,
    sampled_failure_summary$cor_log_kappa_avg,
    oracle_summary$cor_log_kappa_avg[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$cor_log_kappa_avg[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  beta_guard = c(
    sampled_stable_summary$beta_guard_total,
    sampled_failure_summary$beta_guard_total,
    oracle_summary$beta_guard_total[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$beta_guard_total[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  kappa_guard = c(
    sampled_stable_summary$kappa_guard_total,
    sampled_failure_summary$kappa_guard_total,
    oracle_summary$kappa_guard_total[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$kappa_guard_total[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  lambda_guard = c(
    sampled_stable_summary$lambda_guard_total,
    sampled_failure_summary$lambda_guard_total,
    oracle_summary$lambda_guard_total[oracle_summary$diagnostic_label == "oracle_lambda"],
    oracle_summary$lambda_guard_total[oracle_summary$diagnostic_label == "oracle_lambda_phi"]
  ),
  interpretation = c(
    "Ordinary recovery summary for nonfailed S4C sampled-lambda fits.",
    "Problem replicate with sampled lambda; shows latent-process failure.",
    "Fixing lambda restores beta2 but leaves phi/kappa instability.",
    "Fixing lambda and phi restores beta/r and removes numerical guards."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  oracle_detail,
  file.path(out_dir, "scenario4c_oracle_diagnostic_detail.csv"),
  row.names = FALSE
)
write.csv(
  oracle_summary,
  file.path(out_dir, "scenario4c_oracle_diagnostic_summary.csv"),
  row.names = FALSE
)
write.csv(
  ladder,
  file.path(out_dir, "scenario4c_diagnostic_ladder_summary.csv"),
  row.names = FALSE
)

message("Wrote: ", file.path(out_dir, "scenario4c_oracle_diagnostic_detail.csv"))
message("Wrote: ", file.path(out_dir, "scenario4c_oracle_diagnostic_summary.csv"))
message("Wrote: ", file.path(out_dir, "scenario4c_diagnostic_ladder_summary.csv"))
message("Done.")
