# Robust batch rerun + summary for S1 replications (path-safe v3)
# This version does NOT assume the current working directory is the project root.
# You provide project_dir, and the script will locate and source required files.

safe_num <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(x[[1]])
}

safe_mean <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(mean(x, na.rm = TRUE))
}

safe_quant_cover <- function(samples, truth, probs = c(0.025, 0.975)) {
  if (is.null(samples) || length(samples) == 0L || any(dim(samples) == 0L)) return(NA_real_)
  if (is.null(truth) || length(truth) == 0L) return(NA_real_)
  lo <- apply(samples, seq_along(dim(samples))[-1], quantile, probs = probs[1], na.rm = TRUE)
  hi <- apply(samples, seq_along(dim(samples))[-1], quantile, probs = probs[2], na.rm = TRUE)
  mean(truth >= lo & truth <= hi)
}

safe_cor <- function(x, y) {
  if (is.null(x) || is.null(y) || length(x) == 0L || length(y) == 0L) return(NA_real_)
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  suppressWarnings(as.numeric(cor(x, y, use = "pairwise.complete.obs")))
}

extract_col <- function(mat, j) {
  if (is.null(mat)) return(NULL)
  if (is.vector(mat)) {
    if (j == 1L) return(as.numeric(mat)) else return(NULL)
  }
  if (is.matrix(mat) || is.data.frame(mat)) {
    if (ncol(mat) < j) return(NULL)
    return(as.numeric(mat[, j]))
  }
  NULL
}

check_param <- function(samples, truth) {
  if (is.null(samples) || length(samples) == 0L || is.na(truth)) {
    return(list(mean = NA_real_, sd = NA_real_, lo = NA_real_, hi = NA_real_, ok = FALSE))
  }
  pm <- mean(samples)
  psd <- sd(samples)
  ci <- quantile(samples, c(0.025, 0.975), na.rm = TRUE)
  ok <- isTRUE(truth >= ci[1] && truth <= ci[2])
  list(mean = pm, sd = psd, lo = ci[1], hi = ci[2], ok = ok)
}

empty_row <- function(rep_id, fit_file, status = "unknown", error_msg = NA_character_) {
  data.frame(
    rep_id = rep_id,
    fit_file = fit_file,
    status = status,
    error_msg = error_msg,
    beta0_mean = NA_real_,
    beta1_mean = NA_real_,
    beta2_mean = NA_real_,
    gamma1_mean = NA_real_,
    r1_mean = NA_real_,
    delta_mean = NA_real_,
    phi_cor = NA_real_,
    lambda_log_cor = NA_real_,
    lambda_log_rmse = NA_real_,
    lambda_cover = NA_real_,
    phi_accept = NA_real_,
    gamma_accept = NA_real_,
    delta_accept = NA_real_,
    r_accept = NA_real_,
    beta_mean_n_reject = NA_real_,
    beta0_ok = FALSE,
    beta1_ok = FALSE,
    beta2_ok = FALSE,
    gamma1_ok = FALSE,
    r1_ok = FALSE,
    delta_ok = FALSE,
    phi_ok = FALSE,
    lambda_cor_ok = FALSE,
    lambda_cover_ok = FALSE,
    n_pass = 0L,
    n_total = 9L,
    stage3_pass = FALSE,
    stringsAsFactors = FALSE
  )
}

summarise_one_fit <- function(rep_id, dat, fit, fit_file) {
  row <- empty_row(rep_id, fit_file, status = "ok", error_msg = NA_character_)

  s <- fit$samples
  if (is.null(s)) {
    row$status <- "missing_samples"
    row$error_msg <- "fit$samples missing"
    return(row)
  }

  beta0_true <- safe_num(dat$beta0_star_ident)
  beta1_true <- safe_num(dat$beta_star[1])
  beta2_true <- safe_num(dat$beta_star[2])
  gamma1_true <- safe_num(dat$gamma_star[1])
  r1_true <- safe_num(dat$r_star[1])
  delta_true <- safe_num(dat$delta_star)

  beta0 <- check_param(s$beta0, beta0_true)
  beta1 <- check_param(extract_col(s$beta, 1L), beta1_true)
  beta2 <- check_param(extract_col(s$beta, 2L), beta2_true)
  gamma1 <- check_param(extract_col(s$gamma, 1L), gamma1_true)
  r1 <- check_param(extract_col(s$r, 1L), r1_true)
  delta <- check_param(s$delta, delta_true)

  phi_pm <- NULL
  if (!is.null(s$phi)) {
    if (is.matrix(s$phi)) phi_pm <- colMeans(s$phi)
    if (is.data.frame(s$phi)) phi_pm <- colMeans(as.matrix(s$phi))
  }
  phi_cor <- safe_cor(phi_pm, dat$phi_star_ident)

  lt_cor <- NA_real_
  lt_rmse <- NA_real_
  lt_cover <- NA_real_
  if (!is.null(s$lambda_tilde)) {
    lt_arr <- s$lambda_tilde
    if (length(dim(lt_arr)) == 3L) {
      pm <- apply(lt_arr, c(2, 3), mean)
      truth <- log(dat$lambda_tilde_ident)
      lt_cor <- safe_cor(as.vector(log(pm)), as.vector(truth))
      lt_rmse <- sqrt(mean((as.vector(log(pm)) - as.vector(truth))^2))
      lt_cover <- safe_quant_cover(log(lt_arr), truth)
    }
  }

  row$beta0_mean <- beta0$mean
  row$beta1_mean <- beta1$mean
  row$beta2_mean <- beta2$mean
  row$gamma1_mean <- gamma1$mean
  row$r1_mean <- r1$mean
  row$delta_mean <- delta$mean
  row$phi_cor <- phi_cor
  row$lambda_log_cor <- lt_cor
  row$lambda_log_rmse <- lt_rmse
  row$lambda_cover <- lt_cover
  row$phi_accept <- safe_num(fit$diagnostics$phi_accept_rate)
  row$gamma_accept <- safe_mean(fit$diagnostics$gamma_accept_rate)
  row$delta_accept <- safe_num(fit$diagnostics$delta_accept_rate)
  row$r_accept <- safe_mean(fit$diagnostics$r_accept_rate)
  row$beta_mean_n_reject <- safe_num(fit$diagnostics$beta_mean_n_reject)

  row$beta0_ok <- beta0$ok
  row$beta1_ok <- beta1$ok
  row$beta2_ok <- beta2$ok
  row$gamma1_ok <- gamma1$ok
  row$r1_ok <- r1$ok
  row$delta_ok <- delta$ok
  row$phi_ok <- isTRUE(!is.na(phi_cor) && phi_cor > 0.5)
  row$lambda_cor_ok <- isTRUE(!is.na(lt_cor) && lt_cor > 0.3)
  row$lambda_cover_ok <- isTRUE(!is.na(lt_cover) && lt_cover > 0.70)

  checks <- unlist(row[c("beta0_ok", "beta1_ok", "beta2_ok", "gamma1_ok", "r1_ok", "delta_ok", "phi_ok", "lambda_cor_ok", "lambda_cover_ok")])
  row$n_pass <- sum(checks, na.rm = TRUE)
  row$stage3_pass <- isTRUE(row$n_pass >= row$n_total - 1L)
  row
}

resolve_project_path <- function(project_dir = ".") {
  normalizePath(project_dir, winslash = "/", mustWork = TRUE)
}

find_file_in_project <- function(project_dir, filename, subdirs = c("", "R", "R/mcmc", "mcmc", "scripts", "code")) {
  candidates <- unique(file.path(project_dir, subdirs, filename))
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) return(NA_character_)
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

source_required_files <- function(project_dir, verbose = TRUE) {
  required <- c(
    "00_setup.R",
    "mcmc_config.R",
    "mcmc_utils.R",
    "update_kappa.R",
    "update_regression.R",
    "update_icar.R",
    "update_dispersion.R",
    "update_gamma.R",
    "ffbs_lambda.R",
    "update_delta.R",
    "smooth_omega.R",
    "sampler.R"
  )
  paths <- setNames(lapply(required, function(f) find_file_in_project(project_dir, f)), required)
  missing <- names(paths)[vapply(paths, function(x) is.na(x) || !nzchar(x), logical(1))]
  if (length(missing) > 0L) {
    stop(
      "Could not find required project files under project_dir = ", project_dir,
      "\nMissing: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (verbose) {
    cat("\nSourcing project files from:\n")
    for (nm in names(paths)) cat("  ", nm, " <- ", paths[[nm]], "\n", sep = "")
    cat("\n")
  }
  for (nm in required) source(paths[[nm]], local = .GlobalEnv)
  invisible(paths)
}

run_batch_s1_summary <- function(reps = 1:10,
                                 project_dir = ".",
                                 data_dir = "data_revised",
                                 fit_dir = "batch_fits_S1",
                                 out_prefix = "batch_s1_m1_short",
                                 overwrite_existing = TRUE,
                                 verbose = TRUE) {
  project_dir <- resolve_project_path(project_dir)

  # make relative directories project-root relative
  if (!grepl("^(/|[A-Za-z]:)", data_dir)) data_dir <- file.path(project_dir, data_dir)
  if (!grepl("^(/|[A-Za-z]:)", fit_dir)) fit_dir <- file.path(project_dir, fit_dir)
  if (!grepl("^(/|[A-Za-z]:)", out_prefix)) out_prefix <- file.path(project_dir, out_prefix)

  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

  source_required_files(project_dir, verbose = verbose)

  settings_short <- MCMC_SETTINGS
  settings_short$n_iter <- 2000L
  settings_short$n_burnin <- 1000L
  settings_short$n_thin <- 2L

  priors <- MCMC_PRIORS
  spatial <- list(H = H, B_ICAR = B_ICAR, BHB = crossprod(B_ICAR, H %*% B_ICAR))
  constants <- list(A0 = A0, B0 = B0, C0 = C0)

  rows <- vector("list", length(reps))

  for (i in seq_along(reps)) {
    rep_id <- reps[i]
    rr <- sprintf("%02d", rep_id)
    dat_file <- file.path(data_dir, "S1", sprintf("data_rep%s.rds", rr))
    fit_file <- file.path(fit_dir, sprintf("fit_S1_rep%s_short.rds", rr))

    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf("S1 rep %s (%d of %d)\n", rr, i, length(reps)))
      cat("============================================================\n")
    }

    if (!file.exists(dat_file)) {
      warning(sprintf("Data file not found: %s", dat_file))
      rows[[i]] <- empty_row(rep_id, fit_file, status = "missing_data", error_msg = sprintf("Data file not found: %s", dat_file))
      next
    }

    dat <- readRDS(dat_file)

    if (!all(c("beta0_star_ident", "phi_star_ident", "lambda_tilde_ident") %in% names(dat))) {
      rows[[i]] <- empty_row(rep_id, fit_file, status = "missing_identified_truth", error_msg = "Dataset missing identified-truth fields")
      next
    }

    fit <- NULL
    if (file.exists(fit_file) && !overwrite_existing) {
      fit <- readRDS(fit_file)
    } else {
      fit_try <- try(run_mcmc(dat, settings_short, priors, spatial, constants, verbose = 500L), silent = TRUE)
      if (inherits(fit_try, "try-error")) {
        rows[[i]] <- empty_row(rep_id, fit_file, status = "fit_error", error_msg = as.character(fit_try))
        next
      }
      fit <- fit_try
      saveRDS(fit, fit_file)
    }

    rows[[i]] <- summarise_one_fit(rep_id, dat, fit, fit_file)
  }

  summary_df <- do.call(rbind, rows)
  csv_file <- paste0(out_prefix, "_summary.csv")
  txt_file <- paste0(out_prefix, "_summary.txt")
  write.csv(summary_df, csv_file, row.names = FALSE)

  n_ok <- sum(summary_df$status == "ok")
  n_pass <- sum(summary_df$stage3_pass, na.rm = TRUE)
  avg_beta1 <- mean(summary_df$beta1_mean, na.rm = TRUE)
  avg_gamma <- mean(summary_df$gamma1_mean, na.rm = TRUE)
  avg_rmse <- mean(summary_df$lambda_log_rmse, na.rm = TRUE)

  txt <- c(
    "============================================================",
    "Batch rerun summary: S1 ordinary fit (short chain)",
    "============================================================",
    sprintf("Project dir                  : %s", project_dir),
    sprintf("Replications requested       : %d", length(reps)),
    sprintf("Replications with data       : %d", sum(summary_df$status != "missing_data")),
    sprintf("Replications fit successfully: %d", n_ok),
    sprintf("Stage-3 PASS count           : %d", n_pass),
    sprintf("Mean beta1 posterior mean    : %.4f", avg_beta1),
    sprintf("Mean gamma1 posterior mean   : %.4f", avg_gamma),
    sprintf("Mean lambda log-RMSE         : %.4f", avg_rmse),
    "",
    "Status counts:",
    capture.output(print(table(summary_df$status, useNA = "ifany")))
  )
  writeLines(txt, txt_file)

  if (verbose) {
    cat("\n")
    cat(paste(txt, collapse = "\n"), "\n")
    cat(sprintf("\nWrote: %s\n", csv_file))
    cat(sprintf("Wrote: %s\n", txt_file))
  }

  invisible(list(summary = summary_df, csv_file = csv_file, txt_file = txt_file))
}
