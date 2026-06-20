## ============================================================================
## run_s4c_small_r_fixed_gamma_T100.R
##
## Scenario 4C strong-overdispersion stress test with gamma fixed.
##
## Purpose
##   Fit the same dynamic MSSTNB model used in Scenario 3 to Scenario 4C data,
##   where the negative-binomial dispersion truth is reduced from r = 15 to
##   r = 3.  The mean, exposure, covariates, spatial effects, and latent lambda
##   path are kept at the Scenario 3 reference structure.  This isolates strong
##   overdispersion / small-r stress from S4A sparse-count stress and S4B
##   low-exposure stress.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   s4c_small_r_overdispersion_v2.R  (preferred) or s4c_small_r_overdispersion.R
##   data_s4c_overdispersion/S4C_STRONG_OVERDISPERSION_T100/data_repXX.rds
##
## Main outputs
##   output_s4c_small_r_fixed_gamma/
##     S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100/
##
## Notes
##   This script only runs fitting and writes one-row summaries/manifests.
##   Posterior-performance analysis should be handled by a separate
##   run_scenario4c_posterior_performance.R script.
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100"
data_scenario_id <- "S4C_STRONG_OVERDISPERSION_T100"

## Start with short_test.  After rep 1 succeeds, change to "formal".
run_profile <- "formal"

if (identical(run_profile, "short_test")) {
    reps_formal <- 1:1
} else if (identical(run_profile, "formal")) {
    reps_formal <- 1:10
} else {
    stop("Unknown run_profile: ", run_profile, call. = FALSE)
}

s3_core_file <- "s3_dynamic_learned_gamma.R"
s4c_data_file <- if (file.exists(file.path(root_dir, "s4c_small_r_overdispersion_v2.R"))) {
    "s4c_small_r_overdispersion_v2.R"
} else {
    "s4c_small_r_overdispersion.R"
}

## Expected official S4C data-generation constants.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "small_r_overdispersion"
expected_r_reference_truth <- 15
expected_r_stress_truth <- 3

## MCMC settings.  These mirror Scenario 3/S4A/S4B fixed-gamma runs.
if (identical(run_profile, "short_test")) {
    n_iter <- 6000L
    n_burnin <- 1000L
    n_thin <- 5L
} else if (identical(run_profile, "formal")) {
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
}

fixed_gamma_value <- 0.8
gamma_prior <- c(1, 1)

data_dir <- "data_s4c_overdispersion"
output_dir <- "output_s4c_small_r_fixed_gamma"
verbose <- 1000L

## In short_test mode, overwrite by default.  In formal mode, skip existing fits
## unless explicitly set to TRUE.
overwrite_fit <- identical(run_profile, "short_test")

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir_s4c_fit <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s4c_fit <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true_s4c_fit <- function(x, message) {
    if (!isTRUE(x)) stop(message, call. = FALSE)
    invisible(TRUE)
}

cv_s4c_fit <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x, na.rm = TRUE)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) {
        return(NA_real_)
    }
    stats::sd(x, na.rm = TRUE) / mx
}

safe_ratio_s4c_fit <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}

count_stats_s4c_fit <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(yy, probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = cv_s4c_fit(yy),
        variance_to_mean = safe_ratio_s4c_fit(vv, mn),
        max_to_mean = safe_ratio_s4c_fit(max(yy), mn),
        q95_count = as.numeric(qs[6]),
        q99_count = as.numeric(qs[7]),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

kappa_stats_s4c_fit <- function(kappa, prefix = "") {
    kk <- as.numeric(kappa)
    qs <- stats::quantile(kk, probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        kappa_mean = mean(kk),
        kappa_sd = stats::sd(kk),
        kappa_cv = cv_s4c_fit(kk),
        kappa_min = min(kk),
        kappa_q01 = as.numeric(qs[1]),
        kappa_q05 = as.numeric(qs[2]),
        kappa_median = as.numeric(qs[3]),
        kappa_q95 = as.numeric(qs[4]),
        kappa_q99 = as.numeric(qs[5]),
        kappa_max = max(kk),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

s4c_source_data_path <- function(rep_id,
                                 root_arg = root_dir,
                                 data_dir_arg = data_dir,
                                 data_scenario_id_arg = data_scenario_id) {
    file.path(
        root_arg,
        data_dir_arg,
        data_scenario_id_arg,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

summarise_s4c_kappa_recovery <- function(kappa_draws, kappa_truth) {
    if (is.null(kappa_draws) || length(dim(kappa_draws)) != 3L) {
        return(data.frame(
            kappa_rmse = NA_real_,
            log_kappa_rmse = NA_real_,
            cor_kappa = NA_real_,
            cor_log_kappa = NA_real_,
            kappa_truth_cv = NA_real_,
            kappa_post_mean_cv = NA_real_,
            stringsAsFactors = FALSE
        ))
    }
    kappa_mean <- apply(kappa_draws, c(2, 3), mean, na.rm = TRUE)
    kappa_truth <- as.matrix(kappa_truth)
    if (!identical(dim(kappa_mean), dim(kappa_truth))) {
        stop("kappa posterior mean and kappa truth have incompatible dimensions.", call. = FALSE)
    }
    log_km <- log(pmax(kappa_mean, .Machine$double.xmin))
    log_kt <- log(pmax(kappa_truth, .Machine$double.xmin))
    data.frame(
        kappa_rmse = sqrt(mean((kappa_mean - kappa_truth)^2)),
        log_kappa_rmse = sqrt(mean((log_km - log_kt)^2)),
        cor_kappa = suppressWarnings(stats::cor(as.numeric(kappa_mean), as.numeric(kappa_truth))),
        cor_log_kappa = suppressWarnings(stats::cor(as.numeric(log_km), as.numeric(log_kt))),
        kappa_truth_cv = cv_s4c_fit(kappa_truth),
        kappa_post_mean_cv = cv_s4c_fit(kappa_mean),
        stringsAsFactors = FALSE
    )
}

r_recovery_s4c_fit <- function(r_draws, r_truth) {
    if (is.null(r_draws)) {
        return(data.frame(
            r_region_coverage_95 = NA_real_,
            r_region_rmse = NA_real_,
            r_region_mae = NA_real_,
            r_region_mean_sd = NA_real_,
            stringsAsFactors = FALSE
        ))
    }
    r_truth <- as.numeric(r_truth)
    if (length(r_truth) == 1L) r_truth <- rep(r_truth, ncol(r_draws))
    r_mean <- colMeans(r_draws, na.rm = TRUE)
    r_q <- apply(r_draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    r_sd <- apply(r_draws, 2, stats::sd, na.rm = TRUE)
    data.frame(
        r_region_coverage_95 = mean(r_q[1, ] <= r_truth & r_truth <= r_q[3, ]),
        r_region_rmse = sqrt(mean((r_mean - r_truth)^2)),
        r_region_mae = mean(abs(r_mean - r_truth)),
        r_region_mean_sd = mean(r_sd),
        stringsAsFactors = FALSE
    )
}

check_s4c_source_dataset <- function(data_file,
                                     TT_use_arg = 100L,
                                     n1_use_arg = 9L,
                                     expected_stress_type_arg = "small_r_overdispersion",
                                     expected_r_reference_truth_arg = 15,
                                     expected_r_stress_truth_arg = 3,
                                     beta0_ident_abs_limit = 20) {
    assert_file_exists_s4c_fit(data_file, "S4C source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "y_coarse", "y_coarse_reference", "e",
        "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "gamma_star",
        "beta0_star", "beta0_star_ident", "beta_star_ident", "phi_star_ident",
        "kappa", "kappa_reference", "r_star", "r_reference_truth", "r_stress_truth",
        "reference_count_summary", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s4c_fit(
        length(missing) == 0L,
        paste("Dataset is missing required S4C fields:", paste(missing, collapse = ", "))
    )

    assert_true_s4c_fit(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use_arg), as.integer(n1_use_arg))),
        paste0("y_coarse dimension is not TT_use by n1_use. Got ",
               paste(dim(dat$y_coarse), collapse = " x "), ".")
    )

    for (nm in c("y_coarse_reference", "e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident", "kappa", "kappa_reference")) {
        assert_true_s4c_fit(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    assert_true_s4c_fit(
        identical(dat$stress_type, expected_stress_type_arg),
        paste0("stress_type is not ", expected_stress_type_arg, ". Got ", dat$stress_type, ".")
    )
    assert_true_s4c_fit(
        abs(dat$r_reference_truth - expected_r_reference_truth_arg) < 1e-8,
        paste0("r_reference_truth is not ", expected_r_reference_truth_arg,
               ". Got ", dat$r_reference_truth, ".")
    )
    assert_true_s4c_fit(
        abs(dat$r_stress_truth - expected_r_stress_truth_arg) < 1e-8,
        paste0("r_stress_truth is not ", expected_r_stress_truth_arg,
               ". Got ", dat$r_stress_truth, ".")
    )

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$e must be positive finite; dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0) ||
        any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("lambda_tilde and lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$kappa)) || any(dat$kappa <= 0) ||
        any(!is.finite(dat$kappa_reference)) || any(dat$kappa_reference <= 0)) {
        stop("kappa and kappa_reference must be positive and finite.", call. = FALSE)
    }
    assert_true_s4c_fit(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident, ".")
    )

    if (exists("validate_s4c_overdispersion_data", mode = "function")) {
        validate_s4c_overdispersion_data(dat)
    }

    y_stats <- count_stats_s4c_fit(dat$y_coarse, prefix = "")
    y_ref_stats <- count_stats_s4c_fit(dat$y_coarse_reference, prefix = "reference_")
    k_stats <- kappa_stats_s4c_fit(dat$kappa, prefix = "")
    k_ref_stats <- kappa_stats_s4c_fit(dat$kappa_reference, prefix = "reference_")

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat$lambda_tilde_ident, finite = TRUE)

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type,
        r_reference_truth = dat$r_reference_truth,
        r_stress_truth = dat$r_stress_truth,
        r_ratio_to_reference = dat$r_stress_truth / dat$r_reference_truth,
        y_stats = y_stats,
        y_ref_stats = y_ref_stats,
        k_stats = k_stats,
        k_ref_stats = k_ref_stats,
        count_cv_increase = y_stats$count_cv - y_ref_stats$reference_count_cv,
        variance_to_mean_increase = y_stats$variance_to_mean - y_ref_stats$reference_variance_to_mean,
        kappa_cv_increase = k_stats$kappa_cv - k_ref_stats$reference_kappa_cv,
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_ident_min = lambda_ident_range[[1L]],
        lambda_ident_max = lambda_ident_range[[2L]]
    )
}

make_s4c_source_data_manifest <- function(reps,
                                          root_arg = root_dir,
                                          data_dir_arg = data_dir,
                                          data_scenario_id_arg = data_scenario_id) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s4c_source_data_path(
            rep_id = rep_id,
            root_arg = root_arg,
            data_dir_arg = data_dir_arg,
            data_scenario_id_arg = data_scenario_id_arg
        )
        chk <- check_s4c_source_dataset(data_file)
        data.frame(
            scenario_id = data_scenario_id_arg,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            r_reference_truth = chk$r_reference_truth,
            r_stress_truth = chk$r_stress_truth,
            r_ratio_to_reference = chk$r_ratio_to_reference,
            mean_count = chk$y_stats$mean_count,
            reference_mean_count = chk$y_ref_stats$reference_mean_count,
            median_count = chk$y_stats$median_count,
            zero_prop = chk$y_stats$zero_prop,
            reference_zero_prop = chk$y_ref_stats$reference_zero_prop,
            total_count = chk$y_stats$total_count,
            max_count = chk$y_stats$max_count,
            reference_max_count = chk$y_ref_stats$reference_max_count,
            count_sd = chk$y_stats$count_sd,
            count_cv = chk$y_stats$count_cv,
            reference_count_cv = chk$y_ref_stats$reference_count_cv,
            count_cv_increase = chk$count_cv_increase,
            variance_to_mean = chk$y_stats$variance_to_mean,
            reference_variance_to_mean = chk$y_ref_stats$reference_variance_to_mean,
            variance_to_mean_increase = chk$variance_to_mean_increase,
            q95_count = chk$y_stats$q95_count,
            q99_count = chk$y_stats$q99_count,
            kappa_cv = chk$k_stats$kappa_cv,
            reference_kappa_cv = chk$k_ref_stats$reference_kappa_cv,
            kappa_cv_increase = chk$kappa_cv_increase,
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

## ---- S4C fitting functions --------------------------------------------------
fit_s4c_small_r_fixed_gamma_one_rep <- function(rep_id,
                                                scenario_id = "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100",
                                                data_scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
                                                data_dir = "data_s4c_overdispersion",
                                                output_dir = "output_s4c_small_r_fixed_gamma",
                                                root = ".",
                                                settings_override = list(),
                                                priors = MCMC_PRIORS,
                                                spatial = build_s3_spatial(),
                                                gamma_init = NULL,
                                                fixed_gamma_value = 0.8,
                                                gamma_prior = c(1, 1),
                                                verbose = 1000L,
                                                save_result = TRUE,
                                                return_result = TRUE) {
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, data_scenario_id, paste0("data_rep", rr, ".rds"))
    chk <- check_s4c_source_dataset(dat_file)
    dat <- chk$dat

    validate_s3_data(dat)
    if (exists("validate_s4c_overdispersion_data", mode = "function")) {
        validate_s4c_overdispersion_data(dat)
    }

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    gamma_start <- gamma_init %||% fixed_gamma_value
    if (length(gamma_start) == 1L) gamma_start <- rep(gamma_start, dat$n1)

    cat(sprintf("=== Scenario 4C strong-overdispersion fixed-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Stress    : %s; r truth = %.3f (reference %.3f)\n",
                dat$stress_type, dat$r_stress_truth, dat$r_reference_truth))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f, count CV = %.3f\n",
                chk$y_stats$mean_count, chk$y_stats$zero_prop, chk$y_stats$count_cv))
    cat(sprintf("Reference : mean = %.3f, count CV = %.3f\n",
                chk$y_ref_stats$reference_mean_count, chk$y_ref_stats$reference_count_cv))
    cat(sprintf("Fixed     : gamma = %.3f\n", mean(gamma_start)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    reset_s4c_numeric_guards()

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = verbose
    )

    s4c_guard_counts <- get_s4c_numeric_guards()
    fit$diagnostics$s4c_numeric_guards <- s4c_guard_counts
    fit$diagnostics$s4c_beta_guard_count <- s4c_guard_counts$n_beta_guard
    fit$diagnostics$s4c_kappa_guard_count <- s4c_guard_counts$n_kappa_guard
    fit$diagnostics$s4c_lambda_input_guard_count <- s4c_guard_counts$n_lambda_input_guard
    fit$diagnostics$s4c_lambda_output_guard_count <- s4c_guard_counts$n_lambda_output_guard

    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    fit$diagnostics$gamma_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    fit$metadata$model <- "S4C strong-overdispersion dynamic NB-ICAR with fixed gamma"
    fit$metadata$fixed_gamma <- TRUE
    fit$metadata$learned_gamma <- FALSE
    fit$metadata$gamma_fixed_value <- fixed_gamma_value
    fit$metadata$updated_blocks <- setdiff(fit$metadata$updated_blocks, "gamma")
    fit$metadata$disabled_blocks <- unique(c(fit$metadata$disabled_blocks, "gamma", "delta", "omega"))

    summary <- summarise_s3_dynamic_learned_gamma_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id)
    )

    gamma_truth_mean <- mean(dat$gamma_star %||% fixed_gamma_value, na.rm = TRUE)
    gamma_fixed_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    summary$gamma_fixed_in_fit <- TRUE
    summary$gamma_learned_in_fit <- FALSE
    summary$gamma_truth_mean <- gamma_truth_mean
    summary$gamma_mean <- gamma_fixed_mean
    summary$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    summary$gamma_q025 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.025, na.rm = TRUE))
    summary$gamma_q50 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.500, na.rm = TRUE))
    summary$gamma_q975 <- as.numeric(stats::quantile(fit$samples$gamma_common, 0.975, na.rm = TRUE))
    summary$gamma_bias <- gamma_fixed_mean - gamma_truth_mean
    summary$gamma_covered <- as.integer(summary$gamma_q025 <= gamma_truth_mean && gamma_truth_mean <= summary$gamma_q975)
    summary$gamma_accept_rate <- NA_real_
    summary$gamma_proposal_sd_final <- NA_real_

    ## S4C-specific data context.
    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$r_reference_truth <- dat$r_reference_truth %||% NA_real_
    summary$r_stress_truth <- dat$r_stress_truth %||% mean(dat$r_star %||% NA_real_)
    summary$r_ratio_to_reference <- summary$r_stress_truth / summary$r_reference_truth
    summary$reference_mean_count <- chk$y_ref_stats$reference_mean_count
    summary$reference_zero_prop <- chk$y_ref_stats$reference_zero_prop
    summary$observed_mean_count <- chk$y_stats$mean_count
    summary$observed_zero_prop <- chk$y_stats$zero_prop
    summary$observed_total_count <- chk$y_stats$total_count
    summary$observed_max_count <- chk$y_stats$max_count
    summary$count_sd <- chk$y_stats$count_sd
    summary$count_cv <- chk$y_stats$count_cv
    summary$reference_count_cv <- chk$y_ref_stats$reference_count_cv
    summary$count_cv_increase <- chk$count_cv_increase
    summary$variance_to_mean <- chk$y_stats$variance_to_mean
    summary$reference_variance_to_mean <- chk$y_ref_stats$reference_variance_to_mean
    summary$variance_to_mean_increase <- chk$variance_to_mean_increase
    summary$q95_count <- chk$y_stats$q95_count
    summary$q99_count <- chk$y_stats$q99_count
    summary$reference_max_count <- chk$y_ref_stats$reference_max_count
    summary$kappa_truth_cv <- chk$k_stats$kappa_cv
    summary$reference_kappa_cv <- chk$k_ref_stats$reference_kappa_cv
    summary$kappa_cv_increase <- chk$kappa_cv_increase

    r_extra <- r_recovery_s4c_fit(fit$samples$r, dat$r_star)
    for (nm in names(r_extra)) summary[[nm]] <- r_extra[[nm]]

    kappa_extra <- summarise_s4c_kappa_recovery(fit$samples$kappa, dat$kappa)
    for (nm in names(kappa_extra)) summary[[nm]] <- kappa_extra[[nm]]

    summary$s4c_beta_guard_count <- fit$diagnostics$s4c_beta_guard_count
    summary$s4c_kappa_guard_count <- fit$diagnostics$s4c_kappa_guard_count
    summary$s4c_lambda_input_guard_count <- fit$diagnostics$s4c_lambda_input_guard_count
    summary$s4c_lambda_output_guard_count <- fit$diagnostics$s4c_lambda_output_guard_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            stress_type = dat$stress_type,
            r_reference_truth = dat$r_reference_truth,
            r_stress_truth = dat$r_stress_truth,
            model_source = "S3_DYNAMIC_FIXED_GAMMA",
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            fit_file_prefix = "fit_S4C_small_r_fixed_gamma_rep",
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_fixed_value = fixed_gamma_value,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        ensure_dir_s4c_fit(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4C_small_r_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4C_small_r_fixed_gamma_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s4c_small_r_fixed_gamma_batch <- function(reps = 1:10,
                                              scenario_id = "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100",
                                              data_scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
                                              data_dir = "data_s4c_overdispersion",
                                              output_dir = "output_s4c_small_r_fixed_gamma",
                                              root = ".",
                                              settings_override = list(),
                                              gamma_init = NULL,
                                              fixed_gamma_value = 0.8,
                                              gamma_prior = c(1, 1),
                                              verbose = 1000L,
                                              overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    ensure_dir_s4c_fit(out_dir)

    summaries <- list()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4C_small_r_fixed_gamma_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s4c_small_r_fixed_gamma_one_rep(
            rep_id = rep_id,
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            root = root,
            settings_override = settings_override,
            gamma_init = gamma_init,
            fixed_gamma_value = fixed_gamma_value,
            gamma_prior = gamma_prior,
            verbose = verbose,
            save_result = TRUE,
            return_result = TRUE
        )
        summaries[[rr]] <- fit$summary
    }

    summary_all <- do.call(rbind, summaries)
    summary_file <- file.path(out_dir, "summary_S4C_small_r_fixed_gamma_all_reps.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)
    invisible(summary_all)
}

print_s4c_fit_driver_summary <- function(source_manifest, fit_summary, fit_manifest) {
    cat("\n=== Scenario 4C strong-overdispersion fixed-gamma fitting summary ===\n")
    cat("Number of reps: ", nrow(source_manifest), "\n", sep = "")
    cat("Mean count avg: ", round(mean(source_manifest$mean_count, na.rm = TRUE), 4), "\n", sep = "")
    cat("Zero prop avg : ", round(mean(source_manifest$zero_prop, na.rm = TRUE), 4), "\n", sep = "")
    cat("Count CV avg  : ", round(mean(source_manifest$count_cv, na.rm = TRUE), 4), "\n", sep = "")
    cat("Reference CV avg: ", round(mean(source_manifest$reference_count_cv, na.rm = TRUE), 4), "\n", sep = "")
    cat("r truth / reference: ", unique(source_manifest$r_stress_truth), " / ", unique(source_manifest$r_reference_truth), "\n", sep = "")
    cat("Fit files all present: ", all(fit_manifest$fit_exists), "\n", sep = "")

    if (!is.null(fit_summary) && nrow(fit_summary) > 0L) {
        beta_cols <- intersect(c("beta0_mean", "beta1_mean", "beta2_mean", "r_mean", "r_region_coverage_95", "lambda_rmse", "cor_log_lambda", "s4c_kappa_guard_count"), names(fit_summary))
        if (length(beta_cols) > 0L) {
            cat("\nSelected fit-summary columns:\n")
            print(fit_summary[, c("scenario_id", "rep_id", beta_cols), drop = FALSE])
        }
    }
    invisible(TRUE)
}

## ---- pre-flight checks ------------------------------------------------------
assert_file_exists_s4c_fit(file.path(root_dir, s3_core_file), "Scenario 3 core script")
assert_file_exists_s4c_fit(file.path(root_dir, s4c_data_file), "Scenario 4C data-generation script")

source(file.path(root_dir, s4c_data_file))
source_s4c_small_r_overdispersion(root = root_dir, verbose = FALSE)

## ---- S4C numerical guard for xi --------------------------------------------
## S4C has high counts and stronger kappa variation.  We use the same safe xi
## computation as S4A/S4B so overflow is recorded through guards instead of
## crashing the run.
s4c_log_xi_lower <- -40
s4c_log_xi_upper <-  40

compute_s3_xi <- function(e, x1, x2, beta0, beta, phi) {
    if (!is.matrix(e) || !is.matrix(x1) || !is.matrix(x2)) {
        stop("e, x1, and x2 must be matrices.", call. = FALSE)
    }
    if (!identical(dim(e), dim(x1)) || !identical(dim(e), dim(x2))) {
        stop("e, x1, and x2 must have the same dimensions.", call. = FALSE)
    }
    if (length(beta) < 2L || length(phi) != ncol(e)) {
        stop("beta must have length at least 2 and phi must match ncol(e).", call. = FALSE)
    }
    if (any(!is.finite(e)) || any(e <= 0)) {
        stop("Exposure matrix e must be positive and finite.", call. = FALSE)
    }

    TT_now <- nrow(e)
    n1_now <- ncol(e)
    xi <- matrix(NA_real_, TT_now, n1_now)
    for (j in seq_len(n1_now)) {
        eta_j <- beta0 + beta[1] * x1[, j] + beta[2] * x2[, j] + phi[j]
        log_xi_j <- log(e[, j]) + eta_j
        log_xi_j <- pmin(pmax(log_xi_j, s4c_log_xi_lower), s4c_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("S4C safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

cat(sprintf(
    "Using S4C safe compute_s3_xi(): log(xi) clamped to [%.1f, %.1f].\n",
    s4c_log_xi_lower, s4c_log_xi_upper
))

## ---- S4C stabilization guards ----------------------------------------------
.s4c_guard_env <- new.env(parent = emptyenv())
reset_s4c_numeric_guards <- function() {
    .s4c_guard_env$n_beta_guard <- 0L
    .s4c_guard_env$n_kappa_guard <- 0L
    .s4c_guard_env$n_lambda_input_guard <- 0L
    .s4c_guard_env$n_lambda_output_guard <- 0L
    invisible(TRUE)
}
get_s4c_numeric_guards <- function() {
    list(
        n_beta_guard = .s4c_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4c_guard_env$n_kappa_guard %||% 0L,
        n_lambda_input_guard = .s4c_guard_env$n_lambda_input_guard %||% 0L,
        n_lambda_output_guard = .s4c_guard_env$n_lambda_output_guard %||% 0L
    )
}
reset_s4c_numeric_guards()

s4c_beta0_bounds <- c(-30, 10)
s4c_beta_bounds <- c(-5, 5)
s4c_kappa_bounds <- c(1e-10, 1e10)
s4c_lambda_bounds <- c(1e-10, 1e10)

.update_beta_s3_original <- update_beta
update_beta <- function(beta_current, ...) {
    res <- .update_beta_s3_original(beta_current = beta_current, ...)
    bad <- FALSE
    smp <- res$sample
    if (length(smp) < 3L || any(!is.finite(smp))) {
        bad <- TRUE
    } else {
        bad <- smp[1] < s4c_beta0_bounds[1] || smp[1] > s4c_beta0_bounds[2] ||
            any(smp[-1] < s4c_beta_bounds[1] | smp[-1] > s4c_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s4c_guard_env$n_beta_guard <- (.s4c_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s4c_guard_rejected <- TRUE
    } else {
        res$s4c_guard_rejected <- FALSE
    }
    res
}

.update_kappa_s3_original <- update_kappa
update_kappa <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
    y <- as.matrix(y_coarse)
    L <- as.matrix(lambda_tilde)
    X <- as.matrix(xi)
    if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
        stop("S4C safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
    }
    r_vec <- as.numeric(r)
    if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
    if (length(r_vec) != ncol(y)) {
        stop("S4C safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
    }
    R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
    shape <- y + R
    rate <- X * L + R
    guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
    if (any(guard)) {
        .s4c_guard_env$n_kappa_guard <- (.s4c_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s4c_guard_env$n_kappa_guard <- (.s4c_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s4c_kappa_bounds[1]), s4c_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s4c_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

.ffbs_lambda_all_s3_original <- ffbs_lambda_all
ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0, return_diag = TRUE) {
    X <- as.matrix(xi)
    K <- as.matrix(kappa)
    bad_in <- !is.finite(X) | X <= 0 | !is.finite(K) | K <= 0
    if (any(bad_in)) {
        .s4c_guard_env$n_lambda_input_guard <- (.s4c_guard_env$n_lambda_input_guard %||% 0L) + sum(bad_in)
    }
    X <- pmin(pmax(X, 1e-10), 1e10)
    K <- pmin(pmax(K, s4c_kappa_bounds[1]), s4c_kappa_bounds[2])
    out <- .ffbs_lambda_all_s3_original(
        gamma = gamma,
        y_coarse = y_coarse,
        xi = X,
        kappa = K,
        a0 = a0,
        b0 = b0,
        return_diag = return_diag
    )
    L <- out$lambda_tilde
    bad_out <- !is.finite(L) | L <= 0 | L < s4c_lambda_bounds[1] | L > s4c_lambda_bounds[2]
    if (any(bad_out)) {
        .s4c_guard_env$n_lambda_output_guard <- (.s4c_guard_env$n_lambda_output_guard %||% 0L) + sum(bad_out)
        L <- pmin(pmax(L, s4c_lambda_bounds[1]), s4c_lambda_bounds[2])
        out$lambda_tilde <- L
        if (!is.null(out$diag)) {
            out$diag$min_lambda <- min(L)
            out$diag$max_lambda <- max(L)
            out$diag$s4c_lambda_output_guarded <- TRUE
        }
    }
    out
}

cat(sprintf(
    "Using S4C stabilization guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s4c_beta0_bounds[1], s4c_beta0_bounds[2],
    s4c_beta_bounds[1], s4c_beta_bounds[2],
    s4c_log_xi_lower, s4c_log_xi_upper
))

## ---- S4C fixed-gamma update override ---------------------------------------
update_gamma_common_s3 <- function(gamma_current,
                                   lambda_tilde,
                                   y_coarse,
                                   a0 = 10,
                                   gamma_prior = c(1, 1),
                                   proposal_sd = 0.15) {
    gamma_current <- as.numeric(gamma_current)
    gamma_common_current <- if (length(gamma_current) == 1L) gamma_current else mean(gamma_current)
    gamma_common_current <- min(max(gamma_common_current, 1e-12), 1 - 1e-12)
    gamma_new <- rep(gamma_common_current, ncol(lambda_tilde))
    list(
        gamma = gamma_new,
        gamma_common = gamma_common_current,
        accept = FALSE,
        log_alpha = NA_real_,
        proposal_sd = proposal_sd,
        log_target_current = NA_real_,
        log_target_proposal = NA_real_
    )
}
cat(sprintf("Using S4C fixed-gamma override: gamma fixed at %.3f.\n", fixed_gamma_value))

source_data_manifest <- make_s4c_source_data_manifest(
    reps = reps_formal,
    root_arg = root_dir,
    data_dir_arg = data_dir,
    data_scenario_id_arg = data_scenario_id
)

cat("\n=== S4C source-data check ===\n")
print(source_data_manifest[, c(
    "rep_id", "r_stress_truth", "mean_count", "median_count", "zero_prop",
    "total_count", "max_count", "count_cv", "reference_count_cv",
    "variance_to_mean", "reference_variance_to_mean", "q99_count",
    "kappa_cv", "reference_kappa_cv", "beta0_star_ident", "lambda_raw_min", "lambda_raw_max"
), drop = FALSE])

## ---- run Scenario 4C fixed-gamma fit ---------------------------------------
fit_summary <- fit_s4c_small_r_fixed_gamma_batch(
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
    fixed_gamma_value = fixed_gamma_value,
    gamma_prior = gamma_prior,
    verbose = verbose,
    overwrite_existing = overwrite_fit
)

fit_dir_full <- file.path(root_dir, output_dir, scenario_id)
fit_files <- file.path(
    fit_dir_full,
    sprintf("fit_S4C_small_r_fixed_gamma_rep%02d.rds", as.integer(reps_formal))
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
assert_true_s4c_fit(all(fit_manifest$fit_exists), "At least one S4C fit file was not created.")

ensure_dir_s4c_fit(fit_dir_full)
utils::write.csv(source_data_manifest, file.path(fit_dir_full, "s4c_source_data_manifest.csv"), row.names = FALSE)
utils::write.csv(fit_manifest, file.path(fit_dir_full, "s4c_fit_manifest.csv"), row.names = FALSE)
saveRDS(
    list(
        run_profile = run_profile,
        scenario_id = scenario_id,
        data_scenario_id = data_scenario_id,
        reps = reps_formal,
        mcmc = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin),
        fixed_gamma_value = fixed_gamma_value,
        gamma_prior = gamma_prior,
        source_data_manifest = source_data_manifest,
        fit_manifest = fit_manifest,
        fit_summary = fit_summary
    ),
    file.path(fit_dir_full, "run_info_S4C_small_r_fixed_gamma_T100.rds")
)

## ---- final summary ----------------------------------------------------------
print_s4c_fit_driver_summary(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
)

cat("\n=== Main output locations ===\n")
cat("S4C data: ", file.path(data_dir, data_scenario_id), "\n", sep = "")
cat("Fits    : ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("\nScenario 4C strong-overdispersion fixed-gamma T100 fitting finished successfully.\n")

invisible(list(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
))
