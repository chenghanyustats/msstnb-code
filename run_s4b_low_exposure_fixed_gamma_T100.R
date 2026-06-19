## ============================================================================
## run_s4b_low_exposure_fixed_gamma_T100.R
##
## Scenario 4B low/heterogeneous-exposure stress test with gamma fixed.
##
## Purpose
##   Fit the same dynamic MSSTNB model used in Scenario 3 to Scenario 4B
##   low-heterogeneous-exposure data, but fix the common temporal discount
##   parameter gamma at its truth value (0.8).
##
## Why fixed gamma here?
##   S4B is designed to isolate the effect of known exposure being small and
##   heterogeneous.  Gamma learning is not the scientific target of this stress
##   test.  Fixing gamma at truth makes the comparison to S3 and S4A cleaner and
##   lets us focus on recovery of beta, phi, r, and lambda under exposure stress.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   s4b_low_exposure.R
##   data_s4b_low_exposure/S4B_LOW_HETEROGENEOUS_EXPOSURE_T100/data_repXX.rds
##
## Main outputs
##   output_s4b_low_exposure_fixed_gamma/
##     S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_T100/
##
## Notes
##   This script only runs fitting and writes one-row summaries/manifests.
##   Posterior-performance analysis should be handled by a separate
##   run_scenario4b_posterior_performance.R script.
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_T100"
data_scenario_id <- "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100"

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
s4b_data_file <- "s4b_low_exposure.R"

## Expected official S4B data-generation constants.
TT_use <- 100L
n1_use <- 9L
expected_stress_type <- "low_exposure"
expected_exposure_stress_type <- "low_heterogeneous_exposure"
expected_target_mean_multiplier <- 0.05
expected_area_log_sd <- 0.75
expected_time_log_sd <- 0.08
expected_lower_multiplier <- 0.005
expected_upper_multiplier <- 0.25

## MCMC settings.  These mirror Scenario 3/S4A fixed-gamma runs.
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

data_dir <- "data_s4b_low_exposure"
output_dir <- "output_s4b_low_exposure_fixed_gamma"
verbose <- 1000L

## In short_test mode, overwrite by default.  In formal mode, skip existing fits
## unless explicitly set to TRUE.
overwrite_fit <- identical(run_profile, "short_test")

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir_s4b_fit <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s4b_fit <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true_s4b_fit <- function(x, message) {
    if (!isTRUE(x)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

cv_s4b_fit <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0L || !is.finite(mean(x, na.rm = TRUE)) || mean(x, na.rm = TRUE) == 0) {
        return(NA_real_)
    }
    stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

s4b_source_data_path <- function(rep_id,
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

## Area groups by observed low-exposure mean.  Used both for source-data checks
## and posterior diagnostic summaries.
s4b_area_exposure_groups <- function(dat, n_groups = 3L) {
    n_groups <- as.integer(n_groups)
    if (n_groups < 2L) n_groups <- 2L

    area_mean_e <- colMeans(dat$e)
    ord <- order(area_mean_e, decreasing = FALSE)
    group_id <- rep(NA_integer_, length(area_mean_e))
    group_id[ord] <- cut(seq_along(ord), breaks = n_groups, labels = FALSE)

    if (n_groups == 3L) {
        labels <- c("low_exposure", "middle_exposure", "high_exposure")
    } else {
        labels <- paste0("group_", seq_len(n_groups))
    }
    list(group_id = group_id, labels = labels)
}

s4b_exposure_group_count_summary <- function(dat, n_groups = 3L) {
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        y_g <- dat$y_coarse[, jj, drop = FALSE]
        e_g <- dat$e[, jj, drop = FALSE]
        e_ref_g <- dat$e_reference[, jj, drop = FALSE]
        mult_g <- dat$exposure_multiplier[, jj, drop = FALSE]
        data.frame(
            exposure_group = g$labels[gg],
            n_areas = length(jj),
            area_ids = paste(jj, collapse = ","),
            mean_exposure = mean(e_g),
            mean_reference_exposure = mean(e_ref_g),
            mean_multiplier = mean(mult_g),
            mean_count = mean(y_g),
            zero_prop = mean(y_g == 0),
            total_count = sum(y_g),
            max_count = max(y_g),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

s4b_lambda_group_recovery <- function(lambda_draws, lambda_truth, dat, n_groups = 3L) {
    if (is.null(lambda_draws) || length(dim(lambda_draws)) != 3L) {
        return(data.frame())
    }
    lambda_mean <- apply(lambda_draws, c(2, 3), mean, na.rm = TRUE)
    lambda_truth <- as.matrix(lambda_truth)
    if (!identical(dim(lambda_mean), dim(lambda_truth))) {
        stop("lambda posterior mean and lambda truth have incompatible dimensions.", call. = FALSE)
    }
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        lm <- lambda_mean[, jj, drop = FALSE]
        lt <- lambda_truth[, jj, drop = FALSE]
        log_lm <- log(pmax(lm, .Machine$double.xmin))
        log_lt <- log(pmax(lt, .Machine$double.xmin))
        data.frame(
            exposure_group = g$labels[gg],
            lambda_rmse = sqrt(mean((lm - lt)^2)),
            log_lambda_rmse = sqrt(mean((log_lm - log_lt)^2)),
            cor_log_lambda = suppressWarnings(stats::cor(as.numeric(log_lm), as.numeric(log_lt))),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

s4b_phi_group_recovery <- function(phi_draws, phi_truth, dat, n_groups = 3L) {
    if (is.null(phi_draws) || length(dim(phi_draws)) != 2L) {
        return(data.frame())
    }
    phi_mean <- colMeans(phi_draws, na.rm = TRUE)
    phi_truth <- as.numeric(phi_truth)
    if (length(phi_mean) != length(phi_truth)) {
        stop("phi posterior mean and phi truth have incompatible lengths.", call. = FALSE)
    }
    g <- s4b_area_exposure_groups(dat, n_groups = n_groups)
    out <- lapply(seq_along(g$labels), function(gg) {
        jj <- which(g$group_id == gg)
        err <- phi_mean[jj] - phi_truth[jj]
        data.frame(
            exposure_group = g$labels[gg],
            phi_rmse = sqrt(mean(err^2)),
            phi_mae = mean(abs(err)),
            phi_cor = suppressWarnings(stats::cor(phi_mean[jj], phi_truth[jj])),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

prefix_group_metrics <- function(df, prefix) {
    if (is.null(df) || nrow(df) == 0L) return(list())
    out <- list()
    for (ii in seq_len(nrow(df))) {
        grp <- df$exposure_group[ii]
        for (nm in setdiff(names(df), "exposure_group")) {
            out[[paste0(prefix, "_", grp, "_", nm)]] <- df[[nm]][ii]
        }
    }
    out
}

check_s4b_source_dataset <- function(data_file,
                                     TT_use_arg = 100L,
                                     n1_use_arg = 9L,
                                     expected_stress_type_arg = "low_exposure",
                                     expected_exposure_stress_type_arg = "low_heterogeneous_exposure",
                                     expected_target_mean_multiplier_arg = 0.05,
                                     beta0_ident_abs_limit = 20) {
    assert_file_exists_s4b_fit(data_file, "S4B source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "scenario_id", "stress_type", "exposure_stress_type", "y_coarse", "e",
        "e_reference", "exposure_multiplier", "x1", "x2", "lambda_tilde",
        "lambda_tilde_ident", "gamma_star", "beta0_star", "beta0_star_ident",
        "beta_star_ident", "phi_star_ident", "target_mean_exposure_multiplier",
        "realized_mean_exposure_multiplier", "realized_min_exposure_multiplier",
        "realized_max_exposure_multiplier", "mean_count", "zero_prop", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s4b_fit(
        length(missing) == 0L,
        paste("Dataset is missing required S4B fields:", paste(missing, collapse = ", "))
    )

    assert_true_s4b_fit(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use_arg), as.integer(n1_use_arg))),
        paste0(
            "y_coarse dimension is not TT_use by n1_use. Got ",
            paste(dim(dat$y_coarse), collapse = " x "), "."
        )
    )

    for (nm in c("e", "e_reference", "exposure_multiplier", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        assert_true_s4b_fit(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    assert_true_s4b_fit(
        identical(dat$stress_type, expected_stress_type_arg),
        paste0("stress_type is not ", expected_stress_type_arg, ". Got ", dat$stress_type, ".")
    )
    assert_true_s4b_fit(
        identical(dat$exposure_stress_type, expected_exposure_stress_type_arg),
        paste0("exposure_stress_type is not ", expected_exposure_stress_type_arg, ". Got ", dat$exposure_stress_type, ".")
    )
    assert_true_s4b_fit(
        abs(dat$target_mean_exposure_multiplier - expected_target_mean_multiplier_arg) < 1e-8,
        paste0("target_mean_exposure_multiplier is not ", expected_target_mean_multiplier_arg,
               ". Got ", dat$target_mean_exposure_multiplier, ".")
    )

    if (!isTRUE(all.equal(dat$e, dat$e_reference * dat$exposure_multiplier,
                          tolerance = 1e-8, check.attributes = FALSE))) {
        stop("dat$e is not equal to dat$e_reference * dat$exposure_multiplier.", call. = FALSE)
    }
    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$e_reference)) || any(dat$e_reference <= 0) ||
        any(!is.finite(dat$exposure_multiplier)) || any(dat$exposure_multiplier <= 0) ||
        any(dat$exposure_multiplier >= 1)) {
        stop("S4B exposure values/multipliers must be positive, finite, and multipliers must be < 1.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0) ||
        any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("lambda_tilde and lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }
    assert_true_s4b_fit(
        is.finite(dat$beta0_star_ident) && abs(dat$beta0_star_ident) < beta0_ident_abs_limit,
        paste0("beta0_star_ident appears pathological: ", dat$beta0_star_ident, ".")
    )

    ## Run the stricter data-generation validator when available.
    if (exists("validate_s4b_low_exposure_data", mode = "function")) {
        validate_s4b_low_exposure_data(dat)
    }

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat$lambda_tilde_ident, finite = TRUE)
    count_groups <- s4b_exposure_group_count_summary(dat, n_groups = 3L)
    low_g <- count_groups[count_groups$exposure_group == "low_exposure", , drop = FALSE]
    high_g <- count_groups[count_groups$exposure_group == "high_exposure", , drop = FALSE]

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        stress_type = dat$stress_type,
        exposure_stress_type = dat$exposure_stress_type,
        mean_count = mean(dat$y_coarse),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        zero_prop = mean(dat$y_coarse == 0),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        reference_mean_count = dat$reference_count_summary$mean_count %||% NA_real_,
        reference_zero_prop = dat$reference_count_summary$zero_prop %||% NA_real_,
        target_mean_multiplier = dat$target_mean_exposure_multiplier,
        realized_mean_multiplier = mean(dat$exposure_multiplier),
        realized_min_multiplier = min(dat$exposure_multiplier),
        realized_max_multiplier = max(dat$exposure_multiplier),
        reference_mean_exposure = mean(dat$e_reference),
        mean_exposure = mean(dat$e),
        exposure_mean_ratio = mean(dat$e) / mean(dat$e_reference),
        exposure_cv = cv_s4b_fit(dat$e),
        multiplier_cv = cv_s4b_fit(dat$exposure_multiplier),
        area_multiplier_cv = cv_s4b_fit(colMeans(dat$exposure_multiplier)),
        lowest_exposure_group_mean_count = low_g$mean_count %||% NA_real_,
        lowest_exposure_group_zero_prop = low_g$zero_prop %||% NA_real_,
        highest_exposure_group_mean_count = high_g$mean_count %||% NA_real_,
        highest_exposure_group_zero_prop = high_g$zero_prop %||% NA_real_,
        count_exposure_cor_area = suppressWarnings(stats::cor(colMeans(dat$y_coarse), colMeans(dat$e))),
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_ident_min = lambda_ident_range[[1L]],
        lambda_ident_max = lambda_ident_range[[2L]]
    )
}

make_s4b_source_data_manifest <- function(reps,
                                          root_arg = root_dir,
                                          data_dir_arg = data_dir,
                                          data_scenario_id_arg = data_scenario_id) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s4b_source_data_path(
            rep_id = rep_id,
            root_arg = root_arg,
            data_dir_arg = data_dir_arg,
            data_scenario_id_arg = data_scenario_id_arg
        )
        chk <- check_s4b_source_dataset(data_file)
        data.frame(
            scenario_id = data_scenario_id_arg,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            stress_type = chk$stress_type,
            exposure_stress_type = chk$exposure_stress_type,
            mean_count = chk$mean_count,
            median_count = chk$median_count,
            zero_prop = chk$zero_prop,
            total_count = chk$total_count,
            max_count = chk$max_count,
            reference_mean_count = chk$reference_mean_count,
            reference_zero_prop = chk$reference_zero_prop,
            target_mean_multiplier = chk$target_mean_multiplier,
            realized_mean_multiplier = chk$realized_mean_multiplier,
            realized_min_multiplier = chk$realized_min_multiplier,
            realized_max_multiplier = chk$realized_max_multiplier,
            reference_mean_exposure = chk$reference_mean_exposure,
            mean_exposure = chk$mean_exposure,
            exposure_mean_ratio = chk$exposure_mean_ratio,
            exposure_cv = chk$exposure_cv,
            multiplier_cv = chk$multiplier_cv,
            area_multiplier_cv = chk$area_multiplier_cv,
            lowest_exposure_group_mean_count = chk$lowest_exposure_group_mean_count,
            lowest_exposure_group_zero_prop = chk$lowest_exposure_group_zero_prop,
            highest_exposure_group_mean_count = chk$highest_exposure_group_mean_count,
            highest_exposure_group_zero_prop = chk$highest_exposure_group_zero_prop,
            count_exposure_cor_area = chk$count_exposure_cor_area,
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

## ---- S4B fitting functions --------------------------------------------------
fit_s4b_low_exposure_fixed_gamma_one_rep <- function(rep_id,
                                                     scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_T100",
                                                     data_scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
                                                     data_dir = "data_s4b_low_exposure",
                                                     output_dir = "output_s4b_low_exposure_fixed_gamma",
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
    chk <- check_s4b_source_dataset(dat_file)
    dat <- chk$dat

    validate_s3_data(dat)
    if (exists("validate_s4b_low_exposure_data", mode = "function")) {
        validate_s4b_low_exposure_data(dat)
    }

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    gamma_start <- gamma_init %||% fixed_gamma_value
    if (length(gamma_start) == 1L) gamma_start <- rep(gamma_start, dat$n1)

    cat(sprintf("=== Scenario 4B low-exposure fixed-gamma fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Stress    : %s; target exposure multiplier = %.3f\n",
                dat$exposure_stress_type, dat$target_mean_exposure_multiplier))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Exposure  : mean ratio = %.4f, exposure CV = %.3f\n", chk$exposure_mean_ratio, chk$exposure_cv))
    cat(sprintf("Fixed     : gamma = %.3f\n", mean(gamma_start)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    reset_s4b_numeric_guards()

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = verbose
    )

    s4b_guard_counts <- get_s4b_numeric_guards()
    fit$diagnostics$s4b_numeric_guards <- s4b_guard_counts
    fit$diagnostics$s4b_beta_guard_count <- s4b_guard_counts$n_beta_guard
    fit$diagnostics$s4b_kappa_guard_count <- s4b_guard_counts$n_kappa_guard
    fit$diagnostics$s4b_lambda_input_guard_count <- s4b_guard_counts$n_lambda_input_guard
    fit$diagnostics$s4b_lambda_output_guard_count <- s4b_guard_counts$n_lambda_output_guard

    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    fit$diagnostics$gamma_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    fit$metadata$model <- "S4B low-heterogeneous-exposure dynamic NB-ICAR with fixed gamma"
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

    ## S4B-specific data context.
    summary$stress_type <- dat$stress_type %||% NA_character_
    summary$exposure_stress_type <- dat$exposure_stress_type %||% NA_character_
    summary$target_mean_exposure_multiplier <- dat$target_mean_exposure_multiplier %||% NA_real_
    summary$realized_mean_exposure_multiplier <- mean(dat$exposure_multiplier)
    summary$realized_min_exposure_multiplier <- min(dat$exposure_multiplier)
    summary$realized_max_exposure_multiplier <- max(dat$exposure_multiplier)
    summary$reference_mean_count <- chk$reference_mean_count
    summary$reference_zero_prop <- chk$reference_zero_prop
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count
    summary$reference_mean_exposure <- chk$reference_mean_exposure
    summary$mean_exposure <- chk$mean_exposure
    summary$exposure_mean_ratio <- chk$exposure_mean_ratio
    summary$exposure_cv <- chk$exposure_cv
    summary$multiplier_cv <- chk$multiplier_cv
    summary$area_multiplier_cv <- chk$area_multiplier_cv
    summary$lowest_exposure_group_mean_count <- chk$lowest_exposure_group_mean_count
    summary$lowest_exposure_group_zero_prop <- chk$lowest_exposure_group_zero_prop
    summary$highest_exposure_group_mean_count <- chk$highest_exposure_group_mean_count
    summary$highest_exposure_group_zero_prop <- chk$highest_exposure_group_zero_prop
    summary$count_exposure_cor_area <- chk$count_exposure_cor_area

    ## Exposure-group posterior recovery.  These fields will be useful in the
    ## later Scenario 4B posterior-performance script.
    lambda_group <- s4b_lambda_group_recovery(
        lambda_draws = fit$samples$lambda_tilde,
        lambda_truth = dat$lambda_tilde_ident %||% dat$lambda_tilde,
        dat = dat,
        n_groups = 3L
    )
    phi_group <- s4b_phi_group_recovery(
        phi_draws = fit$samples$phi,
        phi_truth = dat$phi_star_ident %||% dat$phi_star,
        dat = dat,
        n_groups = 3L
    )
    count_group <- s4b_exposure_group_count_summary(dat, n_groups = 3L)

    group_metrics <- c(
        prefix_group_metrics(count_group[, c("exposure_group", "n_areas", "mean_exposure", "mean_multiplier", "mean_count", "zero_prop", "total_count"), drop = FALSE], "count"),
        prefix_group_metrics(lambda_group, "lambda"),
        prefix_group_metrics(phi_group, "phi")
    )
    if (length(group_metrics) > 0L) {
        for (nm in names(group_metrics)) summary[[nm]] <- group_metrics[[nm]]
    }

    summary$s4b_beta_guard_count <- fit$diagnostics$s4b_beta_guard_count
    summary$s4b_kappa_guard_count <- fit$diagnostics$s4b_kappa_guard_count
    summary$s4b_lambda_input_guard_count <- fit$diagnostics$s4b_lambda_input_guard_count
    summary$s4b_lambda_output_guard_count <- fit$diagnostics$s4b_lambda_output_guard_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            data_scenario_id = data_scenario_id,
            stress_type = dat$stress_type,
            exposure_stress_type = dat$exposure_stress_type,
            target_mean_exposure_multiplier = dat$target_mean_exposure_multiplier,
            realized_mean_exposure_multiplier = mean(dat$exposure_multiplier),
            model_source = "S3_DYNAMIC_FIXED_GAMMA",
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            fit_file_prefix = "fit_S4B_low_exposure_fixed_gamma_rep",
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_fixed_value = fixed_gamma_value,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        ensure_dir_s4b_fit(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S4B_low_exposure_fixed_gamma_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S4B_low_exposure_fixed_gamma_rep", rr, ".csv"))
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

fit_s4b_low_exposure_fixed_gamma_batch <- function(reps = 1:10,
                                                   scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_T100",
                                                   data_scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
                                                   data_dir = "data_s4b_low_exposure",
                                                   output_dir = "output_s4b_low_exposure_fixed_gamma",
                                                   root = ".",
                                                   settings_override = list(),
                                                   gamma_init = NULL,
                                                   fixed_gamma_value = 0.8,
                                                   gamma_prior = c(1, 1),
                                                   verbose = 1000L,
                                                   overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    ensure_dir_s4b_fit(out_dir)

    summaries <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S4B_low_exposure_fixed_gamma_rep", rr, ".rds"))

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s4b_low_exposure_fixed_gamma_one_rep(
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
    summary_file <- file.path(out_dir, "summary_S4B_low_exposure_fixed_gamma_all_reps.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)

    invisible(summary_all)
}

print_s4b_fit_driver_summary <- function(source_manifest, fit_summary, fit_manifest) {
    cat("\n=== Scenario 4B low-heterogeneous-exposure fixed-gamma fitting summary ===\n")
    cat("Number of reps: ", nrow(source_manifest), "\n", sep = "")
    cat("Mean count avg: ", round(mean(source_manifest$mean_count, na.rm = TRUE), 4), "\n", sep = "")
    cat("Zero prop avg : ", round(mean(source_manifest$zero_prop, na.rm = TRUE), 4), "\n", sep = "")
    cat("Exposure ratio avg: ", round(mean(source_manifest$exposure_mean_ratio, na.rm = TRUE), 5), "\n", sep = "")
    cat("Exposure CV avg   : ", round(mean(source_manifest$exposure_cv, na.rm = TRUE), 4), "\n", sep = "")
    cat("Fit files all present: ", all(fit_manifest$fit_exists), "\n", sep = "")

    if (!is.null(fit_summary) && nrow(fit_summary) > 0L) {
        beta_cols <- intersect(c("beta0_mean", "beta1_mean", "beta2_mean", "r_mean", "lambda_rmse", "cor_log_lambda"), names(fit_summary))
        if (length(beta_cols) > 0L) {
            cat("\nSelected fit-summary columns:\n")
            print(fit_summary[, c("scenario_id", "rep_id", beta_cols), drop = FALSE])
        }
    }
    invisible(TRUE)
}

## ---- pre-flight checks ------------------------------------------------------
assert_file_exists_s4b_fit(file.path(root_dir, s3_core_file), "Scenario 3 core script")
assert_file_exists_s4b_fit(file.path(root_dir, s4b_data_file), "Scenario 4B data-generation script")

source(file.path(root_dir, s4b_data_file))
source_s4b_low_exposure(root = root_dir, verbose = FALSE)

## ---- S4B numerical guard for xi --------------------------------------------
## S4B is less sparse than S4A on average, but heterogeneous exposure can still
## produce low-information areas.  We use the same numerically safe xi computation
## as S4A so any instability is recorded instead of crashing through overflow.
s4b_log_xi_lower <- -40
s4b_log_xi_upper <-  40

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
        log_xi_j <- pmin(pmax(log_xi_j, s4b_log_xi_lower), s4b_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }

    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("S4B safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

cat(sprintf(
    "Using S4B safe compute_s3_xi(): log(xi) clamped to [%.1f, %.1f].\n",
    s4b_log_xi_lower, s4b_log_xi_upper
))

## ---- S4B stabilization guards ----------------------------------------------
.s4b_guard_env <- new.env(parent = emptyenv())
reset_s4b_numeric_guards <- function() {
    .s4b_guard_env$n_beta_guard <- 0L
    .s4b_guard_env$n_kappa_guard <- 0L
    .s4b_guard_env$n_lambda_input_guard <- 0L
    .s4b_guard_env$n_lambda_output_guard <- 0L
    invisible(TRUE)
}
get_s4b_numeric_guards <- function() {
    list(
        n_beta_guard = .s4b_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s4b_guard_env$n_kappa_guard %||% 0L,
        n_lambda_input_guard = .s4b_guard_env$n_lambda_input_guard %||% 0L,
        n_lambda_output_guard = .s4b_guard_env$n_lambda_output_guard %||% 0L
    )
}
reset_s4b_numeric_guards()

s4b_beta0_bounds <- c(-30, 10)
s4b_beta_bounds <- c(-5, 5)
s4b_kappa_bounds <- c(1e-10, 1e10)
s4b_lambda_bounds <- c(1e-10, 1e10)

.update_beta_s3_original <- update_beta
update_beta <- function(beta_current, ...) {
    res <- .update_beta_s3_original(beta_current = beta_current, ...)
    bad <- FALSE
    smp <- res$sample
    if (length(smp) < 3L || any(!is.finite(smp))) {
        bad <- TRUE
    } else {
        bad <- smp[1] < s4b_beta0_bounds[1] || smp[1] > s4b_beta0_bounds[2] ||
            any(smp[-1] < s4b_beta_bounds[1] | smp[-1] > s4b_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s4b_guard_env$n_beta_guard <- (.s4b_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s4b_guard_rejected <- TRUE
    } else {
        res$s4b_guard_rejected <- FALSE
    }
    res
}

.update_kappa_s3_original <- update_kappa
update_kappa <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
    y <- as.matrix(y_coarse)
    L <- as.matrix(lambda_tilde)
    X <- as.matrix(xi)
    if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
        stop("S4B safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
    }
    r_vec <- as.numeric(r)
    if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
    if (length(r_vec) != ncol(y)) {
        stop("S4B safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
    }
    R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
    shape <- y + R
    rate <- X * L + R
    guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
    if (any(guard)) {
        .s4b_guard_env$n_kappa_guard <- (.s4b_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s4b_guard_env$n_kappa_guard <- (.s4b_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s4b_kappa_bounds[1]), s4b_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s4b_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

.ffbs_lambda_all_s3_original <- ffbs_lambda_all
ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0, return_diag = TRUE) {
    X <- as.matrix(xi)
    K <- as.matrix(kappa)
    bad_in <- !is.finite(X) | X <= 0 | !is.finite(K) | K <= 0
    if (any(bad_in)) {
        .s4b_guard_env$n_lambda_input_guard <- (.s4b_guard_env$n_lambda_input_guard %||% 0L) + sum(bad_in)
    }
    X <- pmin(pmax(X, 1e-10), 1e10)
    K <- pmin(pmax(K, s4b_kappa_bounds[1]), s4b_kappa_bounds[2])
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
    bad_out <- !is.finite(L) | L <= 0 | L < s4b_lambda_bounds[1] | L > s4b_lambda_bounds[2]
    if (any(bad_out)) {
        .s4b_guard_env$n_lambda_output_guard <- (.s4b_guard_env$n_lambda_output_guard %||% 0L) + sum(bad_out)
        L <- pmin(pmax(L, s4b_lambda_bounds[1]), s4b_lambda_bounds[2])
        out$lambda_tilde <- L
        if (!is.null(out$diag)) {
            out$diag$min_lambda <- min(L)
            out$diag$max_lambda <- max(L)
            out$diag$s4b_lambda_output_guarded <- TRUE
        }
    }
    out
}

cat(sprintf(
    "Using S4B stabilization guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s4b_beta0_bounds[1], s4b_beta0_bounds[2],
    s4b_beta_bounds[1], s4b_beta_bounds[2],
    s4b_log_xi_lower, s4b_log_xi_upper
))

## ---- S4B fixed-gamma update override ---------------------------------------
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
cat(sprintf("Using S4B fixed-gamma override: gamma fixed at %.3f.\n", fixed_gamma_value))

source_data_manifest <- make_s4b_source_data_manifest(
    reps = reps_formal,
    root_arg = root_dir,
    data_dir_arg = data_dir,
    data_scenario_id_arg = data_scenario_id
)

cat("\n=== S4B source-data check ===\n")
print(source_data_manifest[, c(
    "rep_id", "mean_count", "median_count", "zero_prop", "total_count",
    "max_count", "exposure_mean_ratio", "exposure_cv",
    "lowest_exposure_group_mean_count", "highest_exposure_group_mean_count",
    "beta0_star_ident", "lambda_raw_min", "lambda_raw_max"
), drop = FALSE])

## ---- run Scenario 4B fixed-gamma fit ---------------------------------------
fit_summary <- fit_s4b_low_exposure_fixed_gamma_batch(
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
    sprintf("fit_S4B_low_exposure_fixed_gamma_rep%02d.rds", as.integer(reps_formal))
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
assert_true_s4b_fit(all(fit_manifest$fit_exists), "At least one S4B fit file was not created.")

ensure_dir_s4b_fit(fit_dir_full)
utils::write.csv(
    source_data_manifest,
    file.path(fit_dir_full, "s4b_source_data_manifest.csv"),
    row.names = FALSE
)
utils::write.csv(
    fit_manifest,
    file.path(fit_dir_full, "s4b_fit_manifest.csv"),
    row.names = FALSE
)
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
    file.path(fit_dir_full, "run_info_S4B_low_exposure_fixed_gamma_T100.rds")
)

## ---- final summary ----------------------------------------------------------
print_s4b_fit_driver_summary(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
)

cat("\n=== Main output locations ===\n")
cat("S4B data: ", file.path(data_dir, data_scenario_id), "\n", sep = "")
cat("Fits    : ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("\nScenario 4B low-heterogeneous-exposure fixed-gamma T100 fitting finished successfully.\n")

invisible(list(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
))
