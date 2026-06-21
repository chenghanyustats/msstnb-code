## ============================================================================
## s4d_short_time_series_continuous_x2.R
## Scenario 4D corrected short-time-series / small-T data-generation and calibration
## script for the MSSTNB project.
##
## Scenario 4D purpose:
##   Corrected short temporal horizon stress test corresponding to Scenario 3.
##
## Correction relative to the previous S4D script:
##   The previous S4D data-generation call did not explicitly pass x2_mode, so
##   R's match.arg() selected the first available option, binary_region. That made
##   x2 time-invariant within each region. This corrected script explicitly uses
##   x2_mode = "continuous_time", x2_ar = 0.50, and x2_innov_sd = 0.80.
##
## Design rule:
##   1. Keep the Scenario 3 data-generating mechanism: same mean structure,
##      exposure process, ICAR spatial effect, dynamic lambda process, r truth,
##      and gamma truth.
##   2. Use continuous-time x2 explicitly, matching the intended S3 T100 design.
##   3. Change only the time length T.
##   4. Generate complete short-T datasets directly from the Scenario 3 DGP.
##   5. Summarize whether the generated data remain high/moderate count and do
##      not accidentally become sparse-count, low-exposure, or overdispersion
##      stress regimes.
##
## This script only generates and calibrates data; it does not fit the model.
## ============================================================================

.same_dim_s4dctx <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.safe_ratio_s4dctx <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}

.cv_s4dctx <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) return(NA_real_)
    stats::sd(x) / mx
}

.bind_rows_aligned_s4dctx <- function(x) {
    if (length(x) == 0L) return(data.frame())
    all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    out <- lapply(x, function(df) {
        missing <- setdiff(all_names, names(df))
        for (nm in missing) df[[nm]] <- NA
        df[, all_names, drop = FALSE]
    })
    do.call(rbind, out)
}

.require_file_s4dctx <- function(path) {
    if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
    invisible(path)
}

.source_checked_s4dctx <- function(path, verbose = TRUE) {
    .require_file_s4dctx(path)
    if (isTRUE(verbose)) message("source: ", path)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4dctx <- function(root = ".", s3_script = NULL) {
    if (!is.null(s3_script)) return(s3_script)

    candidates <- c(
        file.path(root, "s3_dynamic_learned_gamma.R"),
        file.path(root, "scripts", "s3_dynamic_learned_gamma.R"),
        file.path(root, "R", "s3_dynamic_learned_gamma.R"),
        file.path(root, "R", "scenarios", "s3_dynamic_learned_gamma.R"),
        file.path(root, "scenarios", "s3_dynamic_learned_gamma.R")
    )
    hits <- candidates[file.exists(candidates)]
    if (length(hits) == 0L) {
        stop(
            "Could not find s3_dynamic_learned_gamma.R. ",
            "Pass its location with s3_script = 'path/to/s3_dynamic_learned_gamma.R'.",
            call. = FALSE
        )
    }
    hits[1L]
}

.require_object_s4dctx <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source Scenario 3 and project dependencies
## -----------------------------------------------------------------------------
source_s4dctx_short_time_continuous_x2_series <- function(root = ".",
                                         s3_script = NULL,
                                         verbose = TRUE) {
    s3_path <- .find_s3_script_s4dctx(root = root, s3_script = s3_script)
    .source_checked_s4dctx(s3_path, verbose = verbose)

    .require_object_s4dctx("source_s3_dynamic_learned_gamma")
    source_s3_dynamic_learned_gamma(root = root, verbose = verbose)

    needed <- c(
        "simulate_s3_dynamic_learned_gamma_one",
        "validate_s3_data",
        "recenter",
        "REP_SEEDS",
        "TT", "N1", "N_CHILDREN"
    )
    missing <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop("After sourcing Scenario 3 dependencies, missing objects: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("Scenario 4D corrected short-time-series continuous-x2 data generator loaded.")
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Summary helpers
## -----------------------------------------------------------------------------
.count_stats_s4dctx <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(
        yy,
        probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
        names = FALSE,
        type = 7
    )
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        nonzero_prop = mean(yy > 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = .cv_s4dctx(yy),
        variance_to_mean = .safe_ratio_s4dctx(vv, mn),
        max_to_mean = .safe_ratio_s4dctx(max(yy), mn),
        q05_count = as.numeric(qs[1]),
        q25_count = as.numeric(qs[2]),
        q50_count = as.numeric(qs[3]),
        q75_count = as.numeric(qs[4]),
        q90_count = as.numeric(qs[5]),
        q95_count = as.numeric(qs[6]),
        q99_count = as.numeric(qs[7]),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

.lambda_stats_s4dctx <- function(lambda_tilde, prefix = "") {
    lam <- as.matrix(lambda_tilde)
    log_lam <- log(pmax(lam, .Machine$double.xmin))
    delta_log <- if (nrow(log_lam) >= 2L) diff(log_lam) else matrix(NA_real_, 0L, ncol(log_lam))
    qs <- stats::quantile(as.numeric(lam), probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    log_qs <- stats::quantile(as.numeric(log_lam), probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                              names = FALSE, type = 7)
    out <- data.frame(
        lambda_raw_min = min(lam),
        lambda_raw_q01 = as.numeric(qs[1]),
        lambda_raw_q05 = as.numeric(qs[2]),
        lambda_raw_median = as.numeric(qs[3]),
        lambda_raw_q95 = as.numeric(qs[4]),
        lambda_raw_q99 = as.numeric(qs[5]),
        lambda_raw_max = max(lam),
        log_lambda_mean = mean(log_lam),
        log_lambda_sd = stats::sd(as.numeric(log_lam)),
        log_lambda_q01 = as.numeric(log_qs[1]),
        log_lambda_q05 = as.numeric(log_qs[2]),
        log_lambda_median = as.numeric(log_qs[3]),
        log_lambda_q95 = as.numeric(log_qs[4]),
        log_lambda_q99 = as.numeric(log_qs[5]),
        delta_log_lambda_mean = if (length(delta_log) > 0L) mean(delta_log) else NA_real_,
        delta_log_lambda_sd = if (length(delta_log) > 0L) stats::sd(as.numeric(delta_log)) else NA_real_,
        delta_log_lambda_abs_mean = if (length(delta_log) > 0L) mean(abs(delta_log)) else NA_real_,
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

.kappa_stats_s4dctx <- function(kappa, prefix = "") {
    kk <- as.numeric(kappa)
    out <- data.frame(
        kappa_mean = mean(kk),
        kappa_sd = stats::sd(kk),
        kappa_cv = .cv_s4dctx(kk),
        kappa_min = min(kk),
        kappa_median = stats::median(kk),
        kappa_max = max(kk),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}


.safe_cor_s4dctx <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (stats::sd(x[ok]) <= .Machine$double.eps ||
        stats::sd(y[ok]) <= .Machine$double.eps) return(NA_real_)
    as.numeric(stats::cor(x[ok], y[ok]))
}

.mean_lag1_by_region_s4dctx <- function(x) {
    xx <- as.matrix(x)
    if (nrow(xx) < 3L) return(NA_real_)
    cors <- apply(xx, 2, function(z) .safe_cor_s4dctx(z[-length(z)], z[-1L]))
    mean(cors, na.rm = TRUE)
}

.within_sd_by_region_s4dctx <- function(x) {
    xx <- as.matrix(x)
    apply(xx, 2, stats::sd)
}

.covariate_geometry_stats_s4dctx <- function(dat, prefix = "") {
    x1 <- as.matrix(dat$x1)
    x2 <- as.matrix(dat$x2)

    x1_within_sd <- .within_sd_by_region_s4dctx(x1)
    x2_within_sd <- .within_sd_by_region_s4dctx(x2)
    x1_bar <- colMeans(x1)
    x2_bar <- colMeans(x2)

    out <- data.frame(
        x2_mode = if (!is.null(dat$x2_mode)) as.character(dat$x2_mode) else NA_character_,
        x2_ar = if (!is.null(dat$x2_ar)) as.numeric(dat$x2_ar) else NA_real_,
        x2_innov_sd = if (!is.null(dat$x2_innov_sd)) as.numeric(dat$x2_innov_sd) else NA_real_,
        mean_x1_within_sd = mean(x1_within_sd, na.rm = TRUE),
        mean_x2_within_sd = mean(x2_within_sd, na.rm = TRUE),
        min_x1_within_sd = min(x1_within_sd, na.rm = TRUE),
        min_x2_within_sd = min(x2_within_sd, na.rm = TRUE),
        max_x1_within_sd = max(x1_within_sd, na.rm = TRUE),
        max_x2_within_sd = max(x2_within_sd, na.rm = TRUE),
        mean_lag1_x1 = .mean_lag1_by_region_s4dctx(x1),
        mean_lag1_x2 = .mean_lag1_by_region_s4dctx(x2),
        between_var_x1 = stats::var(x1_bar),
        between_var_x2 = stats::var(x2_bar),
        within_var_x1 = mean(apply(x1, 2, stats::var), na.rm = TRUE),
        within_var_x2 = mean(apply(x2, 2, stats::var), na.rm = TRUE),
        cor_vec_x1_x2 = .safe_cor_s4dctx(as.vector(x1), as.vector(x2)),
        cor_bar_x1_x2 = .safe_cor_s4dctx(x1_bar, x2_bar),
        stringsAsFactors = FALSE
    )
    out$between_within_ratio_x1 <- .safe_ratio_s4dctx(out$between_var_x1, out$within_var_x1)
    out$between_within_ratio_x2 <- .safe_ratio_s4dctx(out$between_var_x2, out$within_var_x2)
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

## -----------------------------------------------------------------------------
## Validation
## -----------------------------------------------------------------------------
validate_s4dctx_short_time_continuous_x2_data <- function(dat) {
    .require_object_s4dctx("validate_s3_data")
    validate_s3_data(dat)

    required <- c(
        "scenario_id", "stress_type", "stress_description",
        "TT_reference", "TT_short", "TT_ratio_to_reference", "n_cells",
        "r_star", "gamma_star", "lambda_tilde", "phi_star",
        "x2_mode", "x2_ar", "x2_innov_sd"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("S4D short-time dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    if (!identical(dat$stress_type, "short_time_series")) {
        stop("dat$stress_type must be 'short_time_series'.", call. = FALSE)
    }
    if (!identical(as.integer(dat$TT_short), as.integer(dat$TT))) {
        stop("dat$TT_short must equal dat$TT.", call. = FALSE)
    }
    if (!identical(as.integer(dat$n_cells), as.integer(dat$TT) * as.integer(dat$n1))) {
        stop("dat$n_cells is inconsistent with dat$TT and dat$n1.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0)) {
        stop("dat$lambda_tilde must be positive and finite.", call. = FALSE)
    }
    if (!identical(as.character(dat$x2_mode), "continuous_time")) {
        stop("Corrected S4D requires dat$x2_mode == 'continuous_time'.", call. = FALSE)
    }
    x2_within_sd <- apply(as.matrix(dat$x2), 2, stats::sd)
    if (any(!is.finite(x2_within_sd)) || min(x2_within_sd) <= 1e-8) {
        stop(
            "Corrected S4D continuous-time x2 audit failed: at least one region has near-zero x2 temporal SD.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Scenario 4D short-T simulation
## -----------------------------------------------------------------------------
simulate_s4dctx_short_time_continuous_x2_one <- function(seed = 1L,
                                        TT_short = 30,
                                        TT_reference = 100,
                                        n1_use = N1,
                                        n_children_use = N_CHILDREN,
                                        beta0_truth = -1.5,
                                        beta_truth = c(0.5, -0.4),
                                        r_truth = 15,
                                        tau_phi_truth = 2,
                                        gamma_truth = 0.8,
                                        delta_truth = 0.9,
                                        x2_mode = "continuous_time",
                                        x2_ar = 0.50,
                                        x2_innov_sd = 0.80,
                                        scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                        rep_id = NA_integer_,
                                        max_poisson_rate = 1e7,
                                        ...) {
    .require_object_s4dctx("simulate_s3_dynamic_learned_gamma_one")
    .require_object_s4dctx("recenter")

    TT_short <- as.integer(TT_short)
    TT_reference <- as.integer(TT_reference)
    if (!is.finite(TT_short) || TT_short < 2L) {
        stop("TT_short must be an integer >= 2.", call. = FALSE)
    }
    if (!is.finite(TT_reference) || TT_reference <= 0L) {
        stop("TT_reference must be positive.", call. = FALSE)
    }
    if (TT_short >= TT_reference) {
        warning("TT_short is not smaller than TT_reference; S4D is intended for short T.", call. = FALSE)
    }

    x2_mode <- match.arg(x2_mode, choices = "continuous_time")

    dat <- simulate_s3_dynamic_learned_gamma_one(
        seed = seed,
        TT_use = TT_short,
        n1_use = n1_use,
        n_children_use = n_children_use,
        beta0_truth = beta0_truth,
        beta_truth = beta_truth,
        r_truth = r_truth,
        tau_phi_truth = tau_phi_truth,
        gamma_truth = gamma_truth,
        delta_truth = delta_truth,
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        scenario_id = scenario_id,
        rep_id = rep_id,
        max_poisson_rate = max_poisson_rate,
        ...
    )

    rc_truth <- recenter(
        beta0 = beta0_truth,
        phi = dat$phi_star,
        lambda_tilde = dat$lambda_tilde,
        return_diag = TRUE
    )

    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_short_time_series_continuous_x2"
    dat$stress_type <- "short_time_series"
    dat$stress_description <- paste0(
        "Scenario 3 data-generating mechanism with shorter time length T=", TT_short,
        ". This corrected version explicitly uses continuous-time x2, isolating ",
        "short temporal horizon stress from sparse-count, low-exposure, ",
        "strong-overdispersion, and unintended region-level-x2 stresses."
    )
    dat$TT_reference <- TT_reference
    dat$TT_short <- TT_short
    dat$TT_ratio_to_reference <- TT_short / TT_reference
    dat$n_cells <- as.integer(TT_short) * as.integer(dat$n1)
    dat$r_truth <- r_truth
    dat$gamma_truth <- gamma_truth
    dat$delta_truth <- delta_truth
    dat$x2_mode <- x2_mode
    dat$x2_ar <- x2_ar
    dat$x2_innov_sd <- x2_innov_sd
    dat$beta0_star_ident <- rc_truth$beta0
    dat$beta_star_ident <- dat$beta_star
    dat$phi_star_ident <- rc_truth$phi
    dat$lambda_tilde_ident <- rc_truth$lambda_tilde
    dat$lambda_recenter_diag <- rc_truth$diag

    y <- dat$y_coarse
    dat$mean_count <- mean(y)
    dat$median_count <- stats::median(as.numeric(y))
    dat$zero_prop <- mean(y == 0)
    dat$nonzero_prop <- mean(y > 0)
    dat$total_count <- sum(y)
    dat$max_count <- max(y)
    dat$count_sd <- stats::sd(as.numeric(y))
    dat$count_cv <- .cv_s4dctx(y)
    dat$variance_to_mean <- .safe_ratio_s4dctx(stats::var(as.numeric(y)), mean(y))
    dat$count_quantiles <- stats::quantile(
        as.numeric(y),
        probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1),
        names = TRUE
    )

    validate_s4dctx_short_time_continuous_x2_data(dat)
    dat
}

summarise_s4dctx_short_time_continuous_x2_one <- function(dat) {
    validate_s4dctx_short_time_continuous_x2_data(dat)

    y_stats <- .count_stats_s4dctx(dat$y_coarse, prefix = "")
    lam_stats <- .lambda_stats_s4dctx(dat$lambda_tilde, prefix = "")
    lam_ident_stats <- .lambda_stats_s4dctx(dat$lambda_tilde_ident, prefix = "lambda_ident_")
    k_stats <- .kappa_stats_s4dctx(dat$kappa, prefix = "")
    geom_stats <- .covariate_geometry_stats_s4dctx(dat, prefix = "")

    data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        TT_reference = as.integer(dat$TT_reference),
        TT_short = as.integer(dat$TT_short),
        TT_ratio_to_reference = dat$TT_ratio_to_reference,
        n1 = as.integer(dat$n1),
        n_cells = as.integer(dat$n_cells),
        stress_type = dat$stress_type,
        beta0_star = dat$beta0_star,
        beta0_star_ident = dat$beta0_star_ident,
        beta1_star = dat$beta_star[1],
        beta2_star = dat$beta_star[2],
        r_truth_mean = mean(dat$r_star),
        r_truth_min = min(dat$r_star),
        r_truth_max = max(dat$r_star),
        gamma_truth_mean = mean(dat$gamma_star),
        gamma_truth_min = min(dat$gamma_star),
        gamma_truth_max = max(dat$gamma_star),
        delta_truth = dat$delta_star,
        mean_exposure = mean(dat$e),
        min_exposure = min(dat$e),
        max_exposure = max(dat$e),
        exposure_cv = .cv_s4dctx(dat$e),
        x1_sd = stats::sd(as.numeric(dat$x1)),
        x2_sd = stats::sd(as.numeric(dat$x2)),
        geom_stats,
        phi_truth_sd = stats::sd(as.numeric(dat$phi_star)),
        phi_ident_sd = stats::sd(as.numeric(dat$phi_star_ident)),
        coherent = isTRUE(dat$coherent),
        y_stats,
        k_stats,
        lam_stats,
        lam_ident_stats,
        stringsAsFactors = FALSE
    )
}

summarise_s4dctx_short_time_continuous_x2_from_files <- function(files) {
    if (length(files) == 0L) stop("No files supplied.", call. = FALSE)
    out <- lapply(files, function(ff) {
        dat <- readRDS(ff)
        ss <- summarise_s4dctx_short_time_continuous_x2_one(dat)
        ss$file <- ff
        ss
    })
    .bind_rows_aligned_s4dctx(out)
}

summarise_s4dctx_short_time_continuous_x2_from_dir <- function(root = ".",
                                              data_dir = "data_s4d_short_time_continuous_x2",
                                              scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2") {
    in_dir <- file.path(root, data_dir, scenario_id)
    files <- list.files(in_dir, pattern = "^data_rep[0-9]+\\.rds$", full.names = TRUE)
    if (length(files) == 0L) stop("No S4D data files found in: ", in_dir, call. = FALSE)
    summarise_s4dctx_short_time_continuous_x2_from_files(files)
}

simulate_s4dctx_short_time_continuous_x2_batch <- function(reps = 1:10,
                                          data_dir = "data_s4d_short_time_continuous_x2",
                                          scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
                                          root = ".",
                                          seed_base = NULL,
                                          overwrite_existing = TRUE,
                                          verbose = TRUE,
                                          manifest_name = NULL,
                                          ...) {
    out_dir <- file.path(root, data_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    manifest <- list()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        out_file <- file.path(out_dir, paste0("data_rep", rr, ".rds"))

        if (file.exists(out_file) && !isTRUE(overwrite_existing)) {
            if (isTRUE(verbose)) message("Skipping existing file: ", out_file)
            dat <- readRDS(out_file)
        } else {
            seed <- if (!is.null(seed_base)) {
                as.integer(seed_base + rep_id)
            } else if (exists("REP_SEEDS", envir = .GlobalEnv) &&
                       rep_id <= length(REP_SEEDS)) {
                as.integer(REP_SEEDS[rep_id])
            } else {
                as.integer(2026000L + rep_id)
            }

            dat <- simulate_s4dctx_short_time_continuous_x2_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                message(sprintf(
                    paste0(
                        "Saved %s | T=%d n_cells=%d mean_count=%.2f zero_prop=%.3f ",
                        "cv=%.2f v/m=%.2f max=%s log_lambda_sd=%.3f x2_lag1=%.3f"
                    ),
                    out_file,
                    as.integer(dat$TT),
                    as.integer(dat$n_cells),
                    dat$mean_count,
                    dat$zero_prop,
                    dat$count_cv,
                    dat$variance_to_mean,
                    as.character(dat$max_count),
                    stats::sd(as.numeric(log(dat$lambda_tilde))),
                    .mean_lag1_by_region_s4dctx(dat$x2)
                ))
            }
        }

        row <- summarise_s4dctx_short_time_continuous_x2_one(dat)
        row$file <- out_file
        manifest[[rr]] <- row
    }

    manifest_df <- .bind_rows_aligned_s4dctx(manifest)
    if (is.null(manifest_name)) manifest_name <- paste0("manifest_", scenario_id, ".csv")
    manifest_file <- file.path(out_dir, manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message(sprintf(
            paste0(
                "S4D short-T summary | reps=%d T=%s mean_count=%.2f zero_prop=%.3f ",
                "cv=%.2f v/m=%.2f total_count_sum=%s"
            ),
            nrow(manifest_df),
            paste(unique(manifest_df$TT), collapse = ","),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            mean(manifest_df$count_cv),
            mean(manifest_df$variance_to_mean),
            as.character(sum(manifest_df$total_count))
        ))
    }
    invisible(manifest_df)
}

run_s4dctx_short_time_continuous_x2_data_generation <- function(root = ".",
                                               s3_script = NULL,
                                               reps = 1:10,
                                               TT_short = 30,
                                               TT_reference = 100,
                                               data_dir = "data_s4d_short_time_continuous_x2",
                                               scenario_id = NULL,
                                               seed_base = NULL,
                                               overwrite_existing = TRUE,
                                               verbose = TRUE,
                                               ...) {
    source_s4dctx_short_time_continuous_x2_series(root = root, s3_script = s3_script, verbose = verbose)
    if (is.null(scenario_id)) scenario_id <- paste0("S4D_SHORT_TIME_T", as.integer(TT_short), "_CONTINUOUS_X2")

    simulate_s4dctx_short_time_continuous_x2_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_short = TT_short,
        TT_reference = TT_reference,
        ...
    )
}

## -----------------------------------------------------------------------------
## Data-quality checks and calibration helpers
## -----------------------------------------------------------------------------
check_s4dctx_short_time_continuous_x2_data_summary <- function(manifest,
                                              target_mean_count_range = c(80, 350),
                                              target_zero_prop_max = 0.20,
                                              target_count_cv_max = 2.50,
                                              target_variance_to_mean_max = Inf,
                                              target_abs_beta0_ident_max = 20,
                                              target_lambda_raw_max = 10,
                                              minimum_n_cells = 100) {
    if (is.character(manifest)) manifest <- read.csv(manifest)

    required <- c(
        "rep_id", "TT", "TT_reference", "TT_ratio_to_reference", "n1", "n_cells",
        "mean_count", "median_count", "zero_prop", "total_count", "max_count",
        "count_sd", "count_cv", "variance_to_mean", "q95_count", "q99_count",
        "beta0_star_ident", "lambda_raw_median", "lambda_raw_max",
        "log_lambda_sd", "delta_log_lambda_sd", "delta_log_lambda_abs_mean",
        "r_truth_mean", "gamma_truth_mean",
        "x2_mode", "mean_x2_within_sd", "min_x2_within_sd", "mean_lag1_x2"
    )
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0L) {
        stop("Manifest is missing required columns: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    out <- data.frame(
        n_reps = nrow(manifest),
        TT_unique = paste(sort(unique(manifest$TT)), collapse = ","),
        TT_reference_unique = paste(sort(unique(manifest$TT_reference)), collapse = ","),
        TT_ratio_to_reference_avg = mean(manifest$TT_ratio_to_reference),
        n1_unique = paste(sort(unique(manifest$n1)), collapse = ","),
        n_cells_avg = mean(manifest$n_cells),
        n_cells_min = min(manifest$n_cells),
        x2_mode_unique = paste(sort(unique(as.character(manifest$x2_mode))), collapse = ","),
        mean_x2_within_sd_avg = mean(manifest$mean_x2_within_sd, na.rm = TRUE),
        min_x2_within_sd_min = min(manifest$min_x2_within_sd, na.rm = TRUE),
        mean_lag1_x2_avg = mean(manifest$mean_lag1_x2, na.rm = TRUE),
        mean_lag1_x2_min = min(manifest$mean_lag1_x2, na.rm = TRUE),
        mean_lag1_x2_max = max(manifest$mean_lag1_x2, na.rm = TRUE),
        mean_count_avg = mean(manifest$mean_count),
        mean_count_min = min(manifest$mean_count),
        mean_count_max = max(manifest$mean_count),
        median_count_avg = mean(manifest$median_count),
        zero_prop_avg = mean(manifest$zero_prop),
        zero_prop_min = min(manifest$zero_prop),
        zero_prop_max = max(manifest$zero_prop),
        total_count_sum = sum(manifest$total_count),
        total_count_avg = mean(manifest$total_count),
        max_count_avg = mean(manifest$max_count),
        max_count_max = max(manifest$max_count),
        count_sd_avg = mean(manifest$count_sd),
        count_cv_avg = mean(manifest$count_cv),
        count_cv_min = min(manifest$count_cv),
        count_cv_max = max(manifest$count_cv),
        variance_to_mean_avg = mean(manifest$variance_to_mean),
        variance_to_mean_min = min(manifest$variance_to_mean),
        variance_to_mean_max = max(manifest$variance_to_mean),
        q95_count_avg = mean(manifest$q95_count),
        q99_count_avg = mean(manifest$q99_count),
        q99_count_max = max(manifest$q99_count),
        kappa_cv_avg = mean(manifest$kappa_cv),
        beta0_ident_min = min(manifest$beta0_star_ident),
        beta0_ident_median = stats::median(manifest$beta0_star_ident),
        beta0_ident_max = max(manifest$beta0_star_ident),
        max_abs_beta0_ident = max(abs(manifest$beta0_star_ident)),
        lambda_raw_median_avg = mean(manifest$lambda_raw_median),
        lambda_raw_max_max = max(manifest$lambda_raw_max),
        log_lambda_sd_avg = mean(manifest$log_lambda_sd),
        log_lambda_sd_min = min(manifest$log_lambda_sd),
        log_lambda_sd_max = max(manifest$log_lambda_sd),
        delta_log_lambda_sd_avg = mean(manifest$delta_log_lambda_sd, na.rm = TRUE),
        delta_log_lambda_abs_mean_avg = mean(manifest$delta_log_lambda_abs_mean, na.rm = TRUE),
        r_truth_mean_avg = mean(manifest$r_truth_mean),
        gamma_truth_mean_avg = mean(manifest$gamma_truth_mean),
        stringsAsFactors = FALSE
    )

    out$passes_mean_count_target <-
        out$mean_count_avg >= target_mean_count_range[1] &&
        out$mean_count_avg <= target_mean_count_range[2]
    out$passes_zero_prop_target <- out$zero_prop_avg <= target_zero_prop_max
    out$passes_count_cv_target <- out$count_cv_avg <= target_count_cv_max
    out$passes_variance_to_mean_target <- out$variance_to_mean_avg <= target_variance_to_mean_max
    out$passes_identified_scale_target <- out$max_abs_beta0_ident <= target_abs_beta0_ident_max
    out$passes_lambda_raw_target <- out$lambda_raw_max_max <= target_lambda_raw_max
    out$passes_n_cells_target <- out$n_cells_min >= minimum_n_cells
    out$passes_continuous_x2_target <-
        identical(out$x2_mode_unique, "continuous_time") &&
        is.finite(out$min_x2_within_sd_min) && out$min_x2_within_sd_min > 1e-8 &&
        is.finite(out$mean_lag1_x2_avg)

    out$passes_s4d_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_count_cv_target) &&
        isTRUE(out$passes_variance_to_mean_target) &&
        isTRUE(out$passes_identified_scale_target) &&
        isTRUE(out$passes_lambda_raw_target) &&
        isTRUE(out$passes_n_cells_target) &&
        isTRUE(out$passes_continuous_x2_target)

    out
}

calibrate_s4dctx_short_T_grid <- function(root = ".",
                                       s3_script = NULL,
                                       reps = 1:10,
                                       TT_grid = c(80, 50, 30, 20),
                                       TT_reference = 100,
                                       data_dir = "data_s4d_short_time_continuous_x2_calibration",
                                       base_scenario_id = "S4D_SHORT_TIME_CONTINUOUS_X2_CALIB",
                                       overwrite_existing = TRUE,
                                       verbose = TRUE,
                                       ...) {
    source_s4dctx_short_time_continuous_x2_series(root = root, s3_script = s3_script, verbose = verbose)

    out <- list()
    for (TT_now in TT_grid) {
        scenario_id <- paste0(base_scenario_id, "_T", as.integer(TT_now))
        if (isTRUE(verbose)) {
            message("\n--- S4D calibration: TT_short = ", TT_now,
                    " | scenario_id = ", scenario_id, " ---")
        }

        manifest <- simulate_s4dctx_short_time_continuous_x2_batch(
            reps = reps,
            data_dir = data_dir,
            scenario_id = scenario_id,
            root = root,
            overwrite_existing = overwrite_existing,
            verbose = verbose,
            TT_short = as.integer(TT_now),
            TT_reference = TT_reference,
            ...
        )
        chk <- check_s4dctx_short_time_continuous_x2_data_summary(manifest)
        chk$scenario_id <- scenario_id
        chk$TT_short <- as.integer(TT_now)
        out[[paste0("T", as.integer(TT_now))]] <- chk
    }

    calib <- .bind_rows_aligned_s4dctx(out)
    out_dir <- file.path(root, data_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    calib_file <- file.path(out_dir, paste0("calibration_summary_", base_scenario_id, ".csv"))
    write.csv(calib, calib_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved S4D calibration summary: ", calib_file)

    invisible(calib)
}

## Example usage:
##
## source("s4d_short_time_series_continuous_x2.R")
## calib_s4d <- calibrate_s4d_short_T_continuous_x2_grid(
##     root = ".",
##     reps = 1:10,
##     TT_grid = c(80, 50, 30, 20),
##     TT_reference = 100,
##     overwrite_existing = TRUE
## )
## calib_s4d[, c(
##     "TT_short", "n_cells_avg", "mean_count_avg", "zero_prop_avg",
##     "count_cv_avg", "variance_to_mean_avg", "total_count_sum",
##     "x2_mode_unique", "mean_x2_within_sd_avg", "mean_lag1_x2_avg",
##     "max_abs_beta0_ident", "lambda_raw_max_max", "log_lambda_sd_avg",
##     "delta_log_lambda_sd_avg", "passes_continuous_x2_target",
##     "passes_s4d_data_check"
## )]
##
## After choosing the official T, for example T = 30:
## manifest_s4d <- run_s4d_short_time_continuous_x2_data_generation(
##     root = ".",
##     reps = 1:10,
##     TT_short = 30,
##     TT_reference = 100,
##     scenario_id = "S4D_SHORT_TIME_T30_CONTINUOUS_X2",
##     overwrite_existing = TRUE
## )
## check_s4d_short_time_continuous_x2_data_summary(manifest_s4d)


## -----------------------------------------------------------------------------
## Quick audit helper for corrected S4D files
## -----------------------------------------------------------------------------
audit_s4d_short_time_continuous_x2_mode <- function(data_dir,
                                                    reps = 1:10,
                                                    root = ".") {
    out <- lapply(reps, function(r) {
        ff <- file.path(root, data_dir, sprintf("data_rep%02d.rds", as.integer(r)))
        dat <- readRDS(ff)
        x2 <- as.matrix(dat$x2)
        x2_within_sd <- apply(x2, 2, stats::sd)
        data.frame(
            rep_id = as.integer(r),
            x2_mode = if (!is.null(dat$x2_mode)) as.character(dat$x2_mode) else NA_character_,
            dim_x2 = paste(dim(x2), collapse = " x "),
            mean_x2_within_sd = mean(x2_within_sd, na.rm = TRUE),
            min_x2_within_sd = min(x2_within_sd, na.rm = TRUE),
            max_x2_within_sd = max(x2_within_sd, na.rm = TRUE),
            mean_lag1_x2 = .mean_lag1_by_region_s4dctx(x2),
            stringsAsFactors = FALSE
        )
    })
    .bind_rows_aligned_s4dctx(out)
}


## -----------------------------------------------------------------------------
## Clean public aliases
## -----------------------------------------------------------------------------
source_s4d_short_time_continuous_x2 <- source_s4dctx_short_time_continuous_x2_series
simulate_s4d_short_time_continuous_x2_one <- simulate_s4dctx_short_time_continuous_x2_one
simulate_s4d_short_time_continuous_x2_batch <- simulate_s4dctx_short_time_continuous_x2_batch
run_s4d_short_time_continuous_x2_data_generation <- run_s4dctx_short_time_continuous_x2_data_generation
check_s4d_short_time_continuous_x2_data_summary <- check_s4dctx_short_time_continuous_x2_data_summary
calibrate_s4d_short_T_continuous_x2_grid <- calibrate_s4dctx_short_T_grid
calibrate_s4d_short_time_continuous_x2_grid <- calibrate_s4dctx_short_T_grid
