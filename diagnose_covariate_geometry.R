## ============================================================================
## diagnose_covariate_geometry.R
##
## Purpose:
##   Diagnose whether beta2 fragility is driven by covariate geometry rather
##   than by coefficient-index asymmetry.  The script reads existing MSSTNB
##   scenario data files, computes temporal structure, spatial alignment, and
##   covariate-collinearity diagnostics for x1 and x2, and optionally merges
##   replicate-level posterior summaries to relate beta bias to geometry.
##
## Main outputs:
##   analysis_covariate_geometry/tables/
##     scenario_file_manifest.csv
##     covariate_geometry_by_rep.csv
##     covariate_geometry_by_scenario.csv
##     covariate_geometry_long.csv
##     beta_geometry_by_rep.csv                 [if fit summaries exist]
##     beta_geometry_associations.csv           [if fit summaries exist]
##     covariance_diagnostic_interpretation.csv
##
## Intended usage:
##   source("diagnose_covariate_geometry.R")
##
## The object `covariate_geometry_diagnostic` will be created in the workspace.
## ============================================================================

## ---- small utilities --------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

.is_finite_nonzero_sd <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    length(x) >= 2L && is.finite(stats::sd(x)) && stats::sd(x) > 0
}

.safe_num <- function(x) {
    suppressWarnings(as.numeric(x))
}

.safe_mean <- function(x) {
    x <- .safe_num(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
}

.safe_sd <- function(x) {
    x <- .safe_num(x)
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x)
}

.safe_min <- function(x) {
    x <- .safe_num(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    min(x)
}

.safe_max <- function(x) {
    x <- .safe_num(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(x)
}

.safe_cor <- function(x, y) {
    x <- .safe_num(x)
    y <- .safe_num(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    x <- x[ok]
    y <- y[ok]
    if (stats::sd(x) <= 0 || stats::sd(y) <= 0) return(NA_real_)
    as.numeric(stats::cor(x, y))
}

.cv <- function(x) {
    x <- .safe_num(x)
    mu <- mean(x, na.rm = TRUE)
    sig <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(mu) || abs(mu) < .Machine$double.eps) return(NA_real_)
    sig / abs(mu)
}

.bind_rows_aligned <- function(x) {
    x <- x[!vapply(x, is.null, logical(1))]
    if (length(x) == 0L) return(data.frame())
    all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    aligned <- lapply(x, function(df) {
        missing <- setdiff(all_names, names(df))
        for (nm in missing) df[[nm]] <- NA
        df[, all_names, drop = FALSE]
    })
    do.call(rbind, aligned)
}

.parse_rep_id_from_file <- function(path) {
    b <- basename(path)
    out <- sub("^data_rep([0-9]+)\\.rds$", "\\1", b)
    suppressWarnings(as.integer(out))
}

.ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

.first_existing_dir <- function(root, data_dir, scenario_candidates) {
    for (sid in scenario_candidates) {
        p <- file.path(root, data_dir, sid)
        if (dir.exists(p)) return(list(path = p, scenario_id = sid))
    }
    list(path = file.path(root, data_dir, scenario_candidates[1]), scenario_id = scenario_candidates[1])
}

## ---- default scenario registry ---------------------------------------------
## These are the official data/fitting locations used in the current S3/S4 work.
## Missing folders are skipped automatically.
default_covariate_geometry_scenarios <- function(root = ".") {
    make_entry <- function(label, scenario_family, data_dir, data_scenario_candidates,
                           fit_summary_file = NA_character_, subset_label = NA_character_) {
        dd <- .first_existing_dir(root, data_dir, data_scenario_candidates)
        data.frame(
            scenario_label = label,
            scenario_family = scenario_family,
            data_dir = data_dir,
            data_scenario_id = dd$scenario_id,
            data_path = dd$path,
            fit_summary_file = fit_summary_file,
            subset_label = subset_label,
            stringsAsFactors = FALSE
        )
    }

    .bind_rows_aligned(list(
        make_entry(
            label = "S3 control T100",
            scenario_family = "S3",
            data_dir = "data_revised",
            data_scenario_candidates = c("DGP_DYNAMIC_T100", "S2_DYNAMIC_FIXED_GAMMA_T100"),
            fit_summary_file = file.path(
                root,
                "output_s3_control_fixed_gamma_stabilized",
                "S3_CONTROL_FIXED_GAMMA_STABILIZED_T100",
                "summary_S3_control_fixed_gamma_stabilized_all_reps.csv"
            ),
            subset_label = "S3 fixed-gamma control"
        ),
        make_entry(
            label = "S4A sparse counts",
            scenario_family = "S4A",
            data_dir = "data_s4a_sparse_counts",
            data_scenario_candidates = c("S4A_SPARSE_COUNTS_OBS_T100", "S4A_SPARSE_COUNTS_T100"),
            fit_summary_file = file.path(
                root,
                "output_s4a_sparse_counts_fixed_gamma",
                "S4A_SPARSE_COUNTS_OBS_FIXED_GAMMA_T100",
                "summary_S4A_sparse_counts_obs_fixed_gamma_all_reps.csv"
            ),
            subset_label = "S4A fixed-gamma"
        ),
        make_entry(
            label = "S4B low exposure",
            scenario_family = "S4B",
            data_dir = "data_s4b_low_exposure",
            data_scenario_candidates = c("S4B_LOW_HETEROGENEOUS_EXPOSURE_T100"),
            fit_summary_file = file.path(
                root,
                "output_s4b_low_exposure_fixed_gamma",
                "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_T100",
                "summary_S4B_low_exposure_fixed_gamma_all_reps.csv"
            ),
            subset_label = "S4B fixed-gamma"
        ),
        make_entry(
            label = "S4C overdispersion",
            scenario_family = "S4C",
            data_dir = "data_s4c_overdispersion",
            data_scenario_candidates = c("S4C_STRONG_OVERDISPERSION_T100"),
            fit_summary_file = file.path(
                root,
                "output_s4c_small_r_fixed_gamma",
                "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_T100",
                "summary_S4C_small_r_fixed_gamma_all_reps.csv"
            ),
            subset_label = "S4C fixed-gamma"
        ),
        make_entry(
            label = "S4D short time",
            scenario_family = "S4D",
            data_dir = "data_s4d_short_time",
            data_scenario_candidates = c("S4D_SHORT_TIME_T30"),
            fit_summary_file = file.path(
                root,
                "output_s4d_short_time_fixed_gamma",
                "S4D_SHORT_TIME_FIXED_GAMMA_T30",
                "summary_S4D_short_time_fixed_gamma_all_reps.csv"
            ),
            subset_label = "S4D fixed-gamma"
        )
    ))
}

## ---- file discovery ---------------------------------------------------------
scenario_files_from_registry <- function(scenarios, max_reps = Inf, verbose = TRUE) {
    out <- list()
    idx <- 0L
    for (ii in seq_len(nrow(scenarios))) {
        data_path <- scenarios$data_path[ii]
        if (!dir.exists(data_path)) {
            if (isTRUE(verbose)) message("Skipping missing data folder: ", data_path)
            next
        }
        files <- list.files(data_path, pattern = "^data_rep[0-9]+\\.rds$", full.names = TRUE)
        if (length(files) == 0L) {
            if (isTRUE(verbose)) message("No data_rep*.rds files found in: ", data_path)
            next
        }
        reps <- vapply(files, .parse_rep_id_from_file, integer(1))
        ord <- order(reps)
        files <- files[ord]
        reps <- reps[ord]
        if (is.finite(max_reps)) {
            keep <- seq_along(files) <= max_reps
            files <- files[keep]
            reps <- reps[keep]
        }
        for (jj in seq_along(files)) {
            idx <- idx + 1L
            out[[idx]] <- data.frame(
                scenario_label = scenarios$scenario_label[ii],
                scenario_family = scenarios$scenario_family[ii],
                data_scenario_id = scenarios$data_scenario_id[ii],
                rep_id = reps[jj],
                data_file = files[jj],
                fit_summary_file = scenarios$fit_summary_file[ii],
                subset_label = scenarios$subset_label[ii],
                stringsAsFactors = FALSE
            )
        }
    }
    .bind_rows_aligned(out)
}

## ---- covariate geometry computations ---------------------------------------
.extract_phi <- function(dat) {
    if (!is.null(dat$phi_star_ident)) return(list(phi = as.numeric(dat$phi_star_ident), source = "phi_star_ident"))
    if (!is.null(dat$phi_ident)) return(list(phi = as.numeric(dat$phi_ident), source = "phi_ident"))
    if (!is.null(dat$phi_star)) return(list(phi = as.numeric(dat$phi_star), source = "phi_star"))
    if (!is.null(dat$phi)) return(list(phi = as.numeric(dat$phi), source = "phi"))
    stop("Could not find a spatial effect vector in this data object.", call. = FALSE)
}

.extract_phi_raw <- function(dat) {
    if (!is.null(dat$phi_star)) return(as.numeric(dat$phi_star))
    if (!is.null(dat$phi)) return(as.numeric(dat$phi))
    NA_real_
}

.extract_lambda_matrix <- function(dat) {
    lam <- dat$lambda_tilde_ident %||% dat$lambda_tilde %||% dat$lambda_star %||% NULL
    if (is.null(lam)) return(NULL)
    as.matrix(lam)
}

.orient_time_region_matrix <- function(x, n_region, name = "x") {
    x <- as.matrix(x)
    if (ncol(x) == n_region) return(x)
    if (nrow(x) == n_region) return(t(x))
    stop(sprintf(
        "Cannot orient %s: neither nrow nor ncol equals number of regions (%s).",
        name, n_region
    ), call. = FALSE)
}

.lag1_by_region <- function(x) {
    if (nrow(x) < 3L) return(rep(NA_real_, ncol(x)))
    vapply(seq_len(ncol(x)), function(j) {
        .safe_cor(x[-nrow(x), j], x[-1L, j])
    }, numeric(1))
}

.geometry_one_cov <- function(x, phi_ident, phi_raw, log_lambda, covariate_name) {
    x_vec <- as.numeric(x)
    x_bar <- colMeans(x, na.rm = TRUE)
    x_lag1 <- .lag1_by_region(x)
    within_vars <- apply(x, 2, stats::var, na.rm = TRUE)
    between_var <- stats::var(x_bar, na.rm = TRUE)
    within_var <- mean(within_vars, na.rm = TRUE)
    bw_ratio <- if (is.finite(within_var) && within_var > 0) between_var / within_var else NA_real_

    out <- data.frame(
        covariate = covariate_name,
        global_mean = mean(x_vec, na.rm = TRUE),
        global_sd = stats::sd(x_vec, na.rm = TRUE),
        mean_lag1 = .safe_mean(x_lag1),
        sd_lag1 = .safe_sd(x_lag1),
        min_lag1 = .safe_min(x_lag1),
        max_lag1 = .safe_max(x_lag1),
        between_region_var = between_var,
        within_region_var = within_var,
        between_within_ratio = bw_ratio,
        cor_xbar_phi_ident = .safe_cor(x_bar, phi_ident),
        abs_cor_xbar_phi_ident = abs(.safe_cor(x_bar, phi_ident)),
        stringsAsFactors = FALSE
    )

    if (length(phi_raw) == length(phi_ident) && all(is.finite(phi_raw))) {
        out$cor_xbar_phi_raw <- .safe_cor(x_bar, phi_raw)
        out$abs_cor_xbar_phi_raw <- abs(out$cor_xbar_phi_raw)
    } else {
        out$cor_xbar_phi_raw <- NA_real_
        out$abs_cor_xbar_phi_raw <- NA_real_
    }

    if (!is.null(log_lambda) && all(dim(log_lambda) == dim(x))) {
        out$cor_vec_x_log_lambda <- .safe_cor(x_vec, as.numeric(log_lambda))
        out$abs_cor_vec_x_log_lambda <- abs(out$cor_vec_x_log_lambda)

        lambda_bar_region <- colMeans(log_lambda, na.rm = TRUE)
        lambda_bar_time <- rowMeans(log_lambda, na.rm = TRUE)
        out$cor_xbar_lambdabar_region <- .safe_cor(x_bar, lambda_bar_region)
        out$abs_cor_xbar_lambdabar_region <- abs(out$cor_xbar_lambdabar_region)
        out$cor_xtime_lambdabar_time <- .safe_cor(rowMeans(x, na.rm = TRUE), lambda_bar_time)
        out$abs_cor_xtime_lambdabar_time <- abs(out$cor_xtime_lambdabar_time)
    } else {
        out$cor_vec_x_log_lambda <- NA_real_
        out$abs_cor_vec_x_log_lambda <- NA_real_
        out$cor_xbar_lambdabar_region <- NA_real_
        out$abs_cor_xbar_lambdabar_region <- NA_real_
        out$cor_xtime_lambdabar_time <- NA_real_
        out$abs_cor_xtime_lambdabar_time <- NA_real_
    }

    out
}

geometry_from_dataset <- function(data_file,
                                  scenario_label = NA_character_,
                                  scenario_family = NA_character_,
                                  data_scenario_id = NA_character_,
                                  rep_id = NA_integer_) {
    dat <- readRDS(data_file)
    phi_info <- .extract_phi(dat)
    phi_ident <- phi_info$phi
    phi_raw <- .extract_phi_raw(dat)
    n_region <- length(phi_ident)

    if (is.null(dat$x1) || is.null(dat$x2)) {
        stop("Data object must contain x1 and x2.", call. = FALSE)
    }

    x1 <- .orient_time_region_matrix(dat$x1, n_region, "x1")
    x2 <- .orient_time_region_matrix(dat$x2, n_region, "x2")
    TT <- nrow(x1)
    n1 <- ncol(x1)

    lam <- .extract_lambda_matrix(dat)
    if (!is.null(lam)) {
        lam <- .orient_time_region_matrix(lam, n_region, "lambda")
        log_lam <- suppressWarnings(log(pmax(lam, .Machine$double.xmin)))
    } else {
        log_lam <- NULL
    }

    x1g <- .geometry_one_cov(x1, phi_ident, phi_raw, log_lam, "x1")
    x2g <- .geometry_one_cov(x2, phi_ident, phi_raw, log_lam, "x2")

    ## Cross-covariate collinearity.
    cor_vec_x1_x2 <- .safe_cor(as.numeric(x1), as.numeric(x2))
    cor_bar_x1_x2 <- .safe_cor(colMeans(x1, na.rm = TRUE), colMeans(x2, na.rm = TRUE))
    abs_cor_vec_x1_x2 <- abs(cor_vec_x1_x2)
    abs_cor_bar_x1_x2 <- abs(cor_bar_x1_x2)

    y <- dat$y_coarse %||% dat$y %||% NULL
    y_vec <- if (!is.null(y)) as.numeric(y) else numeric(0)

    beta <- dat$beta_star_ident %||% dat$beta_star %||% c(NA_real_, NA_real_)
    beta <- as.numeric(beta)
    beta0 <- dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_

    beta1_truth <- if (length(beta) >= 1L) beta[1] else NA_real_
    beta2_truth <- if (length(beta) >= 2L) beta[2] else NA_real_

    out <- data.frame(
        scenario_label = scenario_label,
        scenario_family = scenario_family,
        data_scenario_id = data_scenario_id,
        rep_id = as.integer(rep_id),
        data_file = data_file,
        TT = TT,
        n1 = n1,
        n_cells = TT * n1,
        phi_source = phi_info$source,
        beta0_truth = as.numeric(beta0),
        beta1_truth = beta1_truth,
        beta2_truth = beta2_truth,
        mean_count = if (length(y_vec) > 0L) mean(y_vec, na.rm = TRUE) else NA_real_,
        zero_prop = if (length(y_vec) > 0L) mean(y_vec == 0, na.rm = TRUE) else NA_real_,
        count_cv = if (length(y_vec) > 0L) .cv(y_vec) else NA_real_,
        variance_to_mean = if (length(y_vec) > 0L) stats::var(y_vec, na.rm = TRUE) / mean(y_vec, na.rm = TRUE) else NA_real_,
        mean_lag1_x1 = x1g$mean_lag1,
        mean_lag1_x2 = x2g$mean_lag1,
        sd_lag1_x1 = x1g$sd_lag1,
        sd_lag1_x2 = x2g$sd_lag1,
        between_region_var_x1 = x1g$between_region_var,
        between_region_var_x2 = x2g$between_region_var,
        within_region_var_x1 = x1g$within_region_var,
        within_region_var_x2 = x2g$within_region_var,
        between_within_ratio_x1 = x1g$between_within_ratio,
        between_within_ratio_x2 = x2g$between_within_ratio,
        cor_x1bar_phi_ident = x1g$cor_xbar_phi_ident,
        cor_x2bar_phi_ident = x2g$cor_xbar_phi_ident,
        abs_cor_x1bar_phi_ident = x1g$abs_cor_xbar_phi_ident,
        abs_cor_x2bar_phi_ident = x2g$abs_cor_xbar_phi_ident,
        cor_x1bar_phi_raw = x1g$cor_xbar_phi_raw,
        cor_x2bar_phi_raw = x2g$cor_xbar_phi_raw,
        abs_cor_x1bar_phi_raw = x1g$abs_cor_xbar_phi_raw,
        abs_cor_x2bar_phi_raw = x2g$abs_cor_xbar_phi_raw,
        cor_vec_x1_log_lambda = x1g$cor_vec_x_log_lambda,
        cor_vec_x2_log_lambda = x2g$cor_vec_x_log_lambda,
        abs_cor_vec_x1_log_lambda = x1g$abs_cor_vec_x_log_lambda,
        abs_cor_vec_x2_log_lambda = x2g$abs_cor_vec_x_log_lambda,
        cor_x1bar_lambdabar_region = x1g$cor_xbar_lambdabar_region,
        cor_x2bar_lambdabar_region = x2g$cor_xbar_lambdabar_region,
        abs_cor_x1bar_lambdabar_region = x1g$abs_cor_xbar_lambdabar_region,
        abs_cor_x2bar_lambdabar_region = x2g$abs_cor_xbar_lambdabar_region,
        cor_xtime1_lambdabar_time = x1g$cor_xtime_lambdabar_time,
        cor_xtime2_lambdabar_time = x2g$cor_xtime_lambdabar_time,
        abs_cor_xtime1_lambdabar_time = x1g$abs_cor_xtime_lambdabar_time,
        abs_cor_xtime2_lambdabar_time = x2g$abs_cor_xtime_lambdabar_time,
        cor_vec_x1_x2 = cor_vec_x1_x2,
        abs_cor_vec_x1_x2 = abs_cor_vec_x1_x2,
        cor_bar_x1_x2 = cor_bar_x1_x2,
        abs_cor_bar_x1_x2 = abs_cor_bar_x1_x2,
        signal_sd_beta1_x1 = if (is.finite(beta1_truth)) stats::sd(as.numeric(beta1_truth * x1), na.rm = TRUE) else NA_real_,
        signal_sd_beta2_x2 = if (is.finite(beta2_truth)) stats::sd(as.numeric(beta2_truth * x2), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
    )

    out$lag1_x1_minus_x2 <- out$mean_lag1_x1 - out$mean_lag1_x2
    out$abs_cor_x2phi_minus_x1phi <- out$abs_cor_x2bar_phi_ident - out$abs_cor_x1bar_phi_ident
    out$bw_ratio_x2_minus_x1 <- out$between_within_ratio_x2 - out$between_within_ratio_x1
    out$signal_sd_beta2_minus_beta1 <- out$signal_sd_beta2_x2 - out$signal_sd_beta1_x1

    out
}

make_covariate_geometry_long <- function(geom) {
    if (nrow(geom) == 0L) return(data.frame())
    rows <- list()
    idx <- 0L
    for (ii in seq_len(nrow(geom))) {
        for (covariate in c("x1", "x2")) {
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
                scenario_label = geom$scenario_label[ii],
                scenario_family = geom$scenario_family[ii],
                data_scenario_id = geom$data_scenario_id[ii],
                rep_id = geom$rep_id[ii],
                covariate = covariate,
                mean_lag1 = geom[[paste0("mean_lag1_", covariate)]][ii],
                between_region_var = geom[[paste0("between_region_var_", covariate)]][ii],
                within_region_var = geom[[paste0("within_region_var_", covariate)]][ii],
                between_within_ratio = geom[[paste0("between_within_ratio_", covariate)]][ii],
                cor_xbar_phi_ident = geom[[paste0("cor_", covariate, "bar_phi_ident")]][ii],
                abs_cor_xbar_phi_ident = geom[[paste0("abs_cor_", covariate, "bar_phi_ident")]][ii],
                cor_vec_x_log_lambda = geom[[paste0("cor_vec_", covariate, "_log_lambda")]][ii],
                abs_cor_vec_x_log_lambda = geom[[paste0("abs_cor_vec_", covariate, "_log_lambda")]][ii],
                signal_sd_beta_x = if (covariate == "x1") geom$signal_sd_beta1_x1[ii] else geom$signal_sd_beta2_x2[ii],
                stringsAsFactors = FALSE
            )
        }
    }
    .bind_rows_aligned(rows)
}

summarise_geometry_by_scenario <- function(geom) {
    if (nrow(geom) == 0L) return(data.frame())
    groups <- split(geom, geom$scenario_label)
    rows <- lapply(names(groups), function(g) {
        d <- groups[[g]]
        data.frame(
            scenario_label = g,
            scenario_family = d$scenario_family[1],
            data_scenario_id = d$data_scenario_id[1],
            n_reps = nrow(d),
            TT_avg = .safe_mean(d$TT),
            n_cells_avg = .safe_mean(d$n_cells),
            mean_count_avg = .safe_mean(d$mean_count),
            zero_prop_avg = .safe_mean(d$zero_prop),
            count_cv_avg = .safe_mean(d$count_cv),
            mean_lag1_x1_avg = .safe_mean(d$mean_lag1_x1),
            mean_lag1_x2_avg = .safe_mean(d$mean_lag1_x2),
            lag1_x1_minus_x2_avg = .safe_mean(d$lag1_x1_minus_x2),
            between_within_ratio_x1_avg = .safe_mean(d$between_within_ratio_x1),
            between_within_ratio_x2_avg = .safe_mean(d$between_within_ratio_x2),
            abs_cor_x1bar_phi_avg = .safe_mean(d$abs_cor_x1bar_phi_ident),
            abs_cor_x2bar_phi_avg = .safe_mean(d$abs_cor_x2bar_phi_ident),
            abs_cor_x2phi_minus_x1phi_avg = .safe_mean(d$abs_cor_x2phi_minus_x1phi),
            abs_cor_vec_x1_log_lambda_avg = .safe_mean(d$abs_cor_vec_x1_log_lambda),
            abs_cor_vec_x2_log_lambda_avg = .safe_mean(d$abs_cor_vec_x2_log_lambda),
            abs_cor_vec_x1_x2_avg = .safe_mean(d$abs_cor_vec_x1_x2),
            abs_cor_bar_x1_x2_avg = .safe_mean(d$abs_cor_bar_x1_x2),
            signal_sd_beta1_x1_avg = .safe_mean(d$signal_sd_beta1_x1),
            signal_sd_beta2_x2_avg = .safe_mean(d$signal_sd_beta2_x2),
            stringsAsFactors = FALSE
        )
    })
    .bind_rows_aligned(rows)
}

## ---- beta summary merge -----------------------------------------------------
.read_fit_summaries_from_registry <- function(scenarios, verbose = TRUE) {
    rows <- list()
    idx <- 0L
    for (ii in seq_len(nrow(scenarios))) {
        ff <- scenarios$fit_summary_file[ii]
        if (!is.finite(match(ff, ff)) || is.na(ff) || !nzchar(ff) || !file.exists(ff)) {
            if (isTRUE(verbose)) message("No fit summary found for ", scenarios$scenario_label[ii], ": ", ff)
            next
        }
        dat <- tryCatch(utils::read.csv(ff, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(dat) || nrow(dat) == 0L) next
        idx <- idx + 1L
        dat$scenario_label <- scenarios$scenario_label[ii]
        dat$scenario_family <- scenarios$scenario_family[ii]
        dat$subset_label <- scenarios$subset_label[ii]
        rows[[idx]] <- dat
    }
    .bind_rows_aligned(rows)
}

make_beta_geometry_by_rep <- function(geom, fit_summaries) {
    if (nrow(geom) == 0L || nrow(fit_summaries) == 0L) return(data.frame())
    if (!"rep_id" %in% names(fit_summaries)) return(data.frame())

    fit_summaries$rep_id <- suppressWarnings(as.integer(as.character(fit_summaries$rep_id)))

    key <- paste(geom$scenario_label, geom$rep_id, sep = "___")
    geom_map <- split(geom, key)

    rows <- list()
    idx <- 0L
    for (ii in seq_len(nrow(fit_summaries))) {
        kk <- paste(fit_summaries$scenario_label[ii], fit_summaries$rep_id[ii], sep = "___")
        if (!kk %in% names(geom_map)) next
        g <- geom_map[[kk]][1, , drop = FALSE]

        for (param in c("beta1", "beta2")) {
            mean_col <- paste0(param, "_mean")
            if (!mean_col %in% names(fit_summaries)) next
            post_mean <- .safe_num(fit_summaries[[mean_col]][ii])
            truth <- if (param == "beta1") g$beta1_truth else g$beta2_truth
            covariate <- if (param == "beta1") "x1" else "x2"
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
                scenario_label = g$scenario_label,
                scenario_family = g$scenario_family,
                rep_id = g$rep_id,
                parameter = param,
                covariate = covariate,
                beta_truth = truth,
                beta_post_mean = post_mean,
                beta_bias = post_mean - truth,
                abs_beta_bias = abs(post_mean - truth),
                mean_lag1 = if (covariate == "x1") g$mean_lag1_x1 else g$mean_lag1_x2,
                between_within_ratio = if (covariate == "x1") g$between_within_ratio_x1 else g$between_within_ratio_x2,
                abs_cor_xbar_phi_ident = if (covariate == "x1") g$abs_cor_x1bar_phi_ident else g$abs_cor_x2bar_phi_ident,
                abs_cor_vec_x_log_lambda = if (covariate == "x1") g$abs_cor_vec_x1_log_lambda else g$abs_cor_vec_x2_log_lambda,
                signal_sd_beta_x = if (covariate == "x1") g$signal_sd_beta1_x1 else g$signal_sd_beta2_x2,
                abs_cor_vec_x1_x2 = g$abs_cor_vec_x1_x2,
                abs_cor_bar_x1_x2 = g$abs_cor_bar_x1_x2,
                stringsAsFactors = FALSE
            )
        }
    }
    .bind_rows_aligned(rows)
}

make_beta_geometry_associations <- function(beta_geom) {
    if (nrow(beta_geom) == 0L) return(data.frame())
    metrics <- c(
        "mean_lag1",
        "between_within_ratio",
        "abs_cor_xbar_phi_ident",
        "abs_cor_vec_x_log_lambda",
        "signal_sd_beta_x",
        "abs_cor_vec_x1_x2",
        "abs_cor_bar_x1_x2"
    )

    groups <- split(beta_geom, paste(beta_geom$scenario_label, beta_geom$parameter, sep = "___"))
    rows <- list()
    idx <- 0L
    for (gname in names(groups)) {
        d <- groups[[gname]]
        for (metric in metrics) {
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
                scenario_label = d$scenario_label[1],
                scenario_family = d$scenario_family[1],
                parameter = d$parameter[1],
                covariate = d$covariate[1],
                metric = metric,
                n = sum(is.finite(.safe_num(d$abs_beta_bias)) & is.finite(.safe_num(d[[metric]]))),
                cor_abs_bias_metric = .safe_cor(d$abs_beta_bias, d[[metric]]),
                cor_signed_bias_metric = .safe_cor(d$beta_bias, d[[metric]]),
                stringsAsFactors = FALSE
            )
        }
    }
    .bind_rows_aligned(rows)
}

make_interpretation_table <- function(scenario_summary) {
    if (nrow(scenario_summary) == 0L) return(data.frame())
    out <- scenario_summary
    out$x2_more_noise_like <- out$mean_lag1_x2_avg < out$mean_lag1_x1_avg
    out$x2_more_phi_aligned <- out$abs_cor_x2bar_phi_avg > out$abs_cor_x1bar_phi_avg
    out$x1_x2_global_collinearity_flag <- out$abs_cor_vec_x1_x2_avg > 0.50
    out$x1_x2_spatial_collinearity_flag <- out$abs_cor_bar_x1_x2_avg > 0.70
    out$recommended_next_check <- ifelse(
        out$x1_x2_global_collinearity_flag | out$x1_x2_spatial_collinearity_flag,
        "Check x1-x2 collinearity before changing covariate design.",
        ifelse(
            out$x2_more_noise_like & out$x2_more_phi_aligned,
            "Both temporal-noise and spatial-alignment mechanisms are plausible; use diagnostic S4E-A before S4E-B.",
            ifelse(
                out$x2_more_noise_like,
                "Temporal-structure mechanism plausible; run matched-temporal x2 diagnostic first.",
                ifelse(
                    out$x2_more_phi_aligned,
                    "Spatial-alignment mechanism plausible; controlled x2-phi alignment is well motivated.",
                    "No strong x2-specific geometry signal; consider coefficient-index swap diagnostic."
                )
            )
        )
    )
    out[, c(
        "scenario_label",
        "n_reps",
        "mean_lag1_x1_avg",
        "mean_lag1_x2_avg",
        "lag1_x1_minus_x2_avg",
        "abs_cor_x1bar_phi_avg",
        "abs_cor_x2bar_phi_avg",
        "abs_cor_x2phi_minus_x1phi_avg",
        "abs_cor_vec_x1_x2_avg",
        "abs_cor_bar_x1_x2_avg",
        "x2_more_noise_like",
        "x2_more_phi_aligned",
        "x1_x2_global_collinearity_flag",
        "x1_x2_spatial_collinearity_flag",
        "recommended_next_check"
    ), drop = FALSE]
}

## ---- main driver ------------------------------------------------------------
run_covariate_geometry_diagnostic <- function(root = ".",
                                              scenarios = NULL,
                                              output_dir = "analysis_covariate_geometry",
                                              max_reps_per_scenario = Inf,
                                              verbose = TRUE) {
    if (is.null(scenarios)) {
        scenarios <- default_covariate_geometry_scenarios(root = root)
    }

    out_dir <- file.path(root, output_dir)
    tab_dir <- file.path(out_dir, "tables")
    .ensure_dir(tab_dir)

    manifest <- scenario_files_from_registry(
        scenarios = scenarios,
        max_reps = max_reps_per_scenario,
        verbose = verbose
    )

    if (nrow(manifest) == 0L) {
        stop("No scenario data files were found. Check data folders or provide a custom scenarios data frame.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("Found ", nrow(manifest), " data files across ", length(unique(manifest$scenario_label)), " scenario(s).")
    }

    geom_list <- vector("list", nrow(manifest))
    for (ii in seq_len(nrow(manifest))) {
        if (isTRUE(verbose)) {
            message(sprintf(
                "[%d/%d] %s rep %s",
                ii, nrow(manifest), manifest$scenario_label[ii], manifest$rep_id[ii]
            ))
        }
        geom_list[[ii]] <- tryCatch(
            geometry_from_dataset(
                data_file = manifest$data_file[ii],
                scenario_label = manifest$scenario_label[ii],
                scenario_family = manifest$scenario_family[ii],
                data_scenario_id = manifest$data_scenario_id[ii],
                rep_id = manifest$rep_id[ii]
            ),
            error = function(e) {
                warning("Failed to process ", manifest$data_file[ii], ": ", conditionMessage(e), call. = FALSE)
                NULL
            }
        )
    }

    geom <- .bind_rows_aligned(geom_list)
    geom_long <- make_covariate_geometry_long(geom)
    geom_summary <- summarise_geometry_by_scenario(geom)
    interpretation <- make_interpretation_table(geom_summary)

    fit_summaries <- .read_fit_summaries_from_registry(scenarios, verbose = verbose)
    beta_geom <- make_beta_geometry_by_rep(geom, fit_summaries)
    beta_assoc <- make_beta_geometry_associations(beta_geom)

    utils::write.csv(manifest, file.path(tab_dir, "scenario_file_manifest.csv"), row.names = FALSE)
    utils::write.csv(geom, file.path(tab_dir, "covariate_geometry_by_rep.csv"), row.names = FALSE)
    utils::write.csv(geom_summary, file.path(tab_dir, "covariate_geometry_by_scenario.csv"), row.names = FALSE)
    utils::write.csv(geom_long, file.path(tab_dir, "covariate_geometry_long.csv"), row.names = FALSE)
    utils::write.csv(interpretation, file.path(tab_dir, "covariance_diagnostic_interpretation.csv"), row.names = FALSE)

    if (nrow(beta_geom) > 0L) {
        utils::write.csv(beta_geom, file.path(tab_dir, "beta_geometry_by_rep.csv"), row.names = FALSE)
        utils::write.csv(beta_assoc, file.path(tab_dir, "beta_geometry_associations.csv"), row.names = FALSE)
    }

    if (isTRUE(verbose)) {
        message("\nSaved covariate-geometry diagnostics to: ", tab_dir)
        message("\nScenario summary:")
        print(geom_summary[, intersect(c(
            "scenario_label", "n_reps", "TT_avg", "mean_count_avg",
            "mean_lag1_x1_avg", "mean_lag1_x2_avg", "lag1_x1_minus_x2_avg",
            "abs_cor_x1bar_phi_avg", "abs_cor_x2bar_phi_avg",
            "abs_cor_vec_x1_x2_avg", "abs_cor_bar_x1_x2_avg"
        ), names(geom_summary)), drop = FALSE])
        if (nrow(beta_assoc) > 0L) {
            message("\nBeta-geometry association rows written: ", nrow(beta_assoc))
        } else {
            message("\nNo beta-geometry association table was written because fit summaries were not found or did not match geometry files.")
        }
    }

    invisible(list(
        scenarios = scenarios,
        manifest = manifest,
        geometry_by_rep = geom,
        geometry_by_scenario = geom_summary,
        geometry_long = geom_long,
        interpretation = interpretation,
        fit_summaries = fit_summaries,
        beta_geometry_by_rep = beta_geom,
        beta_geometry_associations = beta_assoc,
        output_dir = out_dir,
        tables_dir = tab_dir
    ))
}

## ---- automatic run when sourced --------------------------------------------
## Set MSSTNB_COVGEOM_AUTO_RUN=0 before sourcing if you only want the functions.
if (!identical(Sys.getenv("MSSTNB_COVGEOM_AUTO_RUN"), "0")) {
    covariate_geometry_diagnostic <- run_covariate_geometry_diagnostic(
        root = ".",
        output_dir = "analysis_covariate_geometry",
        max_reps_per_scenario = Inf,
        verbose = TRUE
    )
}
