## ============================================================================
## s4c_small_r_overdispersion.R
## Scenario 4C small-r / strong-overdispersion data-generation script
## for the MSSTNB project.
##
## Scenario 4C purpose:
##   Strong overdispersion stress test corresponding to Scenario 3.
##
## Design rule:
##   1. Generate the Scenario 3 reference temporal-spatial latent structure under
##      the usual r_reference_truth, typically r = 15.
##   2. Keep the same exposure, covariates, beta, phi, omega, gamma, and latent
##      lambda_tilde path.
##   3. Replace only the negative-binomial dispersion truth r by a smaller value.
##   4. Re-generate kappa ~ Gamma(r, r), then regenerate observation counts from
##      Poisson(xi * lambda_tilde * kappa).
##
## This isolates observation-level overdispersion from S4A sparse-count stress
## and S4B low-exposure stress.  This script only generates and calibrates data;
## it does not fit the model.
## ============================================================================

.same_dim_s4c <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.cv_s4c <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) return(NA_real_)
    stats::sd(x) / mx
}

.safe_ratio_s4c <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}


.add_increase_columns_s4c <- function(df) {
    if (is.null(df) || !is.data.frame(df)) return(df)
    if (all(c("count_cv", "reference_count_cv") %in% names(df))) {
        df$count_cv_increase <- df$count_cv - df$reference_count_cv
    }
    if (all(c("variance_to_mean", "reference_variance_to_mean") %in% names(df))) {
        df$variance_to_mean_increase <- df$variance_to_mean - df$reference_variance_to_mean
    }
    if (all(c("kappa_cv", "reference_kappa_cv") %in% names(df))) {
        df$kappa_cv_increase <- df$kappa_cv - df$reference_kappa_cv
    }
    df
}

.require_file_s4c <- function(path) {
    if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
    invisible(path)
}

.source_checked_s4c <- function(path, verbose = TRUE) {
    .require_file_s4c(path)
    if (isTRUE(verbose)) message("source: ", path)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4c <- function(root = ".", s3_script = NULL) {
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

.require_object_s4c <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source Scenario 3 and project dependencies
## -----------------------------------------------------------------------------
source_s4c_small_r_overdispersion <- function(root = ".",
                                              s3_script = NULL,
                                              verbose = TRUE) {
    s3_path <- .find_s3_script_s4c(root = root, s3_script = s3_script)
    .source_checked_s4c(s3_path, verbose = verbose)

    .require_object_s4c("source_s3_dynamic_learned_gamma")
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
        message("Scenario 4C small-r overdispersion data generator loaded.")
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Count and overdispersion summaries
## -----------------------------------------------------------------------------
.count_stats_s4c <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(yy, probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        nonzero_prop = mean(yy > 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = .cv_s4c(yy),
        variance_to_mean = .safe_ratio_s4c(vv, mn),
        max_to_mean = .safe_ratio_s4c(max(yy), mn),
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

.kappa_stats_s4c <- function(kappa, prefix = "") {
    kk <- as.numeric(kappa)
    qs <- stats::quantile(kk, probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                          names = FALSE, type = 7)
    out <- data.frame(
        kappa_mean = mean(kk),
        kappa_sd = stats::sd(kk),
        kappa_cv = .cv_s4c(kk),
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

summarise_s4c_overdispersion_counts_one <- function(dat) {
    validate_s4c_overdispersion_data(dat)

    y_stats <- .count_stats_s4c(dat$y_coarse, prefix = "")
    y_ref_stats <- .count_stats_s4c(dat$y_coarse_reference, prefix = "reference_")
    k_stats <- .kappa_stats_s4c(dat$kappa, prefix = "")
    k_ref_stats <- .kappa_stats_s4c(dat$kappa_reference, prefix = "reference_")

    mu <- as.numeric(dat$mu_nb)
    r_vec <- as.numeric(dat$r_star)
    r_by_cell <- rep(r_vec, each = as.integer(dat$TT))
    theoretical_vm <- mean(1 + mu / r_by_cell)

    data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        n1 = as.integer(dat$n1),
        stress_type = dat$stress_type,
        r_reference_truth = dat$r_reference_truth,
        r_stress_truth = dat$r_stress_truth,
        r_stress_min = min(dat$r_star),
        r_stress_max = max(dat$r_star),
        r_ratio_to_reference = dat$r_stress_truth / dat$r_reference_truth,
        beta0_star = dat$beta0_star,
        beta0_star_ident = dat$beta0_star_ident,
        mean_exposure = mean(dat$e),
        min_exposure = min(dat$e),
        max_exposure = max(dat$e),
        gamma_truth_mean = mean(dat$gamma_star),
        delta_truth = dat$delta_star,
        theoretical_variance_to_mean_avg = theoretical_vm,
        lambda_raw_min = min(dat$lambda_tilde),
        lambda_raw_median = stats::median(as.numeric(dat$lambda_tilde)),
        lambda_raw_max = max(dat$lambda_tilde),
        lambda_ident_min = min(dat$lambda_tilde_ident),
        lambda_ident_median = stats::median(as.numeric(dat$lambda_tilde_ident)),
        lambda_ident_max = max(dat$lambda_tilde_ident),
        lambda_ident_log_rm_mean = mean(abs(colMeans(log(dat$lambda_tilde_ident)))),
        coherent = isTRUE(dat$coherent),
        y_stats,
        y_ref_stats,
        k_stats,
        k_ref_stats,
        stringsAsFactors = FALSE
    )
}

summarise_s4c_overdispersion_counts_from_files <- function(files) {
    if (length(files) == 0L) stop("No files supplied.", call. = FALSE)
    out <- lapply(files, function(ff) {
        dat <- readRDS(ff)
        ss <- summarise_s4c_overdispersion_counts_one(dat)
        ss$file <- ff
        ss
    })
    .add_increase_columns_s4c(do.call(rbind, out))
}

summarise_s4c_overdispersion_counts_from_dir <- function(root = ".",
                                                         data_dir = "data_s4c_overdispersion",
                                                         scenario_id = "S4C_STRONG_OVERDISPERSION_T100") {
    in_dir <- file.path(root, data_dir, scenario_id)
    files <- list.files(in_dir, pattern = "^data_rep[0-9]+\\.rds$", full.names = TRUE)
    if (length(files) == 0L) stop("No S4C data files found in: ", in_dir, call. = FALSE)
    summarise_s4c_overdispersion_counts_from_files(files)
}

## -----------------------------------------------------------------------------
## Validation
## -----------------------------------------------------------------------------
validate_s4c_overdispersion_data <- function(dat) {
    .require_object_s4c("validate_s3_data")
    validate_s3_data(dat)

    required <- c(
        "scenario_id", "stress_type", "stress_description",
        "r_reference_truth", "r_stress_truth", "r_star",
        "y_coarse_reference", "kappa_reference", "reference_count_summary"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("S4C overdispersion dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    if (!identical(dat$stress_type, "small_r_overdispersion")) {
        stop("dat$stress_type must be 'small_r_overdispersion'.", call. = FALSE)
    }
    if (!is.finite(dat$r_reference_truth) || dat$r_reference_truth <= 0) {
        stop("dat$r_reference_truth must be positive and finite.", call. = FALSE)
    }
    if (!is.finite(dat$r_stress_truth) || dat$r_stress_truth <= 0) {
        stop("dat$r_stress_truth must be positive and finite.", call. = FALSE)
    }
    if (dat$r_stress_truth >= dat$r_reference_truth) {
        stop("For S4C, r_stress_truth should be smaller than r_reference_truth.", call. = FALSE)
    }
    if (length(dat$r_star) != as.integer(dat$n1) ||
        any(!is.finite(dat$r_star)) || any(dat$r_star <= 0)) {
        stop("dat$r_star must be a positive finite vector of length n1.", call. = FALSE)
    }
    if (!.same_dim_s4c(dat$kappa, dim(dat$y_coarse)) ||
        any(!is.finite(dat$kappa)) || any(dat$kappa <= 0)) {
        stop("dat$kappa must be positive finite matrix with y_coarse dimensions.", call. = FALSE)
    }
    if (!.same_dim_s4c(dat$kappa_reference, dim(dat$y_coarse))) {
        stop("dat$kappa_reference must be a matrix with y_coarse dimensions.", call. = FALSE)
    }

    y_check <- apply(dat$y_fine, c(1, 2), sum)
    if (!isTRUE(all(dat$y_coarse == y_check))) {
        stop("S4C fine counts are not coherent with coarse counts.", call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Scenario 4C small-r overdispersion simulation
## -----------------------------------------------------------------------------
simulate_s4c_overdispersion_one <- function(seed = 1L,
                                            TT_use = TT,
                                            r_stress_truth = 3,
                                            r_reference_truth = 15,
                                            beta0_reference_truth = -1.5,
                                            scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
                                            rep_id = NA_integer_,
                                            max_poisson_rate = 1e9,
                                            ...) {
    .require_object_s4c("simulate_s3_dynamic_learned_gamma_one")
    .require_object_s4c("recenter")

    if (!is.finite(r_reference_truth) || r_reference_truth <= 0) {
        stop("r_reference_truth must be positive and finite.", call. = FALSE)
    }
    if (!is.finite(r_stress_truth) || r_stress_truth <= 0) {
        stop("r_stress_truth must be positive and finite.", call. = FALSE)
    }
    if (r_stress_truth >= r_reference_truth) {
        stop("r_stress_truth should be smaller than r_reference_truth for S4C.", call. = FALSE)
    }

    ## Generate the full Scenario 3 latent structure under the reference r.
    ref_dat <- simulate_s3_dynamic_learned_gamma_one(
        seed = seed,
        TT_use = TT_use,
        beta0_truth = beta0_reference_truth,
        r_truth = r_reference_truth,
        scenario_id = paste0(scenario_id, "_REFERENCE_S3_LATENT"),
        rep_id = rep_id,
        max_poisson_rate = max_poisson_rate,
        ...
    )

    TT_now <- as.integer(ref_dat$TT)
    n1_now <- as.integer(ref_dat$n1)
    n_children_now <- as.integer(ref_dat$n_children)

    r_vec <- rep(as.numeric(r_stress_truth), n1_now)

    ## Re-generate kappa and counts under the stronger overdispersion, while
    ## keeping xi, lambda_tilde, and omega fixed from the reference S3 latent path.
    set.seed(as.integer(seed) + 999000L + as.integer(round(1000 * r_stress_truth)))

    kappa_s4c <- matrix(NA_real_, nrow = TT_now, ncol = n1_now)
    for (j in seq_len(n1_now)) {
        kappa_s4c[, j] <- stats::rgamma(TT_now, shape = r_vec[j], rate = r_vec[j])
    }

    mu_nb_s4c <- ref_dat$xi * ref_dat$lambda_tilde
    poisson_rate_s4c <- mu_nb_s4c * kappa_s4c

    bad <- !is.finite(poisson_rate_s4c) | poisson_rate_s4c < 0 |
        poisson_rate_s4c > max_poisson_rate
    if (any(bad)) {
        idx <- which(bad, arr.ind = TRUE)[1, ]
        stop(sprintf(
            paste0(
                "Bad S4C Poisson rate. First bad cell: t=%d, j=%d, rate=%s, ",
                "lambda=%s, xi=%s, kappa=%s."
            ),
            idx[1], idx[2],
            as.character(poisson_rate_s4c[idx[1], idx[2]]),
            as.character(ref_dat$lambda_tilde[idx[1], idx[2]]),
            as.character(ref_dat$xi[idx[1], idx[2]]),
            as.character(kappa_s4c[idx[1], idx[2]])
        ), call. = FALSE)
    }

    y_coarse_s4c <- matrix(
        stats::rpois(TT_now * n1_now, lambda = as.numeric(poisson_rate_s4c)),
        nrow = TT_now,
        ncol = n1_now
    )
    if (anyNA(y_coarse_s4c)) stop("rpois produced NA values in S4C.", call. = FALSE)

    y_fine_s4c <- array(0L, dim = c(TT_now, n1_now, n_children_now))
    for (t in seq_len(TT_now)) {
        for (j in seq_len(n1_now)) {
            if (y_coarse_s4c[t, j] > 0L) {
                y_fine_s4c[t, j, ] <- as.integer(stats::rmultinom(
                    1L,
                    size = y_coarse_s4c[t, j],
                    prob = ref_dat$omega[t, j, ]
                ))
            }
        }
    }

    coherent <- all(y_coarse_s4c == apply(y_fine_s4c, c(1, 2), sum))
    if (!coherent) stop("S4C fine counts are not coherent.", call. = FALSE)

    ## Identified scale is unchanged because beta0, phi, and lambda_tilde are
    ## unchanged from the reference latent structure.
    rc_truth <- recenter(
        beta0 = beta0_reference_truth,
        phi = ref_dat$phi_star,
        lambda_tilde = ref_dat$lambda_tilde,
        return_diag = TRUE
    )

    dat <- ref_dat
    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_small_r_overdispersion"
    dat$stress_type <- "small_r_overdispersion"
    dat$stress_description <- paste0(
        "Scenario 3 latent structure with observations regenerated under smaller ",
        "negative-binomial r. This isolates strong overdispersion from sparse-count ",
        "or low-exposure stresses."
    )

    dat$y_coarse_reference <- ref_dat$y_coarse
    dat$y_fine_reference <- ref_dat$y_fine
    dat$kappa_reference <- ref_dat$kappa
    dat$poisson_rate_reference <- ref_dat$poisson_rate
    dat$reference_count_summary <- as.list(.count_stats_s4c(ref_dat$y_coarse, prefix = ""))
    dat$reference_kappa_summary <- as.list(.kappa_stats_s4c(ref_dat$kappa, prefix = ""))

    dat$y_coarse <- y_coarse_s4c
    dat$y_fine <- y_fine_s4c
    dat$y_levels <- list(y_coarse_s4c, y_fine_s4c)
    dat$kappa <- kappa_s4c
    dat$xi_reference <- ref_dat$xi
    dat$mu_nb_reference <- ref_dat$mu_nb
    dat$xi <- ref_dat$xi
    dat$mu_nb <- mu_nb_s4c
    dat$poisson_rate <- poisson_rate_s4c

    dat$r_reference_truth <- as.numeric(r_reference_truth)
    dat$r_stress_truth <- as.numeric(r_stress_truth)
    dat$r_star <- r_vec

    dat$beta0_star <- beta0_reference_truth
    dat$beta_star <- ref_dat$beta_star
    dat$phi_star <- ref_dat$phi_star
    dat$beta0_star_ident <- rc_truth$beta0
    dat$beta_star_ident <- ref_dat$beta_star
    dat$phi_star_ident <- rc_truth$phi
    dat$lambda_tilde_ident <- rc_truth$lambda_tilde
    dat$lambda_recenter_diag <- rc_truth$diag

    dat$mean_count <- mean(y_coarse_s4c)
    dat$median_count <- stats::median(as.numeric(y_coarse_s4c))
    dat$zero_prop <- mean(y_coarse_s4c == 0)
    dat$nonzero_prop <- mean(y_coarse_s4c > 0)
    dat$total_count <- sum(y_coarse_s4c)
    dat$max_count <- max(y_coarse_s4c)
    dat$count_quantiles <- stats::quantile(
        as.numeric(y_coarse_s4c),
        probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1),
        names = TRUE
    )
    dat$count_sd <- stats::sd(as.numeric(y_coarse_s4c))
    dat$count_cv <- .cv_s4c(y_coarse_s4c)
    dat$variance_to_mean <- .safe_ratio_s4c(stats::var(as.numeric(y_coarse_s4c)),
                                            mean(y_coarse_s4c))
    dat$coherent <- coherent

    validate_s4c_overdispersion_data(dat)
    dat
}

simulate_s4c_overdispersion_batch <- function(reps = 1:10,
                                               data_dir = "data_s4c_overdispersion",
                                               scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
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

            dat <- simulate_s4c_overdispersion_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                message(sprintf(
                    paste0(
                        "Saved %s | r=%.3f mean_count=%.2f zero_prop=%.3f ",
                        "cv=%.2f v/m=%.2f max=%s"
                    ),
                    out_file,
                    dat$r_stress_truth,
                    dat$mean_count,
                    dat$zero_prop,
                    dat$count_cv,
                    dat$variance_to_mean,
                    as.character(dat$max_count)
                ))
            }
        }

        row <- summarise_s4c_overdispersion_counts_one(dat)
        row$file <- out_file
        manifest[[rr]] <- row
    }

    manifest_df <- .add_increase_columns_s4c(do.call(rbind, manifest))
    if (is.null(manifest_name)) manifest_name <- paste0("manifest_", scenario_id, ".csv")
    manifest_file <- file.path(out_dir, manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message(sprintf(
            paste0(
                "S4C overdispersion summary | reps=%d r=%.3f mean_count=%.2f ",
                "zero_prop=%.3f cv=%.2f v/m=%.2f max_count_max=%s"
            ),
            nrow(manifest_df),
            unique(manifest_df$r_stress_truth)[1],
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            mean(manifest_df$count_cv),
            mean(manifest_df$variance_to_mean),
            as.character(max(manifest_df$max_count))
        ))
    }
    invisible(manifest_df)
}

run_s4c_overdispersion_data_generation <- function(root = ".",
                                                   s3_script = NULL,
                                                   reps = 1:10,
                                                   TT_use = 100,
                                                   r_stress_truth = 3,
                                                   r_reference_truth = 15,
                                                   beta0_reference_truth = -1.5,
                                                   data_dir = "data_s4c_overdispersion",
                                                   scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
                                                   seed_base = NULL,
                                                   overwrite_existing = TRUE,
                                                   verbose = TRUE,
                                                   ...) {
    source_s4c_small_r_overdispersion(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4c_overdispersion_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_use = TT_use,
        r_stress_truth = r_stress_truth,
        r_reference_truth = r_reference_truth,
        beta0_reference_truth = beta0_reference_truth,
        ...
    )
}

## -----------------------------------------------------------------------------
## Data-quality checks and calibration helpers
## -----------------------------------------------------------------------------
check_s4c_overdispersion_data_summary <- function(manifest,
                                                  target_mean_count_range = c(80, 350),
                                                  target_zero_prop_max = 0.20,
                                                  minimum_cv_increase = 0.10,
                                                  minimum_vm_increase = 0.20,
                                                  target_abs_beta0_ident_max = 20,
                                                  max_count_max_limit = Inf) {
    if (is.character(manifest)) manifest <- read.csv(manifest)

    required <- c(
        "rep_id", "TT", "n1", "r_reference_truth", "r_stress_truth",
        "mean_count", "median_count", "zero_prop", "total_count", "max_count",
        "count_sd", "count_cv", "variance_to_mean", "q95_count", "q99_count",
        "reference_mean_count", "reference_count_cv", "reference_variance_to_mean",
        "reference_max_count", "beta0_star_ident", "lambda_raw_median", "lambda_raw_max",
        "kappa_cv", "reference_kappa_cv"
    )
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0L) {
        stop("Manifest is missing required columns: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    out <- data.frame(
        n_reps = nrow(manifest),
        TT_unique = paste(sort(unique(manifest$TT)), collapse = ","),
        n1_unique = paste(sort(unique(manifest$n1)), collapse = ","),
        r_reference_truth_unique = paste(signif(sort(unique(manifest$r_reference_truth)), 6), collapse = ","),
        r_stress_truth_unique = paste(signif(sort(unique(manifest$r_stress_truth)), 6), collapse = ","),
        r_ratio_to_reference_avg = mean(manifest$r_ratio_to_reference),
        mean_count_avg = mean(manifest$mean_count),
        mean_count_min = min(manifest$mean_count),
        mean_count_max = max(manifest$mean_count),
        reference_mean_count_avg = mean(manifest$reference_mean_count),
        median_count_avg = mean(manifest$median_count),
        zero_prop_avg = mean(manifest$zero_prop),
        zero_prop_min = min(manifest$zero_prop),
        zero_prop_max = max(manifest$zero_prop),
        total_count_sum = sum(manifest$total_count),
        max_count_avg = mean(manifest$max_count),
        max_count_max = max(manifest$max_count),
        reference_max_count_avg = mean(manifest$reference_max_count),
        reference_max_count_max = max(manifest$reference_max_count),
        count_sd_avg = mean(manifest$count_sd),
        count_cv_avg = mean(manifest$count_cv),
        reference_count_cv_avg = mean(manifest$reference_count_cv),
        count_cv_increase_avg = mean(manifest$count_cv - manifest$reference_count_cv),
        variance_to_mean_avg = mean(manifest$variance_to_mean),
        reference_variance_to_mean_avg = mean(manifest$reference_variance_to_mean),
        variance_to_mean_increase_avg = mean(manifest$variance_to_mean - manifest$reference_variance_to_mean),
        q95_count_avg = mean(manifest$q95_count),
        q99_count_avg = mean(manifest$q99_count),
        q99_count_max = max(manifest$q99_count),
        kappa_cv_avg = mean(manifest$kappa_cv),
        reference_kappa_cv_avg = mean(manifest$reference_kappa_cv),
        kappa_cv_increase_avg = mean(manifest$kappa_cv - manifest$reference_kappa_cv),
        theoretical_variance_to_mean_avg = mean(manifest$theoretical_variance_to_mean_avg),
        beta0_ident_min = min(manifest$beta0_star_ident),
        beta0_ident_median = stats::median(manifest$beta0_star_ident),
        beta0_ident_max = max(manifest$beta0_star_ident),
        max_abs_beta0_ident = max(abs(manifest$beta0_star_ident)),
        lambda_raw_median_avg = mean(manifest$lambda_raw_median),
        lambda_raw_max_max = max(manifest$lambda_raw_max),
        stringsAsFactors = FALSE
    )

    out$passes_mean_count_target <-
        out$mean_count_avg >= target_mean_count_range[1] &&
        out$mean_count_avg <= target_mean_count_range[2]
    out$passes_zero_prop_target <- out$zero_prop_avg <= target_zero_prop_max
    out$passes_overdispersion_cv_target <- out$count_cv_increase_avg >= minimum_cv_increase
    out$passes_overdispersion_vm_target <- out$variance_to_mean_increase_avg >= minimum_vm_increase
    out$passes_identified_scale_target <- out$max_abs_beta0_ident <= target_abs_beta0_ident_max
    out$passes_max_count_target <- out$max_count_max <= max_count_max_limit

    out$passes_s4c_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_overdispersion_cv_target) &&
        isTRUE(out$passes_overdispersion_vm_target) &&
        isTRUE(out$passes_identified_scale_target) &&
        isTRUE(out$passes_max_count_target)

    out
}

calibrate_s4c_overdispersion_grid <- function(root = ".",
                                              s3_script = NULL,
                                              reps = 1:10,
                                              TT_use = 100,
                                              r_truth_grid = c(10, 7, 5, 3, 2),
                                              r_reference_truth = 15,
                                              beta0_reference_truth = -1.5,
                                              data_dir = "data_s4c_overdispersion_calibration",
                                              base_scenario_id = "S4C_OVERDISPERSION_CALIB_T100",
                                              overwrite_existing = TRUE,
                                              verbose = TRUE,
                                              ...) {
    source_s4c_small_r_overdispersion(root = root, s3_script = s3_script, verbose = verbose)

    out <- list()
    for (rrr in r_truth_grid) {
        tag <- gsub("\\.", "p", sprintf("%.3f", as.numeric(rrr)))
        scenario_id <- paste0(base_scenario_id, "_R", tag)
        if (isTRUE(verbose)) {
            message("\n--- S4C calibration: r_stress_truth = ", rrr,
                    " | scenario_id = ", scenario_id, " ---")
        }

        manifest <- simulate_s4c_overdispersion_batch(
            reps = reps,
            data_dir = data_dir,
            scenario_id = scenario_id,
            root = root,
            overwrite_existing = overwrite_existing,
            verbose = verbose,
            TT_use = TT_use,
            r_stress_truth = as.numeric(rrr),
            r_reference_truth = r_reference_truth,
            beta0_reference_truth = beta0_reference_truth,
            ...
        )
        chk <- check_s4c_overdispersion_data_summary(manifest)
        chk$scenario_id <- scenario_id
        chk$r_stress_truth <- as.numeric(rrr)
        out[[tag]] <- chk
    }

    calib <- do.call(rbind, out)
    out_dir <- file.path(root, data_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    calib_file <- file.path(out_dir, paste0("calibration_summary_", base_scenario_id, ".csv"))
    write.csv(calib, calib_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved S4C calibration summary: ", calib_file)

    invisible(calib)
}

## Example usage:
##
## source("s4c_small_r_overdispersion.R")
## calib_s4c <- calibrate_s4c_overdispersion_grid(
##     root = ".",
##     reps = 1:10,
##     TT_use = 100,
##     r_truth_grid = c(10, 7, 5, 3, 2),
##     r_reference_truth = 15,
##     overwrite_existing = TRUE
## )
## calib_s4c[, c(
##     "r_stress_truth", "mean_count_avg", "zero_prop_avg",
##     "count_cv_avg", "reference_count_cv_avg", "count_cv_increase_avg",
##     "variance_to_mean_avg", "reference_variance_to_mean_avg",
##     "variance_to_mean_increase_avg", "max_count_max", "q99_count_avg",
##     "kappa_cv_avg", "passes_s4c_data_check"
## )]
##
## After choosing the official r, for example r = 3:
## manifest_s4c <- run_s4c_overdispersion_data_generation(
##     root = ".",
##     reps = 1:10,
##     TT_use = 100,
##     r_stress_truth = 3,
##     r_reference_truth = 15,
##     scenario_id = "S4C_STRONG_OVERDISPERSION_T100",
##     overwrite_existing = TRUE
## )
## check_s4c_overdispersion_data_summary(manifest_s4c)



## ============================================================================
## Continuous-time x2 wrapper for revised Scenario 4C
## ============================================================================
## Revised S4C official covariate design:
##   x2_mode     = "continuous_time"
##   x2_ar       = 0.50
##   x2_innov_sd = 0.80
##
## These functions keep the original small-r / strong-overdispersion mechanism
## unchanged and only force the covariate design to match S1--S3 and the revised
## S4A/S4B/S4D workflows.
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.summarise_s4c_x2_design <- function(dat) {
    x2 <- dat$x2
    if (is.null(x2)) {
        return(list(
            x2_found = FALSE,
            x2_mean = NA_real_,
            x2_sd = NA_real_,
            x2_min = NA_real_,
            x2_max = NA_real_,
            x2_binary_like_prop = NA_real_,
            x2_empirical_ar1_mean = NA_real_,
            x2_empirical_ar1_median = NA_real_
        ))
    }

    x2_mat <- as.matrix(x2)
    ar_vals <- rep(NA_real_, ncol(x2_mat))
    if (nrow(x2_mat) >= 3L) {
        for (jj in seq_len(ncol(x2_mat))) {
            ar_vals[jj] <- suppressWarnings(stats::cor(
                x2_mat[-1L, jj],
                x2_mat[-nrow(x2_mat), jj],
                use = "complete.obs"
            ))
        }
    }
    near_zero_one <- abs(as.numeric(x2_mat)) < 1e-10 |
        abs(as.numeric(x2_mat) - 1) < 1e-10

    list(
        x2_found = TRUE,
        x2_mean = mean(x2_mat),
        x2_sd = stats::sd(as.numeric(x2_mat)),
        x2_min = min(x2_mat),
        x2_max = max(x2_mat),
        x2_binary_like_prop = mean(near_zero_one),
        x2_empirical_ar1_mean = mean(ar_vals, na.rm = TRUE),
        x2_empirical_ar1_median = stats::median(ar_vals, na.rm = TRUE)
    )
}

validate_s4c_overdispersion_continuous_x2_data <- function(dat,
                                                           expected_mode = "continuous_time",
                                                           expected_ar = 0.50,
                                                           expected_innov_sd = 0.80,
                                                           strict = TRUE) {
    validate_s4c_overdispersion_data(dat)

    required <- c("x2_mode", "x2_ar", "x2_innov_sd", "x2")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop(
            "S4C continuous-time x2 data is missing fields: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    if (!identical(dat$x2_mode, expected_mode)) {
        stop("dat$x2_mode must be '", expected_mode, "'.", call. = FALSE)
    }
    if (!isTRUE(all.equal(dat$x2_ar, expected_ar))) {
        stop("dat$x2_ar must be ", expected_ar, ".", call. = FALSE)
    }
    if (!isTRUE(all.equal(dat$x2_innov_sd, expected_innov_sd))) {
        stop("dat$x2_innov_sd must be ", expected_innov_sd, ".", call. = FALSE)
    }

    xs <- .summarise_s4c_x2_design(dat)
    if (isTRUE(strict)) {
        if (!isTRUE(xs$x2_found)) {
            stop("x2 could not be found in dat.", call. = FALSE)
        }
        if (!is.finite(xs$x2_sd) || xs$x2_sd < 0.10) {
            stop("x2 empirical SD is too small; x2 may not be continuous.", call. = FALSE)
        }
        if (is.finite(xs$x2_binary_like_prop) && xs$x2_binary_like_prop > 0.25) {
            stop("x2 appears too close to a binary/indicator covariate.", call. = FALSE)
        }
    }

    invisible(TRUE)
}

summarise_s4c_overdispersion_counts_one_continuous_x2 <- function(dat) {
    validate_s4c_overdispersion_continuous_x2_data(dat)
    base <- summarise_s4c_overdispersion_counts_one(dat)
    xs <- .summarise_s4c_x2_design(dat)

    x2_df <- data.frame(
        x2_mode = dat$x2_mode,
        x2_ar = dat$x2_ar,
        x2_innov_sd = dat$x2_innov_sd,
        x2_found = xs$x2_found,
        x2_mean = xs$x2_mean,
        x2_sd = xs$x2_sd,
        x2_min = xs$x2_min,
        x2_max = xs$x2_max,
        x2_binary_like_prop = xs$x2_binary_like_prop,
        x2_empirical_ar1_mean = xs$x2_empirical_ar1_mean,
        x2_empirical_ar1_median = xs$x2_empirical_ar1_median,
        stringsAsFactors = FALSE
    )

    cbind(base, x2_df)
}

simulate_s4c_overdispersion_continuous_x2_batch <- function(reps = 1:10,
                                                            data_dir = "data_s4c_overdispersion_continuous_x2",
                                                            scenario_id = "S4C_STRONG_OVERDISPERSION_CONTINUOUS_X2_T100",
                                                            root = ".",
                                                            seed_base = NULL,
                                                            overwrite_existing = TRUE,
                                                            verbose = TRUE,
                                                            manifest_name = NULL,
                                                            x2_mode = "continuous_time",
                                                            x2_ar = 0.50,
                                                            x2_innov_sd = 0.80,
                                                            ...) {
    if (!identical(x2_mode, "continuous_time")) {
        stop("For revised S4C, x2_mode must be exactly 'continuous_time'.", call. = FALSE)
    }
    if (!isTRUE(all.equal(x2_ar, 0.50))) {
        stop("For revised S4C, x2_ar must be 0.50.", call. = FALSE)
    }
    if (!isTRUE(all.equal(x2_innov_sd, 0.80))) {
        stop("For revised S4C, x2_innov_sd must be 0.80.", call. = FALSE)
    }

    out_dir <- file.path(root, data_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    manifest <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        out_file <- file.path(out_dir, paste0("data_rep", rr, ".rds"))

        if (file.exists(out_file) && !isTRUE(overwrite_existing)) {
            if (isTRUE(verbose)) message("Skipping existing file: ", out_file)
            dat <- readRDS(out_file)
            validate_s4c_overdispersion_continuous_x2_data(dat)
        } else {
            seed <- if (!is.null(seed_base)) {
                as.integer(seed_base + rep_id)
            } else if (exists("REP_SEEDS", envir = .GlobalEnv) && rep_id <= length(REP_SEEDS)) {
                as.integer(REP_SEEDS[rep_id])
            } else {
                as.integer(2026000L + rep_id)
            }

            dat <- simulate_s4c_overdispersion_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                x2_mode = x2_mode,
                x2_ar = x2_ar,
                x2_innov_sd = x2_innov_sd,
                ...
            )

            ## The Scenario 3 generator should create these fields.  We set them
            ## only if they are missing, then validate strictly.
            if (is.null(dat$x2_mode)) dat$x2_mode <- x2_mode
            if (is.null(dat$x2_ar)) dat$x2_ar <- x2_ar
            if (is.null(dat$x2_innov_sd)) dat$x2_innov_sd <- x2_innov_sd

            dat$s4c_continuous_x2_config <- list(
                x2_mode = x2_mode,
                x2_ar = x2_ar,
                x2_innov_sd = x2_innov_sd,
                r_reference_truth = dat$r_reference_truth,
                r_stress_truth = dat$r_stress_truth
            )
            dat$scenario_info <- c(dat$scenario_info %||% list(), list(
                scenario_family = "S4C",
                covariate_design = "continuous_time_x2",
                stress_test = "small_r_overdispersion"
            ))

            validate_s4c_overdispersion_continuous_x2_data(dat)
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                row_tmp <- summarise_s4c_overdispersion_counts_one_continuous_x2(dat)
                message(sprintf(
                    paste0(
                        "Saved %s | mean_count=%.2f zero_prop=%.3f ",
                        "count_cv=%.3f ref_cv=%.3f kappa_cv=%.3f x2_ar_emp=%.3f"
                    ),
                    out_file,
                    row_tmp$mean_count,
                    row_tmp$zero_prop,
                    row_tmp$count_cv,
                    row_tmp$reference_count_cv,
                    row_tmp$kappa_cv,
                    row_tmp$x2_empirical_ar1_mean
                ))
            }
        }

        row <- summarise_s4c_overdispersion_counts_one_continuous_x2(dat)
        row$file <- out_file
        manifest[[rr]] <- row
    }

    manifest_df <- do.call(rbind, manifest)

    if (is.null(manifest_name)) manifest_name <- paste0("manifest_", scenario_id, ".csv")
    manifest_file <- file.path(out_dir, manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message(sprintf(
            paste0(
                "S4C continuous-time x2 summary | reps=%d mean_count=%.2f zero_prop=%.3f ",
                "count_cv=%.3f ref_cv=%.3f vm=%.2f ref_vm=%.2f x2_emp_ar=%.3f"
            ),
            nrow(manifest_df),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            mean(manifest_df$count_cv),
            mean(manifest_df$reference_count_cv),
            mean(manifest_df$variance_to_mean),
            mean(manifest_df$reference_variance_to_mean),
            mean(manifest_df$x2_empirical_ar1_mean)
        ))
    }

    invisible(manifest_df)
}

run_s4c_overdispersion_continuous_x2_data_generation <- function(root = ".",
                                                                 s3_script = NULL,
                                                                 reps = 1:10,
                                                                 TT_use = 100,
                                                                 r_stress_truth = 3,
                                                                 r_reference_truth = 15,
                                                                 beta0_reference_truth = -1.5,
                                                                 data_dir = "data_s4c_overdispersion_continuous_x2",
                                                                 scenario_id = "S4C_STRONG_OVERDISPERSION_CONTINUOUS_X2_T100",
                                                                 seed_base = NULL,
                                                                 overwrite_existing = TRUE,
                                                                 verbose = TRUE,
                                                                 x2_mode = "continuous_time",
                                                                 x2_ar = 0.50,
                                                                 x2_innov_sd = 0.80,
                                                                 ...) {
    source_s4c_small_r_overdispersion(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4c_overdispersion_continuous_x2_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_use = TT_use,
        r_stress_truth = r_stress_truth,
        r_reference_truth = r_reference_truth,
        beta0_reference_truth = beta0_reference_truth,
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        ...
    )
}

check_s4c_overdispersion_continuous_x2_data_summary <- function(manifest,
                                                                target_mean_count_range = c(80, 350),
                                                                target_zero_prop_max = 0.20,
                                                                minimum_cv_increase = 0.10,
                                                                minimum_vm_increase = 0.20,
                                                                target_abs_beta0_ident_max = 20,
                                                                max_count_max_limit = Inf) {
    if (is.character(manifest)) manifest <- read.csv(manifest)

    base_required <- c(
        "rep_id", "TT", "n1", "r_reference_truth", "r_stress_truth",
        "mean_count", "median_count", "zero_prop", "total_count", "max_count",
        "count_sd", "count_cv", "variance_to_mean", "q95_count", "q99_count",
        "reference_mean_count", "reference_count_cv", "reference_variance_to_mean",
        "reference_max_count", "beta0_star_ident", "lambda_raw_median", "lambda_raw_max",
        "kappa_cv", "reference_kappa_cv"
    )
    x2_required <- c("x2_mode", "x2_ar", "x2_innov_sd", "x2_sd", "x2_binary_like_prop", "x2_empirical_ar1_mean")
    missing <- setdiff(c(base_required, x2_required), names(manifest))
    if (length(missing) > 0L) {
        stop("Manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    out <- check_s4c_overdispersion_data_summary(
        manifest,
        target_mean_count_range = target_mean_count_range,
        target_zero_prop_max = target_zero_prop_max,
        minimum_cv_increase = minimum_cv_increase,
        minimum_vm_increase = minimum_vm_increase,
        target_abs_beta0_ident_max = target_abs_beta0_ident_max,
        max_count_max_limit = max_count_max_limit
    )

    out$x2_mode_unique <- paste(sort(unique(manifest$x2_mode)), collapse = ",")
    out$x2_ar_unique <- paste(sort(unique(manifest$x2_ar)), collapse = ",")
    out$x2_innov_sd_unique <- paste(sort(unique(manifest$x2_innov_sd)), collapse = ",")
    out$x2_sd_avg <- mean(manifest$x2_sd)
    out$x2_binary_like_prop_max <- max(manifest$x2_binary_like_prop)
    out$x2_empirical_ar1_mean_avg <- mean(manifest$x2_empirical_ar1_mean)
    out$passes_continuous_x2_target <-
        identical(out$x2_mode_unique, "continuous_time") &&
        identical(out$x2_ar_unique, "0.5") &&
        identical(out$x2_innov_sd_unique, "0.8") &&
        is.finite(out$x2_sd_avg) && out$x2_sd_avg >= 0.10 &&
        is.finite(out$x2_binary_like_prop_max) && out$x2_binary_like_prop_max <= 0.25

    out$passes_s4c_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_overdispersion_cv_target) &&
        isTRUE(out$passes_overdispersion_vm_target) &&
        isTRUE(out$passes_identified_scale_target) &&
        isTRUE(out$passes_max_count_target) &&
        isTRUE(out$passes_continuous_x2_target)

    out
}

## Short aliases matching the revised workflow name.
run_s4c_continuous_x2_data_generation <- run_s4c_overdispersion_continuous_x2_data_generation
check_s4c_continuous_x2_data_summary <- check_s4c_overdispersion_continuous_x2_data_summary

## Example revised usage:
##
## source("s4c_small_r_overdispersion_continuous_x2.R")
## manifest_test <- run_s4c_continuous_x2_data_generation(
##     root = ".",
##     reps = 1,
##     TT_use = 100,
##     r_stress_truth = 3,
##     r_reference_truth = 15,
##     overwrite_existing = TRUE,
##     verbose = TRUE
## )
## check_test <- check_s4c_continuous_x2_data_summary(manifest_test)
## print(check_test)
##
## manifest_s4c <- run_s4c_continuous_x2_data_generation(
##     root = ".",
##     reps = 1:10,
##     TT_use = 100,
##     r_stress_truth = 3,
##     r_reference_truth = 15,
##     overwrite_existing = TRUE,
##     verbose = TRUE
## )
## check_s4c <- check_s4c_continuous_x2_data_summary(manifest_s4c)
## print(check_s4c)
