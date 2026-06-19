## ============================================================================
## s4b_low_exposure.R
## Scenario 4B low-exposure data-generation script for the MSSTNB project.
##
## Scenario 4B purpose:
##   Low and heterogeneous exposure stress test corresponding to Scenario 3.
##
## Design rule:
##   1. Generate the Scenario 3 latent temporal-spatial structure under the
##      reference exposure, covariates, beta, phi, lambda, kappa, gamma, and r.
##   2. Keep the same latent risk path lambda_tilde, covariates, phi, kappa,
##      omega, gamma, delta, and r.
##   3. Replace only the known exposure matrix e by a lower and heterogeneous
##      exposure matrix e_s4b = e_reference * exposure_multiplier.
##   4. Re-generate observation counts from the same latent structure under the
##      reduced exposure.
##
## This isolates the practical effect of low/heterogeneous exposure from the
## dynamic-lambda path-collapse problem.  It is intentionally different from
## lowering beta0, because beta0 changes disease-risk intensity while exposure is
## a known offset/sample-size component.
##
## Scope:
##   This script only generates and checks S4B data.  It does not fit the model.
## ============================================================================

.same_dim_s4b <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.cv_s4b <- function(x) {
    x <- as.numeric(x)
    if (length(x) == 0L || !is.finite(mean(x)) || mean(x) == 0) return(NA_real_)
    stats::sd(x) / mean(x)
}

.require_file_s4b <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked_s4b <- function(path, verbose = TRUE) {
    .require_file_s4b(path)
    if (isTRUE(verbose)) message("source: ", path)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4b <- function(root = ".", s3_script = NULL) {
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

.require_object_s4b <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source Scenario 3 and project dependencies
## -----------------------------------------------------------------------------
source_s4b_low_exposure <- function(root = ".",
                                    s3_script = NULL,
                                    verbose = TRUE) {
    s3_path <- .find_s3_script_s4b(root = root, s3_script = s3_script)
    .source_checked_s4b(s3_path, verbose = verbose)

    .require_object_s4b("source_s3_dynamic_learned_gamma")
    source_s3_dynamic_learned_gamma(root = root, verbose = verbose)

    needed <- c(
        "simulate_s3_dynamic_learned_gamma_one",
        "validate_s3_data",
        "compute_s3_xi",
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
        message("Scenario 4B low-exposure data generator loaded.")
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Exposure multiplier generator
## -----------------------------------------------------------------------------
generate_s4b_exposure_multiplier <- function(TT_use,
                                             n1_use,
                                             target_mean_multiplier = 0.05,
                                             area_log_sd = 0.75,
                                             time_log_sd = 0.08,
                                             lower_multiplier = 0.005,
                                             upper_multiplier = 0.25,
                                             seed = NULL) {
    TT_use <- as.integer(TT_use)
    n1_use <- as.integer(n1_use)

    if (!is.null(seed)) set.seed(as.integer(seed))

    if (!is.finite(target_mean_multiplier) ||
        target_mean_multiplier <= 0 || target_mean_multiplier >= 1) {
        stop("target_mean_multiplier must be in (0, 1).", call. = FALSE)
    }
    if (!is.finite(area_log_sd) || area_log_sd < 0) {
        stop("area_log_sd must be nonnegative and finite.", call. = FALSE)
    }
    if (!is.finite(time_log_sd) || time_log_sd < 0) {
        stop("time_log_sd must be nonnegative and finite.", call. = FALSE)
    }
    if (!is.finite(lower_multiplier) || !is.finite(upper_multiplier) ||
        lower_multiplier <= 0 || upper_multiplier <= lower_multiplier ||
        upper_multiplier >= 1) {
        stop("Require 0 < lower_multiplier < upper_multiplier < 1.", call. = FALSE)
    }

    ## area_raw has approximately unit mean before normalization.
    area_raw <- exp(stats::rnorm(n1_use, mean = -0.5 * area_log_sd^2, sd = area_log_sd))
    time_raw <- matrix(
        exp(stats::rnorm(TT_use * n1_use, mean = -0.5 * time_log_sd^2, sd = time_log_sd)),
        nrow = TT_use,
        ncol = n1_use
    )
    mult <- matrix(rep(area_raw, each = TT_use), nrow = TT_use, ncol = n1_use) * time_raw

    ## Normalize to the desired average multiplier, then clip extreme values.
    mult <- mult / mean(mult) * target_mean_multiplier
    mult <- pmin(pmax(mult, lower_multiplier), upper_multiplier)

    ## Re-normalize after clipping when possible, then clip one more time.
    mult <- mult / mean(mult) * target_mean_multiplier
    mult <- pmin(pmax(mult, lower_multiplier), upper_multiplier)

    if (any(!is.finite(mult)) || any(mult <= 0) || any(mult >= 1)) {
        stop("Generated exposure multipliers are non-finite or outside (0, 1).", call. = FALSE)
    }

    area_mult <- colMeans(mult)
    list(
        multiplier = mult,
        area_multiplier = area_mult,
        target_mean_multiplier = target_mean_multiplier,
        realized_mean_multiplier = mean(mult),
        realized_median_multiplier = stats::median(as.numeric(mult)),
        realized_min_multiplier = min(mult),
        realized_max_multiplier = max(mult),
        realized_cv_multiplier = .cv_s4b(mult),
        area_multiplier_cv = .cv_s4b(area_mult),
        area_log_sd = area_log_sd,
        time_log_sd = time_log_sd,
        lower_multiplier = lower_multiplier,
        upper_multiplier = upper_multiplier
    )
}

## -----------------------------------------------------------------------------
## Summaries
## -----------------------------------------------------------------------------
exposure_group_summary_s4b <- function(dat, n_groups = 3L) {
    validate_s4b_low_exposure_data(dat)
    n_groups <- as.integer(n_groups)
    if (n_groups < 2L) n_groups <- 2L

    area_mean_e <- colMeans(dat$e)
    area_mean_e_ref <- colMeans(dat$e_reference)
    area_mult <- colMeans(dat$exposure_multiplier)

    ord <- order(area_mean_e, decreasing = FALSE)
    group_id <- rep(NA_integer_, length(area_mean_e))
    group_id[ord] <- cut(seq_along(ord), breaks = n_groups, labels = FALSE)

    labels <- paste0("G", seq_len(n_groups), "_", c("lowest", rep("middle", max(n_groups - 2L, 0L)), "highest")[seq_len(n_groups)], "_exposure")
    if (length(labels) != n_groups) labels <- paste0("G", seq_len(n_groups))

    out <- lapply(seq_len(n_groups), function(g) {
        jj <- which(group_id == g)
        y_g <- dat$y_coarse[, jj, drop = FALSE]
        e_g <- dat$e[, jj, drop = FALSE]
        e_ref_g <- dat$e_reference[, jj, drop = FALSE]
        mult_g <- dat$exposure_multiplier[, jj, drop = FALSE]
        data.frame(
            scenario_id = dat$scenario_id,
            rep_id = as.integer(dat$rep_id),
            exposure_group = labels[g],
            n_areas = length(jj),
            area_ids = paste(jj, collapse = ","),
            mean_exposure = mean(e_g),
            median_exposure = stats::median(as.numeric(e_g)),
            mean_reference_exposure = mean(e_ref_g),
            mean_multiplier = mean(mult_g),
            area_multiplier_mean = mean(area_mult[jj]),
            area_multiplier_min = min(area_mult[jj]),
            area_multiplier_max = max(area_mult[jj]),
            mean_count = mean(y_g),
            zero_prop = mean(y_g == 0),
            total_count = sum(y_g),
            max_count = max(y_g),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

summarise_s4b_low_exposure_counts_one <- function(dat) {
    validate_s4b_low_exposure_data(dat)

    y <- as.numeric(dat$y_coarse)
    e <- as.numeric(dat$e)
    e_ref <- as.numeric(dat$e_reference)
    mult <- as.numeric(dat$exposure_multiplier)
    area_mean_e <- colMeans(dat$e)
    area_mean_y <- colMeans(dat$y_coarse)

    low_group <- exposure_group_summary_s4b(dat, n_groups = 3L)
    low_row <- low_group[low_group$exposure_group == "G1_lowest_exposure", , drop = FALSE]
    high_row <- low_group[low_group$exposure_group == "G3_highest_exposure", , drop = FALSE]

    data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        n1 = as.integer(dat$n1),
        stress_type = dat$stress_type,
        exposure_stress_type = dat$exposure_stress_type,
        target_mean_multiplier = dat$target_mean_exposure_multiplier,
        realized_mean_multiplier = mean(mult),
        realized_median_multiplier = stats::median(mult),
        realized_min_multiplier = min(mult),
        realized_max_multiplier = max(mult),
        realized_cv_multiplier = .cv_s4b(mult),
        area_multiplier_cv = .cv_s4b(colMeans(dat$exposure_multiplier)),
        area_log_sd = dat$exposure_area_log_sd,
        time_log_sd = dat$exposure_time_log_sd,
        mean_exposure = mean(e),
        median_exposure = stats::median(e),
        min_exposure = min(e),
        max_exposure = max(e),
        exposure_cv = .cv_s4b(e),
        reference_mean_exposure = mean(e_ref),
        reference_min_exposure = min(e_ref),
        reference_max_exposure = max(e_ref),
        exposure_mean_ratio = mean(e) / mean(e_ref),
        mean_count = mean(y),
        median_count = stats::median(y),
        zero_prop = mean(y == 0),
        nonzero_prop = mean(y > 0),
        total_count = sum(y),
        max_count = max(y),
        q05_count = as.numeric(stats::quantile(y, 0.05, names = FALSE)),
        q25_count = as.numeric(stats::quantile(y, 0.25, names = FALSE)),
        q75_count = as.numeric(stats::quantile(y, 0.75, names = FALSE)),
        q95_count = as.numeric(stats::quantile(y, 0.95, names = FALSE)),
        count_exposure_cor_cell = suppressWarnings(stats::cor(y, e)),
        count_exposure_cor_area = suppressWarnings(stats::cor(area_mean_y, area_mean_e)),
        lowest_exposure_group_mean_count = low_row$mean_count[1],
        lowest_exposure_group_zero_prop = low_row$zero_prop[1],
        highest_exposure_group_mean_count = high_row$mean_count[1],
        highest_exposure_group_zero_prop = high_row$zero_prop[1],
        gamma_truth_mean = mean(dat$gamma_star),
        delta_truth = dat$delta_star,
        r_truth_mean = mean(dat$r_star),
        beta0_star = dat$beta0_star,
        beta0_star_ident = dat$beta0_star_ident,
        lambda_raw_min = min(dat$lambda_tilde),
        lambda_raw_median = stats::median(as.numeric(dat$lambda_tilde)),
        lambda_raw_max = max(dat$lambda_tilde),
        lambda_ident_min = min(dat$lambda_tilde_ident),
        lambda_ident_median = stats::median(as.numeric(dat$lambda_tilde_ident)),
        lambda_ident_max = max(dat$lambda_tilde_ident),
        lambda_ident_log_rm_mean = mean(abs(colMeans(log(dat$lambda_tilde_ident)))),
        reference_mean_count = dat$reference_count_summary$mean_count,
        reference_zero_prop = dat$reference_count_summary$zero_prop,
        coherent = isTRUE(dat$coherent),
        stringsAsFactors = FALSE
    )
}

summarise_s4b_low_exposure_counts_from_files <- function(files) {
    if (length(files) == 0L) stop("No files supplied.", call. = FALSE)
    out <- lapply(files, function(ff) {
        dat <- readRDS(ff)
        ss <- summarise_s4b_low_exposure_counts_one(dat)
        ss$file <- ff
        ss
    })
    do.call(rbind, out)
}

summarise_s4b_low_exposure_counts_from_dir <- function(root = ".",
                                                       data_dir = "data_s4b_low_exposure",
                                                       scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100") {
    in_dir <- file.path(root, data_dir, scenario_id)
    files <- list.files(in_dir, pattern = "^data_rep[0-9]+\\.rds$", full.names = TRUE)
    if (length(files) == 0L) {
        stop("No S4B data files found in: ", in_dir, call. = FALSE)
    }
    summarise_s4b_low_exposure_counts_from_files(files)
}

## -----------------------------------------------------------------------------
## Validation
## -----------------------------------------------------------------------------
validate_s4b_low_exposure_data <- function(dat) {
    .require_object_s4b("validate_s3_data")
    validate_s3_data(dat)

    required <- c(
        "scenario_id", "stress_type", "stress_description",
        "exposure_stress_type", "e_reference", "exposure_multiplier",
        "target_mean_exposure_multiplier", "realized_mean_exposure_multiplier",
        "reference_count_summary", "mean_count", "zero_prop"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("S4B low-exposure dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (!identical(dat$stress_type, "low_exposure")) {
        stop("dat$stress_type must be 'low_exposure'.", call. = FALSE)
    }
    if (!identical(dat$exposure_stress_type, "low_heterogeneous_exposure")) {
        stop("dat$exposure_stress_type must be 'low_heterogeneous_exposure'.", call. = FALSE)
    }

    target_dim <- c(as.integer(dat$TT), as.integer(dat$n1))
    if (!.same_dim_s4b(dat$e_reference, target_dim) ||
        !.same_dim_s4b(dat$exposure_multiplier, target_dim)) {
        stop("e_reference and exposure_multiplier must have dimension TT x n1.", call. = FALSE)
    }

    if (any(!is.finite(dat$e)) || any(dat$e <= 0) ||
        any(!is.finite(dat$e_reference)) || any(dat$e_reference <= 0) ||
        any(!is.finite(dat$exposure_multiplier)) || any(dat$exposure_multiplier <= 0) ||
        any(dat$exposure_multiplier >= 1)) {
        stop("S4B exposure values or multipliers are non-finite or outside the required range.",
             call. = FALSE)
    }

    if (!isTRUE(all.equal(dat$e, dat$e_reference * dat$exposure_multiplier,
                          tolerance = 1e-8, check.attributes = FALSE))) {
        stop("dat$e is not equal to dat$e_reference * dat$exposure_multiplier.", call. = FALSE)
    }

    y_check <- apply(dat$y_fine, c(1, 2), sum)
    if (!isTRUE(all(dat$y_coarse == y_check))) {
        stop("S4B fine counts are not coherent with coarse counts.", call. = FALSE)
    }

    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Scenario 4B low-exposure simulation
## -----------------------------------------------------------------------------
simulate_s4b_low_exposure_one <- function(seed = 1L,
                                          TT_use = TT,
                                          target_mean_multiplier = 0.05,
                                          area_log_sd = 0.75,
                                          time_log_sd = 0.08,
                                          lower_multiplier = 0.005,
                                          upper_multiplier = 0.25,
                                          beta0_truth = -1.5,
                                          scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
                                          rep_id = NA_integer_,
                                          max_poisson_rate = 1e7,
                                          ...) {
    .require_object_s4b("simulate_s3_dynamic_learned_gamma_one")
    .require_object_s4b("compute_s3_xi")
    .require_object_s4b("recenter")

    ## Generate the reference Scenario 3 latent structure.
    ref_dat <- simulate_s3_dynamic_learned_gamma_one(
        seed = seed,
        TT_use = TT_use,
        beta0_truth = beta0_truth,
        scenario_id = paste0(scenario_id, "_REFERENCE_S3_LATENT"),
        rep_id = rep_id,
        max_poisson_rate = max_poisson_rate,
        ...
    )

    TT_now <- as.integer(ref_dat$TT)
    n1_now <- as.integer(ref_dat$n1)
    n_children_now <- as.integer(ref_dat$n_children)

    ## Create low and heterogeneous known exposure.
    mult_obj <- generate_s4b_exposure_multiplier(
        TT_use = TT_now,
        n1_use = n1_now,
        target_mean_multiplier = target_mean_multiplier,
        area_log_sd = area_log_sd,
        time_log_sd = time_log_sd,
        lower_multiplier = lower_multiplier,
        upper_multiplier = upper_multiplier,
        seed = as.integer(seed) + 888888L
    )
    exposure_multiplier <- mult_obj$multiplier
    e_low <- ref_dat$e * exposure_multiplier

    ## Recompute known-offset intensity with the same risk parameters and latent path.
    xi_low <- compute_s3_xi(
        e = e_low,
        x1 = ref_dat$x1,
        x2 = ref_dat$x2,
        beta0 = ref_dat$beta0_star,
        beta = ref_dat$beta_star,
        phi = ref_dat$phi_star
    )

    mu_nb_low <- xi_low * ref_dat$lambda_tilde
    poisson_rate_low <- mu_nb_low * ref_dat$kappa

    bad <- !is.finite(poisson_rate_low) | poisson_rate_low < 0 |
        poisson_rate_low > max_poisson_rate
    if (any(bad)) {
        idx <- which(bad, arr.ind = TRUE)[1, ]
        stop(sprintf(
            paste0(
                "Bad low-exposure Poisson rate in S4B generator. ",
                "First bad cell: t=%d, j=%d, rate=%s."
            ),
            idx[1], idx[2], as.character(poisson_rate_low[idx[1], idx[2]])
        ), call. = FALSE)
    }

    ## Use a separate deterministic seed for low-exposure observations.
    set.seed(as.integer(seed) + 999999L)
    y_coarse_low <- matrix(
        stats::rpois(TT_now * n1_now, lambda = as.numeric(poisson_rate_low)),
        nrow = TT_now,
        ncol = n1_now
    )

    y_fine_low <- array(0L, dim = c(TT_now, n1_now, n_children_now))
    for (t in seq_len(TT_now)) {
        for (j in seq_len(n1_now)) {
            if (y_coarse_low[t, j] > 0L) {
                y_fine_low[t, j, ] <- as.integer(stats::rmultinom(
                    1L,
                    size = y_coarse_low[t, j],
                    prob = ref_dat$omega[t, j, ]
                ))
            }
        }
    }

    coherent <- all(y_coarse_low == apply(y_fine_low, c(1, 2), sum))
    if (!coherent) stop("S4B low-exposure fine counts are not coherent.", call. = FALSE)

    ## Identified truth is unchanged except that the known exposure offset changes.
    rc_truth <- recenter(
        beta0 = ref_dat$beta0_star,
        phi = ref_dat$phi_star,
        lambda_tilde = ref_dat$lambda_tilde,
        return_diag = TRUE
    )

    dat <- ref_dat
    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_low_heterogeneous_exposure"
    dat$stress_type <- "low_exposure"
    dat$exposure_stress_type <- "low_heterogeneous_exposure"
    dat$stress_description <- paste0(
        "Scenario 3 latent structure with known exposure reduced and made heterogeneous; ",
        "this isolates low-exposure observation stress from beta0-shift sparse-count stress."
    )

    dat$y_coarse_reference <- ref_dat$y_coarse
    dat$e_reference <- ref_dat$e
    dat$exposure_multiplier <- exposure_multiplier
    dat$exposure_multiplier_area <- mult_obj$area_multiplier
    dat$target_mean_exposure_multiplier <- target_mean_multiplier
    dat$realized_mean_exposure_multiplier <- mult_obj$realized_mean_multiplier
    dat$realized_median_exposure_multiplier <- mult_obj$realized_median_multiplier
    dat$realized_min_exposure_multiplier <- mult_obj$realized_min_multiplier
    dat$realized_max_exposure_multiplier <- mult_obj$realized_max_multiplier
    dat$realized_cv_exposure_multiplier <- mult_obj$realized_cv_multiplier
    dat$area_multiplier_cv <- mult_obj$area_multiplier_cv
    dat$exposure_area_log_sd <- area_log_sd
    dat$exposure_time_log_sd <- time_log_sd
    dat$exposure_lower_multiplier <- lower_multiplier
    dat$exposure_upper_multiplier <- upper_multiplier

    dat$reference_count_summary <- list(
        mean_count = mean(ref_dat$y_coarse),
        median_count = stats::median(as.numeric(ref_dat$y_coarse)),
        zero_prop = mean(ref_dat$y_coarse == 0),
        total_count = sum(ref_dat$y_coarse),
        max_count = max(ref_dat$y_coarse),
        mean_exposure = mean(ref_dat$e),
        min_exposure = min(ref_dat$e),
        max_exposure = max(ref_dat$e)
    )

    dat$y_coarse <- y_coarse_low
    dat$y_fine <- y_fine_low
    dat$y_levels <- list(y_coarse_low, y_fine_low)
    dat$e <- e_low
    dat$xi_reference <- ref_dat$xi
    dat$mu_nb_reference <- ref_dat$mu_nb
    dat$poisson_rate_reference <- ref_dat$poisson_rate
    dat$xi <- xi_low
    dat$mu_nb <- mu_nb_low
    dat$poisson_rate <- poisson_rate_low

    dat$beta0_star <- ref_dat$beta0_star
    dat$beta_star <- ref_dat$beta_star
    dat$phi_star <- ref_dat$phi_star
    dat$beta0_star_ident <- rc_truth$beta0
    dat$beta_star_ident <- ref_dat$beta_star
    dat$phi_star_ident <- rc_truth$phi
    dat$lambda_tilde_ident <- rc_truth$lambda_tilde
    dat$lambda_recenter_diag <- rc_truth$diag

    dat$mean_count <- mean(y_coarse_low)
    dat$median_count <- stats::median(as.numeric(y_coarse_low))
    dat$zero_prop <- mean(y_coarse_low == 0)
    dat$nonzero_prop <- mean(y_coarse_low > 0)
    dat$total_count <- sum(y_coarse_low)
    dat$max_count <- max(y_coarse_low)
    dat$count_quantiles <- stats::quantile(
        as.numeric(y_coarse_low),
        probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1),
        names = TRUE
    )
    dat$mean_exposure <- mean(e_low)
    dat$median_exposure <- stats::median(as.numeric(e_low))
    dat$min_exposure <- min(e_low)
    dat$max_exposure <- max(e_low)
    dat$exposure_cv <- .cv_s4b(e_low)
    dat$exposure_group_summary <- exposure_group_summary_s4b(dat, n_groups = 3L)
    dat$coherent <- coherent

    validate_s4b_low_exposure_data(dat)
    dat
}

simulate_s4b_low_exposure_batch <- function(reps = 1:10,
                                             data_dir = "data_s4b_low_exposure",
                                             scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
                                             root = ".",
                                             seed_base = NULL,
                                             overwrite_existing = TRUE,
                                             verbose = TRUE,
                                             manifest_name = NULL,
                                             group_manifest_name = NULL,
                                             ...) {
    out_dir <- file.path(root, data_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    manifest <- list()
    group_manifest <- list()

    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        out_file <- file.path(out_dir, paste0("data_rep", rr, ".rds"))

        if (file.exists(out_file) && !isTRUE(overwrite_existing)) {
            if (isTRUE(verbose)) message("Skipping existing file: ", out_file)
            dat <- readRDS(out_file)
        } else {
            seed <- if (!is.null(seed_base)) {
                as.integer(seed_base + rep_id)
            } else if (exists("REP_SEEDS", envir = .GlobalEnv) && rep_id <= length(REP_SEEDS)) {
                as.integer(REP_SEEDS[rep_id])
            } else {
                as.integer(2026000L + rep_id)
            }

            dat <- simulate_s4b_low_exposure_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                message(sprintf(
                    paste0(
                        "Saved %s | mean_count=%.2f zero_prop=%.3f ",
                        "mean_e=%.2f mean_mult=%.4f exposure_cv=%.3f"
                    ),
                    out_file, dat$mean_count, dat$zero_prop,
                    dat$mean_exposure, dat$realized_mean_exposure_multiplier,
                    dat$exposure_cv
                ))
            }
        }

        row <- summarise_s4b_low_exposure_counts_one(dat)
        row$file <- out_file
        manifest[[rr]] <- row

        group_row <- exposure_group_summary_s4b(dat, n_groups = 3L)
        group_row$file <- out_file
        group_manifest[[rr]] <- group_row
    }

    manifest_df <- do.call(rbind, manifest)
    group_manifest_df <- do.call(rbind, group_manifest)

    if (is.null(manifest_name)) manifest_name <- paste0("manifest_", scenario_id, ".csv")
    if (is.null(group_manifest_name)) group_manifest_name <- paste0("manifest_exposure_groups_", scenario_id, ".csv")

    manifest_file <- file.path(out_dir, manifest_name)
    group_manifest_file <- file.path(out_dir, group_manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)
    write.csv(group_manifest_df, group_manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message("Saved exposure-group manifest: ", group_manifest_file)
        message(sprintf(
            paste0(
                "S4B low-exposure summary | reps=%d mean_count=%.2f zero_prop=%.3f ",
                "mean_multiplier=%.4f exposure_cv=%.3f"
            ),
            nrow(manifest_df),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            mean(manifest_df$realized_mean_multiplier),
            mean(manifest_df$exposure_cv)
        ))
    }

    invisible(manifest_df)
}

run_s4b_low_exposure_data_generation <- function(root = ".",
                                                 s3_script = NULL,
                                                 reps = 1:10,
                                                 TT_use = 100,
                                                 target_mean_multiplier = 0.05,
                                                 area_log_sd = 0.75,
                                                 time_log_sd = 0.08,
                                                 lower_multiplier = 0.005,
                                                 upper_multiplier = 0.25,
                                                 beta0_truth = -1.5,
                                                 data_dir = "data_s4b_low_exposure",
                                                 scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
                                                 seed_base = NULL,
                                                 overwrite_existing = TRUE,
                                                 verbose = TRUE,
                                                 ...) {
    source_s4b_low_exposure(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4b_low_exposure_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_use = TT_use,
        target_mean_multiplier = target_mean_multiplier,
        area_log_sd = area_log_sd,
        time_log_sd = time_log_sd,
        lower_multiplier = lower_multiplier,
        upper_multiplier = upper_multiplier,
        beta0_truth = beta0_truth,
        ...
    )
}

## -----------------------------------------------------------------------------
## Data-quality checks and calibration helpers
## -----------------------------------------------------------------------------
check_s4b_low_exposure_data_summary <- function(manifest,
                                                target_mean_count_range = c(5, 20),
                                                target_zero_prop_range = c(0.05, 0.35),
                                                target_mean_multiplier_range = c(0.03, 0.08),
                                                minimum_exposure_cv = 0.40,
                                                target_abs_beta0_ident_max = 20) {
    if (is.character(manifest)) manifest <- read.csv(manifest)

    required <- c(
        "rep_id", "TT", "n1", "target_mean_multiplier", "realized_mean_multiplier",
        "mean_count", "median_count", "zero_prop", "total_count", "max_count",
        "exposure_cv", "realized_cv_multiplier", "area_multiplier_cv",
        "beta0_star_ident", "lambda_raw_median", "lambda_raw_max",
        "lowest_exposure_group_mean_count", "highest_exposure_group_mean_count",
        "lowest_exposure_group_zero_prop", "highest_exposure_group_zero_prop"
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
        target_mean_multiplier_unique = paste(signif(sort(unique(manifest$target_mean_multiplier)), 6), collapse = ","),
        realized_mean_multiplier_avg = mean(manifest$realized_mean_multiplier),
        realized_mean_multiplier_min = min(manifest$realized_mean_multiplier),
        realized_mean_multiplier_max = max(manifest$realized_mean_multiplier),
        mean_count_avg = mean(manifest$mean_count),
        mean_count_min = min(manifest$mean_count),
        mean_count_max = max(manifest$mean_count),
        median_count_avg = mean(manifest$median_count),
        zero_prop_avg = mean(manifest$zero_prop),
        zero_prop_min = min(manifest$zero_prop),
        zero_prop_max = max(manifest$zero_prop),
        total_count_sum = sum(manifest$total_count),
        max_count_max = max(manifest$max_count),
        exposure_cv_avg = mean(manifest$exposure_cv),
        exposure_cv_min = min(manifest$exposure_cv),
        exposure_cv_max = max(manifest$exposure_cv),
        multiplier_cv_avg = mean(manifest$realized_cv_multiplier),
        area_multiplier_cv_avg = mean(manifest$area_multiplier_cv),
        low_exp_mean_count_avg = mean(manifest$lowest_exposure_group_mean_count),
        high_exp_mean_count_avg = mean(manifest$highest_exposure_group_mean_count),
        low_exp_zero_prop_avg = mean(manifest$lowest_exposure_group_zero_prop),
        high_exp_zero_prop_avg = mean(manifest$highest_exposure_group_zero_prop),
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
    out$passes_zero_prop_target <-
        out$zero_prop_avg >= target_zero_prop_range[1] &&
        out$zero_prop_avg <= target_zero_prop_range[2]
    out$passes_multiplier_target <-
        out$realized_mean_multiplier_avg >= target_mean_multiplier_range[1] &&
        out$realized_mean_multiplier_avg <= target_mean_multiplier_range[2]
    out$passes_exposure_heterogeneity_target <-
        out$exposure_cv_avg >= minimum_exposure_cv
    out$passes_identified_scale_target <-
        out$max_abs_beta0_ident <= target_abs_beta0_ident_max

    out$passes_s4b_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_multiplier_target) &&
        isTRUE(out$passes_exposure_heterogeneity_target) &&
        isTRUE(out$passes_identified_scale_target)

    out
}

calibrate_s4b_low_exposure_grid <- function(root = ".",
                                            s3_script = NULL,
                                            reps = 1:10,
                                            TT_use = 100,
                                            target_mean_multipliers = c(0.04, 0.05, 0.06, 0.08),
                                            area_log_sd = 0.75,
                                            time_log_sd = 0.08,
                                            lower_multiplier = 0.005,
                                            upper_multiplier = 0.25,
                                            beta0_truth = -1.5,
                                            data_dir = "data_s4b_low_exposure_calibration",
                                            base_scenario_id = "S4B_LOW_EXPOSURE_CALIB_T100",
                                            overwrite_existing = TRUE,
                                            verbose = TRUE,
                                            ...) {
    source_s4b_low_exposure(root = root, s3_script = s3_script, verbose = verbose)

    out <- list()
    for (mm in target_mean_multipliers) {
        tag <- gsub("\\.", "p", sprintf("%.3f", mm))
        scenario_id <- paste0(base_scenario_id, "_M", tag)
        if (isTRUE(verbose)) {
            message("\n--- S4B calibration: target_mean_multiplier = ", mm,
                    " | scenario_id = ", scenario_id, " ---")
        }

        manifest <- simulate_s4b_low_exposure_batch(
            reps = reps,
            data_dir = data_dir,
            scenario_id = scenario_id,
            root = root,
            overwrite_existing = overwrite_existing,
            verbose = verbose,
            TT_use = TT_use,
            target_mean_multiplier = mm,
            area_log_sd = area_log_sd,
            time_log_sd = time_log_sd,
            lower_multiplier = lower_multiplier,
            upper_multiplier = upper_multiplier,
            beta0_truth = beta0_truth,
            ...
        )
        chk <- check_s4b_low_exposure_data_summary(manifest)
        chk$scenario_id <- scenario_id
        chk$target_mean_multiplier <- mm
        out[[tag]] <- chk
    }

    calib <- do.call(rbind, out)
    out_dir <- file.path(root, data_dir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    calib_file <- file.path(out_dir, paste0("calibration_summary_", base_scenario_id, ".csv"))
    write.csv(calib, calib_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved S4B calibration summary: ", calib_file)

    invisible(calib)
}

## Example usage:
##
## source("s4b_low_exposure.R")
## manifest <- run_s4b_low_exposure_data_generation(
##     root = ".",
##     reps = 1:10,
##     TT_use = 100,
##     target_mean_multiplier = 0.05,
##     area_log_sd = 0.75,
##     time_log_sd = 0.08,
##     scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_T100",
##     overwrite_existing = TRUE
## )
## check_s4b_low_exposure_data_summary(manifest)
##
## Optional calibration grid before choosing official S4B:
## calib <- calibrate_s4b_low_exposure_grid(
##     root = ".",
##     reps = 1:10,
##     target_mean_multipliers = c(0.04, 0.05, 0.06, 0.08)
## )
