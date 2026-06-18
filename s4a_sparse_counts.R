## ============================================================================
## s4a_sparse_counts.R
## Scenario 4A core data-generation script for the MSSTNB project.
##
## Scenario 4A:
##   Sparse-count stress test corresponding to Scenario 3.
##
## Design rule:
##   This file is intentionally a thin wrapper around s3_dynamic_learned_gamma.R.
##   It keeps the Scenario 3 dynamic learned-gamma data-generating mechanism and
##   changes only the baseline log intensity by applying sparse_beta0_shift.
##
##   beta0_truth_s4a = beta0_reference_truth + sparse_beta0_shift
##
##   The calibrated default sparse_beta0_shift is -4.0, so the expected count scale is
##   multiplied by exp(-4.0) = 0.0183 before the dynamic latent risk and kappa
##   components are applied. In the first calibration over 10 replicates, this
##   produced mean_count about 4.24 and zero_prop about 0.335, making it a
##   meaningful sparse-count stress test rather than a mildly lower-count S3.
##
## Scope:
##   This script only generates and checks sparse-count data. It does not fit
##   the model. Fitting should be added only after the generated data are checked.
## ============================================================================

.same_dim_s4a <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.require_file_s4a <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked_s4a <- function(path, verbose = TRUE) {
    .require_file_s4a(path)
    if (isTRUE(verbose)) {
        message("source: ", path)
    }
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4a <- function(root = ".", s3_script = NULL) {
    if (!is.null(s3_script)) {
        return(s3_script)
    }

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

.require_object_s4a <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source Scenario 3 and project dependencies
## -----------------------------------------------------------------------------
source_s4a_sparse_counts <- function(root = ".",
                                     s3_script = NULL,
                                     verbose = TRUE) {
    s3_path <- .find_s3_script_s4a(root = root, s3_script = s3_script)
    .source_checked_s4a(s3_path, verbose = verbose)

    .require_object_s4a("source_s3_dynamic_learned_gamma")
    source_s3_dynamic_learned_gamma(root = root, verbose = verbose)

    needed <- c(
        "simulate_s3_dynamic_learned_gamma_one",
        "validate_s3_data",
        "REP_SEEDS",
        "TT", "N1", "N_CHILDREN"
    )
    missing <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop("After sourcing Scenario 3 dependencies, missing objects: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("Scenario 4A sparse-count data-generation wrapper loaded.")
    }

    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Validators and count summaries
## -----------------------------------------------------------------------------
validate_s4a_data <- function(dat) {
    .require_object_s4a("validate_s3_data")
    validate_s3_data(dat)

    required <- c(
        "scenario_id", "stress_type", "sparse_beta0_shift",
        "beta0_reference_truth", "beta0_sparse_truth",
        "expected_count_multiplier", "mean_count", "zero_prop"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("S4A dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (!identical(dat$stress_type, "sparse_counts")) {
        stop("dat$stress_type must be 'sparse_counts'.", call. = FALSE)
    }
    if (!is.finite(dat$sparse_beta0_shift) || dat$sparse_beta0_shift >= 0) {
        stop("dat$sparse_beta0_shift must be a finite negative value.", call. = FALSE)
    }
    if (!is.finite(dat$expected_count_multiplier) ||
        dat$expected_count_multiplier <= 0 ||
        dat$expected_count_multiplier >= 1) {
        stop("dat$expected_count_multiplier must be in (0, 1).", call. = FALSE)
    }

    invisible(TRUE)
}

summarise_s4a_counts_one <- function(dat) {
    validate_s4a_data(dat)

    y <- as.numeric(dat$y_coarse)
    data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        n1 = as.integer(dat$n1),
        sparse_beta0_shift = dat$sparse_beta0_shift,
        expected_count_multiplier = dat$expected_count_multiplier,
        beta0_reference_truth = dat$beta0_reference_truth,
        beta0_sparse_truth = dat$beta0_sparse_truth,
        beta0_star_ident = dat$beta0_star_ident,
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
        mean_exposure = mean(dat$e),
        min_exposure = min(dat$e),
        max_exposure = max(dat$e),
        gamma_truth_mean = mean(dat$gamma_star),
        delta_truth = dat$delta_star,
        r_truth_mean = mean(dat$r_star),
        lambda_raw_min = min(dat$lambda_tilde),
        lambda_raw_max = max(dat$lambda_tilde),
        lambda_ident_min = min(dat$lambda_tilde_ident),
        lambda_ident_max = max(dat$lambda_tilde_ident),
        lambda_ident_log_rm_mean = mean(abs(colMeans(log(dat$lambda_tilde_ident)))),
        coherent = isTRUE(dat$coherent),
        stringsAsFactors = FALSE
    )
}

summarise_s4a_counts_from_files <- function(files) {
    if (length(files) == 0L) {
        stop("No files supplied.", call. = FALSE)
    }
    out <- lapply(files, function(ff) {
        dat <- readRDS(ff)
        ss <- summarise_s4a_counts_one(dat)
        ss$file <- ff
        ss
    })
    do.call(rbind, out)
}

summarise_s4a_counts_from_dir <- function(root = ".",
                                          data_dir = "data_s4a_sparse_counts",
                                          scenario_id = "S4A_SPARSE_COUNTS_T100") {
    in_dir <- file.path(root, data_dir, scenario_id)
    files <- list.files(in_dir, pattern = "^data_rep[0-9]+\\.rds$",
                        full.names = TRUE)
    if (length(files) == 0L) {
        stop("No S4A data files found in: ", in_dir, call. = FALSE)
    }
    summarise_s4a_counts_from_files(files)
}

make_s4a_file_name <- function(prefix, rep_id, ext = "rds") {
    paste0(prefix, sprintf("%02d", as.integer(rep_id)), ".", ext)
}

## -----------------------------------------------------------------------------
## Scenario 4A simulation
## -----------------------------------------------------------------------------
simulate_s4a_sparse_counts_one <- function(seed = 1L,
                                           sparse_beta0_shift = -4.0,
                                           beta0_reference_truth = -1.5,
                                           scenario_id = "S4A_SPARSE_COUNTS_T100",
                                           rep_id = NA_integer_,
                                           ...) {
    .require_object_s4a("simulate_s3_dynamic_learned_gamma_one")

    if (!is.finite(sparse_beta0_shift) || sparse_beta0_shift >= 0) {
        stop("sparse_beta0_shift must be a finite negative value.", call. = FALSE)
    }
    if (!is.finite(beta0_reference_truth)) {
        stop("beta0_reference_truth must be finite.", call. = FALSE)
    }

    beta0_sparse_truth <- beta0_reference_truth + sparse_beta0_shift

    dat <- simulate_s3_dynamic_learned_gamma_one(
        seed = seed,
        beta0_truth = beta0_sparse_truth,
        scenario_id = scenario_id,
        rep_id = rep_id,
        ...
    )

    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_sparse_counts"
    dat$stress_type <- "sparse_counts"
    dat$stress_description <- "Scenario 3 dynamic learned-gamma DGP with lower baseline log intensity."
    dat$beta0_reference_truth <- beta0_reference_truth
    dat$sparse_beta0_shift <- sparse_beta0_shift
    dat$beta0_sparse_truth <- beta0_sparse_truth
    dat$expected_count_multiplier <- exp(sparse_beta0_shift)

    ## The Scenario 3 generator already stores beta0_star and beta0_star_ident
    ## using the shifted truth and recomputes the identified truth through recenter().
    dat$mean_count <- mean(dat$y_coarse)
    dat$median_count <- stats::median(as.numeric(dat$y_coarse))
    dat$zero_prop <- mean(dat$y_coarse == 0)
    dat$nonzero_prop <- mean(dat$y_coarse > 0)
    dat$total_count <- sum(dat$y_coarse)
    dat$max_count <- max(dat$y_coarse)
    dat$count_quantiles <- stats::quantile(as.numeric(dat$y_coarse),
                                           probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1),
                                           names = TRUE)

    validate_s4a_data(dat)
    dat
}

simulate_s4a_sparse_counts_batch <- function(reps = 1:10,
                                             data_dir = "data_s4a_sparse_counts",
                                             scenario_id = "S4A_SPARSE_COUNTS_T100",
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
            if (isTRUE(verbose)) {
                message("Skipping existing file: ", out_file)
            }
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

            dat <- simulate_s4a_sparse_counts_one(
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
                        "multiplier=%.3f beta0_sparse=%.3f"
                    ),
                    out_file, dat$mean_count, dat$zero_prop,
                    dat$expected_count_multiplier, dat$beta0_sparse_truth
                ))
            }
        }

        row <- summarise_s4a_counts_one(dat)
        row$file <- out_file
        manifest[[rr]] <- row
    }

    manifest_df <- do.call(rbind, manifest)

    if (is.null(manifest_name)) {
        manifest_name <- paste0("manifest_", scenario_id, ".csv")
    }
    manifest_file <- file.path(out_dir, manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message(sprintf(
            "S4A count summary | reps=%d mean_count=%.2f zero_prop=%.3f total_count=%d",
            nrow(manifest_df),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            as.integer(sum(manifest_df$total_count))
        ))
    }

    invisible(manifest_df)
}

## Convenience wrapper for the first data-generation pass.
run_s4a_sparse_counts_data_generation <- function(root = ".",
                                                  s3_script = NULL,
                                                  reps = 1:10,
                                                  TT_use = 100,
                                                  sparse_beta0_shift = -4.0,
                                                  beta0_reference_truth = -1.5,
                                                  data_dir = "data_s4a_sparse_counts",
                                                  scenario_id = "S4A_SPARSE_COUNTS_T100",
                                                  seed_base = NULL,
                                                  overwrite_existing = TRUE,
                                                  verbose = TRUE,
                                                  ...) {
    source_s4a_sparse_counts(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4a_sparse_counts_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_use = TT_use,
        sparse_beta0_shift = sparse_beta0_shift,
        beta0_reference_truth = beta0_reference_truth,
        ...
    )
}


## Official S4A data-level acceptance check.
check_s4a_official_data_summary <- function(manifest,
                                            target_mean_count_range = c(2, 8),
                                            target_zero_prop_range = c(0.20, 0.60)) {
    if (is.character(manifest)) {
        manifest <- read.csv(manifest)
    }

    required <- c(
        "rep_id", "TT", "n1", "sparse_beta0_shift",
        "expected_count_multiplier", "mean_count", "median_count",
        "zero_prop", "total_count", "max_count"
    )
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0L) {
        stop("Manifest is missing required columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    out <- data.frame(
        n_reps = nrow(manifest),
        TT_unique = paste(sort(unique(manifest$TT)), collapse = ","),
        n1_unique = paste(sort(unique(manifest$n1)), collapse = ","),
        sparse_beta0_shift_unique = paste(sort(unique(manifest$sparse_beta0_shift)), collapse = ","),
        expected_count_multiplier_unique = paste(signif(sort(unique(manifest$expected_count_multiplier)), 6), collapse = ","),
        mean_count_avg = mean(manifest$mean_count),
        mean_count_min = min(manifest$mean_count),
        mean_count_max = max(manifest$mean_count),
        median_count_avg = mean(manifest$median_count),
        zero_prop_avg = mean(manifest$zero_prop),
        zero_prop_min = min(manifest$zero_prop),
        zero_prop_max = max(manifest$zero_prop),
        total_count_sum = sum(manifest$total_count),
        max_count_max = max(manifest$max_count),
        stringsAsFactors = FALSE
    )

    out$passes_mean_count_target <-
        out$mean_count_avg >= target_mean_count_range[1] &&
        out$mean_count_avg <= target_mean_count_range[2]
    out$passes_zero_prop_target <-
        out$zero_prop_avg >= target_zero_prop_range[1] &&
        out$zero_prop_avg <= target_zero_prop_range[2]
    out$passes_s4a_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target)

    out
}

## Optional comparison helper. This is for data-level checks only.
compare_s4a_to_s3_count_manifests <- function(s4a_manifest,
                                              s3_manifest,
                                              by = "rep_id") {
    if (is.character(s4a_manifest)) {
        s4a_manifest <- read.csv(s4a_manifest)
    }
    if (is.character(s3_manifest)) {
        s3_manifest <- read.csv(s3_manifest)
    }

    keep_s4a <- c(by, "mean_count", "zero_prop", "total_count")
    keep_s3 <- c(by, "mean_count", "zero_prop")

    missing_s4a <- setdiff(keep_s4a, names(s4a_manifest))
    missing_s3 <- setdiff(keep_s3, names(s3_manifest))
    if (length(missing_s4a) > 0L || length(missing_s3) > 0L) {
        stop(
            "Manifest files do not contain the required count columns. Missing S4A: ",
            paste(missing_s4a, collapse = ", "),
            "; missing S3: ", paste(missing_s3, collapse = ", "),
            call. = FALSE
        )
    }

    s4a <- s4a_manifest[, keep_s4a, drop = FALSE]
    s3 <- s3_manifest[, keep_s3, drop = FALSE]
    names(s4a)[names(s4a) == "mean_count"] <- "mean_count_s4a"
    names(s4a)[names(s4a) == "zero_prop"] <- "zero_prop_s4a"
    names(s4a)[names(s4a) == "total_count"] <- "total_count_s4a"
    names(s3)[names(s3) == "mean_count"] <- "mean_count_s3"
    names(s3)[names(s3) == "zero_prop"] <- "zero_prop_s3"

    out <- merge(s4a, s3, by = by, all = FALSE)
    out$mean_count_ratio_s4a_to_s3 <- out$mean_count_s4a / out$mean_count_s3
    out$zero_prop_increase_s4a_minus_s3 <- out$zero_prop_s4a - out$zero_prop_s3
    out
}
