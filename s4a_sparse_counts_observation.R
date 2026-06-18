## ============================================================================
## s4a_sparse_counts_observation.R
## Scenario 4A clean sparse-count data-generation script for the MSSTNB project.
##
## Scenario 4A clean version:
##   Observation-level sparse-count stress test corresponding to Scenario 3.
##
## Design rule:
##   1. Generate the Scenario 3 latent structure under the reference beta0.
##   2. Keep the same exposure, covariates, phi, lambda_tilde, omega, kappa,
##      gamma, delta, and r.
##   3. Re-generate only the observation counts under a lower baseline intensity:
##
##        beta0_sparse_truth = beta0_reference_truth + sparse_beta0_shift.
##
##   This isolates sparse observed counts from dynamic-lambda path collapse.
##   It is intentionally different from lowering beta0 inside the full sequential
##   S3 generator, where smaller y also feeds back into the gamma-discount
##   filtering shapes and can create pathological lambda scales.
##
## Scope:
##   This script only generates and checks S4A sparse-count data. It does not fit
##   the model.
## ============================================================================

.same_dim_s4a_obs <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.require_file_s4a_obs <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.source_checked_s4a_obs <- function(path, verbose = TRUE) {
    .require_file_s4a_obs(path)
    if (isTRUE(verbose)) {
        message("source: ", path)
    }
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4a_obs <- function(root = ".", s3_script = NULL) {
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

.require_object_s4a_obs <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Source Scenario 3 and project dependencies
## -----------------------------------------------------------------------------
source_s4a_sparse_counts_observation <- function(root = ".",
                                                 s3_script = NULL,
                                                 verbose = TRUE) {
    s3_path <- .find_s3_script_s4a_obs(root = root, s3_script = s3_script)
    .source_checked_s4a_obs(s3_path, verbose = verbose)

    .require_object_s4a_obs("source_s3_dynamic_learned_gamma")
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
        message("Scenario 4A clean observation-level sparse-count generator loaded.")
    }

    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Count and scale summaries
## -----------------------------------------------------------------------------
summarise_s4a_observation_counts_one <- function(dat) {
    validate_s4a_observation_data(dat)

    y <- as.numeric(dat$y_coarse)
    data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        n1 = as.integer(dat$n1),
        stress_type = dat$stress_type,
        sparse_beta0_shift = dat$sparse_beta0_shift,
        expected_count_multiplier = dat$expected_count_multiplier,
        beta0_reference_truth = dat$beta0_reference_truth,
        beta0_sparse_truth = dat$beta0_sparse_truth,
        beta0_star = dat$beta0_star,
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

summarise_s4a_observation_counts_from_files <- function(files) {
    if (length(files) == 0L) {
        stop("No files supplied.", call. = FALSE)
    }
    out <- lapply(files, function(ff) {
        dat <- readRDS(ff)
        ss <- summarise_s4a_observation_counts_one(dat)
        ss$file <- ff
        ss
    })
    do.call(rbind, out)
}

summarise_s4a_observation_counts_from_dir <- function(root = ".",
                                                      data_dir = "data_s4a_sparse_counts",
                                                      scenario_id = "S4A_SPARSE_COUNTS_T100") {
    in_dir <- file.path(root, data_dir, scenario_id)
    files <- list.files(in_dir, pattern = "^data_rep[0-9]+\\.rds$",
                        full.names = TRUE)
    if (length(files) == 0L) {
        stop("No S4A data files found in: ", in_dir, call. = FALSE)
    }
    summarise_s4a_observation_counts_from_files(files)
}

validate_s4a_observation_data <- function(dat) {
    .require_object_s4a_obs("validate_s3_data")
    validate_s3_data(dat)

    required <- c(
        "scenario_id", "stress_type", "stress_description",
        "sparse_beta0_shift", "beta0_reference_truth", "beta0_sparse_truth",
        "expected_count_multiplier", "mean_count", "zero_prop",
        "reference_count_summary"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
        stop("S4A observation dat is missing fields: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    if (!identical(dat$stress_type, "observation_sparse_counts")) {
        stop("dat$stress_type must be 'observation_sparse_counts'.", call. = FALSE)
    }
    if (!is.finite(dat$sparse_beta0_shift) || dat$sparse_beta0_shift >= 0) {
        stop("dat$sparse_beta0_shift must be a finite negative value.", call. = FALSE)
    }
    if (!is.finite(dat$expected_count_multiplier) ||
        dat$expected_count_multiplier <= 0 ||
        dat$expected_count_multiplier >= 1) {
        stop("dat$expected_count_multiplier must be in (0, 1).", call. = FALSE)
    }

    y_check <- apply(dat$y_fine, c(1, 2), sum)
    if (!isTRUE(all(dat$y_coarse == y_check))) {
        stop("S4A fine counts are not coherent with coarse counts.", call. = FALSE)
    }

    invisible(TRUE)
}

## -----------------------------------------------------------------------------
## Scenario 4A clean observation-level sparse simulation
## -----------------------------------------------------------------------------
simulate_s4a_observation_sparse_one <- function(seed = 1L,
                                                TT_use = TT,
                                                sparse_beta0_shift = -4.0,
                                                beta0_reference_truth = -1.5,
                                                scenario_id = "S4A_SPARSE_COUNTS_T100",
                                                rep_id = NA_integer_,
                                                max_poisson_rate = 1e7,
                                                ...) {
    .require_object_s4a_obs("simulate_s3_dynamic_learned_gamma_one")
    .require_object_s4a_obs("compute_s3_xi")
    .require_object_s4a_obs("recenter")

    if (!is.finite(sparse_beta0_shift) || sparse_beta0_shift >= 0) {
        stop("sparse_beta0_shift must be a finite negative value.", call. = FALSE)
    }
    if (!is.finite(beta0_reference_truth)) {
        stop("beta0_reference_truth must be finite.", call. = FALSE)
    }

    beta0_sparse_truth <- beta0_reference_truth + sparse_beta0_shift

    ## Generate the full Scenario 3 latent structure under the reference intensity.
    ref_dat <- simulate_s3_dynamic_learned_gamma_one(
        seed = seed,
        TT_use = TT_use,
        beta0_truth = beta0_reference_truth,
        scenario_id = paste0(scenario_id, "_REFERENCE_S3_LATENT"),
        rep_id = rep_id,
        max_poisson_rate = max_poisson_rate,
        ...
    )

    TT_now <- as.integer(ref_dat$TT)
    n1_now <- as.integer(ref_dat$n1)
    n_children_now <- as.integer(ref_dat$n_children)

    ## New observation intensity with the same latent lambda and kappa paths.
    xi_sparse <- compute_s3_xi(
        e = ref_dat$e,
        x1 = ref_dat$x1,
        x2 = ref_dat$x2,
        beta0 = beta0_sparse_truth,
        beta = ref_dat$beta_star,
        phi = ref_dat$phi_star
    )

    mu_nb_sparse <- xi_sparse * ref_dat$lambda_tilde
    poisson_rate_sparse <- mu_nb_sparse * ref_dat$kappa

    bad <- !is.finite(poisson_rate_sparse) | poisson_rate_sparse < 0 |
        poisson_rate_sparse > max_poisson_rate
    if (any(bad)) {
        idx <- which(bad, arr.ind = TRUE)[1, ]
        stop(sprintf(
            paste0(
                "Bad sparse Poisson rate in S4A observation generator. ",
                "First bad cell: t=%d, j=%d, rate=%s."
            ),
            idx[1], idx[2], as.character(poisson_rate_sparse[idx[1], idx[2]])
        ), call. = FALSE)
    }

    ## Use a separate deterministic seed for the sparse observations so the
    ## reference latent path remains tied to the Scenario 3 replicate seed.
    set.seed(as.integer(seed) + 777777L)
    y_coarse_sparse <- matrix(
        stats::rpois(TT_now * n1_now, lambda = as.numeric(poisson_rate_sparse)),
        nrow = TT_now,
        ncol = n1_now
    )

    y_fine_sparse <- array(0L, dim = c(TT_now, n1_now, n_children_now))
    for (t in seq_len(TT_now)) {
        for (j in seq_len(n1_now)) {
            if (y_coarse_sparse[t, j] > 0L) {
                y_fine_sparse[t, j, ] <- as.integer(stats::rmultinom(
                    1L,
                    size = y_coarse_sparse[t, j],
                    prob = ref_dat$omega[t, j, ]
                ))
            }
        }
    }

    coherent <- all(y_coarse_sparse == apply(y_fine_sparse, c(1, 2), sum))
    if (!coherent) {
        stop("S4A observation sparse fine counts are not coherent.", call. = FALSE)
    }

    ## Identified truth on the sparse beta0 scale, with the stable S3 lambda path.
    rc_truth <- recenter(
        beta0 = beta0_sparse_truth,
        phi = ref_dat$phi_star,
        lambda_tilde = ref_dat$lambda_tilde,
        return_diag = TRUE
    )

    dat <- ref_dat
    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_observation_sparse_counts"
    dat$stress_type <- "observation_sparse_counts"
    dat$stress_description <- paste0(
        "Scenario 3 latent structure with observations regenerated under lower baseline ",
        "log intensity; this isolates sparse observed counts from lambda-path collapse."
    )

    dat$y_coarse_reference <- ref_dat$y_coarse
    dat$reference_count_summary <- list(
        mean_count = mean(ref_dat$y_coarse),
        median_count = stats::median(as.numeric(ref_dat$y_coarse)),
        zero_prop = mean(ref_dat$y_coarse == 0),
        total_count = sum(ref_dat$y_coarse),
        max_count = max(ref_dat$y_coarse)
    )

    dat$y_coarse <- y_coarse_sparse
    dat$y_fine <- y_fine_sparse
    dat$y_levels <- list(y_coarse_sparse, y_fine_sparse)
    dat$xi_reference <- ref_dat$xi
    dat$mu_nb_reference <- ref_dat$mu_nb
    dat$poisson_rate_reference <- ref_dat$poisson_rate
    dat$xi <- xi_sparse
    dat$mu_nb <- mu_nb_sparse
    dat$poisson_rate <- poisson_rate_sparse

    dat$beta0_reference_truth <- beta0_reference_truth
    dat$sparse_beta0_shift <- sparse_beta0_shift
    dat$beta0_sparse_truth <- beta0_sparse_truth
    dat$expected_count_multiplier <- exp(sparse_beta0_shift)

    dat$beta0_star <- beta0_sparse_truth
    dat$beta_star <- ref_dat$beta_star
    dat$phi_star <- ref_dat$phi_star
    dat$beta0_star_ident <- rc_truth$beta0
    dat$beta_star_ident <- ref_dat$beta_star
    dat$phi_star_ident <- rc_truth$phi
    dat$lambda_tilde_ident <- rc_truth$lambda_tilde
    dat$lambda_recenter_diag <- rc_truth$diag

    dat$mean_count <- mean(y_coarse_sparse)
    dat$median_count <- stats::median(as.numeric(y_coarse_sparse))
    dat$zero_prop <- mean(y_coarse_sparse == 0)
    dat$nonzero_prop <- mean(y_coarse_sparse > 0)
    dat$total_count <- sum(y_coarse_sparse)
    dat$max_count <- max(y_coarse_sparse)
    dat$count_quantiles <- stats::quantile(
        as.numeric(y_coarse_sparse),
        probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1),
        names = TRUE
    )
    dat$coherent <- coherent

    validate_s4a_observation_data(dat)
    dat
}

simulate_s4a_observation_sparse_batch <- function(reps = 1:10,
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

            dat <- simulate_s4a_observation_sparse_one(
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
                        "multiplier=%.5f beta0_sparse=%.3f beta0_ident=%.3f"
                    ),
                    out_file, dat$mean_count, dat$zero_prop,
                    dat$expected_count_multiplier, dat$beta0_sparse_truth,
                    dat$beta0_star_ident
                ))
            }
        }

        row <- summarise_s4a_observation_counts_one(dat)
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
            "S4A observation sparse summary | reps=%d mean_count=%.2f zero_prop=%.3f beta0_ident_median=%.3f",
            nrow(manifest_df),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            stats::median(manifest_df$beta0_star_ident)
        ))
    }

    invisible(manifest_df)
}

run_s4a_observation_sparse_data_generation <- function(root = ".",
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
    source_s4a_sparse_counts_observation(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4a_observation_sparse_batch(
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

check_s4a_observation_data_summary <- function(manifest,
                                               target_mean_count_range = c(2, 8),
                                               target_zero_prop_range = c(0.20, 0.60),
                                               target_abs_beta0_ident_max = 20) {
    if (is.character(manifest)) {
        manifest <- read.csv(manifest)
    }

    required <- c(
        "rep_id", "TT", "n1", "sparse_beta0_shift",
        "expected_count_multiplier", "mean_count", "median_count",
        "zero_prop", "total_count", "max_count", "beta0_star_ident",
        "lambda_raw_median", "lambda_raw_max"
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
    out$passes_identified_scale_target <-
        out$max_abs_beta0_ident <= target_abs_beta0_ident_max
    out$passes_s4a_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_identified_scale_target)

    out
}
