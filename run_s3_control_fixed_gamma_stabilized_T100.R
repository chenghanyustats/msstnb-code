## ============================================================================
## run_s3_control_fixed_gamma_stabilized_T100.R
##
## S3 control run using the same fixed-gamma stabilized machinery used for
## S4A sparse-count fixed-gamma diagnostics.
##
## Purpose
##   This is a validation/control script, not the main Scenario 3 analysis.
##   It fits the ordinary Scenario 3 dynamic T100 data with gamma fixed at truth
##   (0.8), while applying the same numerical guards used in the S4A fixed-gamma
##   stress-test fitting script.  The goal is to distinguish a genuine S4A
##   sparse-count failure mode from a general code/sampler problem.
##
## Interpretation
##   If this S3 control run has negligible guard counts and no beta--lambda
##   divergence, while S4A shows instability, that supports the interpretation
##   that the S4A failures are caused by sparse-count stress rather than by a
##   general implementation error.
##
## Required files in project root
##   s3_dynamic_learned_gamma.R
##   data_revised/DGP_DYNAMIC_T100/data_repXX.rds
##   or data_revised/S2_DYNAMIC_FIXED_GAMMA_T100/data_repXX.rds, which will be
##   copied into DGP_DYNAMIC_T100 if missing.
## ============================================================================

## ---- user settings ----------------------------------------------------------
root_dir <- "."

scenario_id <- "S3_CONTROL_FIXED_GAMMA_STABILIZED_T100"
source_data_scenario_id <- "DGP_DYNAMIC_T100"
legacy_source_data_scenario_id <- "S2_DYNAMIC_FIXED_GAMMA_T100"

## Recommended first run: control_3reps.  This uses the same formal MCMC length
## as S4A formal fits but only for reps 1:3.
## Options: "short_test", "control_3reps", "control_10reps", "control_20reps".
run_profile <- "control_3reps"

if (identical(run_profile, "short_test")) {
    reps_formal <- 1:1
    n_iter <- 6000L
    n_burnin <- 1000L
    n_thin <- 5L
} else if (identical(run_profile, "control_3reps")) {
    reps_formal <- 1:3
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
} else if (identical(run_profile, "control_10reps")) {
    reps_formal <- 1:10
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
} else if (identical(run_profile, "control_20reps")) {
    reps_formal <- 1:20
    n_iter <- 40000L
    n_burnin <- 20000L
    n_thin <- 5L
} else {
    stop("Unknown run_profile: ", run_profile, call. = FALSE)
}

s3_core_file <- "s3_dynamic_learned_gamma.R"
copy_legacy_s2_data_if_missing <- TRUE

TT_use <- 100L
n1_use <- 9L
x2_mode <- "continuous_time"
gamma_truth <- 0.80
fixed_gamma_value <- 0.80
gamma_prior <- c(1, 1)

## Output settings.
data_dir <- "data_revised"
output_dir <- "output_s3_control_fixed_gamma_stabilized"
verbose <- 1000L

## For a diagnostic control run, overwrite by default so the results definitely
## reflect this exact stabilized script.
overwrite_fit <- TRUE

## ---- helper functions -------------------------------------------------------
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

ensure_dir_s3ctl <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

assert_file_exists_s3ctl <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(TRUE)
}

assert_true_s3ctl <- function(x, message) {
    if (!isTRUE(x)) {
        stop(message, call. = FALSE)
    }
    invisible(TRUE)
}

s3ctl_data_file <- function(rep_id,
                            root = root_dir,
                            data_dir = data_dir,
                            data_scenario_id = source_data_scenario_id) {
    file.path(
        root,
        data_dir,
        data_scenario_id,
        sprintf("data_rep%02d.rds", as.integer(rep_id))
    )
}

copy_legacy_dynamic_data_if_needed_s3ctl <- function(reps,
                                                     root = root_dir,
                                                     data_dir = data_dir,
                                                     source_data_scenario_id = source_data_scenario_id,
                                                     legacy_source_data_scenario_id = legacy_source_data_scenario_id,
                                                     copy_if_missing = TRUE) {
    target_dir <- file.path(root, data_dir, source_data_scenario_id)
    target_files <- file.path(target_dir, sprintf("data_rep%02d.rds", as.integer(reps)))
    missing_target <- target_files[!file.exists(target_files)]

    if (length(missing_target) == 0L) {
        return(target_files)
    }

    if (!isTRUE(copy_if_missing)) {
        stop(
            "S3 control data files are missing. First missing file: ",
            missing_target[[1L]],
            call. = FALSE
        )
    }

    legacy_dir <- file.path(root, data_dir, legacy_source_data_scenario_id)
    assert_true_s3ctl(
        dir.exists(legacy_dir),
        paste0(
            "Shared DGP data are missing and legacy Scenario 2 data folder does not exist: ",
            legacy_dir,
            ". Run Scenario 2/S3 data generation first."
        )
    )

    ensure_dir_s3ctl(target_dir)
    legacy_files <- file.path(legacy_dir, sprintf("data_rep%02d.rds", as.integer(reps)))
    missing_legacy <- legacy_files[!file.exists(legacy_files)]
    assert_true_s3ctl(
        length(missing_legacy) == 0L,
        paste0("Some legacy dynamic data files are missing. First missing file: ", missing_legacy[[1L]])
    )

    message("Copying existing dynamic DGP data from ", legacy_dir, " to ", target_dir)
    copied <- file.copy(from = legacy_files, to = target_files, overwrite = FALSE)
    assert_true_s3ctl(
        all(copied | file.exists(target_files)),
        "Failed to copy at least one shared DGP data file."
    )

    target_files
}

check_s3ctl_source_dataset <- function(data_file,
                                       TT_use = 100L,
                                       n1_use = 9L,
                                       expected_x2_mode = "continuous_time") {
    assert_file_exists_s3ctl(data_file, "S3 control source dataset")
    dat <- readRDS(data_file)

    required <- c(
        "y_coarse", "e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident",
        "gamma_star", "beta_star", "phi_star", "TT", "n1"
    )
    missing <- required[!vapply(required, function(nm) !is.null(dat[[nm]]), logical(1L))]
    assert_true_s3ctl(
        length(missing) == 0L,
        paste("Dataset is missing required fields:", paste(missing, collapse = ", "))
    )

    assert_true_s3ctl(
        identical(as.integer(dim(dat$y_coarse)), c(as.integer(TT_use), as.integer(n1_use))),
        paste0("y_coarse dimension is not TT_use by n1_use. Got ", paste(dim(dat$y_coarse), collapse = " x "), ".")
    )

    for (nm in c("e", "x1", "x2", "lambda_tilde", "lambda_tilde_ident")) {
        assert_true_s3ctl(
            is.matrix(dat[[nm]]) && identical(dim(dat[[nm]]), dim(dat$y_coarse)),
            paste0("dat$", nm, " must have the same dimension as dat$y_coarse.")
        )
    }

    if (!is.null(dat$x2_mode)) {
        assert_true_s3ctl(
            identical(dat$x2_mode, expected_x2_mode),
            paste0("x2_mode is not ", expected_x2_mode, ". Got ", dat$x2_mode, ".")
        )
    }

    gamma_range <- range(dat$gamma_star, finite = TRUE)
    assert_true_s3ctl(
        all(abs(gamma_range - gamma_truth) < 1e-12),
        paste0("gamma_star is not fixed at gamma_truth. Range is ", paste(gamma_range, collapse = ", "), ".")
    )

    if (any(!is.finite(dat$y_coarse)) || any(dat$y_coarse < 0) ||
        any(abs(dat$y_coarse - round(dat$y_coarse)) > 1e-8)) {
        stop("dat$y_coarse must contain nonnegative integer counts.", call. = FALSE)
    }
    if (any(!is.finite(dat$e)) || any(dat$e <= 0)) {
        stop("dat$e must be positive and finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$x1)) || any(!is.finite(dat$x2))) {
        stop("dat$x1 and dat$x2 must be finite.", call. = FALSE)
    }
    if (any(!is.finite(dat$lambda_tilde)) || any(dat$lambda_tilde <= 0) ||
        any(!is.finite(dat$lambda_tilde_ident)) || any(dat$lambda_tilde_ident <= 0)) {
        stop("lambda_tilde and lambda_tilde_ident must be positive and finite.", call. = FALSE)
    }

    lambda_range <- range(dat$lambda_tilde, finite = TRUE)
    lambda_ident_range <- range(dat$lambda_tilde_ident, finite = TRUE)

    phi_true <- dat$phi_star_ident %||% dat$phi_star
    phi_mat <- matrix(phi_true, nrow = nrow(dat$x2), ncol = ncol(dat$x2), byrow = TRUE)
    x1 <- as.vector(dat$x1)
    x2 <- as.vector(dat$x2)
    loglam <- as.vector(log(dat$lambda_tilde_ident))

    list(
        dat = dat,
        scenario_id = dat$scenario_id %||% NA_character_,
        mean_count = mean(dat$y_coarse),
        median_count = stats::median(as.numeric(dat$y_coarse)),
        zero_prop = mean(dat$y_coarse == 0),
        total_count = sum(dat$y_coarse),
        max_count = max(dat$y_coarse),
        beta0_star_ident = dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_,
        lambda_raw_min = lambda_range[[1L]],
        lambda_raw_max = lambda_range[[2L]],
        lambda_ident_min = lambda_ident_range[[1L]],
        lambda_ident_max = lambda_ident_range[[2L]],
        sd_x1 = stats::sd(x1),
        sd_x2 = stats::sd(x2),
        cor_x1_x2 = suppressWarnings(stats::cor(x1, x2, use = "complete.obs")),
        cor_x2_loglambda = suppressWarnings(stats::cor(x2, loglam, use = "complete.obs")),
        cor_x2_phi = suppressWarnings(stats::cor(x2, as.vector(phi_mat), use = "complete.obs"))
    )
}

make_s3ctl_source_data_manifest <- function(reps,
                                            root = root_dir,
                                            data_dir = data_dir,
                                            data_scenario_id = source_data_scenario_id) {
    out <- lapply(reps, function(rep_id) {
        data_file <- s3ctl_data_file(rep_id, root = root, data_dir = data_dir, data_scenario_id = data_scenario_id)
        chk <- check_s3ctl_source_dataset(data_file)
        data.frame(
            scenario_id = data_scenario_id,
            rep_id = as.integer(rep_id),
            data_file = data_file,
            mean_count = chk$mean_count,
            median_count = chk$median_count,
            zero_prop = chk$zero_prop,
            total_count = chk$total_count,
            max_count = chk$max_count,
            beta0_star_ident = chk$beta0_star_ident,
            lambda_raw_min = chk$lambda_raw_min,
            lambda_raw_max = chk$lambda_raw_max,
            lambda_ident_min = chk$lambda_ident_min,
            lambda_ident_max = chk$lambda_ident_max,
            sd_x1 = chk$sd_x1,
            sd_x2 = chk$sd_x2,
            cor_x1_x2 = chk$cor_x1_x2,
            cor_x2_loglambda = chk$cor_x2_loglambda,
            cor_x2_phi = chk$cor_x2_phi,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

## ---- pre-flight: source Scenario 3 core ------------------------------------
assert_file_exists_s3ctl(file.path(root_dir, s3_core_file), "Scenario 3 core script")

message("S3 control fixed-gamma stabilized run profile: ", run_profile)
message("Replicates: ", paste(reps_formal, collapse = ", "))
message("MCMC: n_iter = ", n_iter, ", n_burnin = ", n_burnin, ", n_thin = ", n_thin)
message("overwrite_fit = ", overwrite_fit)

source_data_files <- copy_legacy_dynamic_data_if_needed_s3ctl(
    reps = reps_formal,
    root = root_dir,
    data_dir = data_dir,
    source_data_scenario_id = source_data_scenario_id,
    legacy_source_data_scenario_id = legacy_source_data_scenario_id,
    copy_if_missing = copy_legacy_s2_data_if_missing
)

source(file.path(root_dir, s3_core_file))
source_s3_dynamic_learned_gamma(root = root_dir)

## ---- same stabilized machinery used in S4A ---------------------------------
## Safe xi computation: work on log scale and clamp log(xi).  This should not
## affect ordinary posterior mass but prevents numerical overflow in extreme
## states.
s3ctl_log_xi_lower <- -40
s3ctl_log_xi_upper <-  40

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
        log_xi_j <- pmin(pmax(log_xi_j, s3ctl_log_xi_lower), s3ctl_log_xi_upper)
        xi[, j] <- exp(log_xi_j)
    }
    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("S3 control safe compute_s3_xi produced nonpositive or nonfinite xi.", call. = FALSE)
    }
    xi
}

cat(sprintf(
    "Using S3-control safe compute_s3_xi(): log(xi) clamped to [%.1f, %.1f].\n",
    s3ctl_log_xi_lower, s3ctl_log_xi_upper
))

.s3ctl_guard_env <- new.env(parent = emptyenv())
reset_s3ctl_numeric_guards <- function() {
    .s3ctl_guard_env$n_beta_guard <- 0L
    .s3ctl_guard_env$n_kappa_guard <- 0L
    .s3ctl_guard_env$n_lambda_input_guard <- 0L
    .s3ctl_guard_env$n_lambda_output_guard <- 0L
    invisible(TRUE)
}
get_s3ctl_numeric_guards <- function() {
    list(
        n_beta_guard = .s3ctl_guard_env$n_beta_guard %||% 0L,
        n_kappa_guard = .s3ctl_guard_env$n_kappa_guard %||% 0L,
        n_lambda_input_guard = .s3ctl_guard_env$n_lambda_input_guard %||% 0L,
        n_lambda_output_guard = .s3ctl_guard_env$n_lambda_output_guard %||% 0L
    )
}
reset_s3ctl_numeric_guards()

s3ctl_beta0_bounds <- c(-30, 10)
s3ctl_beta_bounds <- c(-5, 5)
s3ctl_kappa_bounds <- c(1e-10, 1e10)
s3ctl_lambda_bounds <- c(1e-10, 1e10)

.update_beta_s3_original <- update_beta
update_beta <- function(beta_current, ...) {
    res <- .update_beta_s3_original(beta_current = beta_current, ...)
    smp <- res$sample
    bad <- FALSE
    if (length(smp) < 3L || any(!is.finite(smp))) {
        bad <- TRUE
    } else {
        bad <- smp[1] < s3ctl_beta0_bounds[1] || smp[1] > s3ctl_beta0_bounds[2] ||
            any(smp[-1] < s3ctl_beta_bounds[1] | smp[-1] > s3ctl_beta_bounds[2])
    }
    if (isTRUE(bad)) {
        .s3ctl_guard_env$n_beta_guard <- (.s3ctl_guard_env$n_beta_guard %||% 0L) + 1L
        res$sample <- beta_current
        res$n_reject <- (res$n_reject %||% 0L) + 1L
        res$s3_control_guard_rejected <- TRUE
    } else {
        res$s3_control_guard_rejected <- FALSE
    }
    res
}

.update_kappa_s3_original <- update_kappa
update_kappa <- function(y_coarse, lambda_tilde, xi, r, return_diag = TRUE) {
    y <- as.matrix(y_coarse)
    L <- as.matrix(lambda_tilde)
    X <- as.matrix(xi)
    if (!identical(dim(y), dim(L)) || !identical(dim(y), dim(X))) {
        stop("S3 control safe update_kappa: y, lambda_tilde, and xi must have the same dimensions.", call. = FALSE)
    }
    r_vec <- as.numeric(r)
    if (length(r_vec) == 1L) r_vec <- rep(r_vec, ncol(y))
    if (length(r_vec) != ncol(y)) {
        stop("S3 control safe update_kappa: r must be scalar or length ncol(y).", call. = FALSE)
    }
    R <- matrix(rep(r_vec, each = nrow(y)), nrow = nrow(y), ncol = ncol(y))
    shape <- y + R
    rate <- X * L + R
    guard <- !is.finite(shape) | !is.finite(rate) | shape <= 0 | rate <= 0
    if (any(guard)) {
        .s3ctl_guard_env$n_kappa_guard <- (.s3ctl_guard_env$n_kappa_guard %||% 0L) + sum(guard)
    }
    shape <- pmin(pmax(shape, 1e-10), 1e10)
    rate <- pmin(pmax(rate, 1e-10), 1e10)
    kappa <- matrix(stats::rgamma(length(shape), shape = as.numeric(shape), rate = as.numeric(rate)),
                    nrow = nrow(y), ncol = ncol(y))
    bad_k <- !is.finite(kappa) | kappa <= 0
    if (any(bad_k)) {
        .s3ctl_guard_env$n_kappa_guard <- (.s3ctl_guard_env$n_kappa_guard %||% 0L) + sum(bad_k)
        kappa[bad_k] <- 1
    }
    kappa <- pmin(pmax(kappa, s3ctl_kappa_bounds[1]), s3ctl_kappa_bounds[2])
    diag <- list(
        mean_kappa = mean(kappa),
        min_kappa = min(kappa),
        max_kappa = max(kappa),
        n_guarded = .s3ctl_guard_env$n_kappa_guard %||% 0L
    )
    if (isTRUE(return_diag)) list(kappa = kappa, diag = diag) else kappa
}

.ffbs_lambda_all_s3_original <- ffbs_lambda_all
ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0, return_diag = TRUE) {
    X <- as.matrix(xi)
    K <- as.matrix(kappa)
    bad_in <- !is.finite(X) | X <= 0 | !is.finite(K) | K <= 0
    if (any(bad_in)) {
        .s3ctl_guard_env$n_lambda_input_guard <- (.s3ctl_guard_env$n_lambda_input_guard %||% 0L) + sum(bad_in)
    }
    X <- pmin(pmax(X, 1e-10), 1e10)
    K <- pmin(pmax(K, s3ctl_kappa_bounds[1]), s3ctl_kappa_bounds[2])
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
    bad_out <- !is.finite(L) | L <= 0 | L < s3ctl_lambda_bounds[1] | L > s3ctl_lambda_bounds[2]
    if (any(bad_out)) {
        .s3ctl_guard_env$n_lambda_output_guard <- (.s3ctl_guard_env$n_lambda_output_guard %||% 0L) + sum(bad_out)
        L <- pmin(pmax(L, s3ctl_lambda_bounds[1]), s3ctl_lambda_bounds[2])
        out$lambda_tilde <- L
        if (!is.null(out$diag)) {
            out$diag$min_lambda <- min(L)
            out$diag$max_lambda <- max(L)
            out$diag$s3_control_lambda_output_guarded <- TRUE
        }
    }
    out
}

cat(sprintf(
    "Using S3-control stabilization guards: beta0 in [%.1f, %.1f], beta in [%.1f, %.1f], log(xi) in [%.1f, %.1f].\n",
    s3ctl_beta0_bounds[1], s3ctl_beta0_bounds[2],
    s3ctl_beta_bounds[1], s3ctl_beta_bounds[2],
    s3ctl_log_xi_lower, s3ctl_log_xi_upper
))

## Fixed-gamma override.  This intentionally disables learned-gamma updates so
## the S3 control matches the S4A main fixed-gamma machinery.
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
cat(sprintf("Using S3-control fixed-gamma override: gamma fixed at %.3f.\n", fixed_gamma_value))

## ---- fitting functions ------------------------------------------------------
fit_s3_control_fixed_gamma_stabilized_one_rep <- function(rep_id,
                                                          scenario_id = "S3_CONTROL_FIXED_GAMMA_STABILIZED_T100",
                                                          data_scenario_id = "DGP_DYNAMIC_T100",
                                                          data_dir = "data_revised",
                                                          output_dir = "output_s3_control_fixed_gamma_stabilized",
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
    chk <- check_s3ctl_source_dataset(dat_file)
    dat <- chk$dat

    validate_s3_data(dat)

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    gamma_start <- gamma_init %||% fixed_gamma_value
    if (length(gamma_start) == 1L) gamma_start <- rep(gamma_start, dat$n1)

    cat(sprintf("=== S3 control fixed-gamma stabilized fit: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Counts    : mean = %.3f, zero proportion = %.3f\n", chk$mean_count, chk$zero_prop))
    cat(sprintf("Fixed     : gamma = %.3f\n", mean(gamma_start)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r, kappa, lambda\n")
    cat("Disabled  : gamma, delta, omega updates\n\n")

    reset_s3ctl_numeric_guards()

    fit <- run_s3_dynamic_learned_gamma_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        gamma_init = gamma_start,
        gamma_prior = gamma_prior,
        verbose = verbose
    )

    guard_counts <- get_s3ctl_numeric_guards()
    fit$diagnostics$s3_control_numeric_guards <- guard_counts
    fit$diagnostics$s3_control_beta_guard_count <- guard_counts$n_beta_guard
    fit$diagnostics$s3_control_kappa_guard_count <- guard_counts$n_kappa_guard
    fit$diagnostics$s3_control_lambda_input_guard_count <- guard_counts$n_lambda_input_guard
    fit$diagnostics$s3_control_lambda_output_guard_count <- guard_counts$n_lambda_output_guard

    fit$diagnostics$gamma_accept_rate <- NA_real_
    fit$diagnostics$gamma_proposal_sd_final <- NA_real_
    fit$diagnostics$gamma_sd <- stats::sd(fit$samples$gamma_common, na.rm = TRUE)
    fit$diagnostics$gamma_mean <- mean(fit$samples$gamma_common, na.rm = TRUE)
    fit$metadata$model <- "S3 control dynamic NB-ICAR with fixed gamma and S4A-style stabilization"
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

    ## Data-level context and guard counts.
    summary$source_data_scenario_id <- data_scenario_id
    summary$observed_mean_count <- chk$mean_count
    summary$observed_zero_prop <- chk$zero_prop
    summary$observed_total_count <- chk$total_count
    summary$observed_max_count <- chk$max_count
    summary$s3_control_beta_guard_count <- fit$diagnostics$s3_control_beta_guard_count
    summary$s3_control_kappa_guard_count <- fit$diagnostics$s3_control_kappa_guard_count
    summary$s3_control_lambda_input_guard_count <- fit$diagnostics$s3_control_lambda_input_guard_count
    summary$s3_control_lambda_output_guard_count <- fit$diagnostics$s3_control_lambda_output_guard_count
    summary$numeric_beta_guard_count <- fit$diagnostics$s3_control_beta_guard_count
    summary$numeric_kappa_guard_count <- fit$diagnostics$s3_control_kappa_guard_count
    summary$numeric_lambda_input_guard_count <- fit$diagnostics$s3_control_lambda_input_guard_count
    summary$numeric_lambda_output_guard_count <- fit$diagnostics$s3_control_lambda_output_guard_count

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            source_data_scenario_id = data_scenario_id,
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            fit_file_prefix = "fit_S3_control_fixed_gamma_stabilized_rep",
            gamma_fixed_in_fit = TRUE,
            gamma_learned_in_fit = FALSE,
            gamma_fixed_value = fixed_gamma_value,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        ensure_dir_s3ctl(out_dir)
        fit_file <- file.path(out_dir, paste0("fit_S3_control_fixed_gamma_stabilized_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S3_control_fixed_gamma_stabilized_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        utils::write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) return(fit)
    invisible(NULL)
}

fit_s3_control_fixed_gamma_stabilized_batch <- function(reps = 1:3,
                                                        scenario_id = "S3_CONTROL_FIXED_GAMMA_STABILIZED_T100",
                                                        data_scenario_id = "DGP_DYNAMIC_T100",
                                                        data_dir = "data_revised",
                                                        output_dir = "output_s3_control_fixed_gamma_stabilized",
                                                        root = ".",
                                                        settings_override = list(),
                                                        gamma_init = NULL,
                                                        fixed_gamma_value = 0.8,
                                                        gamma_prior = c(1, 1),
                                                        verbose = 1000L,
                                                        overwrite_existing = TRUE) {
    out_dir <- file.path(root, output_dir, scenario_id)
    ensure_dir_s3ctl(out_dir)

    summaries <- list()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        fit_file <- file.path(out_dir, paste0("fit_S3_control_fixed_gamma_stabilized_rep", rr, ".rds"))
        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            message("Skipping existing fit: ", fit_file)
            fit <- readRDS(fit_file)
            summaries[[rr]] <- fit$summary
            next
        }

        fit <- fit_s3_control_fixed_gamma_stabilized_one_rep(
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
    summary_file <- file.path(out_dir, "summary_S3_control_fixed_gamma_stabilized_all_reps.csv")
    utils::write.csv(summary_all, summary_file, row.names = FALSE)
    message("Saved combined summary: ", summary_file)
    invisible(summary_all)
}

print_s3ctl_fit_driver_summary <- function(source_manifest, fit_summary, fit_manifest) {
    cat("\n=== S3 control fixed-gamma stabilized fitting summary ===\n")
    cat("Number of reps: ", nrow(source_manifest), "\n", sep = "")
    cat("Mean count avg: ", round(mean(source_manifest$mean_count, na.rm = TRUE), 4), "\n", sep = "")
    cat("Zero prop avg : ", round(mean(source_manifest$zero_prop, na.rm = TRUE), 4), "\n", sep = "")
    cat("Fit files all present: ", all(fit_manifest$fit_exists), "\n", sep = "")

    if (!is.null(fit_summary) && nrow(fit_summary) > 0L) {
        keep <- intersect(
            c(
                "scenario_id", "rep_id", "mean_count", "zero_prop",
                "beta0_mean", "beta1_mean", "beta2_mean", "r_mean",
                "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
                "s3_control_beta_guard_count", "s3_control_kappa_guard_count",
                "s3_control_lambda_input_guard_count", "s3_control_lambda_output_guard_count"
            ),
            names(fit_summary)
        )
        cat("\nSelected fit-summary columns:\n")
        print(fit_summary[, keep, drop = FALSE])
    }
    invisible(TRUE)
}

## ---- verify source data and run control fits --------------------------------
source_data_manifest <- make_s3ctl_source_data_manifest(
    reps = reps_formal,
    root = root_dir,
    data_dir = data_dir,
    data_scenario_id = source_data_scenario_id
)

cat("\n=== S3 control source-data check ===\n")
print(source_data_manifest[, c(
    "rep_id", "mean_count", "median_count", "zero_prop", "total_count",
    "max_count", "beta0_star_ident", "lambda_raw_min", "lambda_raw_max"
), drop = FALSE])

fit_summary <- fit_s3_control_fixed_gamma_stabilized_batch(
    reps = reps_formal,
    scenario_id = scenario_id,
    data_scenario_id = source_data_scenario_id,
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
    sprintf("fit_S3_control_fixed_gamma_stabilized_rep%02d.rds", as.integer(reps_formal))
)
fit_manifest <- data.frame(
    scenario_id = scenario_id,
    source_data_scenario_id = source_data_scenario_id,
    rep_id = as.integer(reps_formal),
    data_file = source_data_manifest$data_file,
    fit_file = fit_files,
    fit_exists = file.exists(fit_files),
    stringsAsFactors = FALSE
)
assert_true_s3ctl(all(fit_manifest$fit_exists), "At least one S3 control fit file was not created.")

ensure_dir_s3ctl(fit_dir_full)
utils::write.csv(
    source_data_manifest,
    file.path(fit_dir_full, "s3_control_source_data_manifest.csv"),
    row.names = FALSE
)
utils::write.csv(
    fit_manifest,
    file.path(fit_dir_full, "s3_control_fit_manifest.csv"),
    row.names = FALSE
)
saveRDS(
    list(
        run_profile = run_profile,
        scenario_id = scenario_id,
        source_data_scenario_id = source_data_scenario_id,
        reps = reps_formal,
        mcmc = list(n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin),
        fixed_gamma_value = fixed_gamma_value,
        gamma_prior = gamma_prior,
        source_data_manifest = source_data_manifest,
        fit_manifest = fit_manifest,
        fit_summary = fit_summary
    ),
    file.path(fit_dir_full, "run_info_S3_control_fixed_gamma_stabilized_T100.rds")
)

print_s3ctl_fit_driver_summary(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
)

cat("\n=== Main output locations ===\n")
cat("S3 control data: ", file.path(data_dir, source_data_scenario_id), "\n", sep = "")
cat("Fits           : ", file.path(output_dir, scenario_id), "\n", sep = "")
cat("\nS3 control fixed-gamma stabilized T100 fitting finished successfully.\n")

invisible(list(
    source_manifest = source_data_manifest,
    fit_summary = fit_summary,
    fit_manifest = fit_manifest
))
