## ============================================================================
## clean_AB_old_sampler_ablation.R
##
## Clean A-vs-B diagnostic for the OLD kappa-conditioned MSSTNB sampler.
##
## A = old kappa-conditioned sampler with gamma fixed at truth/default.
## B = old kappa-conditioned sampler with gamma learned by update_gamma().
##
## Intended use: source this file from the MSSTNB project root, after making sure
## the old sampler files are available in R/mcmc/ or another directory supplied
## by old_mcmc_dir.
## ============================================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

.first_nonnull <- function(...) {
    xs <- list(...)
    for (x in xs) if (!is.null(x)) return(x)
    NULL
}

.safe_mean <- function(x) {
    x <- as.numeric(x)
    if (!length(x) || all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.safe_cor <- function(x, y) {
    x <- as.numeric(x); y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    as.numeric(cor(x[ok], y[ok]))
}

.rmse <- function(est, truth) {
    est <- as.numeric(est); truth <- as.numeric(truth)
    ok <- is.finite(est) & is.finite(truth)
    if (!any(ok)) return(NA_real_)
    sqrt(mean((est[ok] - truth[ok])^2))
}

.recenter_truth <- function(beta0, phi, lambda_tilde) {
    if (is.null(beta0) || is.null(phi) || is.null(lambda_tilde)) {
        return(list(beta0 = beta0, phi = phi, lambda_tilde = lambda_tilde))
    }
    s <- colMeans(log(lambda_tilde))
    s_bar <- mean(s)
    lambda_new <- lambda_tilde
    for (j in seq_len(ncol(lambda_new))) {
        lambda_new[, j] <- lambda_new[, j] * exp(-s[j])
    }
    list(
        beta0 = as.numeric(beta0) + s_bar,
        phi = as.numeric(phi) + s - s_bar,
        lambda_tilde = lambda_new
    )
}

## ---- robust truth extraction ------------------------------------------------
extract_truth_clean_AB <- function(dat, default_gamma = 0.8, identify_truth = TRUE) {
    beta0 <- .first_nonnull(
        dat$beta0_star_ident, dat$beta0_ident, dat$beta0_star,
        dat$beta0_true, dat$truth$beta0_star_ident, dat$truth$beta0_true,
        dat$truth$beta0
    )

    beta <- .first_nonnull(
        dat$beta_star_ident, dat$beta_ident, dat$beta_star,
        dat$beta_true, dat$truth$beta_star_ident, dat$truth$beta_true,
        dat$truth$beta
    )

    ## Sometimes beta_true may include the intercept.
    if (!is.null(beta) && length(beta) >= 3L && is.null(beta0)) {
        beta0 <- beta[1]
        beta <- beta[-1]
    }

    phi <- .first_nonnull(
        dat$phi_star_ident, dat$phi_ident, dat$phi_star,
        dat$phi_true, dat$truth$phi_star_ident, dat$truth$phi_true,
        dat$truth$phi
    )

    lambda_tilde <- .first_nonnull(
        dat$lambda_tilde_ident, dat$lambda_tilde_star_ident,
        dat$lambda_tilde_star, dat$lambda_tilde_true,
        dat$truth$lambda_tilde_ident, dat$truth$lambda_tilde_true,
        dat$truth$lambda_tilde
    )

    ## If only raw truth is present and identified truth is not, optionally apply
    ## the same recentering convention used by the sampler.
    has_identified_fields <- !is.null(dat$beta0_star_ident) ||
        !is.null(dat$phi_star_ident) || !is.null(dat$lambda_tilde_ident) ||
        !is.null(dat$truth$beta0_star_ident) || !is.null(dat$truth$phi_star_ident)

    if (identify_truth && !has_identified_fields &&
        !is.null(beta0) && !is.null(phi) && !is.null(lambda_tilde)) {
        rc <- .recenter_truth(beta0, phi, lambda_tilde)
        beta0 <- rc$beta0
        phi <- rc$phi
        lambda_tilde <- rc$lambda_tilde
    }

    gamma <- .first_nonnull(
        dat$gamma_true, dat$gamma_star, dat$gamma,
        dat$truth$gamma_true, dat$truth$gamma_star, dat$truth$gamma
    )
    if (is.null(gamma)) gamma <- default_gamma
    if (length(gamma) == 1L && !is.null(dat$n1)) gamma <- rep(gamma, dat$n1)

    if (is.null(beta0) || length(beta0) == 0L) beta0 <- NA_real_
    if (is.null(beta)  || length(beta)  == 0L) beta  <- c(NA_real_, NA_real_)
    if (length(beta) == 1L) beta <- c(beta, NA_real_)
    if (is.null(phi) || length(phi) == 0L) {
        phi <- if (!is.null(dat$n1)) rep(NA_real_, dat$n1) else NA_real_
    }

    list(
        beta0 = as.numeric(beta0)[1],
        beta = as.numeric(beta),
        phi = as.numeric(phi),
        lambda_tilde = lambda_tilde,
        gamma = as.numeric(gamma)
    )
}

## ---- source old sampler files ---------------------------------------------
source_old_sampler_for_AB <- function(project_dir = ".", old_mcmc_dir = NULL,
                                      source_project_setup = TRUE) {
    project_dir <- normalizePath(project_dir, mustWork = TRUE)
    old_mcmc_dir <- old_mcmc_dir %||% file.path(project_dir, "R", "mcmc")
    old_mcmc_dir <- normalizePath(old_mcmc_dir, mustWork = TRUE)

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(project_dir)

    if (isTRUE(source_project_setup)) {
        helper_file <- file.path(project_dir, "R", "01_helpers.R")
        setup_file  <- file.path(project_dir, "R", "00_setup.R")
        if (file.exists(helper_file)) source(helper_file)
        if (file.exists(setup_file)) source(setup_file)
    }

    needed <- c(
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
    for (f in needed) {
        path <- file.path(old_mcmc_dir, f)
        if (!file.exists(path)) stop("Cannot find old sampler file: ", path)
        source(path)
    }

    invisible(TRUE)
}

build_spatial_from_globals_AB <- function() {
    if (!exists("H", envir = .GlobalEnv)) stop("Object H not found. Source R/00_setup.R first.")
    if (!exists("B_ICAR", envir = .GlobalEnv)) stop("Object B_ICAR not found. Source R/00_setup.R first.")
    H_obj <- get("H", envir = .GlobalEnv)
    B_obj <- get("B_ICAR", envir = .GlobalEnv)
    BHB <- crossprod(B_obj, H_obj %*% B_obj)
    list(H = H_obj, B_ICAR = B_obj, R_BHB = chol(BHB))
}

build_constants_from_globals_AB <- function() {
    list(
        A0 = get("A0", envir = .GlobalEnv),
        B0 = get("B0", envir = .GlobalEnv),
        C0 = get("C0", envir = .GlobalEnv)
    )
}

## ---- one iteration with optional fixed gamma -------------------------------
run_one_iteration_old_AB <- function(state, dat, settings, priors, spatial,
                                     constants, fixed_gamma = NULL) {
    diag <- list()

    if (settings$include_nb) {
        state$kappa <- update_kappa(dat$y_coarse, state$lambda_tilde,
                                    state$xi, state$r)
    }

    if (settings$include_covariates) {
        beta_result <- update_beta(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            kappa = state$kappa, lambda_tilde = state$lambda_tilde,
            phi = state$phi, priors = priors
        )
        state$beta0 <- beta_result$sample[1]
        state$beta  <- beta_result$sample[2:3]
        diag$beta_ess_reject <- beta_result$n_reject
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    if (settings$include_icar) {
        phi_result <- update_phi(
            u_current = state$u, B = spatial$B_ICAR,
            R_BHB = spatial$R_BHB, tau_phi = state$tau_phi,
            y_coarse = dat$y_coarse, e = dat$e,
            x1 = dat$x1, x2 = dat$x2,
            beta0 = state$beta0, beta = state$beta,
            kappa = state$kappa, lambda_tilde = state$lambda_tilde
        )
        state$u   <- phi_result$u
        state$phi <- phi_result$phi
        diag$phi_ess_reject <- phi_result$n_reject
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    if (settings$include_icar) {
        state$tau_phi <- update_tau_phi(state$phi, spatial$H,
                                        dat$n1, priors)
    }

    if (settings$include_nb) {
        r_result <- update_r(state$r, state$kappa, priors,
                             mh_sd = settings$mh_sd_log_r)
        state$r <- r_result$r
        diag$r_accept <- r_result$accept
    }

    if (!is.null(fixed_gamma)) {
        state$gamma <- fixed_gamma
        diag$gamma_accept <- rep(FALSE, length(state$gamma))
    } else {
        gamma_result <- update_gamma(
            state$gamma, dat$y_coarse, state$xi, state$kappa,
            a0 = constants$A0, b0 = constants$B0,
            priors = priors, mh_sd = settings$mh_sd_gamma
        )
        state$gamma <- gamma_result$gamma
        diag$gamma_accept <- gamma_result$accept
    }

    state$lambda_tilde <- ffbs_lambda_all(
        state$gamma, dat$y_coarse, state$xi, state$kappa,
        a0 = constants$A0, b0 = constants$B0
    )

    rc <- recenter(state$beta0, state$phi, state$lambda_tilde)
    state$beta0        <- rc$beta0
    state$phi          <- rc$phi
    state$lambda_tilde <- rc$lambda_tilde

    if (settings$include_icar) {
        state$u <- as.numeric(t(spatial$B_ICAR) %*% state$phi)
    }

    state$xi <- compute_xi_mcmc(dat$e, dat$x1, dat$x2,
                                state$beta0, state$beta, state$phi)

    delta_result <- update_delta(
        state$delta, dat$y_fine, c0 = constants$C0,
        priors = priors, mh_sd = settings$mh_sd_delta
    )
    state$delta <- delta_result$delta
    diag$delta_accept <- delta_result$accept

    state$omega <- smooth_omega_all(state$delta, dat$y_fine,
                                    c0 = constants$C0)

    diag$loglik <- compute_loglik(dat$y_coarse, state$xi,
                                  state$lambda_tilde, state$kappa)

    list(state = state, diag = diag)
}

run_mcmc_old_AB <- function(dat, settings, priors, spatial, constants,
                            fixed_gamma = NULL, verbose = 1000L) {
    n_iter   <- settings$n_iter
    n_burnin <- settings$n_burnin
    n_thin   <- settings$n_thin
    n_stored <- as.integer((n_iter - n_burnin) / n_thin)

    TT <- dat$TT
    n1 <- dat$n1
    K  <- dat$n_children
    p_beta <- length(dat$beta_star %||% dat$beta_true %||% c(NA_real_, NA_real_))
    if (p_beta >= 3L) p_beta <- p_beta - 1L
    if (!is.finite(p_beta) || p_beta < 1L) p_beta <- 2L

    samples <- list(
        beta0        = numeric(n_stored),
        beta         = matrix(NA_real_, n_stored, p_beta),
        phi          = matrix(NA_real_, n_stored, n1),
        tau_phi      = numeric(n_stored),
        r            = matrix(NA_real_, n_stored, n1),
        gamma        = matrix(NA_real_, n_stored, n1),
        delta        = numeric(n_stored),
        lambda_tilde = array(NA_real_, dim = c(n_stored, TT, n1)),
        kappa        = array(NA_real_, dim = c(n_stored, TT, n1)),
        omega        = array(NA_real_, dim = c(n_stored, TT, n1, K)),
        loglik       = numeric(n_stored)
    )

    diag_all <- list(
        loglik_trace       = numeric(n_iter),
        beta_ess_reject    = numeric(n_iter),
        phi_ess_reject     = numeric(n_iter),
        r_accept_rate      = numeric(n1),
        gamma_accept_rate  = numeric(n1),
        delta_accept_count = 0L,
        fixed_gamma        = fixed_gamma
    )

    state <- initialise_state(dat, settings, priors, spatial, constants)
    if (!is.null(fixed_gamma)) state$gamma <- fixed_gamma

    start_time <- proc.time()
    store_idx <- 0L

    for (iter in seq_len(n_iter)) {
        result <- run_one_iteration_old_AB(
            state, dat, settings, priors, spatial, constants,
            fixed_gamma = fixed_gamma
        )
        state <- result$state
        diag  <- result$diag

        diag_all$loglik_trace[iter] <- diag$loglik
        if (!is.null(diag$beta_ess_reject)) diag_all$beta_ess_reject[iter] <- diag$beta_ess_reject
        if (!is.null(diag$phi_ess_reject)) diag_all$phi_ess_reject[iter] <- diag$phi_ess_reject
        if (!is.null(diag$r_accept)) diag_all$r_accept_rate <- diag_all$r_accept_rate + diag$r_accept
        if (!is.null(diag$gamma_accept)) diag_all$gamma_accept_rate <- diag_all$gamma_accept_rate + diag$gamma_accept
        if (!is.null(diag$delta_accept) && isTRUE(diag$delta_accept)) {
            diag_all$delta_accept_count <- diag_all$delta_accept_count + 1L
        }

        if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
            store_idx <- store_idx + 1L
            samples$beta0[store_idx]             <- state$beta0
            samples$beta[store_idx, ]            <- state$beta
            samples$phi[store_idx, ]             <- state$phi
            samples$tau_phi[store_idx]           <- state$tau_phi
            samples$r[store_idx, ]               <- state$r
            samples$gamma[store_idx, ]           <- state$gamma
            samples$delta[store_idx]             <- state$delta
            samples$lambda_tilde[store_idx, , ]  <- state$lambda_tilde
            samples$kappa[store_idx, , ]         <- state$kappa
            samples$omega[store_idx, , , ]       <- state$omega
            samples$loglik[store_idx]            <- diag$loglik
        }

        if (verbose > 0 && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            cat(sprintf("  iter %5d/%d [%.0fs] loglik=%.1f beta0=%.3f gamma[1]=%.3f\n",
                        iter, n_iter, elapsed, diag$loglik,
                        state$beta0, state$gamma[1]))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]
    diag_all$r_accept_rate     <- diag_all$r_accept_rate / n_iter
    diag_all$gamma_accept_rate <- diag_all$gamma_accept_rate / n_iter
    diag_all$delta_accept_rate <- diag_all$delta_accept_count / n_iter
    diag_all$elapsed_sec       <- elapsed_total

    if (verbose > 0) {
        cat(sprintf("Done. %d iterations in %.1f sec; stored %d draws.\n",
                    n_iter, elapsed_total, n_stored))
        cat(sprintf("Acceptance rates: gamma=%.3f delta=%.3f r=%.3f\n",
                    mean(diag_all$gamma_accept_rate),
                    diag_all$delta_accept_rate,
                    mean(diag_all$r_accept_rate)))
    }

    list(
        samples = samples,
        diagnostics = diag_all,
        settings = settings,
        priors = priors,
        method = "old_kappa_conditioned",
        n_stored = n_stored
    )
}

## ---- run A/B fits ----------------------------------------------------------
fit_old_AB_one <- function(data_file, variant = c("A_old_fixed_gamma", "B_old_learned_gamma"),
                           out_file = NULL, n_iter = 10000L, n_burnin = 5000L,
                           n_thin = 5L, gamma_fixed = NULL,
                           default_gamma = 0.8, verbose = 1000L) {
    variant <- match.arg(variant)
    dat <- readRDS(data_file)

    settings <- MCMC_SETTINGS
    settings$n_iter   <- as.integer(n_iter)
    settings$n_burnin <- as.integer(n_burnin)
    settings$n_thin   <- as.integer(n_thin)
    settings$include_nb <- TRUE
    settings$include_icar <- TRUE
    settings$include_covariates <- TRUE

    priors <- MCMC_PRIORS
    spatial <- build_spatial_from_globals_AB()
    constants <- build_constants_from_globals_AB()

    fixed_gamma_vec <- NULL
    if (variant == "A_old_fixed_gamma") {
        truth <- extract_truth_clean_AB(dat, default_gamma = default_gamma)
        fixed_gamma_vec <- gamma_fixed %||% truth$gamma
        if (length(fixed_gamma_vec) == 1L) fixed_gamma_vec <- rep(fixed_gamma_vec, dat$n1)
        if (length(fixed_gamma_vec) != dat$n1) {
            stop("fixed gamma must have length 1 or n1. Got length ", length(fixed_gamma_vec))
        }
    }

    cat("\n=== Running ", variant, " ===\n", sep = "")
    cat("Data: ", data_file, "\n", sep = "")
    if (!is.null(fixed_gamma_vec)) {
        cat("Fixed gamma mean: ", sprintf("%.4f", mean(fixed_gamma_vec)), "\n", sep = "")
    }

    fit <- run_mcmc_old_AB(
        dat = dat, settings = settings, priors = priors,
        spatial = spatial, constants = constants,
        fixed_gamma = fixed_gamma_vec,
        verbose = verbose
    )
    fit$variant <- variant
    fit$data_file <- data_file
    fit$fixed_gamma <- fixed_gamma_vec

    if (!is.null(out_file)) {
        dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
        saveRDS(fit, out_file)
        cat("Saved: ", out_file, "\n", sep = "")
    }
    fit
}

run_clean_AB_old_sampler <- function(project_dir = ".",
                                     old_mcmc_dir = NULL,
                                     scenario_id = "S3",
                                     reps = 1:3,
                                     data_dir = "data_revised",
                                     output_dir = "analysis_sampler_ablation_clean/fits_AB",
                                     n_iter = 10000L,
                                     n_burnin = 5000L,
                                     n_thin = 5L,
                                     default_gamma = 0.8,
                                     overwrite = FALSE,
                                     verbose = 1000L) {
    source_old_sampler_for_AB(project_dir = project_dir, old_mcmc_dir = old_mcmc_dir)

    manifest <- data.frame()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        data_file <- file.path(project_dir, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
        if (!file.exists(data_file)) stop("Data file not found: ", data_file)

        for (variant in c("A_old_fixed_gamma", "B_old_learned_gamma")) {
            out_file <- file.path(project_dir, output_dir, scenario_id,
                                  paste0(variant, "_rep", rr, ".rds"))
            status <- "skipped_existing"
            if (overwrite || !file.exists(out_file)) {
                fit_old_AB_one(
                    data_file = data_file, variant = variant, out_file = out_file,
                    n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
                    default_gamma = default_gamma, verbose = verbose
                )
                status <- "completed"
            }
            manifest <- rbind(manifest, data.frame(
                scenario_id = scenario_id,
                rep = as.integer(rep_id),
                variant = variant,
                data_file = data_file,
                fit_file = out_file,
                status = status,
                stringsAsFactors = FALSE
            ))
        }
    }

    manifest_file <- file.path(project_dir, output_dir, scenario_id, "AB_manifest.csv")
    dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(manifest, manifest_file, row.names = FALSE)
    cat("\nManifest saved: ", manifest_file, "\n", sep = "")
    invisible(manifest)
}

## ---- summarize A/B fits ----------------------------------------------------
summarize_fit_old_AB <- function(fit_file, data_file, default_gamma = 0.8) {
    fit <- readRDS(fit_file)
    dat <- readRDS(data_file)
    truth <- extract_truth_clean_AB(dat, default_gamma = default_gamma)
    s <- fit$samples

    beta0_mean <- .safe_mean(s$beta0)
    beta_mean <- if (!is.null(s$beta)) colMeans(s$beta, na.rm = TRUE) else c(NA_real_, NA_real_)
    beta_truth_full <- c(truth$beta0, truth$beta)
    beta_est_full <- c(beta0_mean, beta_mean)
    beta_rmse <- .rmse(beta_est_full, beta_truth_full)

    phi_mean <- if (!is.null(s$phi)) colMeans(s$phi, na.rm = TRUE) else rep(NA_real_, dat$n1)
    phi_rmse <- .rmse(phi_mean, truth$phi)
    phi_cor  <- .safe_cor(phi_mean, truth$phi)

    log_lambda_rmse <- NA_real_
    log_lambda_cor <- NA_real_
    delta_log_lambda_cor <- NA_real_
    if (!is.null(s$lambda_tilde) && !is.null(truth$lambda_tilde)) {
        log_lambda_draws <- log(s$lambda_tilde)
        log_lambda_mean <- apply(log_lambda_draws, c(2, 3), mean, na.rm = TRUE)
        log_lambda_truth <- log(truth$lambda_tilde)
        log_lambda_rmse <- .rmse(log_lambda_mean, log_lambda_truth)
        log_lambda_cor <- .safe_cor(log_lambda_mean, log_lambda_truth)
        if (nrow(log_lambda_mean) >= 2L) {
            delta_log_lambda_cor <- .safe_cor(
                diff(log_lambda_mean), diff(log_lambda_truth)
            )
        }
    }

    gamma_mean_vec <- if (!is.null(s$gamma)) colMeans(s$gamma, na.rm = TRUE) else NA_real_
    gamma_truth_vec <- truth$gamma
    if (length(gamma_truth_vec) == 1L && length(gamma_mean_vec) > 1L) {
        gamma_truth_vec <- rep(gamma_truth_vec, length(gamma_mean_vec))
    }

    data.frame(
        variant = fit$variant %||% NA_character_,
        fit_file = fit_file,
        data_file = data_file,
        beta0_mean = beta0_mean,
        beta1_mean = beta_mean[1] %||% NA_real_,
        beta2_mean = beta_mean[2] %||% NA_real_,
        beta0_truth = truth$beta0 %||% NA_real_,
        beta1_truth = truth$beta[1] %||% NA_real_,
        beta2_truth = truth$beta[2] %||% NA_real_,
        beta0_bias = beta0_mean - (truth$beta0 %||% NA_real_),
        beta1_bias = (beta_mean[1] %||% NA_real_) - (truth$beta[1] %||% NA_real_),
        beta2_bias = (beta_mean[2] %||% NA_real_) - (truth$beta[2] %||% NA_real_),
        beta_rmse = beta_rmse,
        phi_rmse = phi_rmse,
        phi_cor = phi_cor,
        log_lambda_rmse = log_lambda_rmse,
        log_lambda_cor = log_lambda_cor,
        delta_log_lambda_cor = delta_log_lambda_cor,
        gamma_mean = .safe_mean(gamma_mean_vec),
        gamma_truth = .safe_mean(gamma_truth_vec),
        gamma_bias = .safe_mean(gamma_mean_vec - gamma_truth_vec),
        gamma_accept_mean = .safe_mean(fit$diagnostics$gamma_accept_rate),
        r_accept_mean = .safe_mean(fit$diagnostics$r_accept_rate),
        delta_accept = fit$diagnostics$delta_accept_rate %||% NA_real_,
        elapsed_sec = fit$diagnostics$elapsed_sec %||% NA_real_,
        max_log_lambda_mean = if (!is.null(s$lambda_tilde)) max(apply(log(s$lambda_tilde), c(2, 3), mean, na.rm = TRUE), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
    )
}

summarize_clean_AB_old_sampler <- function(project_dir = ".",
                                           scenario_id = "S3",
                                           output_dir = "analysis_sampler_ablation_clean/fits_AB",
                                           default_gamma = 0.8) {
    manifest_file <- file.path(project_dir, output_dir, scenario_id, "AB_manifest.csv")
    if (!file.exists(manifest_file)) stop("Manifest not found: ", manifest_file)
    manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
    rows <- list()
    for (i in seq_len(nrow(manifest))) {
        if (!file.exists(manifest$fit_file[i])) next
        rows[[length(rows) + 1L]] <- summarize_fit_old_AB(
            fit_file = manifest$fit_file[i],
            data_file = manifest$data_file[i],
            default_gamma = default_gamma
        )
        rows[[length(rows)]]$scenario_id <- manifest$scenario_id[i]
        rows[[length(rows)]]$rep <- manifest$rep[i]
    }
    summary <- do.call(rbind, rows)
    summary <- summary[order(summary$rep, summary$variant), ]

    table_dir <- file.path(project_dir, "analysis_sampler_ablation_clean", "tables")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    out_summary <- file.path(table_dir, paste0("AB_old_sampler_summary_", scenario_id, ".csv"))
    write.csv(summary, out_summary, row.names = FALSE)
    cat("Summary saved: ", out_summary, "\n", sep = "")

    metrics <- c("beta_rmse", "phi_rmse", "phi_cor", "log_lambda_rmse",
                 "log_lambda_cor", "delta_log_lambda_cor", "gamma_mean", "gamma_bias")
    contrast_rows <- list()
    reps <- sort(unique(summary$rep))
    for (m in metrics) {
        for (rr in reps) {
            a <- summary[summary$rep == rr & summary$variant == "A_old_fixed_gamma", m]
            b <- summary[summary$rep == rr & summary$variant == "B_old_learned_gamma", m]
            if (length(a) == 1L && length(b) == 1L) {
                contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
                    rep = rr, metric = m, A = a, B = b, B_minus_A = b - a
                )
            }
        }
    }
    contrasts_long <- if (length(contrast_rows)) do.call(rbind, contrast_rows) else data.frame()
    contrast_summary <- data.frame()
    if (nrow(contrasts_long)) {
        for (m in unique(contrasts_long$metric)) {
            d <- contrasts_long$B_minus_A[contrasts_long$metric == m]
            contrast_summary <- rbind(contrast_summary, data.frame(
                metric = m,
                n_pairs = sum(is.finite(d)),
                mean_B_minus_A = mean(d, na.rm = TRUE),
                sd_B_minus_A = sd(d, na.rm = TRUE),
                mcse_B_minus_A = sd(d, na.rm = TRUE) / sqrt(sum(is.finite(d))),
                stringsAsFactors = FALSE
            ))
        }
    }

    out_contrasts_long <- file.path(table_dir, paste0("AB_old_sampler_paired_contrasts_long_", scenario_id, ".csv"))
    out_contrasts <- file.path(table_dir, paste0("AB_old_sampler_paired_contrasts_", scenario_id, ".csv"))
    write.csv(contrasts_long, out_contrasts_long, row.names = FALSE)
    write.csv(contrast_summary, out_contrasts, row.names = FALSE)
    cat("Paired contrasts saved: ", out_contrasts, "\n", sep = "")

    invisible(list(summary = summary, contrasts_long = contrasts_long,
                   contrast_summary = contrast_summary))
}

## ---- example ----------------------------------------------------------------
## From project root:
## source("clean_AB_old_sampler_ablation.R")
## run_clean_AB_old_sampler(
##     project_dir = ".",
##     old_mcmc_dir = "R/mcmc",       # or a folder containing the old sampler files
##     scenario_id = "S3",
##     reps = 1:3,
##     data_dir = "data_revised",
##     n_iter = 10000,
##     n_burnin = 5000,
##     n_thin = 5,
##     overwrite = TRUE,
##     verbose = 1000
## )
## out <- summarize_clean_AB_old_sampler(project_dir = ".", scenario_id = "S3")
