## ============================================================================
## s2_oracle_lambda_diagnostic.R
##
## Oracle-lambda diagnostic for Scenario 2: dynamic residual risk with fixed gamma.
##
## Purpose:
##   Use the same Scenario 2 data, but fix lambda_tilde at the true identified
##   lambda path during fitting.  Then update only beta, phi, tau_phi, and r.
##   This diagnostic separates two possibilities:
##
##   1. If beta2 recovers under oracle lambda, then the sampled dynamic lambda
##      block is absorbing part of the x2 signal.
##   2. If beta2 still does not recover under oracle lambda, then the issue is
##      more likely an x2/DGP mismatch or a beta-update implementation problem.
##
## Dependencies:
##   This script reuses helper functions from
##   s2_dynamic_fixed_gamma_v2_with_regression_dispersion_patch1.R, especially:
##       source_s2_dynamic_fixed_gamma()
##       validate_s2_data()
##       build_s2_spatial()
##       validate_s2_fit_objects()
##       initialise_s2_state()
##       update_beta_s2()
##       update_phi_s2()
##       update_r_s2()
##       adapt_sd_s2()
##       compute_s2_xi()
##       compute_s2_loglik_nb()
##       make_s2_file_name()
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.require_file_s2_oracle <- function(path) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

source_s2_oracle_lambda_diagnostic <- function(
    root = ".",
    s2_script = "s2_dynamic_fixed_gamma_v2_with_regression_dispersion_patch1.R",
    verbose = TRUE
) {
    ## Source the Scenario 2 local-updater script if needed.  The file may be in
    ## the current working directory or in root.
    if (!exists("source_s2_dynamic_fixed_gamma", envir = .GlobalEnv)) {
        candidates <- c(s2_script, file.path(root, s2_script))
        candidates <- unique(candidates)
        hit <- candidates[file.exists(candidates)][1]
        if (is.na(hit)) {
            stop(
                "Could not find the Scenario 2 support script. Please place ",
                s2_script, " in the working directory or project root, or source it before this diagnostic.",
                call. = FALSE
            )
        }
        if (isTRUE(verbose)) {
            message("source oracle support: ", hit)
        }
        source(hit, local = .GlobalEnv)
    }

    source_s2_dynamic_fixed_gamma(root = root, verbose = verbose)

    needed <- c(
        "MCMC_SETTINGS", "MCMC_PRIORS",
        "validate_s2_data", "build_s2_spatial", "validate_s2_fit_objects",
        "initialise_s2_state", "update_beta_s2", "update_phi_s2",
        "update_r_s2", "update_tau_phi", "adapt_sd_s2",
        "compute_s2_xi", "compute_s2_loglik_nb", "make_s2_file_name"
    )
    missing <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop(
            "After sourcing Scenario 2 dependencies, missing objects: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    if (isTRUE(verbose)) {
        message("Scenario 2 oracle-lambda diagnostic loaded.")
    }
    invisible(TRUE)
}

get_s2_oracle_lambda_truth <- function(dat, lambda_source = c("identified", "raw")) {
    lambda_source <- match.arg(lambda_source)
    if (lambda_source == "identified") {
        lambda_fixed <- dat$lambda_tilde_ident %||% dat$lambda_tilde
    } else {
        lambda_fixed <- dat$lambda_tilde
    }
    if (is.null(lambda_fixed)) {
        stop("No lambda_tilde or lambda_tilde_ident found in dat.", call. = FALSE)
    }
    if (!is.matrix(lambda_fixed) || !identical(dim(lambda_fixed), dim(dat$y_coarse))) {
        stop("Oracle lambda must be a matrix with the same dimension as dat$y_coarse.", call. = FALSE)
    }
    if (any(!is.finite(lambda_fixed)) || any(lambda_fixed <= 0)) {
        stop("Oracle lambda contains nonpositive or nonfinite values.", call. = FALSE)
    }
    lambda_fixed
}

run_s2_oracle_lambda_mcmc <- function(
    dat,
    settings = MCMC_SETTINGS,
    priors = MCMC_PRIORS,
    spatial = build_s2_spatial(),
    lambda_fixed = get_s2_oracle_lambda_truth(dat, lambda_source = "identified"),
    gamma_fixed = dat$gamma_star %||% 0.8,
    verbose = 1000L
) {
    validate_s2_fit_objects(dat, priors, spatial)

    if (!is.matrix(lambda_fixed) || !identical(dim(lambda_fixed), dim(dat$y_coarse))) {
        stop("lambda_fixed must have the same dimensions as dat$y_coarse.", call. = FALSE)
    }
    if (any(!is.finite(lambda_fixed)) || any(lambda_fixed <= 0)) {
        stop("lambda_fixed must be positive and finite.", call. = FALSE)
    }

    n_iter <- as.integer(settings$n_iter %||% 20000L)
    n_burnin <- as.integer(settings$n_burnin %||% 10000L)
    n_thin <- as.integer(settings$n_thin %||% 5L)
    adapt_interval <- as.integer(settings$adapt_interval %||% 50L)

    if (n_iter <= 0L || n_burnin < 0L || n_burnin >= n_iter || n_thin <= 0L) {
        stop("Invalid MCMC iteration settings.", call. = FALSE)
    }

    TT_use <- as.integer(dat$TT)
    n1 <- as.integer(dat$n1)
    p <- length(priors$beta_mean)
    n_stored <- as.integer(floor((n_iter - n_burnin) / n_thin))

    if (length(gamma_fixed) == 1L) {
        gamma_fixed <- rep(gamma_fixed, n1)
    }

    samples <- list(
        beta0 = numeric(n_stored),
        beta = matrix(NA_real_, nrow = n_stored, ncol = p),
        phi = matrix(NA_real_, nrow = n_stored, ncol = n1),
        tau_phi = numeric(n_stored),
        r = matrix(NA_real_, nrow = n_stored, ncol = n1),
        gamma = matrix(rep(gamma_fixed, each = n_stored), nrow = n_stored, ncol = n1),
        loglik_nb = numeric(n_stored)
    )

    diagnostics <- list(
        loglik_nb_trace = rep(NA_real_, n_iter),
        beta_n_reject = rep(NA_real_, n_iter),
        phi_accept_trace = rep(FALSE, n_iter),
        r_accept_trace = matrix(FALSE, nrow = n_iter, ncol = n1),
        phi_log_alpha = rep(NA_real_, n_iter),
        r_log_alpha = matrix(NA_real_, nrow = n_iter, ncol = n1),
        lambda_min_trace = rep(min(lambda_fixed), n_iter),
        lambda_max_trace = rep(max(lambda_fixed), n_iter),
        fixed_gamma = gamma_fixed
    )

    state <- initialise_s2_state(dat, settings, priors, spatial, gamma_fixed)
    state$lambda_tilde <- lambda_fixed
    state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                              state$beta0, state$beta, state$phi)

    phi_accept_window <- 0L
    r_accept_window <- rep(0L, n1)
    store_idx <- 0L
    start_time <- proc.time()

    for (iter in seq_len(n_iter)) {
        ## 1. beta update using marginal NB likelihood, conditional on oracle lambda.
        beta_result <- update_beta_s2(
            beta_current = c(state$beta0, state$beta),
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = lambda_fixed,
            phi = state$phi,
            priors = priors,
            ess_base = state$ess_base,
            r = state$r,
            use_preconditioned = settings$use_preconditioned_beta %||% TRUE
        )
        state$beta0 <- beta_result$sample[1]
        state$beta <- beta_result$sample[2:(p + 1L)]
        state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                                  state$beta0, state$beta, state$phi)

        ## 2. phi update using marginal NB likelihood, conditional on oracle lambda.
        phi_result <- update_phi_s2(
            u_current = state$u,
            B = spatial$B_ICAR,
            BHB = spatial$BHB,
            tau_phi = state$tau_phi,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            lambda_tilde = lambda_fixed,
            r = state$r,
            proposal_sd = state$phi_proposal_sd
        )
        state$u <- phi_result$u
        state$phi <- phi_result$phi
        state$xi <- compute_s2_xi(dat$e, dat$x1, dat$x2,
                                  state$beta0, state$beta, state$phi)

        ## 3. tau_phi update.
        state$tau_phi <- update_tau_phi(
            phi = state$phi,
            H = spatial$H,
            n1 = n1,
            priors = priors
        )

        ## 4. r update using marginal NB likelihood, conditional on oracle lambda.
        r_result <- update_r_s2(
            r_current = state$r,
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            lambda_tilde = lambda_fixed,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            priors = priors,
            mh_sd = state$r_proposal_sd,
            method = "marginal_nb",
            return_diag = TRUE
        )
        state$r <- r_result$r

        loglik_nb <- compute_s2_loglik_nb(
            y_coarse = dat$y_coarse,
            e = dat$e,
            x1 = dat$x1,
            x2 = dat$x2,
            beta0 = state$beta0,
            beta = state$beta,
            phi = state$phi,
            r = state$r,
            lambda_tilde = lambda_fixed
        )

        diagnostics$loglik_nb_trace[iter] <- loglik_nb
        diagnostics$beta_n_reject[iter] <- beta_result$n_reject %||% NA_real_
        diagnostics$phi_accept_trace[iter] <- isTRUE(phi_result$accept)
        diagnostics$r_accept_trace[iter, ] <- r_result$accept
        diagnostics$phi_log_alpha[iter] <- phi_result$log_alpha %||% NA_real_
        diagnostics$r_log_alpha[iter, ] <- r_result$diag$log_alpha %||% rep(NA_real_, n1)

        phi_accept_window <- phi_accept_window + as.integer(isTRUE(phi_result$accept))
        r_accept_window <- r_accept_window + as.integer(r_result$accept)

        if (iter <= n_burnin && iter %% adapt_interval == 0L) {
            state$phi_proposal_sd <- adapt_sd_s2(
                current_sd = state$phi_proposal_sd,
                n_accept = phi_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$phi_target_accept %||% 0.25
            )
            state$r_proposal_sd <- adapt_sd_s2(
                current_sd = state$r_proposal_sd,
                n_accept = r_accept_window,
                n_trials = adapt_interval,
                target_rate = settings$r_target_accept %||% 0.30
            )
            phi_accept_window <- 0L
            r_accept_window <- rep(0L, n1)
        }

        if (iter > n_burnin && (iter - n_burnin) %% n_thin == 0L) {
            store_idx <- store_idx + 1L
            samples$beta0[store_idx] <- state$beta0
            samples$beta[store_idx, ] <- state$beta
            samples$phi[store_idx, ] <- state$phi
            samples$tau_phi[store_idx] <- state$tau_phi
            samples$r[store_idx, ] <- state$r
            samples$loglik_nb[store_idx] <- loglik_nb
        }

        if (verbose > 0L && iter %% verbose == 0L) {
            elapsed <- (proc.time() - start_time)[3]
            i0 <- max(1L, iter - 99L)
            phi_rate <- mean(diagnostics$phi_accept_trace[i0:iter])
            r_rate <- mean(diagnostics$r_accept_trace[i0:iter, , drop = FALSE])
            beta_rej <- mean(diagnostics$beta_n_reject[i0:iter], na.rm = TRUE)
            cat(sprintf(
                paste0(
                    "iter %5d/%d [%.0fs] oracle loglik_nb=%.1f beta0=%.3f ",
                    "beta=(%.3f, %.3f) r_mean=%.2f lambda_truth=[%.3f, %.3f] ",
                    "| phi_acc=%.2f r_acc=%.2f beta_rej=%.1f\n"
                ),
                iter, n_iter, elapsed, loglik_nb, state$beta0, state$beta[1],
                state$beta[2], mean(state$r), min(lambda_fixed), max(lambda_fixed),
                phi_rate, r_rate, beta_rej
            ))
        }
    }

    elapsed_total <- (proc.time() - start_time)[3]
    diagnostics$elapsed_sec <- elapsed_total
    diagnostics$phi_accept_rate <- mean(diagnostics$phi_accept_trace)
    diagnostics$r_accept_rate <- colMeans(diagnostics$r_accept_trace)
    diagnostics$beta_mean_n_reject <- mean(diagnostics$beta_n_reject, na.rm = TRUE)
    diagnostics$phi_proposal_sd_final <- state$phi_proposal_sd
    diagnostics$r_proposal_sd_final <- state$r_proposal_sd
    diagnostics$gamma_fixed_mean <- mean(gamma_fixed)

    list(
        samples = samples,
        diagnostics = diagnostics,
        final_state = state,
        settings = settings,
        priors = priors,
        spatial = list(H = spatial$H, B_ICAR = spatial$B_ICAR, BHB = spatial$BHB),
        metadata = list(
            model = "S2 oracle-lambda NB-ICAR diagnostic",
            fixed_lambda = TRUE,
            dynamic_lambda_in_truth = TRUE,
            lambda_fixed_to_truth = TRUE,
            fixed_gamma = TRUE,
            gamma_fixed_value = gamma_fixed,
            updated_blocks = c("beta", "phi", "tau_phi", "r"),
            disabled_blocks = c("lambda", "kappa", "gamma", "delta", "omega"),
            uses_marginal_nb_for_beta_phi_r = TRUE,
            uses_kappa_for_lambda_ffbs = FALSE,
            uses_recentered_identified_parameterization = TRUE
        ),
        n_stored = n_stored,
        lambda_fixed = lambda_fixed
    )
}

summarise_s2_oracle_lambda_fit <- function(fit, dat, scenario_id = NULL, rep_id = NULL) {
    s <- fit$samples

    beta0_q <- as.numeric(stats::quantile(s$beta0, c(0.025, 0.5, 0.975), na.rm = TRUE))
    beta_q <- apply(s$beta, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    phi_mean <- colMeans(s$phi)
    r_mean_by_region <- colMeans(s$r)
    r_q_by_region <- apply(s$r, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    beta0_true <- dat$beta0_star_ident %||% dat$beta0_star %||% NA_real_
    beta_true <- dat$beta_star_ident %||% dat$beta_star %||% rep(NA_real_, 2)
    phi_true <- dat$phi_star_ident %||% dat$phi_star %||% rep(NA_real_, dat$n1)
    r_true <- dat$r_star %||% rep(NA_real_, dat$n1)
    if (length(r_true) == 1L) {
        r_true <- rep(r_true, dat$n1)
    }

    lambda_fixed <- fit$lambda_fixed %||% get_s2_oracle_lambda_truth(dat, lambda_source = "identified")
    phi_rmse <- sqrt(mean((phi_mean - phi_true)^2))
    phi_cor <- suppressWarnings(stats::cor(phi_mean, phi_true))

    data.frame(
        scenario_id = scenario_id %||% dat$scenario_id %||% "S2_ORACLE_LAMBDA",
        rep_id = rep_id %||% dat$rep_id %||% NA_integer_,
        TT = dat$TT,
        n1 = dat$n1,
        mean_count = mean(dat$y_coarse),
        zero_prop = mean(dat$y_coarse == 0),
        dynamic_lambda_in_truth = TRUE,
        oracle_lambda_in_fit = TRUE,
        lambda_sampled_in_fit = FALSE,
        gamma_fixed_in_fit = isTRUE(fit$metadata$fixed_gamma),
        gamma_truth_mean = mean(dat$gamma_star),
        gamma_fixed_mean = mean(fit$metadata$gamma_fixed_value),
        lambda_truth_min = min(lambda_fixed),
        lambda_truth_max = max(lambda_fixed),

        beta0_true = beta0_true,
        beta0_mean = mean(s$beta0),
        beta0_sd = stats::sd(s$beta0),
        beta0_q025 = beta0_q[1],
        beta0_q50 = beta0_q[2],
        beta0_q975 = beta0_q[3],
        beta0_bias = mean(s$beta0) - beta0_true,
        beta0_covered = as.integer(beta0_q[1] <= beta0_true && beta0_true <= beta0_q[3]),

        beta1_true = beta_true[1],
        beta1_mean = mean(s$beta[, 1]),
        beta1_sd = stats::sd(s$beta[, 1]),
        beta1_q025 = beta_q[1, 1],
        beta1_q50 = beta_q[2, 1],
        beta1_q975 = beta_q[3, 1],
        beta1_bias = mean(s$beta[, 1]) - beta_true[1],
        beta1_covered = as.integer(beta_q[1, 1] <= beta_true[1] && beta_true[1] <= beta_q[3, 1]),

        beta2_true = beta_true[2],
        beta2_mean = mean(s$beta[, 2]),
        beta2_sd = stats::sd(s$beta[, 2]),
        beta2_q025 = beta_q[1, 2],
        beta2_q50 = beta_q[2, 2],
        beta2_q975 = beta_q[3, 2],
        beta2_bias = mean(s$beta[, 2]) - beta_true[2],
        beta2_covered = as.integer(beta_q[1, 2] <= beta_true[2] && beta_true[2] <= beta_q[3, 2]),

        phi_rmse = phi_rmse,
        phi_cor = phi_cor,
        phi_mean_bias = mean(phi_mean - phi_true),

        r_true_mean = mean(r_true),
        r_mean = mean(r_mean_by_region),
        r_bias = mean(r_mean_by_region) - mean(r_true),
        r_q025_mean = mean(r_q_by_region[1, ]),
        r_q50_mean = mean(r_q_by_region[2, ]),
        r_q975_mean = mean(r_q_by_region[3, ]),

        tau_phi_mean = mean(s$tau_phi),
        loglik_nb_mean = mean(s$loglik_nb),
        phi_accept_rate = fit$diagnostics$phi_accept_rate,
        r_accept_rate_mean = mean(fit$diagnostics$r_accept_rate),
        beta_mean_n_reject = fit$diagnostics$beta_mean_n_reject,
        elapsed_sec = fit$diagnostics$elapsed_sec,
        stringsAsFactors = FALSE
    )
}

fit_s2_oracle_lambda_one_rep <- function(
    rep_id,
    scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
    data_dir = "data_revised",
    output_dir = "output_s2_oracle_lambda_diagnostic",
    root = ".",
    settings_override = list(),
    priors = MCMC_PRIORS,
    spatial = build_s2_spatial(),
    lambda_source = c("identified", "raw"),
    gamma_fixed = NULL,
    verbose = 1000L,
    save_result = TRUE,
    return_result = TRUE
) {
    lambda_source <- match.arg(lambda_source)
    rr <- sprintf("%02d", as.integer(rep_id))
    dat_file <- file.path(root, data_dir, scenario_id, paste0("data_rep", rr, ".rds"))
    .require_file_s2_oracle(dat_file)
    dat <- readRDS(dat_file)
    validate_s2_data(dat)

    settings <- MCMC_SETTINGS
    if (length(settings_override) > 0L) {
        for (nm in names(settings_override)) {
            settings[[nm]] <- settings_override[[nm]]
        }
    }

    lambda_fixed <- get_s2_oracle_lambda_truth(dat, lambda_source = lambda_source)
    gamma_use <- gamma_fixed %||% dat$gamma_star

    cat(sprintf("=== Scenario 2 oracle-lambda diagnostic: rep %s ===\n", rr))
    cat(sprintf("Data file : %s\n", dat_file))
    cat(sprintf("Fixed     : lambda = %s truth, gamma = %.3f on average\n", lambda_source, mean(gamma_use)))
    cat("Estimated : beta0, beta1, beta2, phi, tau_phi, r\n")
    cat("Disabled  : lambda, kappa, gamma, delta, omega updates\n\n")

    fit <- run_s2_oracle_lambda_mcmc(
        dat = dat,
        settings = settings,
        priors = priors,
        spatial = spatial,
        lambda_fixed = lambda_fixed,
        gamma_fixed = gamma_use,
        verbose = verbose
    )

    summary <- summarise_s2_oracle_lambda_fit(
        fit = fit,
        dat = dat,
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id)
    )

    fit$summary <- summary
    fit$metadata <- c(
        fit$metadata,
        list(
            scenario_id = scenario_id,
            rep_id = as.integer(rep_id),
            data_file = dat_file,
            output_dir = output_dir,
            lambda_source = lambda_source,
            run_time = Sys.time()
        )
    )

    if (isTRUE(save_result)) {
        out_dir <- file.path(root, output_dir, scenario_id)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        fit_file <- file.path(out_dir, paste0("fit_S2_oracle_lambda_rep", rr, ".rds"))
        csv_file <- file.path(out_dir, paste0("summary_S2_oracle_lambda_rep", rr, ".csv"))
        saveRDS(fit, fit_file)
        write.csv(summary, csv_file, row.names = FALSE)
        cat(sprintf("Saved fit    : %s\n", fit_file))
        cat(sprintf("Saved summary: %s\n\n", csv_file))
    }

    if (isTRUE(return_result)) {
        return(fit)
    }
    invisible(NULL)
}

sanity_check_s2_oracle_lambda <- function(
    root = ".",
    rep_id = 1L,
    scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
    data_dir = "data_revised",
    output_dir = "output_s2_oracle_lambda_diagnostic",
    lambda_source = c("identified", "raw"),
    n_iter = 10000L,
    n_burnin = 5000L,
    n_thin = 5L,
    verbose = 1000L
) {
    lambda_source <- match.arg(lambda_source)
    fit_s2_oracle_lambda_one_rep(
        rep_id = rep_id,
        scenario_id = scenario_id,
        data_dir = data_dir,
        output_dir = output_dir,
        root = root,
        settings_override = list(
            n_iter = as.integer(n_iter),
            n_burnin = as.integer(n_burnin),
            n_thin = as.integer(n_thin)
        ),
        lambda_source = lambda_source,
        verbose = verbose,
        save_result = TRUE,
        return_result = TRUE
    )
}

compare_s2_sampled_vs_oracle_beta <- function(
    sampled_summary_file,
    oracle_summary_file
) {
    .require_file_s2_oracle(sampled_summary_file)
    .require_file_s2_oracle(oracle_summary_file)
    sampled <- read.csv(sampled_summary_file)
    oracle <- read.csv(oracle_summary_file)

    data.frame(
        fit = c("sampled_lambda", "oracle_lambda"),
        beta0_mean = c(sampled$beta0_mean[1], oracle$beta0_mean[1]),
        beta0_covered = c(sampled$beta0_covered[1], oracle$beta0_covered[1]),
        beta1_mean = c(sampled$beta1_mean[1], oracle$beta1_mean[1]),
        beta1_covered = c(sampled$beta1_covered[1], oracle$beta1_covered[1]),
        beta2_mean = c(sampled$beta2_mean[1], oracle$beta2_mean[1]),
        beta2_covered = c(sampled$beta2_covered[1], oracle$beta2_covered[1]),
        r_mean = c(sampled$r_mean[1], oracle$r_mean[1]),
        loglik_nb_mean = c(sampled$loglik_nb_mean[1], oracle$loglik_nb_mean[1]),
        stringsAsFactors = FALSE
    )
}
