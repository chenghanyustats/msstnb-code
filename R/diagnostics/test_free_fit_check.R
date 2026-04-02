## ==========================================================================
## test_mcmc_step6_pathfixed.R
## Complete testing script with identified-space posterior checks
## Path-fixed version for revised DGP output
## ==========================================================================

cat("\n============================================================\n")
cat("  MSSTNB MCMC Sampler — Test Suite (identified truth, path-fixed)\n")
cat("============================================================\n\n")

## User-configurable paths -------------------------------------------------
DAT_FILE_CANDIDATES <- c(
    "data_revised/S1/data_rep01.rds",
    "./data_revised/S1/data_rep01.rds",
    "data/S1/data_rep01.rds",
    "./data/S1/data_rep01.rds"
)
FIT_FILE_CANDIDATES <- c(
    "fit_revised_S1_rep01_short.rds",
    "./fit_revised_S1_rep01_short.rds",
    "fit.Rds",
    "./fit.Rds"
)

find_first_existing <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit) == 0L) return(NA_character_)
    hit[1]
}

DAT_FILE <- find_first_existing(DAT_FILE_CANDIDATES)
FIT_FILE <- find_first_existing(FIT_FILE_CANDIDATES)

## Helpers -----------------------------------------------------------------
require_identified_truth <- function(dat) {
    needed <- c("beta0_star_ident", "phi_star_ident", "lambda_tilde_ident")
    missing <- needed[!needed %in% names(dat)]
    if (length(missing) > 0L) {
        stop(
            "Dataset is missing identified-truth fields: ",
            paste(missing, collapse = ", "),
            ". Re-generate data with the revised DGP before running this test."
        )
    }

    list(
        beta0        = dat$beta0_star_ident,
        beta         = dat$beta_star,
        phi          = dat$phi_star_ident,
        lambda_tilde = dat$lambda_tilde_ident,
        gamma        = dat$gamma_star,
        r            = dat$r_star,
        delta        = dat$delta_star
    )
}

check_param <- function(name, true_val, samples) {
    pm  <- mean(samples)
    psd <- stats::sd(samples)
    ci  <- stats::quantile(samples, c(0.025, 0.975))
    covered <- (true_val >= ci[1]) && (true_val <= ci[2])
    status  <- ifelse(covered, "OK", "MISS")
    cat(sprintf(
        "  %-16s %10.3f %10.3f %10.3f [%6.2f,%6.2f]  %s\n",
        name, true_val, pm, psd, ci[1], ci[2], status
    ))
    covered
}

summarise_latent_recovery <- function(samples_arr, truth_mat, label = "lambda_tilde") {
    eps <- .Machine$double.xmin

    post_mean <- apply(samples_arr, c(2, 3), mean)
    ci_lo     <- apply(samples_arr, c(2, 3), stats::quantile, 0.025)
    ci_hi     <- apply(samples_arr, c(2, 3), stats::quantile, 0.975)

    truth_vec <- as.vector(truth_mat)
    mean_vec  <- as.vector(post_mean)
    lo_vec    <- as.vector(ci_lo)
    hi_vec    <- as.vector(ci_hi)

    log_cor   <- stats::cor(log(pmax(mean_vec, eps)), log(pmax(truth_vec, eps)))
    log_rmse  <- sqrt(mean((log(pmax(mean_vec, eps)) - log(pmax(truth_vec, eps)))^2))
    coverage  <- mean(truth_vec >= lo_vec & truth_vec <= hi_vec)

    cat(sprintf("  %-16s %10s %10.3f\n", paste0("cor(log ", label, ")"), "", log_cor))
    cat(sprintf("  %-16s %10s %10.3f\n", paste0("RMSE(log ", label, ")"), "", log_rmse))
    cat(sprintf("  %-16s %10s %10.3f\n", paste0(label, " cover"), "", coverage))

    list(
        post_mean = post_mean,
        ci_lo     = ci_lo,
        ci_hi     = ci_hi,
        log_cor   = log_cor,
        log_rmse  = log_rmse,
        coverage  = coverage
    )
}

## STAGE 1 -----------------------------------------------------------------
cat("--- Stage 1: Source all files ---\n")
tryCatch({
    source("R/00_setup.R")
    source("R/01_helpers.R")
    source("R/mcmc/mcmc_config.R")
    source("R/mcmc/mcmc_utils.R")
    source("R/mcmc/update_kappa.R")
    source("R/mcmc/update_regression.R")
    source("R/mcmc/update_icar.R")
    source("R/mcmc/update_dispersion.R")
    source("R/mcmc/update_gamma.R")
    source("R/mcmc/ffbs_lambda.R")
    source("R/mcmc/update_delta.R")
    source("R/mcmc/smooth_omega.R")
    source("R/mcmc/sampler.R")
    source("R/mcmc/run_single_fit.R")
    cat("  All files sourced.\n  STAGE 1: PASS\n\n")
}, error = function(e) {
    cat("  STAGE 1: FAIL —", conditionMessage(e), "\n")
    stop("Fix source errors before proceeding.")
})

## STAGE 2 -----------------------------------------------------------------
cat("--- Stage 2: Single iteration ---\n")
if (is.na(DAT_FILE)) stop("No data file found. Expected one of: ", paste(DAT_FILE_CANDIDATES, collapse = ", "))
dat <- readRDS(DAT_FILE)
truth <- require_identified_truth(dat)

cat("  Loaded:", DAT_FILE, "\n")
cat(
    "  Identified truth: beta0=", round(truth$beta0, 4),
    " beta=", paste(round(truth$beta, 4), collapse = ","),
    " r[1]=", round(truth$r[1], 4),
    " gamma[1]=", round(truth$gamma[1], 4),
    "\n"
)

settings  <- method_settings("M1")
priors    <- MCMC_PRIORS
constants <- list(A0 = A0, B0 = B0, C0 = C0)
BHB       <- crossprod(B_ICAR, H %*% B_ICAR)
spatial   <- list(H = H, B_ICAR = B_ICAR, BHB = BHB)

state <- initialise_state(dat, settings, priors, spatial, constants)

tryCatch({
    result <- run_one_iteration(state, dat, settings, priors, spatial, constants)
    s2 <- result$state
    d2 <- result$diag

    cat("  One iteration completed.\n")
    cat("  loglik:", round(d2$loglik, 1), "\n")
    cat("  beta0:", round(state$beta0, 4), "->", round(s2$beta0, 4), "\n")
    cat("  gamma[1]:", round(state$gamma[1], 4), "->", round(s2$gamma[1], 4), "\n")

    stopifnot(is.finite(d2$loglik))
    stopifnot(all(s2$lambda_tilde > 0))
    stopifnot(all(s2$kappa > 0))
    stopifnot(all(s2$gamma > 0 & s2$gamma < 1))
    stopifnot(abs(sum(s2$phi)) < 1e-8)

    cat("  All constraints OK.\n  STAGE 2: PASS\n\n")
}, error = function(e) {
    cat("  STAGE 2: FAIL —", conditionMessage(e), "\n")
    stop("Fix single-iteration errors before proceeding.")
})

## STAGE 3 -----------------------------------------------------------------
cat("--- Stage 3: Short MCMC run (2000 iter, M1, S1 rep01) ---\n")
cat("  This may take 2-5 minutes...\n\n")

settings_short          <- method_settings("M1")
settings_short$n_iter   <- 2000L
settings_short$n_burnin <- 1000L
settings_short$n_thin   <- 2L

if (!is.na(FIT_FILE)) {
    cat("  Using existing fit:", FIT_FILE, "\n")
    fit <- readRDS(FIT_FILE)
} else {
    fit <- run_mcmc(dat, settings_short, priors, spatial, constants, verbose = 500L)
}

s   <- fit$samples

cat("\n--- Posterior summary vs identified true values ---\n")
cat(sprintf(
    "  %-16s %10s %10s %10s %15s  %s\n",
    "Parameter", "True", "PostMean", "PostSD", "95% CI", "Status"
))

results <- c()
results["beta0"]  <- check_param("beta0",    truth$beta0,    s$beta0)
results["beta1"]  <- check_param("beta1",    truth$beta[1],  s$beta[, 1])
results["beta2"]  <- check_param("beta2",    truth$beta[2],  s$beta[, 2])
results["gamma1"] <- check_param("gamma[1]", truth$gamma[1], s$gamma[, 1])
results["r1"]     <- check_param("r[1]",     truth$r[1],     s$r[, 1])
results["delta"]  <- check_param("delta",    truth$delta,    s$delta)

phi_pm  <- colMeans(s$phi)
phi_cor <- stats::cor(phi_pm, truth$phi)
cat(sprintf("  %-16s %10s %10.3f\n", "cor(phi_ident)", "", phi_cor))
results["phi_cor"] <- (phi_cor > 0.5)

cat("\n--- Latent-path recovery vs identified truth ---\n")
lt_diag <- summarise_latent_recovery(s$lambda_tilde, truth$lambda_tilde, label = "lambda_tilde")
results["lambda_cor"]   <- (lt_diag$log_cor > 0.5)
results["lambda_cover"] <- (lt_diag$coverage > 0.70)

cat("\n--- Acceptance / tuning ---\n")
cat(sprintf("  phi MH:   %.3f  (target 0.25)\n", fit$diagnostics$phi_accept_rate))
cat(sprintf("  gamma MH: %.3f  (target 0.30)\n", mean(fit$diagnostics$gamma_accept_rate)))
cat(sprintf("  delta MH: %.3f  (target 0.30)\n", fit$diagnostics$delta_accept_rate))
cat(sprintf("  r MH:     %.3f  (target 0.30)\n", mean(fit$diagnostics$r_accept_rate)))
cat(sprintf("  beta ESS: %.1f mean rejections/step\n", fit$diagnostics$beta_mean_n_reject))

n_pass <- sum(results)
n_total <- length(results)
cat(sprintf("\n  Posterior checks: %d/%d passed\n", n_pass, n_total))

if (n_pass >= n_total - 1L) {
    cat("  STAGE 3: PASS")
    if (n_pass < n_total) cat(" (1 miss OK with 2000 iter)")
    cat("\n\n")
} else {
    cat("  STAGE 3: FAIL\n\n")
}
