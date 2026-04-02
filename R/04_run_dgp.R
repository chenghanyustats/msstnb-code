## ==========================================================================
## 04_run_dgp.R
## Generate all simulated datasets: 5 scenarios × 50 replications
##
## STEP 4 REVISION:
## Within each replication, S1--S4 share the same baseline phi_star.
## S5 uses that same phi_star and adds a separate phi_after after the
## changepoint.
## ==========================================================================

#' Generate and save all simulated datasets
#'
#' @param reps       Integer vector of replication indices (default 1:R_TOTAL)
#' @param scenarios  Character vector of scenario IDs (default "S1"–"S5")
#' @param data_dir   Base directory for saving (default "data")
#' @param verbose    Print progress?
run_all_dgp <- function(reps      = seq_len(R_TOTAL),
                        scenarios = paste0("S", 1:5),
                        data_dir  = "data",
                        verbose   = TRUE) {

    ## ---- structural constants (from 00_setup.R) ----
    constants <- list(
        TT         = TT,
        N1         = N1,
        N_CHILDREN = N_CHILDREN,
        A0         = A0,
        B0         = B0,
        C0         = C0
    )

    ## ---- ensure output directories exist ----
    for (sid in scenarios) {
        dir.create(file.path(data_dir, sid), recursive = TRUE, showWarnings = FALSE)
    }
    dir.create(file.path(data_dir, "shared"), recursive = TRUE, showWarnings = FALSE)

    ## ---- loop over replications ----
    n_total <- length(reps) * length(scenarios)
    counter <- 0L

    for (r in reps) {
        rr <- sprintf("%02d", r)

        ## Set replication seed (shared across scenarios)
        set.seed(REP_SEEDS[r])

        ## Generate shared inputs (exposures, covariates) — same for all scenarios
        inputs <- generate_inputs(TT = TT, n1 = N1)
        shared_file <- file.path(data_dir, "shared",
                                 paste0("shared_inputs_rep", rr, ".rds"))
        saveRDS(inputs, shared_file)

        ## Generate scenario parameters. Under the revised scenario builder,
        ## S1--S4 share the same baseline phi_star within this replication.
        ## S5 reuses that same baseline phi_star and adds a separate phi_after.
        all_params <- all_scenarios(H = H, B = B_ICAR, n1 = N1)

        if (verbose) {
            phi_ref <- all_params[["S1"]]$phi
            same_s1_s4 <- all(vapply(c("S2", "S3", "S4"), function(sid) {
                isTRUE(all.equal(all_params[[sid]]$phi, phi_ref, tolerance = 0))
            }, logical(1)))
            cat(sprintf("Replication %s: shared phi across S1--S4 = %s\n",
                        rr, ifelse(same_s1_s4, "TRUE", "FALSE")))
        }

        for (sid in scenarios) {
            counter <- counter + 1L
            params  <- all_params[[sid]]

            ## Simulate dataset
            dat <- simulate_one_dataset(inputs    = inputs,
                                        params    = params,
                                        constants = constants)

            ## Save
            out_file <- file.path(data_dir, sid,
                                  paste0("data_rep", rr, ".rds"))
            saveRDS(dat, out_file)

            if (verbose) {
                cat(sprintf("[%3d/%d]  %s  rep %s  → %s  (y range: %d–%d)\n",
                            counter, n_total, sid, rr, out_file,
                            min(dat$y_coarse), max(dat$y_coarse)))
            }
        }
    }

    if (verbose) cat("\nDone. Total files:", counter, "\n")
    invisible(counter)
}


#' Quick sanity check: generate 1 replication of S1, print summary
#'
#' @return The simulated dataset (invisibly)
sanity_check <- function() {
    cat("Running sanity check (S1, rep 01)...\n\n")

    set.seed(REP_SEEDS[1])
    inputs <- generate_inputs(TT = TT, n1 = N1)

    all_params <- all_scenarios(H = H, B = B_ICAR, n1 = N1)
    params_s1 <- all_params[["S1"]]
    params_s4 <- all_params[["S4"]]

    print_scenario(params_s1)
    cat(sprintf("Shared phi between S1 and S4: %s\n\n",
                ifelse(isTRUE(all.equal(params_s1$phi, params_s4$phi, tolerance = 0)),
                       "TRUE", "FALSE")))

    constants <- list(TT = TT, N1 = N1, N_CHILDREN = N_CHILDREN,
                      A0 = A0, B0 = B0, C0 = C0)
    dat <- simulate_one_dataset(inputs = inputs, params = params_s1,
                                constants = constants)

    cat("--- Data summary ---\n")
    cat(sprintf("  y_coarse: %d × %d, range [%d, %d], mean %.1f\n",
                nrow(dat$y_coarse), ncol(dat$y_coarse),
                min(dat$y_coarse), max(dat$y_coarse),
                mean(dat$y_coarse)))
    cat(sprintf("  y_fine:   %d × %d × %d\n",
                dim(dat$y_fine)[1], dim(dat$y_fine)[2], dim(dat$y_fine)[3]))
    cat(sprintf("  λ̃ range:  [%.3f, %.3f], mean %.3f\n",
                min(dat$lambda_tilde), max(dat$lambda_tilde),
                mean(dat$lambda_tilde)))
    cat(sprintf("  κ range:  [%.3f, %.3f], mean %.3f\n",
                min(dat$kappa), max(dat$kappa), mean(dat$kappa)))
    cat(sprintf("  Tree coherent: %s\n", dat$coherent))
    cat("\n")

    invisible(dat)
}
