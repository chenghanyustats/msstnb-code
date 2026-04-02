# ============================================================
# Run frozen-a gamma diagnostic on an existing ordinary fit
#
# Purpose:
#   Connect the current ordinary fit to frozen_a_gamma_diagnostic.R
#   with minimal manual editing.
#
# What this script does:
#   1. Reads an existing ordinary fit (.rds)
#   2. Reads the matching dataset (.rds)
#   3. Sources frozen_a_gamma_diagnostic.R
#   4. Runs the diagnostic
#   5. Saves all outputs using a user-specified prefix
#
# Usage example:
#
# source("run_frozen_a_on_existing_fit.R")
#
# out <- run_frozen_a_on_existing_fit(
#     data_file   = "data_revised/S1/data_rep01.rds",
#     fit_file    = "fit_revised_S1_rep01_short.rds",
#     script_file = "frozen_a_gamma_diagnostic.R",
#     out_prefix  = "diag_frozen_a_S1_rep01",
#     gamma_ref_mode = "truth",
#     lambda_a0   = 10
# )
#
# ============================================================

run_frozen_a_on_existing_fit <- function(
    data_file,
    fit_file,
    script_file = "frozen_a_gamma_diagnostic.R",
    out_prefix = "diag_frozen_a",
    gamma_ref_mode = "truth",
    gamma_ref_value = NULL,
    lambda_a0 = 10,
    verbose = TRUE
) {
    # --------------------------------------------------------
    # Basic checks
    # --------------------------------------------------------
    if (!file.exists(data_file)) {
        stop("Data file not found: ", data_file)
    }
    if (!file.exists(fit_file)) {
        stop("Fit file not found: ", fit_file)
    }
    if (!file.exists(script_file)) {
        stop("Diagnostic script not found: ", script_file,
             "\nProvide the full path to frozen_a_gamma_diagnostic.R.")
    }

    # --------------------------------------------------------
    # Read data / fit for sanity checks before sourcing
    # --------------------------------------------------------
    dat <- readRDS(data_file)
    fit <- readRDS(fit_file)

    # Check that this looks like the ordinary fit / data pair
    if (is.null(dat$gamma_star)) {
        stop("The dataset does not contain gamma_star. ",
             "This diagnostic expects a simulated dataset with truth gamma.")
    }

    if (is.null(fit$samples) && is.null(fit$draws)) {
        stop("The fit object does not appear to contain posterior samples/draws.")
    }

    # --------------------------------------------------------
    # Source the diagnostic script
    # --------------------------------------------------------
    source(script_file, local = FALSE)

    if (!exists("run_frozen_a_gamma_diagnostic", mode = "function")) {
        stop("After sourcing ", script_file,
             ", function run_frozen_a_gamma_diagnostic() was not found.")
    }

    # --------------------------------------------------------
    # Run the diagnostic
    # --------------------------------------------------------
    if (verbose) {
        cat("============================================================\n")
        cat("Running frozen-a gamma diagnostic on existing ordinary fit\n")
        cat("============================================================\n")
        cat("Data file   :", data_file, "\n")
        cat("Fit file    :", fit_file, "\n")
        cat("Script file :", script_file, "\n")
        cat("Out prefix  :", out_prefix, "\n")
        cat("gamma_ref_mode :", gamma_ref_mode, "\n")
        if (!is.null(gamma_ref_value)) {
            cat("gamma_ref_value:", gamma_ref_value, "\n")
        }
        cat("lambda_a0      :", lambda_a0, "\n\n")
    }

    args <- list(
        data_file = data_file,
        fit_file = fit_file,
        out_prefix = out_prefix,
        gamma_ref_mode = gamma_ref_mode,
        lambda_a0 = lambda_a0
    )

    if (!is.null(gamma_ref_value)) {
        args$gamma_ref_value <- gamma_ref_value
    }

    out <- do.call(run_frozen_a_gamma_diagnostic, args)

    if (verbose) {
        cat("\n============================================================\n")
        cat("Frozen-a diagnostic finished\n")
        cat("============================================================\n")
        cat("Expected output prefix:", out_prefix, "\n")
        cat("Typical files:\n")
        cat("  ", out_prefix, ".png\n", sep = "")
        cat("  ", out_prefix, ".pdf\n", sep = "")
        cat("  ", out_prefix, "_summary.csv\n", sep = "")
        cat("  ", out_prefix, "_diagnostic.rds\n", sep = "")
    }

    invisible(out)
}


# ------------------------------------------------------------
# Optional convenience wrapper for the current ordinary fit
# ------------------------------------------------------------
run_frozen_a_on_current_ordinary_fit <- function(
    data_file   = "data_revised/S1/data_rep01.rds",
    fit_file    = "fit_revised_S1_rep01_short.rds",
    script_file = "frozen_a_gamma_diagnostic.R",
    out_prefix  = "diag_frozen_a_S1_rep01",
    gamma_ref_mode = "truth",
    gamma_ref_value = NULL,
    lambda_a0   = 10,
    verbose = TRUE
) {
    run_frozen_a_on_existing_fit(
        data_file = data_file,
        fit_file = fit_file,
        script_file = script_file,
        out_prefix = out_prefix,
        gamma_ref_mode = gamma_ref_mode,
        gamma_ref_value = gamma_ref_value,
        lambda_a0 = lambda_a0,
        verbose = verbose
    )
}
