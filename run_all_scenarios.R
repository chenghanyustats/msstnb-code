# ============================================================
# run_all_scenarios.R
# Master runner for MSSTNB simulation scenarios
#
# Purpose:
#   1. Run selected scenarios from one entry point.
#   2. Keep Scenario 1 and Scenario 2 output structure consistent.
#   3. Optionally rerun posterior performance summaries, MCMC figures,
#      and final validation checks without manually sourcing many scripts.
#
# Recommended usage:
#
#   source("run_all_scenarios.R")
#
#   # Run Scenario 1 and Scenario 2 from scratch
#   run_all_scenarios(
#       scenarios = c("S1", "S2"),
#       clean = TRUE,
#       run_mcmc = TRUE,
#       run_postprocess = TRUE,
#       run_validation = TRUE
#   )
#
#   # Only regenerate figures/tables and validation, without rerunning MCMC
#   run_all_scenarios(
#       scenarios = c("S1", "S2"),
#       clean = FALSE,
#       run_mcmc = FALSE,
#       run_postprocess = TRUE,
#       run_validation = TRUE
#   )
# ============================================================


# ------------------------------------------------------------
# Utility helpers
# ------------------------------------------------------------

.source_required <- function(file) {
    if (!file.exists(file)) {
        stop("Required file not found: ", file, call. = FALSE)
    }
    source(file)
    invisible(TRUE)
}

.source_optional <- function(file) {
    if (file.exists(file)) {
        source(file)
        return(TRUE)
    }
    FALSE
}

.clean_path <- function(path, clean = FALSE) {
    if (isTRUE(clean) && dir.exists(path)) {
        message("Removing existing directory: ", path)
        unlink(path, recursive = TRUE, force = TRUE)
    }
    invisible(TRUE)
}

.ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(TRUE)
}


# ------------------------------------------------------------
# Scenario 1
# ------------------------------------------------------------

run_scenario1_all <- function(
    root = ".",
    clean = FALSE,
    run_mcmc = TRUE,
    run_postprocess = TRUE,
    run_extra_mcmc = TRUE
) {
    scenario_id <- "S1_BASELINE_VALIDATION_T100"
    data_dir <- "data_revised"
    output_dir <- "output_s1_baseline_only"
    analysis_dir <- "analysis_s1_baseline_only"

    message("\n============================================================")
    message("Running Scenario 1: baseline validation T100")
    message("============================================================")

    if (isTRUE(clean)) {
        .clean_path(file.path(data_dir, scenario_id), clean = TRUE)
        .clean_path(file.path(output_dir, scenario_id), clean = TRUE)
        .clean_path(file.path(analysis_dir, scenario_id), clean = TRUE)
    }

    # Source the Scenario 1 core if available.
    # Some project versions use s1_baseline.R; some use s1_baseline_only.R.
    if (file.exists("s1_baseline_only.R")) {
        .source_required("s1_baseline_only.R")
        if (exists("source_s1_baseline_only")) {
            source_s1_baseline_only(root = root)
        }
    } else if (file.exists("s1_baseline.R")) {
        .source_required("s1_baseline.R")
        if (exists("source_s1_baseline")) {
            source_s1_baseline(root = root)
        }
    } else {
        stop("Cannot find Scenario 1 core script: s1_baseline_only.R or s1_baseline.R", call. = FALSE)
    }

    if (isTRUE(run_mcmc)) {
        # This runner should generate data, fit all reps, and usually run
        # posterior performance and extra diagnostics depending on its own contents.
        if (file.exists("run_s1_baseline_only_T100.R")) {
            .source_required("run_s1_baseline_only_T100.R")
        } else {
            .source_required("run_s1_baseline_validation_T100_v2.R")
        }
    }

    if (isTRUE(run_postprocess)) {
        .source_required("run_scenario1_posterior_performance.R")

        out_s1 <- run_scenario1_posterior_performance(
            root = root,
            reps = 1:20,
            scenario_id = scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            analysis_dir = analysis_dir,
            plot_format = "png",
            reps_to_show = c(1L, 5L, 10L, 20L)
        )

        s1_rds_file <- file.path(
            analysis_dir,
            scenario_id,
            "scenario1_posterior_performance_results.rds"
        )

        saveRDS(out_s1, s1_rds_file)
        message("Saved Scenario 1 posterior performance RDS: ", s1_rds_file)
    }

    if (isTRUE(run_extra_mcmc)) {
        if (file.exists("run_scenario1_extra_mcmc_diagnostics.R")) {
            .source_required("run_scenario1_extra_mcmc_diagnostics.R")

            run_scenario1_extra_mcmc_diagnostics(
                root = root,
                reps = 1:20,
                scenario_id = scenario_id,
                output_dir = output_dir,
                analysis_dir = analysis_dir,
                plot_format = "png",
                reps_to_show = c(1L, 5L, 10L, 20L),
                selected_phi_regions = c(1L, 5L, 9L),
                selected_r_regions = c(1L, 5L, 9L)
            )
        } else {
            warning("run_scenario1_extra_mcmc_diagnostics.R not found. Skipping Scenario 1 MCMC figures.")
        }
    }

    invisible(TRUE)
}


# ------------------------------------------------------------
# Scenario 2
# ------------------------------------------------------------

run_scenario2_all <- function(
    root = ".",
    clean = FALSE,
    run_mcmc = TRUE,
    run_postprocess = TRUE,
    run_extra_mcmc = TRUE
) {
    scenario_id <- "S2_DYNAMIC_FIXED_GAMMA_T100"
    data_dir <- "data_revised"
    output_dir <- "output_s2_dynamic_fixed_gamma"
    analysis_dir <- "analysis_s2_dynamic_fixed_gamma"

    message("\n============================================================")
    message("Running Scenario 2: dynamic residual risk with fixed gamma T100")
    message("============================================================")

    if (isTRUE(clean)) {
        .clean_path(file.path(data_dir, scenario_id), clean = TRUE)
        .clean_path(file.path(output_dir, scenario_id), clean = TRUE)
        .clean_path(file.path(analysis_dir, scenario_id), clean = TRUE)
    }

    .source_required("s2_dynamic_fixed_gamma.R")
    if (exists("source_s2_dynamic_fixed_gamma")) {
        source_s2_dynamic_fixed_gamma(root = root)
    }

    if (isTRUE(run_mcmc)) {
        .source_required("run_s2_dynamic_fixed_gamma_T100.R")
    }

    if (isTRUE(run_postprocess)) {
        .source_required("run_scenario2_posterior_performance.R")

        out_s2 <- run_scenario2_posterior_performance(
            root = root,
            reps = 1:20,
            scenario_id = scenario_id,
            data_dir = data_dir,
            output_dir = output_dir,
            analysis_dir = analysis_dir,
            plot_format = "png",
            reps_to_show = c(1L, 5L, 10L, 20L),
            selected_regions = c(1L, 5L, 9L)
        )

        s2_rds_file <- file.path(
            analysis_dir,
            scenario_id,
            "scenario2_posterior_performance_results.rds"
        )

        saveRDS(out_s2, s2_rds_file)
        message("Saved Scenario 2 posterior performance RDS: ", s2_rds_file)
    }

    if (isTRUE(run_extra_mcmc)) {
        .source_required("run_scenario2_extra_mcmc_diagnostics.R")

        run_scenario2_extra_mcmc_diagnostics(
            root = root,
            reps = 1:20,
            scenario_id = scenario_id,
            output_dir = output_dir,
            analysis_dir = analysis_dir,
            plot_format = "png",
            reps_to_show = c(1L, 5L, 10L, 20L),
            selected_phi_regions = c(1L, 5L, 9L),
            selected_r_regions = c(1L, 5L, 9L),
            selected_lambda_regions = c(1L, 5L, 9L),
            selected_lambda_times = c(25L, 50L, 75L)
        )
    }

    invisible(TRUE)
}


# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------

run_all_scenarios <- function(
    scenarios = c("S1", "S2"),
    root = ".",
    clean = FALSE,
    run_mcmc = TRUE,
    run_postprocess = TRUE,
    run_extra_mcmc = TRUE,
    run_validation = TRUE,
    stop_on_validation_fail = FALSE
) {
    scenarios <- toupper(scenarios)

    valid_scenarios <- c("S1", "S2", "S3", "S4")
    unknown <- setdiff(scenarios, valid_scenarios)

    if (length(unknown) > 0L) {
        stop(
            "Unknown scenario(s): ",
            paste(unknown, collapse = ", "),
            call. = FALSE
        )
    }

    started_at <- Sys.time()
    message("Started run_all_scenarios at: ", started_at)

    if ("S1" %in% scenarios) {
        run_scenario1_all(
            root = root,
            clean = clean,
            run_mcmc = run_mcmc,
            run_postprocess = run_postprocess,
            run_extra_mcmc = run_extra_mcmc
        )
    }

    if ("S2" %in% scenarios) {
        run_scenario2_all(
            root = root,
            clean = clean,
            run_mcmc = run_mcmc,
            run_postprocess = run_postprocess,
            run_extra_mcmc = run_extra_mcmc
        )
    }

    if ("S3" %in% scenarios) {
        stop("Scenario 3 runner has not been implemented yet.", call. = FALSE)
    }

    if ("S4" %in% scenarios) {
        stop("Scenario 4 runner has not been implemented yet.", call. = FALSE)
    }

    validation <- NULL

    if (isTRUE(run_validation)) {
        if (!file.exists("check_results.R")) {
            warning("check_results.R not found. Skipping validation.")
        } else {
            .source_required("check_results.R")
            validation <- check_results(
                root = root,
                stop_on_fail = stop_on_validation_fail
            )
        }
    }

    finished_at <- Sys.time()
    message("Finished run_all_scenarios at: ", finished_at)
    message("Elapsed time: ", round(difftime(finished_at, started_at, units = "mins"), 2), " minutes")

    invisible(list(
        scenarios = scenarios,
        started_at = started_at,
        finished_at = finished_at,
        validation = validation
    ))
}
