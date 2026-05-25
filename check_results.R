# check_results.R
# Unified validation script for MSSTNB simulation scenarios.
# Current support: Scenario 1 baseline validation and Scenario 2 dynamic fixed-gamma.
# Version: 2026-05-25-clean-beta-counts

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

.add_check <- function(checks, scenario, category, check, status, detail = "") {
    row <- data.frame(
        scenario = scenario,
        category = category,
        check = check,
        status = status,
        detail = as.character(detail),
        stringsAsFactors = FALSE
    )
    rbind(checks, row)
}

.status_from_logical <- function(x, warn_if_false = FALSE) {
    if (isTRUE(x)) {
        "PASS"
    } else if (isTRUE(warn_if_false)) {
        "WARN"
    } else {
        "FAIL"
    }
}

.file_nonempty <- function(path) {
    file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

.dir_nonempty <- function(path) {
    dir.exists(path) && length(list.files(path, all.files = FALSE, recursive = FALSE)) > 0L
}

.read_csv_safe <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    read.csv(path, stringsAsFactors = FALSE)
}

.standardize_beta_parameter <- function(x) {
    x <- as.character(x)
    # Handles names such as beta1.beta1 and beta2.beta2 produced by named vectors.
    ifelse(grepl("beta0", x), "beta0",
        ifelse(grepl("beta1", x), "beta1",
            ifelse(grepl("beta2", x), "beta2", x)
        )
    )
}

.get_beta_mean_col <- function(beta_df) {
    if ("post_mean" %in% names(beta_df)) {
        "post_mean"
    } else if ("mean" %in% names(beta_df)) {
        "mean"
    } else {
        NA_character_
    }
}

.get_beta_sd_col <- function(beta_df) {
    if ("post_sd" %in% names(beta_df)) {
        "post_sd"
    } else if ("sd" %in% names(beta_df)) {
        "sd"
    } else {
        NA_character_
    }
}

.summarise_beta_table <- function(beta_df) {
    if (is.null(beta_df)) {
        return(NULL)
    }
    beta_df$parameter_clean <- .standardize_beta_parameter(beta_df$parameter)
    mean_col <- .get_beta_mean_col(beta_df)
    sd_col <- .get_beta_sd_col(beta_df)
    if (is.na(mean_col) || is.na(sd_col)) {
        return(NULL)
    }
    do.call(
        rbind,
        lapply(split(beta_df, beta_df$parameter_clean), function(d) {
            data.frame(
                parameter = unique(d$parameter_clean)[1],
                n = nrow(d),
                truth_mean = mean(d$truth, na.rm = TRUE),
                mean_of_post_means = mean(d[[mean_col]], na.rm = TRUE),
                mean_bias = mean(d$bias, na.rm = TRUE),
                rmse_bias = sqrt(mean(d$bias^2, na.rm = TRUE)),
                coverage = mean(d$covered, na.rm = TRUE),
                mean_post_sd = mean(d[[sd_col]], na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        })
    )
}

.check_beta_counts <- function(beta_df, expected_reps) {
    if (is.null(beta_df) || !("parameter" %in% names(beta_df))) {
        return(list(ok = FALSE, detail = "missing beta table or parameter column"))
    }
    param <- .standardize_beta_parameter(beta_df$parameter)
    counts <- as.data.frame(table(param), stringsAsFactors = FALSE)
    names(counts) <- c("parameter", "n")

    required <- data.frame(
        parameter = c("beta0", "beta1", "beta2"),
        stringsAsFactors = FALSE
    )
    merged <- merge(required, counts, by = "parameter", all.x = TRUE)
    merged$n[is.na(merged$n)] <- 0L

    ok <- all(merged$n == expected_reps)
    list(
        ok = ok,
        detail = paste(
            paste0(merged$parameter, "=", merged$n),
            collapse = ", "
        )
    )
}

.check_dimension_table <- function(diag_df, beta_df, expected_TT, expected_n1, expected_reps) {
    checks <- data.frame()

    checks <- .add_check(
        checks, "", "dimensions", "TT equals expected value",
        .status_from_logical(!is.null(diag_df) && length(unique(diag_df$TT)) == 1L && unique(diag_df$TT) == expected_TT),
        if (!is.null(diag_df) && "TT" %in% names(diag_df)) paste(unique(diag_df$TT), collapse = ", ") else "missing TT"
    )

    checks <- .add_check(
        checks, "", "dimensions", "n1 equals expected value",
        .status_from_logical(!is.null(diag_df) && length(unique(diag_df$n1)) == 1L && unique(diag_df$n1) == expected_n1),
        if (!is.null(diag_df) && "n1" %in% names(diag_df)) paste(unique(diag_df$n1), collapse = ", ") else "missing n1"
    )

    checks <- .add_check(
        checks, "", "dimensions", "posterior diagnostics has expected number of reps",
        .status_from_logical(!is.null(diag_df) && nrow(diag_df) == expected_reps),
        if (!is.null(diag_df)) paste0("nrow=", nrow(diag_df)) else "missing diagnostics"
    )

    beta_count <- .check_beta_counts(beta_df, expected_reps)
    checks <- .add_check(
        checks, "", "dimensions", "beta0/beta1/beta2 each has expected reps",
        .status_from_logical(beta_count$ok),
        beta_count$detail
    )

    checks
}

.check_beta_reasonable <- function(beta_summary, scenario_name,
                                   max_abs_mean_bias = 0.05,
                                   max_rmse_bias = 0.08,
                                   min_coverage = 0.80,
                                   max_coverage = 1.00) {
    checks <- data.frame()
    if (is.null(beta_summary)) {
        return(.add_check(checks, scenario_name, "beta", "beta summary can be computed", "FAIL", "summary is NULL"))
    }

    for (p in c("beta0", "beta1", "beta2")) {
        d <- beta_summary[beta_summary$parameter == p, , drop = FALSE]
        if (nrow(d) != 1L) {
            checks <- .add_check(checks, scenario_name, "beta", paste0(p, " summary row exists"), "FAIL", "missing or duplicated")
            next
        }

        checks <- .add_check(
            checks, scenario_name, "beta", paste0(p, " mean bias is small"),
            .status_from_logical(abs(d$mean_bias) <= max_abs_mean_bias, warn_if_false = TRUE),
            sprintf("mean_bias=%.4f", d$mean_bias)
        )

        checks <- .add_check(
            checks, scenario_name, "beta", paste0(p, " RMSE bias is small"),
            .status_from_logical(d$rmse_bias <= max_rmse_bias, warn_if_false = TRUE),
            sprintf("rmse_bias=%.4f", d$rmse_bias)
        )

        checks <- .add_check(
            checks, scenario_name, "beta", paste0(p, " coverage is reasonable"),
            .status_from_logical(d$coverage >= min_coverage && d$coverage <= max_coverage, warn_if_false = TRUE),
            sprintf("coverage=%.3f", d$coverage)
        )
    }

    checks
}

.check_s2_lambda_reasonable <- function(lambda_df) {
    checks <- data.frame()
    if (is.null(lambda_df)) {
        return(.add_check(checks, "Scenario 2", "lambda", "lambda recovery table exists", "FAIL", "missing"))
    }

    required <- c("log_lambda_rmse", "log_lambda_coverage_95", "cor_log_lambda", "cor_delta_log_lambda")
    missing <- setdiff(required, names(lambda_df))
    if (length(missing) > 0L) {
        return(.add_check(checks, "Scenario 2", "lambda", "lambda recovery columns exist", "FAIL", paste(missing, collapse = ", ")))
    }

    mean_rmse <- mean(lambda_df$log_lambda_rmse, na.rm = TRUE)
    mean_cov <- mean(lambda_df$log_lambda_coverage_95, na.rm = TRUE)
    mean_cor <- mean(lambda_df$cor_log_lambda, na.rm = TRUE)
    mean_cor_delta <- mean(lambda_df$cor_delta_log_lambda, na.rm = TRUE)

    checks <- .add_check(
        checks, "Scenario 2", "lambda", "mean log-lambda RMSE is reasonable",
        .status_from_logical(mean_rmse <= 0.10, warn_if_false = TRUE),
        sprintf("mean=%.4f", mean_rmse)
    )

    checks <- .add_check(
        checks, "Scenario 2", "lambda", "mean log-lambda coverage is near nominal",
        .status_from_logical(mean_cov >= 0.85 && mean_cov <= 1.00, warn_if_false = TRUE),
        sprintf("mean=%.4f", mean_cov)
    )

    checks <- .add_check(
        checks, "Scenario 2", "lambda", "mean cor(log lambda) is positive and substantial",
        .status_from_logical(mean_cor >= 0.50, warn_if_false = TRUE),
        sprintf("mean=%.4f", mean_cor)
    )

    checks <- .add_check(
        checks, "Scenario 2", "lambda", "mean cor(delta log lambda) is positive",
        .status_from_logical(mean_cor_delta > 0, warn_if_false = TRUE),
        sprintf("mean=%.4f", mean_cor_delta)
    )

    checks
}

.check_s2_confounding_reasonable <- function(conf_df) {
    checks <- data.frame()
    if (is.null(conf_df)) {
        return(.add_check(checks, "Scenario 2", "confounding", "confounding table exists", "FAIL", "missing"))
    }

    required <- c("cor_x1_x2", "cor_x2_loglambda", "cor_x2_phi")
    missing <- setdiff(required, names(conf_df))
    if (length(missing) > 0L) {
        return(.add_check(checks, "Scenario 2", "confounding", "confounding columns exist", "FAIL", paste(missing, collapse = ", ")))
    }

    for (col in required) {
        max_abs <- max(abs(conf_df[[col]]), na.rm = TRUE)
        mean_abs <- mean(abs(conf_df[[col]]), na.rm = TRUE)
        checks <- .add_check(
            checks, "Scenario 2", "confounding", paste0(col, " is small"),
            .status_from_logical(max_abs <= 0.30, warn_if_false = TRUE),
            sprintf("mean_abs=%.4f, max_abs=%.4f", mean_abs, max_abs)
        )
    }

    checks
}

.check_s2_data_mode <- function(root, scenario_id, data_dir, expected_TT, expected_n1) {
    checks <- data.frame()
    data_file <- file.path(root, data_dir, scenario_id, "data_rep01.rds")
    if (!file.exists(data_file)) {
        return(.add_check(checks, "Scenario 2", "data", "data_rep01.rds exists", "FAIL", data_file))
    }

    dat <- readRDS(data_file)
    mode <- dat$x2_mode %||% NA_character_
    checks <- .add_check(
        checks, "Scenario 2", "data", "x2_mode is continuous_time",
        .status_from_logical(identical(as.character(mode), "continuous_time")),
        paste0("x2_mode=", as.character(mode))
    )

    checks <- .add_check(
        checks, "Scenario 2", "data", "x2 dimension is expected",
        .status_from_logical(!is.null(dat$x2) && all(dim(dat$x2) == c(expected_TT, expected_n1))),
        if (!is.null(dat$x2)) paste(dim(dat$x2), collapse = " x ") else "missing x2"
    )

    if (!is.null(dat$x2)) {
        col_sds <- apply(dat$x2, 2, sd)
        checks <- .add_check(
            checks, "Scenario 2", "data", "x2 varies over time within regions",
            .status_from_logical(all(col_sds > 0.05), warn_if_false = TRUE),
            paste(sprintf("%.3f", col_sds), collapse = ", ")
        )
    }

    checks
}

check_scenario1_results <- function(root = ".",
                                    scenario_id = "S1_BASELINE_VALIDATION_T100",
                                    analysis_dir = "analysis_s1_baseline_only",
                                    expected_TT = 100L,
                                    expected_n1 = 9L,
                                    expected_reps = 20L) {
    scenario_name <- "Scenario 1"
    scenario_root <- file.path(root, analysis_dir, scenario_id)
    tables_dir <- file.path(scenario_root, "tables")
    figures_dir <- file.path(scenario_root, "figures")
    figures_mcmc_dir <- file.path(scenario_root, "figures_mcmc")

    checks <- data.frame()

    required_tables <- c(
        "posterior_performance_diagnostics.csv",
        "posterior_beta_summary.csv",
        "posterior_phi_summary.csv",
        "posterior_r_summary.csv"
    )

    for (f in required_tables) {
        checks <- .add_check(
            checks, scenario_name, "files", paste0("table exists: ", f),
            .status_from_logical(.file_nonempty(file.path(tables_dir, f))),
            file.path(tables_dir, f)
        )
    }

    checks <- .add_check(checks, scenario_name, "files", "figures directory is nonempty", .status_from_logical(.dir_nonempty(figures_dir), warn_if_false = TRUE), figures_dir)
    checks <- .add_check(checks, scenario_name, "files", "figures_mcmc directory is nonempty", .status_from_logical(.dir_nonempty(figures_mcmc_dir), warn_if_false = TRUE), figures_mcmc_dir)
    checks <- .add_check(checks, scenario_name, "files", "posterior performance RDS exists", .status_from_logical(.file_nonempty(file.path(scenario_root, "scenario1_posterior_performance_results.rds")), warn_if_false = TRUE), file.path(scenario_root, "scenario1_posterior_performance_results.rds"))

    diag_df <- .read_csv_safe(file.path(tables_dir, "posterior_performance_diagnostics.csv"))
    beta_df <- .read_csv_safe(file.path(tables_dir, "posterior_beta_summary.csv"))

    dim_checks <- .check_dimension_table(diag_df, beta_df, expected_TT, expected_n1, expected_reps)
    dim_checks$scenario <- scenario_name
    checks <- rbind(checks, dim_checks)

    beta_summary <- .summarise_beta_table(beta_df)
    checks <- rbind(checks, .check_beta_reasonable(beta_summary, scenario_name))

    list(
        checks = checks,
        beta_summary = beta_summary
    )
}

check_scenario2_results <- function(root = ".",
                                    scenario_id = "S2_DYNAMIC_FIXED_GAMMA_T100",
                                    data_dir = "data_revised",
                                    analysis_dir = "analysis_s2_dynamic_fixed_gamma",
                                    expected_TT = 100L,
                                    expected_n1 = 9L,
                                    expected_reps = 20L) {
    scenario_name <- "Scenario 2"
    scenario_root <- file.path(root, analysis_dir, scenario_id)
    tables_dir <- file.path(scenario_root, "tables")
    figures_dir <- file.path(scenario_root, "figures")
    figures_mcmc_dir <- file.path(scenario_root, "figures_mcmc")

    checks <- data.frame()

    required_tables <- c(
        "posterior_performance_diagnostics.csv",
        "posterior_beta_summary.csv",
        "posterior_lambda_path_recovery.csv",
        "posterior_phi_summary.csv",
        "posterior_r_summary.csv",
        "scenario2_confounding_diagnostics.csv",
        "scenario2_data_manifest.csv"
    )

    for (f in required_tables) {
        checks <- .add_check(
            checks, scenario_name, "files", paste0("table exists: ", f),
            .status_from_logical(.file_nonempty(file.path(tables_dir, f))),
            file.path(tables_dir, f)
        )
    }

    checks <- .add_check(checks, scenario_name, "files", "figures directory is nonempty", .status_from_logical(.dir_nonempty(figures_dir), warn_if_false = TRUE), figures_dir)
    checks <- .add_check(checks, scenario_name, "files", "figures_mcmc directory is nonempty", .status_from_logical(.dir_nonempty(figures_mcmc_dir), warn_if_false = TRUE), figures_mcmc_dir)
    checks <- .add_check(checks, scenario_name, "files", "posterior performance RDS exists", .status_from_logical(.file_nonempty(file.path(scenario_root, "scenario2_posterior_performance_results.rds")), warn_if_false = TRUE), file.path(scenario_root, "scenario2_posterior_performance_results.rds"))

    diag_df <- .read_csv_safe(file.path(tables_dir, "posterior_performance_diagnostics.csv"))
    beta_df <- .read_csv_safe(file.path(tables_dir, "posterior_beta_summary.csv"))
    lambda_df <- .read_csv_safe(file.path(tables_dir, "posterior_lambda_path_recovery.csv"))
    conf_df <- .read_csv_safe(file.path(tables_dir, "scenario2_confounding_diagnostics.csv"))

    dim_checks <- .check_dimension_table(diag_df, beta_df, expected_TT, expected_n1, expected_reps)
    dim_checks$scenario <- scenario_name
    checks <- rbind(checks, dim_checks)

    beta_summary <- .summarise_beta_table(beta_df)
    checks <- rbind(checks, .check_beta_reasonable(beta_summary, scenario_name))
    checks <- rbind(checks, .check_s2_lambda_reasonable(lambda_df))
    checks <- rbind(checks, .check_s2_confounding_reasonable(conf_df))
    checks <- rbind(checks, .check_s2_data_mode(root, scenario_id, data_dir, expected_TT, expected_n1))

    lambda_summary <- NULL
    if (!is.null(lambda_df)) {
        lambda_summary <- summary(lambda_df[, intersect(c("log_lambda_rmse", "log_lambda_coverage_95", "cor_log_lambda", "cor_delta_log_lambda"), names(lambda_df)), drop = FALSE])
    }

    confounding_summary <- NULL
    if (!is.null(conf_df)) {
        confounding_summary <- summary(conf_df[, intersect(c("cor_x1_x2", "cor_x2_loglambda", "cor_x2_phi"), names(conf_df)), drop = FALSE])
    }

    list(
        checks = checks,
        beta_summary = beta_summary,
        lambda_summary = lambda_summary,
        confounding_summary = confounding_summary
    )
}

check_scenario3_results <- function(...) {
    message("Scenario 3 checks are not implemented yet. Add them once Scenario 3 output structure is finalized.")
    list(checks = data.frame())
}

check_scenario4_results <- function(...) {
    message("Scenario 4 checks are not implemented yet. Add them once Scenario 4 output structure is finalized.")
    list(checks = data.frame())
}

.print_section <- function(title) {
    cat("\n")
    cat(paste(rep("=", nchar(title) + 8L), collapse = ""), "\n", sep = "")
    cat("=== ", title, " ===\n", sep = "")
    cat(paste(rep("=", nchar(title) + 8L), collapse = ""), "\n", sep = "")
}

check_results <- function(root = ".",
                          scenarios = c("scenario1", "scenario2"),
                          stop_on_fail = FALSE,
                          expected_TT = 100L,
                          expected_n1 = 9L,
                          expected_reps = 20L) {
    cat("Using check_results.R version: 2026-05-25-clean-beta-counts\n")

    results <- list()
    all_checks <- data.frame()

    if ("scenario1" %in% scenarios) {
        .print_section("Scenario 1")
        results$scenario1 <- check_scenario1_results(
            root = root,
            expected_TT = expected_TT,
            expected_n1 = expected_n1,
            expected_reps = expected_reps
        )
        print(results$scenario1$checks[, c("scenario", "category", "check", "status", "detail")], row.names = FALSE)
        all_checks <- rbind(all_checks, results$scenario1$checks)
        cat("\nScenario 1 beta summary:\n")
        print(results$scenario1$beta_summary)
    }

    if ("scenario2" %in% scenarios) {
        .print_section("Scenario 2")
        results$scenario2 <- check_scenario2_results(
            root = root,
            expected_TT = expected_TT,
            expected_n1 = expected_n1,
            expected_reps = expected_reps
        )
        print(results$scenario2$checks[, c("scenario", "category", "check", "status", "detail")], row.names = FALSE)
        all_checks <- rbind(all_checks, results$scenario2$checks)
        cat("\nScenario 2 beta summary:\n")
        print(results$scenario2$beta_summary)
        cat("\nScenario 2 lambda summary:\n")
        print(results$scenario2$lambda_summary)
        cat("\nScenario 2 confounding summary:\n")
        print(results$scenario2$confounding_summary)
    }

    n_fail <- sum(all_checks$status == "FAIL", na.rm = TRUE)
    n_warn <- sum(all_checks$status == "WARN", na.rm = TRUE)
    overall_status <- if (n_fail > 0L) {
        "FAIL"
    } else if (n_warn > 0L) {
        "WARN"
    } else {
        "PASS"
    }

    .print_section("Overall validation status")
    cat("Overall status: ", overall_status, "\n", sep = "")
    cat("Number of FAIL: ", n_fail, "\n", sep = "")
    cat("Number of WARN: ", n_warn, "\n", sep = "")

    if (isTRUE(stop_on_fail) && n_fail > 0L) {
        stop("Validation failed. See checks with status == 'FAIL'.", call. = FALSE)
    }

    invisible(list(
        overall_status = overall_status,
        checks = all_checks,
        scenario1 = results$scenario1 %||% NULL,
        scenario2 = results$scenario2 %||% NULL
    ))
}
