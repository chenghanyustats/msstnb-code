# ================================================================
# Scenario 4 overall stress-test summary, continuous-time x2 version (v5)
# ================================================================
#
# Purpose:
#   Collect S4A--S4E posterior-performance outputs into one master
#   cross-scenario summary. This script does not fit any model. It only
#   reads the finalized scenario-level analysis tables and writes overall
#   summary tables/figures for the Scenario 4 stress-test report.
#
# Main entry point:
#   s4_overall <- run_scenario4_overall_summary_continuous_x2(root = ".")
#
# Expected output:
#   analysis_s4_overall_continuous_x2/
#     tables/
#       scenario4_cross_scenario_summary.csv
#       scenario4_stress_calibration_summary.csv
#       scenario4_regression_summary.csv
#       scenario4_latent_spatial_kappa_summary.csv
#       scenario4_numerical_stability_summary.csv
#       scenario4_replicate_level_master.csv
#       scenario4_overall_manifest.csv
#     figures/
#       scenario4_fit_status_counts.png/pdf
#       scenario4_beta2_recovery.png/pdf
#       scenario4_log_lambda_rmse.png/pdf
#       scenario4_guard_counts.png/pdf
#
# Notes:
#   - The primary performance subset is stable_plus_soft_warning when
#     available, then stable_only, then all_sampled_lambda.
#   - This matters mostly for S4A, where numerical-instability replicates
#     should be counted in fit-status summaries but not used as the main
#     stable-performance estimate.
#   - S4A tables are read from output_s4a_sparse_counts_continuous_x2,
#     because the final S4A QMD uses that finalized table location.
#
# ================================================================

.s4overall_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
}

.s4overall_dir_create <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

.s4overall_path <- function(root, ...) {
    file.path(normalizePath(root, winslash = "/", mustWork = FALSE), ...)
}

.s4overall_file_exists <- function(path) {
    !is.na(path) && nzchar(path) && file.exists(path)
}

.s4overall_first_existing <- function(paths) {
    paths <- unique(paths[!is.na(paths) & nzchar(paths)])
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0L) return(NA_character_)
    hits[[1L]]
}

.s4overall_read_optional <- function(path) {
    if (is.na(path) || !file.exists(path)) return(NULL)
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.s4overall_csv_header <- function(path) {
    out <- tryCatch(names(utils::read.csv(path, nrows = 0L, stringsAsFactors = FALSE, check.names = FALSE)), error = function(e) character(0L))
    out
}

.s4overall_rank_recursive_hit <- function(path, label) {
    p <- tolower(path)
    lab <- tolower(label)
    score <- 1000

    if (grepl("performance", lab)) {
        if (grepl("posterior_performance_by_subset|performance_by_subset", p)) score <- min(score, 1)
        if (grepl("performance_aggregate", p)) score <- min(score, 2)
        if (grepl("summary_.*all_reps|all_reps.*summary", p)) score <- min(score, 3)
        if (grepl("posterior_performance_diagnostics|replicate_level", p)) score <- min(score, 4)
    } else if (grepl("status", lab)) {
        if (grepl("fit_status_counts|status_counts", p)) score <- min(score, 1)
        if (grepl("fit_status_by_rep|status_by_rep", p)) score <- min(score, 2)
        if (grepl("summary_.*all_reps|posterior_performance_diagnostics|replicate_level", p)) score <- min(score, 3)
    } else if (grepl("beta", lab)) {
        if (grepl("beta_aggregate|beta_recovery|posterior_beta_aggregate", p)) score <- min(score, 1)
        if (grepl("posterior_beta_summary", p)) score <- min(score, 2)
        if (grepl("summary_.*all_reps|posterior_performance_diagnostics|replicate_level", p)) score <- min(score, 3)
    } else if (grepl("replicate", lab)) {
        if (grepl("replicate_level|fit_status_by_rep|posterior_performance_diagnostics", p)) score <- min(score, 1)
        if (grepl("summary_.*all_reps", p)) score <- min(score, 2)
    } else if (grepl("guard", lab)) {
        if (grepl("guard_summary|numerical_guard", p)) score <- min(score, 1)
        if (grepl("summary_.*all_reps|posterior_performance_diagnostics|replicate_level", p)) score <- min(score, 3)
    } else {
        score <- 10
    }

    # Prefer files in analysis/tables over raw output/fits when otherwise tied,
    # but keep replicate-level summaries usable as fallbacks.
    if (grepl("/analysis_", p)) score <- score - 0.10
    if (grepl("/tables/", p)) score <- score - 0.05
    if (grepl("/fits/", p)) score <- score + 0.10
    score
}

.s4overall_header_matches_label <- function(header, label) {
    if (length(header) == 0L) return(FALSE)
    h <- tolower(header)
    lab <- tolower(label)

    has_any <- function(x) any(x %in% h)
    has_grep <- function(pattern) any(grepl(pattern, h))

    if (grepl("performance", lab)) {
        return(
            has_any(c("subset")) ||
            (has_any(c("rep_id", "rep", "replicate")) && has_grep("beta1|beta_1|log_lambda|lambda_rmse|cor_log_lambda")) ||
            has_grep("beta1_mean_avg|beta2_mean_avg|log_lambda_rmse_avg")
        )
    }
    if (grepl("status", lab)) {
        return(has_any(c("fit_status", "status", "n_reps", "prop_reps")) || has_grep("guard_count|guard_total"))
    }
    if (grepl("beta", lab)) {
        return(has_any(c("parameter")) || has_grep("beta0|beta1|beta2"))
    }
    if (grepl("replicate", lab)) {
        return(has_any(c("rep_id", "rep", "replicate")) || has_grep("beta1_mean|beta2_mean|log_lambda_rmse"))
    }
    if (grepl("guard", lab)) {
        return(has_grep("guard_count|guard_total|beta_guard|kappa_guard|lambda"))
    }
    TRUE
}

.s4overall_recursive_candidates <- function(label) {
    root <- getOption("s4overall_root", default = normalizePath(".", winslash = "/", mustWork = FALSE))
    scenario <- regmatches(label, regexpr("S4[A-E]", label, ignore.case = TRUE))
    if (length(scenario) == 0L || is.na(scenario) || !nzchar(scenario)) return(character(0L))
    scenario <- toupper(scenario)

    all_csv <- tryCatch(
        list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE),
        error = function(e) character(0L)
    )
    if (length(all_csv) == 0L) return(character(0L))

    p_lower <- tolower(all_csv)
    sc_lower <- tolower(scenario)
    keep <- grepl(sc_lower, p_lower) | grepl(paste0("scenario", substring(sc_lower, 3, 3)), p_lower)
    if (scenario == "S4D") {
        keep <- keep | grepl("short.*time|short_time|scenario4d|s4d", p_lower)
        # Strictly avoid silently using the older non-continuous-x2 S4D outputs.
        # Scenario 4 overall must summarize the finalized continuous-time x2 version.
        # If your finalized S4D tables live in a folder without "continuous" in its name,
        # run find_s4d_summary_files() and add that folder explicitly to the S4D config.
        keep <- keep & grepl("continuous|continuous_x2|continuous-time", p_lower)
    }
    cand <- all_csv[keep]
    if (length(cand) == 0L) return(character(0L))

    ok <- vapply(cand, function(z) .s4overall_header_matches_label(.s4overall_csv_header(z), label), logical(1L))
    cand <- cand[ok]
    if (length(cand) == 0L) return(character(0L))

    scores <- vapply(cand, .s4overall_rank_recursive_hit, numeric(1L), label = label)
    cand[order(scores, nchar(cand), cand)]
}

.s4overall_read_required <- function(paths, label, allow_missing = FALSE) {
    paths <- unique(paths[!is.na(paths) & nzchar(paths)])
    path <- .s4overall_first_existing(paths)

    # Recursive fallback: useful for S4D, whose finalized output names may differ
    # from the later S4A/S4B/S4C/S4E naming convention.
    recursive_paths <- character(0L)
    if (is.na(path)) {
        recursive_paths <- .s4overall_recursive_candidates(label)
        path <- .s4overall_first_existing(recursive_paths)
    }

    if (is.na(path)) {
        all_tried <- unique(c(paths, recursive_paths))
        msg <- paste0(
            "Could not find ", label, ". Tried:\n  ",
            paste(all_tried, collapse = "\n  ")
        )
        if (isTRUE(allow_missing)) {
            warning(msg, call. = FALSE)
            return(list(path = NA_character_, data = NULL))
        }
        stop(msg, call. = FALSE)
    }
    list(path = path, data = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
}

.s4overall_pick_col <- function(d, candidates, default = NA_real_) {
    if (is.null(d) || nrow(d) == 0L) return(rep(default, 0L))
    nm <- names(d)
    idx <- match(candidates, nm)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(rep(default, nrow(d)))
    d[[idx[[1L]]]]
}

.s4overall_one_value <- function(d, candidates, default = NA_real_) {
    x <- .s4overall_pick_col(d, candidates, default = default)
    if (length(x) == 0L) return(default)
    x[[1L]]
}

.s4overall_num <- function(x) {
    suppressWarnings(as.numeric(x))
}

.s4overall_safe_sum <- function(x) {
    x <- .s4overall_num(x)
    if (length(x) == 0L || all(is.na(x))) return(NA_real_)
    sum(x, na.rm = TRUE)
}

.s4overall_safe_mean <- function(x) {
    x <- .s4overall_num(x)
    if (length(x) == 0L || all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.s4overall_safe_sd <- function(x) {
    x <- .s4overall_num(x)
    if (length(x) < 2L || all(is.na(x))) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

.s4overall_has_subset <- function(d) {
    !is.null(d) && nrow(d) > 0L && "subset" %in% names(d)
}

.s4overall_choose_subset_row <- function(d, subset_preference = c("stable_plus_soft_warning", "stable_only", "all_sampled_lambda")) {
    if (is.null(d) || nrow(d) == 0L) return(data.frame())
    if (!"subset" %in% names(d)) return(d[1L, , drop = FALSE])
    for (ss in subset_preference) {
        hit <- which(as.character(d$subset) == ss)
        if (length(hit) > 0L) return(d[hit[[1L]], , drop = FALSE])
    }
    d[1L, , drop = FALSE]
}

.s4overall_status_counts <- function(status_counts) {
    if (is.null(status_counts) || nrow(status_counts) == 0L) {
        return(data.frame(
            stable = NA_real_, soft_warning = NA_real_, numerical_instability = NA_real_,
            n_reps = NA_real_, prop_stable = NA_real_, stringsAsFactors = FALSE
        ))
    }

    nm <- names(status_counts)
    status_col <- c("fit_status", "status", "Status", "status_label", "fit_status_label")
    status_col <- status_col[status_col %in% nm]
    if (length(status_col) == 0L) {
        # Sometimes the first column is an unnamed status label after read.csv.
        possible <- which(vapply(status_counts, function(x) is.character(x) || is.factor(x), logical(1L)))
        if (length(possible) > 0L) status_col <- nm[possible[[1L]]]
    } else {
        status_col <- status_col[[1L]]
    }

    count_col <- c("n_reps", "n", "count", "Count", "Freq", "# reps", "X..reps")
    count_col <- count_col[count_col %in% nm]
    if (length(count_col) == 0L) {
        numeric_cols <- which(vapply(status_counts, function(x) is.numeric(x) || is.integer(x), logical(1L)))
        # avoid choosing percent/proportion first when possible
        numeric_names <- nm[numeric_cols]
        numeric_names <- numeric_names[!grepl("prop|percent", numeric_names, ignore.case = TRUE)]
        if (length(numeric_names) > 0L) count_col <- numeric_names[[1L]]
    } else {
        count_col <- count_col[[1L]]
    }

    if (length(status_col) == 0L || length(count_col) == 0L || is.na(status_col) || is.na(count_col)) {
        return(data.frame(
            stable = NA_real_, soft_warning = NA_real_, numerical_instability = NA_real_,
            n_reps = NA_real_, prop_stable = NA_real_, stringsAsFactors = FALSE
        ))
    }

    status_raw <- tolower(trimws(as.character(status_counts[[status_col]])))
    status_raw <- gsub("[ -]+", "_", status_raw)
    status_raw <- gsub("^numerical$", "numerical_instability", status_raw)
    status_raw <- gsub("^unstable$", "numerical_instability", status_raw)
    status_raw <- gsub("^warning$", "soft_warning", status_raw)
    status_raw <- gsub("^softwarning$", "soft_warning", status_raw)

    n <- .s4overall_num(status_counts[[count_col]])
    stable <- sum(n[status_raw == "stable"], na.rm = TRUE)
    soft <- sum(n[status_raw == "soft_warning"], na.rm = TRUE)
    unstable <- sum(n[status_raw == "numerical_instability"], na.rm = TRUE)

    # If there are unrecognized nonempty labels, keep total counts but do not
    # force them into stable; this makes bad status parsing visible.
    n_total <- sum(n, na.rm = TRUE)
    if (!is.finite(n_total) || n_total == 0) n_total <- NA_real_

    data.frame(
        stable = stable,
        soft_warning = soft,
        numerical_instability = unstable,
        n_reps = n_total,
        prop_stable = ifelse(is.na(n_total) || n_total == 0, NA_real_, stable / n_total),
        stringsAsFactors = FALSE
    )
}

.s4overall_add_missing_col <- function(d, col, value = NA) {
    if (is.null(d) || nrow(d) == 0L) return(d)
    if (!col %in% names(d)) d[[col]] <- value
    d
}

.s4overall_bind_rows <- function(x) {
    x <- x[!vapply(x, is.null, logical(1L))]
    x <- x[vapply(x, nrow, integer(1L)) > 0L]
    if (length(x) == 0L) return(data.frame())
    all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    x2 <- lapply(x, function(d) {
        miss <- setdiff(all_names, names(d))
        for (m in miss) d[[m]] <- NA
        d[, all_names, drop = FALSE]
    })
    do.call(rbind, x2)
}

.s4overall_fig_device <- function(path, width = 8, height = 5) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "pdf") {
        grDevices::pdf(path, width = width, height = height)
    } else {
        grDevices::png(path, width = width, height = height, units = "in", res = 300)
    }
}

.s4overall_close_device <- function() {
    grDevices::dev.off()
}

.s4overall_plot_bar <- function(path, names, values, main, ylab, width = 8, height = 5, las = 2) {
    # Robust plotting helper.  Earlier versions failed when a metric was all
    # missing for one source table, or when a metric was negative-only (e.g.,
    # beta2 posterior means).  This helper keeps the script from stopping at
    # the figure stage while preserving the table outputs.
    vals <- .s4overall_num(values)
    labs <- as.character(names)
    if (length(vals) == 0L) {
        vals <- 0
        labs <- "No data"
    }
    if (length(labs) != length(vals)) labs <- seq_along(vals)

    finite <- is.finite(vals)
    plot_vals <- vals
    plot_vals[!finite] <- 0

    y_min <- min(0, plot_vals, na.rm = TRUE)
    y_max <- max(0, plot_vals, na.rm = TRUE)
    if (!is.finite(y_min)) y_min <- 0
    if (!is.finite(y_max)) y_max <- 1
    if (abs(y_max - y_min) < .Machine$double.eps^0.5) {
        pad <- ifelse(abs(y_max) < .Machine$double.eps^0.5, 1, abs(y_max) * 0.15)
        y_min <- y_min - pad
        y_max <- y_max + pad
    } else {
        pad <- 0.10 * (y_max - y_min)
        y_min <- y_min - pad
        y_max <- y_max + pad
    }

    .s4overall_fig_device(path, width = width, height = height)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({ graphics::par(oldpar); .s4overall_close_device() }, add = TRUE)
    graphics::par(mar = c(8, 5, 4, 2) + 0.1)
    graphics::barplot(
        height = plot_vals,
        names.arg = labs,
        main = main,
        ylab = ylab,
        las = las,
        ylim = c(y_min, y_max)
    )
    if (any(!finite)) {
        graphics::mtext("Missing values plotted as 0", side = 1, line = 6.8, cex = 0.7)
    }
}

.s4overall_plot_status <- function(path, stability_summary, width = 8, height = 5) {
    .s4overall_fig_device(path, width = width, height = height)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit({ graphics::par(oldpar); .s4overall_close_device() }, add = TRUE)
    graphics::par(mar = c(7, 5, 4, 2) + 0.1)
    m <- rbind(
        stable = .s4overall_num(stability_summary$stable),
        soft_warning = .s4overall_num(stability_summary$soft_warning),
        numerical_instability = .s4overall_num(stability_summary$numerical_instability)
    )
    m[!is.finite(m)] <- 0
    if (ncol(m) == 0L || nrow(m) == 0L) {
        m <- matrix(0, nrow = 3, ncol = 1,
                    dimnames = list(c("stable", "soft_warning", "numerical_instability"), "No data"))
        labs <- "No data"
    } else {
        labs <- as.character(stability_summary$scenario)
        if (length(labs) != ncol(m)) labs <- seq_len(ncol(m))
    }
    ylim_max <- max(colSums(m), na.rm = TRUE)
    if (!is.finite(ylim_max) || ylim_max <= 0) ylim_max <- 1
    graphics::barplot(
        m,
        beside = FALSE,
        names.arg = labs,
        las = 2,
        ylab = "Number of replicates",
        main = "Scenario 4 fit-status counts",
        ylim = c(0, ylim_max * 1.15)
    )
    graphics::legend("topright", legend = rownames(m), fill = grDevices::gray.colors(3), bty = "n")
}

.s4overall_default_config <- function(root) {
    p <- function(...) .s4overall_path(root, ...)
    list(
        S4A = list(
            scenario = "S4A",
            label = "Sparse counts",
            stress_mechanism = "Sparse-count observation regime",
            interpretation = "Sparse counts are the clearest latent-scale numerical-instability stress; stable replicates recover slopes but unstable replicates expose lambda-scale failure.",
            table_dirs = c(
                p("output_s4a_sparse_counts_continuous_x2", "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100", "tables"),
                p("analysis_s4a_sparse_counts_continuous_x2", "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100", "tables")
            ),
            status_files = c("s4a_fit_status_counts_continuous_x2.csv"),
            perf_files = c("posterior_performance_by_subset_continuous_x2.csv"),
            beta_files = c("posterior_beta_aggregate_continuous_x2.csv", "posterior_beta_summary_continuous_x2.csv"),
            rep_files = c("s4a_replicate_level_summary_continuous_x2.csv"),
            guard_files = c("numerical_guard_summary_by_subset_continuous_x2.csv"),
            special_files = list(
                lambda_subset = c("lambda_recovery_by_subset_continuous_x2.csv")
            )
        ),
        S4B = list(
            scenario = "S4B",
            label = "Low/heterogeneous exposure",
            stress_mechanism = "Low and heterogeneous exposure; local information imbalance",
            interpretation = "Low exposure mainly weakens low-exposure latent-risk recovery while regression slopes remain stable.",
            table_dirs = c(
                p("analysis_s4b_low_exposure_continuous_x2", "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables"),
                p("output_s4b_low_exposure_continuous_x2", "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables")
            ),
            status_files = c("scenario4b_fit_status_counts.csv", "scenario4b_fit_status_counts_continuous_x2.csv"),
            perf_files = c("scenario4b_performance_by_subset.csv", "posterior_performance_by_subset_continuous_x2.csv"),
            beta_files = c("scenario4b_beta_recovery_by_subset.csv", "posterior_beta_aggregate_continuous_x2.csv"),
            rep_files = c("posterior_performance_diagnostics.csv", "scenario4b_fit_status_by_rep.csv"),
            guard_files = c("scenario4b_guard_summary_by_subset.csv", "numerical_guard_summary_by_subset_continuous_x2.csv"),
            special_files = list(
                low_high = c("scenario4b_low_high_exposure_contrast_by_subset.csv"),
                exposure_group = c("scenario4b_exposure_group_recovery_by_subset.csv")
            )
        ),
        S4C = list(
            scenario = "S4C",
            label = "Strong overdispersion",
            stress_mechanism = "Small-r negative-binomial overdispersion",
            interpretation = "Strong overdispersion is manageable when counts remain abundant; r and kappa-level variability are recovered well.",
            table_dirs = c(
                p("analysis_s4c_overdispersion_continuous_x2", "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables"),
                p("output_s4c_small_r_fixed_gamma_continuous_x2", "S4C_STRONG_OVERDISPERSION_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables")
            ),
            status_files = c("s4c_fit_status_counts_continuous_x2.csv"),
            perf_files = c("posterior_performance_by_subset_continuous_x2.csv"),
            beta_files = c("posterior_beta_aggregate_continuous_x2.csv", "posterior_beta_summary_continuous_x2.csv"),
            rep_files = c("s4c_replicate_level_summary_continuous_x2.csv"),
            guard_files = c("numerical_guard_summary_by_subset_continuous_x2.csv"),
            special_files = list(
                r = c("r_recovery_by_subset_continuous_x2.csv"),
                kappa = c("kappa_recovery_by_subset_continuous_x2.csv")
            )
        ),
        S4D = list(
            scenario = "S4D",
            label = "Short time series",
            stress_mechanism = "Short temporal record, T = 30",
            interpretation = "Short T preserves slope recovery but weakens temporal-path recovery relative to T = 100 scenarios.",
            table_dirs = c(
                p("analysis_s4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "tables"),
                p("analysis_s4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_CONTINUOUS_X2_T30", "tables"),
                p("analysis_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "tables"),
                p("analysis_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_SERIES_CONTINUOUS_X2_T30", "tables"),
                p("analysis_scenario4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "tables"),
                p("analysis_scenario4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_CONTINUOUS_X2_T30", "tables"),
                p("output_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "tables"),
                p("output_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_SERIES_CONTINUOUS_X2_T30", "tables"),
                p("output_s4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "tables"),
                p("output_s4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_CONTINUOUS_X2_T30", "tables"),
                p("output_s4d_short_time_series_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "fits"),
                p("output_s4d_short_time_fixed_gamma_continuous_x2", "S4D_SHORT_TIME_SERIES_FIXED_GAMMA_CONTINUOUS_X2_T30", "fits")
            ),
            status_files = c("s4d_fit_status_counts_continuous_x2.csv", "scenario4d_fit_status_counts_continuous_x2.csv", "scenario4d_fit_status_counts.csv", "fit_status_counts.csv", "status_counts.csv"),
            perf_files = c("posterior_performance_by_subset_continuous_x2.csv", "scenario4d_performance_by_subset.csv", "posterior_performance_by_subset.csv", "performance_by_subset.csv", "summary_S4D_short_time_fixed_gamma_continuous_x2_all_reps.csv", "summary_S4D_short_time_series_continuous_x2_all_reps.csv", "summary_S4D_all_reps.csv"),
            beta_files = c("posterior_beta_aggregate_continuous_x2.csv", "scenario4d_beta_recovery_by_subset.csv", "posterior_beta_summary_continuous_x2.csv", "posterior_beta_summary.csv"),
            rep_files = c("s4d_replicate_level_summary_continuous_x2.csv", "scenario4d_replicate_level_summary_continuous_x2.csv", "posterior_performance_diagnostics.csv", "summary_S4D_short_time_fixed_gamma_continuous_x2_all_reps.csv", "summary_S4D_short_time_series_continuous_x2_all_reps.csv", "summary_S4D_all_reps.csv"),
            guard_files = c("numerical_guard_summary_by_subset_continuous_x2.csv", "scenario4d_guard_summary_by_subset.csv", "guard_summary_by_subset.csv"),
            special_files = list()
        ),
        S4E = list(
            scenario = "S4E",
            label = "Spatial/covariate confounding",
            stress_mechanism = "Strong x2-phi spatial/covariate confounding",
            interpretation = "Strong x2-phi confounding is successfully calibrated and remains stable in this high-count fixed-gamma setting.",
            table_dirs = c(
                p("analysis_s4e_spatial_covariate_confounding_continuous_x2", "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables"),
                p("output_s4e_spatial_covariate_confounding_continuous_x2", "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_FIXED_GAMMA_CONTINUOUS_X2_T100", "tables")
            ),
            status_files = c("s4e_fit_status_counts_continuous_x2.csv"),
            perf_files = c("posterior_performance_by_subset_continuous_x2.csv"),
            beta_files = c("posterior_beta_aggregate_continuous_x2.csv", "posterior_beta_summary_continuous_x2.csv"),
            rep_files = c("s4e_replicate_level_summary_continuous_x2.csv"),
            guard_files = c("numerical_guard_summary_by_subset_continuous_x2.csv"),
            special_files = list(
                confounding = c("confounding_calibration_by_subset_continuous_x2.csv"),
                beta_confounding = c("beta_confounding_diagnostics_continuous_x2.csv")
            )
        )
    )
}

.s4overall_candidate_paths <- function(cfg, file_names) {
    # Use an explicit Cartesian product.  Do not use outer(file.path), because
    # file.path is vectorized and can recycle arguments in a way that silently
    # drops many directory/file combinations.
    dirs <- unique(cfg$table_dirs[!is.na(cfg$table_dirs) & nzchar(cfg$table_dirs)])
    files <- unique(file_names[!is.na(file_names) & nzchar(file_names)])
    if (length(dirs) == 0L || length(files) == 0L) return(character(0L))
    grid <- expand.grid(dir = dirs, file = files, stringsAsFactors = FALSE)
    file.path(grid$dir, grid$file)
}


.s4overall_rep_col <- function(d, candidates) {
    nm <- names(d)
    idx <- match(candidates, nm)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(rep(NA_real_, nrow(d)))
    d[[idx[[1L]]]]
}

.s4overall_is_replicate_level <- function(d) {
    if (is.null(d) || nrow(d) <= 1L) return(FALSE)
    nm <- names(d)
    any(c("rep_id", "rep", "replicate") %in% nm) ||
        any(grepl("^rep", nm, ignore.case = TRUE)) &&
        any(grepl("beta1|beta2|log_lambda|lambda_rmse", nm, ignore.case = TRUE))
}

.s4overall_mean_col <- function(d, candidates) .s4overall_safe_mean(.s4overall_rep_col(d, candidates))
.s4overall_sd_col <- function(d, candidates) .s4overall_safe_sd(.s4overall_rep_col(d, candidates))
.s4overall_sum_col <- function(d, candidates) .s4overall_safe_sum(.s4overall_rep_col(d, candidates))
.s4overall_max_col <- function(d, candidates) {
    x <- .s4overall_num(.s4overall_rep_col(d, candidates))
    if (length(x) == 0L || all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
}

.s4overall_aggregate_replicate_to_performance <- function(d, scenario = NA_character_) {
    if (is.null(d) || nrow(d) == 0L) return(d)
    rep_id <- .s4overall_rep_col(d, c("rep_id", "rep", "replicate"))
    if (all(is.na(rep_id))) rep_id <- seq_len(nrow(d))

    status <- NULL
    if ("fit_status" %in% names(d)) status <- as.character(d$fit_status)
    if (is.null(status)) {
        guard_cols <- grep("guard_count$|_guard_count$|guard_total$", names(d), value = TRUE)
        if (length(guard_cols) > 0L) {
            guards <- d[, guard_cols, drop = FALSE]
            guards[] <- lapply(guards, function(x) suppressWarnings(as.numeric(x)))
            status <- ifelse(rowSums(guards, na.rm = TRUE) == 0, "stable", "soft_warning")
        } else {
            status <- rep("stable", nrow(d))
        }
    }

    make_row <- function(row_index, subset_name) {
        dd <- d[row_index, , drop = FALSE]
        data.frame(
            subset = subset_name,
            n_reps = nrow(dd),
            reps = paste(rep_id[row_index], collapse = ","),
            mean_count_avg = .s4overall_mean_col(dd, c("mean_count_avg", "observed_mean_count", "observed_mean_count_avg", "mean_count")),
            zero_prop_avg = .s4overall_mean_col(dd, c("zero_prop_avg", "observed_zero_prop", "observed_zero_prop_avg", "zero_prop")),
            total_count_sum = .s4overall_sum_col(dd, c("total_count_sum", "observed_total_count", "total_count")),
            max_count_max = .s4overall_max_col(dd, c("max_count_max", "max_count", "observed_max_count")),
            beta0_mean_avg = .s4overall_mean_col(dd, c("beta0_mean", "beta0_mean_avg")),
            beta1_mean_avg = .s4overall_mean_col(dd, c("beta1_mean", "beta1_mean_avg")),
            beta1_mean_sd = .s4overall_sd_col(dd, c("beta1_mean", "beta1_mean_avg")),
            beta2_mean_avg = .s4overall_mean_col(dd, c("beta2_mean", "beta2_mean_avg")),
            beta2_mean_sd = .s4overall_sd_col(dd, c("beta2_mean", "beta2_mean_avg")),
            beta0_coverage = .s4overall_mean_col(dd, c("beta0_coverage", "beta0_covered", "beta0_covered_avg", "covered_beta0")),
            beta1_coverage = .s4overall_mean_col(dd, c("beta1_coverage", "beta1_covered", "beta1_covered_avg", "covered_beta1")),
            beta2_coverage = .s4overall_mean_col(dd, c("beta2_coverage", "beta2_covered", "beta2_covered_avg", "covered_beta2")),
            r_truth_avg = .s4overall_mean_col(dd, c("r_truth", "r_truth_avg", "r_stress_truth", "r_reference_truth")),
            r_mean_avg = .s4overall_mean_col(dd, c("r_mean", "r_mean_avg")),
            r_mean_sd = .s4overall_sd_col(dd, c("r_mean", "r_mean_avg")),
            phi_rmse_avg = .s4overall_mean_col(dd, c("phi_rmse", "phi_rmse_avg")),
            phi_cor_avg = .s4overall_mean_col(dd, c("phi_cor", "phi_cor_avg")),
            lambda_rmse_avg = .s4overall_mean_col(dd, c("lambda_rmse", "lambda_rmse_avg")),
            log_lambda_rmse_avg = .s4overall_mean_col(dd, c("log_lambda_rmse", "log_lambda_rmse_avg")),
            cor_log_lambda_avg = .s4overall_mean_col(dd, c("cor_log_lambda", "cor_log_lambda_avg")),
            kappa_truth_cv_avg = .s4overall_mean_col(dd, c("kappa_truth_cv", "kappa_truth_cv_avg", "truth_cv", "truth_CV")),
            kappa_post_mean_cv_avg = .s4overall_mean_col(dd, c("kappa_post_mean_cv", "kappa_post_mean_cv_avg", "post_mean_cv", "post_mean_CV")),
            log_kappa_rmse_avg = .s4overall_mean_col(dd, c("log_kappa_rmse", "log_kappa_rmse_avg")),
            cor_log_kappa_avg = .s4overall_mean_col(dd, c("cor_log_kappa", "cor_log_kappa_avg")),
            beta_guard_total = .s4overall_sum_col(dd, c("beta_guard_total", "beta_guard_count", "s4a_beta_guard_count", "s4b_beta_guard_count", "s4c_beta_guard_count", "s4d_beta_guard_count", "s4e_beta_guard_count")),
            kappa_guard_total = .s4overall_sum_col(dd, c("kappa_guard_total", "kappa_guard_count", "s4a_kappa_guard_count", "s4b_kappa_guard_count", "s4c_kappa_guard_count", "s4d_kappa_guard_count", "s4e_kappa_guard_count")),
            lambda_input_guard_total = .s4overall_sum_col(dd, c("lambda_input_guard_total", "lambda_input_guard_count", "s4a_lambda_input_guard_count", "s4b_lambda_input_guard_count", "s4c_lambda_input_guard_count", "s4d_lambda_input_guard_count", "s4e_lambda_input_guard_count")),
            lambda_output_guard_total = .s4overall_sum_col(dd, c("lambda_output_guard_total", "lambda_output_guard_count", "s4a_lambda_output_guard_count", "s4b_lambda_output_guard_count", "s4c_lambda_output_guard_count", "s4d_lambda_output_guard_count", "s4e_lambda_output_guard_count")),
            stringsAsFactors = FALSE
        )
    }

    all_idx <- seq_len(nrow(d))
    stable_idx <- which(status == "stable")
    stable_soft_idx <- which(status %in% c("stable", "soft_warning"))
    out <- make_row(all_idx, "all_sampled_lambda")
    if (length(stable_idx) > 0L) out <- rbind(out, make_row(stable_idx, "stable_only"))
    if (length(stable_soft_idx) > 0L) out <- rbind(out, make_row(stable_soft_idx, "stable_plus_soft_warning"))
    out
}

.s4overall_status_from_replicate_or_perf <- function(replicate_data = NULL, performance_row = NULL) {
    # Prefer replicate-level fit_status when available.
    if (!is.null(replicate_data) && nrow(replicate_data) > 0L && "fit_status" %in% names(replicate_data)) {
        tab <- as.data.frame(table(as.character(replicate_data$fit_status)), stringsAsFactors = FALSE)
        names(tab) <- c("fit_status", "n_reps")
        tab$prop_reps <- tab$n_reps / sum(tab$n_reps)
        return(tab)
    }

    # If replicate-level guards are present but no status column exists, classify
    # stable when all available guard counts are zero.
    if (!is.null(replicate_data) && nrow(replicate_data) > 0L) {
        guard_cols <- grep("guard_count$|_guard_count$|guard_total$", names(replicate_data), value = TRUE)
        if (length(guard_cols) > 0L) {
            guards <- replicate_data[, guard_cols, drop = FALSE]
            guards[] <- lapply(guards, function(x) suppressWarnings(as.numeric(x)))
            row_guard_sum <- rowSums(guards, na.rm = TRUE)
            status <- ifelse(row_guard_sum == 0, "stable", "soft_warning")
            tab <- as.data.frame(table(status), stringsAsFactors = FALSE)
            names(tab) <- c("fit_status", "n_reps")
            tab$prop_reps <- tab$n_reps / sum(tab$n_reps)
            return(tab)
        }
    }

    # Last fallback: use n_reps from the selected performance row and assume stable.
    if (!is.null(performance_row) && nrow(performance_row) > 0L) {
        n_reps <- .s4overall_one_value(performance_row, c("n_reps", "n_reps_total"), default = NA_real_)
        if (!is.na(n_reps)) {
            return(data.frame(
                fit_status = "stable",
                n_reps = n_reps,
                prop_reps = 1,
                stringsAsFactors = FALSE
            ))
        }
    }

    data.frame(fit_status = character(0L), n_reps = numeric(0L), prop_reps = numeric(0L))
}

.s4overall_read_scenario <- function(cfg, subset_preference, allow_missing = FALSE, verbose = TRUE) {
    .s4overall_msg("Reading ", cfg$scenario, ": ", cfg$label, verbose = verbose)

    status_obj <- .s4overall_read_required(
        .s4overall_candidate_paths(cfg, cfg$status_files),
        paste0(cfg$scenario, " status-count table"),
        allow_missing = TRUE
    )
    perf_obj <- .s4overall_read_required(
        .s4overall_candidate_paths(cfg, cfg$perf_files),
        paste0(cfg$scenario, " performance aggregate table"),
        allow_missing = allow_missing
    )

    if (is.null(perf_obj$data)) return(NULL)

    perf_raw_replicate <- NULL
    if (.s4overall_is_replicate_level(perf_obj$data)) {
        perf_raw_replicate <- perf_obj$data
        perf_obj$data <- .s4overall_aggregate_replicate_to_performance(perf_obj$data, cfg$scenario)
    }

    beta_obj <- .s4overall_read_required(
        .s4overall_candidate_paths(cfg, cfg$beta_files),
        paste0(cfg$scenario, " beta table"),
        allow_missing = TRUE
    )
    rep_obj <- .s4overall_read_required(
        .s4overall_candidate_paths(cfg, cfg$rep_files),
        paste0(cfg$scenario, " replicate-level table"),
        allow_missing = TRUE
    )
    guard_obj <- .s4overall_read_required(
        .s4overall_candidate_paths(cfg, cfg$guard_files),
        paste0(cfg$scenario, " numerical-guard table"),
        allow_missing = TRUE
    )

    if ((is.null(rep_obj$data) || nrow(rep_obj$data) == 0L) && !is.null(perf_raw_replicate)) {
        rep_obj$data <- perf_raw_replicate
        rep_obj$path <- paste0(perf_obj$path, " [replicate-level source used to derive performance aggregate]")
    }

    special <- list()
    if (!is.null(cfg$special_files) && length(cfg$special_files) > 0L) {
        for (nm in names(cfg$special_files)) {
            obj <- .s4overall_read_required(
                .s4overall_candidate_paths(cfg, cfg$special_files[[nm]]),
                paste0(cfg$scenario, " special table: ", nm),
                allow_missing = TRUE
            )
            special[[nm]] <- obj
        }
    }

    perf_row <- .s4overall_choose_subset_row(perf_obj$data, subset_preference)
    primary_subset <- if ("subset" %in% names(perf_row)) as.character(perf_row$subset[[1L]]) else "first_row"

    if (is.null(status_obj$data) || nrow(status_obj$data) == 0L) {
        status_obj$data <- .s4overall_status_from_replicate_or_perf(rep_obj$data, perf_row)
        status_obj$path <- "derived_from_replicate_or_performance"
    }
    status <- .s4overall_status_counts(status_obj$data)

    # Regression summaries may be either aggregate-by-parameter or replicate-level long form.
    beta_long_summary <- .s4overall_make_regression_summary(cfg, beta_obj$data, perf_row, primary_subset)

    list(
        cfg = cfg,
        status = status,
        status_counts = status_obj$data,
        performance = perf_obj$data,
        performance_row = perf_row,
        beta = beta_obj$data,
        beta_summary = beta_long_summary,
        replicate = rep_obj$data,
        guard = guard_obj$data,
        special = special,
        manifest = data.frame(
            scenario = cfg$scenario,
            label = cfg$label,
            status_path = status_obj$path,
            performance_path = perf_obj$path,
            beta_path = beta_obj$path,
            replicate_path = rep_obj$path,
            guard_path = guard_obj$path,
            primary_subset = primary_subset,
            stringsAsFactors = FALSE
        )
    )
}

.s4overall_make_regression_summary <- function(cfg, beta_data, perf_row, primary_subset) {
    scenario <- cfg$scenario
    label <- cfg$label
    pars <- c("beta0", "beta1", "beta2")

    # Helper: fallback directly from performance aggregate row.
    make_from_perf <- function() {
        out <- data.frame()
        for (pp in pars) {
            idx <- sub("beta", "", pp)
            out <- rbind(out, data.frame(
                scenario = scenario,
                label = label,
                subset = primary_subset,
                parameter = pp,
                truth = .s4overall_one_value(perf_row, c(paste0(pp, "_truth_avg"), paste0(pp, "_truth"), paste0("beta", idx, "_truth_avg"), paste0("beta", idx, "_truth"))),
                posterior_mean = .s4overall_one_value(perf_row, c(paste0(pp, "_mean_avg"), paste0(pp, "_mean"), paste0("beta", idx, "_mean_avg"), paste0("beta", idx, "_mean"))),
                across_rep_sd = .s4overall_one_value(perf_row, c(paste0(pp, "_mean_sd"), paste0(pp, "_sd"))),
                mean_bias = .s4overall_one_value(perf_row, c(paste0(pp, "_bias_avg"), paste0(pp, "_bias"))),
                mean_abs_bias = .s4overall_one_value(perf_row, c(paste0(pp, "_abs_bias_avg"), paste0(pp, "_abs_bias"))),
                rmse_bias = NA_real_,
                coverage = .s4overall_one_value(perf_row, c(paste0(pp, "_coverage"), paste0(pp, "_covered_avg"))),
                mean_post_sd = NA_real_,
                stringsAsFactors = FALSE
            ))
        }
        out
    }

    if (!is.null(beta_data) && nrow(beta_data) > 0L && "parameter" %in% names(beta_data)) {
        d0 <- beta_data

        # If this is a replicate-level long table, aggregate instead of taking
        # the first row.  This prevents S4A from being summarized by rep 1 only.
        is_replicate_long <- any(c("rep_id", "rep", "replicate") %in% names(d0)) ||
            (any(duplicated(as.character(d0$parameter))) &&
                !any(c("mean_post_mean", "sd_post_mean_across_reps", "mean_bias", "mean_abs_bias") %in% names(d0)))

        if (is_replicate_long) {
            d <- d0
            if ("subset" %in% names(d)) {
                d <- d[as.character(d$subset) == primary_subset, , drop = FALSE]
            } else if ("fit_status" %in% names(d) && primary_subset == "stable_plus_soft_warning") {
                d <- d[as.character(d$fit_status) %in% c("stable", "soft_warning"), , drop = FALSE]
            } else if ("fit_status" %in% names(d) && primary_subset == "stable_only") {
                d <- d[as.character(d$fit_status) == "stable", , drop = FALSE]
            }
            if (nrow(d) == 0L) d <- d0

            out <- data.frame()
            for (pp in pars) {
                dd <- d[as.character(d$parameter) == pp, , drop = FALSE]
                if (nrow(dd) == 0L) next
                truth <- .s4overall_one_value(dd, c("truth", "mean_truth", "parameter_truth"), default = NA_real_)
                post_col <- .s4overall_pick_col(dd, c("post_mean", "posterior_mean", "mean", "mean_post_mean"), default = NA_real_)
                sd_col <- .s4overall_pick_col(dd, c("post_sd", "posterior_sd", "sd", "mean_post_sd"), default = NA_real_)
                covered_col <- .s4overall_pick_col(dd, c("covered", "coverage", "covered_95", "covered_beta"), default = NA_real_)
                post_mean <- .s4overall_safe_mean(post_col)
                bias <- post_mean - truth
                out <- rbind(out, data.frame(
                    scenario = scenario,
                    label = label,
                    subset = primary_subset,
                    parameter = pp,
                    truth = truth,
                    posterior_mean = post_mean,
                    across_rep_sd = .s4overall_safe_sd(post_col),
                    mean_bias = bias,
                    mean_abs_bias = .s4overall_safe_mean(abs(.s4overall_num(post_col) - truth)),
                    rmse_bias = sqrt(.s4overall_safe_mean((.s4overall_num(post_col) - truth)^2)),
                    coverage = .s4overall_safe_mean(covered_col),
                    mean_post_sd = .s4overall_safe_mean(sd_col),
                    stringsAsFactors = FALSE
                ))
            }
            if (nrow(out) > 0L) return(out)
            return(make_from_perf())
        }

        d <- d0
        if ("subset" %in% names(d)) {
            d <- d[as.character(d$subset) == primary_subset, , drop = FALSE]
            if (nrow(d) == 0L) d <- d0[as.character(d0$subset) %in% c("stable_plus_soft_warning", "stable_only", "all_sampled_lambda"), , drop = FALSE]
        }
        out <- data.frame()
        for (pp in pars) {
            dd <- d[as.character(d$parameter) == pp, , drop = FALSE]
            if (nrow(dd) == 0L) next
            dd <- dd[1L, , drop = FALSE]
            out <- rbind(out, data.frame(
                scenario = scenario,
                label = label,
                subset = primary_subset,
                parameter = pp,
                truth = .s4overall_one_value(dd, c("mean_truth", "truth", "parameter_truth")),
                posterior_mean = .s4overall_one_value(dd, c("mean_post_mean", "post_mean", "posterior_mean", "mean")),
                across_rep_sd = .s4overall_one_value(dd, c("sd_post_mean_across_reps", "post_mean_sd", "sd")),
                mean_bias = .s4overall_one_value(dd, c("mean_bias", "bias", "mean_post_bias")),
                mean_abs_bias = .s4overall_one_value(dd, c("mean_abs_bias", "mean_absolute_bias", "abs_bias")),
                rmse_bias = .s4overall_one_value(dd, c("rmse_bias", "bias_rmse")),
                coverage = .s4overall_one_value(dd, c("coverage", "covered", "coverage_95")),
                mean_post_sd = .s4overall_one_value(dd, c("mean_post_sd", "posterior_sd", "post_sd")),
                stringsAsFactors = FALSE
            ))
        }
        if (nrow(out) > 0L) return(out)
    }

    make_from_perf()
}

.s4overall_special_values <- function(scenario_obj, subset_preference) {
    cfg <- scenario_obj$cfg
    perf <- scenario_obj$performance_row
    special <- scenario_obj$special
    out <- data.frame(
        special_metric_1 = NA_character_, special_value_1 = NA_real_,
        special_metric_2 = NA_character_, special_value_2 = NA_real_,
        special_metric_3 = NA_character_, special_value_3 = NA_real_,
        stringsAsFactors = FALSE
    )

    sc <- cfg$scenario
    if (sc == "S4A") {
        out$special_metric_1 <- "stable_reps"
        out$special_value_1 <- scenario_obj$status$stable
        out$special_metric_2 <- "numerical_instability_reps"
        out$special_value_2 <- scenario_obj$status$numerical_instability
        out$special_metric_3 <- "lambda_output_guard_total"
        out$special_value_3 <- .s4overall_one_value(perf, c("lambda_output_guard_total", "s4a_lambda_output_guard_total"))
    } else if (sc == "S4B") {
        low_high <- NULL
        if (!is.null(special$low_high)) low_high <- .s4overall_choose_subset_row(special$low_high$data, subset_preference)
        source_row <- if (!is.null(low_high) && nrow(low_high) > 0L) low_high else perf
        out$special_metric_1 <- "low_exp_log_lambda_rmse"
        out$special_value_1 <- .s4overall_one_value(source_row, c("low_log_lambda_rmse", "low_exp_log_lambda_rmse_avg", "lambda_low_exposure_log_lambda_rmse_avg"))
        out$special_metric_2 <- "high_exp_log_lambda_rmse"
        out$special_value_2 <- .s4overall_one_value(source_row, c("high_log_lambda_rmse", "high_exp_log_lambda_rmse_avg", "lambda_high_exposure_log_lambda_rmse_avg"))
        out$special_metric_3 <- "low_minus_high_log_lambda_rmse"
        out$special_value_3 <- .s4overall_one_value(source_row, c("low_minus_high_log_lambda_rmse", "low_high_log_lambda_rmse_gap", "low_exp_minus_high_exp_log_lambda_rmse_avg"))
        if (is.na(out$special_value_3) && !is.na(out$special_value_1) && !is.na(out$special_value_2)) {
            out$special_value_3 <- out$special_value_1 - out$special_value_2
        }
    } else if (sc == "S4C") {
        kappa <- NULL
        if (!is.null(special$kappa)) kappa <- .s4overall_choose_subset_row(special$kappa$data, subset_preference)
        source_row <- if (!is.null(kappa) && nrow(kappa) > 0L) kappa else perf
        out$special_metric_1 <- "kappa_truth_cv"
        out$special_value_1 <- .s4overall_one_value(source_row, c("kappa_truth_cv_avg", "truth_CV", "truth_cv"))
        out$special_metric_2 <- "reference_kappa_cv"
        out$special_value_2 <- .s4overall_one_value(source_row, c("reference_kappa_cv_avg", "reference_cv"))
        out$special_metric_3 <- "kappa_cv_increase"
        out$special_value_3 <- .s4overall_one_value(source_row, c("kappa_cv_increase_avg", "cv_increase"))
    } else if (sc == "S4D") {
        out$special_metric_1 <- "time_points"
        out$special_value_1 <- .s4overall_one_value(perf, c("TT", "TT_unique", "time_points", "T"), default = 30)
        out$special_metric_2 <- "short_T_log_lambda_rmse"
        out$special_value_2 <- .s4overall_one_value(perf, c("log_lambda_rmse_avg", "mean_log_lambda_rmse"))
        out$special_metric_3 <- "short_T_cor_log_lambda"
        out$special_value_3 <- .s4overall_one_value(perf, c("cor_log_lambda_avg", "mean_cor_log_lambda"))
    } else if (sc == "S4E") {
        conf <- NULL
        if (!is.null(special$confounding)) conf <- .s4overall_choose_subset_row(special$confounding$data, subset_preference)
        source_row <- if (!is.null(conf) && nrow(conf) > 0L) conf else perf
        out$special_metric_1 <- "abs_cell_cor_x2_phi"
        out$special_value_1 <- .s4overall_one_value(source_row, c("x2_phi_abs_cell_cor_avg", "abs_cell_cor_avg"))
        out$special_metric_2 <- "abs_area_mean_cor_x2_phi"
        out$special_value_2 <- .s4overall_one_value(source_row, c("x2_phi_abs_area_mean_cor_avg", "abs_area_cor_avg"))
        out$special_metric_3 <- "confounded_minus_reference_cell_cor"
        out$special_value_3 <- .s4overall_one_value(source_row, c("confounded_minus_reference_cell_cor", "cell_cor_increase"))
    }
    out
}

.s4overall_make_summary_row <- function(scenario_obj, subset_preference) {
    cfg <- scenario_obj$cfg
    perf <- scenario_obj$performance_row
    status <- scenario_obj$status
    special <- .s4overall_special_values(scenario_obj, subset_preference)

    data.frame(
        scenario = cfg$scenario,
        label = cfg$label,
        stress_mechanism = cfg$stress_mechanism,
        primary_subset = if ("subset" %in% names(perf)) as.character(perf$subset[[1L]]) else "first_row",
        n_reps_total = status$n_reps,
        n_stable = status$stable,
        n_soft_warning = status$soft_warning,
        n_numerical_instability = status$numerical_instability,
        prop_stable = status$prop_stable,
        mean_count_avg = .s4overall_one_value(perf, c("mean_count_avg", "observed_mean_count_avg", "mean_count")),
        zero_prop_avg = .s4overall_one_value(perf, c("zero_prop_avg", "observed_zero_prop_avg", "zero_prop")),
        max_count_max = .s4overall_one_value(perf, c("max_count_max", "max_count")),
        beta1_mean_avg = .s4overall_one_value(perf, c("beta1_mean_avg", "beta1_mean")),
        beta2_mean_avg = .s4overall_one_value(perf, c("beta2_mean_avg", "beta2_mean")),
        beta1_coverage = .s4overall_one_value(perf, c("beta1_coverage", "beta1_covered_avg")),
        beta2_coverage = .s4overall_one_value(perf, c("beta2_coverage", "beta2_covered_avg")),
        r_truth_avg = .s4overall_one_value(perf, c("r_truth_avg", "r_truth")),
        r_mean_avg = .s4overall_one_value(perf, c("r_mean_avg", "r_mean")),
        lambda_rmse_avg = .s4overall_one_value(perf, c("lambda_rmse_avg", "lambda_rmse")),
        log_lambda_rmse_avg = .s4overall_one_value(perf, c("log_lambda_rmse_avg", "log_lambda_rmse")),
        cor_log_lambda_avg = .s4overall_one_value(perf, c("cor_log_lambda_avg", "cor_log_lambda")),
        phi_rmse_avg = .s4overall_one_value(perf, c("phi_rmse_avg", "phi_rmse")),
        phi_cor_avg = .s4overall_one_value(perf, c("phi_cor_avg", "phi_cor")),
        kappa_truth_cv_avg = .s4overall_one_value(perf, c("kappa_truth_cv_avg", "truth_CV", "truth_cv")),
        kappa_post_mean_cv_avg = .s4overall_one_value(perf, c("kappa_post_mean_cv_avg", "post_mean_CV", "post_mean_cv")),
        log_kappa_rmse_avg = .s4overall_one_value(perf, c("log_kappa_rmse_avg", "log_kappa_rmse")),
        cor_log_kappa_avg = .s4overall_one_value(perf, c("cor_log_kappa_avg", "cor_log_kappa")),
        beta_guard_total = .s4overall_one_value(perf, c("beta_guard_total", "s4a_beta_guard_total", "s4b_beta_guard_total", "s4c_beta_guard_total", "s4e_beta_guard_total")),
        kappa_guard_total = .s4overall_one_value(perf, c("kappa_guard_total", "s4a_kappa_guard_total", "s4b_kappa_guard_total", "s4c_kappa_guard_total", "s4e_kappa_guard_total")),
        lambda_input_guard_total = .s4overall_one_value(perf, c("lambda_input_guard_total", "s4a_lambda_input_guard_total", "s4b_lambda_input_guard_total", "s4c_lambda_input_guard_total", "s4e_lambda_input_guard_total")),
        lambda_output_guard_total = .s4overall_one_value(perf, c("lambda_output_guard_total", "s4a_lambda_output_guard_total", "s4b_lambda_output_guard_total", "s4c_lambda_output_guard_total", "s4e_lambda_output_guard_total")),
        special_metric_1 = special$special_metric_1,
        special_value_1 = special$special_value_1,
        special_metric_2 = special$special_metric_2,
        special_value_2 = special$special_value_2,
        special_metric_3 = special$special_metric_3,
        special_value_3 = special$special_value_3,
        main_interpretation = cfg$interpretation,
        stringsAsFactors = FALSE
    )
}

.s4overall_make_stress_calibration <- function(summary_table) {
    data.frame(
        scenario = summary_table$scenario,
        label = summary_table$label,
        stress_mechanism = summary_table$stress_mechanism,
        mean_count_avg = summary_table$mean_count_avg,
        zero_prop_avg = summary_table$zero_prop_avg,
        special_metric_1 = summary_table$special_metric_1,
        special_value_1 = summary_table$special_value_1,
        special_metric_2 = summary_table$special_metric_2,
        special_value_2 = summary_table$special_value_2,
        special_metric_3 = summary_table$special_metric_3,
        special_value_3 = summary_table$special_value_3,
        stringsAsFactors = FALSE
    )
}

.s4overall_make_latent_summary <- function(summary_table) {
    data.frame(
        scenario = summary_table$scenario,
        label = summary_table$label,
        primary_subset = summary_table$primary_subset,
        lambda_rmse_avg = summary_table$lambda_rmse_avg,
        log_lambda_rmse_avg = summary_table$log_lambda_rmse_avg,
        cor_log_lambda_avg = summary_table$cor_log_lambda_avg,
        phi_rmse_avg = summary_table$phi_rmse_avg,
        phi_cor_avg = summary_table$phi_cor_avg,
        kappa_truth_cv_avg = summary_table$kappa_truth_cv_avg,
        kappa_post_mean_cv_avg = summary_table$kappa_post_mean_cv_avg,
        log_kappa_rmse_avg = summary_table$log_kappa_rmse_avg,
        cor_log_kappa_avg = summary_table$cor_log_kappa_avg,
        stringsAsFactors = FALSE
    )
}

.s4overall_make_stability_summary <- function(summary_table) {
    data.frame(
        scenario = summary_table$scenario,
        label = summary_table$label,
        n_reps_total = summary_table$n_reps_total,
        n_stable = summary_table$n_stable,
        n_soft_warning = summary_table$n_soft_warning,
        n_numerical_instability = summary_table$n_numerical_instability,
        prop_stable = summary_table$prop_stable,
        beta_guard_total = summary_table$beta_guard_total,
        kappa_guard_total = summary_table$kappa_guard_total,
        lambda_input_guard_total = summary_table$lambda_input_guard_total,
        lambda_output_guard_total = summary_table$lambda_output_guard_total,
        stringsAsFactors = FALSE
    )
}

.s4overall_make_replicate_master <- function(scenarios) {
    out <- list()
    for (nm in names(scenarios)) {
        obj <- scenarios[[nm]]
        d <- obj$replicate
        if (is.null(d) || nrow(d) == 0L) next
        d$scenario <- obj$cfg$scenario
        d$label <- obj$cfg$label
        out[[nm]] <- d
    }
    .s4overall_bind_rows(out)
}

.s4overall_make_figures <- function(summary_table, fig_dir, formats = c("png", "pdf"), verbose = TRUE) {
    if (is.null(summary_table) || nrow(summary_table) == 0L) return(invisible(NULL))
    .s4overall_dir_create(fig_dir)

    # Fit status stacked counts.
    stability <- .s4overall_make_stability_summary(summary_table)
    for (ext in formats) {
        .s4overall_plot_status(file.path(fig_dir, paste0("scenario4_fit_status_counts.", ext)), stability)
    }

    # Beta2 recovery.
    for (ext in formats) {
        .s4overall_plot_bar(
            file.path(fig_dir, paste0("scenario4_beta2_recovery.", ext)),
            names = summary_table$scenario,
            values = .s4overall_num(summary_table$beta2_mean_avg),
            main = "Scenario 4 beta2 posterior mean by stress test",
            ylab = "Posterior mean of beta2"
        )
    }

    # log-lambda RMSE.
    for (ext in formats) {
        .s4overall_plot_bar(
            file.path(fig_dir, paste0("scenario4_log_lambda_rmse.", ext)),
            names = summary_table$scenario,
            values = .s4overall_num(summary_table$log_lambda_rmse_avg),
            main = "Scenario 4 log-lambda RMSE by stress test",
            ylab = "Mean log-lambda RMSE"
        )
    }

    # Total numerical guards, log scale display.
    total_guards <- .s4overall_num(summary_table$beta_guard_total) +
        .s4overall_num(summary_table$kappa_guard_total) +
        .s4overall_num(summary_table$lambda_input_guard_total) +
        .s4overall_num(summary_table$lambda_output_guard_total)
    total_guards[is.na(total_guards)] <- 0
    for (ext in formats) {
        .s4overall_plot_bar(
            file.path(fig_dir, paste0("scenario4_guard_counts.", ext)),
            names = summary_table$scenario,
            values = log10(total_guards + 1),
            main = "Scenario 4 numerical guard activity",
            ylab = "log10(total guards + 1)"
        )
    }

    .s4overall_msg("Wrote figures to: ", fig_dir, verbose = verbose)
    invisible(TRUE)
}

run_scenario4_overall_summary_continuous_x2 <- function(
    root = ".",
    out_dir = file.path(root, "analysis_s4_overall_continuous_x2"),
    subset_preference = c("stable_plus_soft_warning", "stable_only", "all_sampled_lambda"),
    scenario_config = NULL,
    allow_missing = FALSE,
    make_figures = TRUE,
    figure_formats = c("png", "pdf"),
    verbose = TRUE
) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    options(s4overall_root = root)
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
    table_dir <- file.path(out_dir, "tables")
    fig_dir <- file.path(out_dir, "figures")
    .s4overall_dir_create(table_dir)
    .s4overall_dir_create(fig_dir)

    if (is.null(scenario_config)) scenario_config <- .s4overall_default_config(root)

    scenario_objs <- list()
    for (nm in names(scenario_config)) {
        obj <- .s4overall_read_scenario(
            cfg = scenario_config[[nm]],
            subset_preference = subset_preference,
            allow_missing = allow_missing,
            verbose = verbose
        )
        if (!is.null(obj)) scenario_objs[[nm]] <- obj
    }

    if (length(scenario_objs) == 0L) {
        stop("No scenario-level tables were read successfully.", call. = FALSE)
    }

    cross_summary <- .s4overall_bind_rows(lapply(scenario_objs, .s4overall_make_summary_row, subset_preference = subset_preference))
    # Preserve S4A--S4E order.
    cross_summary$scenario <- factor(cross_summary$scenario, levels = c("S4A", "S4B", "S4C", "S4D", "S4E"), ordered = TRUE)
    cross_summary <- cross_summary[order(cross_summary$scenario), , drop = FALSE]
    cross_summary$scenario <- as.character(cross_summary$scenario)

    stress_summary <- .s4overall_make_stress_calibration(cross_summary)
    regression_summary <- .s4overall_bind_rows(lapply(scenario_objs, function(x) x$beta_summary))
    latent_summary <- .s4overall_make_latent_summary(cross_summary)
    stability_summary <- .s4overall_make_stability_summary(cross_summary)
    replicate_master <- .s4overall_make_replicate_master(scenario_objs)
    manifest <- .s4overall_bind_rows(lapply(scenario_objs, function(x) x$manifest))

    files <- list(
        cross_scenario_summary = file.path(table_dir, "scenario4_cross_scenario_summary.csv"),
        stress_calibration_summary = file.path(table_dir, "scenario4_stress_calibration_summary.csv"),
        regression_summary = file.path(table_dir, "scenario4_regression_summary.csv"),
        latent_spatial_kappa_summary = file.path(table_dir, "scenario4_latent_spatial_kappa_summary.csv"),
        numerical_stability_summary = file.path(table_dir, "scenario4_numerical_stability_summary.csv"),
        replicate_level_master = file.path(table_dir, "scenario4_replicate_level_master.csv"),
        manifest = file.path(table_dir, "scenario4_overall_manifest.csv")
    )

    utils::write.csv(cross_summary, files$cross_scenario_summary, row.names = FALSE)
    utils::write.csv(stress_summary, files$stress_calibration_summary, row.names = FALSE)
    utils::write.csv(regression_summary, files$regression_summary, row.names = FALSE)
    utils::write.csv(latent_summary, files$latent_spatial_kappa_summary, row.names = FALSE)
    utils::write.csv(stability_summary, files$numerical_stability_summary, row.names = FALSE)
    utils::write.csv(replicate_master, files$replicate_level_master, row.names = FALSE)
    utils::write.csv(manifest, files$manifest, row.names = FALSE)

    if (isTRUE(make_figures)) {
        .s4overall_make_figures(cross_summary, fig_dir, formats = figure_formats, verbose = verbose)
    }

    .s4overall_msg("Wrote Scenario 4 overall tables to: ", table_dir, verbose = verbose)

    out <- list(
        cross_scenario_summary = cross_summary,
        stress_calibration_summary = stress_summary,
        regression_summary = regression_summary,
        latent_spatial_kappa_summary = latent_summary,
        numerical_stability_summary = stability_summary,
        replicate_level_master = replicate_master,
        manifest = manifest,
        files = files,
        table_dir = table_dir,
        fig_dir = fig_dir,
        scenario_objects = scenario_objs
    )
    invisible(out)
}

# Helper to locate candidate S4D CSV files when finalized continuous-time x2
# tables are not under the expected folder names.
find_s4d_summary_files <- function(root = ".") {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    all_csv <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(all_csv) == 0L) return(data.frame())
    p <- tolower(all_csv)
    keep <- grepl("s4d|scenario4d|short.*time|short_time", p)
    out <- data.frame(path = all_csv[keep], stringsAsFactors = FALSE)
    if (nrow(out) == 0L) return(out)
    out$has_continuous <- grepl("continuous|continuous_x2|continuous-time", tolower(out$path))
    out$file <- basename(out$path)
    out$header <- vapply(out$path, function(z) paste(.s4overall_csv_header(z), collapse = ", "), character(1L))
    out[order(!out$has_continuous, out$file, out$path), c("has_continuous", "file", "path", "header")]
}

# Convenience auto-run when sourced interactively.
# Comment this block out if you prefer to call the function manually.
if (sys.nframe() == 0L) {
    s4_overall <- run_scenario4_overall_summary_continuous_x2(root = ".")
    print(s4_overall$cross_scenario_summary)
}

# ================================================================
# v6 overrides: guard totals and S4B special metrics
# ================================================================
# These definitions intentionally override earlier v5 definitions.
# They are placed after all original function definitions so that the
# main run function will use the corrected versions at call time.

.s4overall_choose_all_row <- function(d) {
    if (is.null(d) || nrow(d) == 0L) return(NULL)
    if ("subset" %in% names(d)) {
        ss <- tolower(as.character(d$subset))
        idx <- which(ss %in% c("all_sampled_lambda", "all", "all_reps", "all replicates"))
        if (length(idx) > 0L) return(d[idx[[1L]], , drop = FALSE])
    }
    d[1L, , drop = FALSE]
}

.s4overall_zero_if_missing <- function(d, candidates) {
    nm <- names(d)
    idx <- match(candidates, nm)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(0)
    x <- suppressWarnings(as.numeric(d[[idx[[1L]]]][[1L]]))
    if (length(x) == 0L || is.na(x)) return(0)
    x
}

.s4overall_guard_values <- function(scenario_obj) {
    guard_row <- NULL
    if (!is.null(scenario_obj$guard) && nrow(scenario_obj$guard) > 0L) {
        guard_row <- .s4overall_choose_all_row(scenario_obj$guard)
    }

    # Use all-sampled guard row when available. This avoids double-counting
    # subset rows, which can otherwise happen for S4A if the guard table has
    # all / stable / numerical rows.
    if (!is.null(guard_row) && nrow(guard_row) > 0L) {
        return(list(
            beta_guard_total = .s4overall_zero_if_missing(guard_row, c("beta_guard_total", "beta_guard_count", "s4a_beta_guard_count", "s4b_beta_guard_count", "s4c_beta_guard_count", "s4d_beta_guard_count", "s4e_beta_guard_count", "Beta guards", "Beta.guards")),
            kappa_guard_total = .s4overall_zero_if_missing(guard_row, c("kappa_guard_total", "kappa_guard_count", "s4a_kappa_guard_count", "s4b_kappa_guard_count", "s4c_kappa_guard_count", "s4d_kappa_guard_count", "s4e_kappa_guard_count", "Kappa guards", "Kappa.guards")),
            lambda_input_guard_total = .s4overall_zero_if_missing(guard_row, c("lambda_input_guard_total", "lambda_input_guard_count", "s4a_lambda_input_guard_count", "s4b_lambda_input_guard_count", "s4c_lambda_input_guard_count", "s4d_lambda_input_guard_count", "s4e_lambda_input_guard_count", "Lambda input guards", "Lambda.input.guards")),
            lambda_output_guard_total = .s4overall_zero_if_missing(guard_row, c("lambda_output_guard_total", "lambda_output_guard_count", "s4a_lambda_output_guard_count", "s4b_lambda_output_guard_count", "s4c_lambda_output_guard_count", "s4d_lambda_output_guard_count", "s4e_lambda_output_guard_count", "Lambda output guards", "Lambda.output.guards"))
        ))
    }

    # Fallback: use all-sampled performance row, not the primary stable row.
    perf_all <- .s4overall_choose_all_row(scenario_obj$performance)
    if (!is.null(perf_all) && nrow(perf_all) > 0L) {
        return(list(
            beta_guard_total = .s4overall_zero_if_missing(perf_all, c("beta_guard_total", "s4a_beta_guard_total", "s4b_beta_guard_total", "s4c_beta_guard_total", "s4d_beta_guard_total", "s4e_beta_guard_total")),
            kappa_guard_total = .s4overall_zero_if_missing(perf_all, c("kappa_guard_total", "s4a_kappa_guard_total", "s4b_kappa_guard_total", "s4c_kappa_guard_total", "s4d_kappa_guard_total", "s4e_kappa_guard_total")),
            lambda_input_guard_total = .s4overall_zero_if_missing(perf_all, c("lambda_input_guard_total", "s4a_lambda_input_guard_total", "s4b_lambda_input_guard_total", "s4c_lambda_input_guard_total", "s4d_lambda_input_guard_total", "s4e_lambda_input_guard_total")),
            lambda_output_guard_total = .s4overall_zero_if_missing(perf_all, c("lambda_output_guard_total", "s4a_lambda_output_guard_total", "s4b_lambda_output_guard_total", "s4c_lambda_output_guard_total", "s4d_lambda_output_guard_total", "s4e_lambda_output_guard_total"))
        ))
    }

    list(beta_guard_total = 0, kappa_guard_total = 0, lambda_input_guard_total = 0, lambda_output_guard_total = 0)
}

.s4overall_one_value_grep <- function(d, patterns, default = NA_real_) {
    if (is.null(d) || nrow(d) == 0L) return(default)
    nm <- names(d)
    nml <- tolower(nm)
    for (pat in patterns) {
        idx <- grep(pat, nml, perl = TRUE)
        if (length(idx) > 0L) {
            x <- suppressWarnings(as.numeric(d[[idx[[1L]]]][[1L]]))
            if (length(x) > 0L && !is.na(x)) return(x)
        }
    }
    default
}

.s4overall_special_values <- function(scenario_obj, subset_preference) {
    sc <- scenario_obj$cfg$scenario
    perf <- scenario_obj$performance_row
    special <- scenario_obj$special
    guards <- .s4overall_guard_values(scenario_obj)

    out <- list(
        special_metric_1 = NA_character_, special_value_1 = NA_real_,
        special_metric_2 = NA_character_, special_value_2 = NA_real_,
        special_metric_3 = NA_character_, special_value_3 = NA_real_
    )

    if (sc == "S4A") {
        out$special_metric_1 <- "stable_reps"
        out$special_value_1 <- scenario_obj$status$stable
        out$special_metric_2 <- "numerical_instability_reps"
        out$special_value_2 <- scenario_obj$status$numerical_instability
        out$special_metric_3 <- "lambda_output_guard_total"
        out$special_value_3 <- guards$lambda_output_guard_total
    } else if (sc == "S4B") {
        low_high <- NULL
        if (!is.null(special$low_high)) low_high <- .s4overall_choose_subset_row(special$low_high$data, subset_preference)
        source_row <- if (!is.null(low_high) && nrow(low_high) > 0L) low_high else perf
        out$special_metric_1 <- "low_exp_log_lambda_rmse"
        out$special_value_1 <- .s4overall_one_value(source_row, c(
            "low_log_lambda_rmse", "low_log_lambda_rmse_avg",
            "low_exp_log_lambda_rmse", "low_exp_log_lambda_rmse_avg", "low_exp_log_lambda_rmse_mean",
            "low_exposure_log_lambda_rmse", "low_exposure_log_lambda_rmse_avg",
            "lambda_low_exposure_log_lambda_rmse_avg"
        ))
        if (is.na(out$special_value_1)) out$special_value_1 <- .s4overall_one_value_grep(source_row, c("low.*log.*lambda.*rmse", "low.*log.*rmse"))

        out$special_metric_2 <- "high_exp_log_lambda_rmse"
        out$special_value_2 <- .s4overall_one_value(source_row, c(
            "high_log_lambda_rmse", "high_log_lambda_rmse_avg",
            "high_exp_log_lambda_rmse", "high_exp_log_lambda_rmse_avg", "high_exp_log_lambda_rmse_mean",
            "high_exposure_log_lambda_rmse", "high_exposure_log_lambda_rmse_avg",
            "lambda_high_exposure_log_lambda_rmse_avg"
        ))
        if (is.na(out$special_value_2)) out$special_value_2 <- .s4overall_one_value_grep(source_row, c("high.*log.*lambda.*rmse", "high.*log.*rmse"))

        out$special_metric_3 <- "low_minus_high_log_lambda_rmse"
        out$special_value_3 <- .s4overall_one_value(source_row, c(
            "low_minus_high_log_lambda_rmse", "low_minus_high_log_lambda_rmse_avg",
            "low_high_log_lambda_rmse_gap", "low_exp_minus_high_exp_log_lambda_rmse_avg",
            "low_exp_minus_high_exp_log_lambda_rmse", "low_exposure_minus_high_exposure_log_lambda_rmse_avg"
        ))
        if (is.na(out$special_value_3)) out$special_value_3 <- .s4overall_one_value_grep(source_row, c("low.*minus.*high.*log.*lambda.*rmse", "low.*high.*gap"))
        if (is.na(out$special_value_3) && !is.na(out$special_value_1) && !is.na(out$special_value_2)) {
            out$special_value_3 <- out$special_value_1 - out$special_value_2
        }
    } else if (sc == "S4C") {
        kappa <- NULL
        if (!is.null(special$kappa)) kappa <- .s4overall_choose_subset_row(special$kappa$data, subset_preference)
        source_row <- if (!is.null(kappa) && nrow(kappa) > 0L) kappa else perf
        out$special_metric_1 <- "kappa_truth_cv"
        out$special_value_1 <- .s4overall_one_value(source_row, c("kappa_truth_cv_avg", "truth_CV", "truth_cv"))
        out$special_metric_2 <- "reference_kappa_cv"
        out$special_value_2 <- .s4overall_one_value(source_row, c("reference_kappa_cv_avg", "reference_cv"))
        out$special_metric_3 <- "kappa_cv_increase"
        out$special_value_3 <- .s4overall_one_value(source_row, c("kappa_cv_increase_avg", "cv_increase"))
    } else if (sc == "S4D") {
        out$special_metric_1 <- "time_points"
        out$special_value_1 <- .s4overall_one_value(perf, c("TT", "TT_unique", "time_points", "T"), default = 30)
        out$special_metric_2 <- "short_T_log_lambda_rmse"
        out$special_value_2 <- .s4overall_one_value(perf, c("log_lambda_rmse_avg", "mean_log_lambda_rmse"))
        out$special_metric_3 <- "short_T_cor_log_lambda"
        out$special_value_3 <- .s4overall_one_value(perf, c("cor_log_lambda_avg", "mean_cor_log_lambda"))
    } else if (sc == "S4E") {
        conf <- NULL
        if (!is.null(special$confounding)) conf <- .s4overall_choose_subset_row(special$confounding$data, subset_preference)
        source_row <- if (!is.null(conf) && nrow(conf) > 0L) conf else perf
        out$special_metric_1 <- "abs_cell_cor_x2_phi"
        out$special_value_1 <- .s4overall_one_value(source_row, c("x2_phi_abs_cell_cor_avg", "abs_cell_cor_avg"))
        out$special_metric_2 <- "abs_area_mean_cor_x2_phi"
        out$special_value_2 <- .s4overall_one_value(source_row, c("x2_phi_abs_area_mean_cor_avg", "abs_area_cor_avg"))
        out$special_metric_3 <- "confounded_minus_reference_cell_cor"
        out$special_value_3 <- .s4overall_one_value(source_row, c("confounded_minus_reference_cell_cor", "cell_cor_increase"))
    }
    out
}

.s4overall_make_summary_row <- function(scenario_obj, subset_preference) {
    cfg <- scenario_obj$cfg
    perf <- scenario_obj$performance_row
    status <- scenario_obj$status
    special <- .s4overall_special_values(scenario_obj, subset_preference)
    guards <- .s4overall_guard_values(scenario_obj)

    data.frame(
        scenario = cfg$scenario,
        label = cfg$label,
        stress_mechanism = cfg$stress_mechanism,
        primary_subset = if ("subset" %in% names(perf)) as.character(perf$subset[[1L]]) else "first_row",
        n_reps_total = status$n_reps,
        n_stable = status$stable,
        n_soft_warning = status$soft_warning,
        n_numerical_instability = status$numerical_instability,
        prop_stable = status$prop_stable,
        mean_count_avg = .s4overall_one_value(perf, c("mean_count_avg", "observed_mean_count_avg", "mean_count")),
        zero_prop_avg = .s4overall_one_value(perf, c("zero_prop_avg", "observed_zero_prop_avg", "zero_prop")),
        max_count_max = .s4overall_one_value(perf, c("max_count_max", "max_count")),
        beta1_mean_avg = .s4overall_one_value(perf, c("beta1_mean_avg", "beta1_mean")),
        beta2_mean_avg = .s4overall_one_value(perf, c("beta2_mean_avg", "beta2_mean")),
        beta1_coverage = .s4overall_one_value(perf, c("beta1_coverage", "beta1_covered_avg")),
        beta2_coverage = .s4overall_one_value(perf, c("beta2_coverage", "beta2_covered_avg")),
        r_truth_avg = .s4overall_one_value(perf, c("r_truth_avg", "r_truth")),
        r_mean_avg = .s4overall_one_value(perf, c("r_mean_avg", "r_mean")),
        lambda_rmse_avg = .s4overall_one_value(perf, c("lambda_rmse_avg", "lambda_rmse")),
        log_lambda_rmse_avg = .s4overall_one_value(perf, c("log_lambda_rmse_avg", "log_lambda_rmse")),
        cor_log_lambda_avg = .s4overall_one_value(perf, c("cor_log_lambda_avg", "cor_log_lambda")),
        phi_rmse_avg = .s4overall_one_value(perf, c("phi_rmse_avg", "phi_rmse")),
        phi_cor_avg = .s4overall_one_value(perf, c("phi_cor_avg", "phi_cor")),
        kappa_truth_cv_avg = .s4overall_one_value(perf, c("kappa_truth_cv_avg", "truth_CV", "truth_cv")),
        kappa_post_mean_cv_avg = .s4overall_one_value(perf, c("kappa_post_mean_cv_avg", "post_mean_CV", "post_mean_cv")),
        log_kappa_rmse_avg = .s4overall_one_value(perf, c("log_kappa_rmse_avg", "log_kappa_rmse")),
        cor_log_kappa_avg = .s4overall_one_value(perf, c("cor_log_kappa_avg", "cor_log_kappa")),
        beta_guard_total = guards$beta_guard_total,
        kappa_guard_total = guards$kappa_guard_total,
        lambda_input_guard_total = guards$lambda_input_guard_total,
        lambda_output_guard_total = guards$lambda_output_guard_total,
        special_metric_1 = special$special_metric_1,
        special_value_1 = special$special_value_1,
        special_metric_2 = special$special_metric_2,
        special_value_2 = special$special_value_2,
        special_metric_3 = special$special_metric_3,
        special_value_3 = special$special_value_3,
        main_interpretation = cfg$interpretation,
        stringsAsFactors = FALSE
    )
}
