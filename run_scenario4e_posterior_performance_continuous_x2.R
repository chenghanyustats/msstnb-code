## ============================================================================
## VERSION: S4E_CONTINUOUS_X2_POSTERIOR_PERFORMANCE_V1_2026_06_22
## run_scenario4e_posterior_performance_continuous_x2.R
##
## Posterior-performance summaries for Scenario 4E:
## strong spatial/covariate confounding stress test with fixed gamma and
## continuous-time x2.
##
## Typical use from project root:
##   source("run_scenario4e_posterior_performance_continuous_x2.R")
##
## Or programmatically:
##   s4e_perf <- run_scenario4e_posterior_performance(root = ".")
## ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.s4e_norm_path <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

.s4e_read_required_csv <- function(path) {
    if (!file.exists(path)) stop("Required CSV not found: ", path, call. = FALSE)
    utils::read.csv(path, stringsAsFactors = FALSE)
}

.s4e_col <- function(d, nm, default = NA_real_) {
    if (!nm %in% names(d)) return(rep(default, nrow(d)))
    d[[nm]]
}

.s4e_num1 <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    suppressWarnings(as.numeric(x[[1L]]))
}

.s4e_safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.s4e_safe_sd <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (sum(is.finite(x)) <= 1L) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

.s4e_safe_sum <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(0)
    sum(x, na.rm = TRUE)
}

.s4e_safe_median <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    stats::median(x, na.rm = TRUE)
}

.s4e_safe_cor <- function(x, y) {
    x <- suppressWarnings(as.numeric(x))
    y <- suppressWarnings(as.numeric(y))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
    stats::cor(x[ok], y[ok])
}

.s4e_reps_label <- function(d) {
    if (is.null(d) || nrow(d) == 0L || !"rep_id" %in% names(d)) return(NA_character_)
    paste(as.integer(d$rep_id), collapse = ",")
}

.s4e_first_available <- function(d, nms, default = NA_real_) {
    for (nm in nms) if (nm %in% names(d)) return(d[[nm]])
    rep(default, nrow(d))
}

.s4e_fix_aliases <- function(d) {
    ## Generic aliases inherited from the S3 summary helper.
    if (!"mean_count" %in% names(d) && "observed_mean_count" %in% names(d)) d$mean_count <- d$observed_mean_count
    if (!"zero_prop" %in% names(d) && "observed_zero_prop" %in% names(d)) d$zero_prop <- d$observed_zero_prop
    if (!"total_count" %in% names(d) && "observed_total_count" %in% names(d)) d$total_count <- d$observed_total_count
    if (!"max_count" %in% names(d) && "observed_max_count" %in% names(d)) d$max_count <- d$observed_max_count

    ## Truth aliases.
    if (!"r_truth_for_eval" %in% names(d)) {
        d$r_truth_for_eval <- .s4e_first_available(d, c("r_truth", "r_reference_truth", "r_stress_truth"), 15)
    }
    if (!"beta0_truth_for_eval" %in% names(d)) {
        d$beta0_truth_for_eval <- .s4e_first_available(d, c("beta0_true", "beta0_confounding_truth", "beta0_reference_truth"), NA_real_)
    }
    if (!"beta1_truth_for_eval" %in% names(d)) {
        d$beta1_truth_for_eval <- .s4e_first_available(d, c("beta1_true"), 0.5)
    }
    if (!"beta2_truth_for_eval" %in% names(d)) {
        d$beta2_truth_for_eval <- .s4e_first_available(d, c("beta2_true"), -0.4)
    }

    ## Guard columns may be absent in older interim runs.
    for (nm in c(
        "s4e_beta_guard_count", "s4e_kappa_guard_count",
        "s4e_lambda_input_guard_count", "s4e_lambda_output_guard_count"
    )) {
        if (!nm %in% names(d)) d[[nm]] <- 0
    }
    d
}

.s4e_classify_one <- function(row) {
    beta_guard <- .s4e_num1(row$s4e_beta_guard_count, 0)
    kappa_guard <- .s4e_num1(row$s4e_kappa_guard_count, 0)
    lambda_input_guard <- .s4e_num1(row$s4e_lambda_input_guard_count, 0)
    lambda_output_guard <- .s4e_num1(row$s4e_lambda_output_guard_count, 0)

    beta0_mean <- .s4e_num1(row$beta0_mean, NA_real_)
    beta1_mean <- .s4e_num1(row$beta1_mean, NA_real_)
    beta2_mean <- .s4e_num1(row$beta2_mean, NA_real_)
    r_mean <- .s4e_num1(row$r_mean, NA_real_)
    lambda_rmse <- .s4e_num1(row$lambda_rmse, NA_real_)
    log_lambda_rmse <- .s4e_num1(row$log_lambda_rmse, NA_real_)

    severe <-
        (!is.finite(beta0_mean)) || (!is.finite(beta1_mean)) || (!is.finite(beta2_mean)) ||
        (!is.finite(r_mean)) || (!is.finite(lambda_rmse)) || (!is.finite(log_lambda_rmse)) ||
        abs(beta0_mean) > 30 || abs(beta1_mean) > 5 || abs(beta2_mean) > 5 ||
        r_mean <= 0 ||
        lambda_rmse > 10 || log_lambda_rmse > 5 ||
        beta_guard > 100 || kappa_guard > 100 ||
        lambda_input_guard > 100 || lambda_output_guard > 100

    if (isTRUE(severe)) return("numerical_instability")

    warn <- beta_guard > 0 || kappa_guard > 0 || lambda_input_guard > 0 || lambda_output_guard > 0
    if (isTRUE(warn)) return("soft_warning")

    "stable"
}

.s4e_add_fit_status <- function(d) {
    d <- .s4e_fix_aliases(d)
    d$fit_status <- vapply(seq_len(nrow(d)), function(i) .s4e_classify_one(d[i, , drop = FALSE]), character(1))
    d$nonfailed <- d$fit_status %in% c("stable", "soft_warning")
    d$stable_fit <- d$fit_status == "stable"
    d
}

.s4e_beta_fields <- function(d, p) {
    truth_nm <- switch(p,
        beta0 = "beta0_truth_for_eval",
        beta1 = "beta1_truth_for_eval",
        beta2 = "beta2_truth_for_eval"
    )
    truth <- .s4e_col(d, truth_nm)
    meanv <- .s4e_col(d, paste0(p, "_mean"))
    sdv <- .s4e_col(d, paste0(p, "_sd"))
    q025 <- .s4e_col(d, paste0(p, "_q025"))
    q50 <- .s4e_col(d, paste0(p, "_q50"))
    q975 <- .s4e_col(d, paste0(p, "_q975"))
    if (all(is.na(q50)) && paste0(p, "_median") %in% names(d)) q50 <- d[[paste0(p, "_median")]]
    bias <- .s4e_col(d, paste0(p, "_bias"), meanv - truth)
    covered <- .s4e_col(d, paste0(p, "_covered"), as.integer(q025 <= truth & truth <= q975))
    list(truth = truth, meanv = meanv, sdv = sdv, q025 = q025, q50 = q50, q975 = q975, bias = bias, covered = covered)
}

.s4e_make_beta_long <- function(d, scenario_id) {
    params <- c("beta0", "beta1", "beta2")
    out <- lapply(params, function(p) {
        f <- .s4e_beta_fields(d, p)
        data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(d$rep_id),
            parameter = p,
            truth = as.numeric(f$truth),
            post_mean = as.numeric(f$meanv),
            post_sd = as.numeric(f$sdv),
            q025 = as.numeric(f$q025),
            q50 = as.numeric(f$q50),
            q975 = as.numeric(f$q975),
            bias = as.numeric(f$bias),
            covered = as.integer(f$covered),
            fit_status = d$fit_status,
            nonfailed = d$nonfailed,
            stable_fit = d$stable_fit,
            stringsAsFactors = FALSE
        )
    })
    ans <- do.call(rbind, out)
    rownames(ans) <- NULL
    ans
}

.s4e_subset_df <- function(d, subset_name) {
    if (subset_name == "all_sampled_lambda") return(d)
    if (subset_name == "stable_only") return(d[d$fit_status == "stable", , drop = FALSE])
    if (subset_name == "stable_plus_soft_warning") return(d[d$fit_status %in% c("stable", "soft_warning"), , drop = FALSE])
    if (subset_name == "numerical_instability") return(d[d$fit_status == "numerical_instability", , drop = FALSE])
    d[FALSE, , drop = FALSE]
}

.s4e_beta_aggregate <- function(beta_long, subset_name) {
    if (is.null(beta_long) || nrow(beta_long) == 0L) return(data.frame())
    ans <- lapply(split(beta_long, beta_long$parameter), function(x) {
        data.frame(
            subset = subset_name,
            parameter = unique(x$parameter),
            n_reps = length(unique(x$rep_id)),
            mean_truth = .s4e_safe_mean(x$truth),
            mean_post_mean = .s4e_safe_mean(x$post_mean),
            sd_post_mean_across_reps = .s4e_safe_sd(x$post_mean),
            mean_bias = .s4e_safe_mean(x$bias),
            median_bias = .s4e_safe_median(x$bias),
            rmse_bias = sqrt(.s4e_safe_mean(x$bias^2)),
            mean_abs_bias = .s4e_safe_mean(abs(x$bias)),
            coverage = .s4e_safe_mean(x$covered),
            mean_post_sd = .s4e_safe_mean(x$post_sd),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, ans)
    rownames(out) <- NULL
    out
}

.s4e_performance_aggregate <- function(d, subset_name) {
    if (nrow(d) == 0L) {
        return(data.frame(subset = subset_name, n_reps = 0L, reps = NA_character_, stringsAsFactors = FALSE))
    }
    data.frame(
        subset = subset_name,
        n_reps = nrow(d),
        reps = .s4e_reps_label(d),
        mean_count_avg = .s4e_safe_mean(.s4e_col(d, "mean_count")),
        zero_prop_avg = .s4e_safe_mean(.s4e_col(d, "zero_prop")),
        reference_mean_count_avg = .s4e_safe_mean(.s4e_col(d, "reference_mean_count")),
        total_count_sum = .s4e_safe_sum(.s4e_col(d, "total_count")),
        max_count_max = suppressWarnings(max(as.numeric(.s4e_col(d, "max_count")), na.rm = TRUE)),
        x2_confounding_strength_avg = .s4e_safe_mean(.s4e_col(d, "x2_confounding_strength")),
        x2_phi_abs_cell_cor_avg = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_cell_cor")),
        x2_phi_abs_area_mean_cor_avg = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_area_mean_cor")),
        reference_x2_phi_abs_cell_cor_avg = .s4e_safe_mean(.s4e_col(d, "reference_x2_phi_abs_cell_cor")),
        reference_x2_phi_abs_area_mean_cor_avg = .s4e_safe_mean(.s4e_col(d, "reference_x2_phi_abs_area_mean_cor")),
        beta0_mean_avg = .s4e_safe_mean(.s4e_col(d, "beta0_mean")),
        beta1_mean_avg = .s4e_safe_mean(.s4e_col(d, "beta1_mean")),
        beta1_mean_sd = .s4e_safe_sd(.s4e_col(d, "beta1_mean")),
        beta2_mean_avg = .s4e_safe_mean(.s4e_col(d, "beta2_mean")),
        beta2_mean_sd = .s4e_safe_sd(.s4e_col(d, "beta2_mean")),
        beta0_coverage = .s4e_safe_mean(.s4e_col(d, "beta0_covered")),
        beta1_coverage = .s4e_safe_mean(.s4e_col(d, "beta1_covered")),
        beta2_coverage = .s4e_safe_mean(.s4e_col(d, "beta2_covered")),
        r_truth_avg = .s4e_safe_mean(.s4e_col(d, "r_truth_for_eval")),
        r_mean_avg = .s4e_safe_mean(.s4e_col(d, "r_mean")),
        r_mean_sd = .s4e_safe_sd(.s4e_col(d, "r_mean")),
        r_bias_avg = .s4e_safe_mean(.s4e_col(d, "r_mean") - .s4e_col(d, "r_truth_for_eval")),
        r_abs_bias_avg = .s4e_safe_mean(abs(.s4e_col(d, "r_mean") - .s4e_col(d, "r_truth_for_eval"))),
        r_region_coverage_95_avg = .s4e_safe_mean(.s4e_col(d, "r_region_coverage_95")),
        phi_rmse_avg = .s4e_safe_mean(.s4e_col(d, "phi_rmse")),
        phi_cor_avg = .s4e_safe_mean(.s4e_col(d, "phi_cor")),
        lambda_rmse_avg = .s4e_safe_mean(.s4e_col(d, "lambda_rmse")),
        lambda_rmse_median = .s4e_safe_median(.s4e_col(d, "lambda_rmse")),
        log_lambda_rmse_avg = .s4e_safe_mean(.s4e_col(d, "log_lambda_rmse")),
        log_lambda_rmse_median = .s4e_safe_median(.s4e_col(d, "log_lambda_rmse")),
        cor_log_lambda_avg = .s4e_safe_mean(.s4e_col(d, "cor_log_lambda")),
        kappa_truth_cv_avg = .s4e_safe_mean(.s4e_col(d, "kappa_truth_cv")),
        kappa_post_mean_cv_avg = .s4e_safe_mean(.s4e_col(d, "kappa_post_mean_cv")),
        kappa_cv_bias_avg = .s4e_safe_mean(.s4e_col(d, "kappa_post_mean_cv") - .s4e_col(d, "kappa_truth_cv")),
        kappa_rmse_avg = .s4e_safe_mean(.s4e_col(d, "kappa_rmse")),
        log_kappa_rmse_avg = .s4e_safe_mean(.s4e_col(d, "log_kappa_rmse")),
        cor_log_kappa_avg = .s4e_safe_mean(.s4e_col(d, "cor_log_kappa")),
        beta_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_beta_guard_count", 0)),
        kappa_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_kappa_guard_count", 0)),
        lambda_input_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_lambda_input_guard_count", 0)),
        lambda_output_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_lambda_output_guard_count", 0)),
        stringsAsFactors = FALSE
    )
}

.s4e_confounding_aggregate <- function(d, subset_name) {
    if (nrow(d) == 0L) return(data.frame(subset = subset_name, n_reps = 0L, stringsAsFactors = FALSE))
    data.frame(
        subset = subset_name,
        n_reps = nrow(d),
        confounding_strength_avg = .s4e_safe_mean(.s4e_col(d, "x2_confounding_strength")),
        x2_phi_abs_cell_cor_avg = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_cell_cor")),
        x2_phi_abs_cell_cor_min = suppressWarnings(min(as.numeric(.s4e_col(d, "x2_phi_abs_cell_cor")), na.rm = TRUE)),
        x2_phi_abs_cell_cor_max = suppressWarnings(max(as.numeric(.s4e_col(d, "x2_phi_abs_cell_cor")), na.rm = TRUE)),
        x2_phi_abs_area_mean_cor_avg = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_area_mean_cor")),
        x2_phi_abs_area_mean_cor_min = suppressWarnings(min(as.numeric(.s4e_col(d, "x2_phi_abs_area_mean_cor")), na.rm = TRUE)),
        x2_phi_abs_area_mean_cor_max = suppressWarnings(max(as.numeric(.s4e_col(d, "x2_phi_abs_area_mean_cor")), na.rm = TRUE)),
        reference_x2_phi_abs_cell_cor_avg = .s4e_safe_mean(.s4e_col(d, "reference_x2_phi_abs_cell_cor")),
        reference_x2_phi_abs_area_mean_cor_avg = .s4e_safe_mean(.s4e_col(d, "reference_x2_phi_abs_area_mean_cor")),
        confounded_minus_reference_cell_cor = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_cell_cor") - .s4e_col(d, "reference_x2_phi_abs_cell_cor")),
        confounded_minus_reference_area_cor = .s4e_safe_mean(.s4e_col(d, "x2_phi_abs_area_mean_cor") - .s4e_col(d, "reference_x2_phi_abs_area_mean_cor")),
        x2_sd_avg = .s4e_safe_mean(.s4e_col(d, "x2_sd")),
        x2_empirical_ar1_mean_avg = .s4e_safe_mean(.s4e_col(d, "x2_empirical_ar1_mean")),
        x2_binary_like_prop_max = suppressWarnings(max(as.numeric(.s4e_col(d, "x2_binary_like_prop")), na.rm = TRUE)),
        beta0_shift_avg = .s4e_safe_mean(.s4e_col(d, "beta0_shift_to_preserve_reference_mean")),
        stringsAsFactors = FALSE
    )
}

.s4e_r_aggregate <- function(d, subset_name) {
    if (nrow(d) == 0L) return(data.frame(subset = subset_name, n_reps = 0L, stringsAsFactors = FALSE))
    data.frame(
        subset = subset_name,
        n_reps = nrow(d),
        r_truth_avg = .s4e_safe_mean(.s4e_col(d, "r_truth_for_eval")),
        r_mean_avg = .s4e_safe_mean(.s4e_col(d, "r_mean")),
        r_mean_sd = .s4e_safe_sd(.s4e_col(d, "r_mean")),
        r_bias_avg = .s4e_safe_mean(.s4e_col(d, "r_mean") - .s4e_col(d, "r_truth_for_eval")),
        r_abs_bias_avg = .s4e_safe_mean(abs(.s4e_col(d, "r_mean") - .s4e_col(d, "r_truth_for_eval"))),
        r_region_coverage_95_avg = .s4e_safe_mean(.s4e_col(d, "r_region_coverage_95")),
        r_region_rmse_avg = .s4e_safe_mean(.s4e_col(d, "r_region_rmse")),
        r_region_mae_avg = .s4e_safe_mean(.s4e_col(d, "r_region_mae")),
        stringsAsFactors = FALSE
    )
}

.s4e_kappa_aggregate <- function(d, subset_name) {
    if (nrow(d) == 0L) return(data.frame(subset = subset_name, n_reps = 0L, stringsAsFactors = FALSE))
    data.frame(
        subset = subset_name,
        n_reps = nrow(d),
        kappa_truth_cv_avg = .s4e_safe_mean(.s4e_col(d, "kappa_truth_cv")),
        kappa_post_mean_cv_avg = .s4e_safe_mean(.s4e_col(d, "kappa_post_mean_cv")),
        kappa_cv_bias_avg = .s4e_safe_mean(.s4e_col(d, "kappa_post_mean_cv") - .s4e_col(d, "kappa_truth_cv")),
        kappa_rmse_avg = .s4e_safe_mean(.s4e_col(d, "kappa_rmse")),
        log_kappa_rmse_avg = .s4e_safe_mean(.s4e_col(d, "log_kappa_rmse")),
        cor_log_kappa_avg = .s4e_safe_mean(.s4e_col(d, "cor_log_kappa")),
        stringsAsFactors = FALSE
    )
}

.s4e_guard_aggregate <- function(d, subset_name) {
    if (nrow(d) == 0L) return(data.frame(subset = subset_name, n_reps = 0L, stringsAsFactors = FALSE))
    data.frame(
        subset = subset_name,
        n_reps = nrow(d),
        beta_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_beta_guard_count", 0)),
        kappa_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_kappa_guard_count", 0)),
        lambda_input_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_lambda_input_guard_count", 0)),
        lambda_output_guard_total = .s4e_safe_sum(.s4e_col(d, "s4e_lambda_output_guard_count", 0)),
        stringsAsFactors = FALSE
    )
}

.s4e_beta_confounding_diagnostics <- function(beta_long, rep_summary, subset_name) {
    if (nrow(rep_summary) == 0L || nrow(beta_long) == 0L) return(data.frame(subset = subset_name, stringsAsFactors = FALSE))
    x <- rep_summary[, c("rep_id", "x2_phi_abs_cell_cor", "x2_phi_abs_area_mean_cor", "reference_x2_phi_abs_cell_cor", "reference_x2_phi_abs_area_mean_cor"), drop = FALSE]
    out <- lapply(c("beta0", "beta1", "beta2"), function(p) {
        b <- beta_long[beta_long$parameter == p & beta_long$rep_id %in% x$rep_id, , drop = FALSE]
        m <- merge(b[, c("rep_id", "bias", "post_mean")], x, by = "rep_id", all.x = TRUE, sort = FALSE)
        data.frame(
            subset = subset_name,
            parameter = p,
            n_reps = nrow(m),
            mean_bias = .s4e_safe_mean(m$bias),
            mean_abs_bias = .s4e_safe_mean(abs(m$bias)),
            cor_bias_abs_cell_confounding = .s4e_safe_cor(m$bias, m$x2_phi_abs_cell_cor),
            cor_abs_bias_abs_cell_confounding = .s4e_safe_cor(abs(m$bias), m$x2_phi_abs_cell_cor),
            cor_bias_abs_area_confounding = .s4e_safe_cor(m$bias, m$x2_phi_abs_area_mean_cor),
            cor_abs_bias_abs_area_confounding = .s4e_safe_cor(abs(m$bias), m$x2_phi_abs_area_mean_cor),
            stringsAsFactors = FALSE
        )
    })
    ans <- do.call(rbind, out)
    rownames(ans) <- NULL
    ans
}

.s4e_optional_s3_row <- function() {
    ## Reference values from the previously checked S3 fixed-gamma stabilized summary.
    data.frame(
        scenario = "S3 control",
        subset = "control reps",
        n_reps = 3L,
        mean_count_avg = 195.05,
        zero_prop_avg = 0.000,
        x2_phi_abs_cell_cor_avg = NA_real_,
        x2_phi_abs_area_mean_cor_avg = NA_real_,
        beta1_mean_avg = 0.4918,
        beta2_mean_avg = -0.4077,
        beta2_mean_sd = 0.0089,
        r_mean_avg = NA_real_,
        log_lambda_rmse_avg = 0.0472,
        cor_log_lambda_avg = 0.7340,
        stringsAsFactors = FALSE
    )
}

.s4e_make_plots <- function(rep_summary, beta_long, fig_dir) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
    ggplot2 <- asNamespace("ggplot2")

    status_lab <- function(x) {
        out <- as.character(x)
        out[out == "stable"] <- "Stable"
        out[out == "soft_warning"] <- "Soft warning"
        out[out == "numerical_instability"] <- "Numerical"
        out
    }
    rep_summary$status_label <- status_lab(rep_summary$fit_status)
    beta_long$status_label <- status_lab(beta_long$fit_status)

    ## Fit status counts.
    sc <- as.data.frame(table(rep_summary$status_label), stringsAsFactors = FALSE)
    names(sc) <- c("status_label", "n_reps")
    p0 <- ggplot2$ggplot(sc, ggplot2$aes(x = status_label, y = n_reps)) +
        ggplot2$geom_col(width = 0.65) +
        ggplot2$geom_text(ggplot2$aes(label = n_reps), vjust = -0.35, size = 3.5) +
        ggplot2$labs(x = "Fit status", y = "Number of replicates") +
        ggplot2$theme_bw(base_size = 11)
    ggplot2$ggsave(file.path(fig_dir, "s4e_fit_status_counts.pdf"), p0, width = 5.8, height = 3.4)

    ## Confounding calibration by replicate.
    cc <- rbind(
        data.frame(rep_id = rep_summary$rep_id, value = rep_summary$x2_phi_abs_cell_cor, quantity = "Cell-level abs cor"),
        data.frame(rep_id = rep_summary$rep_id, value = rep_summary$x2_phi_abs_area_mean_cor, quantity = "Area-mean abs cor"),
        data.frame(rep_id = rep_summary$rep_id, value = rep_summary$reference_x2_phi_abs_cell_cor, quantity = "Reference cell-level abs cor"),
        data.frame(rep_id = rep_summary$rep_id, value = rep_summary$reference_x2_phi_abs_area_mean_cor, quantity = "Reference area-mean abs cor")
    )
    p1 <- ggplot2$ggplot(cc, ggplot2$aes(x = rep_id, y = value, linetype = quantity)) +
        ggplot2$geom_point(size = 2.0) +
        ggplot2$geom_line(alpha = 0.55) +
        ggplot2$scale_x_continuous(breaks = sort(unique(rep_summary$rep_id))) +
        ggplot2$labs(x = "Replicate", y = "Absolute correlation", linetype = "Quantity") +
        ggplot2$theme_bw(base_size = 11)
    ggplot2$ggsave(file.path(fig_dir, "s4e_confounding_correlations_by_rep.pdf"), p1, width = 7.4, height = 4.0)

    ## beta1/beta2 posterior means.
    b12 <- beta_long[beta_long$parameter %in% c("beta1", "beta2"), , drop = FALSE]
    p2 <- ggplot2$ggplot(b12, ggplot2$aes(x = rep_id, y = post_mean, group = status_label)) +
        ggplot2$geom_point(ggplot2$aes(shape = status_label), size = 2.2) +
        ggplot2$geom_line(alpha = 0.45) +
        ggplot2$geom_hline(data = unique(b12[, c("parameter", "truth")]), ggplot2$aes(yintercept = truth), linetype = "dashed") +
        ggplot2$facet_wrap(~ parameter, scales = "free_y") +
        ggplot2$scale_x_continuous(breaks = sort(unique(rep_summary$rep_id))) +
        ggplot2$labs(x = "Replicate", y = "Posterior mean", shape = "Status") +
        ggplot2$theme_bw(base_size = 11)
    ggplot2$ggsave(file.path(fig_dir, "s4e_beta_means_by_rep.pdf"), p2, width = 7.2, height = 3.6)

    ## beta2 bias vs confounding strength/correlation.
    b2 <- beta_long[beta_long$parameter == "beta2", c("rep_id", "bias"), drop = FALSE]
    b2 <- merge(b2, rep_summary[, c("rep_id", "x2_phi_abs_cell_cor", "x2_phi_abs_area_mean_cor", "fit_status"), drop = FALSE], by = "rep_id", all.x = TRUE, sort = FALSE)
    b2$status_label <- status_lab(b2$fit_status)
    p3 <- ggplot2$ggplot(b2, ggplot2$aes(x = x2_phi_abs_cell_cor, y = bias)) +
        ggplot2$geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2$geom_point(ggplot2$aes(shape = status_label), size = 2.4) +
        ggplot2$labs(x = "Abs cor(x2, phi), cell level", y = "Beta2 posterior-mean bias", shape = "Status") +
        ggplot2$theme_bw(base_size = 11)
    ggplot2$ggsave(file.path(fig_dir, "s4e_beta2_bias_vs_confounding.pdf"), p3, width = 6.4, height = 3.8)

    ## Global log-lambda RMSE.
    p4 <- ggplot2$ggplot(rep_summary, ggplot2$aes(x = rep_id, y = log_lambda_rmse, group = status_label)) +
        ggplot2$geom_point(ggplot2$aes(shape = status_label), size = 2.2) +
        ggplot2$geom_line(alpha = 0.45) +
        ggplot2$scale_x_continuous(breaks = sort(unique(rep_summary$rep_id))) +
        ggplot2$labs(x = "Replicate", y = "Log-lambda RMSE", shape = "Status") +
        ggplot2$theme_bw(base_size = 11)
    ggplot2$ggsave(file.path(fig_dir, "s4e_log_lambda_rmse_by_rep.pdf"), p4, width = 6.5, height = 3.5)

    invisible(TRUE)
}

run_scenario4e_posterior_performance <- function(root = ".",
                                                 make_plots = TRUE,
                                                 verbose = TRUE) {
    scenario_id <- "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_FIXED_GAMMA_CONTINUOUS_X2_T100"
    sampled_path <- file.path(
        root,
        "output_s4e_spatial_covariate_confounding_fixed_gamma_continuous_x2",
        scenario_id,
        "summary_S4E_spatial_covariate_confounding_fixed_gamma_continuous_x2_all_reps.csv"
    )

    out_base <- file.path(
        root,
        "analysis_s4e_spatial_covariate_confounding_continuous_x2",
        scenario_id
    )
    table_dir <- file.path(out_base, "tables")
    fig_dir <- file.path(out_base, "figures")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    sampled <- .s4e_read_required_csv(sampled_path)
    sampled <- .s4e_add_fit_status(sampled)

    status_counts <- as.data.frame(table(sampled$fit_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("fit_status", "n_reps")
    status_counts$prop_reps <- status_counts$n_reps / nrow(sampled)
    status_counts <- status_counts[order(match(status_counts$fit_status, c("numerical_instability", "soft_warning", "stable"))), , drop = FALSE]
    rownames(status_counts) <- NULL

    beta_long <- .s4e_make_beta_long(sampled, scenario_id)

    subset_names <- c("all_sampled_lambda", "stable_only", "stable_plus_soft_warning", "numerical_instability")
    subset_dfs <- setNames(lapply(subset_names, function(s) .s4e_subset_df(sampled, s)), subset_names)
    subset_names_present <- subset_names[vapply(subset_dfs, nrow, integer(1)) > 0]

    performance_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) .s4e_performance_aggregate(subset_dfs[[s]], s)))
    confounding_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) .s4e_confounding_aggregate(subset_dfs[[s]], s)))
    r_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) .s4e_r_aggregate(subset_dfs[[s]], s)))
    kappa_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) .s4e_kappa_aggregate(subset_dfs[[s]], s)))
    guard_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) .s4e_guard_aggregate(subset_dfs[[s]], s)))
    beta_aggregate <- do.call(rbind, lapply(subset_names_present, function(s) {
        reps <- unique(subset_dfs[[s]]$rep_id)
        .s4e_beta_aggregate(beta_long[beta_long$rep_id %in% reps, , drop = FALSE], s)
    }))
    beta_confounding_diagnostics <- do.call(rbind, lapply(subset_names_present, function(s) {
        reps <- unique(subset_dfs[[s]]$rep_id)
        .s4e_beta_confounding_diagnostics(
            beta_long[beta_long$rep_id %in% reps, , drop = FALSE],
            subset_dfs[[s]],
            s
        )
    }))

    replicate_summary <- sampled[, intersect(c(
        "scenario_id", "data_scenario_id", "rep_id", "fit_status",
        "mean_count", "zero_prop", "reference_mean_count", "observed_total_count", "max_count",
        "x2_confounding_strength", "x2_residual_weight", "preserve_reference_mean",
        "x2_phi_cell_cor", "x2_phi_abs_cell_cor", "x2_phi_area_mean_cor", "x2_phi_abs_area_mean_cor",
        "reference_x2_phi_abs_cell_cor", "reference_x2_phi_abs_area_mean_cor",
        "beta0_reference_truth", "beta0_confounding_truth", "beta0_shift_to_preserve_reference_mean",
        "beta0_mean", "beta1_mean", "beta2_mean", "r_mean", "r_truth_for_eval",
        "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "kappa_truth_cv", "kappa_post_mean_cv", "kappa_rmse", "log_kappa_rmse", "cor_log_kappa",
        "phi_rmse", "phi_cor",
        "s4e_beta_guard_count", "s4e_kappa_guard_count", "s4e_lambda_input_guard_count", "s4e_lambda_output_guard_count",
        "x2_mode", "x2_ar", "x2_innov_sd", "x2_sd", "x2_empirical_ar1_mean", "x2_binary_like_prop"
    ), names(sampled)), drop = FALSE]

    stable_plus <- performance_aggregate[performance_aggregate$subset == "stable_plus_soft_warning", , drop = FALSE]
    scenario_comparison <- rbind(
        .s4e_optional_s3_row(),
        data.frame(
            scenario = "S4E spatial/covariate confounding",
            subset = "stable+soft",
            n_reps = stable_plus$n_reps,
            mean_count_avg = stable_plus$mean_count_avg,
            zero_prop_avg = stable_plus$zero_prop_avg,
            x2_phi_abs_cell_cor_avg = stable_plus$x2_phi_abs_cell_cor_avg,
            x2_phi_abs_area_mean_cor_avg = stable_plus$x2_phi_abs_area_mean_cor_avg,
            beta1_mean_avg = stable_plus$beta1_mean_avg,
            beta2_mean_avg = stable_plus$beta2_mean_avg,
            beta2_mean_sd = stable_plus$beta2_mean_sd,
            r_mean_avg = stable_plus$r_mean_avg,
            log_lambda_rmse_avg = stable_plus$log_lambda_rmse_avg,
            cor_log_lambda_avg = stable_plus$cor_log_lambda_avg,
            stringsAsFactors = FALSE
        )
    )

    manifest <- data.frame(
        item = c("sampled_summary", "table_dir", "figure_dir", "scenario_id"),
        value = c(.s4e_norm_path(sampled_path), .s4e_norm_path(table_dir), .s4e_norm_path(fig_dir), scenario_id),
        stringsAsFactors = FALSE
    )

    ## Write tables.
    utils::write.csv(status_counts, file.path(table_dir, "s4e_fit_status_counts_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(replicate_summary, file.path(table_dir, "s4e_replicate_level_summary_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(beta_long, file.path(table_dir, "posterior_beta_summary_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(beta_aggregate, file.path(table_dir, "posterior_beta_aggregate_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(beta_confounding_diagnostics, file.path(table_dir, "beta_confounding_diagnostics_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(performance_aggregate, file.path(table_dir, "posterior_performance_by_subset_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(confounding_aggregate, file.path(table_dir, "confounding_calibration_by_subset_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(r_aggregate, file.path(table_dir, "r_recovery_by_subset_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(kappa_aggregate, file.path(table_dir, "kappa_recovery_by_subset_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(guard_aggregate, file.path(table_dir, "numerical_guard_summary_by_subset_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(scenario_comparison, file.path(table_dir, "scenario_comparison_s3_s4e_continuous_x2.csv"), row.names = FALSE)
    utils::write.csv(manifest, file.path(table_dir, "scenario4e_posterior_performance_manifest_continuous_x2.csv"), row.names = FALSE)

    if (isTRUE(make_plots)) .s4e_make_plots(sampled, beta_long, fig_dir)

    out <- list(
        status_counts = status_counts,
        replicate_summary = replicate_summary,
        beta_summary = beta_long,
        beta_aggregate = beta_aggregate,
        beta_confounding_diagnostics = beta_confounding_diagnostics,
        performance_aggregate = performance_aggregate,
        confounding_aggregate = confounding_aggregate,
        r_aggregate = r_aggregate,
        kappa_aggregate = kappa_aggregate,
        guard_aggregate = guard_aggregate,
        scenario_comparison = scenario_comparison,
        manifest = manifest,
        table_dir = table_dir,
        fig_dir = fig_dir,
        sampled_summary_path = sampled_path
    )

    saveRDS(out, file.path(out_base, "scenario4e_posterior_performance_results.rds"))

    if (isTRUE(verbose)) {
        cat("\n================ Scenario 4E posterior performance ================\n")
        cat("Main fit: strong spatial/covariate confounding, fixed gamma, sampled lambda, continuous-time x2\n")
        cat("Fit status counts:\n")
        print(status_counts, row.names = FALSE)
        cat("\nConfounding calibration:\n")
        print(confounding_aggregate[, intersect(c(
            "subset", "n_reps", "confounding_strength_avg",
            "x2_phi_abs_cell_cor_avg", "x2_phi_abs_area_mean_cor_avg",
            "reference_x2_phi_abs_cell_cor_avg", "reference_x2_phi_abs_area_mean_cor_avg",
            "x2_empirical_ar1_mean_avg", "x2_binary_like_prop_max"
        ), names(confounding_aggregate)), drop = FALSE], row.names = FALSE)
        cat("\nPerformance by subset:\n")
        print(performance_aggregate[, intersect(c(
            "subset", "n_reps", "beta1_mean_avg", "beta2_mean_avg", "beta2_mean_sd",
            "beta1_coverage", "beta2_coverage", "r_mean_avg",
            "log_lambda_rmse_avg", "cor_log_lambda_avg", "phi_rmse_avg", "phi_cor_avg",
            "beta_guard_total", "kappa_guard_total", "lambda_output_guard_total"
        ), names(performance_aggregate)), drop = FALSE], row.names = FALSE)
        cat("\nBeta aggregate:\n")
        print(beta_aggregate, row.names = FALSE)
        cat("\nBeta-confounding diagnostics:\n")
        print(beta_confounding_diagnostics, row.names = FALSE)
        cat("====================================================================\n\n")
        cat("Saved Scenario 4E posterior performance RDS: ", .s4e_norm_path(file.path(out_base, "scenario4e_posterior_performance_results.rds")), "\n", sep = "")
        cat("Tables: ", .s4e_norm_path(table_dir), "\n", sep = "")
        cat("Figures: ", .s4e_norm_path(fig_dir), "\n", sep = "")
    }

    out
}

summarize_scenario4e_posterior_performance <- run_scenario4e_posterior_performance

## Auto-run when sourced interactively from project root.
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
    s4e_perf <- run_scenario4e_posterior_performance(root = ".", make_plots = TRUE, verbose = TRUE)
    assign("s4e_perf", s4e_perf, envir = .GlobalEnv)
}
