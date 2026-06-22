## ============================================================================
## VERSION: S4B_CONTINUOUS_X2_POSTERIOR_PERFORMANCE_V1_2026_06_22
## run_scenario4b_posterior_performance_continuous_x2.R
##
## Posterior-performance and failure-mode summaries for Scenario 4B:
## low and heterogeneous exposure stress test with fixed gamma and continuous-time x2.
##
## Design goal
##   Scenario 4B is not just another sparse-count scenario.  It is designed to
##   evaluate whether low and heterogeneous known exposure weakens recovery,
##   especially in low-exposure areas.  This script therefore summarizes:
##     1. fit-status classification: stable / soft_warning / numerical_instability;
##     2. all, stable-only, stable+soft-warning, and unstable subsets;
##     3. regression recovery, especially beta1 vs beta2;
##     4. global lambda and phi recovery;
##     5. low-exposure vs high-exposure recovery diagnostics;
##     6. comparison to S3 control and S4A when those summaries are available;
##     7. oracle-lambda and oracle-lambda+phi diagnostic ladder when available.
##
## Typical use from project root:
##   source("run_scenario4b_posterior_performance_continuous_x2.R")
##
## Or programmatically:
##   out <- summarize_scenario4b_posterior_performance(root = ".")
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s4b_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s4b_read_required_csv <- function(path) {
    if (!file.exists(path)) {
        stop("Required CSV not found: ", path, call. = FALSE)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
}

.s4b_read_optional_csv <- function(path) {
    if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else NULL
}

.s4b_num <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    suppressWarnings(as.numeric(x[[1L]]))
}

.s4b_col <- function(d, nm, default = NA_real_) {
    if (!nm %in% names(d)) return(rep(default, nrow(d)))
    d[[nm]]
}

.s4b_safe_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

.s4b_safe_sd <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (sum(is.finite(x)) <= 1L) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
}

.s4b_safe_sum <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!any(is.finite(x))) return(0)
    sum(x, na.rm = TRUE)
}

.s4b_safe_cor <- function(x, y) {
    out <- suppressWarnings(stats::cor(as.vector(x), as.vector(y), use = "complete.obs"))
    if (!is.finite(out)) NA_real_ else out
}

.s4b_safe_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(rep(NA_real_, length(probs)))
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))
}

.s4b_fix_summary_aliases <- function(d) {
    ## Make the script robust to either source-data-manifest names or fit-summary names.
    alias_pairs <- list(
        count_low_exposure_mean_count = "lowest_exposure_group_mean_count",
        count_low_exposure_zero_prop = "lowest_exposure_group_zero_prop",
        count_high_exposure_mean_count = "highest_exposure_group_mean_count",
        count_high_exposure_zero_prop = "highest_exposure_group_zero_prop",
        target_mean_multiplier = "target_mean_exposure_multiplier",
        realized_mean_multiplier = "realized_mean_exposure_multiplier"
    )
    for (to in names(alias_pairs)) {
        from <- alias_pairs[[to]]
        if (!to %in% names(d) && from %in% names(d)) d[[to]] <- d[[from]]
    }
    if (!"mean_count" %in% names(d) && "observed_mean_count" %in% names(d)) d$mean_count <- d$observed_mean_count
    if (!"zero_prop" %in% names(d) && "observed_zero_prop" %in% names(d)) d$zero_prop <- d$observed_zero_prop
    if (!"total_count" %in% names(d) && "observed_total_count" %in% names(d)) d$total_count <- d$observed_total_count
    if (!"max_count" %in% names(d) && "observed_max_count" %in% names(d)) d$max_count <- d$observed_max_count
    d
}

.s4b_classify_one <- function(row) {
    beta_guard <- .s4b_num(row$s4b_beta_guard_count, 0)
    kappa_guard <- .s4b_num(row$s4b_kappa_guard_count, 0)
    lambda_in_guard <- .s4b_num(row$s4b_lambda_input_guard_count, 0)
    lambda_out_guard <- .s4b_num(row$s4b_lambda_output_guard_count, 0)
    beta0_mean <- .s4b_num(row$beta0_mean, NA_real_)
    lambda_rmse <- .s4b_num(row$lambda_rmse, NA_real_)
    log_lambda_rmse <- .s4b_num(row$log_lambda_rmse, NA_real_)
    phi_rmse <- .s4b_num(row$phi_rmse, NA_real_)

    severe <- (!is.finite(beta0_mean)) ||
        (!is.finite(lambda_rmse)) ||
        (!is.finite(log_lambda_rmse)) ||
        abs(beta0_mean) > 30 ||
        lambda_rmse > 10 ||
        log_lambda_rmse > 5 ||
        (is.finite(phi_rmse) && phi_rmse > 100) ||
        kappa_guard > 0 ||
        lambda_in_guard > 0 ||
        lambda_out_guard > 1000 ||
        beta_guard > 1000

    if (isTRUE(severe)) return("numerical_instability")

    warn <- beta_guard > 0 || lambda_out_guard > 0
    if (isTRUE(warn)) return("soft_warning")

    "stable"
}

.s4b_add_fit_status <- function(d) {
    d <- .s4b_fix_summary_aliases(d)
    if (!"s4b_beta_guard_count" %in% names(d)) d$s4b_beta_guard_count <- 0
    if (!"s4b_kappa_guard_count" %in% names(d)) d$s4b_kappa_guard_count <- 0
    if (!"s4b_lambda_input_guard_count" %in% names(d)) d$s4b_lambda_input_guard_count <- 0
    if (!"s4b_lambda_output_guard_count" %in% names(d)) d$s4b_lambda_output_guard_count <- 0

    d$fit_status <- vapply(seq_len(nrow(d)), function(i) .s4b_classify_one(d[i, , drop = FALSE]), character(1))
    d$nonfailed <- d$fit_status %in% c("stable", "soft_warning")
    d$stable_fit <- d$fit_status == "stable"
    d
}

.s4b_make_beta_long <- function(d, scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100") {
    params <- c("beta0", "beta1", "beta2")
    out <- list()
    for (p in params) {
        truth <- .s4b_col(d, paste0(p, "_true"))
        meanv <- .s4b_col(d, paste0(p, "_mean"))
        sdv <- .s4b_col(d, paste0(p, "_sd"))
        q025 <- .s4b_col(d, paste0(p, "_q025"))
        q50 <- .s4b_col(d, paste0(p, "_q50"))
        q975 <- .s4b_col(d, paste0(p, "_q975"))
        if (all(is.na(q50)) && paste0(p, "_median") %in% names(d)) q50 <- d[[paste0(p, "_median")]]
        bias <- .s4b_col(d, paste0(p, "_bias"), meanv - truth)
        covered <- .s4b_col(d, paste0(p, "_covered"), as.integer(q025 <= truth & truth <= q975))
        out[[p]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(d$rep_id),
            parameter = p,
            truth = as.numeric(truth),
            post_mean = as.numeric(meanv),
            post_sd = as.numeric(sdv),
            q025 = as.numeric(q025),
            q50 = as.numeric(q50),
            q975 = as.numeric(q975),
            bias = as.numeric(bias),
            covered = as.integer(covered),
            covered_95 = as.integer(covered),
            fit_status = d$fit_status,
            nonfailed = d$nonfailed,
            stable_fit = d$stable_fit,
            stringsAsFactors = FALSE
        )
    }
    ans <- do.call(rbind, out)
    rownames(ans) <- NULL
    ans
}

.s4b_make_lambda_summary <- function(d, scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100") {
    cols <- c(
        "rep_id", "mean_count", "zero_prop", "exposure_mean_ratio", "exposure_cv",
        "lambda_rmse", "lambda_mae", "lambda_coverage_95", "log_lambda_rmse",
        "log_lambda_mae", "log_lambda_coverage_95", "cor_lambda", "cor_log_lambda",
        "rmse_delta_log_lambda", "mae_delta_log_lambda", "cor_delta_log_lambda",
        "fit_status", "nonfailed", "stable_fit"
    )
    keep <- intersect(cols, names(d))
    out <- d[, keep, drop = FALSE]
    out$scenario_id <- scenario_id
    out <- out[, c("scenario_id", setdiff(names(out), "scenario_id")), drop = FALSE]
    out
}

.s4b_make_exposure_group_long <- function(d, scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100") {
    groups <- c("low_exposure", "middle_exposure", "high_exposure")
    out <- list()
    for (g in groups) {
        prefix <- paste0("_", g, "_")
        count_mean <- .s4b_col(d, paste0("count_", g, "_mean_count"))
        count_zero <- .s4b_col(d, paste0("count_", g, "_zero_prop"))
        count_n <- .s4b_col(d, paste0("count_", g, "_n_areas"))
        count_mult <- .s4b_col(d, paste0("count_", g, "_mean_multiplier"))
        lam_rmse <- .s4b_col(d, paste0("lambda_", g, "_lambda_rmse"))
        log_lam_rmse <- .s4b_col(d, paste0("lambda_", g, "_log_lambda_rmse"))
        cor_log_lam <- .s4b_col(d, paste0("lambda_", g, "_cor_log_lambda"))
        phi_rmse <- .s4b_col(d, paste0("phi_", g, "_phi_rmse"))
        phi_mae <- .s4b_col(d, paste0("phi_", g, "_phi_mae"))
        phi_cor <- .s4b_col(d, paste0("phi_", g, "_phi_cor"))
        ## Backward-compatible aliases for low/high.
        if (g == "low_exposure") {
            if (all(is.na(count_mean))) count_mean <- .s4b_col(d, "lowest_exposure_group_mean_count")
            if (all(is.na(count_zero))) count_zero <- .s4b_col(d, "lowest_exposure_group_zero_prop")
        }
        if (g == "high_exposure") {
            if (all(is.na(count_mean))) count_mean <- .s4b_col(d, "highest_exposure_group_mean_count")
            if (all(is.na(count_zero))) count_zero <- .s4b_col(d, "highest_exposure_group_zero_prop")
        }
        out[[g]] <- data.frame(
            scenario_id = scenario_id,
            rep_id = as.integer(d$rep_id),
            exposure_group = g,
            n_areas = as.numeric(count_n),
            mean_multiplier = as.numeric(count_mult),
            mean_count = as.numeric(count_mean),
            zero_prop = as.numeric(count_zero),
            lambda_rmse = as.numeric(lam_rmse),
            log_lambda_rmse = as.numeric(log_lam_rmse),
            cor_log_lambda = as.numeric(cor_log_lam),
            phi_rmse = as.numeric(phi_rmse),
            phi_mae = as.numeric(phi_mae),
            phi_cor = as.numeric(phi_cor),
            fit_status = d$fit_status,
            nonfailed = d$nonfailed,
            stable_fit = d$stable_fit,
            stringsAsFactors = FALSE
        )
    }
    ans <- do.call(rbind, out)
    rownames(ans) <- NULL
    ans
}

.s4b_parameter_aggregate <- function(beta_df, subset_name = "all") {
    if (is.null(beta_df) || nrow(beta_df) == 0L) return(data.frame())
    out <- do.call(rbind, lapply(split(beta_df, beta_df$parameter), function(d) {
        data.frame(
            subset = subset_name,
            parameter = unique(d$parameter),
            n_reps = length(unique(d$rep_id)),
            mean_truth = .s4b_safe_mean(d$truth),
            mean_post_mean = .s4b_safe_mean(d$post_mean),
            sd_post_mean_across_reps = .s4b_safe_sd(d$post_mean),
            mean_bias = .s4b_safe_mean(d$bias),
            median_bias = stats::median(d$bias, na.rm = TRUE),
            rmse_bias = sqrt(mean(d$bias^2, na.rm = TRUE)),
            mean_abs_bias = .s4b_safe_mean(abs(d$bias)),
            coverage = .s4b_safe_mean(d$covered),
            mean_post_sd = .s4b_safe_mean(d$post_sd),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

.s4b_diagnostics_aggregate <- function(diag_df, subset_name = "all") {
    if (is.null(diag_df) || nrow(diag_df) == 0L) return(data.frame())
    data.frame(
        subset = subset_name,
        n_reps = nrow(diag_df),
        mean_count_avg = .s4b_safe_mean(diag_df$mean_count),
        zero_prop_avg = .s4b_safe_mean(diag_df$zero_prop),
        exposure_mean_ratio_avg = .s4b_safe_mean(.s4b_col(diag_df, "exposure_mean_ratio")),
        exposure_cv_avg = .s4b_safe_mean(.s4b_col(diag_df, "exposure_cv")),
        low_exp_mean_count_avg = .s4b_safe_mean(.s4b_col(diag_df, "count_low_exposure_mean_count")),
        high_exp_mean_count_avg = .s4b_safe_mean(.s4b_col(diag_df, "count_high_exposure_mean_count")),
        low_exp_zero_prop_avg = .s4b_safe_mean(.s4b_col(diag_df, "count_low_exposure_zero_prop")),
        high_exp_zero_prop_avg = .s4b_safe_mean(.s4b_col(diag_df, "count_high_exposure_zero_prop")),
        beta1_mean_avg = .s4b_safe_mean(diag_df$beta1_mean),
        beta1_mean_sd = .s4b_safe_sd(diag_df$beta1_mean),
        beta2_mean_avg = .s4b_safe_mean(diag_df$beta2_mean),
        beta2_mean_sd = .s4b_safe_sd(diag_df$beta2_mean),
        beta0_coverage = .s4b_safe_mean(.s4b_col(diag_df, "beta0_covered")),
        beta1_coverage = .s4b_safe_mean(.s4b_col(diag_df, "beta1_covered")),
        beta2_coverage = .s4b_safe_mean(.s4b_col(diag_df, "beta2_covered")),
        r_mean_avg = .s4b_safe_mean(.s4b_col(diag_df, "r_mean")),
        phi_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "phi_rmse")),
        phi_cor_avg = .s4b_safe_mean(.s4b_col(diag_df, "phi_cor")),
        lambda_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "lambda_rmse")),
        log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "log_lambda_rmse")),
        cor_log_lambda_avg = .s4b_safe_mean(.s4b_col(diag_df, "cor_log_lambda")),
        low_exp_log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "lambda_low_exposure_log_lambda_rmse")),
        high_exp_log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "lambda_high_exposure_log_lambda_rmse")),
        low_minus_high_log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "lambda_low_exposure_log_lambda_rmse") - .s4b_col(diag_df, "lambda_high_exposure_log_lambda_rmse")),
        low_exp_phi_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "phi_low_exposure_phi_rmse")),
        high_exp_phi_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "phi_high_exposure_phi_rmse")),
        low_minus_high_phi_rmse_avg = .s4b_safe_mean(.s4b_col(diag_df, "phi_low_exposure_phi_rmse") - .s4b_col(diag_df, "phi_high_exposure_phi_rmse")),
        beta_guard_total = .s4b_safe_sum(.s4b_col(diag_df, "s4b_beta_guard_count", 0)),
        kappa_guard_total = .s4b_safe_sum(.s4b_col(diag_df, "s4b_kappa_guard_count", 0)),
        lambda_output_guard_total = .s4b_safe_sum(.s4b_col(diag_df, "s4b_lambda_output_guard_count", 0)),
        stringsAsFactors = FALSE
    )
}

.s4b_exposure_group_aggregate <- function(group_df, subset_name = "all") {
    if (is.null(group_df) || nrow(group_df) == 0L) return(data.frame())
    out <- do.call(rbind, lapply(split(group_df, group_df$exposure_group), function(d) {
        data.frame(
            subset = subset_name,
            exposure_group = unique(d$exposure_group),
            n_reps = length(unique(d$rep_id)),
            mean_count_avg = .s4b_safe_mean(d$mean_count),
            zero_prop_avg = .s4b_safe_mean(d$zero_prop),
            lambda_rmse_avg = .s4b_safe_mean(d$lambda_rmse),
            log_lambda_rmse_avg = .s4b_safe_mean(d$log_lambda_rmse),
            cor_log_lambda_avg = .s4b_safe_mean(d$cor_log_lambda),
            phi_rmse_avg = .s4b_safe_mean(d$phi_rmse),
            phi_mae_avg = .s4b_safe_mean(d$phi_mae),
            phi_cor_avg = .s4b_safe_mean(d$phi_cor),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

.s4b_low_high_contrast_by_rep <- function(d) {
    data.frame(
        rep_id = as.integer(d$rep_id),
        fit_status = d$fit_status,
        nonfailed = d$nonfailed,
        stable_fit = d$stable_fit,
        mean_count = d$mean_count,
        zero_prop = d$zero_prop,
        low_mean_count = .s4b_col(d, "count_low_exposure_mean_count"),
        high_mean_count = .s4b_col(d, "count_high_exposure_mean_count"),
        low_zero_prop = .s4b_col(d, "count_low_exposure_zero_prop"),
        high_zero_prop = .s4b_col(d, "count_high_exposure_zero_prop"),
        low_log_lambda_rmse = .s4b_col(d, "lambda_low_exposure_log_lambda_rmse"),
        high_log_lambda_rmse = .s4b_col(d, "lambda_high_exposure_log_lambda_rmse"),
        low_minus_high_log_lambda_rmse = .s4b_col(d, "lambda_low_exposure_log_lambda_rmse") - .s4b_col(d, "lambda_high_exposure_log_lambda_rmse"),
        low_phi_rmse = .s4b_col(d, "phi_low_exposure_phi_rmse"),
        high_phi_rmse = .s4b_col(d, "phi_high_exposure_phi_rmse"),
        low_minus_high_phi_rmse = .s4b_col(d, "phi_low_exposure_phi_rmse") - .s4b_col(d, "phi_high_exposure_phi_rmse"),
        stringsAsFactors = FALSE
    )
}

.s4b_low_high_contrast_aggregate <- function(contrast_df, subset_name = "all") {
    if (is.null(contrast_df) || nrow(contrast_df) == 0L) return(data.frame())
    data.frame(
        subset = subset_name,
        n_reps = nrow(contrast_df),
        low_mean_count_avg = .s4b_safe_mean(contrast_df$low_mean_count),
        high_mean_count_avg = .s4b_safe_mean(contrast_df$high_mean_count),
        low_zero_prop_avg = .s4b_safe_mean(contrast_df$low_zero_prop),
        high_zero_prop_avg = .s4b_safe_mean(contrast_df$high_zero_prop),
        low_log_lambda_rmse_avg = .s4b_safe_mean(contrast_df$low_log_lambda_rmse),
        high_log_lambda_rmse_avg = .s4b_safe_mean(contrast_df$high_log_lambda_rmse),
        low_minus_high_log_lambda_rmse_avg = .s4b_safe_mean(contrast_df$low_minus_high_log_lambda_rmse),
        low_phi_rmse_avg = .s4b_safe_mean(contrast_df$low_phi_rmse),
        high_phi_rmse_avg = .s4b_safe_mean(contrast_df$high_phi_rmse),
        low_minus_high_phi_rmse_avg = .s4b_safe_mean(contrast_df$low_minus_high_phi_rmse),
        stringsAsFactors = FALSE
    )
}

.s4b_build_comparison_table <- function(s4b_diag_agg,
                                        s4b_status_counts,
                                        s3_control_summary = NULL,
                                        s4a_performance_subset = NULL,
                                        s4a_status_counts = NULL,
                             oracle_lambda_summary = NULL,
                             oracle_lamphi_summary = NULL,
                             oracle_lambda_detail = NULL,
                             oracle_lamphi_detail = NULL) {
    rows <- list()

    if (!is.null(s3_control_summary) && nrow(s3_control_summary) > 0L) {
        rows[["s3"]] <- data.frame(
            scenario = "S3 control fixed gamma",
            subset = "all_control_reps",
            n_reps = nrow(s3_control_summary),
            mean_count_avg = .s4b_safe_mean(s3_control_summary$mean_count),
            zero_prop_avg = .s4b_safe_mean(s3_control_summary$zero_prop),
            beta1_mean_avg = .s4b_safe_mean(s3_control_summary$beta1_mean),
            beta2_mean_avg = .s4b_safe_mean(s3_control_summary$beta2_mean),
            beta2_mean_sd = .s4b_safe_sd(s3_control_summary$beta2_mean),
            lambda_rmse_avg = .s4b_safe_mean(s3_control_summary$lambda_rmse),
            log_lambda_rmse_avg = .s4b_safe_mean(s3_control_summary$log_lambda_rmse),
            cor_log_lambda_avg = .s4b_safe_mean(s3_control_summary$cor_log_lambda),
            numerical_instability_rate = 0,
            conclusion = "Non-sparse reference setting; same stabilized machinery should be stable.",
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(s4a_performance_subset) && nrow(s4a_performance_subset) > 0L) {
        s4a_row <- s4a_performance_subset[s4a_performance_subset$subset %in% c("stable_plus_soft_warning", "stable_plus_soft", "stable"), , drop = FALSE]
        if (nrow(s4a_row) > 1L) {
            pref <- match(c("stable_plus_soft_warning", "stable_plus_soft", "stable"), s4a_row$subset)
            pref <- pref[is.finite(pref)][1L]
            if (!is.na(pref)) s4a_row <- s4a_row[pref, , drop = FALSE]
        }
        if (nrow(s4a_row) == 0L) s4a_row <- s4a_performance_subset[1L, , drop = FALSE]
        inst_rate <- NA_real_
        if (!is.null(s4a_status_counts) && nrow(s4a_status_counts) > 0L && "fit_status" %in% names(s4a_status_counts)) {
            idx <- s4a_status_counts$fit_status == "numerical_instability"
            if (any(idx)) inst_rate <- sum(s4a_status_counts$n_reps[idx], na.rm = TRUE) / sum(s4a_status_counts$n_reps, na.rm = TRUE)
        }
        rows[["s4a"]] <- data.frame(
            scenario = "S4A sparse counts fixed gamma, continuous x2",
            subset = as.character(s4a_row$subset[1]),
            n_reps = as.integer(s4a_row$n_reps[1]),
            mean_count_avg = .s4b_num(s4a_row$mean_count_avg),
            zero_prop_avg = .s4b_num(s4a_row$zero_prop_avg),
            beta1_mean_avg = .s4b_num(s4a_row$beta1_mean_avg),
            beta2_mean_avg = .s4b_num(s4a_row$beta2_mean_avg),
            beta2_mean_sd = .s4b_num(s4a_row$beta2_between_rep_sd, .s4b_num(s4a_row$beta2_mean_sd)),
            lambda_rmse_avg = .s4b_num(s4a_row$lambda_rmse_avg, .s4b_num(s4a_row$lambda_rmse_mean)),
            log_lambda_rmse_avg = .s4b_num(s4a_row$log_lambda_rmse_avg, .s4b_num(s4a_row$log_lambda_rmse_mean)),
            cor_log_lambda_avg = .s4b_num(s4a_row$cor_log_lambda_avg, .s4b_num(s4a_row$cor_log_lambda_mean)),
            numerical_instability_rate = inst_rate,
            conclusion = "Global sparse-count stress under the continuous-time x2 design; useful reference but not the same stress source.",
            stringsAsFactors = FALSE
        )
    }

    s4b_row <- s4b_diag_agg[s4b_diag_agg$subset == "stable_plus_soft_warning", , drop = FALSE]
    if (nrow(s4b_row) == 0L) s4b_row <- s4b_diag_agg[1L, , drop = FALSE]
    inst_rate <- NA_real_
    idx <- s4b_status_counts$fit_status == "numerical_instability"
    if (any(idx)) inst_rate <- sum(s4b_status_counts$n_reps[idx], na.rm = TRUE) / sum(s4b_status_counts$n_reps, na.rm = TRUE)
    rows[["s4b"]] <- data.frame(
        scenario = "S4B low heterogeneous exposure fixed gamma, continuous x2",
        subset = as.character(s4b_row$subset[1]),
        n_reps = as.integer(s4b_row$n_reps[1]),
        mean_count_avg = .s4b_num(s4b_row$mean_count_avg),
        zero_prop_avg = .s4b_num(s4b_row$zero_prop_avg),
        beta1_mean_avg = .s4b_num(s4b_row$beta1_mean_avg),
        beta2_mean_avg = .s4b_num(s4b_row$beta2_mean_avg),
        beta2_mean_sd = .s4b_num(s4b_row$beta2_mean_sd),
        lambda_rmse_avg = .s4b_num(s4b_row$lambda_rmse_avg),
        log_lambda_rmse_avg = .s4b_num(s4b_row$log_lambda_rmse_avg),
        cor_log_lambda_avg = .s4b_num(s4b_row$cor_log_lambda_avg),
        numerical_instability_rate = inst_rate,
        conclusion = "Exposure imbalance stress; low-exposure areas show weaker latent-risk recovery and difficult reps may fail.",
        stringsAsFactors = FALSE
    )

    ans <- do.call(rbind, rows)
    rownames(ans) <- NULL
    ans
}

.s4b_open_plot_device <- function(file, width = 8, height = 5, res = 180) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
        grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
        grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
        stop("Unsupported plot extension: ", ext, call. = FALSE)
    }
}

.s4b_close_plot_device <- function() grDevices::dev.off()

.s4b_make_plots <- function(diag, contrast, fig_dir, plot_format = c("pdf", "png")) {
    plot_format <- match.arg(plot_format)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    ext <- plot_format
    files <- list(
        beta2_by_status = file.path(fig_dir, paste0("s4b_beta2_mean_by_fit_status.", ext)),
        log_lambda_rmse_by_status = file.path(fig_dir, paste0("s4b_log_lambda_rmse_by_fit_status.", ext)),
        low_high_lambda_rmse = file.path(fig_dir, paste0("s4b_low_high_log_lambda_rmse_by_rep.", ext)),
        guard_counts = file.path(fig_dir, paste0("s4b_guard_counts_by_rep.", ext))
    )
    pch_map <- c(stable = 19, soft_warning = 17, numerical_instability = 4)
    pch_vals <- pch_map[diag$fit_status]
    pch_vals[is.na(pch_vals)] <- 1

    .s4b_open_plot_device(files$beta2_by_status, width = 8, height = 5)
    graphics::plot(diag$rep_id, diag$beta2_mean,
                   pch = pch_vals, xlab = "Replicate", ylab = "Posterior mean of beta2",
                   main = "S4B beta2 recovery by fit status")
    if ("beta2_true" %in% names(diag)) graphics::abline(h = unique(diag$beta2_true)[1L], lty = 2)
    graphics::legend("topright", legend = names(pch_map), pch = pch_map, bty = "n")
    .s4b_close_plot_device()

    .s4b_open_plot_device(files$log_lambda_rmse_by_status, width = 8, height = 5)
    y <- pmin(diag$log_lambda_rmse, 10)
    graphics::plot(diag$rep_id, y,
                   pch = pch_vals, xlab = "Replicate", ylab = "log-lambda RMSE (capped at 10)",
                   main = "S4B latent temporal-risk recovery")
    graphics::legend("topright", legend = names(pch_map), pch = pch_map, bty = "n")
    .s4b_close_plot_device()

    .s4b_open_plot_device(files$low_high_lambda_rmse, width = 8, height = 5)
    yy_low <- pmin(contrast$low_log_lambda_rmse, 10)
    yy_high <- pmin(contrast$high_log_lambda_rmse, 10)
    ylim <- range(c(yy_low, yy_high), na.rm = TRUE)
    graphics::plot(contrast$rep_id, yy_low, pch = 19, ylim = ylim,
                   xlab = "Replicate", ylab = "log-lambda RMSE (capped at 10)",
                   main = "S4B low- vs high-exposure lambda recovery")
    graphics::points(contrast$rep_id, yy_high, pch = 17)
    graphics::legend("topright", legend = c("Low exposure", "High exposure"), pch = c(19, 17), bty = "n")
    .s4b_close_plot_device()

    .s4b_open_plot_device(files$guard_counts, width = 8, height = 5)
    total_guards <- .s4b_col(diag, "s4b_beta_guard_count", 0) +
        .s4b_col(diag, "s4b_kappa_guard_count", 0) +
        .s4b_col(diag, "s4b_lambda_output_guard_count", 0)
    graphics::plot(diag$rep_id, log10(total_guards + 1), type = "h", lwd = 3,
                   xlab = "Replicate", ylab = "log10(total guards + 1)",
                   main = "S4B numerical guard activity")
    graphics::points(diag$rep_id, log10(total_guards + 1), pch = pch_vals)
    .s4b_close_plot_device()

    invisible(files)
}


.s4b_summarize_oracle_diag <- function(d, label, conclusion) {
    if (is.null(d) || nrow(d) == 0L) return(data.frame())
    data.frame(
        diagnostic_step = label,
        n_reps = nrow(d),
        mean_count_avg = .s4b_safe_mean(.s4b_col(d, "mean_count")),
        zero_prop_avg = .s4b_safe_mean(.s4b_col(d, "zero_prop")),
        exposure_mean_ratio_avg = .s4b_safe_mean(.s4b_col(d, "exposure_mean_ratio")),
        exposure_cv_avg = .s4b_safe_mean(.s4b_col(d, "exposure_cv")),
        beta0_mean_avg = .s4b_safe_mean(.s4b_col(d, "beta0_mean")),
        beta0_coverage = .s4b_safe_mean(.s4b_col(d, "beta0_covered")),
        beta1_mean_avg = .s4b_safe_mean(.s4b_col(d, "beta1_mean")),
        beta1_coverage = .s4b_safe_mean(.s4b_col(d, "beta1_covered")),
        beta2_mean_avg = .s4b_safe_mean(.s4b_col(d, "beta2_mean")),
        beta2_mean_sd = .s4b_safe_sd(.s4b_col(d, "beta2_mean")),
        beta2_coverage = .s4b_safe_mean(.s4b_col(d, "beta2_covered")),
        r_mean_avg = .s4b_safe_mean(.s4b_col(d, "r_mean")),
        phi_rmse_avg = .s4b_safe_mean(.s4b_col(d, "phi_rmse")),
        phi_cor_avg = .s4b_safe_mean(.s4b_col(d, "phi_cor")),
        lambda_rmse_avg = .s4b_safe_mean(.s4b_col(d, "lambda_rmse")),
        log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(d, "log_lambda_rmse")),
        cor_log_lambda_avg = .s4b_safe_mean(.s4b_col(d, "cor_log_lambda")),
        low_exp_log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(d, "lambda_low_exposure_log_lambda_rmse")),
        high_exp_log_lambda_rmse_avg = .s4b_safe_mean(.s4b_col(d, "lambda_high_exposure_log_lambda_rmse")),
        low_exp_phi_rmse_avg = .s4b_safe_mean(.s4b_col(d, "phi_low_exposure_phi_rmse")),
        high_exp_phi_rmse_avg = .s4b_safe_mean(.s4b_col(d, "phi_high_exposure_phi_rmse")),
        beta_guard_total = .s4b_safe_sum(if ("s4b_oracle_beta_guard_count" %in% names(d)) d$s4b_oracle_beta_guard_count else if ("s4b_oracle_lamphi_beta_guard_count" %in% names(d)) d$s4b_oracle_lamphi_beta_guard_count else rep(0, nrow(d))),
        kappa_guard_total = .s4b_safe_sum(if ("s4b_oracle_kappa_guard_count" %in% names(d)) d$s4b_oracle_kappa_guard_count else if ("s4b_oracle_lamphi_kappa_guard_count" %in% names(d)) d$s4b_oracle_lamphi_kappa_guard_count else rep(0, nrow(d))),
        conclusion = conclusion,
        stringsAsFactors = FALSE
    )
}

.s4b_build_diagnostic_ladder <- function(diag_agg,
                                         status_counts,
                                         comparison = NULL,
                                         oracle_lambda_summary = NULL,
                                         oracle_lamphi_summary = NULL) {
    rows <- list()

    main <- diag_agg[diag_agg$subset == "stable_plus_soft_warning", , drop = FALSE]
    if (nrow(main) == 0L) main <- diag_agg[1L, , drop = FALSE]
    n_stable <- sum(status_counts$n_reps[status_counts$fit_status == "stable"], na.rm = TRUE)
    n_soft <- sum(status_counts$n_reps[status_counts$fit_status == "soft_warning"], na.rm = TRUE)
    n_inst <- sum(status_counts$n_reps[status_counts$fit_status == "numerical_instability"], na.rm = TRUE)
    rows[["s4b"]] <- data.frame(
        diagnostic_step = "S4B sampled lambda, fixed gamma (stable+soft)",
        n_reps = as.integer(main$n_reps[1]),
        n_stable = n_stable,
        n_soft_warning = n_soft,
        n_numerical_instability = n_inst,
        beta1_mean_avg = .s4b_num(main$beta1_mean_avg),
        beta2_mean_avg = .s4b_num(main$beta2_mean_avg),
        beta2_mean_sd = .s4b_num(main$beta2_mean_sd),
        lambda_rmse_avg = .s4b_num(main$lambda_rmse_avg),
        log_lambda_rmse_avg = .s4b_num(main$log_lambda_rmse_avg),
        cor_log_lambda_avg = .s4b_num(main$cor_log_lambda_avg),
        low_exp_log_lambda_rmse_avg = .s4b_num(main$low_exp_log_lambda_rmse_avg),
        high_exp_log_lambda_rmse_avg = .s4b_num(main$high_exp_log_lambda_rmse_avg),
        beta_guard_total = .s4b_num(main$beta_guard_total, 0),
        kappa_guard_total = .s4b_num(main$kappa_guard_total, 0),
        lambda_output_guard_total = .s4b_num(main$lambda_output_guard_total, 0),
        conclusion = "Main exposure-stress fit summarized on stable+soft reps; failure rate is reported separately.",
        stringsAsFactors = FALSE
    )

    if (!is.null(comparison) && nrow(comparison) > 0L) {
        s3 <- comparison[comparison$scenario == "S3 control fixed gamma", , drop = FALSE]
        if (nrow(s3) > 0L) {
            rows[["s3"]] <- data.frame(
                diagnostic_step = "S3 control, fixed gamma stabilized",
                n_reps = as.integer(s3$n_reps[1]),
                n_stable = as.integer(s3$n_reps[1]),
                n_soft_warning = 0L,
                n_numerical_instability = 0L,
                beta1_mean_avg = .s4b_num(s3$beta1_mean_avg),
                beta2_mean_avg = .s4b_num(s3$beta2_mean_avg),
                beta2_mean_sd = .s4b_num(s3$beta2_mean_sd),
                lambda_rmse_avg = .s4b_num(s3$lambda_rmse_avg),
                log_lambda_rmse_avg = .s4b_num(s3$log_lambda_rmse_avg),
                cor_log_lambda_avg = .s4b_num(s3$cor_log_lambda_avg),
                low_exp_log_lambda_rmse_avg = NA_real_,
                high_exp_log_lambda_rmse_avg = NA_real_,
                beta_guard_total = 0,
                kappa_guard_total = 0,
                lambda_output_guard_total = 0,
                conclusion = "Same machinery stable in non-exposure-stressed reference setting.",
                stringsAsFactors = FALSE
            )
        }
        s4a <- comparison[comparison$scenario == "S4A sparse counts fixed gamma", , drop = FALSE]
        if (nrow(s4a) > 0L) {
            rows[["s4a"]] <- data.frame(
                diagnostic_step = "S4A sparse counts, fixed gamma",
                n_reps = as.integer(s4a$n_reps[1]),
                n_stable = NA_integer_,
                n_soft_warning = NA_integer_,
                n_numerical_instability = NA_integer_,
                beta1_mean_avg = .s4b_num(s4a$beta1_mean_avg),
                beta2_mean_avg = .s4b_num(s4a$beta2_mean_avg),
                beta2_mean_sd = .s4b_num(s4a$beta2_mean_sd),
                lambda_rmse_avg = .s4b_num(s4a$lambda_rmse_avg),
                log_lambda_rmse_avg = .s4b_num(s4a$log_lambda_rmse_avg),
                cor_log_lambda_avg = .s4b_num(s4a$cor_log_lambda_avg),
                low_exp_log_lambda_rmse_avg = NA_real_,
                high_exp_log_lambda_rmse_avg = NA_real_,
                beta_guard_total = NA_real_,
                kappa_guard_total = NA_real_,
                lambda_output_guard_total = NA_real_,
                conclusion = "Global sparse-count stress; useful reference for S4B but not the same stress source.",
                stringsAsFactors = FALSE
            )
        }
    }

    if (!is.null(oracle_lambda_summary) && nrow(oracle_lambda_summary) > 0L) {
        ol <- oracle_lambda_summary[1L, , drop = FALSE]
        rows[["oracle_lambda"]] <- data.frame(
            diagnostic_step = ol$diagnostic_step[1],
            n_reps = as.integer(ol$n_reps[1]),
            n_stable = NA_integer_,
            n_soft_warning = NA_integer_,
            n_numerical_instability = NA_integer_,
            beta1_mean_avg = .s4b_num(ol$beta1_mean_avg),
            beta2_mean_avg = .s4b_num(ol$beta2_mean_avg),
            beta2_mean_sd = .s4b_num(ol$beta2_mean_sd),
            lambda_rmse_avg = .s4b_num(ol$lambda_rmse_avg),
            log_lambda_rmse_avg = .s4b_num(ol$log_lambda_rmse_avg),
            cor_log_lambda_avg = .s4b_num(ol$cor_log_lambda_avg),
            low_exp_log_lambda_rmse_avg = .s4b_num(ol$low_exp_log_lambda_rmse_avg),
            high_exp_log_lambda_rmse_avg = .s4b_num(ol$high_exp_log_lambda_rmse_avg),
            beta_guard_total = .s4b_num(ol$beta_guard_total, 0),
            kappa_guard_total = .s4b_num(ol$kappa_guard_total, 0),
            lambda_output_guard_total = 0,
            conclusion = ol$conclusion[1],
            stringsAsFactors = FALSE
        )
    }

    if (!is.null(oracle_lamphi_summary) && nrow(oracle_lamphi_summary) > 0L) {
        op <- oracle_lamphi_summary[1L, , drop = FALSE]
        rows[["oracle_lamphi"]] <- data.frame(
            diagnostic_step = op$diagnostic_step[1],
            n_reps = as.integer(op$n_reps[1]),
            n_stable = NA_integer_,
            n_soft_warning = NA_integer_,
            n_numerical_instability = NA_integer_,
            beta1_mean_avg = .s4b_num(op$beta1_mean_avg),
            beta2_mean_avg = .s4b_num(op$beta2_mean_avg),
            beta2_mean_sd = .s4b_num(op$beta2_mean_sd),
            lambda_rmse_avg = .s4b_num(op$lambda_rmse_avg),
            log_lambda_rmse_avg = .s4b_num(op$log_lambda_rmse_avg),
            cor_log_lambda_avg = .s4b_num(op$cor_log_lambda_avg),
            low_exp_log_lambda_rmse_avg = .s4b_num(op$low_exp_log_lambda_rmse_avg),
            high_exp_log_lambda_rmse_avg = .s4b_num(op$high_exp_log_lambda_rmse_avg),
            beta_guard_total = .s4b_num(op$beta_guard_total, 0),
            kappa_guard_total = .s4b_num(op$kappa_guard_total, 0),
            lambda_output_guard_total = 0,
            conclusion = op$conclusion[1],
            stringsAsFactors = FALSE
        )
    }

    ans <- do.call(rbind, rows)
    rownames(ans) <- NULL
    ans
}

.s4b_save_tables <- function(diag,
                             beta_long,
                             lambda_summary,
                             group_long,
                             contrast,
                             out_dir,
                             s3_control_summary = NULL,
                             s4a_performance_subset = NULL,
                             s4a_status_counts = NULL,
                             oracle_lambda_summary = NULL,
                             oracle_lamphi_summary = NULL,
                             oracle_lambda_detail = NULL,
                             oracle_lamphi_detail = NULL) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    files <- list(
        beta = file.path(out_dir, "posterior_beta_summary.csv"),
        lambda = file.path(out_dir, "posterior_lambda_path_recovery.csv"),
        latent = file.path(out_dir, "posterior_latent_risk_performance.csv"),
        diagnostics = file.path(out_dir, "posterior_performance_diagnostics.csv"),
        fit_status = file.path(out_dir, "scenario4b_fit_status_by_rep.csv"),
        fit_status_counts = file.path(out_dir, "scenario4b_fit_status_counts.csv"),
        beta_aggregate = file.path(out_dir, "scenario4b_beta_recovery_by_subset.csv"),
        performance_aggregate = file.path(out_dir, "scenario4b_performance_by_subset.csv"),
        exposure_group_by_rep = file.path(out_dir, "scenario4b_exposure_group_recovery_by_rep.csv"),
        exposure_group_aggregate = file.path(out_dir, "scenario4b_exposure_group_recovery_by_subset.csv"),
        low_high_contrast_by_rep = file.path(out_dir, "scenario4b_low_high_exposure_contrast_by_rep.csv"),
        low_high_contrast_aggregate = file.path(out_dir, "scenario4b_low_high_exposure_contrast_by_subset.csv"),
        guard_summary = file.path(out_dir, "scenario4b_guard_summary_by_subset.csv"),
        comparison = file.path(out_dir, "scenario4b_comparison_to_s3_s4a_summary.csv"),
        oracle_lambda = file.path(out_dir, "scenario4b_oracle_lambda_diagnostic_summary.csv"),
        oracle_lambda_phi = file.path(out_dir, "scenario4b_oracle_lambda_phi_diagnostic_summary.csv"),
        oracle_lambda_detail = file.path(out_dir, "scenario4b_oracle_lambda_diagnostic_detail.csv"),
        oracle_lambda_phi_detail = file.path(out_dir, "scenario4b_oracle_lambda_phi_diagnostic_detail.csv"),
        diagnostic_ladder = file.path(out_dir, "scenario4b_diagnostic_ladder_summary.csv")
    )

    utils::write.csv(beta_long, files$beta, row.names = FALSE)
    utils::write.csv(lambda_summary, files$lambda, row.names = FALSE)
    utils::write.csv(lambda_summary, files$latent, row.names = FALSE)
    utils::write.csv(diag, files$diagnostics, row.names = FALSE)

    status_cols <- intersect(c(
        "scenario_id", "rep_id", "fit_status", "nonfailed", "stable_fit",
        "mean_count", "zero_prop", "exposure_mean_ratio", "exposure_cv",
        "count_low_exposure_mean_count", "count_low_exposure_zero_prop",
        "count_high_exposure_mean_count", "count_high_exposure_zero_prop",
        "beta0_mean", "beta1_mean", "beta2_mean", "r_mean",
        "phi_rmse", "phi_cor", "lambda_rmse", "log_lambda_rmse", "cor_log_lambda",
        "lambda_low_exposure_log_lambda_rmse", "lambda_high_exposure_log_lambda_rmse",
        "phi_low_exposure_phi_rmse", "phi_high_exposure_phi_rmse",
        "s4b_beta_guard_count", "s4b_kappa_guard_count",
        "s4b_lambda_input_guard_count", "s4b_lambda_output_guard_count"
    ), names(diag))
    utils::write.csv(diag[, status_cols, drop = FALSE], files$fit_status, row.names = FALSE)

    status_counts <- as.data.frame(table(diag$fit_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("fit_status", "n_reps")
    status_counts$prop_reps <- status_counts$n_reps / sum(status_counts$n_reps)
    utils::write.csv(status_counts, files$fit_status_counts, row.names = FALSE)

    beta_agg <- do.call(rbind, list(
        .s4b_parameter_aggregate(beta_long, "all_sampled_lambda"),
        .s4b_parameter_aggregate(beta_long[beta_long$fit_status == "stable", , drop = FALSE], "stable_only"),
        .s4b_parameter_aggregate(beta_long[beta_long$fit_status %in% c("stable", "soft_warning"), , drop = FALSE], "stable_plus_soft_warning"),
        .s4b_parameter_aggregate(beta_long[beta_long$fit_status == "numerical_instability", , drop = FALSE], "numerical_instability")
    ))
    utils::write.csv(beta_agg, files$beta_aggregate, row.names = FALSE)

    diag_agg <- do.call(rbind, list(
        .s4b_diagnostics_aggregate(diag, "all_sampled_lambda"),
        .s4b_diagnostics_aggregate(diag[diag$fit_status == "stable", , drop = FALSE], "stable_only"),
        .s4b_diagnostics_aggregate(diag[diag$fit_status %in% c("stable", "soft_warning"), , drop = FALSE], "stable_plus_soft_warning"),
        .s4b_diagnostics_aggregate(diag[diag$fit_status == "numerical_instability", , drop = FALSE], "numerical_instability")
    ))
    utils::write.csv(diag_agg, files$performance_aggregate, row.names = FALSE)

    group_agg <- do.call(rbind, list(
        .s4b_exposure_group_aggregate(group_long, "all_sampled_lambda"),
        .s4b_exposure_group_aggregate(group_long[group_long$fit_status == "stable", , drop = FALSE], "stable_only"),
        .s4b_exposure_group_aggregate(group_long[group_long$fit_status %in% c("stable", "soft_warning"), , drop = FALSE], "stable_plus_soft_warning"),
        .s4b_exposure_group_aggregate(group_long[group_long$fit_status == "numerical_instability", , drop = FALSE], "numerical_instability")
    ))
    utils::write.csv(group_long, files$exposure_group_by_rep, row.names = FALSE)
    utils::write.csv(group_agg, files$exposure_group_aggregate, row.names = FALSE)

    contrast_agg <- do.call(rbind, list(
        .s4b_low_high_contrast_aggregate(contrast, "all_sampled_lambda"),
        .s4b_low_high_contrast_aggregate(contrast[contrast$fit_status == "stable", , drop = FALSE], "stable_only"),
        .s4b_low_high_contrast_aggregate(contrast[contrast$fit_status %in% c("stable", "soft_warning"), , drop = FALSE], "stable_plus_soft_warning"),
        .s4b_low_high_contrast_aggregate(contrast[contrast$fit_status == "numerical_instability", , drop = FALSE], "numerical_instability")
    ))
    utils::write.csv(contrast, files$low_high_contrast_by_rep, row.names = FALSE)
    utils::write.csv(contrast_agg, files$low_high_contrast_aggregate, row.names = FALSE)

    guard_summary <- diag_agg[, intersect(c(
        "subset", "n_reps", "beta_guard_total", "kappa_guard_total", "lambda_output_guard_total"
    ), names(diag_agg)), drop = FALSE]
    utils::write.csv(guard_summary, files$guard_summary, row.names = FALSE)

    comp <- .s4b_build_comparison_table(
        s4b_diag_agg = diag_agg,
        s4b_status_counts = status_counts,
        s3_control_summary = s3_control_summary,
        s4a_performance_subset = s4a_performance_subset,
        s4a_status_counts = s4a_status_counts,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        oracle_lambda_detail = oracle_lambda_detail,
        oracle_lamphi_detail = oracle_lamphi_detail
    )
    utils::write.csv(comp, files$comparison, row.names = FALSE)

    if (!is.null(oracle_lambda_detail) && nrow(oracle_lambda_detail) > 0L) {
        utils::write.csv(oracle_lambda_detail, files$oracle_lambda_detail, row.names = FALSE)
    }
    if (!is.null(oracle_lamphi_detail) && nrow(oracle_lamphi_detail) > 0L) {
        utils::write.csv(oracle_lamphi_detail, files$oracle_lambda_phi_detail, row.names = FALSE)
    }
    if (!is.null(oracle_lambda_summary) && nrow(oracle_lambda_summary) > 0L) {
        utils::write.csv(oracle_lambda_summary, files$oracle_lambda, row.names = FALSE)
    }
    if (!is.null(oracle_lamphi_summary) && nrow(oracle_lamphi_summary) > 0L) {
        utils::write.csv(oracle_lamphi_summary, files$oracle_lambda_phi, row.names = FALSE)
    }

    ladder <- .s4b_build_diagnostic_ladder(
        diag_agg = diag_agg,
        status_counts = status_counts,
        comparison = comp,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary
    )
    utils::write.csv(ladder, files$diagnostic_ladder, row.names = FALSE)

    invisible(list(files = files, status_counts = status_counts, beta_agg = beta_agg,
                   diag_agg = diag_agg, group_agg = group_agg, contrast_agg = contrast_agg,
                   comparison = comp, oracle_lambda_summary = oracle_lambda_summary,
                   oracle_lamphi_summary = oracle_lamphi_summary, diagnostic_ladder = ladder))
}

.s4b_console_summary <- function(diag, saved) {
    cat("\n================ Scenario 4B posterior performance ================\n")
    cat("Main fit: low/heterogeneous exposure, fixed gamma, sampled lambda\n")
    cat("Fit status counts:\n")
    print(saved$status_counts, row.names = FALSE)
    cat("\nData and exposure stress:\n")
    cat(sprintf("  mean count average : %.4f\n", .s4b_safe_mean(diag$mean_count)))
    cat(sprintf("  zero proportion avg: %.4f\n", .s4b_safe_mean(diag$zero_prop)))
    cat(sprintf("  exposure ratio avg : %.5f\n", .s4b_safe_mean(diag$exposure_mean_ratio)))
    cat(sprintf("  exposure CV avg    : %.4f\n", .s4b_safe_mean(diag$exposure_cv)))
    cat("\nPerformance by subset:\n")
    print(saved$diag_agg[, intersect(c(
        "subset", "n_reps", "beta1_mean_avg", "beta2_mean_avg", "beta2_mean_sd",
        "log_lambda_rmse_avg", "cor_log_lambda_avg", "low_exp_log_lambda_rmse_avg",
        "high_exp_log_lambda_rmse_avg", "beta_guard_total", "kappa_guard_total",
        "lambda_output_guard_total"
    ), names(saved$diag_agg)), drop = FALSE], row.names = FALSE)
    cat("\nLow/high exposure contrast by subset:\n")
    print(saved$contrast_agg, row.names = FALSE)
    if (!is.null(saved$diagnostic_ladder) && nrow(saved$diagnostic_ladder) > 0L) {
        cat("\nDiagnostic ladder summary:\n")
        print(saved$diagnostic_ladder[, intersect(c("diagnostic_step", "n_reps", "beta1_mean_avg", "beta2_mean_avg", "beta2_mean_sd", "log_lambda_rmse_avg", "cor_log_lambda_avg", "kappa_guard_total"), names(saved$diagnostic_ladder)), drop = FALSE], row.names = FALSE)
    }
    cat("====================================================================\n\n")
}

summarize_scenario4b_posterior_performance <- function(root = ".",
                                                       summary_file = "output_s4b_low_exposure_continuous_x2/S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100/summary_S4B_low_exposure_fixed_gamma_continuous_x2_all_reps.csv",
                                                       analysis_dir = "analysis_s4b_low_exposure_continuous_x2/S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100",
                                                       scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100",
                                                       s3_control_summary_file = "output_s3_control_fixed_gamma_stabilized/S3_CONTROL_FIXED_GAMMA_STABILIZED_T100/summary_S3_control_fixed_gamma_stabilized_all_reps.csv",
                                                       s4a_performance_subset_file = "output_s4a_sparse_counts_continuous_x2/S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100/tables/posterior_performance_by_subset_continuous_x2.csv",
                                                       s4a_status_counts_file = "output_s4a_sparse_counts_continuous_x2/S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100/tables/s4a_fit_status_counts_continuous_x2.csv",
                                                       oracle_lambda_summary_file = "output_s4b_oracle_lambda_fixed_gamma_continuous_x2/S4B_ORACLE_LAMBDA_FIXED_GAMMA_CONTINUOUS_X2_DIAGNOSTIC_T100/summary_S4B_oracle_lambda_fixed_gamma_continuous_x2_all_reps.csv",
                                                       oracle_lamphi_summary_file = "output_s4b_oracle_lambda_phi_fixed_gamma_continuous_x2/S4B_ORACLE_LAMBDA_PHI_FIXED_GAMMA_CONTINUOUS_X2_DIAGNOSTIC_T100/summary_S4B_oracle_lambda_phi_fixed_gamma_continuous_x2_all_reps.csv",
                                                       make_plots = TRUE,
                                                       plot_format = c("pdf", "png"),
                                                       verbose = TRUE) {
    plot_format <- match.arg(plot_format)
    root <- .s4b_norm_path(root)
    summary_path <- if (grepl("^/", summary_file)) summary_file else file.path(root, summary_file)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)

    tab_dir <- file.path(analysis_dir, "tables")
    fig_dir <- file.path(analysis_dir, "figures")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    raw <- .s4b_read_required_csv(summary_path)
    diag <- .s4b_add_fit_status(raw)
    if (!"scenario_id" %in% names(diag)) diag$scenario_id <- scenario_id

    beta_long <- .s4b_make_beta_long(diag, scenario_id = scenario_id)
    lambda_summary <- .s4b_make_lambda_summary(diag, scenario_id = scenario_id)
    group_long <- .s4b_make_exposure_group_long(diag, scenario_id = scenario_id)
    contrast <- .s4b_low_high_contrast_by_rep(diag)

    s3_control_summary <- .s4b_read_optional_csv(if (grepl("^/", s3_control_summary_file)) s3_control_summary_file else file.path(root, s3_control_summary_file))
    s4a_performance_subset <- .s4b_read_optional_csv(if (grepl("^/", s4a_performance_subset_file)) s4a_performance_subset_file else file.path(root, s4a_performance_subset_file))
    s4a_status_counts <- .s4b_read_optional_csv(if (grepl("^/", s4a_status_counts_file)) s4a_status_counts_file else file.path(root, s4a_status_counts_file))
    oracle_lambda_detail <- .s4b_read_optional_csv(if (grepl("^/", oracle_lambda_summary_file)) oracle_lambda_summary_file else file.path(root, oracle_lambda_summary_file))
    oracle_lamphi_detail <- .s4b_read_optional_csv(if (grepl("^/", oracle_lamphi_summary_file)) oracle_lamphi_summary_file else file.path(root, oracle_lamphi_summary_file))
    oracle_lambda_summary <- .s4b_summarize_oracle_diag(
        oracle_lambda_detail,
        label = "Oracle lambda, estimate beta/phi/r",
        conclusion = "Fixing lambda restores beta2, but beta0/phi/kappa can remain unstable."
    )
    oracle_lamphi_summary <- .s4b_summarize_oracle_diag(
        oracle_lamphi_detail,
        label = "Oracle lambda + oracle phi, estimate beta/r",
        conclusion = "Fixing both lambda and phi restores beta/r and eliminates numerical guards."
    )

    saved <- .s4b_save_tables(
        diag = diag,
        beta_long = beta_long,
        lambda_summary = lambda_summary,
        group_long = group_long,
        contrast = contrast,
        out_dir = tab_dir,
        s3_control_summary = s3_control_summary,
        s4a_performance_subset = s4a_performance_subset,
        s4a_status_counts = s4a_status_counts,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        oracle_lambda_detail = oracle_lambda_detail,
        oracle_lamphi_detail = oracle_lamphi_detail
    )

    plot_files <- list()
    if (isTRUE(make_plots)) {
        plot_files <- .s4b_make_plots(diag = diag, contrast = contrast, fig_dir = fig_dir, plot_format = plot_format)
    }

    rds_file <- file.path(analysis_dir, "scenario4b_posterior_performance_results.rds")
    saveRDS(list(
        diagnostics = diag,
        beta = beta_long,
        lambda = lambda_summary,
        exposure_group = group_long,
        low_high_contrast = contrast,
        aggregates = saved,
        s3_control_summary = s3_control_summary,
        s4a_performance_subset = s4a_performance_subset,
        s4a_status_counts = s4a_status_counts,
        oracle_lambda_detail = oracle_lambda_detail,
        oracle_lamphi_detail = oracle_lamphi_detail,
        oracle_lambda_summary = oracle_lambda_summary,
        oracle_lamphi_summary = oracle_lamphi_summary,
        diagnostic_ladder = saved$diagnostic_ladder
    ), rds_file)

    if (isTRUE(verbose)) {
        .s4b_console_summary(diag, saved)
        cat("Saved Scenario 4B posterior performance RDS: ", rds_file, "\n", sep = "")
        cat("Tables: ", tab_dir, "\n", sep = "")
        if (isTRUE(make_plots)) cat("Figures: ", fig_dir, "\n", sep = "")
    }

    invisible(list(
        diagnostics = diag,
        beta_summary = beta_long,
        lambda_summary = lambda_summary,
        exposure_group_summary = group_long,
        low_high_contrast = contrast,
        tables = saved$files,
        status_counts = saved$status_counts,
        beta_aggregate = saved$beta_agg,
        performance_aggregate = saved$diag_agg,
        exposure_group_aggregate = saved$group_agg,
        contrast_aggregate = saved$contrast_agg,
        comparison = saved$comparison,
        oracle_lambda_summary = saved$oracle_lambda_summary,
        oracle_lamphi_summary = saved$oracle_lamphi_summary,
        diagnostic_ladder = saved$diagnostic_ladder,
        plots = plot_files,
        rds_file = rds_file,
        analysis_dir = analysis_dir,
        tables_dir = tab_dir,
        figures_dir = fig_dir
    ))
}

run_scenario4b_posterior_performance <- function(root = ".",
                                                 scenario_id = "S4B_LOW_HETEROGENEOUS_EXPOSURE_FIXED_GAMMA_CONTINUOUS_X2_T100",
                                                 output_dir = "output_s4b_low_exposure_continuous_x2",
                                                 analysis_dir = "analysis_s4b_low_exposure_continuous_x2",
                                                 make_plots = TRUE,
                                                 plot_format = c("pdf", "png"),
                                                 verbose = TRUE,
                                                 ...) {
    plot_format <- match.arg(plot_format)
    summary_file <- file.path(output_dir, scenario_id, "summary_S4B_low_exposure_fixed_gamma_continuous_x2_all_reps.csv")
    out_dir <- file.path(analysis_dir, scenario_id)
    summarize_scenario4b_posterior_performance(
        root = root,
        summary_file = summary_file,
        analysis_dir = out_dir,
        scenario_id = scenario_id,
        make_plots = make_plots,
        plot_format = plot_format,
        verbose = verbose,
        ...
    )
}

## Auto-run when sourced in a normal interactive workflow.  Set environment
## variable S4B_PERF_AUTO_RUN=false before sourcing if you only want functions.
if (!identical(tolower(Sys.getenv("S4B_PERF_AUTO_RUN", "true")), "false")) {
    run_scenario4b_posterior_performance(root = ".")
}
