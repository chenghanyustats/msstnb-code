# ============================================================
# External marginal predictive gamma target vs internal audit curve
#
# Purpose:
#   Given one gamma audit record saved from update_gamma_audit_version_v2.R,
#   independently re-compute the same marginal predictive gamma target
#   outside the sampler, then compare:
#       - internal curve used by sampler
#       - external independently re-computed curve
#
# Main function:
#   compare_gamma_internal_external(...)
#
# Typical use:
#
# source("compare_gamma_internal_external.R")
#
# out <- compare_gamma_internal_external(
#     audit_rds   = "gamma_audit_fit_S1_rep01_audit_record.rds",
#     record_index = 1L,
#     project_dir = "/Users/chenghanyu/Dropbox/academia/research/msstnb_new/simulation",
#     out_prefix  = "gamma_internal_external_rep01"
# )
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.safe_scalar <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(default)
    x[1L]
}

discrete_weighted_quantile <- function(x, w, probs) {
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1L))
}

summarize_curve_from_logtarget <- function(grid, log_target) {
    finite_idx <- is.finite(log_target)
    if (!any(finite_idx)) {
        return(list(mode = NA_real_, mean = NA_real_, q025 = NA_real_, q975 = NA_real_))
    }
    mlt <- max(log_target[finite_idx])
    w <- rep(0, length(log_target))
    w[finite_idx] <- exp(log_target[finite_idx] - mlt)
    w <- w / sum(w)
    qs <- discrete_weighted_quantile(grid, w, c(0.025, 0.975))
    list(
        mode = grid[which.max(w)],
        mean = sum(grid * w),
        q025 = qs[1L],
        q975 = qs[2L],
        weights = w
    )
}

find_file <- function(filename, roots) {
    for (root in roots) {
        cand <- file.path(root, filename)
        if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
    stop("Could not find required file: ", filename,
         "\nSearched in:\n", paste(roots, collapse = "\n"))
}

external_log_marginal_gamma <- function(gamma_j, y_j, zeta_j, a0, b0, return_parts = FALSE) {
    TT <- length(y_j)

    if (!is.finite(gamma_j) || gamma_j <= 0 || gamma_j >= 1) {
        if (!return_parts) return(-Inf)
        return(list(
            log_ml = -Inf,
            alpha = rep(NA_real_, TT),
            beta = rep(NA_real_, TT),
            log_pred = rep(NA_real_, TT),
            a_prev = rep(NA_real_, TT),
            b_prev = rep(NA_real_, TT)
        ))
    }

    log_ml <- 0
    at <- a0
    bt <- b0

    if (return_parts) {
        alpha_vec <- beta_vec <- log_pred_vec <- rep(NA_real_, TT)
        a_prev_vec <- b_prev_vec <- rep(NA_real_, TT)
    }

    for (t in seq_len(TT)) {
        alpha <- gamma_j * at
        beta  <- gamma_j * bt
        y_t   <- y_j[t]
        z_t   <- zeta_j[t]

        if (!is.finite(z_t) || z_t <= 0 || !is.finite(alpha) || !is.finite(beta) ||
            alpha <= 0 || beta <= 0) {
            if (!return_parts) return(-Inf)
            return(list(
                log_ml = -Inf,
                alpha = alpha_vec,
                beta = beta_vec,
                log_pred = log_pred_vec,
                a_prev = a_prev_vec,
                b_prev = b_prev_vec
            ))
        }

        log_pred <- lgamma(alpha + y_t) - lgamma(alpha) - lgamma(y_t + 1) +
            alpha * log(beta) + y_t * log(z_t) -
            (alpha + y_t) * log(beta + z_t)

        log_ml <- log_ml + log_pred

        if (return_parts) {
            a_prev_vec[t] <- at
            b_prev_vec[t] <- bt
            alpha_vec[t] <- alpha
            beta_vec[t] <- beta
            log_pred_vec[t] <- log_pred
        }

        at <- alpha + y_t
        bt <- beta  + z_t
    }

    if (!return_parts) return(log_ml)

    list(
        log_ml = log_ml,
        alpha = alpha_vec,
        beta = beta_vec,
        log_pred = log_pred_vec,
        a_prev = a_prev_vec,
        b_prev = b_prev_vec
    )
}

external_gamma_log_target <- function(gamma_j, y_j, zeta_j, a0, b0, prior_a, prior_b,
                                      return_parts = FALSE) {
    ml <- external_log_marginal_gamma(
        gamma_j = gamma_j, y_j = y_j, zeta_j = zeta_j,
        a0 = a0, b0 = b0, return_parts = return_parts
    )

    if (!return_parts) {
        if (!is.finite(ml)) return(-Inf)
        log_prior <- dbeta(gamma_j, prior_a, prior_b, log = TRUE)
        log_jac   <- log(gamma_j) + log(1 - gamma_j)
        return(ml + log_prior + log_jac)
    }

    log_ml <- ml$log_ml
    log_prior <- if (is.finite(log_ml)) dbeta(gamma_j, prior_a, prior_b, log = TRUE) else -Inf
    log_jac   <- if (is.finite(log_ml)) log(gamma_j) + log(1 - gamma_j) else -Inf
    total <- log_ml + log_prior + log_jac

    c(ml, list(
        log_prior = log_prior,
        log_jac = log_jac,
        log_target = total
    ))
}

build_external_curve <- function(grid, y_j, zeta_j, a0, b0, prior_a, prior_b) {
    out <- lapply(grid, function(gm) {
        parts <- external_gamma_log_target(
            gamma_j = gm, y_j = y_j, zeta_j = zeta_j,
            a0 = a0, b0 = b0, prior_a = prior_a, prior_b = prior_b,
            return_parts = TRUE
        )
        data.frame(
            gamma = gm,
            log_target = .safe_scalar(parts$log_target),
            log_ml = .safe_scalar(parts$log_ml),
            log_prior = .safe_scalar(parts$log_prior),
            log_jac = .safe_scalar(parts$log_jac)
        )
    })
    do.call(rbind, out)
}

compare_gamma_internal_external <- function(
    audit_rds,
    record_index = 1L,
    project_dir = NULL,
    prior_a = NULL,
    prior_b = NULL,
    out_prefix = "gamma_internal_external_compare",
    save_png = TRUE,
    save_pdf = TRUE,
    verbose = TRUE
) {
    if (!file.exists(audit_rds)) {
        stop("Audit record file not found: ", audit_rds)
    }

    recs <- readRDS(audit_rds)
    if (!is.list(recs) || length(recs) < record_index) {
        stop("Audit record file does not contain record_index = ", record_index)
    }
    rec <- recs[[record_index]]

    a0 <- .safe_scalar(rec$a_prev_current[1L])
    b0 <- .safe_scalar(rec$b_prev_current[1L])

    if (!is.finite(a0) || !is.finite(b0)) {
        stop("Could not infer a0/b0 from audit record.")
    }

    if (is.null(prior_a) || is.null(prior_b)) {
        if (is.null(project_dir)) {
            stop("You must supply either (prior_a, prior_b) or project_dir.")
        }
        roots <- c(
            project_dir,
            file.path(project_dir, "R"),
            file.path(project_dir, "R", "mcmc"),
            file.path(project_dir, "mcmc"),
            file.path(project_dir, "scripts"),
            file.path(project_dir, "code"),
            file.path(project_dir, "R_gpt"),
            file.path(project_dir, "R_gpt", "mcmc")
        )
        f_cfg <- find_file("mcmc_config.R", roots)
        source(f_cfg, local = FALSE)
        prior_a <- prior_a %||% .safe_scalar(MCMC_PRIORS$gamma_a)
        prior_b <- prior_b %||% .safe_scalar(MCMC_PRIORS$gamma_b)
    }

    if (!is.finite(prior_a) || !is.finite(prior_b)) {
        stop("Could not determine gamma prior parameters.")
    }

    grid <- rec$internal_curve$gamma
    y_j <- as.numeric(rec$y)
    zeta_j <- as.numeric(rec$zeta)

    external_curve <- build_external_curve(
        grid = grid, y_j = y_j, zeta_j = zeta_j,
        a0 = a0, b0 = b0, prior_a = prior_a, prior_b = prior_b
    )

    internal_curve <- rec$internal_curve

    diff_raw <- internal_curve$log_target - external_curve$log_target
    finite_idx <- is.finite(diff_raw)

    max_abs_diff <- if (any(finite_idx)) max(abs(diff_raw[finite_idx])) else NA_real_
    mean_diff <- if (any(finite_idx)) mean(diff_raw[finite_idx]) else NA_real_
    diff_centered <- diff_raw - mean_diff
    max_abs_diff_centered <- if (any(finite_idx)) max(abs(diff_centered[finite_idx])) else NA_real_
    cor_logtarget <- if (sum(finite_idx) >= 2L) suppressWarnings(cor(internal_curve$log_target[finite_idx],
                                                                     external_curve$log_target[finite_idx])) else NA_real_

    ext_current <- external_gamma_log_target(
        gamma_j = rec$gamma_current, y_j = y_j, zeta_j = zeta_j,
        a0 = a0, b0 = b0, prior_a = prior_a, prior_b = prior_b,
        return_parts = TRUE
    )
    ext_proposal <- external_gamma_log_target(
        gamma_j = rec$gamma_proposal, y_j = y_j, zeta_j = zeta_j,
        a0 = a0, b0 = b0, prior_a = prior_a, prior_b = prior_b,
        return_parts = TRUE
    )

    current_diff <- rec$log_target_current - ext_current$log_target
    proposal_diff <- rec$log_target_proposal - ext_proposal$log_target
    log_alpha_internal <- rec$log_target_proposal - rec$log_target_current
    log_alpha_external <- ext_proposal$log_target - ext_current$log_target
    log_alpha_diff <- log_alpha_internal - log_alpha_external

    int_sum <- summarize_curve_from_logtarget(internal_curve$gamma, internal_curve$log_target)
    ext_sum <- summarize_curve_from_logtarget(external_curve$gamma, external_curve$log_target)

    curve_df <- data.frame(
        gamma = grid,
        internal_log_target = internal_curve$log_target,
        internal_log_ml = internal_curve$log_ml,
        internal_log_prior = internal_curve$log_prior,
        internal_log_jac = internal_curve$log_jac,
        external_log_target = external_curve$log_target,
        external_log_ml = external_curve$log_ml,
        external_log_prior = external_curve$log_prior,
        external_log_jac = external_curve$log_jac,
        diff_raw = diff_raw,
        diff_centered = diff_centered
    )
    utils::write.csv(curve_df, paste0(out_prefix, "_curve.csv"), row.names = FALSE)

    lines <- c(
        "============================================================",
        "Gamma internal vs external marginal predictive target",
        "============================================================",
        paste("Audit record file        :", audit_rds),
        paste("Record index             :", record_index),
        paste("Recorded iter            :", rec$iter),
        paste("Recorded region          :", rec$region),
        paste("a0, b0                   :", round(a0, 6), ", ", round(b0, 6), sep = ""),
        paste("prior_a, prior_b         :", round(prior_a, 6), ", ", round(prior_b, 6), sep = ""),
        "",
        "Internal curve summary",
        paste("  mode                   :", round(int_sum$mode, 6)),
        paste("  mean                   :", round(int_sum$mean, 6)),
        paste("  95% CI                 : [", round(int_sum$q025, 6), ", ", round(int_sum$q975, 6), "]", sep = ""),
        "",
        "External curve summary",
        paste("  mode                   :", round(ext_sum$mode, 6)),
        paste("  mean                   :", round(ext_sum$mean, 6)),
        paste("  95% CI                 : [", round(ext_sum$q025, 6), ", ", round(ext_sum$q975, 6), "]", sep = ""),
        "",
        "Pointwise curve comparison",
        paste("  cor(log_target)        :", round(cor_logtarget, 6)),
        paste("  max abs diff           :", round(max_abs_diff, 10)),
        paste("  mean diff              :", round(mean_diff, 10)),
        paste("  max abs diff centered  :", round(max_abs_diff_centered, 10)),
        "",
        "Current/proposal direct check",
        paste("  gamma_current          :", round(rec$gamma_current, 6)),
        paste("  gamma_proposal         :", round(rec$gamma_proposal, 6)),
        paste("  internal current       :", round(rec$log_target_current, 10)),
        paste("  external current       :", round(ext_current$log_target, 10)),
        paste("  diff current           :", round(current_diff, 10)),
        paste("  internal proposal      :", round(rec$log_target_proposal, 10)),
        paste("  external proposal      :", round(ext_proposal$log_target, 10)),
        paste("  diff proposal          :", round(proposal_diff, 10)),
        paste("  internal log alpha     :", round(log_alpha_internal, 10)),
        paste("  external log alpha     :", round(log_alpha_external, 10)),
        paste("  diff log alpha         :", round(log_alpha_diff, 10)),
        paste("  accepted               :", rec$accepted)
    )
    writeLines(lines, con = paste0(out_prefix, "_summary.txt"))

    if (save_png) {
        png(paste0(out_prefix, ".png"), width = 1400, height = 1100, res = 140)
        par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
        plot(grid, internal_curve$log_target, type = "l", lwd = 2,
             main = "Internal vs external log target",
             xlab = expression(gamma), ylab = "log target")
        lines(grid, external_curve$log_target, lwd = 2, lty = 2)
        abline(v = rec$gamma_current, lty = 3)
        abline(v = rec$gamma_proposal, lty = 3)
        legend("topright", legend = c("internal", "external"),
               lty = c(1, 2), lwd = 2, bty = "n")

        plot(grid, diff_raw, type = "l", lwd = 2,
             main = "Raw difference: internal - external",
             xlab = expression(gamma), ylab = "difference")
        abline(h = 0, lty = 2)

        plot(grid, diff_centered, type = "l", lwd = 2,
             main = "Centered difference",
             xlab = expression(gamma), ylab = "difference after removing mean shift")
        abline(h = 0, lty = 2)

        plot(grid, int_sum$weights, type = "l", lwd = 2,
             main = "Normalized posterior mass on grid",
             xlab = expression(gamma), ylab = "mass")
        lines(grid, ext_sum$weights, lwd = 2, lty = 2)
        abline(v = rec$gamma_current, lty = 3)
        abline(v = rec$gamma_proposal, lty = 3)
        legend("topright", legend = c("internal", "external"),
               lty = c(1, 2), lwd = 2, bty = "n")
        dev.off()
    }

    if (save_pdf) {
        pdf(paste0(out_prefix, ".pdf"), width = 10, height = 8)
        par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
        plot(grid, internal_curve$log_target, type = "l", lwd = 2,
             main = "Internal vs external log target",
             xlab = expression(gamma), ylab = "log target")
        lines(grid, external_curve$log_target, lwd = 2, lty = 2)
        abline(v = rec$gamma_current, lty = 3)
        abline(v = rec$gamma_proposal, lty = 3)
        legend("topright", legend = c("internal", "external"),
               lty = c(1, 2), lwd = 2, bty = "n")

        plot(grid, diff_raw, type = "l", lwd = 2,
             main = "Raw difference: internal - external",
             xlab = expression(gamma), ylab = "difference")
        abline(h = 0, lty = 2)

        plot(grid, diff_centered, type = "l", lwd = 2,
             main = "Centered difference",
             xlab = expression(gamma), ylab = "difference after removing mean shift")
        abline(h = 0, lty = 2)

        plot(grid, int_sum$weights, type = "l", lwd = 2,
             main = "Normalized posterior mass on grid",
             xlab = expression(gamma), ylab = "mass")
        lines(grid, ext_sum$weights, lwd = 2, lty = 2)
        abline(v = rec$gamma_current, lty = 3)
        abline(v = rec$gamma_proposal, lty = 3)
        legend("topright", legend = c("internal", "external"),
               lty = c(1, 2), lwd = 2, bty = "n")
        dev.off()
    }

    if (verbose) {
        cat(paste(lines, collapse = "\n"), "\n")
        cat("\nWrote:\n")
        cat("  ", paste0(out_prefix, "_curve.csv"), "\n", sep = "")
        cat("  ", paste0(out_prefix, "_summary.txt"), "\n", sep = "")
        cat("  ", paste0(out_prefix, ".png"), "\n", sep = "")
        cat("  ", paste0(out_prefix, ".pdf"), "\n", sep = "")
    }

    invisible(list(
        record = rec,
        internal_curve = internal_curve,
        external_curve = external_curve,
        summary_lines = lines,
        curve_df = curve_df,
        metrics = list(
            cor_logtarget = cor_logtarget,
            max_abs_diff = max_abs_diff,
            mean_diff = mean_diff,
            max_abs_diff_centered = max_abs_diff_centered,
            current_diff = current_diff,
            proposal_diff = proposal_diff,
            log_alpha_diff = log_alpha_diff,
            internal_mode = int_sum$mode,
            internal_mean = int_sum$mean,
            external_mode = ext_sum$mode,
            external_mean = ext_sum$mean
        )
    ))
}

compare_gamma_internal_external_current <- function(
    audit_rds = "gamma_audit_fit_S1_rep01_audit_record.rds",
    record_index = 1L,
    project_dir = "/Users/chenghanyu/Dropbox/academia/research/msstnb_new/simulation",
    out_prefix = "gamma_internal_external_rep01",
    save_png = TRUE,
    save_pdf = TRUE,
    verbose = TRUE
) {
    compare_gamma_internal_external(
        audit_rds = audit_rds,
        record_index = record_index,
        project_dir = project_dir,
        out_prefix = out_prefix,
        save_png = save_png,
        save_pdf = save_pdf,
        verbose = verbose
    )
}
