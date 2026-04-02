
# ============================================================================
# kappa_diagnostic.R
# Diagnose whether free kappa is absorbing mean structure in the MSSTNB sampler
#
# Usage:
#   source("kappa_diagnostic.R")
#   out <- diagnose_kappa(
#       data_file = "data_revised/S1/data_rep01.rds",
#       fit_file  = "fit_revised_S1_rep01_short.rds",
#       out_prefix = "kappa_diag_S1_rep01"
#   )
#
# Outputs:
#   <out_prefix>_summary.txt
#   <out_prefix>_cell_level.csv
#   <out_prefix>_region_level.csv
#   <out_prefix>_draw_level.csv
#   <out_prefix>_plots.pdf
# ============================================================================

find_component_recursive <- function(x, candidates, path = "fit") {
    if (!is.list(x)) return(NULL)

    nms <- names(x)
    if (!is.null(nms)) {
        for (nm in candidates) {
            if (nm %in% nms) {
                return(list(value = x[[nm]], path = paste0(path, "$", nm)))
            }
        }
        for (nm in nms) {
            out <- find_component_recursive(x[[nm]], candidates, path = paste0(path, "$", nm))
            if (!is.null(out)) return(out)
        }
    } else {
        for (i in seq_along(x)) {
            out <- find_component_recursive(x[[i]], candidates, path = paste0(path, "[[", i, "]]"))
            if (!is.null(out)) return(out)
        }
    }

    NULL
}

extract_kappa_array <- function(fit, TT, n1) {
    res <- find_component_recursive(
        fit,
        candidates = c("kappa", "kappa_draws", "kappa_samples", "draw_kappa")
    )
    if (is.null(res)) {
        stop("Could not find kappa draws inside fit object.")
    }

    x <- res$value
    src_path <- res$path

    # Case 1: list of TT x n1 matrices
    if (is.list(x) && length(x) > 0 && is.matrix(x[[1]]) &&
        nrow(x[[1]]) == TT && ncol(x[[1]]) == n1) {
        S <- length(x)
        arr <- array(NA_real_, dim = c(S, TT, n1))
        for (s in seq_len(S)) arr[s, , ] <- x[[s]]
        return(list(arr = arr, path = src_path))
    }

    # Case 2: array
    if (is.array(x)) {
        d <- dim(x)
        if (length(d) == 3L) {
            # Try to identify TT and n1 positions
            idx_TT <- which(d == TT)
            idx_n1 <- which(d == n1)
            if (length(idx_TT) >= 1L && length(idx_n1) >= 1L) {
                # pick first non-overlapping TT and n1 axes
                used <- integer(0)
                ax_TT <- idx_TT[1]
                used <- c(used, ax_TT)
                ax_n1 <- idx_n1[which(!(idx_n1 %in% used))[1]]
                if (is.na(ax_n1)) stop("Could not identify kappa array axes.")
                ax_S <- setdiff(1:3, c(ax_TT, ax_n1))
                arr <- aperm(x, c(ax_S, ax_TT, ax_n1))
                return(list(arr = arr, path = src_path))
            }
        }
    }

    # Case 3: matrix S x (TT*n1) or (TT*n1) x S
    if (is.matrix(x)) {
        d <- dim(x)
        if (d[2] == TT * n1) {
            S <- d[1]
            arr <- array(as.numeric(t(x)), dim = c(TT, n1, S))
            arr <- aperm(arr, c(3, 1, 2))
            return(list(arr = arr, path = src_path))
        }
        if (d[1] == TT * n1) {
            S <- d[2]
            arr <- array(as.numeric(x), dim = c(TT, n1, S))
            arr <- aperm(arr, c(3, 1, 2))
            return(list(arr = arr, path = src_path))
        }
    }

    stop("Found kappa object at ", src_path, " but could not coerce it to draws array [S, TT, n1].")
}

safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    cor(x[ok], y[ok])
}

diagnose_kappa <- function(data_file = "data_revised/S1/data_rep01.rds",
                           fit_file  = "fit_revised_S1_rep01_short.rds",
                           out_prefix = "kappa_diag") {

    dat <- readRDS(data_file)
    fit <- readRDS(fit_file)

    TT <- dat$TT
    n1 <- dat$n1

    ext <- extract_kappa_array(fit, TT = TT, n1 = n1)
    kappa_arr <- ext$arr
    kappa_src <- ext$path

    S <- dim(kappa_arr)[1]

    # Cell-level summaries
    kappa_mean <- apply(kappa_arr, c(2, 3), mean)
    kappa_sd   <- apply(kappa_arr, c(2, 3), sd)
    logk_mean  <- apply(log(kappa_arr), c(2, 3), mean)
    logk_sd    <- apply(log(kappa_arr), c(2, 3), sd)

    # Draw-level summaries
    draw_mean_kappa <- apply(kappa_arr, 1, mean)
    draw_mean_logk  <- apply(log(kappa_arr), 1, mean)

    # Truth and covariates
    kappa_true <- if ("kappa" %in% names(dat)) dat$kappa else NULL
    x1 <- dat$x1
    x2 <- dat$x2
    xi <- dat$xi
    phi <- dat$phi_star_ident
    if (is.null(phi) && "phi_star" %in% names(dat)) phi <- dat$phi_star

    phi_mat <- matrix(rep(phi, each = TT), nrow = TT, ncol = n1)

    vec_post_mean_kappa <- as.vector(kappa_mean)
    vec_post_mean_logk  <- as.vector(logk_mean)
    vec_x1   <- as.vector(x1)
    vec_x2   <- as.vector(x2)
    vec_phi  <- as.vector(phi_mat)
    vec_logxi <- as.vector(log(xi))

    cor_x1   <- safe_cor(vec_post_mean_logk, vec_x1)
    cor_x2   <- safe_cor(vec_post_mean_logk, vec_x2)
    cor_phi  <- safe_cor(vec_post_mean_logk, vec_phi)
    cor_logxi <- safe_cor(vec_post_mean_logk, vec_logxi)

    lm_post <- lm(vec_post_mean_logk ~ vec_x1 + vec_x2 + vec_phi + vec_logxi)
    r2_post <- summary(lm_post)$r.squared

    truth_cor_logk <- NA_real_
    truth_rmse_logk <- NA_real_
    truth_r2 <- NA_real_

    if (!is.null(kappa_true)) {
        vec_truth_logk <- as.vector(log(kappa_true))
        truth_cor_logk <- safe_cor(vec_post_mean_logk, vec_truth_logk)
        truth_rmse_logk <- sqrt(mean((vec_post_mean_logk - vec_truth_logk)^2))
        lm_truth <- lm(vec_truth_logk ~ vec_x1 + vec_x2 + vec_phi + vec_logxi)
        truth_r2 <- summary(lm_truth)$r.squared
    }

    # Region-level summaries
    region_df <- data.frame(
        region = seq_len(n1),
        post_mean_kappa = colMeans(kappa_mean),
        post_sd_kappa   = apply(kappa_mean, 2, sd),
        post_mean_logk  = colMeans(logk_mean),
        post_sd_logk    = apply(logk_mean, 2, sd),
        x1_mean         = colMeans(x1),
        x2_mean         = colMeans(x2),
        logxi_mean      = colMeans(log(xi)),
        phi_ident       = phi
    )
    if (!is.null(kappa_true)) {
        region_df$truth_mean_kappa <- colMeans(kappa_true)
        region_df$truth_mean_logk  <- colMeans(log(kappa_true))
    }

    cell_df <- data.frame(
        time = rep(seq_len(TT), times = n1),
        region = rep(seq_len(n1), each = TT),
        post_mean_kappa = as.vector(kappa_mean),
        post_sd_kappa   = as.vector(kappa_sd),
        post_mean_logk  = as.vector(logk_mean),
        post_sd_logk    = as.vector(logk_sd),
        x1 = vec_x1,
        x2 = vec_x2,
        phi_ident = vec_phi,
        xi = as.vector(xi),
        logxi = vec_logxi
    )
    if (!is.null(kappa_true)) {
        cell_df$truth_kappa <- as.vector(kappa_true)
        cell_df$truth_logk  <- as.vector(log(kappa_true))
    }

    draw_df <- data.frame(
        draw = seq_len(S),
        mean_kappa = draw_mean_kappa,
        mean_logk  = draw_mean_logk
    )

    # Summary text
    lines <- c(
        "============================================================",
        "  Kappa Diagnostic Report",
        "============================================================",
        paste("Data file :", data_file),
        paste("Fit file  :", fit_file),
        paste("Kappa source in fit :", kappa_src),
        "",
        "--- Key interpretation note ---",
        "For kappa ~ Gamma(r, r), E(kappa)=1, but E(log kappa) is generally NOT 0.",
        "So inspect mean(kappa) relative to 1, and inspect whether posterior mean log(kappa)",
        "is correlated with x1, x2, phi, or log(xi), which would suggest kappa is",
        "absorbing mean structure.",
        "",
        "--- Draw-level summaries ---",
        sprintf("Number of saved draws          : %d", S),
        sprintf("Posterior mean of mean(kappa)  : %.4f", mean(draw_mean_kappa)),
        sprintf("Posterior SD of mean(kappa)    : %.4f", sd(draw_mean_kappa)),
        sprintf("Posterior mean of mean(log k)) : %.4f", mean(draw_mean_logk)),
        sprintf("Posterior SD of mean(log k))   : %.4f", sd(draw_mean_logk)),
        "",
        "--- Cell-level structure diagnostics (posterior mean log kappa) ---",
        sprintf("cor(post_mean_logk, x1)        : %.3f", cor_x1),
        sprintf("cor(post_mean_logk, x2)        : %.3f", cor_x2),
        sprintf("cor(post_mean_logk, phi_ident) : %.3f", cor_phi),
        sprintf("cor(post_mean_logk, log xi)    : %.3f", cor_logxi),
        sprintf("R^2[post_mean_logk ~ x1+x2+phi+logxi] : %.3f", r2_post)
    )

    if (!is.null(kappa_true)) {
        lines <- c(
            lines,
            "",
            "--- Truth comparison ---",
            sprintf("cor(post_mean_logk, truth_logk): %.3f", truth_cor_logk),
            sprintf("RMSE(post_mean_logk, truth_logk): %.3f", truth_rmse_logk),
            sprintf("R^2[truth_logk ~ x1+x2+phi+logxi]    : %.3f", truth_r2)
        )
    }

    lines <- c(
        lines,
        "",
        "--- Heuristic reading ---",
        "Large |cor(post_mean_logk, x1/x2/phi/logxi)| or a large R^2 suggests kappa",
        "is absorbing structured signal instead of behaving like mean-one observation noise.",
        "If this is much stronger than in truth, that strongly implicates the kappa block."
    )

    writeLines(lines, con = paste0(out_prefix, "_summary.txt"))
    write.csv(cell_df,   paste0(out_prefix, "_cell_level.csv"),   row.names = FALSE)
    write.csv(region_df, paste0(out_prefix, "_region_level.csv"), row.names = FALSE)
    write.csv(draw_df,   paste0(out_prefix, "_draw_level.csv"),   row.names = FALSE)

    # Plots
    pdf(paste0(out_prefix, "_plots.pdf"), width = 11, height = 8.5)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

    hist(draw_mean_kappa, breaks = 30, main = "Draw-level mean(kappa)",
         xlab = "mean kappa across all cells")
    abline(v = 1, lty = 2)

    plot(vec_logxi, vec_post_mean_logk, pch = 16, cex = 0.5,
         xlab = "log(xi)", ylab = "posterior mean log(kappa)",
         main = sprintf("post mean log(kappa) vs log(xi)\ncor = %.3f", cor_logxi))
    abline(lm(vec_post_mean_logk ~ vec_logxi), lty = 2)

    plot(vec_x1, vec_post_mean_logk, pch = 16, cex = 0.5,
         xlab = "x1", ylab = "posterior mean log(kappa)",
         main = sprintf("post mean log(kappa) vs x1\ncor = %.3f", cor_x1))
    abline(lm(vec_post_mean_logk ~ vec_x1), lty = 2)

    plot(vec_x2, vec_post_mean_logk, pch = 16, cex = 0.5,
         xlab = "x2", ylab = "posterior mean log(kappa)",
         main = sprintf("post mean log(kappa) vs x2\ncor = %.3f", cor_x2))
    abline(lm(vec_post_mean_logk ~ vec_x2), lty = 2)

    plot(vec_phi, vec_post_mean_logk, pch = 16, cex = 0.5,
         xlab = "phi_ident", ylab = "posterior mean log(kappa)",
         main = sprintf("post mean log(kappa) vs phi\ncor = %.3f", cor_phi))
    abline(lm(vec_post_mean_logk ~ vec_phi), lty = 2)

    if (!is.null(kappa_true)) {
        plot(as.vector(log(kappa_true)), vec_post_mean_logk, pch = 16, cex = 0.5,
             xlab = "truth log(kappa)", ylab = "posterior mean log(kappa)",
             main = sprintf("truth vs posterior mean log(kappa)\ncor = %.3f", truth_cor_logk))
        abline(0, 1, lty = 2)
    } else {
        plot(region_df$region, region_df$post_mean_kappa, type = "b",
             xlab = "region", ylab = "time-avg posterior mean kappa",
             main = "Region-wise posterior mean kappa")
        abline(h = 1, lty = 2)
    }

    par(op)
    dev.off()

    cat(paste(lines, collapse = "\n"), "\n")
    invisible(list(
        summary_lines = lines,
        cell_df = cell_df,
        region_df = region_df,
        draw_df = draw_df
    ))
}
