## ==========================================================================
## update_gamma.R  -- AUDIT VERSION v2
##
## Fixes:
##   - robust handling of audit iteration matching
##   - internal call counter so sampler does NOT need to pass iter index
## ==========================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.safe_scalar <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0L) return(default)
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(default)
    x[1L]
}

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))

log_marginal_gamma <- function(gamma_j, y_j, zeta_j, a0, b0, return_parts = FALSE) {
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

gamma_internal_log_target <- function(gamma_j, y_j, zeta_j, a0, b0, priors,
                                      return_parts = FALSE) {
    ml <- log_marginal_gamma(gamma_j, y_j, zeta_j, a0, b0, return_parts = return_parts)

    if (!return_parts) {
        if (!is.finite(ml)) return(-Inf)
        log_prior <- dbeta(gamma_j, priors$gamma_a, priors$gamma_b, log = TRUE)
        log_jac   <- log(gamma_j) + log(1 - gamma_j)
        return(ml + log_prior + log_jac)
    }

    log_ml <- ml$log_ml
    log_prior <- if (is.finite(log_ml)) dbeta(gamma_j, priors$gamma_a, priors$gamma_b, log = TRUE) else -Inf
    log_jac   <- if (is.finite(log_ml)) log(gamma_j) + log(1 - gamma_j) else -Inf
    total <- log_ml + log_prior + log_jac

    c(ml, list(
        log_prior = log_prior,
        log_jac = log_jac,
        log_target = total
    ))
}

build_gamma_internal_curve <- function(grid, y_j, zeta_j, a0, b0, priors) {
    out <- lapply(grid, function(gm) {
        parts <- gamma_internal_log_target(
            gamma_j = gm, y_j = y_j, zeta_j = zeta_j,
            a0 = a0, b0 = b0, priors = priors,
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

get_gamma_audit_spec <- function() {
    spec <- getOption("MSSTNB_GAMMA_AUDIT", NULL)
    if (is.null(spec) || !is.list(spec)) return(NULL)
    if (!isTRUE(spec$enabled)) return(NULL)
    spec
}

store_gamma_audit_record <- function(record, spec) {
    env <- spec$env %||% .GlobalEnv
    if (is.null(env$.gamma_audit_records)) {
        env$.gamma_audit_records <- list()
    }
    env$.gamma_audit_records[[length(env$.gamma_audit_records) + 1L]] <- record
    invisible(NULL)
}

.next_gamma_audit_iter <- function() {
    ctr <- getOption("MSSTNB_GAMMA_AUDIT_COUNTER", 0L)
    ctr <- as.integer(ctr) + 1L
    options(MSSTNB_GAMMA_AUDIT_COUNTER = ctr)
    ctr
}

.should_audit_iter <- function(current_iter, audit_iter) {
    if (is.null(audit_iter)) return(TRUE)
    if (length(audit_iter) == 0L) return(TRUE)
    if (length(current_iter) == 0L || is.na(current_iter)) return(FALSE)
    isTRUE(any(current_iter == audit_iter))
}

.should_audit_region <- function(j, audit_regions) {
    if (is.null(audit_regions)) return(TRUE)
    if (length(audit_regions) == 0L) return(TRUE)
    isTRUE(any(j == audit_regions))
}

update_gamma <- function(gamma_current, y_coarse, xi, kappa,
                         a0, b0, priors, mh_sd, ...) {
    n1 <- ncol(y_coarse)
    gamma_new <- gamma_current
    accept <- logical(n1)

    audit_spec <- get_gamma_audit_spec()
    current_iter <- .next_gamma_audit_iter()
    audit_iter <- audit_spec$iter %||% NULL
    audit_regions <- audit_spec$region %||% NULL
    audit_grid <- audit_spec$grid %||% seq(0.01, 0.999, length.out = 200L)
    label <- audit_spec$label %||% "gamma_audit"

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]

        logit_current  <- logit(gamma_current[j])
        logit_proposal <- logit_current + rnorm(
            1, mean = 0,
            sd = if (length(mh_sd) > 1L) mh_sd[j] else mh_sd
        )
        gamma_proposal <- expit(logit_proposal)

        parts_current <- gamma_internal_log_target(
            gamma_j = gamma_current[j], y_j = y_coarse[, j], zeta_j = zeta_j,
            a0 = a0, b0 = b0, priors = priors, return_parts = TRUE
        )
        parts_proposal <- gamma_internal_log_target(
            gamma_j = gamma_proposal, y_j = y_coarse[, j], zeta_j = zeta_j,
            a0 = a0, b0 = b0, priors = priors, return_parts = TRUE
        )

        log_alpha <- parts_proposal$log_target - parts_current$log_target

        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
            gamma_new[j] <- gamma_proposal
            accept[j] <- TRUE
        }

        do_audit_iter <- .should_audit_iter(current_iter, audit_iter)
        do_audit_region <- .should_audit_region(j, audit_regions)

        if (!is.null(audit_spec) && do_audit_iter && do_audit_region) {
            curve_df <- build_gamma_internal_curve(
                grid = audit_grid,
                y_j = y_coarse[, j],
                zeta_j = zeta_j,
                a0 = a0, b0 = b0,
                priors = priors
            )

            lt <- curve_df$log_target
            finite_idx <- is.finite(lt)
            if (any(finite_idx)) {
                mlt <- max(lt[finite_idx])
                w <- rep(0, length(lt))
                w[finite_idx] <- exp(lt[finite_idx] - mlt)
                w <- w / sum(w)
                curve_mode <- curve_df$gamma[which.max(w)]
                curve_mean <- sum(curve_df$gamma * w)
            } else {
                curve_mode <- NA_real_
                curve_mean <- NA_real_
            }

            rec <- list(
                label = label,
                iter = current_iter,
                region = j,
                gamma_current = gamma_current[j],
                gamma_proposal = gamma_proposal,
                log_target_current = parts_current$log_target,
                log_target_proposal = parts_proposal$log_target,
                log_ml_current = parts_current$log_ml,
                log_ml_proposal = parts_proposal$log_ml,
                log_prior_current = parts_current$log_prior,
                log_prior_proposal = parts_proposal$log_prior,
                log_jac_current = parts_current$log_jac,
                log_jac_proposal = parts_proposal$log_jac,
                log_alpha = log_alpha,
                accepted = accept[j],
                zeta = zeta_j,
                y = y_coarse[, j],
                a_prev_current = parts_current$a_prev,
                b_prev_current = parts_current$b_prev,
                alpha_current = parts_current$alpha,
                beta_current = parts_current$beta,
                log_pred_current = parts_current$log_pred,
                a_prev_proposal = parts_proposal$a_prev,
                b_prev_proposal = parts_proposal$b_prev,
                alpha_proposal = parts_proposal$alpha,
                beta_proposal = parts_proposal$beta,
                log_pred_proposal = parts_proposal$log_pred,
                internal_curve = curve_df,
                internal_curve_mode = curve_mode,
                internal_curve_mean = curve_mean,
                mh_sd_used = if (length(mh_sd) > 1L) mh_sd[j] else mh_sd
            )
            store_gamma_audit_record(rec, audit_spec)
        }
    }

    list(gamma = gamma_new, accept = accept)
}
