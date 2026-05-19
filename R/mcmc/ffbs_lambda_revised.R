## ==========================================================================
## ffbs_lambda_revised.R
## Robust Forward Filter Backward Sampler for the coarsest residual risk
## lambda_tilde, with deterministic re-centering to an identified space.
##
## Model block:
##     y_{t,j} | lambda_tilde_{t,j}, zeta_{t,j}
##         ~ Poisson(zeta_{t,j} * lambda_tilde_{t,j}),
##     zeta_{t,j} = xi_{t,j} * kappa_{t,j}.
##
## Forward filtering recursion:
##     a_{t,j} = gamma_j * a_{t-1,j} + y_{t,j},
##     b_{t,j} = gamma_j * b_{t-1,j} + zeta_{t,j}.
##
## Backward kernel:
##     lambda_tilde_{T,j} ~ Gamma(a_{T,j}, b_{T,j}),
##     lambda_tilde_{t-1,j} = U_{t-1,j} + gamma_j * lambda_tilde_{t,j},
##     U_{t-1,j} ~ Gamma((1 - gamma_j) * a_{t-1,j}, b_{t-1,j}).
## ===========================================================================

.lambda_floor <- function() .Machine$double.xmin

.is_count_vector <- function(x) {
    is.numeric(x) && all(is.finite(x)) && all(x >= 0) &&
        all(abs(x - round(x)) < sqrt(.Machine$double.eps))
}

.check_positive_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(name, " must be a positive finite scalar.", call. = FALSE)
    }
    invisible(TRUE)
}

.check_same_dim <- function(x, y, x_name, y_name) {
    if (is.null(dim(x)) || is.null(dim(y)) || !identical(dim(x), dim(y))) {
        stop(x_name, " and ", y_name, " must have identical dimensions.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Forward filter for one region
#'
#' @param gamma_j Discount factor for region j. Must lie in (0, 1).
#' @param y_j Length-T nonnegative integer count vector.
#' @param zeta_j Length-T positive effective offset vector.
#' @param a0,b0 Positive initial Gamma shape and rate.
#' @param min_shape Lower numerical bound for Gamma shapes.
#' @param min_rate Lower numerical bound for Gamma rates.
#' @return A list with filtered shapes and rates.
lambda_forward_filter_j <- function(gamma_j, y_j, zeta_j, a0, b0,
                                    min_shape = 1e-12,
                                    min_rate = 1e-12) {
    if (!is.numeric(gamma_j) || length(gamma_j) != 1L ||
        !is.finite(gamma_j) || gamma_j <= 0 || gamma_j >= 1) {
        stop("gamma_j must be a finite scalar strictly between 0 and 1.",
             call. = FALSE)
    }
    if (!.is_count_vector(y_j)) {
        stop("y_j must be a finite nonnegative integer-valued numeric vector.",
             call. = FALSE)
    }
    if (!is.numeric(zeta_j) || length(zeta_j) != length(y_j) ||
        any(!is.finite(zeta_j)) || any(zeta_j <= 0)) {
        stop("zeta_j must be a positive finite vector with the same length as y_j.",
             call. = FALSE)
    }
    .check_positive_scalar(a0, "a0")
    .check_positive_scalar(b0, "b0")
    .check_positive_scalar(min_shape, "min_shape")
    .check_positive_scalar(min_rate, "min_rate")

    TT <- length(y_j)
    if (TT < 1L) {
        stop("y_j must have length at least 1.", call. = FALSE)
    }

    a_filt <- numeric(TT)
    b_filt <- numeric(TT)

    a_prev <- a0
    b_prev <- b0

    for (t in seq_len(TT)) {
        a_pred <- max(gamma_j * a_prev, min_shape)
        b_pred <- max(gamma_j * b_prev, min_rate)

        a_filt[t] <- max(a_pred + y_j[t], min_shape)
        b_filt[t] <- max(b_pred + zeta_j[t], min_rate)

        a_prev <- a_filt[t]
        b_prev <- b_filt[t]
    }

    list(a_filt = a_filt, b_filt = b_filt)
}

#' Forward Filter Backward Sampler for lambda_tilde_{1:T,j}
#'
#' @param gamma_j Discount factor for region j.
#' @param y_j Length-T count vector for region j.
#' @param zeta_j Length-T effective offset vector, zeta = xi * kappa.
#' @param a0,b0 Initial Gamma filter shape and rate.
#' @param lambda_floor Positive lower bound used to avoid exact zero draws.
#' @param return_filter If TRUE, return the sampled path and filtering arrays.
#' @return Either a length-T vector, or a list with lambda_path, a_filt, b_filt.
ffbs_lambda_j <- function(gamma_j, y_j, zeta_j, a0, b0,
                          lambda_floor = .lambda_floor(),
                          min_shape = 1e-12,
                          min_rate = 1e-12,
                          return_filter = FALSE) {
    .check_positive_scalar(lambda_floor, "lambda_floor")

    filt <- lambda_forward_filter_j(
        gamma_j = gamma_j,
        y_j = y_j,
        zeta_j = zeta_j,
        a0 = a0,
        b0 = b0,
        min_shape = min_shape,
        min_rate = min_rate
    )

    a_filt <- filt$a_filt
    b_filt <- filt$b_filt
    TT <- length(y_j)

    lambda_path <- numeric(TT)

    terminal <- rgamma(1L, shape = a_filt[TT], rate = b_filt[TT])
    lambda_path[TT] <- max(terminal, lambda_floor)

    if (TT >= 2L) {
        for (t in TT:2L) {
            shape_U <- max((1 - gamma_j) * a_filt[t - 1L], min_shape)
            rate_U  <- max(b_filt[t - 1L], min_rate)
            U <- rgamma(1L, shape = shape_U, rate = rate_U)

            lambda_path[t - 1L] <- max(U + gamma_j * lambda_path[t], lambda_floor)
        }
    }

    if (!return_filter) {
        return(lambda_path)
    }

    list(
        lambda_path = lambda_path,
        a_filt = a_filt,
        b_filt = b_filt
    )
}

#' Run FFBS for all coarsest regions
#'
#' @param gamma Length-n1 vector of discount factors.
#' @param y_coarse T by n1 count matrix.
#' @param xi T by n1 positive baseline mean factor matrix.
#' @param kappa T by n1 positive current kappa matrix.
#' @param a0,b0 Initial Gamma filter shape and rate.
#' @param return_diag If TRUE, return sampled lambda and filtering diagnostics.
#' @return Either a T by n1 lambda_tilde matrix, or a list with diagnostics.
ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0,
                            lambda_floor = .lambda_floor(),
                            min_shape = 1e-12,
                            min_rate = 1e-12,
                            return_diag = FALSE) {
    if (is.null(dim(y_coarse)) || length(dim(y_coarse)) != 2L) {
        stop("y_coarse must be a matrix.", call. = FALSE)
    }
    if (!.is_count_vector(as.vector(y_coarse))) {
        stop("y_coarse must contain finite nonnegative integer-valued counts.",
             call. = FALSE)
    }
    .check_same_dim(xi, y_coarse, "xi", "y_coarse")
    .check_same_dim(kappa, y_coarse, "kappa", "y_coarse")

    if (any(!is.finite(xi)) || any(xi <= 0)) {
        stop("xi must contain positive finite values.", call. = FALSE)
    }
    if (any(!is.finite(kappa)) || any(kappa <= 0)) {
        stop("kappa must contain positive finite values.", call. = FALSE)
    }

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    if (!is.numeric(gamma) || length(gamma) != n1 ||
        any(!is.finite(gamma)) || any(gamma <= 0) || any(gamma >= 1)) {
        stop("gamma must be a finite length-n1 vector with entries in (0, 1).",
             call. = FALSE)
    }
    .check_positive_scalar(a0, "a0")
    .check_positive_scalar(b0, "b0")

    lambda_tilde_new <- matrix(NA_real_, TT, n1)

    if (return_diag) {
        a_filt <- matrix(NA_real_, TT, n1)
        b_filt <- matrix(NA_real_, TT, n1)
    }

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]

        out_j <- ffbs_lambda_j(
            gamma_j = gamma[j],
            y_j = y_coarse[, j],
            zeta_j = zeta_j,
            a0 = a0,
            b0 = b0,
            lambda_floor = lambda_floor,
            min_shape = min_shape,
            min_rate = min_rate,
            return_filter = return_diag
        )

        if (return_diag) {
            lambda_tilde_new[, j] <- out_j$lambda_path
            a_filt[, j] <- out_j$a_filt
            b_filt[, j] <- out_j$b_filt
        } else {
            lambda_tilde_new[, j] <- out_j
        }
    }

    if (!return_diag) {
        return(lambda_tilde_new)
    }

    list(
        lambda_tilde = lambda_tilde_new,
        a_filt = a_filt,
        b_filt = b_filt,
        diag = list(
            mean_lambda = mean(lambda_tilde_new),
            min_lambda = min(lambda_tilde_new),
            max_lambda = max(lambda_tilde_new),
            n_lambda_floor = sum(lambda_tilde_new <= lambda_floor),
            mean_a_filt = mean(a_filt),
            mean_b_filt = mean(b_filt)
        )
    )
}

#' Deterministic re-centering to identified space
#'
#' Constraints imposed by this reparameterization:
#'     mean_t log(lambda_tilde_{t,j}) = 0 for every region j,
#'     mean_j phi_j = 0.
#'
#' The transformation leaves
#'     exp(beta0 + phi_j) * lambda_tilde_{t,j}
#' exactly unchanged up to floating-point error.
#'
#' @param beta0 Current intercept.
#' @param phi Length-n1 spatial random effect vector.
#' @param lambda_tilde T by n1 residual risk matrix.
#' @param phi_center Numeric target mean for phi, default 0.
#' @param return_diag If TRUE, include shift diagnostics and invariance check.
#' @return List with beta0, phi, lambda_tilde, and optionally diagnostics.
recenter <- function(beta0, phi, lambda_tilde,
                     lambda_floor = .lambda_floor(),
                     phi_center = 0,
                     return_diag = FALSE) {
    if (!is.numeric(beta0) || length(beta0) != 1L || !is.finite(beta0)) {
        stop("beta0 must be a finite scalar.", call. = FALSE)
    }
    if (!is.numeric(phi) || any(!is.finite(phi))) {
        stop("phi must be a finite numeric vector.", call. = FALSE)
    }
    if (is.null(dim(lambda_tilde)) || length(dim(lambda_tilde)) != 2L) {
        stop("lambda_tilde must be a matrix.", call. = FALSE)
    }
    if (length(phi) != ncol(lambda_tilde)) {
        stop("length(phi) must equal ncol(lambda_tilde).", call. = FALSE)
    }
    if (any(!is.finite(lambda_tilde)) || any(lambda_tilde < 0)) {
        stop("lambda_tilde must contain finite nonnegative values.", call. = FALSE)
    }
    .check_positive_scalar(lambda_floor, "lambda_floor")
    if (!is.numeric(phi_center) || length(phi_center) != 1L || !is.finite(phi_center)) {
        stop("phi_center must be a finite scalar.", call. = FALSE)
    }

    lambda_old <- pmax(lambda_tilde, lambda_floor)
    log_lambda_old <- log(lambda_old)

    s <- colMeans(log_lambda_old)

    ## The global shift is chosen so that the re-centered phi has mean phi_center.
    ## If mean(phi) is already zero and phi_center = 0, this reduces to mean(s).
    s_bar <- mean(phi + s) - phi_center

    lambda_new <- sweep(lambda_old, 2L, exp(s), FUN = "/")
    lambda_new <- pmax(lambda_new, lambda_floor)

    phi_new <- phi + s - s_bar
    beta0_new <- beta0 + s_bar

    if (!return_diag) {
        return(list(
            beta0 = beta0_new,
            phi = phi_new,
            lambda_tilde = lambda_new
        ))
    }

    log_core_old <- outer(rep(1, nrow(lambda_old)), beta0 + phi, "+") +
        log(lambda_old)
    log_core_new <- outer(rep(1, nrow(lambda_new)), beta0_new + phi_new, "+") +
        log(lambda_new)

    list(
        beta0 = beta0_new,
        phi = phi_new,
        lambda_tilde = lambda_new,
        diag = list(
            region_log_lambda_shift = s,
            global_shift = s_bar,
            mean_phi_before = mean(phi),
            mean_phi_after = mean(phi_new),
            max_abs_region_log_lambda_mean_after = max(abs(colMeans(log(lambda_new)))),
            max_abs_log_core_difference = max(abs(log_core_new - log_core_old)),
            n_lambda_floor_before = sum(lambda_tilde <= lambda_floor),
            n_lambda_floor_after = sum(lambda_new <= lambda_floor)
        )
    )
}

# Backward-compatible aliases, in case an older sampler uses these names.
update_lambda_ffbs <- ffbs_lambda_all
recenter_lambda_phi_beta0 <- recenter
