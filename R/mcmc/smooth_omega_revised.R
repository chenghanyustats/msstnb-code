## -----------------------------------------------------------------------------
## smooth_omega_revised.R
##
## Smoothed posterior sampling for dynamic multiscale splitting probabilities.
##
## Model for one parent node j with K children:
##     y_t | y_{t+}, omega_t ~ Multinomial(y_{t+}, omega_t),
##     omega_t | D_{t-1}, delta ~ Dirichlet(delta * c_{t-1}),
##     c_t = delta * c_{t-1} + y_t.
##
## This file implements a Dirichlet forward filter and backward smoother.
## It intentionally does not rely on external rdirichlet() implementations.
## -----------------------------------------------------------------------------

## ---- small numerical helpers ------------------------------------------------

.check_positive_scalar <- function(x, name) {
    if (length(x) != 1L || !is.finite(x) || x <= 0) {
        stop(name, " must be a finite positive scalar.")
    }
    invisible(TRUE)
}

.check_probability_scalar <- function(x, name, eps = 1e-12) {
    if (length(x) != 1L || !is.finite(x) || x <= eps || x >= 1 - eps) {
        stop(name, " must lie strictly in (", eps, ", ", 1 - eps, ").")
    }
    invisible(TRUE)
}

.check_count_matrix <- function(y, name = "y_kids") {
    if (!is.matrix(y)) {
        stop(name, " must be a matrix.")
    }
    if (any(!is.finite(y))) {
        stop(name, " must contain only finite values.")
    }
    if (any(y < 0)) {
        stop(name, " must contain nonnegative counts.")
    }
    if (any(abs(y - round(y)) > 1e-8)) {
        stop(name, " must be integer-like counts.")
    }
    invisible(TRUE)
}

## Draw one Dirichlet vector.  This avoids depending on gtools or MCMCpack.
.rdirichlet1 <- function(alpha, min_conc = 1e-12) {
    if (any(!is.finite(alpha)) || any(alpha <= 0)) {
        stop("Dirichlet concentration parameters must be positive and finite.")
    }
    alpha <- pmax(alpha, min_conc)
    g <- rgamma(length(alpha), shape = alpha, rate = 1)

    ## Extremely small concentration parameters can occasionally underflow to all
    ## zeros.  If that happens, fall back to a one-hot draw according to alpha.
    s <- sum(g)
    if (!is.finite(s) || s <= 0) {
        out <- rep(0, length(alpha))
        out[sample.int(length(alpha), size = 1, prob = alpha)] <- 1
        return(out)
    }

    g / s
}

## ---- forward filtering ------------------------------------------------------

omega_forward_filter_j <- function(delta, y_kids, c0,
                                   min_conc = 1e-12,
                                   check_inputs = TRUE) {
    if (check_inputs) {
        .check_probability_scalar(delta, "delta")
        .check_count_matrix(y_kids, "y_kids")
        if (length(c0) != ncol(y_kids)) {
            stop("length(c0) must equal ncol(y_kids).")
        }
        if (any(!is.finite(c0)) || any(c0 <= 0)) {
            stop("c0 must contain positive finite values.")
        }
        .check_positive_scalar(min_conc, "min_conc")
    }

    TT <- nrow(y_kids)
    K <- ncol(y_kids)

    c_filt <- matrix(NA_real_, nrow = TT, ncol = K)
    c_prev <- pmax(as.numeric(c0), min_conc)

    for (t in seq_len(TT)) {
        c_prior <- pmax(delta * c_prev, min_conc)
        c_filt[t, ] <- c_prior + y_kids[t, ]
        c_prev <- c_filt[t, ]
    }

    c_filt
}

## ---- backward smoothing for one parent node ---------------------------------

smooth_omega_j <- function(delta, y_kids, c0,
                           min_conc = 1e-12,
                           check_inputs = TRUE,
                           return_diag = FALSE) {
    if (check_inputs) {
        .check_probability_scalar(delta, "delta")
        .check_count_matrix(y_kids, "y_kids")
        if (length(c0) != ncol(y_kids)) {
            stop("length(c0) must equal ncol(y_kids).")
        }
        if (any(!is.finite(c0)) || any(c0 <= 0)) {
            stop("c0 must contain positive finite values.")
        }
        .check_positive_scalar(min_conc, "min_conc")
    }

    TT <- nrow(y_kids)
    K <- ncol(y_kids)

    c_filt <- omega_forward_filter_j(
        delta = delta,
        y_kids = y_kids,
        c0 = c0,
        min_conc = min_conc,
        check_inputs = FALSE
    )

    omega_path <- matrix(NA_real_, nrow = TT, ncol = K)
    S_draws <- rep(NA_real_, TT)

    ## Terminal smoothed draw.
    omega_path[TT, ] <- .rdirichlet1(c_filt[TT, ], min_conc = min_conc)

    ## Backward smoother.  Guard TT = 1 because in R, 1:2 is c(1, 2), not empty.
    if (TT >= 2L) {
        for (t in TT:2) {
            c_prev <- pmax(c_filt[t - 1L, ], min_conc)
            c_prev_plus <- sum(c_prev)

            shape_s1 <- pmax(delta * c_prev_plus, min_conc)
            shape_s2 <- pmax((1 - delta) * c_prev_plus, min_conc)

            S_t <- rbeta(1, shape1 = shape_s1, shape2 = shape_s2)
            omega_resid <- .rdirichlet1((1 - delta) * c_prev,
                                        min_conc = min_conc)

            omega_path[t - 1L, ] <- (1 - S_t) * omega_resid +
                S_t * omega_path[t, ]

            ## Numerical cleanup: keep exact simplex normalization.
            omega_path[t - 1L, ] <- pmax(omega_path[t - 1L, ], 0)
            omega_path[t - 1L, ] <- omega_path[t - 1L, ] /
                sum(omega_path[t - 1L, ])

            S_draws[t] <- S_t
        }
    }

    row_sum_error <- max(abs(rowSums(omega_path) - 1))

    if (!return_diag) {
        return(omega_path)
    }

    list(
        omega = omega_path,
        c_filt = c_filt,
        diag = list(
            row_sum_error = row_sum_error,
            min_omega = min(omega_path),
            max_omega = max(omega_path),
            mean_omega = mean(omega_path),
            S_draws = S_draws,
            mean_S = mean(S_draws, na.rm = TRUE),
            min_concentration = min(c_filt),
            max_concentration = max(c_filt)
        )
    )
}

## ---- smoothing for all parent nodes -----------------------------------------

smooth_omega_all <- function(delta, y_fine, c0,
                             min_conc = 1e-12,
                             check_inputs = TRUE,
                             return_diag = FALSE) {
    if (check_inputs) {
        .check_probability_scalar(delta, "delta")
        if (length(dim(y_fine)) != 3L) {
            stop("y_fine must be a 3-dimensional array with dim = c(TT, n1, K).")
        }
        if (any(!is.finite(y_fine))) {
            stop("y_fine must contain only finite values.")
        }
        if (any(y_fine < 0)) {
            stop("y_fine must contain nonnegative counts.")
        }
        if (any(abs(y_fine - round(y_fine)) > 1e-8)) {
            stop("y_fine must be integer-like counts.")
        }
        if (length(c0) != dim(y_fine)[3]) {
            stop("length(c0) must equal dim(y_fine)[3].")
        }
        if (any(!is.finite(c0)) || any(c0 <= 0)) {
            stop("c0 must contain positive finite values.")
        }
        .check_positive_scalar(min_conc, "min_conc")
    }

    TT <- dim(y_fine)[1]
    n1 <- dim(y_fine)[2]
    K <- dim(y_fine)[3]

    omega_new <- array(NA_real_, dim = c(TT, n1, K))

    if (return_diag) {
        c_filt_all <- array(NA_real_, dim = c(TT, n1, K))
        row_sum_error_j <- rep(NA_real_, n1)
        mean_S_j <- rep(NA_real_, n1)
        min_conc_j <- rep(NA_real_, n1)
        max_conc_j <- rep(NA_real_, n1)
    }

    for (j in seq_len(n1)) {
        y_kids_j <- matrix(y_fine[, j, ], nrow = TT, ncol = K)

        out_j <- smooth_omega_j(
            delta = delta,
            y_kids = y_kids_j,
            c0 = c0,
            min_conc = min_conc,
            check_inputs = FALSE,
            return_diag = return_diag
        )

        if (return_diag) {
            omega_new[, j, ] <- out_j$omega
            c_filt_all[, j, ] <- out_j$c_filt
            row_sum_error_j[j] <- out_j$diag$row_sum_error
            mean_S_j[j] <- out_j$diag$mean_S
            min_conc_j[j] <- out_j$diag$min_concentration
            max_conc_j[j] <- out_j$diag$max_concentration
        } else {
            omega_new[, j, ] <- out_j
        }
    }

    ## Global simplex diagnostic across all t, j.
    omega_sums <- apply(omega_new, c(1, 2), sum)
    max_row_sum_error <- max(abs(omega_sums - 1))

    if (!return_diag) {
        return(omega_new)
    }

    list(
        omega = omega_new,
        c_filt = c_filt_all,
        diag = list(
            max_row_sum_error = max_row_sum_error,
            row_sum_error_by_node = row_sum_error_j,
            mean_S_by_node = mean_S_j,
            mean_S_overall = mean(mean_S_j, na.rm = TRUE),
            min_omega = min(omega_new),
            max_omega = max(omega_new),
            mean_omega = mean(omega_new),
            min_concentration_by_node = min_conc_j,
            max_concentration_by_node = max_conc_j,
            min_concentration_overall = min(min_conc_j),
            max_concentration_overall = max(max_conc_j)
        )
    )
}

## ---- backward-compatible aliases -------------------------------------------

update_omega_smoother <- smooth_omega_all
smooth_omega <- smooth_omega_all
