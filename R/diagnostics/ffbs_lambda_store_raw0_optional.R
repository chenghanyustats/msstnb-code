## ==========================================================================
## ffbs_lambda.R
## FFBS for coarsest residual risk lambda_tilde, with optional raw 0:T storage
##
## New in this version:
##   - ffbs_lambda_j_fullpath() returns the raw path lambda_{0:T}
##   - ffbs_lambda_all_with_raw0() returns both lambda_{1:T} and lambda_{0:T}
##   - original ffbs_lambda_all() is preserved for backward compatibility
## ==========================================================================

ffbs_lambda_j <- function(gamma_j, y_j, zeta_j, a0, b0) {
    TT <- length(y_j)

    a_filt <- numeric(TT)
    b_filt <- numeric(TT)

    a_filt[1] <- gamma_j * a0 + y_j[1]
    b_filt[1] <- gamma_j * b0 + zeta_j[1]

    if (TT >= 2L) {
        for (t in 2:TT) {
            a_filt[t] <- gamma_j * a_filt[t - 1] + y_j[t]
            b_filt[t] <- gamma_j * b_filt[t - 1] + zeta_j[t]
        }
    }

    lambda_path <- numeric(TT)
    LAMBDA_FLOOR <- .Machine$double.xmin

    lambda_path[TT] <- max(rgamma(1, shape = a_filt[TT], rate = b_filt[TT]),
                           LAMBDA_FLOOR)

    if (TT >= 2L) {
        for (t in TT:2) {
            shape_U <- (1 - gamma_j) * a_filt[t - 1]
            rate_U  <- b_filt[t - 1]
            U <- rgamma(1, shape = shape_U, rate = rate_U)
            lambda_path[t - 1] <- U + gamma_j * lambda_path[t]
        }
    }

    lambda_path
}

ffbs_lambda_j_fullpath <- function(gamma_j, y_j, zeta_j, a0, b0) {
    TT <- length(y_j)

    a_filt <- numeric(TT)
    b_filt <- numeric(TT)

    a_filt[1] <- gamma_j * a0 + y_j[1]
    b_filt[1] <- gamma_j * b0 + zeta_j[1]

    if (TT >= 2L) {
        for (t in 2:TT) {
            a_filt[t] <- gamma_j * a_filt[t - 1] + y_j[t]
            b_filt[t] <- gamma_j * b_filt[t - 1] + zeta_j[t]
        }
    }

    lambda_0T <- numeric(TT + 1L)
    LAMBDA_FLOOR <- .Machine$double.xmin

    ## terminal draw for lambda_T
    lambda_0T[TT + 1L] <- max(rgamma(1, shape = a_filt[TT], rate = b_filt[TT]),
                              LAMBDA_FLOOR)

    ## backward draws for lambda_{T-1},...,lambda_1
    if (TT >= 2L) {
        for (t in TT:2) {
            shape_U <- (1 - gamma_j) * a_filt[t - 1]
            rate_U  <- b_filt[t - 1]
            U <- rgamma(1, shape = shape_U, rate = rate_U)
            lambda_0T[t] <- U + gamma_j * lambda_0T[t + 1L]
        }
    } else {
        ## nothing to do; lambda_1 is already terminal draw
    }

    ## one more backward step to sample raw lambda_0
    ## uses the t=0 filtered state (a0, b0)
    U0 <- rgamma(1, shape = (1 - gamma_j) * a0, rate = b0)
    lambda_0T[1L] <- max(U0 + gamma_j * lambda_0T[2L], LAMBDA_FLOOR)

    ## numerical safety
    lambda_0T <- pmax(lambda_0T, LAMBDA_FLOOR)

    lambda_0T
}

ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0) {
    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    lambda_tilde_new <- matrix(NA_real_, TT, n1)

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]
        lambda_tilde_new[, j] <- ffbs_lambda_j(gamma[j], y_coarse[, j], zeta_j, a0, b0)
    }

    lambda_tilde_new
}

ffbs_lambda_all_with_raw0 <- function(gamma, y_coarse, xi, kappa, a0, b0) {
    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    lambda_tilde_new <- matrix(NA_real_, TT, n1)
    lambda_raw_0T    <- array(NA_real_, dim = c(TT + 1L, n1))

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]
        lam0T <- ffbs_lambda_j_fullpath(gamma[j], y_coarse[, j], zeta_j, a0, b0)
        lambda_raw_0T[, j] <- lam0T
        lambda_tilde_new[, j] <- lam0T[2:(TT + 1L)]
    }

    list(lambda_tilde = lambda_tilde_new, lambda_raw_0T = lambda_raw_0T)
}

recenter <- function(beta0, phi, lambda_tilde) {
    TT <- nrow(lambda_tilde)
    n1 <- ncol(lambda_tilde)

    lambda_tilde <- pmax(lambda_tilde, .Machine$double.xmin)
    s <- colMeans(log(lambda_tilde))
    s_bar <- mean(s)

    for (j in seq_len(n1)) {
        lambda_tilde[, j] <- lambda_tilde[, j] * exp(-s[j])
    }

    phi_new <- phi + s - s_bar
    beta0_new <- beta0 + s_bar

    list(beta0 = beta0_new, phi = phi_new, lambda_tilde = lambda_tilde)
}
