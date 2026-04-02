## ==========================================================================
## ffbs_lambda.R
## Forward Filter Backward Sampler for the coarsest residual risk λ̃
## and the deterministic re-centering step
##
## Paper references:
##   Algorithm 2, Proposition 4.4 — backward kernel for λ̃
##   Eqs. 38      — Gamma filtering recursion (forward pass)
##   Eqs. 24–27   — deterministic re-centering to identified space C
## ==========================================================================


#' Forward Filter Backward Sampler for λ̃_{1:T,1j} (Algorithm 2)
#'
#' Given the accepted γ_j and current conditioning variables, this function:
#'   1. Runs the Gamma filter forward (Eq. 38) to get (a_t, b_t) for t=1,...,T
#'   2. Draws λ̃_T from the filtered Ga(a_T, b_T)
#'   3. Draws λ̃_{t-1} | λ̃_t backwards using Prop. 4.4:
#'        λ̃_{t-1} = U + γ λ̃_t,  U ~ Ga((1-γ) a_{t-1}, b_{t-1})
#'
#' @param gamma_j  Accepted discount factor for region j
#' @param y_j      Length-T vector of coarsest counts for region j
#' @param zeta_j   Length-T vector of effective offsets ζ = ξ·κ for region j
#' @param a0       Initial Gamma filter shape
#' @param b0       Initial Gamma filter rate
#' @return         Length-T vector of sampled λ̃_{1:T,1j}
ffbs_lambda_j <- function(gamma_j, y_j, zeta_j, a0, b0) {

    TT <- length(y_j)

    ## ---- Forward filter (Eq. 38) ----
    ## Store filtered (a_t, b_t) for t = 1,...,T
    a_filt <- numeric(TT)
    b_filt <- numeric(TT)

    ## t = 1: prior from t=0 state
    a_filt[1] <- gamma_j * a0 + y_j[1]
    b_filt[1] <- gamma_j * b0 + zeta_j[1]

    ## t = 2,...,T
    for (t in 2:TT) {
        a_filt[t] <- gamma_j * a_filt[t - 1] + y_j[t]
        b_filt[t] <- gamma_j * b_filt[t - 1] + zeta_j[t]
    }

    ## ---- Backward sampler (Proposition 4.4) ----
    lambda_path <- numeric(TT)

    ## Numerical floor: rgamma can return exactly 0.0 in double precision
    ## when shape is very small or rate is very large. A zero λ̃ value
    ## causes log(λ̃) = -Inf in the re-centering step, which cascades to
    ## β₀ = -Inf. We clamp at the smallest positive normalised double.
    LAMBDA_FLOOR <- .Machine$double.xmin   # ≈ 2.2e-308

    ## Draw terminal state: λ̃_T ~ Ga(a_T, b_T)
    lambda_path[TT] <- max(rgamma(1, shape = a_filt[TT], rate = b_filt[TT]),
                           LAMBDA_FLOOR)

    ## Backward steps: for t = T, T-1, ..., 2, compute λ̃_{t-1}
    for (t in TT:2) {
        ## Filtered values at time t-1: a_filt[t-1], b_filt[t-1]
        ## (At t=2 this is a_filt[1], b_filt[1] — the filtered state at time 1)

        ## U ~ Ga((1-γ) a_{t-1}, b_{t-1})   (Prop. 4.4)
        shape_U <- (1 - gamma_j) * a_filt[t - 1]
        rate_U  <- b_filt[t - 1]
        U <- rgamma(1, shape = shape_U, rate = rate_U)
        ## U = 0 is fine here: lambda_path[t-1] = 0 + γ λ̃_t > 0

        ## λ̃_{t-1} = U + γ λ̃_t
        lambda_path[t - 1] <- U + gamma_j * lambda_path[t]
    }

    return(lambda_path)
}


#' Run FFBS for all coarsest regions
#'
#' @param gamma        Length-n1 vector of accepted discount factors
#' @param y_coarse     T × n1 count matrix
#' @param xi           T × n1 effective offset matrix (Eq. 29)
#' @param kappa        T × n1 current κ matrix
#' @param a0, b0       Initial Gamma filter parameters
#' @return             T × n1 matrix of sampled λ̃
ffbs_lambda_all <- function(gamma, y_coarse, xi, kappa, a0, b0) {

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    lambda_tilde_new <- matrix(NA_real_, TT, n1)

    for (j in seq_len(n1)) {
        zeta_j <- xi[, j] * kappa[, j]    # ζ = ξ·κ  (Eq. 36)
        lambda_tilde_new[, j] <- ffbs_lambda_j(gamma[j], y_coarse[, j],
                                               zeta_j, a0, b0)
    }

    return(lambda_tilde_new)
}


#' Deterministic re-centering (Eqs. 24–27)
#'
#' Projects (β₀, φ, λ̃) onto the identified space C defined in Eq. 28:
#'   Σ_j φ_j = 0  AND  (1/T) Σ_t log λ̃_{t,j} = 0 for all j
#'
#' This resolves the level-shift equivalence class (Section 3.6).
#'
#' @param beta0        Current intercept
#' @param phi          Length-n1 current ICAR vector
#' @param lambda_tilde T × n1 current residual risk matrix
#' @return             List with: beta0, phi, lambda_tilde (all re-centered)
recenter <- function(beta0, phi, lambda_tilde) {

    TT <- nrow(lambda_tilde)
    n1 <- ncol(lambda_tilde)

    ## Safety: clamp λ̃ away from zero before taking log.
    ## rgamma can return exactly 0.0 in double precision for extreme
    ## parameter combinations, which would make log(λ̃) = -Inf and
    ## cascade to β₀ = -Inf.  The floor is astronomically small and
    ## has no effect on inference.
    lambda_tilde <- pmax(lambda_tilde, .Machine$double.xmin)

    ## Eq. 24: temporal mean of log λ̃ for each region
    s <- colMeans(log(lambda_tilde))   # length-n1 vector

    ## Global mean of the shifts
    s_bar <- mean(s)

    ## Eq. 25: re-center λ̃ (multiplicative on original scale)
    for (j in seq_len(n1)) {
        lambda_tilde[, j] <- lambda_tilde[, j] * exp(-s[j])
    }

    ## Eq. 26: absorb region-specific shifts into φ
    phi_new <- phi + s - s_bar

    ## Eq. 27: absorb global shift into β₀
    beta0_new <- beta0 + s_bar

    return(list(
        beta0        = beta0_new,
        phi          = phi_new,
        lambda_tilde = lambda_tilde
    ))
}
