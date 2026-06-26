## ==========================================================================
## smooth_omega.R
## Dirichlet filtering and backward smoothing for split proportions ω
##
## Paper references:
##   Eq. 45   — Dirichlet filtering: c_t = δ c_{t-1} + y_children
##   Eqs. 46–48 — backward smoother:
##     S_t ~ Beta(δ c_{t-1,+}, (1-δ) c_{t-1,+})
##     ω̃_{t-1} ~ Dir((1-δ) c_{t-1})
##     ω_{t-1} = (1-S_t) ω̃_{t-1} + S_t ω_t
##
## The proof of correctness (Appendix A.8) uses the Gamma representation
## of the Dirichlet: if independent Gamma variables are normalised,
## they form a Dirichlet. The backward step decomposes the filtered
## Dirichlet at time t-1 into a convex combination of the smoothed
## value at time t and an independent residual.
## ==========================================================================


#' Dirichlet forward filter + backward sampler for one internal node (l,j)
#'
#' @param delta    Accepted split discount factor
#' @param y_kids   T × K matrix of child counts for node (l,j)
#' @param c0       Length-K initial Dirichlet concentration
#' @return         T × K matrix of smoothed split proportions ω_{1:T,lj}
smooth_omega_j <- function(delta, y_kids, c0) {

    TT <- nrow(y_kids)
    K  <- ncol(y_kids)

    ## ---- Forward filter (Eq. 45) ----
    ## Store filtered concentrations for t = 1,...,T
    ## c_filt[t, ] = filtered concentration at time t
    c_filt <- matrix(NA_real_, TT, K)

    ## t = 1: prior from t=0 state
    c_filt[1, ] <- delta * c0 + y_kids[1, ]

    ## t = 2,...,T
    for (t in 2:TT) {
        c_filt[t, ] <- delta * c_filt[t - 1, ] + y_kids[t, ]
    }

    ## ---- Backward sampler (Eqs. 46–48) ----
    omega_path <- matrix(NA_real_, TT, K)

    ## Terminal draw: ω_T ~ Dir(c_T)
    omega_path[TT, ] <- rdirichlet(1, c_filt[TT, ])

    ## Backward steps: for t = T, T-1, ..., 2, compute ω_{t-1}
    for (t in TT:2) {

        ## Filtered concentration at time t-1
        c_prev <- c_filt[t - 1, ]

        c_prev_plus <- sum(c_prev)    # c_{t-1,+}

        ## Eq. 46: mixing weight
        S_t <- rbeta(1,
                     shape1 = delta * c_prev_plus,
                     shape2 = (1 - delta) * c_prev_plus)

        ## Eq. 47: independent residual Dirichlet
        omega_tilde <- rdirichlet(1, (1 - delta) * c_prev)

        ## Eq. 48: convex combination
        omega_path[t - 1, ] <- (1 - S_t) * omega_tilde + S_t * omega_path[t, ]
    }

    return(omega_path)
}


#' Run Dirichlet smoothing for all internal nodes
#'
#' @param delta    Accepted split discount factor (scalar)
#' @param y_fine   T × n1 × K array of child counts
#' @param c0       Length-K initial Dirichlet concentration
#' @return         T × n1 × K array of smoothed split proportions
smooth_omega_all <- function(delta, y_fine, c0) {

    TT <- dim(y_fine)[1]
    n1 <- dim(y_fine)[2]
    K  <- dim(y_fine)[3]

    omega_new <- array(NA_real_, dim = c(TT, n1, K))

    for (j in seq_len(n1)) {
        omega_new[, j, ] <- smooth_omega_j(delta, y_fine[, j, ], c0)
    }

    return(omega_new)
}
