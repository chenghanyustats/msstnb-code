## ==========================================================================
## update_kappa.R
## Conjugate update for the NB random effects κ_{t,1j}
##
## Paper reference: Eq. 30
##   κ_{t,1j} | y, λ̃, β₀, β, φ, r ~ Ga(y_{t,1j} + r_{1j},  r_{1j} + ξ_{t,1j} λ̃_{t,1j})
##
## where ξ_{t,1j} = e_{t,1j} exp(β₀ + x'β + φ_j)   (Eq. 29)
## ==========================================================================


#' Update all κ_{t,1j} from their conjugate full conditionals
#'
#' @param y_coarse     T × n1 integer matrix of observed counts
#' @param lambda_tilde T × n1 matrix of current residual risks
#' @param xi           T × n1 matrix of effective offsets (Eq. 29)
#' @param r            Length-n1 vector of current dispersion parameters
#' @return             T × n1 matrix of updated κ values
update_kappa <- function(y_coarse, lambda_tilde, xi, r) {

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    kappa_new <- matrix(NA_real_, TT, n1)

    for (j in seq_len(n1)) {
        for (t in seq_len(TT)) {

            ## Full conditional:  Ga(y + r,  r + ξ·λ̃)
            shape_post <- y_coarse[t, j] + r[j]
            rate_post  <- r[j] + xi[t, j] * lambda_tilde[t, j]

            kappa_new[t, j] <- rgamma(1, shape = shape_post, rate = rate_post)
        }
    }

    return(kappa_new)
}
