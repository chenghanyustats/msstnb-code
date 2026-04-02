## ======================================================================
## update_kappa.R
## Stable conjugate update for the NB random effects kappa_{t,j}
##
## Full conditional (unchanged mathematically):
##   kappa_{t,j} | y, lambda_tilde, xi, r
##     ~ Gamma(y_{t,j} + r_j, r_j + xi_{t,j} lambda_tilde_{t,j})
##
## This rewrite is mainly for numerical stability and vectorization.
## The larger sampler fix is to STOP conditioning beta/phi on sampled kappa
## when updating baseline mean structure; that is handled in the revised
## update_regression / update_icar / sampler files.
## ======================================================================

update_kappa <- function(y_coarse, lambda_tilde, xi, r,
                         min_rate = 1e-12,
                         return_diag = FALSE) {

    TT <- nrow(y_coarse)
    n1 <- ncol(y_coarse)

    if (length(r) != n1) stop("length(r) must equal ncol(y_coarse)")

    ## Poisson special case: if r is infinite, kappa is degenerate at 1.
    if (all(is.infinite(r))) {
        kappa_new <- matrix(1, TT, n1)
        if (return_diag) {
            return(list(kappa = kappa_new,
                        diag = list(mean_kappa = 1,
                                    min_rate_used = NA_real_)))
        }
        return(kappa_new)
    }

    ## Shape/rate matrices
    shape_mat <- sweep(y_coarse, 2, r, FUN = "+")
    rate_mat  <- sweep(xi * lambda_tilde, 2, r, FUN = "+")
    rate_mat  <- pmax(rate_mat, min_rate)

    kappa_new <- matrix(
        rgamma(length(shape_mat),
               shape = as.vector(shape_mat),
               rate  = as.vector(rate_mat)),
        nrow = TT, ncol = n1
    )

    if (!return_diag) return(kappa_new)

    list(kappa = kappa_new,
         diag = list(
             mean_kappa    = mean(kappa_new),
             sd_kappa      = sd(as.vector(kappa_new)),
             min_rate_used = min(rate_mat),
             mean_rate     = mean(rate_mat)
         ))
}
