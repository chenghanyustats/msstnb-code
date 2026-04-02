## ==========================================================================
## mcmc_config.R
## MCMC configuration: priors, iteration settings, and method flags
##
## DESIGN: All tunable MCMC settings live in MCMC_PRIORS and MCMC_SETTINGS.
## Method variants (M1–M4) are controlled by flags in MCMC_SETTINGS.
## ==========================================================================


## ======================================================================
## PRIORS — change here to adjust prior distributions
## ======================================================================
MCMC_PRIORS <- list(

    ## ---- regression (Eq. 17) ----
    beta0_mean  = 0,          # prior mean for β₀
    beta0_sd    = 10,         # prior sd for β₀
    beta_mean   = c(0, 0),    # prior mean for (β₁, β₂)
    beta_sd     = c(5, 5),    # prior sd for (β₁, β₂) — diagonal

    ## ---- ICAR precision (Eq. 19) ----
    tau_phi_shape = 1,        # Ga(a_φ, b_φ)
    tau_phi_rate  = 0.01,

    ## ---- NB dispersion (Eq. 20, Remark 3.4) ----
    r_shape = 1,              # Ga(a_r, b_r) with a_r ≥ 1
    r_rate  = 0.1,

    ## ---- discount factors (Eq. 20) ----
    gamma_a = 1,              # Beta(a_γ, b_γ) — final adopted default
    gamma_b = 1,
    delta_a = 5,              # Beta(a_δ, b_δ)
    delta_b = 1
)


## ======================================================================
## MCMC SETTINGS — iteration control, tuning, method flags
## ======================================================================
MCMC_SETTINGS <- list(

    ## ---- iterations ----
    n_iter   = 20000L,        # total MCMC iterations
    n_burnin = 10000L,        # burn-in (not stored)
    n_thin   = 5L,            # thinning interval

    ## NOTE: MH proposal SDs are initialised in sampler.R (initialise_state)
    ## and adapted automatically during burn-in. No manual tuning needed.

    ## ---- method flags ----
    ## Set these to select which model variant to fit:
    ##   M1 (full):       include_nb=T, include_icar=T, include_covariates=T
    ##   M2 (no NB):      include_nb=F, include_icar=T, include_covariates=T
    ##   M3 (no ICAR):    include_nb=T, include_icar=F, include_covariates=T
    ##   M4 (FF2017):     include_nb=F, include_icar=F, include_covariates=F
    include_nb         = TRUE,
    include_icar       = TRUE,
    include_covariates = TRUE
)


## ======================================================================
## Convenience: build a method label from flags
## ======================================================================
method_label <- function(settings = MCMC_SETTINGS) {
    if (settings$include_nb && settings$include_icar && settings$include_covariates) {
        return("M1_full")
    } else if (!settings$include_nb && settings$include_icar && settings$include_covariates) {
        return("M2_no_nb")
    } else if (settings$include_nb && !settings$include_icar && settings$include_covariates) {
        return("M3_no_icar")
    } else if (!settings$include_nb && !settings$include_icar && !settings$include_covariates) {
        return("M4_ff2017")
    } else {
        return("custom")
    }
}
