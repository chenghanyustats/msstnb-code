# MCMC Sampler for the MSSTNB Model

## Files (in Algorithm 3 step order)

| File | Algorithm 3 step | Paper reference | What it does |
|------|-----------------|-----------------|-------------|
| `mcmc_config.R` | — | Eqs. 17, 19, 20 | Prior hyperparameters, MCMC settings, method flags |
| `mcmc_utils.R` | — | Murray et al. (2010) | ESS algorithm, logit/expit, log-sum-exp |
| `update_kappa.R` | Step 2 | Eq. 30 | Conjugate κ update: Ga(y+r, r+ξλ̃) |
| `update_regression.R` | Step 3 | Eq. 31 | ESS for (β₀, β) with Gaussian prior |
| `update_icar.R` | Steps 4–5 | Eqs. 32–34 | ESS for φ (reduced coords) + conjugate τ_φ |
| `update_dispersion.R` | Step 6 | Eq. 35 | MH for r on log scale |
| `update_gamma.R` | Step 7 | Algorithm 1, Eqs. 40–43 | MH for γ on logit scale, marginal likelihood |
| `ffbs_lambda.R` | Steps 8–9 | Algorithm 2, Eqs. 24–27 | FFBS for λ̃ + re-centering |
| `update_delta.R` | Step 10 | Eqs. 49–50 | MH for δ on logit scale |
| `smooth_omega.R` | Step 11 | Eqs. 45–48 | Dirichlet FFBS for ω |
| `sampler.R` | Full loop | Algorithm 3 | `run_one_iteration()` + `run_mcmc()` |
| `run_single_fit.R` | — | — | Load data, configure method, run, save |

## Method variants (controlled by flags)

```
M1 (full):    include_nb=T  include_icar=T  include_covariates=T
M2 (no NB):   include_nb=F  include_icar=T  include_covariates=T
M3 (no ICAR): include_nb=T  include_icar=F  include_covariates=T
M4 (FF2017):  include_nb=F  include_icar=F  include_covariates=F
```

## Quick start

```r
# From the simulation/ directory:
source("R/mcmc/run_single_fit.R")   # sources everything
source_mcmc()

# Sanity check (2000 iter, ~1 min)
fit <- sanity_check_mcmc()

# Full fit (20000 iter)
fit <- fit_one("S1", rep_id = 1, method = "M1")

# Access results
str(fit$samples$beta0)           # 2000 posterior draws
str(fit$diagnostics)             # acceptance rates, trace, timing
```

## Output structure

Each `fit_*.rds` file contains:
- `samples`: thinned post-burn-in posterior draws for all parameters
- `diagnostics`: log-likelihood trace, acceptance rates, ESS rejections, timing
- `settings`: the MCMC settings used
- `priors`: the prior hyperparameters used
- `method`: string label ("M1_full", etc.)
- `n_stored`: number of stored samples
