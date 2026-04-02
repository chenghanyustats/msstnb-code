# Final adopted MSSTNB ordinary sampler package

This folder contains the current final adopted version of the ordinary MSSTNB sampler,
after the debugging sequence completed in this project.

## Final adopted defaults
- Ordinary sampler uses the **kappa-collapsed** baseline updates for beta and phi.
- Recommended default gamma prior is **Beta(1,1)**.
- Core production sampler does **not** store raw lambda paths for gamma diagnostics.
- Optional diagnostic scripts are included separately in `R/diagnostics/`.

## Core scripts required to run the model
Located in `R/` and `R/mcmc/`:
- `R/00_setup.R`
- `R/01_helpers.R`
- `R/mcmc/mcmc_config.R`
- `R/mcmc/mcmc_utils.R`
- `R/mcmc/update_kappa.R`
- `R/mcmc/update_regression.R`
- `R/mcmc/update_icar.R`
- `R/mcmc/update_dispersion.R`
- `R/mcmc/update_gamma.R`
- `R/mcmc/ffbs_lambda.R`
- `R/mcmc/update_delta.R`
- `R/mcmc/smooth_omega.R`
- `R/mcmc/sampler.R`
- `R/mcmc/run_single_fit.R`

## Recommended default settings
In `R/mcmc/mcmc_config.R`:
- `gamma_a = 1`
- `gamma_b = 1`
- `delta_a = 5`
- `delta_b = 1`
- other defaults left unchanged from the working project version

## Core file roles
- `update_kappa.R`: conjugate kappa update, numerically cleaned up
- `update_regression.R`: beta update using kappa-collapsed marginal likelihood
- `update_icar.R`: phi update using kappa-collapsed marginal likelihood
- `sampler.R`: ordinary sampler with revised block ordering; baseline updates occur before kappa draw
- `update_gamma.R`: ordinary marginal predictive gamma update (logic verified via internal/external audit)

## Optional diagnostics
Located in `R/diagnostics/`:
- `test_free_fit_check.R`
- `batch_rerun_s1_summary.R`
- `s1_longchains_beta11.R`
- `kappa_diagnostic.R`
- `rerun_gamma_prior_beta11.R`
- advanced gamma diagnostics:
  - `ffbs_lambda_store_raw0_optional.R`
  - `sampler_store_raw_lambda_for_gamma_optional.R`
  - `frozen_a_gamma_diagnostic.R`
  - `run_frozen_a_on_existing_fit.R`
  - `update_gamma_audit_version.R`
  - `run_minimal_gamma_audit.R`
  - `compare_gamma_internal_external.R`

## Recommended usage path
1. Use the core scripts in `R/` and `R/mcmc/` for production ordinary fits.
2. Use `R/mcmc/run_single_fit.R` as the standard entry point.
3. Only switch to the optional diagnostic scripts if you need additional auditing.
