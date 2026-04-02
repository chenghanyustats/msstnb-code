source("R/mcmc/run_single_fit.R")
source_mcmc()

fit <- fit_one(
    scenario_id = "S1",
    rep_id = 1,
    method = "M1",
    data_dir = "data",
    output_dir = "output",
    verbose = 500L
)

source("R/mcmc/analysis.R")

plot_mcmc_one(fit, "beta0", 1)
plot_mcmc_one(fit, "gamma", 1)
plot_mcmc_one(fit, "beta", 1)
plot_mcmc_one(fit, "beta", 2)


data_obj <- readRDS("data/S1/data_rep01.rds")

data_obj$beta0_star_ident
data_obj$beta_star

plot(fit$samples$beta0, type = "l", main = "M1 beta0")
plot(fit$samples$beta[,1], type = "l", main = "M1 beta1")
plot(fit$samples$beta[,2], type = "l", main = "M1 beta2")
plot(fit$samples$gamma[,1], type = "l", main = "M1 gamma1")

par(mfrow = c(2, 2))
plot(fit$samples$phi[,1], type = "l", main = "M1 phi1")
plot(fit$samples$phi[,2], type = "l", main = "M1 phi2")
acf(fit$samples$phi[,1], main = "ACF phi1")
acf(fit$samples$phi[,2], main = "ACF phi2")


source("R/mcmc/check.R")

data_obj <- readRDS("data/S1/data_rep01.rds")

res <- posterior_summary_against_truth(fit, data_obj)
short_table(res)

lam1 <- lambda_region_summary(fit, data_obj, j = 1)
head(lam1)
mean(lam1$covered)

par(mfrow = c(3, 2))
plot_gamma_r_key(
    fit = fit,
    data_obj = data_obj,
    gamma_idx = c(4, 6, 9),
    r_idx = c(4, 3, 2)
)


#
# fit_M3 <- fit_one(
#     scenario_id = "S1",
#     rep_id = 1,
#     method = "M3",
#     data_dir = "data",
#     output_dir = "output",
#     verbose = 500L
# )
#
# par(mfrow = c(2, 2))
# plot(fit_M3$samples$beta0, type = "l", main = "M3 beta0")
# plot(fit_M3$samples$beta[,1], type = "l", main = "M3 beta1")
# plot(fit_M3$samples$beta[,2], type = "l", main = "M3 beta2")
# plot(fit_M3$samples$gamma[,1], type = "l", main = "M3 gamma1")
