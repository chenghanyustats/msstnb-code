source("R/00_setup.R")
source("R/01_helpers.R")
source("R/02_simulate_data.R")
source("R/03_scenarios.R")
source("R/04_run_dgp.R")


dat <- sanity_check()
export_setup_csv()
export_scenario_csv()

run_all_dgp(
    reps = 1:10,
    scenarios = "S1",
    data_dir = "data",
    verbose = TRUE
)



# run_all_dgp(reps = 1:3, scenarios = c("S1", "S2"))
# run_all_dgp()