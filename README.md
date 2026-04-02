# MSSTNB Simulation Study

## Purpose
Generate simulated data for 5 scenarios × 50 replications to validate the
Dynamic Multiscale Spatiotemporal Negative Binomial (MSSTNB) model.
Targeting JASA Theory & Methods.

## Quick start

```r
source("R/00_setup.R")        # load packages, build spatial structure
source("R/01_helpers.R")       # helper functions
source("R/02_simulate_data.R") # main DGP function
source("R/03_scenarios.R")     # scenario configurations
source("R/04_run_dgp.R")       # generate all data (or run interactively)
```

## Folder structure

```
├── R/
│   ├── 00_setup.R            # packages, spatial graph, shared constants
│   ├── 01_helpers.R          # generate_icar(), rdirichlet(), etc.
│   ├── 02_simulate_data.R    # simulate_one_dataset(): the full DGP
│   ├── 03_scenarios.R        # scenario_config(): returns list of 5 scenarios
│   └── 04_run_dgp.R          # run_all_dgp(): loop scenarios × reps, save
├── data/
│   ├── shared/               # covariates & exposures (shared across scenarios)
│   │   └── shared_inputs_rep{RR}.rds
│   ├── S1/ ... S5/           # simulated datasets per scenario
│   │   └── data_rep{RR}.rds  # one file per replication
│   └── (created by 04_run_dgp.R)
├── dictionaries/
│   ├── param_dictionary.csv          # master parameter dictionary
│   ├── data_dictionary.csv           # what's inside each .rds file
│   └── scenario_summary.csv          # true values by scenario
└── README.md
```

## Scenario overview

| ID | Name               | Key manipulation     | r*   | γ*  | τ*_φ | δ*   |
|----|--------------------|----------------------|------|-----|------|------|
| S1 | Baseline           | None                 | 5    | 0.9 | 3    | 0.95 |
| S2 | Heavy OD           | Small r*             | 0.5  | 0.9 | 3    | 0.95 |
| S3 | Near Poisson       | Large r*             | 100  | 0.9 | 3    | 0.95 |
| S4 | Fast dynamics      | Small γ*             | 5    | 0.7 | 3    | 0.95 |
| S5 | Changepoint        | φ shifts at t = 31   | 5    | 0.9 | 3    | 0.95 |

## Naming conventions

- `rep` or `RR`: replication index (01–50), zero-padded to 2 digits
- `S1`–`S5`: scenario ID
- `*_star` suffix: true (DGP) parameter value
- `*_hat` suffix: posterior estimate (used later in fitting)
