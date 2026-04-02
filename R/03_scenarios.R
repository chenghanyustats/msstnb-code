## ==========================================================================
## 03_scenarios.R
## Scenario configurations for the MSSTNB simulation study
##
## STEP 4 REVISION:
## For a given replication, S1--S4 share the SAME baseline phi_star so that
## scenario comparisons isolate the intended parameter changes (r, gamma, etc.).
## Only S5 (changepoint) generates an additional phi_after.
## ==========================================================================


## ======================================================================
## SINGLE SOURCE OF TRUTH: change values here, everything else follows
## ======================================================================
DEFAULTS <- list(
    beta0   = -1,
    beta    = c(0.5, -0.4),
    r       = 10,         # scalar; rep()'d to n1 inside build_scenario_params
    gamma   = 0.8,        # baseline gamma
    tau_phi = 2,
    delta   = 0.9
)

## Each scenario: name, overrides (only fields that differ from DEFAULTS),
## and the DGP type.
SCENARIO_SPEC <- list(
    S1 = list(
        name      = "Baseline",
        overrides = list(),
        dgp_type  = "standard"
    ),
    S2 = list(
        name      = "Heavy overdispersion",
        overrides = list(r = 0.5),
        dgp_type  = "standard"
    ),
    S3 = list(
        name      = "Near Poisson",
        overrides = list(r = 100),
        dgp_type  = "standard"
    ),
    S4 = list(
        name      = "Fast dynamics",
        overrides = list(gamma = 0.55),
        dgp_type  = "standard"
    ),
    S5 = list(
        name      = "Changepoint",
        overrides = list(),
        dgp_type  = "changepoint",
        change_t  = 30L
    )
)


## ======================================================================
## Shared phi generation within a replication
## ======================================================================

#' Generate the baseline ICAR field shared by S1--S4 within one replication
#'
#' @param tau_phi ICAR precision scalar
#' @param H       Graph Laplacian (from 00_setup.R)
#' @param B       ICAR basis (from 00_setup.R)
#' @return        Numeric vector phi_star of length n1
build_shared_phi <- function(tau_phi, H, B) {
    generate_icar(tau_phi = tau_phi, H = H, B = B)
}


## ======================================================================
## Build resolved parameter list for one scenario + one replication
## ======================================================================

#' @param scenario_id  "S1"..."S5"
#' @param H            Graph Laplacian (from 00_setup.R)
#' @param B            ICAR basis (from 00_setup.R)
#' @param n1           Number of coarsest regions
#' @param phi_shared   Optional shared baseline phi_star for this replication.
#'                     If supplied, S1--S4 use this exact field.
#' @return             Named list of all true parameters (scalars expanded)
build_scenario_params <- function(scenario_id, H, B, n1, phi_shared = NULL) {

    spec <- SCENARIO_SPEC[[scenario_id]]
    if (is.null(spec)) stop("Unknown scenario_id: ", scenario_id)

    ## ---- merge defaults + overrides ----
    resolved <- DEFAULTS
    for (nm in names(spec$overrides)) {
        resolved[[nm]] <- spec$overrides[[nm]]
    }

    ## ---- expand scalars to vectors of length n1 ----
    resolved$r     <- rep(resolved$r,     length.out = n1)
    resolved$gamma <- rep(resolved$gamma, length.out = n1)

    ## ---- baseline ICAR field ----
    ## S1--S4 should share the same phi_star within a replication when
    ## phi_shared is supplied. S5 also uses the same baseline phi_star and then
    ## adds phi_after for the changepoint.
    if (is.null(phi_shared)) {
        phi_star <- generate_icar(tau_phi = resolved$tau_phi, H = H, B = B)
    } else {
        phi_star <- as.numeric(phi_shared)
        if (length(phi_star) != n1) {
            stop("phi_shared must have length n1.")
        }
        phi_star <- phi_star - mean(phi_star)
    }

    phi_after <- NULL
    change_t  <- NULL

    if (spec$dgp_type == "changepoint") {
        phi_raw   <- generate_icar(tau_phi = resolved$tau_phi, H = H, B = B)
        phi_after <- 0.5 * phi_star + 0.5 * phi_raw
        phi_after <- phi_after - mean(phi_after)
        change_t  <- spec$change_t
    }

    ## ---- assemble ----
    list(
        scenario_id = scenario_id,
        name        = spec$name,
        dgp_type    = spec$dgp_type,
        beta0       = resolved$beta0,
        beta        = resolved$beta,
        phi         = phi_star,
        phi_after   = phi_after,
        change_t    = change_t,
        tau_phi     = resolved$tau_phi,
        r           = resolved$r,
        gamma       = resolved$gamma,
        delta       = resolved$delta
    )
}


#' Build all 5 scenario parameter sets for one replication
#'
#' S1--S4 share the same baseline phi_star. S5 uses that same phi_star plus a
#' newly generated phi_after for the post-changepoint period.
all_scenarios <- function(H, B, n1) {
    ids <- names(SCENARIO_SPEC)

    phi_shared <- build_shared_phi(tau_phi = DEFAULTS$tau_phi, H = H, B = B)

    out <- lapply(ids, function(sid) {
        build_scenario_params(
            scenario_id = sid,
            H = H,
            B = B,
            n1 = n1,
            phi_shared = phi_shared
        )
    })
    names(out) <- ids
    out
}


## ======================================================================
## Auto-export: scenario_summary.csv that stays in sync with DEFAULTS
## ======================================================================

#' Write dictionaries/scenario_summary.csv from DEFAULTS + SCENARIO_SPEC
#'
#' Call this after changing any value in DEFAULTS or SCENARIO_SPEC.
#'
#' @param out_dir  Directory for the CSV (default "dictionary")
export_scenario_csv <- function(out_dir = "dictionary") {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    rows <- lapply(names(SCENARIO_SPEC), function(sid) {
        spec     <- SCENARIO_SPEC[[sid]]
        resolved <- DEFAULTS
        for (nm in names(spec$overrides)) {
            resolved[[nm]] <- spec$overrides[[nm]]
        }
        data.frame(
            scenario_id   = sid,
            name          = spec$name,
            beta0         = resolved$beta0,
            beta1         = resolved$beta[1],
            beta2         = resolved$beta[2],
            r_star        = resolved$r,
            gamma_star    = resolved$gamma,
            tau_phi_star  = resolved$tau_phi,
            delta_star    = resolved$delta,
            dgp_type      = spec$dgp_type,
            change_t      = ifelse(is.null(spec$change_t), NA, spec$change_t),
            phi_design    = ifelse(sid %in% c("S1", "S2", "S3", "S4"),
                                   "shared_within_rep", "shared_baseline_plus_phi_after"),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)

    out_path <- file.path(out_dir, "scenario_summary.csv")
    write.csv(df, out_path, row.names = FALSE)
    cat("Wrote", out_path, "\n")
    invisible(df)
}


## ======================================================================
## Pretty-print
## ======================================================================

print_scenario <- function(params) {
    cat(sprintf("=== %s: %s ===\n", params$scenario_id, params$name))
    cat(sprintf("  beta0 = %.2f,  beta = (%.2f, %.2f)\n",
                params$beta0, params$beta[1], params$beta[2]))
    cat(sprintf("  tau_phi = %.1f,  delta = %.2f\n",
                params$tau_phi, params$delta))
    cat(sprintf("  r*     = %s\n", paste(unique(params$r), collapse = ", ")))
    cat(sprintf("  gamma* = %s\n", paste(unique(params$gamma), collapse = ", ")))
    cat(sprintf("  phi* range: [%.3f, %.3f], sum = %.1e\n",
                min(params$phi), max(params$phi), sum(params$phi)))
    if (!is.null(params$change_t)) {
        cat(sprintf("  Changepoint at t = %d\n", params$change_t))
        cat(sprintf("  phi_after range: [%.3f, %.3f]\n",
                    min(params$phi_after), max(params$phi_after)))
    }
    cat("\n")
}


#' Print all scenarios with resolved values (quick review)
print_all_scenarios <- function() {
    cat("========== DEFAULTS ==========\n")
    cat(sprintf("  beta0 = %.2f,  beta = (%s)\n",
                DEFAULTS$beta0,
                paste(DEFAULTS$beta, collapse = ", ")))
    cat(sprintf("  r = %s,  gamma = %s,  tau_phi = %s,  delta = %s\n",
                DEFAULTS$r, DEFAULTS$gamma, DEFAULTS$tau_phi, DEFAULTS$delta))
    cat("==============================\n\n")

    for (sid in names(SCENARIO_SPEC)) {
        spec <- SCENARIO_SPEC[[sid]]
        ov   <- spec$overrides
        if (length(ov) == 0) {
            diff_str <- "(no overrides)"
        } else {
            diff_str <- paste(names(ov), "=", ov, collapse = ", ")
        }
        shared_str <- if (sid %in% c("S1", "S2", "S3", "S4")) {
            "shared phi_star within replication"
        } else {
            "shared phi_star + new phi_after"
        }
        cat(sprintf("  %s  %-25s  %s  |  %s\n",
                    sid, spec$name, diff_str, shared_str))
    }
    cat("\n")
}
