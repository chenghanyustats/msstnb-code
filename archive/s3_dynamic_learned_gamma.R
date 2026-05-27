# ============================================================
# Scenario 3: dynamic latent risk with learned gamma
# ============================================================
# Purpose:
#   Fit the full dynamic MSSTNB model to the same dynamic data used by
#   Scenario 2, but learn gamma inside the MCMC instead of fixing gamma
#   at the truth.
#
# Main comparison:
#   Scenario 2: dynamic DGP, dynamic fit, gamma fixed at truth.
#   Scenario 3: dynamic DGP, dynamic fit, gamma learned by MCMC.
#
# This script is intentionally a thin driver. It should use the package
# level code in R/ and R/mcmc/ rather than duplicating sampler internals.
# ============================================================

# ------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------
.s3_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) {
        cat(..., "\n")
    }
}

.s3_assert <- function(cond, msg) {
    if (!isTRUE(cond)) {
        stop(msg, call. = FALSE)
    }
}

.s3_norm_path <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.s3_find_existing_file <- function(candidates, label = "file") {
    candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0L) {
        stop(
            "Could not find ", label, ". Tried:\n",
            paste0("  ", candidates, collapse = "\n"),
            call. = FALSE
        )
    }
    hit[[1L]]
}

.s3_source_if_exists <- function(path, verbose = TRUE) {
    if (file.exists(path)) {
        .s3_msg("Sourcing: ", path, verbose = verbose)
        source(path, local = .GlobalEnv)
        return(TRUE)
    }
    FALSE
}

.s3_source_project_code <- function(root = ".", verbose = TRUE) {
    root <- .s3_norm_path(root)
    r_dir <- file.path(root, "R")
    mcmc_dir <- file.path(root, "R", "mcmc")

    .s3_assert(dir.exists(r_dir), paste0("Cannot find R directory: ", r_dir))

    preferred_first <- c(
        file.path(root, "00_setup.R"),
        file.path(r_dir, "00_setup.R"),
        file.path(root, "mcmc_config.R"),
        file.path(r_dir, "mcmc_config.R")
    )

    sourced <- character(0L)
    for (f in preferred_first) {
        if (.s3_source_if_exists(f, verbose = verbose)) {
            sourced <- c(sourced, .s3_norm_path(f))
        }
    }

    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    mcmc_files <- character(0L)
    if (dir.exists(mcmc_dir)) {
        mcmc_files <- list.files(mcmc_dir, pattern = "\\.R$", full.names = TRUE)
    }

    all_files <- unique(c(r_files, mcmc_files))
    all_files <- all_files[!(.s3_norm_path(all_files) %in% sourced)]

    # Source non setup files in alphabetical order. This is robust when files
    # are function definitions rather than scripts with side effects.
    all_files <- sort(all_files)
    for (f in all_files) {
        .s3_source_if_exists(f, verbose = verbose)
    }

    invisible(TRUE)
}

.s3_infer_TT_from_data <- function(data_obj) {
    if (!is.null(data_obj$y_levels) && length(data_obj$y_levels) >= 1L) {
        return(nrow(data_obj$y_levels[[1L]]))
    }
    if (!is.null(data_obj$y)) {
        return(nrow(data_obj$y))
    }
    if (!is.null(data_obj$y1)) {
        return(nrow(data_obj$y1))
    }
    NA_integer_
}

.s3_find_data_file <- function(
    data_dir,
    rep_id,
    TT = NULL,
    data_scenario = "S2",
    data_file_template = NULL
) {
    rep2 <- sprintf("%02d", rep_id)
    rep4 <- sprintf("%04d", rep_id)

    if (!is.null(data_file_template)) {
        cand <- sprintf(data_file_template, rep_id)
        return(.s3_find_existing_file(cand, label = "Scenario 3 data file"))
    }

    candidates <- c(
        file.path(data_dir, data_scenario, sprintf("data_rep%s.rds", rep2)),
        file.path(data_dir, data_scenario, sprintf("data_rep%s_T%s.rds", rep2, TT)),
        file.path(data_dir, data_scenario, sprintf("rep%s.rds", rep4)),
        file.path(data_dir, data_scenario, sprintf("rep%s_T%s.rds", rep4, TT)),
        file.path(data_dir, data_scenario, sprintf("rep%s.rds", rep2)),
        file.path(data_dir, data_scenario, sprintf("data_rep%s.rds", rep4)),
        file.path(data_dir, sprintf("%s_data_rep%s.rds", data_scenario, rep2)),
        file.path(data_dir, sprintf("data_rep%s.rds", rep2)),
        file.path(data_dir, sprintf("rep%s.rds", rep4)),
        file.path(data_dir, sprintf("rep%s.rds", rep2))
    )

    .s3_find_existing_file(candidates, label = "Scenario 3 data file")
}

.s3_extract_data_and_truth <- function(ds) {
    data_obj <- if (!is.null(ds$data)) ds$data else ds
    truth <- if (!is.null(ds$truth)) ds$truth else NULL
    list(data = data_obj, truth = truth)
}

.s3_make_default_hyper <- function(data_obj, truth = NULL, gamma_prior = c(1, 1)) {
    # Prefer existing project constructors when available.
    if (exists("make_msstnb_hyperparameters", mode = "function")) {
        out <- tryCatch(
            make_msstnb_hyperparameters(data = data_obj, truth = truth),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }
    if (exists("make_hyper", mode = "function")) {
        out <- tryCatch(
            make_hyper(data = data_obj, truth = truth),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }
    if (exists("default_hyper", mode = "function")) {
        out <- tryCatch(
            default_hyper(data = data_obj, truth = truth),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }

    # Minimal fallback. The sampler may require additional fields depending on
    # the current package implementation. If so, use hyper_override below.
    list(
        gamma_prior = gamma_prior,
        gamma_a = gamma_prior[[1L]],
        gamma_b = gamma_prior[[2L]]
    )
}

.s3_make_default_mcmc <- function(n_iter, burn, thin, save_lambda = TRUE) {
    if (exists("make_msstnb_mcmc_config", mode = "function")) {
        out <- tryCatch(
            make_msstnb_mcmc_config(n_iter = n_iter, burn = burn, thin = thin),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }
    if (exists("make_mcmc_config", mode = "function")) {
        out <- tryCatch(
            make_mcmc_config(n_iter = n_iter, burn = burn, thin = thin),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }
    if (exists("default_mcmc", mode = "function")) {
        out <- tryCatch(
            default_mcmc(n_iter = n_iter, burn = burn, thin = thin),
            error = function(e) NULL
        )
        if (!is.null(out)) return(out)
    }

    list(
        n_iter = as.integer(n_iter),
        burn = as.integer(burn),
        thin = as.integer(thin),
        save_lambda = isTRUE(save_lambda),
        save_lambda_tilde = isTRUE(save_lambda)
    )
}

.s3_find_sampler <- function() {
    candidates <- c(
        "msstnb_algo2_mcmc",
        "msstnb_mcmc",
        "run_msstnb_mcmc",
        "sampler_msstnb",
        "run_sampler"
    )
    hit <- candidates[vapply(candidates, exists, logical(1L), mode = "function")]
    if (length(hit) == 0L) {
        stop(
            "Cannot find the MCMC sampler function. Expected one of: ",
            paste(candidates, collapse = ", "),
            call. = FALSE
        )
    }
    get(hit[[1L]], mode = "function")
}

.s3_set_learned_gamma_controls <- function(hyper, mcmc, gamma_prior = c(1, 1)) {
    # Store controls in both hyper and mcmc because different versions of the
    # sampler have used different locations for these flags.
    hyper$fixed_gamma <- NULL
    hyper$gamma_fixed <- NULL
    hyper$learn_gamma <- TRUE
    hyper$update_gamma <- TRUE
    hyper$gamma_mode <- "learned"
    hyper$gamma_prior <- gamma_prior
    hyper$gamma_a <- gamma_prior[[1L]]
    hyper$gamma_b <- gamma_prior[[2L]]

    mcmc$fixed_gamma <- NULL
    mcmc$gamma_fixed <- NULL
    mcmc$learn_gamma <- TRUE
    mcmc$update_gamma <- TRUE
    mcmc$gamma_mode <- "learned"
    mcmc$save_gamma <- TRUE
    mcmc$save_lambda <- TRUE
    mcmc$save_lambda_tilde <- TRUE

    list(hyper = hyper, mcmc = mcmc)
}

.s3_merge_list <- function(x, y) {
    if (is.null(y)) return(x)
    .s3_assert(is.list(y), "Override object must be a list.")
    for (nm in names(y)) {
        x[[nm]] <- y[[nm]]
    }
    x
}

.s3_call_sampler <- function(sampler, data_obj, hyper, mcmc, init, seed, verbose) {
    fml <- names(formals(sampler))
    args <- list()

    if ("data" %in% fml) args$data <- data_obj else if ("data_obj" %in% fml) args$data_obj <- data_obj
    if ("hyper" %in% fml) args$hyper <- hyper
    if ("mcmc" %in% fml) args$mcmc <- mcmc
    if ("init" %in% fml) args$init <- init
    if ("seed" %in% fml) args$seed <- seed
    if ("verbose" %in% fml) args$verbose <- verbose

    do.call(sampler, args)
}

# ------------------------------------------------------------
# Main Scenario 3 driver
# ------------------------------------------------------------
source_s3_dynamic_learned_gamma <- function(
    root = ".",
    reps = 1:20,
    TT = NULL,
    data_dir = "data_revised",
    data_scenario = "S2",
    data_file_template = NULL,
    fit_dir = "fits_s3_dynamic_learned_gamma",
    analysis_dir = "analysis_s3_dynamic_learned_gamma",
    n_iter = 20000L,
    burn = 5000L,
    thin = 10L,
    gamma_prior = c(1, 1),
    seed_base = 300000L,
    hyper_override = NULL,
    mcmc_override = NULL,
    init = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
    scenario_name <- "S3_DYNAMIC_LEARNED_GAMMA"

    root <- .s3_norm_path(root)
    data_dir <- if (grepl("^/", data_dir)) data_dir else file.path(root, data_dir)
    fit_dir <- if (grepl("^/", fit_dir)) fit_dir else file.path(root, fit_dir)
    analysis_dir <- if (grepl("^/", analysis_dir)) analysis_dir else file.path(root, analysis_dir)

    dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

    .s3_source_project_code(root = root, verbose = verbose)
    sampler <- .s3_find_sampler()

    manifest <- data.frame(
        scenario = character(0L),
        rep_id = integer(0L),
        TT = integer(0L),
        data_file = character(0L),
        fit_file = character(0L),
        status = character(0L),
        stringsAsFactors = FALSE
    )

    fits <- vector("list", length(reps))
    names(fits) <- paste0("rep", sprintf("%04d", reps))

    for (ii in seq_along(reps)) {
        rep_id <- as.integer(reps[[ii]])
        .s3_msg("", verbose = verbose)
        .s3_msg("==============================", verbose = verbose)
        .s3_msg("Scenario 3 learned gamma", verbose = verbose)
        .s3_msg("rep = ", rep_id, verbose = verbose)
        .s3_msg("==============================", verbose = verbose)

        data_file <- .s3_find_data_file(
            data_dir = data_dir,
            rep_id = rep_id,
            TT = TT,
            data_scenario = data_scenario,
            data_file_template = data_file_template
        )

        ds <- readRDS(data_file)
        extracted <- .s3_extract_data_and_truth(ds)
        data_obj <- extracted$data
        truth <- extracted$truth
        TT_i <- .s3_infer_TT_from_data(data_obj)
        if (!is.null(TT) && !is.na(TT_i) && TT_i != TT) {
            warning("Requested TT = ", TT, " but data file has TT = ", TT_i, ". Using data TT.")
        }
        TT_label <- ifelse(is.na(TT_i), ifelse(is.null(TT), NA_integer_, TT), TT_i)

        fit_file <- file.path(
            fit_dir,
            sprintf("s3_dynamic_learned_gamma_rep%02d_T%s.rds", rep_id, TT_label)
        )

        if (file.exists(fit_file) && !isTRUE(overwrite)) {
            .s3_msg("Existing fit found. Reading: ", fit_file, verbose = verbose)
            fit <- readRDS(fit_file)
            fits[[ii]] <- fit
            manifest <- rbind(
                manifest,
                data.frame(
                    scenario = scenario_name,
                    rep_id = rep_id,
                    TT = TT_label,
                    data_file = data_file,
                    fit_file = fit_file,
                    status = "existing",
                    stringsAsFactors = FALSE
                )
            )
            next
        }

        hyper <- .s3_make_default_hyper(data_obj = data_obj, truth = truth, gamma_prior = gamma_prior)
        mcmc <- .s3_make_default_mcmc(n_iter = n_iter, burn = burn, thin = thin, save_lambda = TRUE)
        controls <- .s3_set_learned_gamma_controls(hyper, mcmc, gamma_prior = gamma_prior)
        hyper <- .s3_merge_list(controls$hyper, hyper_override)
        mcmc <- .s3_merge_list(controls$mcmc, mcmc_override)

        hyper$scenario_name <- scenario_name
        mcmc$scenario_name <- scenario_name

        seed <- as.integer(seed_base + rep_id)
        set.seed(seed)

        fit <- .s3_call_sampler(
            sampler = sampler,
            data_obj = data_obj,
            hyper = hyper,
            mcmc = mcmc,
            init = init,
            seed = seed,
            verbose = verbose
        )

        fit$scenario_name <- scenario_name
        fit$rep_id <- rep_id
        fit$TT <- TT_label
        fit$data_file <- data_file
        fit$truth <- truth
        fit$scenario3_controls <- list(
            gamma_prior = gamma_prior,
            fixed_gamma = NULL,
            learn_gamma = TRUE,
            update_gamma = TRUE,
            n_iter = n_iter,
            burn = burn,
            thin = thin,
            seed = seed
        )

        saveRDS(fit, fit_file)
        fits[[ii]] <- fit

        manifest <- rbind(
            manifest,
            data.frame(
                scenario = scenario_name,
                rep_id = rep_id,
                TT = TT_label,
                data_file = data_file,
                fit_file = fit_file,
                status = "new",
                stringsAsFactors = FALSE
            )
        )
    }

    manifest_file <- file.path(analysis_dir, "scenario3_fit_manifest.csv")
    write.csv(manifest, manifest_file, row.names = FALSE)
    .s3_msg("Saved manifest: ", manifest_file, verbose = verbose)

    invisible(list(
        scenario = scenario_name,
        fits = fits,
        manifest = manifest,
        manifest_file = manifest_file
    ))
}
