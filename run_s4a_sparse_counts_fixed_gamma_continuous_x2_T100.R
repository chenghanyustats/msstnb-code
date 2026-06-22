## ============================================================================
## run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100.R
## MSSTNB Scenario 4A fixed-gamma fitting script
##
## Scenario:
##   S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100
##
## Purpose:
##   Fit the dynamic-lambda MSSTNB model to the revised Scenario 4A sparse-count
##   stress-test data generated with the continuous-time x2 covariate design:
##
##       x2_mode     = "continuous_time"
##       x2_ar       = 0.50
##       x2_innov_sd = 0.80
##
##   Gamma is fixed at the truth, gamma = 0.80, matching the fixed-gamma
##   stress-test convention used for the revised Scenario 4 workflow.
##
## Recommended use from the project root:
##
##   source("run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100.R")
##
##   ## Smoke test, one replicate first:
##   out_test <- run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100(
##       root = ".",
##       reps = 1,
##       overwrite_existing = TRUE,
##       verbose = TRUE
##   )
##
##   ## Formal run:
##   out <- run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100(
##       root = ".",
##       reps = 1:10,
##       overwrite_existing = FALSE,
##       verbose = TRUE
##   )
##
## Notes:
##   1. This script is path-safe and does not assume the working directory except
##      through root = ".".
##   2. It first tries to reuse an existing project-level fixed-gamma fitting
##      function.  If your project uses a different function name, pass it with
##      fit_fun = "your_function_name" or fit_engine_script = "path/to/script.R".
##   3. The script saves every replicate fit as an RDS file and writes a manifest
##      plus status-count table under output_s4a_sparse_counts_continuous_x2/.
## ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.s4a_fit_dir_create <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(path)
}

.s4a_fit_msg <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
    invisible(NULL)
}

.s4a_fit_require_file <- function(path, label = "file") {
    if (!file.exists(path)) {
        stop("Required ", label, " not found: ", path, call. = FALSE)
    }
    invisible(path)
}

.s4a_fit_safe_source <- function(path, verbose = TRUE) {
    .s4a_fit_require_file(path, label = "R script")
    .s4a_fit_msg("source: ", path, verbose = verbose)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.s4a_fit_find_first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0L) return(NULL)
    hits[[1L]]
}

.s4a_fit_source_optional <- function(path, verbose = TRUE) {
    if (!is.null(path) && file.exists(path)) {
        .s4a_fit_safe_source(path, verbose = verbose)
        return(TRUE)
    }
    FALSE
}

.s4a_fit_source_project_dependencies <- function(root = ".",
                                                 s3_script = NULL,
                                                 s4a_data_script = NULL,
                                                 fit_engine_script = NULL,
                                                 verbose = TRUE) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)

    ## Source S4A data script when available.  This gives access to validation
    ## helpers such as validate_s4a_continuous_x2_design().
    if (is.null(s4a_data_script)) {
        s4a_data_script <- .s4a_fit_find_first_existing(c(
            file.path(root, "s4a_sparse_counts_continuous_x2.R"),
            file.path(root, "scripts", "s4a_sparse_counts_continuous_x2.R"),
            file.path(root, "R", "s4a_sparse_counts_continuous_x2.R"),
            file.path(root, "R", "scenarios", "s4a_sparse_counts_continuous_x2.R")
        ))
    } else {
        s4a_data_script <- file.path(root, s4a_data_script)
    }
    .s4a_fit_source_optional(s4a_data_script, verbose = verbose)

    ## Source Scenario 3 script when available.  In the current project this also
    ## loads the MCMC dependencies used by the dynamic-lambda model.
    if (is.null(s3_script)) {
        s3_script <- .s4a_fit_find_first_existing(c(
            file.path(root, "s3_dynamic_learned_gamma.R"),
            file.path(root, "scripts", "s3_dynamic_learned_gamma.R"),
            file.path(root, "R", "s3_dynamic_learned_gamma.R"),
            file.path(root, "R", "scenarios", "s3_dynamic_learned_gamma.R")
        ))
    } else {
        s3_script <- file.path(root, s3_script)
    }
    .s4a_fit_source_optional(s3_script, verbose = verbose)

    if (exists("source_s3_dynamic_learned_gamma", mode = "function", inherits = TRUE)) {
        get("source_s3_dynamic_learned_gamma", mode = "function", inherits = TRUE)(
            root = root,
            verbose = verbose
        )
    }

    ## Source Scenario 2 fixed-gamma script when present.  This is usually the
    ## closest fitting engine for S4A/S4D fixed-gamma runs.
    s2_script <- .s4a_fit_find_first_existing(c(
        file.path(root, "s2_dynamic_fixed_gamma.R"),
        file.path(root, "scripts", "s2_dynamic_fixed_gamma.R"),
        file.path(root, "R", "s2_dynamic_fixed_gamma.R"),
        file.path(root, "R", "scenarios", "s2_dynamic_fixed_gamma.R")
    ))
    .s4a_fit_source_optional(s2_script, verbose = verbose)

    if (exists("source_s2_dynamic_fixed_gamma", mode = "function", inherits = TRUE)) {
        get("source_s2_dynamic_fixed_gamma", mode = "function", inherits = TRUE)(
            root = root,
            verbose = verbose
        )
    }

    ## Optional explicit fitting-engine script.  Use this when the project-level
    ## engine lives in a run script with a different name.
    if (!is.null(fit_engine_script)) {
        fit_engine_script <- file.path(root, fit_engine_script)
        .s4a_fit_safe_source(fit_engine_script, verbose = verbose)
    }

    invisible(TRUE)
}

.s4a_fit_candidate_fit_functions <- function() {
    c(
        ## Most likely project-level one-replicate fixed-gamma engines.
        "fit_s4a_sparse_counts_fixed_gamma_continuous_x2_one",
        "fit_s4a_sparse_counts_fixed_gamma_one",
        "fit_s4a_fixed_gamma_one",
        "run_s4a_sparse_counts_fixed_gamma_continuous_x2_one",
        "run_s4a_sparse_counts_fixed_gamma_one",
        "run_s4a_fixed_gamma_one",

        ## S4D/S2 dynamic fixed-gamma engines, reusable for S4A data.
        "fit_s4d_short_time_fixed_gamma_continuous_x2_one",
        "run_s4d_short_time_fixed_gamma_continuous_x2_one",
        "fit_s2_dynamic_fixed_gamma_one",
        "run_s2_dynamic_fixed_gamma_one",
        "fit_one_s2_dynamic_fixed_gamma",
        "run_one_s2_dynamic_fixed_gamma",
        "fit_dynamic_fixed_gamma_one",
        "run_dynamic_fixed_gamma_one",
        "fit_msstnb_fixed_gamma_one",
        "run_msstnb_fixed_gamma_one",
        "fit_one_replicate_dynamic_fixed_gamma",
        "fit_one_replicate_s2_dynamic_fixed_gamma",

        ## Older matched-gamma engine name from the gamma-grid diagnostics.
        "fit_one_replicate_matched_gamma_omega_recenter_gamma_fixed_truth_delta_truth",
        "run_one_replicate_matched_gamma_omega_recenter_gamma_fixed_truth_delta_truth",

        ## Generic sampler fallbacks.  These are intentionally last.
        "fit_one_replicate",
        "run_one_replicate",
        "fit_ms_stnb_one",
        "run_ms_stnb_one",
        "sampler",
        "run_sampler"
    )
}

.s4a_fit_list_available_engines <- function() {
    objs <- ls(envir = .GlobalEnv)
    funs <- objs[vapply(objs, function(x) exists(x, mode = "function", envir = .GlobalEnv), logical(1))]
    patterns <- c("s4a", "s4d", "s2", "dynamic", "fixed", "gamma", "sampler", "fit", "mcmc")
    keep <- Reduce(`|`, lapply(patterns, function(p) grepl(p, funs, ignore.case = TRUE)))
    sort(funs[keep])
}

.s4a_fit_resolve_engine <- function(fit_fun = NULL, verbose = TRUE) {
    if (is.function(fit_fun)) {
        .s4a_fit_msg("Using supplied fit_fun function object.", verbose = verbose)
        return(fit_fun)
    }

    if (is.character(fit_fun) && length(fit_fun) == 1L) {
        if (!exists(fit_fun, mode = "function", inherits = TRUE)) {
            stop("Requested fit_fun = '", fit_fun, "' was not found in the R session.", call. = FALSE)
        }
        .s4a_fit_msg("Using requested fit function: ", fit_fun, verbose = verbose)
        return(get(fit_fun, mode = "function", inherits = TRUE))
    }

    candidates <- .s4a_fit_candidate_fit_functions()
    hits <- candidates[vapply(candidates, exists, logical(1), mode = "function", inherits = TRUE)]
    if (length(hits) > 0L) {
        .s4a_fit_msg("Using fixed-gamma fitting engine: ", hits[[1L]], verbose = verbose)
        return(get(hits[[1L]], mode = "function", inherits = TRUE))
    }

    available <- .s4a_fit_list_available_engines()
    stop(
        "Could not find a usable fixed-gamma fitting engine.\n",
        "Pass fit_fun = 'your_function_name' or source it with fit_engine_script = 'path/to/script.R'.\n",
        "Candidate names searched:\n  - ", paste(candidates, collapse = "\n  - "), "\n\n",
        "Potentially relevant functions currently available:\n  - ",
        if (length(available) > 0L) paste(available, collapse = "\n  - ") else "<none>",
        call. = FALSE
    )
}

.s4a_fit_call_with_supported_args <- function(fun, args_full, args_conservative) {
    fml <- names(formals(fun))
    has_dots <- "..." %in% fml

    build_args <- function(args) {
        if (has_dots) return(args)
        args[names(args) %in% fml]
    }

    attempts <- list(
        full = build_args(args_full),
        conservative = build_args(args_conservative)
    )

    last_err <- NULL
    for (nm in names(attempts)) {
        aa <- attempts[[nm]]
        if (!has_dots && length(aa) == 0L) next
        res <- try(do.call(fun, aa), silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
        last_err <- res
    }

    formal_msg <- if (length(fml) > 0L) paste(fml, collapse = ", ") else "<none>"
    stop(
        "The selected fixed-gamma fitting engine was found but could not be called.\n",
        "Formal arguments of selected engine: ", formal_msg, "\n",
        "Last error:\n", as.character(last_err),
        call. = FALSE
    )
}

.s4a_fit_get_rep_file <- function(data_dir, rep_id) {
    candidates <- c(
        file.path(data_dir, sprintf("data_rep%02d.rds", as.integer(rep_id))),
        file.path(data_dir, sprintf("data_rep%03d.rds", as.integer(rep_id))),
        file.path(data_dir, sprintf("data_rep%04d.rds", as.integer(rep_id))),
        file.path(data_dir, sprintf("rep%02d.rds", as.integer(rep_id))),
        file.path(data_dir, sprintf("rep%03d.rds", as.integer(rep_id))),
        file.path(data_dir, sprintf("rep%04d.rds", as.integer(rep_id)))
    )
    hit <- .s4a_fit_find_first_existing(candidates)
    if (is.null(hit)) {
        stop(
            "Could not find data file for rep ", rep_id, " in ", data_dir, ".\n",
            "Searched:\n  - ", paste(candidates, collapse = "\n  - "),
            call. = FALSE
        )
    }
    hit
}

.s4a_fit_extract_scalar <- function(x, names, default = NA_real_) {
    for (nm in names) {
        if (!is.null(x[[nm]]) && length(x[[nm]]) >= 1L && is.finite(as.numeric(x[[nm]][1L]))) {
            return(as.numeric(x[[nm]][1L]))
        }
    }
    default
}

.s4a_fit_summarise_data <- function(dat, rep_id, data_file) {
    y <- dat$y_coarse %||% dat$Y %||% dat$y
    y_vec <- as.numeric(y)
    y_vec <- y_vec[is.finite(y_vec)]

    data.frame(
        rep_id = as.integer(rep_id),
        data_file = data_file,
        scenario_id = as.character(dat$scenario_id %||% NA_character_),
        TT = as.integer(dat$TT %||% if (!is.null(dim(y))) dim(y)[1L] else NA_integer_),
        n1 = as.integer(dat$n1 %||% if (!is.null(dim(y))) dim(y)[2L] else NA_integer_),
        stress_type = as.character(dat$stress_type %||% NA_character_),
        x2_mode = as.character(dat$x2_mode %||% NA_character_),
        x2_ar = .s4a_fit_extract_scalar(dat, "x2_ar"),
        x2_innov_sd = .s4a_fit_extract_scalar(dat, "x2_innov_sd"),
        sparse_beta0_shift = .s4a_fit_extract_scalar(dat, "sparse_beta0_shift"),
        beta0_truth_ident = .s4a_fit_extract_scalar(dat, c("beta0_star_ident", "beta0_ident")),
        beta1_truth = if (!is.null(dat$beta_star)) as.numeric(dat$beta_star[1L]) else NA_real_,
        beta2_truth = if (!is.null(dat$beta_star)) as.numeric(dat$beta_star[2L]) else NA_real_,
        r_truth = if (!is.null(dat$r_star)) mean(as.numeric(dat$r_star)) else NA_real_,
        gamma_truth = if (!is.null(dat$gamma_star)) mean(as.numeric(dat$gamma_star)) else NA_real_,
        mean_count = if (length(y_vec)) mean(y_vec) else NA_real_,
        median_count = if (length(y_vec)) stats::median(y_vec) else NA_real_,
        zero_prop = if (length(y_vec)) mean(y_vec == 0) else NA_real_,
        total_count = if (length(y_vec)) sum(y_vec) else NA_real_,
        max_count = if (length(y_vec)) max(y_vec) else NA_real_,
        lambda_raw_median = if (!is.null(dat$lambda_tilde)) stats::median(as.numeric(dat$lambda_tilde)) else NA_real_,
        lambda_raw_max = if (!is.null(dat$lambda_tilde)) max(as.numeric(dat$lambda_tilde)) else NA_real_,
        stringsAsFactors = FALSE
    )
}

.s4a_fit_validate_data <- function(dat,
                                   expected_scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                                   expected_TT = 100L,
                                   expected_n1 = 9L,
                                   expected_x2_mode = "continuous_time",
                                   expected_x2_ar = 0.50,
                                   expected_x2_innov_sd = 0.80,
                                   strict = TRUE) {
    if (isTRUE(strict)) {
        if (!is.null(dat$scenario_id) && !identical(as.character(dat$scenario_id), expected_scenario_id)) {
            stop("Unexpected scenario_id: ", as.character(dat$scenario_id),
                 "; expected ", expected_scenario_id, call. = FALSE)
        }
        if (!is.null(dat$TT) && as.integer(dat$TT) != as.integer(expected_TT)) {
            stop("Unexpected TT: ", dat$TT, "; expected ", expected_TT, call. = FALSE)
        }
        if (!is.null(dat$n1) && as.integer(dat$n1) != as.integer(expected_n1)) {
            stop("Unexpected n1: ", dat$n1, "; expected ", expected_n1, call. = FALSE)
        }
        if (!is.null(dat$x2_mode) && !identical(dat$x2_mode, expected_x2_mode)) {
            stop("Unexpected x2_mode: ", dat$x2_mode, "; expected ", expected_x2_mode, call. = FALSE)
        }
        if (!is.null(dat$x2_ar) && !isTRUE(all.equal(dat$x2_ar, expected_x2_ar))) {
            stop("Unexpected x2_ar: ", dat$x2_ar, "; expected ", expected_x2_ar, call. = FALSE)
        }
        if (!is.null(dat$x2_innov_sd) && !isTRUE(all.equal(dat$x2_innov_sd, expected_x2_innov_sd))) {
            stop("Unexpected x2_innov_sd: ", dat$x2_innov_sd,
                 "; expected ", expected_x2_innov_sd, call. = FALSE)
        }
        if (!is.null(dat$stress_type) && !identical(dat$stress_type, "observation_sparse_counts")) {
            stop("Unexpected stress_type: ", dat$stress_type,
                 "; expected observation_sparse_counts", call. = FALSE)
        }
    }

    if (exists("validate_s4a_continuous_x2_design", mode = "function", inherits = TRUE)) {
        validate_s4a_continuous_x2_design(
            dat,
            expected_mode = expected_x2_mode,
            expected_ar = expected_x2_ar,
            expected_innov_sd = expected_x2_innov_sd,
            strict = strict
        )
    }

    if (exists("validate_s4a_observation_data", mode = "function", inherits = TRUE)) {
        validate_s4a_observation_data(dat)
    }

    invisible(TRUE)
}

.s4a_fit_extract_draw_count <- function(fit) {
    if (is.null(fit)) return(NA_integer_)
    draw_containers <- list(
        fit$draws,
        fit$samples,
        fit$mcmc,
        fit$posterior_draws
    )
    for (obj in draw_containers) {
        if (is.null(obj)) next
        if (is.data.frame(obj) || is.matrix(obj)) return(nrow(obj))
        if (is.list(obj) && length(obj) > 0L) {
            first <- obj[[1L]]
            d <- dim(first)
            if (!is.null(d) && length(d) >= 1L) return(d[1L])
            if (is.atomic(first)) return(length(first))
        }
    }
    NA_integer_
}

.s4a_fit_extract_fit_status <- function(fit) {
    if (is.null(fit)) return("completed_no_return_object")
    for (nm in c("fit_status", "status", "convergence_status")) {
        if (!is.null(fit[[nm]]) && length(fit[[nm]]) >= 1L) {
            return(as.character(fit[[nm]][1L]))
        }
    }
    "completed"
}

.s4a_fit_one_rep <- function(rep_id,
                             data_file,
                             fit_file,
                             fit_fun,
                             gamma_fixed = 0.80,
                             delta_fixed = NULL,
                             n_iter = 12000L,
                             burn = 2000L,
                             thin = 10L,
                             seed = 440500L,
                             scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
                             output_dir = dirname(fit_file),
                             strict_data_check = TRUE,
                             save_failed_fit = TRUE,
                             verbose = TRUE,
                             extra_fit_args = list()) {
    dat <- readRDS(data_file)

    .s4a_fit_validate_data(
        dat,
        expected_scenario_id = scenario_id,
        expected_TT = 100L,
        expected_n1 = 9L,
        expected_x2_mode = "continuous_time",
        expected_x2_ar = 0.50,
        expected_x2_innov_sd = 0.80,
        strict = strict_data_check
    )

    ## Stamp fitting metadata onto the in-memory object.  This does not change
    ## the saved data file; it only makes the fixed-gamma setting explicit for
    ## engines that inspect the data object.
    dat$fit_scenario_id <- scenario_id
    dat$gamma_fixed <- gamma_fixed
    dat$gamma_is_fixed <- TRUE
    dat$learn_gamma <- FALSE
    dat$estimate_gamma <- FALSE
    dat$update_gamma <- FALSE

    rep_seed <- as.integer(seed) + as.integer(rep_id)
    set.seed(rep_seed)

    args_core <- list(
        dat = dat,
        data = dat,
        dataset = dat,
        data_obj = dat,
        data_file = data_file,
        input_file = data_file,
        rep = as.integer(rep_id),
        rep_id = as.integer(rep_id),
        replicate = as.integer(rep_id),
        scenario_id = scenario_id,
        scenario_name = scenario_id,
        scenario = scenario_id,
        gamma_fixed = gamma_fixed,
        fixed_gamma = gamma_fixed,
        gamma_truth = gamma_fixed,
        gamma_common = gamma_fixed,
        fix_gamma = TRUE,
        learn_gamma = FALSE,
        estimate_gamma = FALSE,
        update_gamma = FALSE,
        delta_fixed = delta_fixed,
        fixed_delta = delta_fixed,
        learn_delta = FALSE,
        estimate_delta = FALSE,
        update_delta = FALSE,
        use_fine_counts = FALSE,
        use_omega = FALSE,
        fit_omega = FALSE,
        n_iter = as.integer(n_iter),
        iter = as.integer(n_iter),
        n_mcmc = as.integer(n_iter),
        burn = as.integer(burn),
        n_burn = as.integer(burn),
        burn_in = as.integer(burn),
        thin = as.integer(thin),
        seed = rep_seed,
        fit_seed = rep_seed,
        output_dir = output_dir,
        out_dir = output_dir,
        fit_file = fit_file,
        out_file = fit_file,
        output_file = fit_file,
        save_path = fit_file,
        verbose = verbose
    )

    args_full <- utils::modifyList(args_core, extra_fit_args)
    args_conservative <- utils::modifyList(list(
        dat = dat,
        data = dat,
        data_file = data_file,
        rep = as.integer(rep_id),
        rep_id = as.integer(rep_id),
        gamma_fixed = gamma_fixed,
        n_iter = as.integer(n_iter),
        burn = as.integer(burn),
        thin = as.integer(thin),
        seed = rep_seed,
        output_dir = output_dir,
        fit_file = fit_file,
        verbose = verbose
    ), extra_fit_args)

    t0 <- proc.time()[["elapsed"]]
    fit <- .s4a_fit_call_with_supported_args(
        fun = fit_fun,
        args_full = args_full,
        args_conservative = args_conservative
    )
    elapsed <- proc.time()[["elapsed"]] - t0

    if (is.null(fit) && file.exists(fit_file)) {
        fit <- readRDS(fit_file)
    }

    fit_record <- list(
        scenario_id = scenario_id,
        rep_id = as.integer(rep_id),
        data_file = data_file,
        fit_file = fit_file,
        gamma_fixed = gamma_fixed,
        n_iter = as.integer(n_iter),
        burn = as.integer(burn),
        thin = as.integer(thin),
        seed = rep_seed,
        elapsed_sec = elapsed,
        engine_name = NA_character_,
        fit_status = .s4a_fit_extract_fit_status(fit),
        n_saved_draws = .s4a_fit_extract_draw_count(fit),
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )

    if (is.list(fit)) {
        fit$run_metadata <- utils::modifyList(fit$run_metadata %||% list(), fit_record)
    } else {
        fit <- list(
            fit_object = fit,
            run_metadata = fit_record
        )
    }

    if (!file.exists(fit_file) || !is.null(fit)) {
        saveRDS(fit, fit_file)
    }

    cbind(
        .s4a_fit_summarise_data(dat, rep_id = rep_id, data_file = data_file),
        data.frame(
            fit_file = fit_file,
            gamma_fixed = gamma_fixed,
            n_iter = as.integer(n_iter),
            burn = as.integer(burn),
            thin = as.integer(thin),
            seed = rep_seed,
            elapsed_sec = elapsed,
            fit_status = fit_record$fit_status,
            n_saved_draws = fit_record$n_saved_draws,
            stringsAsFactors = FALSE
        )
    )
}

run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100 <- function(
    root = ".",
    reps = 1:10,
    scenario_id = "S4A_SPARSE_COUNTS_CONTINUOUS_X2_T100",
    data_dir = file.path(root, "data_s4a_sparse_counts_continuous_x2", scenario_id),
    output_dir = file.path(root, "output_s4a_sparse_counts_continuous_x2", scenario_id),
    fit_subdir = "fits",
    tables_subdir = "tables",
    gamma_fixed = 0.80,
    delta_fixed = NULL,
    n_iter = 12000L,
    burn = 2000L,
    thin = 10L,
    seed = 440500L,
    fit_fun = NULL,
    fit_engine_script = NULL,
    s3_script = NULL,
    s4a_data_script = NULL,
    extra_fit_args = list(),
    overwrite_existing = FALSE,
    strict_data_check = TRUE,
    continue_on_error = TRUE,
    verbose = TRUE
) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
    fit_dir <- file.path(output_dir, fit_subdir)
    table_dir <- file.path(output_dir, tables_subdir)

    .s4a_fit_dir_create(output_dir)
    .s4a_fit_dir_create(fit_dir)
    .s4a_fit_dir_create(table_dir)

    .s4a_fit_msg("\n=== MSSTNB S4A fixed-gamma fitting: continuous-time x2 ===", verbose = verbose)
    .s4a_fit_msg("scenario_id: ", scenario_id, verbose = verbose)
    .s4a_fit_msg("root:        ", root, verbose = verbose)
    .s4a_fit_msg("data_dir:    ", data_dir, verbose = verbose)
    .s4a_fit_msg("output_dir:  ", output_dir, verbose = verbose)
    .s4a_fit_msg("fit_dir:     ", fit_dir, verbose = verbose)
    .s4a_fit_msg("gamma_fixed: ", gamma_fixed, verbose = verbose)
    .s4a_fit_msg("n_iter/burn/thin: ", n_iter, "/", burn, "/", thin, verbose = verbose)

    .s4a_fit_source_project_dependencies(
        root = root,
        s3_script = s3_script,
        s4a_data_script = s4a_data_script,
        fit_engine_script = fit_engine_script,
        verbose = verbose
    )

    fit_engine <- .s4a_fit_resolve_engine(fit_fun = fit_fun, verbose = verbose)
    engine_name <- deparse(substitute(fit_engine))

    manifest_list <- vector("list", length(reps))

    for (ii in seq_along(reps)) {
        rep_id <- as.integer(reps[[ii]])
        data_file <- .s4a_fit_get_rep_file(data_dir = data_dir, rep_id = rep_id)
        fit_file <- file.path(fit_dir, sprintf("fit_rep%02d.rds", rep_id))

        .s4a_fit_msg("\n--- S4A fixed-gamma fit rep ", rep_id, " ---", verbose = verbose)
        .s4a_fit_msg("data: ", data_file, verbose = verbose)
        .s4a_fit_msg("fit:  ", fit_file, verbose = verbose)

        if (file.exists(fit_file) && !isTRUE(overwrite_existing)) {
            .s4a_fit_msg("Skipping existing fit because overwrite_existing = FALSE.", verbose = verbose)
            dat <- readRDS(data_file)
            one <- cbind(
                .s4a_fit_summarise_data(dat, rep_id = rep_id, data_file = data_file),
                data.frame(
                    fit_file = fit_file,
                    gamma_fixed = gamma_fixed,
                    n_iter = as.integer(n_iter),
                    burn = as.integer(burn),
                    thin = as.integer(thin),
                    seed = as.integer(seed) + rep_id,
                    elapsed_sec = NA_real_,
                    fit_status = "skipped_existing",
                    n_saved_draws = NA_integer_,
                    stringsAsFactors = FALSE
                )
            )
            manifest_list[[ii]] <- one
            next
        }

        one <- tryCatch(
            .s4a_fit_one_rep(
                rep_id = rep_id,
                data_file = data_file,
                fit_file = fit_file,
                fit_fun = fit_engine,
                gamma_fixed = gamma_fixed,
                delta_fixed = delta_fixed,
                n_iter = n_iter,
                burn = burn,
                thin = thin,
                seed = seed,
                scenario_id = scenario_id,
                output_dir = fit_dir,
                strict_data_check = strict_data_check,
                verbose = verbose,
                extra_fit_args = extra_fit_args
            ),
            error = function(e) {
                if (!isTRUE(continue_on_error)) stop(e)
                .s4a_fit_msg("ERROR in rep ", rep_id, ": ", conditionMessage(e), verbose = TRUE)
                dat <- tryCatch(readRDS(data_file), error = function(e2) NULL)
                if (is.null(dat)) {
                    base <- data.frame(
                        rep_id = rep_id,
                        data_file = data_file,
                        scenario_id = scenario_id,
                        TT = NA_integer_,
                        n1 = NA_integer_,
                        stress_type = NA_character_,
                        x2_mode = NA_character_,
                        x2_ar = NA_real_,
                        x2_innov_sd = NA_real_,
                        sparse_beta0_shift = NA_real_,
                        beta0_truth_ident = NA_real_,
                        beta1_truth = NA_real_,
                        beta2_truth = NA_real_,
                        r_truth = NA_real_,
                        gamma_truth = NA_real_,
                        mean_count = NA_real_,
                        median_count = NA_real_,
                        zero_prop = NA_real_,
                        total_count = NA_real_,
                        max_count = NA_real_,
                        lambda_raw_median = NA_real_,
                        lambda_raw_max = NA_real_,
                        stringsAsFactors = FALSE
                    )
                } else {
                    base <- .s4a_fit_summarise_data(dat, rep_id = rep_id, data_file = data_file)
                }
                cbind(
                    base,
                    data.frame(
                        fit_file = fit_file,
                        gamma_fixed = gamma_fixed,
                        n_iter = as.integer(n_iter),
                        burn = as.integer(burn),
                        thin = as.integer(thin),
                        seed = as.integer(seed) + rep_id,
                        elapsed_sec = NA_real_,
                        fit_status = paste0("error: ", conditionMessage(e)),
                        n_saved_draws = NA_integer_,
                        stringsAsFactors = FALSE
                    )
                )
            }
        )

        ## Add the actual engine name where possible.
        one$engine <- if (!is.null(fit_fun) && is.character(fit_fun)) fit_fun else NA_character_
        manifest_list[[ii]] <- one

        manifest_so_far <- do.call(rbind, manifest_list[!vapply(manifest_list, is.null, logical(1))])
        utils::write.csv(
            manifest_so_far,
            file.path(output_dir, paste0("manifest_", scenario_id, "_fits.csv")),
            row.names = FALSE
        )
    }

    manifest <- do.call(rbind, manifest_list)
    manifest_file <- file.path(output_dir, paste0("manifest_", scenario_id, "_fits.csv"))
    utils::write.csv(manifest, manifest_file, row.names = FALSE)

    status_counts <- as.data.frame(table(manifest$fit_status), stringsAsFactors = FALSE)
    names(status_counts) <- c("fit_status", "n_reps")
    status_counts$prop <- status_counts$n_reps / sum(status_counts$n_reps)
    utils::write.csv(
        status_counts,
        file.path(table_dir, "fit_status_counts.csv"),
        row.names = FALSE
    )

    compact_summary <- data.frame(
        scenario_id = scenario_id,
        n_reps = nrow(manifest),
        gamma_fixed = gamma_fixed,
        n_completed = sum(manifest$fit_status %in% c("completed", "stable", "ok"), na.rm = TRUE),
        n_skipped_existing = sum(manifest$fit_status == "skipped_existing", na.rm = TRUE),
        n_error = sum(grepl("^error:", manifest$fit_status), na.rm = TRUE),
        mean_count_avg = mean(manifest$mean_count, na.rm = TRUE),
        zero_prop_avg = mean(manifest$zero_prop, na.rm = TRUE),
        total_count_sum = sum(manifest$total_count, na.rm = TRUE),
        median_elapsed_sec = stats::median(manifest$elapsed_sec, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    utils::write.csv(
        compact_summary,
        file.path(table_dir, "s4a_fixed_gamma_fit_compact_summary.csv"),
        row.names = FALSE
    )

    .s4a_fit_msg("\nSaved fit manifest: ", manifest_file, verbose = verbose)
    .s4a_fit_msg("Saved status counts: ", file.path(table_dir, "fit_status_counts.csv"), verbose = verbose)

    invisible(list(
        manifest = manifest,
        status_counts = status_counts,
        compact_summary = compact_summary,
        output_dir = output_dir,
        fit_dir = fit_dir,
        table_dir = table_dir
    ))
}

## Convenience alias with a shorter name.
run_s4a_fixed_gamma_continuous_x2_T100 <- run_s4a_sparse_counts_fixed_gamma_continuous_x2_T100

## End of file.
