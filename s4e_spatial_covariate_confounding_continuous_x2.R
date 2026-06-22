# ============================================================================
# MSSTNB Scenario 4E: strong spatial/covariate confounding stress test
# with continuous-time x2
# ============================================================================
# Purpose:
#   Generate Scenario 4E data using the same continuous-time x2 covariate design
#   used in S1--S3 and the revised S4A--S4D workflows, but deliberately align
#   the spatial pattern of x2 with the spatial random effect phi.
#
# Stress mechanism (Option A):
#   x2[t, j] = continuous-time AR(1) temporal covariate + spatial component
#              aligned with standardized phi[j].
#
#   The resulting x2 remains continuous and temporally varying, but its area-
#   mean pattern is strongly correlated with phi.  This isolates spatial/
#   covariate confounding from sparse counts, low exposure, strong
#   overdispersion, or short time series.
#
# Official revised design defaults:
#   x2_mode     = "continuous_time"
#   x2_ar       = 0.50
#   x2_innov_sd = 0.80
#   T           = 100
#   r           = 15
#   gamma       = 0.8 in fitting
#
# This script only generates and calibrates data; it does not fit the model.
# ============================================================================

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.same_dim_s4e <- function(x, target_dim) {
    d <- dim(x)
    !is.null(d) && length(d) == length(target_dim) &&
        all(as.integer(d) == as.integer(target_dim))
}

.cv_s4e <- function(x) {
    x <- as.numeric(x)
    mx <- mean(x)
    if (length(x) == 0L || !is.finite(mx) || abs(mx) < .Machine$double.eps) return(NA_real_)
    stats::sd(x) / mx
}

.safe_ratio_s4e <- function(num, den) {
    if (!is.finite(num) || !is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    num / den
}

.safe_cor_s4e <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    if (stats::sd(x[ok]) <= .Machine$double.eps || stats::sd(y[ok]) <= .Machine$double.eps) {
        return(NA_real_)
    }
    as.numeric(stats::cor(x[ok], y[ok]))
}

.standardize_s4e <- function(x) {
    x <- as.numeric(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= .Machine$double.eps) {
        stop("Cannot standardize vector with zero or non-finite SD.", call. = FALSE)
    }
    (x - mean(x)) / s
}

.require_file_s4e <- function(path) {
    if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
    invisible(path)
}

.source_checked_s4e <- function(path, verbose = TRUE) {
    .require_file_s4e(path)
    if (isTRUE(verbose)) message("source: ", path)
    source(path, local = .GlobalEnv)
    invisible(TRUE)
}

.find_s3_script_s4e <- function(root = ".", s3_script = NULL) {
    if (!is.null(s3_script)) return(s3_script)

    candidates <- c(
        file.path(root, "s3_dynamic_learned_gamma.R"),
        file.path(root, "scripts", "s3_dynamic_learned_gamma.R"),
        file.path(root, "R", "s3_dynamic_learned_gamma.R"),
        file.path(root, "R", "scenarios", "s3_dynamic_learned_gamma.R"),
        file.path(root, "scenarios", "s3_dynamic_learned_gamma.R")
    )
    hits <- candidates[file.exists(candidates)]
    if (length(hits) == 0L) {
        stop(
            "Could not find s3_dynamic_learned_gamma.R. ",
            "Pass its location with s3_script = 'path/to/s3_dynamic_learned_gamma.R'.",
            call. = FALSE
        )
    }
    hits[1L]
}

.require_object_s4e <- function(name) {
    if (!exists(name, envir = .GlobalEnv)) {
        stop("Required object not found after sourcing Scenario 3 script: ", name,
             call. = FALSE)
    }
    invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Source Scenario 3 and project dependencies
# -----------------------------------------------------------------------------
source_s4e_spatial_covariate_confounding <- function(root = ".",
                                                     s3_script = NULL,
                                                     verbose = TRUE) {
    s3_path <- .find_s3_script_s4e(root = root, s3_script = s3_script)
    .source_checked_s4e(s3_path, verbose = verbose)

    .require_object_s4e("source_s3_dynamic_learned_gamma")
    source_s3_dynamic_learned_gamma(root = root, verbose = verbose)

    needed <- c(
        "simulate_s3_dynamic_learned_gamma_one",
        "validate_s3_data",
        "compute_s3_xi",
        "recenter",
        "REP_SEEDS",
        "TT", "N1", "N_CHILDREN"
    )
    missing <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
    if (length(missing) > 0L) {
        stop("After sourcing Scenario 3 dependencies, missing objects: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    if (isTRUE(verbose)) {
        message("Scenario 4E spatial/covariate confounding data generator loaded.")
    }
    invisible(TRUE)
}

.call_s3_generator_s4e <- function(...) {
    .require_object_s4e("simulate_s3_dynamic_learned_gamma_one")
    fun <- get("simulate_s3_dynamic_learned_gamma_one", envir = .GlobalEnv)
    args <- list(...)
    fml <- names(formals(fun))

    if ("..." %in% fml) {
        return(do.call(fun, args))
    }

    args_to_pass <- args[names(args) %in% fml]
    do.call(fun, args_to_pass)
}

# -----------------------------------------------------------------------------
# Basic summaries
# -----------------------------------------------------------------------------
.count_stats_s4e <- function(y, prefix = "") {
    yy <- as.numeric(y)
    mn <- mean(yy)
    vv <- stats::var(yy)
    qs <- stats::quantile(
        yy,
        probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
        names = FALSE,
        type = 7
    )
    out <- data.frame(
        mean_count = mn,
        median_count = stats::median(yy),
        zero_prop = mean(yy == 0),
        nonzero_prop = mean(yy > 0),
        total_count = sum(yy),
        max_count = max(yy),
        count_sd = stats::sd(yy),
        count_var = vv,
        count_cv = .cv_s4e(yy),
        variance_to_mean = .safe_ratio_s4e(vv, mn),
        max_to_mean = .safe_ratio_s4e(max(yy), mn),
        q05_count = as.numeric(qs[1]),
        q25_count = as.numeric(qs[2]),
        q50_count = as.numeric(qs[3]),
        q75_count = as.numeric(qs[4]),
        q90_count = as.numeric(qs[5]),
        q95_count = as.numeric(qs[6]),
        q99_count = as.numeric(qs[7]),
        stringsAsFactors = FALSE
    )
    if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
    out
}

.summarise_s4e_x2_design <- function(dat) {
    x2 <- dat$x2
    if (is.null(x2)) {
        return(list(
            x2_found = FALSE,
            x2_mean = NA_real_,
            x2_sd = NA_real_,
            x2_min = NA_real_,
            x2_max = NA_real_,
            x2_binary_like_prop = NA_real_,
            x2_empirical_ar1_mean = NA_real_,
            x2_empirical_ar1_median = NA_real_
        ))
    }

    x2_mat <- as.matrix(x2)
    ar_vals <- rep(NA_real_, ncol(x2_mat))
    if (nrow(x2_mat) >= 3L) {
        for (jj in seq_len(ncol(x2_mat))) {
            ar_vals[jj] <- .safe_cor_s4e(
                x2_mat[-1L, jj],
                x2_mat[-nrow(x2_mat), jj]
            )
        }
    }

    near_zero_one <- abs(as.numeric(x2_mat)) < 1e-10 |
        abs(as.numeric(x2_mat) - 1) < 1e-10

    list(
        x2_found = TRUE,
        x2_mean = mean(x2_mat),
        x2_sd = stats::sd(as.numeric(x2_mat)),
        x2_min = min(x2_mat),
        x2_max = max(x2_mat),
        x2_binary_like_prop = mean(near_zero_one),
        x2_empirical_ar1_mean = mean(ar_vals, na.rm = TRUE),
        x2_empirical_ar1_median = stats::median(ar_vals, na.rm = TRUE)
    )
}

.get_phi_for_confounding_s4e <- function(dat) {
    if (!is.null(dat$phi_star_ident)) return(as.numeric(dat$phi_star_ident))
    .require_object_s4e("recenter")
    rc <- recenter(
        beta0 = dat$beta0_star,
        phi = dat$phi_star,
        lambda_tilde = dat$lambda_tilde,
        return_diag = TRUE
    )
    as.numeric(rc$phi)
}

.make_confounded_x2_s4e <- function(base_x2,
                                    phi,
                                    confounding_strength = 1.25,
                                    residual_weight = 1.00,
                                    target_sd = 1.00) {
    base <- as.matrix(base_x2)
    TT_now <- nrow(base)
    n1_now <- ncol(base)

    if (length(phi) != n1_now) {
        stop("Length of phi does not match the number of coarse areas.", call. = FALSE)
    }
    if (!is.finite(confounding_strength) || confounding_strength <= 0) {
        stop("confounding_strength must be positive and finite.", call. = FALSE)
    }
    if (!is.finite(residual_weight) || residual_weight < 0) {
        stop("residual_weight must be nonnegative and finite.", call. = FALSE)
    }
    if (!is.finite(target_sd) || target_sd <= 0) {
        stop("target_sd must be positive and finite.", call. = FALSE)
    }

    base_z <- matrix(.standardize_s4e(base), nrow = TT_now, ncol = n1_now)
    phi_z <- .standardize_s4e(phi)
    spatial_component <- matrix(rep(phi_z, each = TT_now), nrow = TT_now, ncol = n1_now)

    raw <- residual_weight * base_z + confounding_strength * spatial_component
    out <- matrix(.standardize_s4e(raw) * target_sd, nrow = TT_now, ncol = n1_now)
    colnames(out) <- colnames(base)
    rownames(out) <- rownames(base)
    out
}

.confounding_stats_s4e <- function(x2, phi) {
    x2_mat <- as.matrix(x2)
    TT_now <- nrow(x2_mat)
    n1_now <- ncol(x2_mat)
    phi <- as.numeric(phi)
    phi_z <- .standardize_s4e(phi)

    if (length(phi_z) != n1_now) {
        stop("Length of phi does not match x2 columns in confounding summary.", call. = FALSE)
    }

    phi_cell <- rep(phi_z, each = TT_now)
    x2_vec <- as.numeric(x2_mat)
    x2_area_mean <- colMeans(x2_mat)
    x2_area_sd <- apply(x2_mat, 2, stats::sd)

    data.frame(
        x2_phi_cell_cor = .safe_cor_s4e(x2_vec, phi_cell),
        x2_phi_abs_cell_cor = abs(.safe_cor_s4e(x2_vec, phi_cell)),
        x2_phi_area_mean_cor = .safe_cor_s4e(x2_area_mean, phi_z),
        x2_phi_abs_area_mean_cor = abs(.safe_cor_s4e(x2_area_mean, phi_z)),
        x2_area_mean_sd = stats::sd(x2_area_mean),
        x2_area_sd_mean = mean(x2_area_sd),
        x2_area_sd_min = min(x2_area_sd),
        x2_area_sd_max = max(x2_area_sd),
        stringsAsFactors = FALSE
    )
}

# -----------------------------------------------------------------------------
# Data validation
# -----------------------------------------------------------------------------
validate_s4e_confounding_data <- function(dat,
                                          expected_scenario_prefix = "S4E",
                                          require_coherence = TRUE) {
    required <- c(
        "scenario_id", "rep_id", "TT", "n1", "n_children",
        "y_coarse", "y_fine", "x1", "x2", "e",
        "xi", "mu_nb", "poisson_rate", "lambda_tilde", "kappa",
        "omega", "beta0_star", "beta_star", "phi_star", "r_star",
        "x2_mode", "x2_ar", "x2_innov_sd",
        "x2_phi_cell_cor", "x2_phi_area_mean_cor"
    )
    missing <- required[!vapply(required, function(z) !is.null(dat[[z]]), logical(1))]
    if (length(missing) > 0L) {
        stop("S4E data is missing fields: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    if (!startsWith(dat$scenario_id, expected_scenario_prefix)) {
        stop("scenario_id should start with ", expected_scenario_prefix, ".", call. = FALSE)
    }

    TT_now <- as.integer(dat$TT)
    n1_now <- as.integer(dat$n1)
    n_children_now <- as.integer(dat$n_children)

    if (!.same_dim_s4e(dat$y_coarse, c(TT_now, n1_now))) {
        stop("y_coarse has wrong dimension.", call. = FALSE)
    }
    if (!.same_dim_s4e(dat$y_fine, c(TT_now, n1_now, n_children_now))) {
        stop("y_fine has wrong dimension.", call. = FALSE)
    }
    for (nm in c("x1", "x2", "e", "xi", "mu_nb", "poisson_rate", "lambda_tilde", "kappa")) {
        if (!.same_dim_s4e(dat[[nm]], c(TT_now, n1_now))) {
            stop(nm, " has wrong dimension.", call. = FALSE)
        }
    }
    if (!.same_dim_s4e(dat$omega, c(TT_now, n1_now, n_children_now))) {
        stop("omega has wrong dimension.", call. = FALSE)
    }

    if (isTRUE(require_coherence)) {
        coherent <- all(dat$y_coarse == apply(dat$y_fine, c(1, 2), sum))
        if (!coherent) stop("Fine counts are not coherent with coarse counts.", call. = FALSE)
    }

    if (!identical(dat$x2_mode, "continuous_time")) {
        stop("S4E revised data requires x2_mode = 'continuous_time'.", call. = FALSE)
    }
    if (!isTRUE(all.equal(dat$x2_ar, 0.50))) {
        stop("S4E revised data requires x2_ar = 0.50.", call. = FALSE)
    }
    if (!isTRUE(all.equal(dat$x2_innov_sd, 0.80))) {
        stop("S4E revised data requires x2_innov_sd = 0.80.", call. = FALSE)
    }

    xs <- .summarise_s4e_x2_design(dat)
    if (!isTRUE(xs$x2_found) || !is.finite(xs$x2_sd) || xs$x2_sd < 0.10) {
        stop("x2 is missing or not sufficiently variable.", call. = FALSE)
    }
    if (is.finite(xs$x2_binary_like_prop) && xs$x2_binary_like_prop > 0.25) {
        stop("x2 appears too close to a binary/indicator covariate.", call. = FALSE)
    }

    invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Simulate one S4E replicate
# -----------------------------------------------------------------------------
simulate_s4e_confounding_one <- function(seed = 1L,
                                         TT_use = TT,
                                         beta0_reference_truth = -1.5,
                                         r_truth = 15,
                                         confounding_strength = 1.25,
                                         residual_weight = 1.00,
                                         preserve_reference_mean = TRUE,
                                         scenario_id = "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_CONTINUOUS_X2_T100",
                                         rep_id = NA_integer_,
                                         max_poisson_rate = 1e9,
                                         x2_mode = "continuous_time",
                                         x2_ar = 0.50,
                                         x2_innov_sd = 0.80,
                                         ...) {
    .require_object_s4e("compute_s3_xi")
    .require_object_s4e("recenter")

    if (!identical(x2_mode, "continuous_time")) {
        stop("For revised S4E, x2_mode must be exactly 'continuous_time'.", call. = FALSE)
    }
    if (!isTRUE(all.equal(x2_ar, 0.50))) {
        stop("For revised S4E, x2_ar must be 0.50.", call. = FALSE)
    }
    if (!isTRUE(all.equal(x2_innov_sd, 0.80))) {
        stop("For revised S4E, x2_innov_sd must be 0.80.", call. = FALSE)
    }

    ref_dat <- .call_s3_generator_s4e(
        seed = seed,
        TT_use = TT_use,
        beta0_truth = beta0_reference_truth,
        r_truth = r_truth,
        scenario_id = paste0(scenario_id, "_REFERENCE_S3_LATENT"),
        rep_id = rep_id,
        max_poisson_rate = max_poisson_rate,
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        ...
    )

    TT_now <- as.integer(ref_dat$TT)
    n1_now <- as.integer(ref_dat$n1)
    n_children_now <- as.integer(ref_dat$n_children)

    phi_for_confounding <- .get_phi_for_confounding_s4e(ref_dat)
    x2_reference <- ref_dat$x2
    x2_confounded <- .make_confounded_x2_s4e(
        base_x2 = x2_reference,
        phi = phi_for_confounding,
        confounding_strength = confounding_strength,
        residual_weight = residual_weight,
        target_sd = 1.00
    )

    beta0_prelim <- beta0_reference_truth
    xi_prelim <- compute_s3_xi(
        e = ref_dat$e,
        x1 = ref_dat$x1,
        x2 = x2_confounded,
        beta0 = beta0_prelim,
        beta = ref_dat$beta_star,
        phi = ref_dat$phi_star
    )
    mu_prelim <- xi_prelim * ref_dat$lambda_tilde

    target_mu_mean <- if (!is.null(ref_dat$mu_nb)) {
        mean(ref_dat$mu_nb)
    } else {
        mean(ref_dat$xi * ref_dat$lambda_tilde)
    }
    prelim_mu_mean <- mean(mu_prelim)
    beta0_shift <- if (isTRUE(preserve_reference_mean)) {
        if (!is.finite(target_mu_mean) || target_mu_mean <= 0 ||
            !is.finite(prelim_mu_mean) || prelim_mu_mean <= 0) {
            stop("Cannot compute beta0 shift to preserve reference mean.", call. = FALSE)
        }
        log(target_mu_mean / prelim_mu_mean)
    } else {
        0
    }
    beta0_confounding_truth <- beta0_prelim + beta0_shift

    xi_confounded <- compute_s3_xi(
        e = ref_dat$e,
        x1 = ref_dat$x1,
        x2 = x2_confounded,
        beta0 = beta0_confounding_truth,
        beta = ref_dat$beta_star,
        phi = ref_dat$phi_star
    )
    mu_nb_confounded <- xi_confounded * ref_dat$lambda_tilde
    poisson_rate_confounded <- mu_nb_confounded * ref_dat$kappa

    bad <- !is.finite(poisson_rate_confounded) | poisson_rate_confounded < 0 |
        poisson_rate_confounded > max_poisson_rate
    if (any(bad)) {
        idx <- which(bad, arr.ind = TRUE)[1, ]
        stop(sprintf(
            paste0(
                "Bad S4E Poisson rate. First bad cell: t=%d, j=%d, rate=%s, ",
                "xi=%s, lambda=%s, kappa=%s."
            ),
            idx[1], idx[2],
            as.character(poisson_rate_confounded[idx[1], idx[2]]),
            as.character(xi_confounded[idx[1], idx[2]]),
            as.character(ref_dat$lambda_tilde[idx[1], idx[2]]),
            as.character(ref_dat$kappa[idx[1], idx[2]])
        ), call. = FALSE)
    }

    set.seed(as.integer(seed) + 888888L + as.integer(round(1000 * confounding_strength)))
    y_coarse_confounded <- matrix(
        stats::rpois(TT_now * n1_now, lambda = as.numeric(poisson_rate_confounded)),
        nrow = TT_now,
        ncol = n1_now
    )
    if (anyNA(y_coarse_confounded)) stop("rpois produced NA values in S4E.", call. = FALSE)

    y_fine_confounded <- array(0L, dim = c(TT_now, n1_now, n_children_now))
    for (t in seq_len(TT_now)) {
        for (j in seq_len(n1_now)) {
            if (y_coarse_confounded[t, j] > 0L) {
                y_fine_confounded[t, j, ] <- as.integer(stats::rmultinom(
                    1L,
                    size = y_coarse_confounded[t, j],
                    prob = ref_dat$omega[t, j, ]
                ))
            }
        }
    }

    coherent <- all(y_coarse_confounded == apply(y_fine_confounded, c(1, 2), sum))
    if (!coherent) stop("S4E fine counts are not coherent.", call. = FALSE)

    rc_truth <- recenter(
        beta0 = beta0_confounding_truth,
        phi = ref_dat$phi_star,
        lambda_tilde = ref_dat$lambda_tilde,
        return_diag = TRUE
    )

    conf_stats <- .confounding_stats_s4e(x2_confounded, rc_truth$phi)
    ref_conf_stats <- .confounding_stats_s4e(x2_reference, rc_truth$phi)
    x2_stats <- .summarise_s4e_x2_design(list(x2 = x2_confounded))

    dat <- ref_dat
    dat$scenario_id <- scenario_id
    dat$reference_scenario_id <- "S3_DYNAMIC_LEARNED_GAMMA"
    dat$data_type <- "dynamic_lambda_learned_gamma_spatial_covariate_confounding"
    dat$stress_type <- "spatial_covariate_confounding"
    dat$confounding_type <- "x2_area_mean_aligned_with_phi"
    dat$stress_description <- paste0(
        "Scenario 3 latent structure with continuous-time x2 replaced by a ",
        "temporally varying covariate whose spatial mean is strongly aligned ",
        "with phi. This isolates spatial/covariate confounding."
    )

    dat$y_coarse_reference <- ref_dat$y_coarse
    dat$y_fine_reference <- ref_dat$y_fine
    dat$x2_reference <- x2_reference
    dat$xi_reference <- ref_dat$xi
    dat$mu_nb_reference <- ref_dat$mu_nb
    dat$poisson_rate_reference <- ref_dat$poisson_rate
    dat$reference_count_summary <- as.list(.count_stats_s4e(ref_dat$y_coarse, prefix = ""))

    dat$y_coarse <- y_coarse_confounded
    dat$y_fine <- y_fine_confounded
    dat$y_levels <- list(y_coarse_confounded, y_fine_confounded)
    dat$x2 <- x2_confounded
    dat$xi <- xi_confounded
    dat$mu_nb <- mu_nb_confounded
    dat$poisson_rate <- poisson_rate_confounded

    dat$beta0_reference_truth <- beta0_reference_truth
    dat$beta0_confounding_truth <- beta0_confounding_truth
    dat$beta0_shift_to_preserve_reference_mean <- beta0_shift
    dat$preserve_reference_mean <- isTRUE(preserve_reference_mean)
    dat$target_mu_mean_reference <- target_mu_mean
    dat$mu_mean_pre_beta0_shift <- prelim_mu_mean
    dat$mu_mean_after_beta0_shift <- mean(mu_nb_confounded)

    dat$beta0_star <- beta0_confounding_truth
    dat$beta_star <- ref_dat$beta_star
    dat$phi_star <- ref_dat$phi_star
    dat$beta0_star_ident <- rc_truth$beta0
    dat$beta_star_ident <- ref_dat$beta_star
    dat$phi_star_ident <- rc_truth$phi
    dat$lambda_tilde_ident <- rc_truth$lambda_tilde
    dat$lambda_recenter_diag <- rc_truth$diag

    dat$r_reference_truth <- as.numeric(r_truth)
    dat$r_star <- rep(as.numeric(r_truth), n1_now)
    dat$x2_mode <- x2_mode
    dat$x2_ar <- x2_ar
    dat$x2_innov_sd <- x2_innov_sd
    dat$covariate_design <- list(
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        description = "continuous-time AR(1) x2 with spatial mean aligned with phi"
    )

    dat$x2_confounding_strength <- as.numeric(confounding_strength)
    dat$x2_residual_weight <- as.numeric(residual_weight)
    dat$x2_phi_cell_cor <- conf_stats$x2_phi_cell_cor
    dat$x2_phi_abs_cell_cor <- conf_stats$x2_phi_abs_cell_cor
    dat$x2_phi_area_mean_cor <- conf_stats$x2_phi_area_mean_cor
    dat$x2_phi_abs_area_mean_cor <- conf_stats$x2_phi_abs_area_mean_cor
    dat$x2_area_mean_sd <- conf_stats$x2_area_mean_sd
    dat$x2_area_sd_mean <- conf_stats$x2_area_sd_mean
    dat$x2_area_sd_min <- conf_stats$x2_area_sd_min
    dat$x2_area_sd_max <- conf_stats$x2_area_sd_max
    dat$reference_x2_phi_cell_cor <- ref_conf_stats$x2_phi_cell_cor
    dat$reference_x2_phi_area_mean_cor <- ref_conf_stats$x2_phi_area_mean_cor
    dat$reference_x2_phi_abs_cell_cor <- ref_conf_stats$x2_phi_abs_cell_cor
    dat$reference_x2_phi_abs_area_mean_cor <- ref_conf_stats$x2_phi_abs_area_mean_cor

    dat$x2_mean <- x2_stats$x2_mean
    dat$x2_sd <- x2_stats$x2_sd
    dat$x2_min <- x2_stats$x2_min
    dat$x2_max <- x2_stats$x2_max
    dat$x2_binary_like_prop <- x2_stats$x2_binary_like_prop
    dat$x2_empirical_ar1_mean <- x2_stats$x2_empirical_ar1_mean
    dat$x2_empirical_ar1_median <- x2_stats$x2_empirical_ar1_median

    dat$s4e_continuous_x2_config <- list(
        scenario_id = scenario_id,
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        confounding_strength = confounding_strength,
        residual_weight = residual_weight,
        preserve_reference_mean = preserve_reference_mean,
        beta0_reference_truth = beta0_reference_truth,
        beta0_confounding_truth = beta0_confounding_truth
    )
    dat$scenario_info <- c(dat$scenario_info %||% list(), list(
        scenario_family = "S4E",
        covariate_design = "continuous_time_x2_spatially_confounded",
        stress_test = "spatial_covariate_confounding"
    ))

    dat$mean_count <- mean(y_coarse_confounded)
    dat$median_count <- stats::median(as.numeric(y_coarse_confounded))
    dat$zero_prop <- mean(y_coarse_confounded == 0)
    dat$nonzero_prop <- mean(y_coarse_confounded > 0)
    dat$total_count <- sum(y_coarse_confounded)
    dat$max_count <- max(y_coarse_confounded)
    dat$count_quantiles <- stats::quantile(
        as.numeric(y_coarse_confounded),
        probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1),
        names = TRUE
    )
    dat$count_sd <- stats::sd(as.numeric(y_coarse_confounded))
    dat$count_cv <- .cv_s4e(y_coarse_confounded)
    dat$variance_to_mean <- .safe_ratio_s4e(
        stats::var(as.numeric(y_coarse_confounded)),
        mean(y_coarse_confounded)
    )
    dat$coherent <- coherent

    validate_s4e_confounding_data(dat)
    dat
}

# -----------------------------------------------------------------------------
# Summarise one replicate
# -----------------------------------------------------------------------------
summarise_s4e_confounding_one <- function(dat) {
    validate_s4e_confounding_data(dat)

    y_stats <- .count_stats_s4e(dat$y_coarse, prefix = "")
    y_ref_stats <- .count_stats_s4e(dat$y_coarse_reference, prefix = "reference_")
    xs <- .summarise_s4e_x2_design(dat)

    main <- data.frame(
        scenario_id = dat$scenario_id,
        rep_id = as.integer(dat$rep_id),
        TT = as.integer(dat$TT),
        n1 = as.integer(dat$n1),
        stress_type = dat$stress_type,
        confounding_type = dat$confounding_type,
        x2_confounding_strength = dat$x2_confounding_strength,
        x2_residual_weight = dat$x2_residual_weight,
        preserve_reference_mean = dat$preserve_reference_mean,
        beta0_reference_truth = dat$beta0_reference_truth,
        beta0_confounding_truth = dat$beta0_confounding_truth,
        beta0_shift_to_preserve_reference_mean = dat$beta0_shift_to_preserve_reference_mean,
        beta0_star_ident = dat$beta0_star_ident,
        beta1_truth = dat$beta_star[1],
        beta2_truth = dat$beta_star[2],
        r_truth = mean(dat$r_star),
        mean_exposure = mean(dat$e),
        min_exposure = min(dat$e),
        max_exposure = max(dat$e),
        gamma_truth_mean = mean(dat$gamma_star),
        delta_truth = dat$delta_star,
        lambda_raw_min = min(dat$lambda_tilde),
        lambda_raw_median = stats::median(as.numeric(dat$lambda_tilde)),
        lambda_raw_max = max(dat$lambda_tilde),
        lambda_ident_min = min(dat$lambda_tilde_ident),
        lambda_ident_median = stats::median(as.numeric(dat$lambda_tilde_ident)),
        lambda_ident_max = max(dat$lambda_tilde_ident),
        x2_mode = dat$x2_mode,
        x2_ar = dat$x2_ar,
        x2_innov_sd = dat$x2_innov_sd,
        x2_mean = xs$x2_mean,
        x2_sd = xs$x2_sd,
        x2_min = xs$x2_min,
        x2_max = xs$x2_max,
        x2_binary_like_prop = xs$x2_binary_like_prop,
        x2_empirical_ar1_mean = xs$x2_empirical_ar1_mean,
        x2_empirical_ar1_median = xs$x2_empirical_ar1_median,
        x2_phi_cell_cor = dat$x2_phi_cell_cor,
        x2_phi_abs_cell_cor = dat$x2_phi_abs_cell_cor,
        x2_phi_area_mean_cor = dat$x2_phi_area_mean_cor,
        x2_phi_abs_area_mean_cor = dat$x2_phi_abs_area_mean_cor,
        reference_x2_phi_cell_cor = dat$reference_x2_phi_cell_cor,
        reference_x2_phi_abs_cell_cor = dat$reference_x2_phi_abs_cell_cor,
        reference_x2_phi_area_mean_cor = dat$reference_x2_phi_area_mean_cor,
        reference_x2_phi_abs_area_mean_cor = dat$reference_x2_phi_abs_area_mean_cor,
        x2_area_mean_sd = dat$x2_area_mean_sd,
        x2_area_sd_mean = dat$x2_area_sd_mean,
        x2_area_sd_min = dat$x2_area_sd_min,
        x2_area_sd_max = dat$x2_area_sd_max,
        mu_mean_after_beta0_shift = dat$mu_mean_after_beta0_shift,
        target_mu_mean_reference = dat$target_mu_mean_reference,
        stringsAsFactors = FALSE
    )
    cbind(main, y_stats, y_ref_stats)
}

# -----------------------------------------------------------------------------
# Batch generation
# -----------------------------------------------------------------------------
simulate_s4e_confounding_batch <- function(reps = 1:10,
                                           data_dir = "data_s4e_spatial_covariate_confounding_continuous_x2",
                                           scenario_id = "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_CONTINUOUS_X2_T100",
                                           root = ".",
                                           seed_base = NULL,
                                           overwrite_existing = TRUE,
                                           verbose = TRUE,
                                           manifest_name = NULL,
                                           ...) {
    out_dir <- file.path(root, data_dir, scenario_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    manifest <- list()
    for (rep_id in reps) {
        rr <- sprintf("%02d", as.integer(rep_id))
        out_file <- file.path(out_dir, paste0("data_rep", rr, ".rds"))

        if (file.exists(out_file) && !isTRUE(overwrite_existing)) {
            if (isTRUE(verbose)) message("Skipping existing file: ", out_file)
            dat <- readRDS(out_file)
            validate_s4e_confounding_data(dat)
        } else {
            seed <- if (!is.null(seed_base)) {
                as.integer(seed_base + rep_id)
            } else if (exists("REP_SEEDS", envir = .GlobalEnv) && rep_id <= length(REP_SEEDS)) {
                as.integer(REP_SEEDS[rep_id])
            } else {
                as.integer(2026000L + rep_id)
            }

            dat <- simulate_s4e_confounding_one(
                seed = seed,
                scenario_id = scenario_id,
                rep_id = as.integer(rep_id),
                ...
            )
            saveRDS(dat, out_file)
            if (isTRUE(verbose)) {
                row_tmp <- summarise_s4e_confounding_one(dat)
                message(sprintf(
                    paste0(
                        "Saved %s | mean_count=%.2f zero_prop=%.3f ",
                        "cell_cor=%.3f area_cor=%.3f x2_ar_emp=%.3f beta0_shift=%.3f"
                    ),
                    out_file,
                    row_tmp$mean_count,
                    row_tmp$zero_prop,
                    row_tmp$x2_phi_cell_cor,
                    row_tmp$x2_phi_area_mean_cor,
                    row_tmp$x2_empirical_ar1_mean,
                    row_tmp$beta0_shift_to_preserve_reference_mean
                ))
            }
        }

        row <- summarise_s4e_confounding_one(dat)
        row$file <- out_file
        manifest[[rr]] <- row
    }

    manifest_df <- do.call(rbind, manifest)
    if (is.null(manifest_name)) manifest_name <- paste0("manifest_", scenario_id, ".csv")
    manifest_file <- file.path(out_dir, manifest_name)
    write.csv(manifest_df, manifest_file, row.names = FALSE)

    if (isTRUE(verbose)) {
        message("Saved manifest: ", manifest_file)
        message(sprintf(
            paste0(
                "S4E confounding summary | reps=%d mean_count=%.2f zero_prop=%.3f ",
                "cell_cor=%.3f area_cor=%.3f x2_emp_ar=%.3f"
            ),
            nrow(manifest_df),
            mean(manifest_df$mean_count),
            mean(manifest_df$zero_prop),
            mean(manifest_df$x2_phi_cell_cor),
            mean(manifest_df$x2_phi_area_mean_cor),
            mean(manifest_df$x2_empirical_ar1_mean)
        ))
    }

    invisible(manifest_df)
}

run_s4e_spatial_covariate_confounding_continuous_x2_data_generation <- function(root = ".",
                                                                                s3_script = NULL,
                                                                                reps = 1:10,
                                                                                TT_use = 100,
                                                                                beta0_reference_truth = -1.5,
                                                                                r_truth = 15,
                                                                                confounding_strength = 1.25,
                                                                                residual_weight = 1.00,
                                                                                preserve_reference_mean = TRUE,
                                                                                data_dir = "data_s4e_spatial_covariate_confounding_continuous_x2",
                                                                                scenario_id = "S4E_STRONG_SPATIAL_COVARIATE_CONFOUNDING_CONTINUOUS_X2_T100",
                                                                                seed_base = NULL,
                                                                                overwrite_existing = TRUE,
                                                                                verbose = TRUE,
                                                                                x2_mode = "continuous_time",
                                                                                x2_ar = 0.50,
                                                                                x2_innov_sd = 0.80,
                                                                                ...) {
    source_s4e_spatial_covariate_confounding(root = root, s3_script = s3_script, verbose = verbose)

    simulate_s4e_confounding_batch(
        reps = reps,
        data_dir = data_dir,
        scenario_id = scenario_id,
        root = root,
        seed_base = seed_base,
        overwrite_existing = overwrite_existing,
        verbose = verbose,
        TT_use = TT_use,
        beta0_reference_truth = beta0_reference_truth,
        r_truth = r_truth,
        confounding_strength = confounding_strength,
        residual_weight = residual_weight,
        preserve_reference_mean = preserve_reference_mean,
        x2_mode = x2_mode,
        x2_ar = x2_ar,
        x2_innov_sd = x2_innov_sd,
        ...
    )
}

# -----------------------------------------------------------------------------
# Data-quality checks
# -----------------------------------------------------------------------------
check_s4e_spatial_covariate_confounding_continuous_x2_data_summary <- function(manifest,
                                                                               target_mean_count_range = c(80, 350),
                                                                               target_zero_prop_max = 0.10,
                                                                               minimum_abs_cell_cor = 0.50,
                                                                               minimum_abs_area_mean_cor = 0.80,
                                                                               target_abs_beta0_ident_max = 20,
                                                                               max_count_max_limit = Inf) {
    if (is.character(manifest)) manifest <- read.csv(manifest)

    required <- c(
        "rep_id", "TT", "n1", "mean_count", "median_count", "zero_prop",
        "total_count", "max_count", "count_cv", "variance_to_mean",
        "reference_mean_count", "reference_zero_prop", "beta0_star_ident",
        "beta0_confounding_truth", "beta0_shift_to_preserve_reference_mean",
        "x2_mode", "x2_ar", "x2_innov_sd", "x2_sd", "x2_binary_like_prop",
        "x2_empirical_ar1_mean", "x2_phi_cell_cor", "x2_phi_abs_cell_cor",
        "x2_phi_area_mean_cor", "x2_phi_abs_area_mean_cor",
        "reference_x2_phi_abs_cell_cor", "reference_x2_phi_abs_area_mean_cor",
        "lambda_raw_median", "lambda_raw_max"
    )
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0L) {
        stop("Manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    out <- data.frame(
        n_reps = nrow(manifest),
        TT_unique = paste(sort(unique(manifest$TT)), collapse = ","),
        n1_unique = paste(sort(unique(manifest$n1)), collapse = ","),
        confounding_strength_unique = paste(signif(sort(unique(manifest$x2_confounding_strength)), 6), collapse = ","),
        residual_weight_unique = paste(signif(sort(unique(manifest$x2_residual_weight)), 6), collapse = ","),
        preserve_reference_mean_unique = paste(sort(unique(as.character(manifest$preserve_reference_mean))), collapse = ","),
        mean_count_avg = mean(manifest$mean_count),
        mean_count_min = min(manifest$mean_count),
        mean_count_max = max(manifest$mean_count),
        reference_mean_count_avg = mean(manifest$reference_mean_count),
        median_count_avg = mean(manifest$median_count),
        zero_prop_avg = mean(manifest$zero_prop),
        zero_prop_min = min(manifest$zero_prop),
        zero_prop_max = max(manifest$zero_prop),
        reference_zero_prop_avg = mean(manifest$reference_zero_prop),
        total_count_sum = sum(manifest$total_count),
        max_count_avg = mean(manifest$max_count),
        max_count_max = max(manifest$max_count),
        count_cv_avg = mean(manifest$count_cv),
        variance_to_mean_avg = mean(manifest$variance_to_mean),
        beta0_confounding_truth_avg = mean(manifest$beta0_confounding_truth),
        beta0_shift_avg = mean(manifest$beta0_shift_to_preserve_reference_mean),
        beta0_ident_min = min(manifest$beta0_star_ident),
        beta0_ident_median = stats::median(manifest$beta0_star_ident),
        beta0_ident_max = max(manifest$beta0_star_ident),
        max_abs_beta0_ident = max(abs(manifest$beta0_star_ident)),
        x2_phi_cell_cor_avg = mean(manifest$x2_phi_cell_cor),
        x2_phi_cell_cor_min = min(manifest$x2_phi_cell_cor),
        x2_phi_cell_cor_max = max(manifest$x2_phi_cell_cor),
        x2_phi_abs_cell_cor_avg = mean(manifest$x2_phi_abs_cell_cor),
        x2_phi_abs_cell_cor_min = min(manifest$x2_phi_abs_cell_cor),
        x2_phi_area_mean_cor_avg = mean(manifest$x2_phi_area_mean_cor),
        x2_phi_area_mean_cor_min = min(manifest$x2_phi_area_mean_cor),
        x2_phi_area_mean_cor_max = max(manifest$x2_phi_area_mean_cor),
        x2_phi_abs_area_mean_cor_avg = mean(manifest$x2_phi_abs_area_mean_cor),
        x2_phi_abs_area_mean_cor_min = min(manifest$x2_phi_abs_area_mean_cor),
        reference_x2_phi_abs_cell_cor_avg = mean(manifest$reference_x2_phi_abs_cell_cor),
        reference_x2_phi_abs_area_mean_cor_avg = mean(manifest$reference_x2_phi_abs_area_mean_cor),
        x2_mode_unique = paste(sort(unique(manifest$x2_mode)), collapse = ","),
        x2_ar_unique = paste(sort(unique(manifest$x2_ar)), collapse = ","),
        x2_innov_sd_unique = paste(sort(unique(manifest$x2_innov_sd)), collapse = ","),
        x2_sd_avg = mean(manifest$x2_sd),
        x2_binary_like_prop_max = max(manifest$x2_binary_like_prop),
        x2_empirical_ar1_mean_avg = mean(manifest$x2_empirical_ar1_mean),
        lambda_raw_median_avg = mean(manifest$lambda_raw_median),
        lambda_raw_max_max = max(manifest$lambda_raw_max),
        stringsAsFactors = FALSE
    )

    out$passes_mean_count_target <-
        out$mean_count_avg >= target_mean_count_range[1] &&
        out$mean_count_avg <= target_mean_count_range[2]
    out$passes_zero_prop_target <- out$zero_prop_avg <= target_zero_prop_max
    out$passes_confounding_cell_target <- out$x2_phi_abs_cell_cor_avg >= minimum_abs_cell_cor
    out$passes_confounding_area_target <- out$x2_phi_abs_area_mean_cor_avg >= minimum_abs_area_mean_cor
    out$passes_confounding_increase_target <-
        out$x2_phi_abs_cell_cor_avg > out$reference_x2_phi_abs_cell_cor_avg &&
        out$x2_phi_abs_area_mean_cor_avg > out$reference_x2_phi_abs_area_mean_cor_avg
    out$passes_identified_scale_target <- out$max_abs_beta0_ident <= target_abs_beta0_ident_max
    out$passes_max_count_target <- out$max_count_max <= max_count_max_limit
    out$passes_continuous_x2_target <-
        identical(out$x2_mode_unique, "continuous_time") &&
        identical(out$x2_ar_unique, "0.5") &&
        identical(out$x2_innov_sd_unique, "0.8") &&
        is.finite(out$x2_sd_avg) && out$x2_sd_avg >= 0.10 &&
        is.finite(out$x2_binary_like_prop_max) && out$x2_binary_like_prop_max <= 0.25

    out$passes_s4e_data_check <-
        isTRUE(out$passes_mean_count_target) &&
        isTRUE(out$passes_zero_prop_target) &&
        isTRUE(out$passes_confounding_cell_target) &&
        isTRUE(out$passes_confounding_area_target) &&
        isTRUE(out$passes_confounding_increase_target) &&
        isTRUE(out$passes_identified_scale_target) &&
        isTRUE(out$passes_max_count_target) &&
        isTRUE(out$passes_continuous_x2_target)

    out
}

# Short aliases matching the revised workflow name.
run_s4e_continuous_x2_data_generation <- run_s4e_spatial_covariate_confounding_continuous_x2_data_generation
check_s4e_continuous_x2_data_summary <- check_s4e_spatial_covariate_confounding_continuous_x2_data_summary

# -----------------------------------------------------------------------------
# Example revised usage
# -----------------------------------------------------------------------------
# source("s4e_spatial_covariate_confounding_continuous_x2.R")
#
# manifest_test <- run_s4e_continuous_x2_data_generation(
#     root = ".",
#     reps = 1,
#     TT_use = 100,
#     confounding_strength = 1.25,
#     residual_weight = 1.00,
#     preserve_reference_mean = TRUE,
#     overwrite_existing = TRUE,
#     verbose = TRUE
# )
#
# check_test <- check_s4e_continuous_x2_data_summary(manifest_test)
# print(check_test)
#
# manifest_s4e <- run_s4e_continuous_x2_data_generation(
#     root = ".",
#     reps = 1:10,
#     TT_use = 100,
#     confounding_strength = 1.25,
#     residual_weight = 1.00,
#     preserve_reference_mean = TRUE,
#     overwrite_existing = TRUE,
#     verbose = TRUE
# )
#
# check_s4e <- check_s4e_continuous_x2_data_summary(manifest_s4e)
# print(check_s4e)
