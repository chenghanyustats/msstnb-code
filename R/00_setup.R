## ==========================================================================
## 00_setup.R
## Packages, spatial graph, and shared constants for the MSSTNB simulation
##
## DESIGN: All user-adjustable values live in CONFIG. Everything below is
## derived from CONFIG and should not be edited directly.
## ==========================================================================


## ======================================================================
## SINGLE SOURCE OF TRUTH — edit here only
## ======================================================================
CONFIG <- list(

    ## ---- simulation scale ----
    R_TOTAL     = 50L,       # replications per scenario
    MASTER_SEED = 2026L,     # master RNG seed

    ## ---- spatial grid ----
    NROW_GRID   = 3L,        # grid rows    → N1 = NROW * NCOL
    NCOL_GRID   = 3L,        # grid columns

    ## ---- tree ----
    N_CHILDREN  = 3L,        # children per coarsest node
    L_TREE      = 2L,        # tree depth

    ## ---- temporal ----
    TT          = 60L,       # number of time points

    ## ---- covariates ----
    P_COV       = 2L,        # number of covariates

    ## ---- filter initial values ----
    A0          = 1.0,       # Gamma filter initial shape
    B0          = 1.0,       # Gamma filter initial rate (prior mean λ̃ = 1)
    C0          = c(5, 3, 2) # Dirichlet filter initial concentration (length = N_CHILDREN)
)


## ======================================================================
## PACKAGES
## ======================================================================
library(Matrix)
if (requireNamespace("MCMCpack", quietly = TRUE)) library(MCMCpack)


## ======================================================================
## DERIVED QUANTITIES — do not edit; change CONFIG above instead
## ======================================================================

## ---- unpack CONFIG into global constants (for convenience) ----
R_TOTAL     <- CONFIG$R_TOTAL
MASTER_SEED <- CONFIG$MASTER_SEED
NROW_GRID   <- CONFIG$NROW_GRID
NCOL_GRID   <- CONFIG$NCOL_GRID
N_CHILDREN  <- CONFIG$N_CHILDREN
L_TREE      <- CONFIG$L_TREE
TT          <- CONFIG$TT
P_COV       <- CONFIG$P_COV
A0          <- CONFIG$A0
B0          <- CONFIG$B0
C0          <- CONFIG$C0

## ---- derived spatial ----
N1 <- NROW_GRID * NCOL_GRID              # number of coarsest regions
N2 <- N1 * N_CHILDREN                    # number of level-2 nodes

COORDS <- expand.grid(row = 1:NROW_GRID, col = 1:NCOL_GRID)

# Queen adjacency
W <- matrix(0L, N1, N1)
for (i in 1:N1) {
    for (j in 1:N1) {
        if (i != j &&
            abs(COORDS$row[i] - COORDS$row[j]) <= 1 &&
            abs(COORDS$col[i] - COORDS$col[j]) <= 1) {
            W[i, j] <- 1L
        }
    }
}

# Graph Laplacian
H <- diag(rowSums(W)) - W
H_sparse <- Matrix(H, sparse = TRUE)

# Eigen-decomposition (precompute; used by generate_icar)
H_EIG  <- eigen(H, symmetric = TRUE)
B_ICAR <- H_EIG$vectors[, 1:(N1 - 1), drop = FALSE]

## ---- derived tree ----
CHILDREN <- lapply(1:N1, function(j) {
    (N_CHILDREN * (j - 1) + 1):(N_CHILDREN * j)
})

## ---- derived seeds ----
set.seed(MASTER_SEED)
REP_SEEDS <- sample.int(1e6, size = R_TOTAL)


## ======================================================================
## AUTO-EXPORT: write CONFIG to dictionaries/setup_config.csv
## ======================================================================

#' Write dictionaries/setup_config.csv from CONFIG
#'
#' Call after changing any value in CONFIG.
#'
#' @param out_dir Directory for the CSV (default "dictionaries")
export_setup_csv <- function(out_dir = "dictionary") {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    rows <- lapply(names(CONFIG), function(nm) {
        val <- CONFIG[[nm]]
        data.frame(
            name        = nm,
            value       = paste(val, collapse = ", "),
            type        = class(val)[1],
            length      = length(val),
            description = switch(nm,
                R_TOTAL     = "Replications per scenario",
                MASTER_SEED = "Master RNG seed for reproducibility",
                NROW_GRID   = "Grid rows (coarsest level)",
                NCOL_GRID   = "Grid columns (coarsest level)",
                N_CHILDREN  = "Children per coarsest node",
                L_TREE      = "Tree depth (number of levels)",
                TT          = "Number of time points",
                P_COV       = "Number of covariates",
                A0          = "Gamma filter initial shape",
                B0          = "Gamma filter initial rate",
                C0          = "Dirichlet filter initial concentration",
                ""
            ),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)

    # Append derived quantities
    derived <- data.frame(
        name        = c("N1", "N2"),
        value       = c(N1, N2),
        type        = c("integer", "integer"),
        length      = c(1, 1),
        description = c(
            paste0("Coarsest regions (", NROW_GRID, " x ", NCOL_GRID, ")"),
            paste0("Level-2 nodes (", N1, " x ", N_CHILDREN, ")")
        ),
        stringsAsFactors = FALSE
    )
    df <- rbind(df, derived)

    out_path <- file.path(out_dir, "setup_config.csv")
    write.csv(df, out_path, row.names = FALSE)
    cat("Wrote", out_path, "\n")
    invisible(df)
}


## ======================================================================
## Print summary
## ======================================================================
cat(sprintf("Setup complete: N1 = %d, N2 = %d, T = %d, R = %d\n",
            N1, N2, TT, R_TOTAL))
