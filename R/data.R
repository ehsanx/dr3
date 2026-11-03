## =============================================================================
## data.R -- dataset loaders for the tutorials. dr3 ships no data: each loader
## pulls an openly available dataset from another package (NHEFS via causaldata,
## RHC via Hmisc) and returns the common structure
##   list(name, Y, A, X (data.frame), family_y, outcome_type, truth)
## where `truth` is the checkable ATE (NA for the real observational datasets).
## =============================================================================

.dr_need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required for this loader. Install it with install.packages('%s').",
                 pkg, pkg), call. = FALSE)
}

.dr_impute_median <- function(df) {
  for (j in seq_along(df)) {
    x <- df[[j]]; if (anyNA(x)) df[[j]][is.na(x)] <- median(x, na.rm = TRUE)
  }
  df
}

#' Load a benchmark dataset
#'
#' Unified accessor for the tutorial datasets. dr3 ships no data; each dataset is
#' pulled from another package.
#'
#' @param dataset `"nhefs"` (via \pkg{causaldata}) or `"rhc"` (via \pkg{Hmisc}).
#' @param ... Passed to the specific loader (e.g. `outcome=` for `"rhc"`).
#' @return A list with `name`, `Y`, `A`, `X`, `family_y`, `outcome_type`, `truth`.
#' @seealso [dr_nhefs()], [dr_rhc()]
#' @examples
#' \dontrun{ d <- dr_data("nhefs") }
#' @export
dr_data <- function(dataset = c("nhefs", "rhc"), ...) {
  dataset <- match.arg(dataset)
  switch(dataset, nhefs = dr_nhefs(...), rhc = dr_rhc(...))
}

#' NHEFS: smoking cessation and weight change (real, observational)
#'
#' The canonical Hernan & Robins *What If* example: quitting smoking (`qsmk`)
#' and 1971-1982 weight change (`wt82_71`). No ground-truth ATE (observational);
#' the literature benchmark is about 3.4-3.5 kg. Loaded from \pkg{causaldata}.
#'
#' @return A list with `name`, `Y`, `A`, `X`, `family_y`, `outcome_type`, `truth`.
#' @examples
#' \dontrun{ str(dr_nhefs()) }
#' @export
dr_nhefs <- function() {
  .dr_need("causaldata")
  d <- causaldata::nhefs
  covs <- c("sex", "age", "race", "education", "smokeintensity",
            "smokeyrs", "exercise", "active", "wt71")
  keep <- complete.cases(d[, c("qsmk", "wt82_71", covs)])
  d <- d[keep, ]
  list(name = "NHEFS (smoking cessation -> weight change)",
       Y = d$wt82_71, A = as.integer(d$qsmk),
       X = .dr_impute_median(as.data.frame(lapply(d[, covs], as.numeric))),
       family_y = stats::gaussian(), outcome_type = "continuous (kg)",
       truth = NA_real_)
}

#' RHC: right-heart catheterization and 30-day mortality (real, observational)
#'
#' Connors et al. (1996): does right-heart catheterization within 24h of ICU
#' admission (`swang1`) affect mortality? Binary outcome, no ground-truth ATE.
#' [Hmisc::getHdata()] downloads Vanderbilt's R-format export of the dataset
#' (`rhc.sav`) from the public Vanderbilt biostatistics repository; neither dr3 nor
#' Hmisc bundles it. Requires \pkg{Hmisc} and an internet connection.
#'
#' The paper analyses Vanderbilt's CSV export of the same dataset (`rhc.csv`). The
#' two exports differ slightly: four factors have different reference levels, and
#' numeric covariates differ by at most about 5e-6. Estimates from this loader can
#' therefore differ slightly from the paper's.
#'
#' @param outcome `"dth30"` (30-day mortality, default) or `"death"`.
#' @return A list with the common dataset structure.
#' @examples
#' \dontrun{ d <- dr_rhc() }
#' @export
dr_rhc <- function(outcome = c("dth30", "death")) {
  .dr_need("Hmisc")
  outcome <- match.arg(outcome)
  ## getHdata() assigns `rhc` into the global environment; capture it and restore
  ## the global environment so the loader has no side effects.
  ge  <- globalenv()
  had <- exists("rhc", envir = ge, inherits = FALSE)
  old <- if (had) get("rhc", envir = ge, inherits = FALSE)
  on.exit(if (had) assign("rhc", old, envir = ge)
          else suppressWarnings(rm(list = "rhc", envir = ge)), add = TRUE)
  suppressMessages(Hmisc::getHdata("rhc"))
  d <- get("rhc", envir = ge, inherits = FALSE)

  A <- as.integer(d$swang1 == "RHC")
  Y <- as.integer(d[[outcome]] == "Yes")
  drop <- c("swang1", "death", "dth30", "sadmdte", "dschdte", "dthdte", "lstctdte",
            "t3d30", "ptid", "cat2", "adld3p", "urin1")  # treatment, outcomes, dates/ids, high-NA
  covs <- droplevels(d[, setdiff(names(d), drop), drop = FALSE])
  X <- as.data.frame(model.matrix(~ . - 1, data = covs))  # one-hot the factors
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  list(name = paste0("RHC (right heart cath -> ", outcome, ")"),
       Y = Y, A = A, X = X,
       family_y = stats::binomial(), outcome_type = "binary (mortality risk difference)",
       truth = NA_real_)
}
