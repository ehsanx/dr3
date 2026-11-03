## =============================================================================
## panel.R -- run the full estimator x nuisance-strategy grid in one call, and
## give the result a first-class `dr_panel` object with print/summary/plot.
## =============================================================================

#' Run the full doubly-robust estimator panel
#'
#' Fits the shared Super Learner nuisances under all three schemes (full-sample,
#' cross-fit, double cross-fit) and runs every estimator on them, returning the
#' eight-row panel that is the organizing idea of *Three Routes to One Answer*:
#' the singly-robust baselines (g-computation, IPTW) plus the six doubly-robust
#' estimators (AIPW, TMLE, DML, CV-TMLE, DC-AIPW, DC-TMLE). The estimators share
#' one Super Learner library, so differences between them come from the
#' combination rule and the fitting scheme (and, for split-based estimators, from
#' the random split).
#'
#' The full-sample fits draw Super Learner's internal cross-validation folds from
#' the current random-number state; call [set.seed()] before `dr_panel()` for
#' reproducible results. `seed` sets the folds of the cross-fit and
#' double-cross-fit schemes (see [dr_folds()]).
#'
#' @inheritParams dr_nuisance
#' @param y_range Optional outcome range for the TMLE scaling; inferred as
#'   `c(0, 1)` automatically when `Y` is binary.
#' @param truth Optional known ATE. When supplied, a `bias` column is filled in
#'   (useful on semi-synthetic benchmarks such as IHDP or Twins).
#' @return A [data.frame] of class `dr_panel` with one row per
#'   estimator x nuisance strategy (`estimator`, `nuisance`, `estimate`, `se`,
#'   `ci_lo`, `ci_hi`, `bias`), carrying `truth`, `n`, and `call` attributes.
#' @seealso [dr_nuisance()], [dr_aipw()], [dr_tmle()], [plot.dr_panel()]
#' @examples
#' \dontrun{
#' d <- dr_nhefs()
#' set.seed(1)
#' fit <- dr_panel(d$Y, d$A, d$X)
#' fit
#' plot(fit)
#' }
#' @export
dr_panel <- function(Y, A, X, family_y = stats::gaussian(), y_range = NULL,
                     V = 5L, Vd = 3L, SL.library = default_dr_library, seed = 1L,
                     truth = NA_real_) {
  binary_y <- isTRUE(all.equal(sort(unique(Y)), c(0, 1)))
  if (binary_y) y_range <- c(0, 1)
  nu_F <- .nuisance_full(Y, A, X, SL.library, family_y)
  nu_X <- .nuisance_crossfit(Y, A, X, V,  SL.library, family_y, seed = seed)
  nu_D <- .nuisance_double(Y, A, X, Vd, SL.library, family_y, seed = seed)

  mk <- function(res, label, strat) data.frame(
    estimator = label, nuisance = strat,
    estimate = res$estimate, se = res$se,
    ci_lo = res$ci[1], ci_hi = res$ci[2],
    bias = res$estimate - truth, stringsAsFactors = FALSE)

  out <- rbind(
    mk(dr_gcomp(Y, A, nu_F), "g-computation",   "full"),
    mk(dr_iptw (Y, A, nu_F), "IPTW",            "full"),
    mk(dr_aipw (Y, A, nu_F), "AIPW",            "full"),
    mk(dr_tmle (Y, A, nu_F, y_range), "TMLE",   "full"),
    mk(dr_aipw (Y, A, nu_X), "AIPW = DML",      "cross-fit"),
    mk(dr_tmle (Y, A, nu_X, y_range), "TMLE (CV-TMLE)",  "cross-fit"),
    mk(dr_aipw (Y, A, nu_D), "AIPW (DC-AIPW)",  "double-cross-fit"),
    mk(dr_tmle (Y, A, nu_D, y_range), "TMLE (DC-TMLE)",  "double-cross-fit")
  )
  attr(out, "truth") <- truth
  attr(out, "n")     <- length(Y)
  attr(out, "call")  <- match.call()
  class(out) <- c("dr_panel", "data.frame")
  out
}

#' @export
print.dr_panel <- function(x, digits = 3, ...) {
  truth <- attr(x, "truth"); n <- attr(x, "n")
  cat(sprintf("dr3 estimator panel  (n = %s%s)\n", n,
              if (!is.na(truth)) sprintf(", truth = %.3f", truth) else ""))
  df <- as.data.frame(x)
  df$estimate <- round(df$estimate, digits)
  ci <- ifelse(is.na(df$ci_lo), "",
               sprintf("(%.*f, %.*f)", digits, df$ci_lo, digits, df$ci_hi))
  show <- data.frame(estimator = df$estimator, nuisance = df$nuisance,
                     estimate = df$estimate, `95% CI` = ci, check.names = FALSE)
  if (!is.na(truth)) show$bias <- round(df$bias, digits)
  print(show, row.names = FALSE, ...)
  invisible(x)
}

#' @export
summary.dr_panel <- function(object, ...) {
  df <- as.data.frame(object)
  dr <- df[grepl("AIPW|TMLE|DML", df$estimator), "estimate"]
  cat(sprintf("Doubly-robust estimates span %.3f (%.3f to %.3f) across %d routes.\n",
              diff(range(dr)), min(dr), max(dr), length(dr)))
  print(object, ...)
  invisible(object)
}

#' Forest plot of a dr3 estimator panel
#'
#' A base-graphics forest plot: point estimates with 95% Wald intervals for the
#' doubly-robust estimators (the singly-robust plug-ins have no analytic interval),
#' with an optional reference line at the known `truth`.
#'
#' @param x A `dr_panel` object from [dr_panel()].
#' @param ref Optional reference value drawn as a vertical line; defaults to the
#'   panel's `truth` attribute when available.
#' @param xlab Axis label.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.dr_panel <- function(x, ref = attr(x, "truth"), xlab = "ATE estimate", ...) {
  df <- as.data.frame(x)
  k  <- nrow(df); ypos <- rev(seq_len(k))
  xr <- range(c(df$estimate, df$ci_lo, df$ci_hi, ref), na.rm = TRUE)
  pad <- diff(xr) * 0.08; xr <- xr + c(-pad, pad)
  cols <- ifelse(grepl("TMLE", df$estimator), "#D55E00",
          ifelse(grepl("AIPW|DML", df$estimator), "#0072B2", "grey45"))
  op <- graphics::par(mar = c(4.2, 11, 1, 1)); on.exit(graphics::par(op))
  plot(NA, xlim = xr, ylim = c(0.5, k + 0.5), yaxt = "n",
       xlab = xlab, ylab = "", ...)
  if (length(ref) && !is.na(ref)) graphics::abline(v = ref, lty = 2, col = "grey60")
  graphics::segments(df$ci_lo, ypos, df$ci_hi, ypos, col = cols, lwd = 2)
  graphics::points(df$estimate, ypos, pch = 19, col = cols)
  graphics::axis(2, at = ypos, labels = df$estimator, las = 1, tick = FALSE)
  invisible(x)
}
