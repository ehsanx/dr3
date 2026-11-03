## =============================================================================
## nuisances.R -- shared Super Learner nuisance estimation
## g(X) = P(A = 1 | X) and Q_a(X) = E[Y | A = a, X], fit with ONE library so that
## differences between estimators come from the combination rule and the fitting scheme
## (and, for split-based estimators, from the split), not from separate implementations.
## =============================================================================

#' Default Super Learner library
#'
#' The Super Learner library shared by the propensity (`g`) and outcome (`Q`)
#' nuisance regressions. Sharing one library is what makes the estimators in
#' [dr_panel()] comparable: differences between them then come from the
#' combination rule and the fitting scheme (and, for split-based estimators, from
#' the random split), not from separate implementations.
#'
#' @format A character vector of \pkg{SuperLearner} wrapper names.
#' @export
default_dr_library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")

.dr_bound <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' Construct cross-fitting folds
#'
#' @param n Integer sample size.
#' @param V Number of folds.
#' @param seed RNG seed for reproducible fold assignment.
#' @details Calls [set.seed()] with `seed`, so it also resets the global
#'   random-number state.
#' @return An integer vector of length `n` with fold ids in `1:V`.
#' @examples
#' table(dr_folds(20, V = 5))
#' @export
dr_folds <- function(n, V = 5L, seed = 1L) {
  set.seed(seed)
  sample(rep(seq_len(V), length.out = n))
}

## Fit one Super Learner and return a predictor closure (keeps calling code readable).
## `env = asNamespace("SuperLearner")` makes the learner and screening ("All")
## wrappers resolvable when dr3 is only imported, not attached (SuperLearner looks
## them up in `env`, which defaults to the caller's frame).
.dr_fit_sl <- function(Y, X, family, SL.library) {
  fit <- SuperLearner::SuperLearner(Y = Y, X = X, family = family,
                                    SL.library = SL.library,
                                    cvControl = list(V = 5L),
                                    env = asNamespace("SuperLearner"))
  function(newX) as.numeric(predict(fit, newX, onlySL = TRUE)$pred)
}

## Every strategy returns the SAME object so the downstream estimators are identical:
##   pi = g(X),  Q1 = E[Y|A=1,X],  Q0 = E[Y|A=0,X],  QA = Q at the observed A.

## (F) full-sample: fit on all data, predict in-sample  -> classic AIPW / TMLE.
.nuisance_full <- function(Y, A, X, SL.library = default_dr_library,
                           family_y = stats::gaussian(), gbound = c(.025, .975)) {
  Xdf <- as.data.frame(X)
  g   <- .dr_fit_sl(A, Xdf, stats::binomial(), SL.library)
  q   <- .dr_fit_sl(Y, cbind(A = A, Xdf), family_y, SL.library)
  pi  <- .dr_bound(g(Xdf), gbound[1], gbound[2])
  Q1  <- q(cbind(A = 1, Xdf)); Q0 <- q(cbind(A = 0, Xdf))
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0))
}

## (X) cross-fit (DML / CV-TMLE): V folds; g and Q trained on the same V-1 folds,
## predicted on the held-out fold => out-of-fold nuisances for every observation.
.nuisance_crossfit <- function(Y, A, X, V = 5L, SL.library = default_dr_library,
                               family_y = stats::gaussian(), gbound = c(.025, .975),
                               seed = 1L) {
  Xdf <- as.data.frame(X); n <- length(Y); fold <- dr_folds(n, V, seed)
  pi <- Q1 <- Q0 <- rep(NA_real_, n)
  for (v in seq_len(V)) {
    tr <- fold != v; te <- fold == v
    g <- .dr_fit_sl(A[tr], Xdf[tr, , drop = FALSE], stats::binomial(), SL.library)
    q <- .dr_fit_sl(Y[tr], cbind(A = A[tr], Xdf[tr, , drop = FALSE]), family_y, SL.library)
    pi[te] <- g(Xdf[te, , drop = FALSE])
    Q1[te] <- q(cbind(A = 1, Xdf[te, , drop = FALSE]))
    Q0[te] <- q(cbind(A = 0, Xdf[te, , drop = FALSE]))
  }
  pi <- .dr_bound(pi, gbound[1], gbound[2])
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0), fold = fold)
}

## (D) double cross-fit (DC-AIPW / DC-TMLE; Newey & Robins 2018; Zivich & Breskin 2021):
## split into V = 3 folds and use one fixed rotation of roles: observations in fold v are
## scored with g fit on fold v+1 alone and Q fit on fold v+2 alone (indices mod 3). Every
## scored observation is excluded from both training sets, and its g and Q never share
## training data.
.nuisance_double <- function(Y, A, X, V = 3L, SL.library = default_dr_library,
                             family_y = stats::gaussian(), gbound = c(.025, .975),
                             seed = 1L) {
  stopifnot(V == 3L)
  Xdf <- as.data.frame(X); n <- length(Y); fold <- dr_folds(n, V, seed)
  pi <- Q1 <- Q0 <- rep(NA_real_, n)
  for (v in seq_len(V)) {
    te   <- fold == v
    tr_g <- fold == (v %% V) + 1L            # the next fold trains g
    tr_q <- fold == ((v + 1L) %% V) + 1L     # the fold after that trains Q
    g <- .dr_fit_sl(A[tr_g], Xdf[tr_g, , drop = FALSE], stats::binomial(), SL.library)
    q <- .dr_fit_sl(Y[tr_q], cbind(A = A[tr_q], Xdf[tr_q, , drop = FALSE]),
                    family_y, SL.library)
    pi[te] <- g(Xdf[te, , drop = FALSE])
    Q1[te] <- q(cbind(A = 1, Xdf[te, , drop = FALSE]))
    Q0[te] <- q(cbind(A = 0, Xdf[te, , drop = FALSE]))
  }
  pi <- .dr_bound(pi, gbound[1], gbound[2])
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0), fold = fold)
}

#' Estimate the shared nuisances (propensity `g` and outcome `Q`)
#'
#' Fits the propensity score `g(X) = P(A = 1 | X)` and the outcome regression
#' `Q_a(X) = E[Y | A = a, X]` with a shared Super Learner library, under one of
#' three fitting schemes. Every scheme returns the same object shape, so the
#' downstream estimators ([dr_aipw()], [dr_tmle()], ...) are identical regardless
#' of how the nuisances were fit.
#'
#' @param Y Numeric outcome vector.
#' @param A Binary treatment vector (0/1).
#' @param X Covariates (data frame or matrix).
#' @param strategy Nuisance-fitting scheme: `"full"` (in-sample, classic
#'   AIPW/TMLE), `"crossfit"` (out-of-fold, DML / CV-TMLE), or `"double"`
#'   (double cross-fit, DC-AIPW / DC-TMLE). The double cross-fit scores fold `v`
#'   with `g` fit on fold `v + 1` alone and `Q` fit on fold `v + 2` alone (mod 3),
#'   so each scored observation is excluded from both fits and `g` and `Q` never
#'   share training data (Zivich & Breskin 2021).
#' @param V Number of folds for `"crossfit"`.
#' @param Vd Number of folds for `"double"`; must be 3 (the rotation uses one
#'   fold to score, one to fit `g`, and one to fit `Q`).
#' @param SL.library Super Learner library shared by `g` and `Q`
#'   (default [default_dr_library]).
#' @param family_y Outcome family: [stats::gaussian()] (continuous) or
#'   [stats::binomial()] (binary).
#' @param gbound Truncation bounds for the estimated propensity.
#' @param seed RNG seed for fold construction.
#' @details The `"full"` scheme sets no seed: Super Learner's internal
#'   cross-validation then draws from the current random-number state, so call
#'   [set.seed()] first for a reproducible full-sample fit. The cross-fit schemes
#'   build their folds with [dr_folds()], which calls `set.seed(seed)`.
#' @return A list with `pi` (bounded propensity), `Q1`, `Q0`, `QA`, and, for the
#'   cross-fit schemes, the `fold` assignment.
#' @seealso [dr_aipw()], [dr_tmle()], [dr_panel()]
#' @examples
#' \dontrun{
#' d  <- dr_nhefs()
#' nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "crossfit")
#' dr_aipw(d$Y, d$A, nu)   # = DML
#' }
#' @export
dr_nuisance <- function(Y, A, X,
                        strategy = c("full", "crossfit", "double"),
                        V = 5L, Vd = 3L,
                        SL.library = default_dr_library,
                        family_y = stats::gaussian(),
                        gbound = c(0.025, 0.975), seed = 1L) {
  strategy <- match.arg(strategy)
  switch(strategy,
    full     = .nuisance_full(Y, A, X, SL.library, family_y, gbound),
    crossfit = .nuisance_crossfit(Y, A, X, V,  SL.library, family_y, gbound, seed),
    double   = .nuisance_double(Y, A, X, Vd, SL.library, family_y, gbound, seed))
}
