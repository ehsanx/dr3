## =============================================================================
## estimators.R -- ATE estimators that are AGNOSTIC to how the nuisances were fit.
## Each takes (Y, A, nu), where `nu` comes from dr_nuisance(). Pair an estimator
## with a nuisance strategy to get a named method:
##   dr_aipw + full      -> AIPW (one-step)      dr_tmle + full      -> TMLE
##   dr_aipw + crossfit  -> DML                  dr_tmle + crossfit  -> CV-TMLE
##   dr_aipw + double    -> DC-AIPW              dr_tmle + double    -> DC-TMLE
## =============================================================================

#' g-computation (outcome-regression plug-in)
#'
#' Singly robust: consistent only if the outcome model `Q` is correct.
#' @param Y Numeric outcome vector.
#' @param A Binary treatment vector (0/1).
#' @param nu A nuisance object from [dr_nuisance()].
#' @return A list with `estimate`, `se`, `ci`, `method`.
#' @examples
#' \dontrun{ dr_gcomp(d$Y, d$A, dr_nuisance(d$Y, d$A, d$X, "full")) }
#' @export
dr_gcomp <- function(Y, A, nu) {
  psi <- mean(nu$Q1 - nu$Q0)
  list(estimate = psi, se = NA_real_, ci = c(NA, NA), method = "g-computation")
}

#' Inverse-probability-of-treatment weighting (Hajek / stabilized)
#'
#' Singly robust: consistent only if the propensity model `g` is correct.
#' @inheritParams dr_gcomp
#' @return A list with `estimate`, `se`, `ci`, `method`.
#' @export
dr_iptw <- function(Y, A, nu) {
  pi <- nu$pi
  m1 <- sum(A * Y / pi) / sum(A / pi)
  m0 <- sum((1 - A) * Y / (1 - pi)) / sum((1 - A) / (1 - pi))
  list(estimate = m1 - m0, se = NA_real_, ci = c(NA, NA), method = "IPTW (Hajek)")
}

#' AIPW / one-step estimator (= DML with cross-fit nuisances)
#'
#' Doubly robust: consistent if EITHER `g` or `Q` is correct. The estimate is the
#' plug-in plus the mean of the efficient influence function; the standard error
#' comes from the influence curve (no bootstrap). With `strategy = "crossfit"`
#' nuisances this is exactly double machine learning (DML).
#' @inheritParams dr_gcomp
#' @return A list with `estimate`, `se`, `ci`, `ic` (influence curve), `method`.
#' @references Robins, Rotnitzky & Zhao (1994); Chernozhukov et al. (2018).
#' @export
dr_aipw <- function(Y, A, nu) {
  pi <- nu$pi; Q1 <- nu$Q1; Q0 <- nu$Q0
  psi_i <- (Q1 - Q0) + (A / pi) * (Y - Q1) - ((1 - A) / (1 - pi)) * (Y - Q0)
  est <- mean(psi_i)
  ic  <- psi_i - est
  se  <- sqrt(var(ic) / length(Y))
  list(estimate = est, se = se, ci = est + c(-1, 1) * 1.96 * se,
       ic = ic, method = "AIPW / one-step")
}

#' Targeted maximum likelihood estimator (TMLE)
#'
#' Doubly robust. Fluctuates the outcome fit `Q` through a logistic submodel with
#' the clever covariate `H(A, X)`, then plugs in. Continuous outcomes are scaled
#' to `[0, 1]` (Gruber & van der Laan) so the logistic update is valid; binary
#' outcomes use `y_range = c(0, 1)` (the scaling is then the identity). With
#' `strategy = "crossfit"` nuisances this is CV-TMLE.
#' @inheritParams dr_gcomp
#' @param y_range Optional outcome range `c(min, max)` used for scaling. Defaults
#'   to `range(Y)`; pass `c(0, 1)` for binary outcomes. A zero range is an error,
#'   since the outcome cannot then be rescaled.
#' @return A list with `estimate`, `se`, `ci`, `ic`, `eps` (fluctuation), `method`.
#' @references van der Laan & Rubin (2006); Gruber & van der Laan (2012).
#' @export
dr_tmle <- function(Y, A, nu, y_range = NULL) {
  pi <- nu$pi
  if (is.null(y_range)) y_range <- range(Y)
  a <- y_range[1]; b <- y_range[2]; span <- b - a
  if (!(span > 0)) stop("dr_tmle: the outcome range is zero, so it cannot be rescaled to [0, 1]")
  Ys  <- (Y      - a) / span
  Q1s <- .dr_bound((nu$Q1 - a) / span, 1e-5, 1 - 1e-5)
  Q0s <- .dr_bound((nu$Q0 - a) / span, 1e-5, 1 - 1e-5)
  QAs <- ifelse(A == 1, Q1s, Q0s)
  H  <- A / pi - (1 - A) / (1 - pi)
  H1 <- 1 / pi; H0 <- -1 / (1 - pi)
  eps <- coef(suppressWarnings(
    glm(Ys ~ -1 + H, offset = qlogis(QAs), family = binomial())))
  Q1s_star <- plogis(qlogis(Q1s) + eps * H1)
  Q0s_star <- plogis(qlogis(Q0s) + eps * H0)
  QAs_star <- plogis(qlogis(QAs) + eps * H)
  Q1_star <- Q1s_star * span + a; Q0_star <- Q0s_star * span + a
  est <- mean(Q1_star - Q0_star)
  ic <- span * H * (Ys - QAs_star) + (Q1_star - Q0_star) - est
  se <- sqrt(var(ic) / length(Y))
  list(estimate = est, se = se, ci = est + c(-1, 1) * 1.96 * se,
       ic = ic, eps = eps, method = "TMLE")
}
