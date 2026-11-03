#' @keywords internal
"_PACKAGE"

## dr3: doubly-robust ATE estimators, three routes to one answer.
##
## AIPW, TMLE, and double machine learning all target the SAME efficient
## influence function for the average treatment effect
##   D(O) = [A/g(X) - (1-A)/(1-g(X))] (Y - Q_A(X)) + (Q_1(X) - Q_0(X)) - psi
## and differ only in how they solve it (average the score vs. targeting) and how
## the nuisances g, Q are fit (full-sample, cross-fit, double cross-fit). dr3 is a
## transparent reference implementation that runs the whole family on ONE shared
## Super Learner library.

#' @importFrom stats binomial coef complete.cases gaussian glm median
#'   model.matrix plogis predict qlogis rbinom rnorm var
#' @importFrom SuperLearner SuperLearner
#' @importFrom graphics abline axis par points segments
NULL
