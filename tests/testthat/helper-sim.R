# A small, well-behaved confounded dataset with a known ATE, used across tests.
# Continuous outcome, true ATE = 2. Cheap enough to fit SuperLearner repeatedly.
sim_data <- function(n = 400, seed = 1L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rbinom(n, 1, 0.5)
  ps <- plogis(0.6 * x1 - 0.5 * x2 + 0.3 * x3)
  A  <- rbinom(n, 1, ps)
  Y  <- 1 + 2 * A + x1 + 0.5 * x2 - x3 + rnorm(n)   # true ATE = 2
  list(Y = Y, A = A, X = data.frame(x1 = x1, x2 = x2, x3 = x3), truth = 2)
}

# A binary-outcome version (risk difference), true-ish ATE via potential outcomes.
sim_binary <- function(n = 500, seed = 2L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ps <- plogis(0.5 * x1 - 0.4 * x2)
  A  <- rbinom(n, 1, ps)
  p1 <- plogis(-0.2 + 0.8 + x1 * 0.5)
  p0 <- plogis(-0.2 + x1 * 0.5)
  y1 <- rbinom(n, 1, p1); y0 <- rbinom(n, 1, p0)
  Y  <- A * y1 + (1 - A) * y0
  list(Y = Y, A = A, X = data.frame(x1 = x1, x2 = x2), truth = mean(p1 - p0))
}

# A fast, dependency-light library for tests (no glmnet/ranger/earth needed).
test_lib <- c("SL.mean", "SL.glm")

# mirrors the internal [0,1] guard used by dr_tmle
.safe01 <- function(x) pmin(pmax(x, 1e-5), 1 - 1e-5)
