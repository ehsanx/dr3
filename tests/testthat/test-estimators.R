test_that("AIPW estimate equals the mean of its efficient influence function", {
  d  <- sim_data()
  nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "full", SL.library = test_lib)
  res <- dr_aipw(d$Y, d$A, nu)
  eif <- with(nu, (Q1 - Q0) + (d$A / pi) * (d$Y - Q1) -
                  ((1 - d$A) / (1 - pi)) * (d$Y - Q0))
  expect_equal(res$estimate, mean(eif))
  # the influence curve is centred and yields the reported SE
  expect_equal(mean(res$ic), 0, tolerance = 1e-8)
  expect_equal(res$se, sqrt(var(res$ic) / length(d$Y)))
})

test_that("TMLE solves its efficient score (clever-covariate equation ~ 0)", {
  d  <- sim_data()
  nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "full", SL.library = test_lib)
  res <- dr_tmle(d$Y, d$A, nu)
  # after fluctuation, mean( H(A,X) * (Y - Q*_A) ) must be ~ 0 by construction
  pi <- nu$pi; H <- d$A / pi - (1 - d$A) / (1 - pi)
  span <- diff(range(d$Y)); a <- min(d$Y)
  QAs  <- ifelse(d$A == 1, .safe01((nu$Q1 - a) / span), .safe01((nu$Q0 - a) / span))
  QAs_star <- plogis(qlogis(QAs) + res$eps * H)
  score <- mean(H * ((d$Y - a) / span - QAs_star))
  expect_lt(abs(score), 1e-6)
})

test_that("TMLE stops on a zero outcome range", {
  nu <- list(pi = rep(0.5, 4), Q1 = rep(1, 4), Q0 = rep(1, 4))
  expect_error(dr_tmle(rep(1, 4), c(0, 1, 0, 1), nu), "outcome range is zero")
})

test_that("both doubly-robust estimators recover a known ATE", {
  d  <- sim_data()
  nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "full", SL.library = test_lib)
  expect_equal(dr_aipw(d$Y, d$A, nu)$estimate, d$truth, tolerance = 0.4)
  expect_equal(dr_tmle(d$Y, d$A, nu)$estimate, d$truth, tolerance = 0.4)
})
