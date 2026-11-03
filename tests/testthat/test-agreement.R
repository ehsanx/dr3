# Cross-package agreement: dr3's hand-coded estimators should match an established
# package given the same data and Super Learner library. Exact agreement is not
# expected (fold schemes and internal truncation differ), so the tolerances are
# loose; the detailed reconciliation is in the software-bridge vignette.

test_that("dr_tmle (full) agrees with the tmle package on a shared library", {
  skip_on_cran()
  skip_if_not_installed("tmle")
  d  <- sim_data(n = 600)
  nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "full", SL.library = test_lib)
  ours <- dr_tmle(d$Y, d$A, nu)$estimate

  set.seed(1)
  pkg <- tmle::tmle(Y = d$Y, A = d$A, W = as.data.frame(d$X),
                    Q.SL.library = test_lib, g.SL.library = test_lib,
                    family = "gaussian")
  theirs <- pkg$estimates$ATE$psi

  expect_equal(sign(ours - 1), sign(theirs - 1))   # same direction vs. null
  expect_lt(abs(ours - theirs), 0.4)               # loose numeric agreement
})

test_that("dr_aipw (cross-fit = DML) agrees with the AIPW package", {
  skip_on_cran()
  skip_if_not_installed("AIPW")
  skip_if_not_installed("SuperLearner")
  withr::local_package("SuperLearner")   # AIPW resolves its SL wrappers on the search path
  d  <- sim_data(n = 800)
  nu <- dr_nuisance(d$Y, d$A, d$X, strategy = "crossfit", SL.library = test_lib)
  ours <- dr_aipw(d$Y, d$A, nu)$estimate

  set.seed(1)
  pkg <- AIPW::AIPW$new(Y = d$Y, A = d$A, W = as.data.frame(d$X),
                        Q.SL.library = test_lib, g.SL.library = test_lib,
                        k_split = 5, verbose = FALSE)$fit()$summary()
  theirs <- as.numeric(pkg$estimates$RD["Estimate"])

  expect_true(is.finite(theirs))
  expect_lt(abs(ours - theirs), 0.5)
})
