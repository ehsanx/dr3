test_that("dr_panel returns the full eight-row grid with the expected shape", {
  d   <- sim_data()
  fit <- dr_panel(d$Y, d$A, d$X, SL.library = test_lib, truth = d$truth)
  expect_s3_class(fit, "dr_panel")
  expect_identical(nrow(fit), 8L)
  expect_true(all(c("estimator", "nuisance", "estimate", "se",
                    "ci_lo", "ci_hi", "bias") %in% names(fit)))
  expect_true(all(is.finite(fit$estimate)))
  expect_identical(attr(fit, "n"), length(d$Y))
  # the six doubly-robust estimates should cluster near the truth
  dr <- fit$estimate[grepl("AIPW|TMLE|DML", fit$estimator)]
  expect_true(all(abs(dr - d$truth) < 0.6))
  # singly-robust plug-ins carry no analytic interval
  expect_true(all(is.na(fit$ci_lo[fit$estimator %in% c("g-computation", "IPTW")])))
})

test_that("S3 methods run without error", {
  d   <- sim_data()
  fit <- dr_panel(d$Y, d$A, d$X, SL.library = test_lib, truth = d$truth)
  expect_output(print(fit))
  expect_output(summary(fit))
  tmpf <- tempfile(fileext = ".png"); grDevices::png(tmpf)
  expect_invisible(plot(fit))
  grDevices::dev.off(); unlink(tmpf)
})

test_that("binary outcomes run and stay on the risk-difference scale", {
  d   <- sim_binary()
  fit <- dr_panel(d$Y, d$A, d$X, family_y = binomial(), SL.library = test_lib,
                  truth = d$truth)
  expect_identical(nrow(fit), 8L)
  expect_true(all(fit$estimate > -1 & fit$estimate < 1))
})
