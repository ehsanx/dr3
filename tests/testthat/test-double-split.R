# The double cross-fit's splitting rule, checked without fitting a model: a recording stub
# replaces the Super Learner fit and logs which rows train each g and Q fit, and which rows
# each fit scores. Column x1 carries the row id.

record_double_split <- function(n = 30L) {
  X <- data.frame(x1 = seq_len(n), x2 = stats::rnorm(n))
  A <- rep(0:1, length.out = n); Y <- stats::rnorm(n)
  log <- new.env(); log$fits <- list(); log$preds <- list()
  stub <- function(Y, X, family, SL.library, ...) {
    id <- length(log$fits) + 1L
    kind <- if ("A" %in% names(X)) "Q" else "g"
    log$fits[[id]] <- list(kind = kind, train = X$x1)
    function(newX) {
      log$preds[[length(log$preds) + 1L]] <- list(fit = id, kind = kind, rows = newX$x1)
      rep(if (kind == "g") 0.5 else 0, nrow(newX))
    }
  }
  local_mocked_bindings(.dr_fit_sl = stub, .env = parent.frame())
  nu <- dr_nuisance(Y, A, X, strategy = "double", seed = 1)
  list(n = n, nu = nu, fits = log$fits, preds = log$preds)
}

test_that("each fold is scored by g and Q trained on disjoint other folds", {
  r <- record_double_split()
  g_pred <- Filter(function(p) p$kind == "g", r$preds)
  q_pred <- Filter(function(p) p$kind == "Q", r$preds)
  expect_length(g_pred, 3L)
  expect_length(q_pred, 6L)                          # Q scores each fold under A = 1 and A = 0
  for (p in g_pred) {
    tr_g <- r$fits[[p$fit]]$train
    qs <- Filter(function(q) setequal(q$rows, p$rows), q_pred)
    expect_length(qs, 2L)
    expect_identical(qs[[1]]$fit, qs[[2]]$fit)
    tr_q <- r$fits[[qs[[1]]$fit]]$train
    expect_length(intersect(p$rows, tr_g), 0L)       # scored rows excluded from both fits
    expect_length(intersect(p$rows, tr_q), 0L)
    expect_length(intersect(tr_g, tr_q), 0L)         # g and Q never share training rows
  }
})

test_that("every observation is scored exactly once", {
  r <- record_double_split()
  scored <- unlist(lapply(Filter(function(p) p$kind == "g", r$preds), `[[`, "rows"))
  expect_setequal(scored, seq_len(r$n))
  expect_equal(anyDuplicated(scored), 0L)
  expect_equal(r$nu$pi, rep(0.5, r$n))
})

test_that("the double cross-fit requires three folds", {
  d <- sim_data(n = 60)
  expect_error(dr_nuisance(d$Y, d$A, d$X, strategy = "double", Vd = 4L, SL.library = test_lib))
})
