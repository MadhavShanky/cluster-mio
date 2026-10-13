# Regression / sanity tests. The planted-truth recovery and the determinism check
# guard the property that matters for a release: a given (data, lambda, seed) maps
# to a fixed selected set and objective. If packaging or an "efficiency" refactor
# ever changes the numerics, these fail loudly.

make_data <- function(K = 30, nk = 20, deviating = c(5, 22), eff = c(5, -5), seed = 1) {
  set.seed(seed)
  n <- K * nk
  cluster <- factor(rep(seq_len(K), each = nk))
  X <- matrix(rnorm(n * 2), n, 2)
  g <- numeric(K); g[deviating] <- eff
  y <- as.numeric(X %*% c(1, -1)) + g[as.integer(cluster)] + rnorm(n)
  list(X = X, y = y, cluster = cluster, deviating = as.character(deviating))
}

test_that("symmetric fit recovers the planted deviating clusters", {
  d <- make_data()
  fit <- scs_fit(d$X, d$y, d$cluster, lambda = 2)
  expect_s3_class(fit, "scs")
  expect_setequal(fit$selected, d$deviating)
  expect_true(is.finite(fit$obj))
})

test_that("fits are deterministic given the seed", {
  d <- make_data()
  f1 <- scs_fit(d$X, d$y, d$cluster, lambda = 2, seed = 7)
  f2 <- scs_fit(d$X, d$y, d$cluster, lambda = 2, seed = 7)
  expect_identical(f1$selected, f2$selected)
  expect_equal(f1$obj, f2$obj)
  expect_equal(unname(f1$gamma), unname(f2$gamma))
})

test_that("asymmetric tail budgets respect the directional caps", {
  d <- make_data()
  fit <- scs_fit(d$X, d$y, d$cluster, lambda = c(worst = 1, best = 1))
  expect_lte(length(fit$selected_pos), 1L)
  expect_lte(length(fit$selected_neg), 1L)
})

test_that("objective is nonincreasing along the budget path", {
  d <- make_data()
  p <- scs_path(d$X, d$y, d$cluster, lambda_max = 4)
  expect_length(p$obj, 5L)
  # larger L0 budget cannot worsen the (minimized) selection objective beyond tolerance
  expect_true(all(diff(p$obj) <= 1e-8))
})

test_that("complementary-pairs stability validates its design arguments", {
  d <- make_data()
  fit <- scs_fit(d$X, d$y, d$cluster, lambda = 2)
  expect_error(stability(fit, B = 3L, pairing = "complementary"), "even B")
  expect_error(stability(fit, B = 10L, frac = 0.4, pairing = "complementary"), "frac=0.5")
  st <- stability(fit, B = 20L, pairing = "complementary")
  expect_s3_class(st, "scs_stability")
  expect_true(is.finite(st$expected_false_flags_bound))
})
