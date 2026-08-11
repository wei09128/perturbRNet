test_that("counterdirection separates cancellation from inactivity", {
  fingerprint <- matrix(
    c(0, 1, -1),
    ncol = 1,
    dimnames = list(c("g0", "g1", "g2"), "perturbation")
  )
  prior <- data.frame(
    source = c("A", "A", "B", "B"),
    target = c("g1", "g0", "g2", "g0"),
    sign = c(1, 1, 1, -1)
  )

  fit <- prn_build(fingerprint, prior, min_targets = 1, lambda = 0.1)
  result <- prn_counterdirection(fit)

  expect_true(all(result$counterdirection >= -1e-10))
  expect_true(all(result$balance_fraction >= -1e-10))
  expect_true(all(result$balance_fraction <= 1 + 1e-10))
})

