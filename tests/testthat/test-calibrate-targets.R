test_that("degree-preserving swaps preserve both margins and source sign counts", {
  net <- data.frame(
    source = rep(c("A", "B", "C"), each = 3),
    target = c("g1", "g2", "g3", "g3", "g4", "g5", "g1", "g5", "g6"),
    sign = c(1, 1, -1, 1, -1, 1, -1, 1, 1)
  )
  set.seed(11)
  out <- perturbRNet:::.prn_degree_preserving_swap(net, 500)$network
  expect_equal(sort(table(out$source)), sort(table(net$source)))
  expect_equal(sort(table(out$target)), sort(table(net$target)))
  before <- aggregate(sign ~ source, net, function(x) paste(sort(x), collapse = ","))
  after <- aggregate(sign ~ source, out, function(x) paste(sort(x), collapse = ","))
  expect_equal(before[order(before$source), ], after[order(after$source), ])
  expect_false(anyDuplicated(out[c("source", "target")]) > 0)
})

test_that("target calibration returns calibrated target-comparison rows", {
  fingerprint <- matrix(
    c(0.2, -0.1, 0.4, -0.3, 0.1, -0.2,
      -0.1, 0.2, -0.3, 0.4, -0.2, 0.1),
    nrow = 6,
    dimnames = list(paste0("g", 1:6), c("c1", "c2"))
  )
  prior <- data.frame(
    source = rep(c("A", "B", "C"), each = 3),
    target = c("g1", "g2", "g3", "g3", "g4", "g5", "g1", "g5", "g6"),
    sign = c(1, 1, -1, 1, -1, 1, -1, 1, 1)
  )
  fit <- prn_build(fingerprint, prior, min_targets = 1, lambda = 0.1)
  result <- prn_calibrate_targets(
    fit, fingerprint, prior, n_null = 3, swaps_per_edge = 10,
    min_targets = 1, lambda = 0.1, seed = 7
  )
  expect_true(all(c(
    "incoming_degree", "null_mean", "null_sd", "degree_conditioned_z",
    "null_percentile", "empirical_p_upper"
  ) %in% names(result$calibration)))
  expect_equal(nrow(result$calibration), nrow(prn_counterdirection(fit)))
  expect_true(all(result$calibration$empirical_p_upper >= 1 / 4))
  expect_true(all(result$calibration$empirical_p_upper <= 1))
})
