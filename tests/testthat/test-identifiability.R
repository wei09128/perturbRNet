test_that("identifiability detects separable and confounded regulators", {
  fingerprint <- matrix(
    c(0, 0, 1, 1, -1, -1),
    ncol = 1,
    dimnames = list(paste0("G", 1:6), "perturbation")
  )
  prior <- data.frame(
    source = c("A", "B", "A", "A", "B", "B"),
    target = c("G1", "G1", "G3", "G4", "G5", "G6"),
    sign = c(1, -1, 1, 1, 1, 1)
  )
  fit <- prn_build(fingerprint, prior, min_targets = 1, lambda = 0.1)
  diagnostic <- prn_identifiability(fit)

  expect_equal(diagnostic$global$rank, 2)
  expect_true(all(diagnostic$regulator$unique_information > 0))
  expect_true(all(diagnostic$pair$pair_separation > 0))
  expect_true(all(c("global", "regulator", "pair", "target") %in% names(diagnostic)))
})

test_that("duplicate regulator footprints are flagged as non-identifiable", {
  W <- cbind(A = c(1, 1, -1), B = c(1, 1, -1))
  object <- list(weight_matrix = W)
  class(object) <- "prn_fit"

  local_mocked_bindings(
    prn_edge_table = function(object) data.frame(
      comparison = "x", target = "G1", regulator = c("A", "B"),
      inferred_contribution = c(1, -1)
    ),
    .package = "perturbRNet"
  )
  diagnostic <- prn_identifiability(object)

  expect_lt(diagnostic$global$rank_fraction, 1)
  expect_equal(diagnostic$pair$pair_separation, 0, tolerance = 1e-12)
  expect_equal(diagnostic$target$identifiability, "low")
})
