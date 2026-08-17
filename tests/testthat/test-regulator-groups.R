make_group_test_fit <- function(W, activity = seq_len(ncol(W))) {
  regulators <- colnames(W)
  genes <- rownames(W)
  object <- list(
    weight_matrix = W,
    comparisons = "perturbation",
    fits = list(perturbation = list(
      activity = stats::setNames(as.numeric(activity), regulators)
    ))
  )
  class(object) <- "prn"
  object
}

test_that("identical footprints are grouped", {
  W <- cbind(A = c(1, 1, 0), B = c(1, 1, 0), C = c(0, 0, 1))
  rownames(W) <- paste0("G", 1:3)
  result <- prn_group_regulators(make_group_test_fit(W, c(2, 3, 4)))

  ab <- result$membership[result$membership$regulator %in% c("A", "B"), ]
  expect_equal(length(unique(ab$group)), 1)
  expect_true(all(ab$orientation_to_anchor == 1))
  expect_false(all(ab$individual_activity_identifiable))
  expect_equal(
    result$activities$identifiable_group_activity[
      result$activities$group == ab$group[1]
    ],
    5
  )
})

test_that("opposite footprints are grouped with reversed orientation", {
  W <- cbind(A = c(1, 1, 0), B = c(-1, -1, 0), C = c(0, 0, 1))
  rownames(W) <- paste0("G", 1:3)
  result <- prn_group_regulators(make_group_test_fit(W, c(5, 2, 0)))

  ab <- result$membership[result$membership$regulator %in% c("A", "B"), ]
  expect_equal(length(unique(ab$group)), 1)
  expect_equal(ab$orientation_to_anchor[match(c("A", "B"), ab$regulator)], c(1, -1))
  expect_equal(
    result$activities$identifiable_group_activity[
      result$activities$group == ab$group[1]
    ],
    3
  )
})

test_that("distinct footprints remain separate", {
  W <- diag(3)
  colnames(W) <- c("A", "B", "C")
  rownames(W) <- paste0("G", 1:3)
  result <- prn_group_regulators(make_group_test_fit(W))

  expect_equal(nrow(result$groups), 3)
  expect_true(all(result$membership$individual_activity_identifiable))
  expect_false(any(result$pairs$same_group))
})
