test_that("mostly positive regulons do not hide a counterbalanced motif", {
  set.seed(101)
  genes <- paste0("g", seq_len(200))
  cancel <- genes[1:20]

  prior <- rbind(
    data.frame(source = "A", target = cancel, sign = 1),
    data.frame(source = "B", target = cancel, sign = -1),
    data.frame(source = "A", target = genes[21:100], sign = 1),
    data.frame(source = "B", target = genes[101:180], sign = 1)
  )

  response <- setNames(rep(0, length(genes)), genes)
  response[21:180] <- 1
  response <- response + rnorm(length(response), sd = 0.05)
  fingerprint <- matrix(
    response,
    ncol = 1,
    dimnames = list(genes, "perturbation")
  )

  fit <- prn_build(
    fingerprint,
    prior,
    min_targets = 10,
    lambda = 0.01
  )
  result <- prn_counterdirection(fit)

  cancel_score <- mean(result$counterdirection[result$target %in% cancel])
  other_score <- mean(result$counterdirection[!result$target %in% cancel])

  expect_gt(cancel_score, 1)
  expect_equal(other_score, 0, tolerance = 1e-8)
  expect_true(is.finite(fit$fits$perturbation$intercept))
  expect_identical(colnames(fit$weight_matrix), c("A", "B"))
  expect_true(all(is.finite(prn_edge_table(fit)$inferred_contribution)))
})
