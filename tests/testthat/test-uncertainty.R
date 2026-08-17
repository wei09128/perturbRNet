test_that("certain signs reproduce a stable ensemble", {
  set.seed(8)
  genes <- paste0("g", seq_len(120))
  cancel <- genes[1:20]
  prior <- rbind(
    data.frame(source = "A", target = cancel, sign = 1, confidence = 1),
    data.frame(source = "B", target = cancel, sign = -1, confidence = 1),
    data.frame(source = "A", target = genes[21:70], sign = 1, confidence = 1),
    data.frame(source = "B", target = genes[71:120], sign = 1, confidence = 1)
  )
  response <- setNames(rep(0, length(genes)), genes)
  response[21:120] <- 1
  response <- response + rnorm(length(response), sd = 0.03)
  fingerprint <- matrix(response, ncol = 1,
                        dimnames = list(genes, "perturbation"))

  fit <- prn_build_uncertain(
    fingerprint,
    prior,
    n_draws = 10,
    min_targets = 5,
    lambda = 0.01,
    seed = 99
  )
  summary <- fit$target_summary
  cancel_summary <- summary[summary$target %in% cancel, ]

  expect_s3_class(fit, "prn_uncertain")
  expect_true(all(cancel_summary$mean_counterdirection > 1))
  expect_true(all(cancel_summary$sd_counterdirection < 1e-10))
  expect_true(all(cancel_summary$probability_counterdirection == 1))
  expect_true(all(is.finite(fit$activity_summary$mean_activity)))
  expect_equal(
    fit$target_summary$mean_counterdirection,
    fit$target_summary$supported_mean_counterdirection,
    tolerance = 1e-10
  )
  expect_equal(
    aggregate(model_support_weight ~ comparison, fit$draw_support, sum)$model_support_weight,
    1,
    tolerance = 1e-10
  )
})

test_that("uncertain signs produce measurable instability", {
  genes <- paste0("g", seq_len(120))
  cancel <- genes[1:20]
  prior <- rbind(
    data.frame(source = "A", target = cancel, sign = 1, confidence = 0.55),
    data.frame(source = "B", target = cancel, sign = -1, confidence = 0.55),
    data.frame(source = "A", target = genes[21:70], sign = 1, confidence = 0.55),
    data.frame(source = "B", target = genes[71:120], sign = 1, confidence = 0.55)
  )
  response <- setNames(rep(0, length(genes)), genes)
  response[21:120] <- 1
  fingerprint <- matrix(response, ncol = 1,
                        dimnames = list(genes, "perturbation"))

  fit <- prn_build_uncertain(
    fingerprint,
    prior,
    n_draws = 20,
    min_targets = 5,
    lambda = 0.1,
    seed = 7
  )

  expect_true(any(fit$target_summary$sd_counterdirection > 0))
  expect_true(all(fit$target_summary$top_candidate_stability >= 0))
  expect_true(all(fit$target_summary$top_candidate_stability <= 1))
  expect_true(all(fit$target_summary$supported_top_candidate_stability >= 0))
  expect_true(all(fit$target_summary$supported_top_candidate_stability <= 1))
})

test_that("expected-sign weighting is zero at confidence one half", {
  genes <- paste0("g", seq_len(80))
  prior <- rbind(
    data.frame(source = "A", target = genes[1:40], sign = 1, confidence = 0.5),
    data.frame(source = "B", target = genes[41:80], sign = 1, confidence = 1)
  )
  response <- matrix(
    seq(-1, 1, length.out = length(genes)),
    ncol = 1,
    dimnames = list(genes, "perturbation")
  )

  # Source A is removed from the deterministic expected-sign prior because
  # c=0.5 implies E[S]=0. Source B remains identifiable.
  fit <- prn_build_uncertain(
    response,
    prior,
    n_draws = 5,
    min_targets = 5,
    lambda = 0.1,
    seed = 5
  )

  expect_identical(colnames(fit$weighted_fit$weight_matrix), "B")
  expect_equal(sum(fit$draw_support$model_support_weight), 1)
})
