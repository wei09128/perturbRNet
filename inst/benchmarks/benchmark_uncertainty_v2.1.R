#!/usr/bin/env Rscript

# ============================================================
# perturbRNet Benchmark 2.1: uncertainty-aware sign inference
# ============================================================
# Usage:
#   Rscript benchmark_uncertainty_v2.R quick
#   Rscript benchmark_uncertainty_v2.R full
#
# Each simulated dataset is analyzed using exactly the same observed response
# and reported prior by four scoring strategies:
#   1. binary reported signs;
#   2. expected-sign weights, s*(2c-1);
#   3. sign-ensemble mean counterdirection;
#   4. sign-ensemble top-candidate stability.
#   5. data-supported ensemble mean;
#   6. data-supported top-candidate stability.
#
# Confidence is assigned before sign errors are sampled. It controls the
# probability that a reported sign is correct, but never reveals which
# particular edges are wrong.
# ============================================================

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1]) else "quick"
if (!mode %in% c("quick", "full")) {
  stop("Mode must be 'quick' or 'full'.")
}

output_dir <- paste0("benchmark_uncertainty_", mode)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seed_base <- 20260812L
n_genes <- 800L
n_regulators <- 12L
n_counter_targets <- 80L
random_auprc <- n_counter_targets / n_genes

auc_roc <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok])
  score <- score[ok]
  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

auc_pr <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok])
  score <- score[ok]
  n_pos <- sum(truth)
  if (!n_pos) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  precision <- cumsum(y) / seq_along(y)
  mean(precision[y])
}

top_k_accuracy <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok])
  score <- score[ok]
  k <- sum(truth)
  if (!k) return(NA_real_)
  mean(truth[order(score, decreasing = TRUE)[seq_len(k)]])
}

make_truth <- function(overlap_fraction, cancellation_strength, seed) {
  set.seed(seed)
  genes <- paste0("G", seq_len(n_genes))
  regulators <- paste0("R", seq_len(n_regulators))
  counter_targets <- genes[seq_len(n_counter_targets)]
  available <- setdiff(genes, counter_targets)

  shared_n <- round(180 * overlap_fraction)
  shared_targets <- if (shared_n > 0) sample(available, shared_n) else character()
  remaining <- setdiff(available, shared_targets)
  unique_n <- max(30L, 180L - shared_n)
  r1_unique <- sample(remaining, unique_n)
  r2_unique <- sample(setdiff(remaining, r1_unique), unique_n)

  prior <- rbind(
    data.frame(source = "R1", target = counter_targets, true_sign = 1),
    data.frame(source = "R2", target = counter_targets, true_sign = -1),
    data.frame(source = "R1", target = shared_targets, true_sign = 1),
    data.frame(source = "R2", target = shared_targets, true_sign = 1),
    data.frame(source = "R1", target = r1_unique, true_sign = 1),
    data.frame(source = "R2", target = r2_unique, true_sign = 1)
  )

  common_pool_n <- round(220 * overlap_fraction)
  common_pool <- if (common_pool_n > 0) sample(available, common_pool_n) else character()
  for (r in regulators[-c(1, 2)]) {
    common_take <- min(length(common_pool), round(100 * overlap_fraction))
    shared <- if (common_take > 0) sample(common_pool, common_take) else character()
    unique_targets <- sample(setdiff(available, shared), 100L - length(shared))
    prior <- rbind(
      prior,
      data.frame(
        source = r,
        target = c(shared, unique_targets),
        true_sign = 1
      )
    )
  }
  prior <- unique(prior)

  W <- matrix(
    0,
    nrow = length(genes),
    ncol = length(regulators),
    dimnames = list(genes, regulators)
  )
  W[cbind(match(prior$target, genes), match(prior$source, regulators))] <-
    prior$true_sign

  activity <- c(
    R1 = cancellation_strength,
    R2 = cancellation_strength,
    stats::setNames(
      seq(0.30, 0.90, length.out = n_regulators - 2L),
      regulators[-c(1, 2)]
    )
  )
  contribution <- sweep(W, 2, activity, `*`)
  positive <- rowSums(pmax(contribution, 0))
  negative <- rowSums(pmax(-contribution, 0))
  net_response <- positive - negative
  true_counter <- positive + negative - abs(net_response)

  list(
    genes = genes,
    prior = prior,
    activity = activity,
    net_response = net_response,
    true_counter = true_counter,
    counter_targets = counter_targets
  )
}

make_reported_prior <- function(
    truth_prior,
    mean_error,
    calibration,
    seed,
    concentration = 40) {

  set.seed(seed)
  # Heterogeneous edge-specific error probabilities with the requested mean.
  alpha <- max(mean_error * concentration, 0.05)
  beta <- max((1 - mean_error) * concentration, 0.05)
  true_error_probability <- stats::rbeta(nrow(truth_prior), alpha, beta)
  is_wrong <- stats::runif(nrow(truth_prior)) < true_error_probability
  reported_sign <- truth_prior$true_sign * ifelse(is_wrong, -1, 1)

  reported_error_probability <- switch(
    calibration,
    calibrated = true_error_probability,
    overconfident = true_error_probability * 0.35,
    underconfident = pmin(0.49, true_error_probability * 1.75),
    stop("Unknown calibration setting.")
  )

  data.frame(
    source = truth_prior$source,
    target = truth_prior$target,
    sign = reported_sign,
    confidence = 1 - reported_error_probability,
    is_wrong = is_wrong,
    true_error_probability = true_error_probability,
    stringsAsFactors = FALSE
  )
}

score_method <- function(
    method,
    genes,
    truth_label,
    truth_counter,
    score,
    scenario_fields) {

  score <- score[match(genes, names(score))]
  ok <- is.finite(score)
  data.frame(
    scenario_fields,
    method = method,
    auroc = auc_roc(truth_label, score),
    auprc = auc_pr(truth_label, score),
    top_k_accuracy = top_k_accuracy(truth_label, score),
    counter_spearman = suppressWarnings(stats::cor(
      truth_counter,
      score,
      method = "spearman",
      use = "complete.obs"
    )),
    target_coverage = mean(ok),
    counter_target_coverage = mean(ok[truth_label]),
    stringsAsFactors = FALSE
  )
}

run_one <- function(
    mean_error,
    calibration,
    donor_n,
    noise_sd,
    overlap_fraction,
    cancellation_strength,
    replicate_id,
    n_draws) {

  scenario_text <- paste(
    mean_error, calibration, donor_n, noise_sd, overlap_fraction,
    cancellation_strength, replicate_id,
    sep = "_"
  )
  seed <- seed_base + sum(utf8ToInt(scenario_text)) + replicate_id * 1009L
  truth <- make_truth(overlap_fraction, cancellation_strength, seed)
  reported_prior <- make_reported_prior(
    truth$prior,
    mean_error,
    calibration,
    seed + 1L
  )

  set.seed(seed + 2L)
  observed_sd <- noise_sd * sqrt(3 / donor_n)
  observed <- truth$net_response + rnorm(n_genes, sd = observed_sd)
  fingerprint <- matrix(
    observed,
    ncol = 1,
    dimnames = list(truth$genes, "perturbation")
  )
  label <- truth$genes %in% truth$counter_targets

  scenario_fields <- list(
    mean_error = mean_error,
    realized_error = mean(reported_prior$is_wrong),
    calibration = calibration,
    donor_n = donor_n,
    noise_sd = noise_sd,
    overlap_fraction = overlap_fraction,
    cancellation_strength = cancellation_strength,
    replicate_id = replicate_id
  )

  binary_fit <- prn_build(
    fingerprint,
    reported_prior,
    source_col = "source",
    target_col = "target",
    sign_col = "sign",
    min_targets = 10,
    evidence_level = "reported_binary_sign"
  )
  binary_target <- prn_counterdirection(binary_fit)
  binary_score <- setNames(
    binary_target$counterdirection,
    binary_target$target
  )

  uncertain_fit <- prn_build_uncertain(
    fingerprint,
    reported_prior,
    source_col = "source",
    target_col = "target",
    sign_col = "sign",
    confidence_col = "confidence",
    n_draws = n_draws,
    min_targets = 10,
    counter_threshold = 0,
    top_fraction = random_auprc,
    seed = seed + 3L,
    keep_draws = FALSE
  )

  weighted_target <- prn_counterdirection(uncertain_fit$weighted_fit)
  weighted_score <- setNames(
    weighted_target$counterdirection,
    weighted_target$target
  )
  ensemble_mean <- setNames(
    uncertain_fit$target_summary$mean_counterdirection,
    uncertain_fit$target_summary$target
  )
  ensemble_stability <- setNames(
    uncertain_fit$target_summary$top_candidate_stability,
    uncertain_fit$target_summary$target
  )
  supported_mean <- setNames(
    uncertain_fit$target_summary$supported_mean_counterdirection,
    uncertain_fit$target_summary$target
  )
  supported_stability <- setNames(
    uncertain_fit$target_summary$supported_top_candidate_stability,
    uncertain_fit$target_summary$target
  )

  do.call(rbind, list(
    score_method(
      "binary", truth$genes, label, truth$true_counter,
      binary_score, scenario_fields
    ),
    score_method(
      "expected_sign", truth$genes, label, truth$true_counter,
      weighted_score, scenario_fields
    ),
    score_method(
      "ensemble_mean", truth$genes, label, truth$true_counter,
      ensemble_mean, scenario_fields
    ),
    score_method(
      "ensemble_stability", truth$genes, label, truth$true_counter,
      ensemble_stability, scenario_fields
    ),
    score_method(
      "supported_mean", truth$genes, label, truth$true_counter,
      supported_mean, scenario_fields
    ),
    score_method(
      "supported_stability", truth$genes, label, truth$true_counter,
      supported_stability, scenario_fields
    )
  ))
}

if (mode == "quick") {
  n_draws <- 20L
  grid <- expand.grid(
    mean_error = c(0.01, 0.05, 0.15),
    calibration = c("calibrated", "overconfident"),
    donor_n = c(3, 8),
    noise_sd = c(0.25, 0.75),
    overlap_fraction = c(0.10, 0.60),
    cancellation_strength = c(0.50, 1.00),
    replicate_id = seq_len(3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  n_draws <- 100L
  grid <- expand.grid(
    mean_error = c(0.005, 0.01, 0.05, 0.10, 0.15),
    calibration = c("calibrated", "overconfident", "underconfident"),
    donor_n = c(3, 6, 12),
    noise_sd = c(0.10, 0.40, 0.80),
    overlap_fraction = c(0.05, 0.35, 0.70),
    cancellation_strength = c(0.30, 0.60, 1.00),
    replicate_id = seq_len(10),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

message("Mode: ", mode)
message("Scenarios: ", nrow(grid))
message("Ensemble draws per scenario: ", n_draws)
message("Total PRN fits approximately: ", nrow(grid) * (n_draws + 2L))

result_list <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  if (i == 1L || i %% 10L == 0L || i == nrow(grid)) {
    message("Scenario ", i, " / ", nrow(grid))
  }
  result_list[[i]] <- tryCatch(
    do.call(run_one, c(as.list(grid[i, ]), list(n_draws = n_draws))),
    error = function(e) {
      data.frame(
        mean_error = grid$mean_error[i],
        realized_error = NA_real_,
        calibration = grid$calibration[i],
        donor_n = grid$donor_n[i],
        noise_sd = grid$noise_sd[i],
        overlap_fraction = grid$overlap_fraction[i],
        cancellation_strength = grid$cancellation_strength[i],
        replicate_id = grid$replicate_id[i],
        method = "FAILED",
        auroc = NA_real_,
        auprc = NA_real_,
        top_k_accuracy = NA_real_,
        counter_spearman = NA_real_,
        target_coverage = NA_real_,
        counter_target_coverage = NA_real_,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}
results <- do.call(rbind, result_list)

write.csv(
  results,
  file.path(output_dir, "benchmark2_raw_results.csv"),
  row.names = FALSE
)

ok <- results$method != "FAILED"
valid <- results[ok, , drop = FALSE]

summary_by_method <- aggregate(
  valid[c("auroc", "auprc", "top_k_accuracy", "counter_spearman")],
  by = valid[c("method", "mean_error", "calibration")],
  FUN = function(x) mean(x, na.rm = TRUE)
)
write.csv(
  summary_by_method,
  file.path(output_dir, "benchmark2_method_summary.csv"),
  row.names = FALSE
)

# Paired AUPRC improvement relative to binary within every scenario.
scenario_columns <- c(
  "mean_error", "calibration", "donor_n", "noise_sd",
  "overlap_fraction", "cancellation_strength", "replicate_id"
)
binary <- valid[valid$method == "binary", c(scenario_columns, "auprc")]
names(binary)[ncol(binary)] <- "binary_auprc"
paired <- merge(valid, binary, by = scenario_columns, all.x = TRUE)
paired$auprc_gain_vs_binary <- paired$auprc - paired$binary_auprc
write.csv(
  paired,
  file.path(output_dir, "benchmark2_paired_results.csv"),
  row.names = FALSE
)

gain_summary <- aggregate(
  paired$auprc_gain_vs_binary,
  by = paired[c("method", "mean_error", "calibration")],
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    win_rate = mean(x > 0, na.rm = TRUE),
    n = sum(is.finite(x))
  )
)
gain_values <- if (is.matrix(gain_summary$x)) {
  gain_summary$x
} else {
  do.call(rbind, gain_summary$x)
}
gain_summary <- cbind(
  gain_summary[c("method", "mean_error", "calibration")],
  as.data.frame(gain_values)
)
write.csv(
  gain_summary,
  file.path(output_dir, "benchmark2_gain_summary.csv"),
  row.names = FALSE
)

pdf(file.path(output_dir, "benchmark2_overview.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (cal in unique(valid$calibration)) {
  plot_data <- summary_by_method[summary_by_method$calibration == cal, ]
  methods <- unique(plot_data$method)
  colors <- seq_along(methods)
  plot(
    range(plot_data$mean_error),
    range(c(random_auprc, plot_data$auprc), na.rm = TRUE),
    type = "n",
    xlab = "Mean sign-error probability",
    ylab = "Mean AUPRC",
    main = paste("Calibration:", cal)
  )
  abline(h = random_auprc, lty = 2, col = "grey50")
  for (j in seq_along(methods)) {
    dat <- plot_data[plot_data$method == methods[j], ]
    dat <- dat[order(dat$mean_error), ]
    lines(dat$mean_error, dat$auprc, type = "b", col = colors[j], pch = j)
  }
  legend("topright", legend = methods, col = colors, pch = seq_along(methods),
         lty = 1, cex = 0.7)
}
boxplot(
  auprc_gain_vs_binary ~ method,
  data = paired[paired$method != "binary", ],
  las = 2,
  ylab = "AUPRC gain versus binary",
  main = "Paired method improvement"
)
abline(h = 0, lty = 2, col = "grey50")
plot.new()
text(0.5, 0.65, paste("Scenarios:", nrow(grid)))
text(0.5, 0.50, paste("Ensemble draws:", n_draws))
text(0.5, 0.35, paste("Failed scenarios:", sum(!ok)))
dev.off()

overall <- aggregate(
  valid[c("auroc", "auprc", "top_k_accuracy", "counter_spearman")],
  by = list(method = valid$method),
  FUN = function(x) mean(x, na.rm = TRUE)
)

summary_lines <- c(
  paste0("perturbRNet Benchmark 2: ", mode),
  paste0("Scenarios: ", nrow(grid)),
  paste0("Ensemble draws per scenario: ", n_draws),
  paste0("Failed scenarios: ", sum(!ok)),
  paste0("Random AUPRC baseline: ", random_auprc),
  "",
  capture.output(print(overall, row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "BENCHMARK2_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
