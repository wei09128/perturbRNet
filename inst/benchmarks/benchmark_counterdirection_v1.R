#!/usr/bin/env Rscript

# ============================================================
# perturbRNet Benchmark 1: counterdirection recovery
# ============================================================
# Usage:
#   Rscript benchmark_counterdirection_v1.R quick
#   Rscript benchmark_counterdirection_v1.R full
#
# Tests whether PRN distinguishes deliberately counterbalanced targets
# from inactive-looking targets as the experiment and prior deteriorate.
# ============================================================

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1]) else "quick"
if (!mode %in% c("quick", "full")) {
  stop("Mode must be 'quick' or 'full'.")
}

output_dir <- paste0("benchmark_counterdirection_", mode)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seed_base <- 20260811L
n_genes <- 800L
n_regulators <- 12L
n_counter_targets <- 80L

# -----------------------------
# Evaluation metrics
# -----------------------------

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

# -----------------------------
# Known signed response network
# -----------------------------

make_truth <- function(overlap_fraction, cancellation_strength, seed) {
  set.seed(seed)

  genes <- paste0("G", seq_len(n_genes))
  regulators <- paste0("R", seq_len(n_regulators))
  counter_targets <- genes[seq_len(n_counter_targets)]
  available <- setdiff(genes, counter_targets)

  # R1/R2 share increasing numbers of same-sign targets. This makes their
  # regulons increasingly collinear, while the counter targets retain
  # opposite edge signs.
  shared_n <- round(180 * overlap_fraction)
  shared_targets <- if (shared_n > 0) sample(available, shared_n) else character()
  remaining <- setdiff(available, shared_targets)
  unique_n <- max(30L, 180L - shared_n)
  r1_unique <- sample(remaining, unique_n)
  r2_unique <- sample(setdiff(remaining, r1_unique), unique_n)

  prior <- rbind(
    data.frame(source = "R1", target = counter_targets, sign = 1),
    data.frame(source = "R2", target = counter_targets, sign = -1),
    data.frame(source = "R1", target = shared_targets, sign = 1),
    data.frame(source = "R2", target = shared_targets, sign = 1),
    data.frame(source = "R1", target = r1_unique, sign = 1),
    data.frame(source = "R2", target = r2_unique, sign = 1)
  )

  # Other regulators have positive activities and positive edges. They can
  # overlap one another but cannot create counterdirection by construction.
  common_pool_n <- round(220 * overlap_fraction)
  common_pool <- if (common_pool_n > 0) sample(available, common_pool_n) else character()
  for (r in regulators[-c(1, 2)]) {
    common_take <- min(length(common_pool), round(100 * overlap_fraction))
    shared <- if (common_take > 0) sample(common_pool, common_take) else character()
    unique <- sample(setdiff(available, shared), 100L - length(shared))
    prior <- rbind(
      prior,
      data.frame(source = r, target = c(shared, unique), sign = 1)
    )
  }

  prior <- unique(prior)

  W <- matrix(
    0,
    nrow = length(genes),
    ncol = length(regulators),
    dimnames = list(genes, regulators)
  )
  W[cbind(match(prior$target, genes), match(prior$source, regulators))] <- prior$sign

  activity <- rep(0.65, n_regulators)
  names(activity) <- regulators
  activity[c("R1", "R2")] <- cancellation_strength

  contribution <- sweep(W, 2, activity, `*`)
  positive <- rowSums(pmax(contribution, 0))
  negative <- rowSums(pmax(-contribution, 0))
  true_net <- positive - negative
  true_gross <- positive + negative
  true_counter <- true_gross - abs(true_net)

  list(
    genes = genes,
    regulators = regulators,
    prior = prior,
    activity = activity,
    net_response = true_net,
    true_counter = true_counter,
    counter_targets = counter_targets
  )
}

corrupt_prior <- function(prior, drop_fraction, flip_fraction, seed) {
  set.seed(seed)
  out <- prior

  if (drop_fraction > 0) {
    keep <- stats::runif(nrow(out)) >= drop_fraction
    out <- out[keep, , drop = FALSE]
  }
  if (flip_fraction > 0 && nrow(out)) {
    flip <- stats::runif(nrow(out)) < flip_fraction
    out$sign[flip] <- -out$sign[flip]
  }
  out
}

run_one <- function(
    donor_n,
    noise_sd,
    drop_fraction,
    flip_fraction,
    overlap_fraction,
    cancellation_strength,
    replicate_id) {

  condition_id <- paste(
    donor_n, noise_sd, drop_fraction, flip_fraction,
    overlap_fraction, cancellation_strength, replicate_id,
    sep = "_"
  )
  seed <- seed_base + sum(utf8ToInt(condition_id)) + replicate_id * 1009L

  truth <- make_truth(overlap_fraction, cancellation_strength, seed)
  observed_sd <- noise_sd * sqrt(3 / donor_n)
  set.seed(seed + 1L)
  observed <- truth$net_response + rnorm(n_genes, sd = observed_sd)
  fingerprint <- matrix(
    observed,
    ncol = 1,
    dimnames = list(truth$genes, "perturbation")
  )

  observed_prior <- corrupt_prior(
    truth$prior,
    drop_fraction,
    flip_fraction,
    seed + 2L
  )

  fit <- tryCatch(
    prn_build(
      fingerprint = fingerprint,
      prior = observed_prior,
      min_targets = 10L,
      lambda = NULL,
      evidence_level = "simulation_corrupted_prior"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      donor_n, noise_sd, drop_fraction, flip_fraction,
      overlap_fraction, cancellation_strength, replicate_id,
      status = "failed", error = conditionMessage(fit),
      auroc = NA_real_, auprc = NA_real_, top_k_accuracy = NA_real_,
      counter_spearman = NA_real_, regulator_spearman = NA_real_,
      target_coverage = NA_real_, counter_target_coverage = NA_real_,
      mean_positive_score = NA_real_, mean_negative_score = NA_real_,
      lambda = NA_real_, variance_explained = NA_real_
    ))
  }

  inferred <- prn_counterdirection(fit)
  inferred <- inferred[match(truth$genes, inferred$target), ]
  label <- truth$genes %in% truth$counter_targets

  activity_est <- fit$fits[["perturbation"]]$activity
  names(activity_est) <- colnames(fit$weight_matrix)
  common_regulators <- intersect(names(truth$activity), names(activity_est))
  true_column_norm <- sqrt(table(truth$prior$source))
  expected_normalized_activity <-
    truth$activity[names(true_column_norm)] * true_column_norm
  regulator_cor <- if (length(common_regulators) >= 3) {
    suppressWarnings(stats::cor(
      expected_normalized_activity[common_regulators],
      activity_est[common_regulators],
      method = "spearman",
      use = "complete.obs"
    ))
  } else {
    NA_real_
  }

  diagnostics <- prn_summary(fit)

  data.frame(
    donor_n, noise_sd, drop_fraction, flip_fraction,
    overlap_fraction, cancellation_strength, replicate_id,
    status = "ok", error = "",
    auroc = auc_roc(label, inferred$counterdirection),
    auprc = auc_pr(label, inferred$counterdirection),
    top_k_accuracy = top_k_accuracy(label, inferred$counterdirection),
    counter_spearman = suppressWarnings(stats::cor(
      truth$true_counter,
      inferred$counterdirection,
      method = "spearman",
      use = "complete.obs"
    )),
    regulator_spearman = regulator_cor,
    target_coverage = mean(is.finite(inferred$counterdirection)),
    counter_target_coverage = mean(is.finite(inferred$counterdirection[label])),
    mean_positive_score = mean(inferred$counterdirection[label]),
    mean_negative_score = mean(inferred$counterdirection[!label]),
    lambda = diagnostics$lambda,
    variance_explained = diagnostics$variance_explained
  )
}

# -----------------------------
# Benchmark grid
# -----------------------------

if (mode == "quick") {
  grid <- expand.grid(
    donor_n = c(3, 8),
    noise_sd = c(0.25, 0.75),
    drop_fraction = c(0, 0.25),
    flip_fraction = c(0, 0.15),
    overlap_fraction = c(0.10, 0.60),
    cancellation_strength = c(0.50, 1.00),
    replicate_id = seq_len(5),
    KEEP.OUT.ATTRS = FALSE
  )
} else {
  grid <- expand.grid(
    donor_n = c(3, 6, 12),
    noise_sd = c(0.10, 0.35, 0.70, 1.00),
    drop_fraction = c(0, 0.15, 0.35),
    flip_fraction = c(0, 0.10, 0.25),
    overlap_fraction = c(0.05, 0.35, 0.70),
    cancellation_strength = c(0.25, 0.60, 1.00),
    replicate_id = seq_len(10),
    KEEP.OUT.ATTRS = FALSE
  )
}

message("Mode: ", mode)
message("Benchmark fits: ", nrow(grid))

results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  if (i == 1L || i %% 50L == 0L || i == nrow(grid)) {
    message("Fit ", i, " / ", nrow(grid))
  }
  results[[i]] <- do.call(run_one, as.list(grid[i, ]))
}
results <- do.call(rbind, results)

write.csv(
  results,
  file.path(output_dir, "benchmark_raw_results.csv"),
  row.names = FALSE
)

# Marginal summaries: average each factor level across all other factors.
metric_names <- c(
  "auroc", "auprc", "top_k_accuracy", "counter_spearman",
  "regulator_spearman", "target_coverage", "counter_target_coverage",
  "variance_explained"
)
factor_names <- c(
  "donor_n", "noise_sd", "drop_fraction", "flip_fraction",
  "overlap_fraction", "cancellation_strength"
)

ok_results <- results[results$status == "ok", , drop = FALSE]
marginal <- do.call(rbind, lapply(factor_names, function(factor_name) {
  pieces <- lapply(metric_names, function(metric_name) {
    aggregate(
      ok_results[[metric_name]],
      by = list(level = ok_results[[factor_name]]),
      FUN = function(x) c(
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        n = sum(is.finite(x))
      )
    ) -> tmp
    values <- if (is.matrix(tmp$x)) tmp$x else do.call(rbind, tmp$x)
    data.frame(
      factor = factor_name,
      level = tmp$level,
      metric = metric_name,
      mean = values[, "mean"],
      sd = values[, "sd"],
      n = values[, "n"]
    )
  })
  do.call(rbind, pieces)
}))

write.csv(
  marginal,
  file.path(output_dir, "benchmark_marginal_summary.csv"),
  row.names = FALSE
)

# Compact base-R diagnostic figure; no visualization package required.
pdf(file.path(output_dir, "benchmark_overview.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

boxplot(auprc ~ noise_sd, data = ok_results,
        xlab = "Noise SD", ylab = "AUPRC", main = "Noise robustness")
abline(h = n_counter_targets / n_genes, lty = 2, col = "grey50")

boxplot(auprc ~ drop_fraction, data = ok_results,
        xlab = "Dropped prior edges", ylab = "AUPRC",
        main = "Missing-prior robustness")

boxplot(auprc ~ flip_fraction, data = ok_results,
        xlab = "Sign-flipped prior edges", ylab = "AUPRC",
        main = "Wrong-sign robustness")

boxplot(auprc ~ overlap_fraction, data = ok_results,
        xlab = "Regulon overlap", ylab = "AUPRC",
        main = "Collinearity robustness")

dev.off()

summary_lines <- c(
  paste0("perturbRNet counterdirection benchmark: ", mode),
  paste0("Total fits: ", nrow(results)),
  paste0("Successful fits: ", sum(results$status == "ok")),
  paste0("Failed fits: ", sum(results$status != "ok")),
  paste0("Mean AUROC: ", round(mean(ok_results$auroc, na.rm = TRUE), 4)),
  paste0("Mean AUPRC: ", round(mean(ok_results$auprc, na.rm = TRUE), 4)),
  paste0("Random AUPRC baseline: ", round(n_counter_targets / n_genes, 4)),
  paste0("Mean top-k accuracy: ", round(mean(ok_results$top_k_accuracy, na.rm = TRUE), 4)),
  paste0("Mean counterdirection Spearman: ",
         round(mean(ok_results$counter_spearman, na.rm = TRUE), 4)),
  paste0("Mean regulator-activity Spearman: ",
         round(mean(ok_results$regulator_spearman, na.rm = TRUE), 4))
)

writeLines(summary_lines, file.path(output_dir, "BENCHMARK_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
