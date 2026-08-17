#!/usr/bin/env Rscript

# perturbRNet Benchmark 4.1: identifiability versus detectability
# Usage: Rscript benchmark_identifiability_v4.1.R quick|full

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1L]) else "quick"
if (!mode %in% c("quick", "full")) stop("Mode must be quick or full.")

output_dir <- paste0("benchmark_identifiability_v4.1_", mode)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
seed_base <- 20260813L
n_genes <- 400L
n_counter <- 40L

auc_pr <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep]); score <- score[keep]
  if (!length(score) || !sum(truth)) return(NA_real_)
  order_index <- order(score, decreasing = TRUE)
  truth <- truth[order_index]
  precision <- cumsum(truth) / seq_along(truth)
  mean(precision[truth])
}

auc_roc <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep]); score <- score[keep]
  positive_n <- sum(truth); negative_n <- sum(!truth)
  if (!positive_n || !negative_n) return(NA_real_)
  score_rank <- rank(score, ties.method = "average")
  (sum(score_rank[truth]) - positive_n * (positive_n + 1) / 2) /
    (positive_n * negative_n)
}

simulate_one <- function(
    unique_targets,
    cancellation_strength,
    noise_sd,
    missing_fraction,
    sign_error,
    replicate_id) {

  scenario_text <- paste(
    unique_targets, cancellation_strength, noise_sd,
    missing_fraction, sign_error, replicate_id, sep = "_"
  )
  seed <- seed_base + sum(utf8ToInt(scenario_text)) + 1009L * replicate_id
  set.seed(seed)

  genes <- paste0("G", seq_len(n_genes))
  counter_targets <- genes[seq_len(n_counter)]
  quiet_opposing_targets <- genes[n_counter + seq_len(n_counter)]
  available <- setdiff(genes, c(counter_targets, quiet_opposing_targets))

  a_unique <- if (unique_targets) sample(available, unique_targets) else character()
  available <- setdiff(available, a_unique)
  b_unique <- if (unique_targets) sample(available, unique_targets) else character()
  available <- setdiff(available, b_unique)
  c_unique <- if (unique_targets) sample(available, unique_targets) else character()
  available <- setdiff(available, c_unique)
  d_unique <- if (unique_targets) sample(available, unique_targets) else character()
  available <- setdiff(available, d_unique)

  prior <- rbind(
    data.frame(source = "A", target = counter_targets, sign = 1),
    data.frame(source = "B", target = counter_targets, sign = -1),
    data.frame(source = "A", target = a_unique, sign = 1),
    data.frame(source = "B", target = b_unique, sign = 1),
    # C and D have the same opposing topology, but zero true activity. These
    # are the hard negatives: 0 + 0 rather than active +s + (-s).
    data.frame(source = "C", target = quiet_opposing_targets, sign = 1),
    data.frame(source = "D", target = quiet_opposing_targets, sign = -1),
    data.frame(source = "C", target = c_unique, sign = 1),
    data.frame(source = "D", target = d_unique, sign = 1)
  )

  # Inactive regulators provide many quiet decoys within the fitted network.
  for (regulator in paste0("E", 1:4)) {
    decoy_n <- min(40L, length(available))
    decoy <- sample(available, decoy_n)
    available <- setdiff(available, decoy)
    prior <- rbind(
      prior,
      data.frame(source = regulator, target = decoy, sign = sample(c(-1, 1), decoy_n, TRUE))
    )
  }
  rownames(prior) <- NULL

  true_prior <- prior
  if (sign_error > 0) {
    flip_n <- round(nrow(prior) * sign_error)
    if (flip_n) {
      flip <- sample(seq_len(nrow(prior)), flip_n)
      prior$sign[flip] <- -prior$sign[flip]
    }
  }

  regulators <- unique(true_prior$source)
  W <- matrix(0, nrow = n_genes, ncol = length(regulators),
              dimnames = list(genes, regulators))
  W[cbind(match(true_prior$target, genes), match(true_prior$source, regulators))] <-
    true_prior$sign
  activity <- stats::setNames(rep(0, length(regulators)), regulators)
  activity[c("A", "B")] <- cancellation_strength
  true_response <- drop(W %*% activity)
  observed_response <- true_response + stats::rnorm(n_genes, 0, noise_sd)

  missing_n <- round(n_genes * missing_fraction)
  retained <- if (missing_n) setdiff(genes, sample(genes, missing_n)) else genes
  fingerprint <- matrix(
    observed_response[match(retained, genes)],
    ncol = 1,
    dimnames = list(retained, "perturbation")
  )

  fit <- prn_build(fingerprint, prior, min_targets = 1, lambda = 0.1,
                   evidence_level = "simulation_known_prior")
  counter <- prn_counterdirection(fit)
  score <- stats::setNames(counter$counterdirection, counter$target)
  score <- score[match(genes, names(score))]
  truth <- genes %in% counter_targets
  diagnostic <- prn_identifiability(fit)

  ab <- diagnostic$pair[
    (diagnostic$pair$regulator_1 == "A" & diagnostic$pair$regulator_2 == "B") |
      (diagnostic$pair$regulator_1 == "B" & diagnostic$pair$regulator_2 == "A"),
    , drop = FALSE
  ]
  regulator_ab <- diagnostic$regulator[
    diagnostic$regulator$regulator %in% c("A", "B"), , drop = FALSE
  ]
  counter_id <- diagnostic$target[
    diagnostic$target$target %in% counter_targets &
      diagnostic$target$comparison == "perturbation", , drop = FALSE
  ]

  data.frame(
    unique_targets = unique_targets,
    cancellation_strength = cancellation_strength,
    noise_sd = noise_sd,
    missing_fraction = missing_fraction,
    sign_error = sign_error,
    replicate_id = replicate_id,
    auroc = auc_roc(truth, score),
    auprc = auc_pr(truth, score),
    counter_spearman = suppressWarnings(stats::cor(
      as.numeric(truth), score, method = "spearman", use = "complete.obs"
    )),
    pair_separation = if (nrow(ab)) ab$pair_separation[1L] else NA_real_,
    minimum_regulator_uniqueness = if (nrow(regulator_ab)) {
      min(regulator_ab$unique_information, na.rm = TRUE)
    } else NA_real_,
    rank_fraction = diagnostic$global$rank_fraction,
    effective_rank_fraction = diagnostic$global$effective_rank_fraction,
    median_counter_target_separation = if (
      nrow(counter_id) && any(is.finite(counter_id$weighted_pair_separation))
    ) stats::median(
      counter_id$weighted_pair_separation[is.finite(counter_id$weighted_pair_separation)]
    ) else NA_real_,
    target_coverage = mean(is.finite(score)),
    counter_target_coverage = mean(is.finite(score[truth])),
    error = NA_character_
  )
}

if (mode == "quick") {
  grid <- expand.grid(
    unique_targets = c(1, 2, 8, 30),
    cancellation_strength = c(0.05, 0.15, 0.50),
    noise_sd = c(0.05, 0.25),
    missing_fraction = c(0, 0.30),
    sign_error = c(0, 0.05),
    replicate_id = seq_len(5),
    KEEP.OUT.ATTRS = FALSE
  )
} else {
  grid <- expand.grid(
    unique_targets = c(1, 2, 5, 10, 20, 40, 60),
    cancellation_strength = c(0.02, 0.05, 0.10, 0.20, 0.50),
    noise_sd = c(0.02, 0.05, 0.15, 0.30, 0.60),
    missing_fraction = c(0, 0.15, 0.30, 0.50),
    sign_error = c(0, 0.01, 0.05, 0.15),
    replicate_id = seq_len(10),
    KEEP.OUT.ATTRS = FALSE
  )
}

message("Mode: ", mode)
message("Scenarios: ", nrow(grid))
pieces <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  if (i == 1L || i %% 50L == 0L || i == nrow(grid)) {
    message("Scenario ", i, " / ", nrow(grid))
  }
  pieces[[i]] <- tryCatch(
    do.call(simulate_one, as.list(grid[i, ])),
    error = function(e) data.frame(
      grid[i, ], auroc = NA_real_, auprc = NA_real_, counter_spearman = NA_real_,
      pair_separation = NA_real_, minimum_regulator_uniqueness = NA_real_,
      rank_fraction = NA_real_, effective_rank_fraction = NA_real_,
      median_counter_target_separation = NA_real_, target_coverage = NA_real_,
      counter_target_coverage = NA_real_, error = conditionMessage(e)
    )
  )
}
results <- do.call(rbind, pieces)
write.csv(results, file.path(output_dir, "benchmark4_raw_results.csv"), row.names = FALSE)

valid <- results[is.finite(results$auprc), , drop = FALSE]
metric_names <- c(
  "pair_separation", "minimum_regulator_uniqueness", "rank_fraction",
  "effective_rank_fraction", "median_counter_target_separation",
  "cancellation_strength", "noise_sd", "missing_fraction", "sign_error"
)
association <- do.call(rbind, lapply(metric_names, function(metric) {
  keep <- is.finite(valid[[metric]]) & is.finite(valid$auprc)
  data.frame(
    predictor = metric,
    spearman_with_auprc = suppressWarnings(stats::cor(
      valid[[metric]][keep], valid$auprc[keep], method = "spearman"
    )),
    n = sum(keep)
  )
}))
write.csv(association, file.path(output_dir, "benchmark4_associations.csv"), row.names = FALSE)

by_geometry_signal <- aggregate(
  valid$auprc,
  by = list(
    unique_targets = valid$unique_targets,
    cancellation_strength = valid$cancellation_strength
  ),
  FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x))
)
values <- if (is.matrix(by_geometry_signal$x)) {
  by_geometry_signal$x
} else do.call(rbind, by_geometry_signal$x)
by_geometry_signal <- data.frame(
  by_geometry_signal[c("unique_targets", "cancellation_strength")],
  mean_auprc = values[, "mean"], sd_auprc = values[, "sd"], n = values[, "n"]
)
write.csv(by_geometry_signal,
          file.path(output_dir, "benchmark4_geometry_by_signal.csv"), row.names = FALSE)

joint <- stats::lm(
  auprc ~ pair_separation + cancellation_strength + noise_sd +
    missing_fraction + sign_error,
  data = valid
)
capture.output(summary(joint),
               file = file.path(output_dir, "benchmark4_joint_model.txt"))

summary_lines <- c(
  paste0("perturbRNet Benchmark 4.1: ", mode),
  paste0("Scenarios: ", nrow(grid)),
  paste0("Successful scenarios: ", nrow(valid)),
  paste0("Failed scenarios: ", nrow(grid) - nrow(valid)),
  paste0("Mean AUPRC: ", round(mean(valid$auprc), 4)),
  paste0("Random AUPRC baseline: ", n_counter / n_genes),
  "",
  "Spearman associations with AUPRC:",
  capture.output(print(association, row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "BENCHMARK4_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
