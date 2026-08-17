#!/usr/bin/env Rscript

# ============================================================
# perturbRNet Benchmark 3: cell counts to counterdirection
# ============================================================
# Usage:
#   Rscript benchmark_end_to_end_counts_v3.R quick
#   Rscript benchmark_end_to_end_counts_v3.R full
#
# Tests the complete analysis path:
#   single-cell counts -> donor pseudobulk -> paired edgeR contrast
#   -> perturbation fingerprint -> PRN -> counterdirection
# ============================================================

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1]) else "quick"
if (!mode %in% c("quick", "full")) stop("Mode must be quick or full.")

output_dir <- paste0("benchmark_end_to_end_", mode)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed_base <- 20260812L
n_genes <- 600L
n_regulators <- 10L
n_counter <- 60L

auc_roc <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok]); score <- score[ok]
  n1 <- sum(truth); n0 <- sum(!truth)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auc_pr <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok]); score <- score[ok]
  if (!sum(truth)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- truth[ord]
  precision <- cumsum(y) / seq_along(y)
  mean(precision[y])
}

top_k_accuracy <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.logical(truth[ok]); score <- score[ok]
  k <- sum(truth)
  if (!k) return(NA_real_)
  mean(truth[order(score, decreasing = TRUE)[seq_len(k)]])
}

make_network <- function(cancellation_strength, seed) {
  set.seed(seed)
  genes <- paste0("G", seq_len(n_genes))
  regulators <- paste0("R", seq_len(n_regulators))
  counter_targets <- genes[seq_len(n_counter)]
  available <- setdiff(genes, counter_targets)

  # R1 and R2 oppose one another at counter_targets but have distinct
  # reinforcing footprints elsewhere, which makes their activity estimable.
  r1_unique <- sample(available, 100)
  r2_unique <- sample(setdiff(available, r1_unique), 100)
  prior <- rbind(
    data.frame(source = "R1", target = counter_targets, sign = 1),
    data.frame(source = "R2", target = counter_targets, sign = -1),
    data.frame(source = "R1", target = r1_unique, sign = 1),
    data.frame(source = "R2", target = r2_unique, sign = 1)
  )

  # R3-R7 are active reinforcing programs. R8-R10 are inactive programs,
  # creating unchanged decoy targets that DE magnitude alone cannot separate
  # from actively counterbalanced genes.
  for (r in regulators[-c(1, 2)]) {
    targets <- sample(available, 80)
    prior <- rbind(prior, data.frame(source = r, target = targets, sign = 1))
  }
  prior <- unique(prior)

  W <- matrix(0, nrow = n_genes, ncol = n_regulators,
              dimnames = list(genes, regulators))
  W[cbind(match(prior$target, genes), match(prior$source, regulators))] <-
    prior$sign

  activity <- c(
    R1 = cancellation_strength,
    R2 = cancellation_strength,
    R3 = 0.20, R4 = 0.28, R5 = 0.36, R6 = 0.44, R7 = 0.52,
    R8 = 0, R9 = 0, R10 = 0
  )
  contribution <- sweep(W, 2, activity, `*`)
  positive <- rowSums(pmax(contribution, 0))
  negative <- rowSums(pmax(-contribution, 0))
  true_effect <- positive - negative
  true_counter <- positive + negative - abs(true_effect)

  list(
    genes = genes,
    prior = prior,
    true_effect = true_effect,
    true_counter = true_counter,
    counter_targets = counter_targets
  )
}

simulate_cells <- function(
    truth,
    donor_n,
    cells_per_group,
    biological_dispersion,
    seed) {

  set.seed(seed)
  genes <- truth$genes

  # Gene baseline abundance varies over roughly two orders of magnitude.
  baseline <- exp(stats::rnorm(n_genes, log(2.5), 0.9))
  baseline <- pmax(baseline, 0.15)

  # Donor-by-gene random effect is shared between control and treatment,
  # producing a genuinely paired design.
  donor_effect <- matrix(
    stats::rnorm(n_genes * donor_n, 0, biological_dispersion),
    nrow = n_genes,
    ncol = donor_n
  )

  metadata <- vector("list", donor_n * 2L)
  count_parts <- vector("list", donor_n * 2L)
  part <- 0L

  for (donor in seq_len(donor_n)) {
    for (condition in c("control", "treated")) {
      part <- part + 1L
      treatment_effect <- if (condition == "treated") truth$true_effect else 0
      gene_mean <- baseline *
        exp(donor_effect[, donor]) *
        (2 ^ treatment_effect)

      # Cell-specific library-size factors mimic sequencing-depth variation.
      library_factor <- stats::rlnorm(cells_per_group, 0, 0.25)
      mu <- outer(gene_mean, library_factor)

      # Negative-binomial technical/biological cell-count variation.
      counts <- matrix(
        stats::rnbinom(
          length(mu),
          mu = as.numeric(mu),
          size = 1 / 0.35
        ),
        nrow = n_genes,
        ncol = cells_per_group,
        dimnames = list(genes, NULL)
      )
      count_parts[[part]] <- counts
      metadata[[part]] <- data.frame(
        donor = paste0("D", donor),
        cell_type = "simulated_cell",
        condition = condition,
        stringsAsFactors = FALSE,
        row.names = NULL
      )[rep(1, cells_per_group), , drop = FALSE]
    }
  }

  counts <- do.call(cbind, count_parts)
  metadata <- do.call(rbind, metadata)
  rownames(metadata) <- NULL
  colnames(counts) <- paste0("C", seq_len(ncol(counts)))

  list(counts = counts, metadata = metadata)
}

evaluate_scores <- function(
    truth,
    score,
    method,
    scenario,
    estimated_effect = NULL) {

  score <- score[match(truth$genes, names(score))]
  label <- truth$genes %in% truth$counter_targets
  result <- data.frame(
    scenario,
    method = method,
    auroc = auc_roc(label, score),
    auprc = auc_pr(label, score),
    top_k_accuracy = top_k_accuracy(label, score),
    counter_spearman = suppressWarnings(stats::cor(
      truth$true_counter, score,
      method = "spearman", use = "complete.obs"
    )),
    target_coverage = mean(is.finite(score)),
    counter_target_coverage = mean(is.finite(score[label])),
    stringsAsFactors = FALSE
  )
  if (!is.null(estimated_effect)) {
    estimated_effect <- estimated_effect[
      match(truth$genes, names(estimated_effect))
    ]
    result$effect_spearman <- suppressWarnings(stats::cor(
      truth$true_effect, estimated_effect,
      method = "spearman", use = "complete.obs"
    ))
    result$effect_rmse <- sqrt(mean(
      (truth$true_effect - estimated_effect) ^ 2,
      na.rm = TRUE
    ))
  } else {
    result$effect_spearman <- NA_real_
    result$effect_rmse <- NA_real_
  }
  result
}

run_one <- function(
    donor_n,
    cells_per_group,
    biological_dispersion,
    cancellation_strength,
    replicate_id) {

  scenario_text <- paste(
    donor_n, cells_per_group, biological_dispersion,
    cancellation_strength, replicate_id,
    sep = "_"
  )
  seed <- seed_base + sum(utf8ToInt(scenario_text)) + replicate_id * 1009L
  truth <- make_network(cancellation_strength, seed)
  sim <- simulate_cells(
    truth, donor_n, cells_per_group, biological_dispersion, seed + 1L
  )

  pb <- prn_aggregate_pseudobulk(
    counts = sim$counts,
    metadata = sim$metadata,
    sample_col = "donor",
    cell_type_col = "cell_type",
    condition_cols = "condition",
    min_cells = 20L
  )
  pb$metadata$donor <- factor(pb$metadata$donor)
  pb$metadata$condition <- stats::relevel(
    factor(pb$metadata$condition), ref = "control"
  )

  de <- prn_fit_contrasts(
    counts = pb$counts,
    metadata = pb$metadata,
    design_formula = ~ donor + condition,
    contrasts = list(treated_vs_control = "conditiontreated")
  )
  fingerprint <- prn_make_fingerprint(de)

  estimated_fit <- prn_build(
    fingerprint,
    truth$prior,
    min_targets = 10,
    evidence_level = "estimated_pseudobulk_response"
  )
  estimated_target <- prn_counterdirection(estimated_fit)
  estimated_score <- setNames(
    estimated_target$counterdirection,
    estimated_target$target
  )

  # Oracle PRN receives the exact simulated log2 effect. Its performance is
  # the ceiling against which response-estimation loss is measured.
  oracle_fp <- matrix(
    truth$true_effect,
    ncol = 1,
    dimnames = list(truth$genes, "treated_vs_control")
  )
  oracle_fit <- prn_build(
    oracle_fp,
    truth$prior,
    min_targets = 10,
    evidence_level = "oracle_response"
  )
  oracle_target <- prn_counterdirection(oracle_fit)
  oracle_score <- setNames(oracle_target$counterdirection, oracle_target$target)

  estimated_effect <- setNames(
    de$effect[de$comparison == "treated_vs_control"],
    de$gene_id[de$comparison == "treated_vs_control"]
  )
  quiet_score <- -abs(estimated_effect)

  scenario <- list(
    donor_n = donor_n,
    cells_per_group = cells_per_group,
    biological_dispersion = biological_dispersion,
    cancellation_strength = cancellation_strength,
    replicate_id = replicate_id,
    pseudobulk_samples = ncol(pb$counts),
    tested_genes = nrow(fingerprint)
  )

  do.call(rbind, list(
    evaluate_scores(
      truth, oracle_score, "oracle_prn", scenario
    ),
    evaluate_scores(
      truth, estimated_score, "estimated_prn", scenario, estimated_effect
    ),
    evaluate_scores(
      truth, quiet_score, "small_abs_logfc", scenario, estimated_effect
    )
  ))
}

if (mode == "quick") {
  grid <- expand.grid(
    donor_n = c(3, 8),
    cells_per_group = c(50, 200),
    biological_dispersion = c(0.15, 0.50),
    cancellation_strength = c(0.30, 0.70),
    replicate_id = seq_len(3),
    KEEP.OUT.ATTRS = FALSE
  )
} else {
  grid <- expand.grid(
    donor_n = c(3, 6, 12),
    cells_per_group = c(30, 100, 300),
    biological_dispersion = c(0.10, 0.35, 0.70),
    cancellation_strength = c(0.20, 0.45, 0.80),
    replicate_id = seq_len(10),
    KEEP.OUT.ATTRS = FALSE
  )
}

message("Mode: ", mode)
message("Scenarios: ", nrow(grid))

pieces <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  if (i == 1L || i %% 10L == 0L || i == nrow(grid)) {
    message("Scenario ", i, " / ", nrow(grid))
  }
  pieces[[i]] <- tryCatch(
    do.call(run_one, as.list(grid[i, ])),
    error = function(e) data.frame(
      donor_n = grid$donor_n[i],
      cells_per_group = grid$cells_per_group[i],
      biological_dispersion = grid$biological_dispersion[i],
      cancellation_strength = grid$cancellation_strength[i],
      replicate_id = grid$replicate_id[i],
      pseudobulk_samples = NA_integer_, tested_genes = NA_integer_,
      method = "FAILED", auroc = NA_real_, auprc = NA_real_,
      top_k_accuracy = NA_real_, counter_spearman = NA_real_,
      target_coverage = NA_real_, counter_target_coverage = NA_real_,
      effect_spearman = NA_real_, effect_rmse = NA_real_,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  )
}
results <- do.call(rbind, pieces)

write.csv(results, file.path(output_dir, "benchmark3_raw_results.csv"),
          row.names = FALSE)
valid <- results[results$method != "FAILED", , drop = FALSE]

overall <- aggregate(
  valid[c("auroc", "auprc", "top_k_accuracy", "counter_spearman",
          "target_coverage", "counter_target_coverage",
          "effect_spearman", "effect_rmse")],
  by = list(method = valid$method),
  FUN = function(x) mean(x, na.rm = TRUE)
)
write.csv(overall, file.path(output_dir, "benchmark3_method_summary.csv"),
          row.names = FALSE)

factor_names <- c(
  "donor_n", "cells_per_group", "biological_dispersion",
  "cancellation_strength"
)
marginal <- do.call(rbind, lapply(factor_names, function(factor_name) {
  do.call(rbind, lapply(unique(valid$method), function(method_name) {
    dat <- valid[valid$method == method_name, ]
    tmp <- aggregate(
      dat$auprc,
      by = list(level = dat[[factor_name]]),
      FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x))
    )
    values <- if (is.matrix(tmp$x)) tmp$x else do.call(rbind, tmp$x)
    data.frame(
      factor = factor_name,
      level = tmp$level,
      method = method_name,
      mean_auprc = values[, "mean"],
      sd_auprc = values[, "sd"],
      n = values[, "n"]
    )
  }))
}))
write.csv(marginal, file.path(output_dir, "benchmark3_marginal_summary.csv"),
          row.names = FALSE)

pdf(file.path(output_dir, "benchmark3_overview.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (factor_name in factor_names) {
  dat <- marginal[marginal$factor == factor_name, ]
  methods <- unique(dat$method)
  plot(range(dat$level), range(dat$mean_auprc), type = "n",
       xlab = factor_name, ylab = "Mean AUPRC", main = factor_name)
  for (j in seq_along(methods)) {
    x <- dat[dat$method == methods[j], ]
    x <- x[order(x$level), ]
    lines(x$level, x$mean_auprc, type = "b", col = j, pch = j)
  }
  legend("bottomright", legend = methods, col = seq_along(methods),
         pch = seq_along(methods), lty = 1, cex = 0.7)
}
dev.off()

summary_lines <- c(
  paste0("perturbRNet Benchmark 3: ", mode),
  paste0("Scenarios: ", nrow(grid)),
  paste0("Successful scenarios: ", sum(results$method == "oracle_prn")),
  paste0("Failed scenarios: ", sum(results$method == "FAILED")),
  "",
  capture.output(print(overall, row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "BENCHMARK3_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
