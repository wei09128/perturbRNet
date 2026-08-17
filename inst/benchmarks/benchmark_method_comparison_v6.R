#!/usr/bin/env Rscript

# perturbRNet method-comparison experiment
# Strict endpoint: counterbalanced versus inactive genes among quiet genes.
# Usage: Rscript benchmark_method_comparison_v6.R quick|full

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1L]) else "quick"
if (!mode %in% c("quick", "full")) stop("Mode must be quick or full.")
output_dir <- paste0("benchmark_method_comparison_", mode)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_per_state <- 50L
n_background <- 300L
seed_base <- 20260815L

auc_roc <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep]); score <- score[keep]
  n1 <- sum(truth); n0 <- sum(!truth)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auc_pr <- function(truth, score) {
  keep <- !is.na(truth) & is.finite(score)
  truth <- as.logical(truth[keep]); score <- score[keep]
  if (!length(score) || !sum(truth)) return(NA_real_)
  truth <- truth[order(score, decreasing = TRUE)]
  precision <- cumsum(truth) / seq_along(truth)
  mean(precision[truth])
}

counter_from_activity <- function(W, activity) {
  contribution <- sweep(W, 2L, activity, "*")
  positive <- rowSums(pmax(contribution, 0))
  negative <- rowSums(pmax(-contribution, 0))
  positive + negative - abs(positive - negative)
}

fit_univariate <- function(W, y) {
  denominator <- colSums(W^2)
  estimate <- as.numeric(crossprod(W, y)) / denominator
  estimate[!is.finite(estimate)] <- 0
  stats::setNames(estimate, colnames(W))
}

fit_ols <- function(W, y) {
  estimate <- tryCatch(
    qr.solve(cbind(Intercept = 1, W), y, tol = 1e-7)[-1L],
    error = function(e) rep(NA_real_, ncol(W))
  )
  stats::setNames(as.numeric(estimate), colnames(W))
}

make_scenario <- function(
    signal,
    noise_sd,
    unique_targets,
    overlap_fraction,
    missing_fraction,
    sign_error,
    replicate_id) {

  scenario_text <- paste(
    signal, noise_sd, unique_targets, overlap_fraction,
    missing_fraction, sign_error, replicate_id, sep = "_"
  )
  set.seed(seed_base + sum(utf8ToInt(scenario_text)) + 1009L * replicate_id)

  state <- rep(
    c("inactive", "positive", "negative", "counterbalanced", "background"),
    c(rep(n_per_state, 4L), n_background)
  )
  genes <- paste0("G", seq_along(state))
  names(state) <- genes
  regulators <- c("A", "B", "C", "D", paste0("N", 1:6))
  W_true <- matrix(0, nrow = length(genes), ncol = length(regulators),
                   dimnames = list(genes, regulators))

  counter <- genes[state == "counterbalanced"]
  inactive <- genes[state == "inactive"]
  positive <- genes[state == "positive"]
  negative <- genes[state == "negative"]
  background <- genes[state == "background"]

  # A and B are active. C and D are inactive topology-matched controls.
  W_true[counter, "A"] <- 1; W_true[counter, "B"] <- -1
  W_true[inactive, "C"] <- 1; W_true[inactive, "D"] <- -1
  W_true[positive, "A"] <- 1
  W_true[negative, "B"] <- 1

  available <- background
  for (regulator in regulators[1:4]) {
    targets <- sample(available, unique_targets)
    W_true[targets, regulator] <- 1
    available <- setdiff(available, targets)
  }
  for (regulator in regulators[5:10]) {
    targets <- sample(background, min(40L, length(background)))
    W_true[targets, regulator] <- sample(c(-1, 1), length(targets), TRUE)
  }

  # Add controlled footprint overlap without changing the four state labels.
  overlap_n <- round(length(background) * overlap_fraction)
  if (overlap_n > 0) {
    overlap_genes <- sample(background, overlap_n)
    W_true[overlap_genes, "A"] <- 1
    W_true[overlap_genes, "B"] <- sample(c(-1, 1), overlap_n, TRUE)
  }

  activity_true <- stats::setNames(rep(0, length(regulators)), regulators)
  activity_true[c("A", "B")] <- signal
  y_true <- drop(W_true %*% activity_true)
  y <- y_true + stats::rnorm(length(y_true), 0, noise_sd)

  W_observed <- W_true
  nonzero <- which(W_observed != 0, arr.ind = TRUE)
  if (missing_fraction > 0 && nrow(nonzero)) {
    drop_n <- round(nrow(nonzero) * missing_fraction)
    if (drop_n) {
      selected <- nonzero[sample(seq_len(nrow(nonzero)), drop_n), , drop = FALSE]
      W_observed[selected] <- 0
    }
  }
  nonzero <- which(W_observed != 0, arr.ind = TRUE)
  if (sign_error > 0 && nrow(nonzero)) {
    flip_n <- round(nrow(nonzero) * sign_error)
    if (flip_n) {
      selected <- nonzero[sample(seq_len(nrow(nonzero)), flip_n), , drop = FALSE]
      W_observed[selected] <- -W_observed[selected]
    }
  }

  keep_regulator <- colSums(W_observed != 0) > 0
  W_fit <- W_observed[, keep_regulator, drop = FALSE]
  observed_norm <- sqrt(colSums(W_fit^2))
  W_fit <- sweep(W_fit, 2L, observed_norm, "/")
  # If W is normalized to W / ||W_r||, the coefficient representing the same
  # per-edge activity is a_r * ||W_r||.
  activity_reference <- activity_true[colnames(W_fit)] * observed_norm

  prior_index <- which(W_observed != 0, arr.ind = TRUE)
  prior <- data.frame(
    source = colnames(W_observed)[prior_index[, 2L]],
    target = rownames(W_observed)[prior_index[, 1L]],
    sign = W_observed[prior_index],
    stringsAsFactors = FALSE
  )
  fingerprint <- matrix(y, ncol = 1, dimnames = list(genes, "perturbation"))
  ridge_fit <- prn_build(fingerprint, prior, min_targets = 1L)
  ridge <- prn_counterdirection(ridge_fit)
  ridge_score <- stats::setNames(ridge$counterdirection, ridge$target)

  activities <- list(
    univariate = fit_univariate(W_fit, y),
    multivariate_ols = fit_ols(W_fit, y),
    oracle_activity = activity_reference
  )
  scores <- list(
    abs_logfc_quietness = -abs(y),
    ridge_prn = ridge_score[genes]
  )
  for (method in names(activities)) {
    scores[[paste0(method, "_plus_counterdirection")]] <-
      counter_from_activity(W_fit, activities[[method]])
  }

  quiet <- state %in% c("inactive", "counterbalanced")
  truth <- state[quiet] == "counterbalanced"
  scenario <- data.frame(
    signal = signal, noise_sd = noise_sd, unique_targets = unique_targets,
    overlap_fraction = overlap_fraction, missing_fraction = missing_fraction,
    sign_error = sign_error, replicate_id = replicate_id
  )
  result <- do.call(rbind, lapply(names(scores), function(method) {
    score <- unname(scores[[method]])[quiet]
    data.frame(
      scenario,
      method = method,
      auroc = auc_roc(truth, score),
      auprc = auc_pr(truth, score),
      stringsAsFactors = FALSE
    )
  }))

  activity_result <- do.call(rbind, lapply(names(activities), function(method) {
    estimate <- activities[[method]]
    common <- intersect(names(estimate), names(activity_reference))
    data.frame(
      scenario,
      method = method,
      activity_spearman = suppressWarnings(stats::cor(
        activity_reference[common], estimate[common], method = "spearman"
      )),
      activity_rmse = sqrt(mean((activity_reference[common] - estimate[common])^2))
    )
  }))
  list(classification = result, activity = activity_result)
}

if (mode == "quick") {
  grid <- expand.grid(
    signal = c(0.10, 0.50), noise_sd = c(0.10, 0.40),
    unique_targets = c(5, 30), overlap_fraction = c(0, 0.20),
    missing_fraction = c(0, 0.25), sign_error = c(0, 0.05),
    replicate_id = seq_len(4), KEEP.OUT.ATTRS = FALSE
  )
} else {
  grid <- expand.grid(
    signal = c(0.05, 0.10, 0.25, 0.50, 1),
    noise_sd = c(0.05, 0.10, 0.25, 0.50),
    unique_targets = c(2, 5, 15, 30, 60),
    overlap_fraction = c(0, 0.10, 0.30),
    missing_fraction = c(0, 0.15, 0.35),
    sign_error = c(0, 0.01, 0.05, 0.15),
    replicate_id = seq_len(10), KEEP.OUT.ATTRS = FALSE
  )
}

message("Mode: ", mode)
message("Scenarios: ", nrow(grid))
pieces <- vector("list", nrow(grid)); failures <- list()
for (i in seq_len(nrow(grid))) {
  if (i == 1L || i %% 25L == 0L || i == nrow(grid)) {
    message("Scenario ", i, " / ", nrow(grid))
  }
  pieces[[i]] <- tryCatch(
    do.call(make_scenario, as.list(grid[i, ])),
    error = function(e) {
      failures[[length(failures) + 1L]] <<- data.frame(
        grid[i, ], error = conditionMessage(e)
      )
      NULL
    }
  )
}
pieces <- pieces[!vapply(pieces, is.null, logical(1))]
classification <- do.call(rbind, lapply(pieces, `[[`, "classification"))
activity <- do.call(rbind, lapply(pieces, `[[`, "activity"))
write.csv(classification, file.path(output_dir, "method_comparison_classification.csv"),
          row.names = FALSE)
write.csv(activity, file.path(output_dir, "method_comparison_activity.csv"),
          row.names = FALSE)
if (length(failures)) write.csv(do.call(rbind, failures),
                                file.path(output_dir, "failures.csv"), row.names = FALSE)

summary_table <- aggregate(
  classification[c("auroc", "auprc")],
  by = list(method = classification$method),
  FUN = function(x) mean(x[is.finite(x)])
)
write.csv(summary_table, file.path(output_dir, "method_comparison_summary.csv"),
          row.names = FALSE)
summary_lines <- c(
  paste0("perturbRNet method comparison: ", mode),
  paste0("Scenarios: ", nrow(grid)),
  paste0("Successful scenarios: ", length(pieces)),
  paste0("Failed scenarios: ", length(failures)),
  "Primary endpoint: counterbalanced versus inactive among quiet genes",
  "Random AUPRC baseline: 0.5",
  "",
  capture.output(print(summary_table, row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "METHOD_COMPARISON_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
