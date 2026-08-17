#!/usr/bin/env Rscript

# perturbRNet Benchmark 5: real-data destroyed-information controls
#
# Usage:
#   Rscript benchmark_realdata_nulls_v5.R fingerprint.csv prior.csv output_dir 100
#
# fingerprint.csv: first column = gene; remaining numeric columns = comparisons
# prior.csv: columns source, target, sign

suppressPackageStartupMessages(library(perturbRNet))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: benchmark_realdata_nulls_v5.R fingerprint.csv prior.csv output_dir [n_null]")
}
fingerprint_file <- args[1L]
prior_file <- args[2L]
output_dir <- args[3L]
n_null <- if (length(args) >= 4L) as.integer(args[4L]) else 100L
if (!is.finite(n_null) || n_null < 10L) stop("n_null must be at least 10.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_fingerprint <- function(path) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(dat) < 2L) stop("Fingerprint needs a gene column and >=1 comparison.")
  genes <- as.character(dat[[1L]])
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Fingerprint gene names must be nonempty and unique.")
  }
  x <- as.matrix(dat[-1L])
  storage.mode(x) <- "double"
  rownames(x) <- genes
  if (is.null(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("Fingerprint comparison columns require names.")
  }
  x
}

fingerprint <- read_fingerprint(fingerprint_file)
prior <- read.csv(prior_file, stringsAsFactors = FALSE)
required <- c("source", "target", "sign")
if (!all(required %in% names(prior))) {
  stop("Prior lacks columns: ", paste(setdiff(required, names(prior)), collapse = ", "))
}
prior <- prior[is.finite(as.numeric(prior$sign)) & as.numeric(prior$sign) != 0, required]
prior$source <- as.character(prior$source)
prior$target <- as.character(prior$target)
prior$sign <- as.numeric(prior$sign)
prior <- prior[prior$target %in% rownames(fingerprint), , drop = FALSE]
if (!nrow(prior)) stop("No prior targets overlap fingerprint genes.")

summarize_fit <- function(fit, null_type, draw_id) {
  counter <- prn_counterdirection(fit)
  diagnostic <- prn_identifiability(fit)
  cvalue <- counter$counterdirection[is.finite(counter$counterdirection)]
  balance <- counter$balance_fraction[is.finite(counter$balance_fraction)]
  fit_summary <- prn_summary(fit)
  variance <- fit_summary$variance_explained[is.finite(fit_summary$variance_explained)]
  data.frame(
    null_type = null_type,
    draw_id = draw_id,
    counterdirection_mean = if (length(cvalue)) mean(cvalue) else NA_real_,
    counterdirection_energy = if (length(cvalue)) mean(cvalue^2) else NA_real_,
    counterdirection_q95 = if (length(cvalue)) {
      unname(stats::quantile(cvalue, 0.95, names = FALSE))
    } else NA_real_,
    balanced_fraction = if (length(balance)) mean(balance >= 0.8) else NA_real_,
    variance_explained = if (length(variance)) mean(variance) else NA_real_,
    rank_fraction = diagnostic$global$rank_fraction,
    effective_rank_fraction = diagnostic$global$effective_rank_fraction,
    stringsAsFactors = FALSE
  )
}

fit_once <- function(fp, net, null_type, draw_id) {
  tryCatch({
    fit <- prn_build(fp, net, min_targets = 10L,
                     evidence_level = paste0("realdata_", null_type))
    list(summary = summarize_fit(fit, null_type, draw_id), error = NULL)
  }, error = function(e) list(summary = NULL, error = conditionMessage(e)))
}

permute_response <- function(fp) {
  out <- fp
  for (j in seq_len(ncol(out))) out[, j] <- sample(out[, j])
  out
}

shuffle_signs <- function(net) {
  out <- net
  out$sign <- ave(
    out$sign,
    out$source,
    FUN = function(x) x[sample.int(length(x))]
  )
  out
}

rewire_targets <- function(net, genes) {
  pieces <- split(net, net$source)
  pieces <- lapply(pieces, function(x) {
    # Preserve regulator degree and its sign vector while destroying targets.
    x$target <- sample(genes, nrow(x), replace = nrow(x) > length(genes))
    x
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

set.seed(20260814L)
records <- list()
errors <- list()
index <- 1L

observed <- fit_once(fingerprint, prior, "observed", 0L)
if (is.null(observed$summary)) stop("Observed fit failed: ", observed$error)
records[[index]] <- observed$summary
index <- index + 1L

null_types <- c("response_permutation", "sign_shuffle", "target_rewire")
for (draw in seq_len(n_null)) {
  if (draw == 1L || draw %% 10L == 0L || draw == n_null) {
    message("Null draw ", draw, " / ", n_null)
  }
  candidates <- list(
    response_permutation = list(fp = permute_response(fingerprint), net = prior),
    sign_shuffle = list(fp = fingerprint, net = shuffle_signs(prior)),
    target_rewire = list(fp = fingerprint, net = rewire_targets(prior, rownames(fingerprint)))
  )
  for (type in null_types) {
    result <- fit_once(candidates[[type]]$fp, candidates[[type]]$net, type, draw)
    if (is.null(result$summary)) {
      errors[[length(errors) + 1L]] <- data.frame(
        null_type = type, draw_id = draw, error = result$error
      )
    } else {
      records[[index]] <- result$summary
      index <- index + 1L
    }
  }
}

raw <- do.call(rbind, records)
write.csv(raw, file.path(output_dir, "benchmark5_raw_results.csv"), row.names = FALSE)
if (length(errors)) {
  write.csv(do.call(rbind, errors), file.path(output_dir, "benchmark5_errors.csv"),
            row.names = FALSE)
}

metrics <- c(
  "counterdirection_mean", "counterdirection_energy", "counterdirection_q95",
  "balanced_fraction", "variance_explained"
)
observed_row <- raw[raw$null_type == "observed", , drop = FALSE]
comparison <- do.call(rbind, lapply(null_types, function(type) {
  null <- raw[raw$null_type == type, , drop = FALSE]
  do.call(rbind, lapply(metrics, function(metric) {
    values <- null[[metric]][is.finite(null[[metric]])]
    observed_value <- observed_row[[metric]][1L]
    data.frame(
      null_type = type,
      metric = metric,
      observed = observed_value,
      null_mean = if (length(values)) mean(values) else NA_real_,
      null_sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      empirical_p_upper = if (length(values) && is.finite(observed_value)) {
        (1 + sum(values >= observed_value)) / (1 + length(values))
      } else NA_real_,
      standardized_difference = if (
        length(values) > 1L && is.finite(observed_value) && stats::sd(values) > 0
      ) (observed_value - mean(values)) / stats::sd(values) else NA_real_,
      successful_nulls = length(values)
    )
  }))
}))
write.csv(comparison, file.path(output_dir, "benchmark5_null_comparison.csv"),
          row.names = FALSE)

summary_lines <- c(
  "perturbRNet Benchmark 5: real-data destroyed-information controls",
  paste0("Genes: ", nrow(fingerprint)),
  paste0("Comparisons: ", ncol(fingerprint)),
  paste0("Prior edges retained: ", nrow(prior)),
  paste0("Null draws requested per type: ", n_null),
  paste0("Failed null fits: ", length(errors)),
  "",
  capture.output(print(comparison, row.names = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "BENCHMARK5_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
