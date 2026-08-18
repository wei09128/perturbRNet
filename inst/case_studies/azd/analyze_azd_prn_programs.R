#!/usr/bin/env Rscript

# Paper 1 bounded biological-program analysis for perturbRNet restoration states.
#
# Inputs:
#   restoration_classification_all_genes.csv
#
# Outputs:
#   program_state_counts.csv
#   matched_pathway_enrichment_all.csv
#   top_programs.csv
#   candidate_targets_descriptive.csv
#   program_enrichment_heatmap.{png,pdf}
#   PROGRAM_ANALYSIS_SUMMARY.txt
#
# The formal test uses the full PRN-covered universe and samples null target
# sets matched exactly on treatment-specific incoming-degree and expression
# strata. Sets with fewer than 10 mapped targets are descriptive only.

suppressPackageStartupMessages({
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required. Install with install.packages('msigdbr').")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
  }
})

args <- commandArgs(trailingOnly = TRUE)

option <- function(name, default = NULL) {
  hit <- match(name, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", name)
  args[[hit + 1L]]
}

input_file <- option(
  "--input",
  "/mnt/c/Wei/AZD/brain/prn_analysis/restoration_classification_final/restoration_classification_all_genes.csv"
)
output_directory <- option(
  "--out",
  "/mnt/c/Wei/AZD/brain/prn_analysis/restoration_program_analysis"
)
permutations <- as.integer(option("--permutations", "10000"))
seed <- as.integer(option("--seed", "20260828"))
minimum_formal_targets <- as.integer(option("--minimum-formal-targets", "10"))
degree_bins <- as.integer(option("--degree-bins", "5"))
expression_bins <- as.integer(option("--expression-bins", "5"))

if (!file.exists(input_file)) stop("Input not found: ", input_file)
if (!is.finite(permutations) || permutations < 100L) stop("Use at least 100 permutations.")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

x <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

required_common <- c("SYMBOL", "logCPM")
missing_common <- setdiff(required_common, names(x))
if (length(missing_common)) {
  stop("Missing required columns: ", paste(missing_common, collapse = ", "))
}

qbin <- function(values, bins) {
  values <- as.numeric(values)
  result <- rep(NA_integer_, length(values))
  ok <- is.finite(values)
  if (!any(ok)) return(result)
  ranks <- rank(values[ok], ties.method = "average")
  result[ok] <- pmin(bins, pmax(1L, ceiling(ranks / sum(ok) * bins)))
  result
}

get_msig <- function(collection, subcollection) {
  formal <- names(formals(msigdbr::msigdbr))
  collection_arg <- if ("collection" %in% formal) "collection" else "category"
  subcollection_arg <- if ("subcollection" %in% formal) "subcollection" else "subcategory"

  call_msig <- function(selected_collection, native_mouse) {
    call_args <- list(species = "Mus musculus")
    if (native_mouse && "db_species" %in% formal) call_args$db_species <- "MM"
    call_args[[collection_arg]] <- selected_collection
    call_args[[subcollection_arg]] <- subcollection
    do.call(msigdbr::msigdbr, call_args)
  }

  # Mouse-native MSigDB prefixes parallel the human collections:
  # human C5 ontology -> mouse M5; human C2 curated -> mouse M2.
  native_collection <- unname(c(C5 = "M5", C2 = "M2")[[collection]])
  if ("db_species" %in% formal && !is.null(native_collection)) {
    native_result <- tryCatch(
      call_msig(native_collection, native_mouse = TRUE),
      error = function(e) e
    )
    if (!inherits(native_result, "error")) return(native_result)
    warning(
      "Mouse-native collection ", native_collection, "/", subcollection,
      " was unavailable; falling back to human ", collection,
      " with mouse ortholog mapping. Details: ", conditionMessage(native_result)
    )
  }

  call_msig(collection, native_mouse = FALSE)
}

message("Loading mouse GO Biological Process and Reactome gene sets...")
go <- get_msig("C5", "GO:BP")
go$database <- "GO_BP"
reactome <- tryCatch(
  {
    z <- get_msig("C2", "CP:REACTOME")
    z$database <- "REACTOME"
    z
  },
  error = function(e) {
    warning("Reactome collection unavailable in this msigdbr version: ", conditionMessage(e))
    NULL
  }
)
gene_sets <- rbind(go, reactome)

entrez_column <- intersect(c("ncbi_gene", "entrez_gene"), names(gene_sets))[1]
if (is.na(entrez_column)) stop("msigdbr output does not contain an Entrez identifier column.")
gene_sets$entrez <- as.character(gene_sets[[entrez_column]])
gene_sets$pathway <- as.character(gene_sets$gs_name)
gene_sets <- unique(gene_sets[c("database", "pathway", "entrez")])
gene_sets <- gene_sets[nzchar(gene_sets$entrez) & !is.na(gene_sets$entrez), ]

treatment_specs <- list(
  AZD = list(
    id = "azd_final_target",
    degree = "azd_final_incoming_degree",
    restored = "azd_expression_restored",
    calibrated = "azd_calibrated",
    direct = "azd_direct_restoration",
    strong = "azd_strong_counterbalanced_restoration",
    low = "azd_low_amplitude_counterbalanced_restoration"
  ),
  shATM = list(
    id = "shatm_final_target",
    degree = "shatm_final_incoming_degree",
    restored = "shatm_expression_restored",
    calibrated = "shatm_calibrated",
    direct = "shatm_direct_restoration",
    strong = "shatm_strong_counterbalanced_restoration",
    low = "shatm_low_amplitude_counterbalanced_restoration"
  )
)

all_results <- list()
count_rows <- list()
candidate_rows <- list()

matched_enrichment <- function(universe, selected, treatment, state, B, seed_value) {
  universe <- universe[!duplicated(universe$entrez), ]
  selected_ids <- intersect(unique(selected), universe$entrez)
  selected_index <- match(selected_ids, universe$entrez)
  selected_index <- selected_index[!is.na(selected_index)]

  mapping <- gene_sets[gene_sets$entrez %in% universe$entrez, ]
  sizes <- table(interaction(mapping$database, mapping$pathway, drop = TRUE))
  keep_keys <- names(sizes)[sizes >= 10L & sizes <= 500L]
  mapping$key <- interaction(mapping$database, mapping$pathway, drop = TRUE)
  mapping <- mapping[as.character(mapping$key) %in% keep_keys, ]
  mapping <- unique(mapping[c("database", "pathway", "entrez")])

  pathways <- unique(mapping[c("database", "pathway")])
  pathways$key <- paste(pathways$database, pathways$pathway, sep = "\r")
  mapping$key2 <- paste(mapping$database, mapping$pathway, sep = "\r")
  gene_index <- match(mapping$entrez, universe$entrez)
  pathway_index <- match(mapping$key2, pathways$key)
  incidence <- Matrix::sparseMatrix(
    i = gene_index,
    j = pathway_index,
    x = 1,
    dims = c(nrow(universe), nrow(pathways))
  )
  incidence[incidence > 0] <- 1

  observed <- as.numeric(Matrix::colSums(incidence[selected_index, , drop = FALSE]))
  pathway_size <- as.numeric(Matrix::colSums(incidence))
  formal <- length(selected_index) >= minimum_formal_targets

  result <- data.frame(
    treatment = treatment,
    restoration_state = state,
    formal_test = formal,
    selected_targets = length(selected_index),
    database = pathways$database,
    pathway = pathways$pathway,
    pathway_size_in_universe = pathway_size,
    observed_overlap = observed,
    null_mean_overlap = NA_real_,
    null_sd_overlap = NA_real_,
    enrichment_ratio = NA_real_,
    empirical_p_upper = NA_real_,
    BH_FDR = NA_real_,
    BH_FDR = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!formal) return(result)

  strata <- split(seq_len(nrow(universe)), universe$stratum)
  selected_strata <- table(universe$stratum[selected_index])
  selected_strata <- selected_strata[selected_strata > 0]
  pool_sizes <- vapply(
    names(selected_strata),
    function(stratum_name) length(strata[[stratum_name]]),
    integer(1)
  )
  missing_pool <- names(selected_strata)[
    pool_sizes < as.integer(selected_strata)
  ]
  if (length(missing_pool)) stop("Insufficient matched pool in strata: ", paste(missing_pool, collapse = ", "))

  set.seed(seed_value)
  exceed <- numeric(nrow(pathways))
  sum_null <- numeric(nrow(pathways))
  sumsq_null <- numeric(nrow(pathways))
  for (b in seq_len(B)) {
    draw <- unlist(lapply(names(selected_strata), function(s) {
      sample(strata[[s]], size = selected_strata[[s]], replace = FALSE)
    }), use.names = FALSE)
    overlap <- as.numeric(Matrix::colSums(incidence[draw, , drop = FALSE]))
    exceed <- exceed + (overlap >= observed)
    sum_null <- sum_null + overlap
    sumsq_null <- sumsq_null + overlap^2
    if (b %% max(100L, floor(B / 10L)) == 0L) {
      message("  ", treatment, " / ", state, ": ", b, " / ", B)
    }
  }
  null_mean <- sum_null / B
  null_var <- pmax(0, (sumsq_null - B * null_mean^2) / max(1, B - 1L))
  result$null_mean_overlap <- null_mean
  result$null_sd_overlap <- sqrt(null_var)
  result$enrichment_ratio <- ifelse(null_mean > 0, observed / null_mean, NA_real_)
  result$empirical_p_upper <- (exceed + 1) / (B + 1)
  result$BH_FDR <- p.adjust(result$empirical_p_upper, method = "BH")
  result
}

for (treatment in names(treatment_specs)) {
  spec <- treatment_specs[[treatment]]
  missing <- setdiff(unlist(spec), names(x))
  if (length(missing)) stop(treatment, " columns missing: ", paste(missing, collapse = ", "))

  z <- data.frame(
    entrez = as.character(x[[spec$id]]),
    symbol = as.character(x$SYMBOL),
    logCPM = as.numeric(x$logCPM),
    degree = as.numeric(x[[spec$degree]]),
    restored = as.logical(x[[spec$restored]]),
    calibrated = as.logical(x[[spec$calibrated]]),
    direct = as.logical(x[[spec$direct]]),
    strong = as.logical(x[[spec$strong]]),
    low = as.logical(x[[spec$low]]),
    stringsAsFactors = FALSE
  )
  for (nm in c("restored", "calibrated", "direct", "strong", "low")) z[[nm]][is.na(z[[nm]])] <- FALSE
  z <- z[nzchar(z$entrez) & !is.na(z$entrez) & is.finite(z$degree) & is.finite(z$logCPM), ]
  z <- z[!duplicated(z$entrez), ]
  z$degree_bin <- qbin(log1p(z$degree), degree_bins)
  z$expression_bin <- qbin(z$logCPM, expression_bins)
  z$stratum <- paste(z$degree_bin, z$expression_bin, sep = "_")

  states <- list(
    direct_restoration = z$entrez[z$restored & z$direct],
    mechanism_unresolved = z$entrez[z$restored & !z$direct & !z$strong & !z$low],
    calibrated_restoration = z$entrez[z$restored & z$calibrated],
    counterbalanced_candidates = z$entrez[z$restored & (z$strong | z$low)]
  )

  count_rows[[treatment]] <- data.frame(
    treatment = treatment,
    restoration_state = names(states),
    targets = vapply(states, length, integer(1)),
    formal_test = vapply(states, length, integer(1)) >= minimum_formal_targets,
    stringsAsFactors = FALSE
  )

  candidates <- z[z$restored & (z$strong | z$low), c("entrez", "symbol", "strong", "low")]
  candidates$treatment <- treatment
  candidates$interpretation <- ifelse(candidates$strong, "strong candidate", "low-amplitude exploratory")
  candidate_rows[[treatment]] <- candidates[c("treatment", "entrez", "symbol", "interpretation")]

  for (i in seq_along(states)) {
    state <- names(states)[i]
    message("Testing ", treatment, " / ", state, " (n=", length(states[[i]]), ")")
    key <- paste(treatment, state, sep = "__")
    all_results[[key]] <- matched_enrichment(
      universe = z,
      selected = states[[i]],
      treatment = treatment,
      state = state,
      B = permutations,
      seed_value = seed + match(treatment, names(treatment_specs)) * 100L + i
    )
  }
}

counts <- do.call(rbind, count_rows)
candidates <- do.call(rbind, candidate_rows)
results <- do.call(rbind, all_results)
row.names(counts) <- row.names(candidates) <- row.names(results) <- NULL

write.csv(counts, file.path(output_directory, "program_state_counts.csv"), row.names = FALSE)
write.csv(candidates, file.path(output_directory, "candidate_targets_descriptive.csv"), row.names = FALSE)
write.csv(results, file.path(output_directory, "matched_pathway_enrichment_all.csv"), row.names = FALSE)

reportable <- results[
  results$formal_test & results$observed_overlap >= 3 & is.finite(results$empirical_p_upper),
]
reportable <- reportable[order(
  reportable$treatment,
  reportable$restoration_state,
  reportable$BH_FDR,
  -reportable$enrichment_ratio,
  reportable$pathway
), ]
top <- do.call(rbind, lapply(split(reportable, interaction(reportable$treatment, reportable$restoration_state)), head, 15))
if (is.null(top)) top <- reportable[FALSE, ]
row.names(top) <- NULL
write.csv(top, file.path(output_directory, "top_programs.csv"), row.names = FALSE)

if (nrow(top)) {
  plot_data <- top
  plot_data$label <- gsub("^(GOBP_|REACTOME_)", "", plot_data$pathway)
  plot_data$label <- gsub("_", " ", plot_data$label)
  plot_data$label <- factor(plot_data$label, levels = rev(unique(plot_data$label)))
  plot_data$score <- -log10(pmax(plot_data$empirical_p_upper, 1 / (permutations + 1)))
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = score, y = label, fill = enrichment_ratio)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_grid(treatment ~ restoration_state, scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey80") +
    ggplot2::labs(
      x = expression(-log[10]("matched empirical P")),
      y = NULL,
      fill = "Observed /\nnull overlap",
      title = "Biological programs associated with perturbRNet restoration states",
      subtitle = "Null sets matched on incoming degree and expression"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right"
    )
  ggplot2::ggsave(file.path(output_directory, "program_enrichment_heatmap.png"), p, width = 11, height = 8, dpi = 400)
  ggplot2::ggsave(file.path(output_directory, "program_enrichment_heatmap.pdf"), p, width = 11, height = 8)
}

significant <- reportable[is.finite(reportable$BH_FDR) & reportable$BH_FDR <= 0.10, ]
summary_lines <- c(
  "perturbRNet restoration-state program analysis",
  paste0("Input: ", normalizePath(input_file)),
  paste0("Permutations per formal state: ", permutations),
  paste0("Matched strata: ", degree_bins, " degree bins x ", expression_bins, " expression bins"),
  paste0("Minimum targets for formal enrichment: ", minimum_formal_targets),
  "",
  "State counts:",
  paste(capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
  "",
  paste0("Reportable pathway rows (observed overlap >=3): ", nrow(reportable)),
  paste0("Pathways at BH FDR <=0.10: ", nrow(significant)),
  "",
  "Interpretation boundary:",
  "Candidate sets below the minimum target count are descriptive only. Enrichment identifies state-associated programs, not causal pathways or uniquely identifiable regulators."
)
writeLines(summary_lines, file.path(output_directory, "PROGRAM_ANALYSIS_SUMMARY.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nOutput directory: ", normalizePath(output_directory), "\n", sep = "")
