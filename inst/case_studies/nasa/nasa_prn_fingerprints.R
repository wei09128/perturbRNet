#!/usr/bin/env Rscript

# ============================================================
# NASA GCR brain RNA-seq: PRN fingerprint construction
# ============================================================
# Two cohorts are two INDEPENDENT experiments (different tissue,
# different age, different timepoint, different exposure-route
# availability). Each gets its OWN DGEList, its OWN dispersion
# estimate, its OWN glmQLFit, and its OWN contrast set -- exactly
# mirroring azd_prn_fingerprints.R's discipline, applied twice.
# They are NEVER combined into one design matrix, and their
# fingerprints are written to two SEPARATE files.
#
# Young (August): cortex, shielded-only, 4 groups (ctrl/LD/MD/HD)
#   -> all 6 pairwise contrasts (same logic as AZD's 4-group set)
# Old (November): whole brain, Rad1(direct)/Rad2(shielded), MD/HD
#   only, sharing one Control -> 5 groups -> 8 well-motivated
#   pairwise contrasts (drops the 2 that cross both dose AND
#   route simultaneously, which aren't a clean single question)
#
# Metadata are read from the shared samples.csv and subset by
# Tissue inside run_cohort().
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
})

# ------------------------------------------------------------
# 0. Paths -- EDIT THESE to match your actual files
# ------------------------------------------------------------

young_count_file <- "/mnt/c/Wei/NASA/Novogene/row_counts_aug.csv"
old_count_file   <- "/mnt/c/Wei/NASA/Novogene/row_counts_nov.csv"

meta_file <- "/mnt/c/Wei/NASA/pipeline_config/samples.csv"

out_dir <- "/mnt/c/Wei/NASA/prn_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annotation_columns <- c("ID", "ENSEMBL", "ENTREZID", "GENENAME", "SYMBOL")

# ------------------------------------------------------------
# Shared helper: one cohort's fingerprint, start to finish.
# Called once per cohort with that cohort's own inputs. Nothing
# here is shared across calls -- each call gets a fresh DGEList,
# fresh dispersion estimate, fresh model fit.
# ------------------------------------------------------------

run_cohort <- function(
    cohort_label,
    count_file,
    meta_file,
    tissue_value,
    group_levels,
    build_group,
    contrast_defs,
    out_dir) {

  cat("\n============================================================\n")
  cat("Cohort:", cohort_label, "\n")
  cat("============================================================\n")

  raw  <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)

  stopifnot(all(annotation_columns %in% names(raw)))
  stopifnot(all(c("Sample", "Tissue", "Group") %in% names(meta)))

  meta <- meta[
    meta$Tissue == tissue_value &
      !is.na(meta$Group) &
      nzchar(trimws(meta$Group)),
    , drop = FALSE
  ]
  meta$Sample <- trimws(meta$Sample)
  meta$Group <- trimws(meta$Group)

  if (!nrow(meta)) {
    stop(cohort_label, ": no usable metadata rows for Tissue = ", tissue_value)
  }
  stopifnot(!anyDuplicated(meta$Sample))

  cat("Count file:", count_file, "\n")
  cat("Metadata file:", meta_file, "\n")
  cat("Metadata tissue:", tissue_value, "\n")
  cat("Raw dimensions:", nrow(raw), "rows x", ncol(raw), "columns\n")
  cat("Retained metadata samples:", nrow(meta), "\n")

  # Keep only samples this cohort's metadata actually lists.
  # This is where any unexplained columns (e.g. stray M01-M12
  # entries not present in samples.csv) get silently excluded --
  # printed below so it is never silent in practice.
  available <- intersect(meta$Sample, names(raw))
  missing_from_counts <- setdiff(meta$Sample, names(raw))
  extra_in_counts <- setdiff(
    setdiff(names(raw), annotation_columns), meta$Sample
  )
  if (length(missing_from_counts)) {
    stop(
      cohort_label, ": samples.csv lists samples absent from counts: ",
      paste(missing_from_counts, collapse = ", ")
    )
  }
  if (length(extra_in_counts)) {
    cat(
      "NOTE:", cohort_label, "count file has columns not in samples.csv",
      "(excluded from this fit):\n  ",
      paste(extra_in_counts, collapse = ", "), "\n"
    )
  }

  meta <- meta[match(available, meta$Sample), , drop = FALSE]
  annotations <- raw[, annotation_columns, drop = FALSE]
  counts <- data.matrix(raw[, meta$Sample, drop = FALSE])

  stopifnot(
    identical(colnames(counts), meta$Sample),
    all(is.finite(counts)),
    all(counts >= 0)
  )

  gene_ids <- as.character(annotations$ENSEMBL)
  bad_gene_id <- is.na(gene_ids) | !nzchar(trimws(gene_ids))
  gene_ids[bad_gene_id] <- as.character(annotations$ID[bad_gene_id])
  gene_ids <- make.unique(gene_ids)
  if (length(gene_ids) != nrow(counts)) {
    stop(
      cohort_label, ": gene_ids has length ", length(gene_ids),
      " but counts has ", nrow(counts), " rows -- the counts file's row ",
      "count no longer matches its own annotation columns. This usually ",
      "means the counts CSV was hand-edited (e.g. a column deleted with a ",
      "text editor) rather than a sample simply being removed from ",
      "samples.csv. Re-check ", count_file, " directly (dim() should be ",
      "the same as when it was first inspected) rather than editing it in ",
      "place -- drop unwanted samples from samples.csv only."
    )
  }
  rownames(counts) <- gene_ids
  rownames(annotations) <- gene_ids

  # ---------------- design ----------------
  group_values <- build_group(meta)
  if (is.null(group_values) || length(group_values) != nrow(meta)) {
    stop(
      cohort_label, ": group builder returned ", length(group_values),
      " values for ", nrow(meta), " samples."
    )
  }
  group_values <- as.character(group_values)
  unexpected_groups <- setdiff(unique(group_values), group_levels)
  if (length(unexpected_groups)) {
    stop(
      cohort_label, ": unexpected groups: ",
      paste(unexpected_groups, collapse = ", "),
      ". Expected: ", paste(group_levels, collapse = ", ")
    )
  }
  group <- factor(group_values, levels = group_levels)
  if (anyNA(group)) stop(cohort_label, ": missing group assignments.")

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  rownames(design) <- meta$Sample

  cat("Group sizes:\n")
  print(table(group))
  cat("\nDesign rank:", qr(design)$rank, "/", ncol(design), "\n")
  stopifnot(qr(design)$rank == ncol(design))

  write.csv(
    cbind(meta, design),
    file.path(out_dir, paste0(cohort_label, "_design_matrix.csv")),
    row.names = FALSE
  )

  # ---------------- edgeR model (fit ONCE per cohort) ----------------
  dge <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(dge, design = design)
  cat("\nGenes before filtering:", nrow(dge), "\n")
  cat("Genes after filtering: ", sum(keep), "\n")

  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)

  pdf(
    file.path(out_dir, paste0(cohort_label, "_edgeR_qc.pdf")),
    width = 8, height = 7
  )
  plotMDS(
    dge, labels = meta$Sample, col = as.integer(group),
    main = paste0(cohort_label, ": MDS")
  )
  legend(
    "topright", legend = levels(group),
    col = seq_along(levels(group)), pch = 16, bty = "n"
  )
  plotBCV(dge, main = paste0(cohort_label, ": biological coefficient of variation"))
  dev.off()

  # ---------------- contrasts ----------------
  contrast_matrix <- do.call(makeContrasts, c(contrast_defs, list(levels = design)))
  write.csv(
    contrast_matrix,
    file.path(out_dir, paste0(cohort_label, "_contrast_matrix.csv"))
  )

  fingerprint <- data.frame(
    gene_id = rownames(dge),
    annotations[rownames(dge), c("ENSEMBL", "SYMBOL", "GENENAME")],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  summary_rows <- list()

  for (contrast_name in colnames(contrast_matrix)) {
    test <- glmQLFTest(fit, contrast = contrast_matrix[, contrast_name])
    tab <- topTags(test, n = Inf, sort.by = "none")$table

    result <- data.frame(
      gene_id = rownames(tab),
      annotations[rownames(tab), c("ENSEMBL", "SYMBOL", "GENENAME")],
      comparison = contrast_name, tab,
      row.names = NULL, check.names = FALSE
    )
    write.csv(
      result,
      file.path(out_dir, paste0("DE_", cohort_label, "_", contrast_name, ".csv")),
      row.names = FALSE
    )

    fingerprint[[contrast_name]] <- tab[fingerprint$gene_id, "logFC"]

    summary_rows[[contrast_name]] <- data.frame(
      comparison = contrast_name,
      tested_genes = nrow(tab),
      FDR_0.05 = sum(tab$FDR < 0.05),
      FDR_0.10 = sum(tab$FDR < 0.10),
      nominal_0.05 = sum(tab$PValue < 0.05),
      up_FDR_0.05 = sum(tab$FDR < 0.05 & tab$logFC > 0),
      down_FDR_0.05 = sum(tab$FDR < 0.05 & tab$logFC < 0)
    )
  }

  write.csv(
    fingerprint,
    file.path(out_dir, paste0(cohort_label, "_prn_fingerprint.csv")),
    row.names = FALSE
  )
  de_summary <- do.call(rbind, summary_rows)
  write.csv(
    de_summary,
    file.path(out_dir, paste0(cohort_label, "_DE_summary.csv")),
    row.names = FALSE
  )

  cat("\n", cohort_label, "differential-expression summary:\n")
  print(de_summary, row.names = FALSE)

  fingerprint
}

# ------------------------------------------------------------
# 1. YOUNG cohort (August, cortex, shielded only)
#    4 groups -> all 6 pairwise contrasts, same structure as AZD.
# ------------------------------------------------------------

young_group_levels <- c("ctrl", "LD", "MD", "HD")

young_build_group <- function(meta) {
  meta$Group
}

young_contrasts <- list(
  LD_vs_ctrl = "LD - ctrl",
  MD_vs_ctrl = "MD - ctrl",
  HD_vs_ctrl = "HD - ctrl",
  MD_vs_LD   = "MD - LD",
  HD_vs_LD   = "HD - LD",
  HD_vs_MD   = "HD - MD"
)

young_fingerprint <- run_cohort(
  cohort_label   = "young",
  count_file     = young_count_file,
  meta_file      = meta_file,
  tissue_value   = "Cortex_Young",
  group_levels   = young_group_levels,
  build_group    = young_build_group,
  contrast_defs  = young_contrasts,
  out_dir        = out_dir
)

# ------------------------------------------------------------
# 2. OLD cohort (November, whole brain, Rad1-direct / Rad2-shielded)
#    Control + 4 route-x-dose groups -> 5 groups -> 8 contrasts.
#    Drops the 2 combinations that cross BOTH dose and route at
#    once (Rad1_MD vs Rad2_HD, Rad1_HD vs Rad2_MD) -- not a clean
#    single biological question.
# ------------------------------------------------------------

old_group_levels <- c("Control", "Rad1_MD", "Rad1_HD", "Rad2_MD", "Rad2_HD")

old_build_group <- function(meta) {
  gsub(" +", "_", meta$Group)
}

old_contrasts <- list(
  Rad1_MD_vs_Control = "Rad1_MD - Control",
  Rad1_HD_vs_Control = "Rad1_HD - Control",
  Rad2_MD_vs_Control = "Rad2_MD - Control",
  Rad2_HD_vs_Control = "Rad2_HD - Control",
  Rad1_HD_vs_Rad1_MD = "Rad1_HD - Rad1_MD",
  Rad2_HD_vs_Rad2_MD = "Rad2_HD - Rad2_MD",
  Rad1_vs_Rad2_at_MD = "Rad1_MD - Rad2_MD",
  Rad1_vs_Rad2_at_HD = "Rad1_HD - Rad2_HD"
)

old_fingerprint <- run_cohort(
  cohort_label   = "old",
  count_file     = old_count_file,
  meta_file      = meta_file,
  tissue_value   = "WholeBrain_Old",
  group_levels   = old_group_levels,
  build_group    = old_build_group,
  contrast_defs  = old_contrasts,
  out_dir        = out_dir
)

# ------------------------------------------------------------
# 3. Done. Two separate fingerprints, two separate output sets.
#    Do NOT merge young_fingerprint and old_fingerprint into one
#    file/matrix -- they come from different edgeR fits on
#    different tissue from different mice at different ages, and
#    PRN should be run on each as its own fingerprint.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Done. Two independent fingerprints written:\n")
cat("  ", file.path(out_dir, "young_prn_fingerprint.csv"),
    "(", ncol(young_fingerprint) - 4, "contrasts )\n")
cat("  ", file.path(out_dir, "old_prn_fingerprint.csv"),
    "(", ncol(old_fingerprint) - 4, "contrasts )\n")
cat("These are NOT combined. Run prn_build() on each separately.\n")
cat("Output directory:", out_dir, "\n")
