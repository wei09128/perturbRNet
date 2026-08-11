#' Fit arbitrary donor-aware pseudobulk contrasts with edgeR
#'
#' @param counts Gene-by-pseudobulk count matrix.
#' @param metadata Sample metadata aligned with columns of `counts`.
#' @param design_formula R formula, for example `~ 0 + condition + batch`.
#' @param contrasts Named list of numeric contrast vectors or contrast strings.
#' @param min_count Passed to edgeR filtering through `filterByExpr`; retained
#'   for API stability and currently unused.
#' @return Long data frame of gene-level contrast estimates and uncertainty.
#' @export
prn_fit_contrasts <- function(
    counts,
    metadata,
    design_formula,
    contrasts,
    min_count = NULL) {

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required.", call. = FALSE)
  }
  if (ncol(counts) != nrow(metadata)) {
    stop("Counts and metadata are not aligned.", call. = FALSE)
  }
  if (is.null(names(contrasts)) || any(names(contrasts) == "")) {
    stop("contrasts must be a named list.", call. = FALSE)
  }

  design <- stats::model.matrix(design_formula, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("Design matrix is not full rank.", call. = FALSE)
  }

  y <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  resolve_contrast <- function(x) {
    if (is.character(x) && length(x) == 1L) {
      return(limma::makeContrasts(contrasts = x, levels = design)[, 1L])
    }
    x <- as.numeric(x)
    if (length(x) != ncol(design)) {
      stop("A numeric contrast has the wrong length.", call. = FALSE)
    }
    x
  }

  pieces <- lapply(names(contrasts), function(label) {
    contrast <- resolve_contrast(contrasts[[label]])
    test <- edgeR::glmQLFTest(fit, contrast = contrast)
    tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
    data.frame(
      gene_id = rownames(tab),
      comparison = label,
      effect = tab$logFC,
      statistic = tab$F,
      p_value = tab$PValue,
      fdr = tab$FDR,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, pieces)
}

#' Convert differential results into a perturbation fingerprint matrix
#'
#' @param results Long data frame.
#' @param gene_col,comparison_col,effect_col Column names.
#' @return Gene-by-comparison numeric matrix.
#' @export
prn_make_fingerprint <- function(
    results,
    gene_col = "gene_id",
    comparison_col = "comparison",
    effect_col = "effect") {

  needed <- c(gene_col, comparison_col, effect_col)
  absent <- setdiff(needed, names(results))
  if (length(absent)) {
    stop("Missing result columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  dat <- results[, needed, drop = FALSE]
  names(dat) <- c("gene", "comparison", "effect")
  if (anyDuplicated(dat[c("gene", "comparison")])) {
    stop("Each gene-comparison combination must be unique.", call. = FALSE)
  }

  genes <- unique(as.character(dat$gene))
  comparisons <- unique(as.character(dat$comparison))
  out <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(comparisons),
    dimnames = list(genes, comparisons)
  )
  out[cbind(match(dat$gene, genes), match(dat$comparison, comparisons))] <-
    as.numeric(dat$effect)
  out
}
