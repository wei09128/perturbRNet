#' Aggregate cells into donor-aware pseudobulk samples
#'
#' @param counts Gene-by-cell count matrix. Dense or sparse.
#' @param metadata Data frame with one row per cell in the same order as
#'   columns of `counts`.
#' @param sample_col Biological replicate/donor column.
#' @param cell_type_col Cell-type or cell-state column.
#' @param condition_cols Columns defining perturbation and context.
#' @param min_cells Minimum cells required in one pseudobulk sample.
#' @return A list containing pseudobulk counts and sample metadata.
#' @export
prn_aggregate_pseudobulk <- function(
    counts,
    metadata,
    sample_col,
    cell_type_col,
    condition_cols,
    min_cells = 20L) {

  if (ncol(counts) != nrow(metadata)) {
    stop("ncol(counts) must equal nrow(metadata).", call. = FALSE)
  }

  needed <- unique(c(sample_col, cell_type_col, condition_cols))
  absent <- setdiff(needed, names(metadata))
  if (length(absent)) {
    stop("Missing metadata columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(rownames(counts))) {
    stop("counts must have gene row names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(counts))) {
    stop("counts contains duplicated gene row names.", call. = FALSE)
  }

  key_data <- metadata[, needed, drop = FALSE]
  if (anyNA(key_data)) {
    stop("Aggregation metadata contains missing values.", call. = FALSE)
  }

  group_key <- do.call(
    paste,
    c(lapply(key_data, as.character), sep = "\r")
  )
  group_factor <- factor(group_key, levels = unique(group_key))
  group_index <- as.integer(group_factor)
  cell_counts <- tabulate(group_index, nbins = nlevels(group_factor))
  keep_group <- cell_counts >= as.integer(min_cells)

  if (!any(keep_group)) {
    stop("No pseudobulk sample passed min_cells.", call. = FALSE)
  }

  design_cells <- Matrix::sparseMatrix(
    i = seq_len(ncol(counts)),
    j = group_index,
    x = 1,
    dims = c(ncol(counts), nlevels(group_factor))
  )
  pb_counts <- counts %*% design_cells
  pb_counts <- pb_counts[, keep_group, drop = FALSE]

  first_cell <- match(levels(group_factor), group_key)
  pb_meta <- key_data[first_cell, , drop = FALSE]
  pb_meta$n_cells <- cell_counts
  pb_meta <- pb_meta[keep_group, , drop = FALSE]
  rownames(pb_meta) <- paste0("PB", seq_len(nrow(pb_meta)))
  colnames(pb_counts) <- rownames(pb_meta)

  list(counts = pb_counts, metadata = pb_meta)
}

