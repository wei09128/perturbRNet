.prn_zscore <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

#' Construct a signed Perturbation Response Network
#'
#' @param fingerprint Gene-by-comparison response matrix.
#' @param prior Directed signed edge table.
#' @param source_col,target_col,sign_col Prior column names.
#' @param min_targets Minimum observed targets per source.
#' @param lambda Ridge penalty. `NULL` selects it by generalized
#'   cross-validation independently for each comparison.
#' @param evidence_level Label attached to returned regulatory edges.
#' @return A `prn` object.
#' @export
prn_build <- function(
    fingerprint,
    prior,
    source_col = "source",
    target_col = "target",
    sign_col = "sign",
    min_targets = 10L,
    lambda = NULL,
    evidence_level = "prior_informed") {

  fingerprint <- as.matrix(fingerprint)
  storage.mode(fingerprint) <- "double"
  if (is.null(rownames(fingerprint)) || is.null(colnames(fingerprint))) {
    stop("fingerprint requires gene and comparison names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(fingerprint))) {
    stop("fingerprint has duplicated genes.", call. = FALSE)
  }

  needed <- c(source_col, target_col, sign_col)
  absent <- setdiff(needed, names(prior))
  if (length(absent)) {
    stop("Missing prior columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  net <- data.frame(
    source = as.character(prior[[source_col]]),
    target = as.character(prior[[target_col]]),
    weight = as.numeric(prior[[sign_col]]),
    stringsAsFactors = FALSE
  )
  net <- net[
    is.finite(net$weight) & net$weight != 0 &
      net$target %in% rownames(fingerprint),
    , drop = FALSE
  ]

  # Collapse duplicate evidence; conflicting signs attenuate or cancel.
  key <- paste(net$source, net$target, sep = "\r")
  split_weight <- split(net$weight, key)
  collapsed <- data.frame(
    key = names(split_weight),
    weight = vapply(split_weight, mean, numeric(1)),
    stringsAsFactors = FALSE
  )
  key_parts <- strsplit(collapsed$key, "\r", fixed = TRUE)
  collapsed$source <- vapply(key_parts, `[`, character(1), 1L)
  collapsed$target <- vapply(key_parts, `[`, character(1), 2L)
  net <- collapsed[collapsed$weight != 0, c("source", "target", "weight")]

  target_count <- table(net$source)
  keep_source <- names(target_count[target_count >= as.integer(min_targets)])
  net <- net[net$source %in% keep_source, , drop = FALSE]
  if (!nrow(net)) stop("No prior edges passed filtering.", call. = FALSE)

  genes <- sort(unique(net$target))
  sources <- sort(unique(net$source))
  W <- Matrix::sparseMatrix(
    i = match(net$target, genes),
    j = match(net$source, sources),
    x = net$weight,
    dims = c(length(genes), length(sources)),
    dimnames = list(genes, sources)
  )
  norms <- sqrt(Matrix::colSums(W ^ 2))
  Wn <- W %*% Matrix::Diagonal(x = 1 / norms)
  # Matrix multiplication by an unnamed Diagonal can drop regulator names.
  # Preserve them because downstream activity and edge tables align by name.
  dimnames(Wn) <- dimnames(W)
  lambda_grid <- 10 ^ seq(-4, 4, length.out = 81)

  fits <- lapply(colnames(fingerprint), function(comparison) {
    y <- .prn_zscore(fingerprint[genes, comparison])
    ok <- is.finite(y)
    X <- Wn[ok, , drop = FALSE]
    yy <- y[ok]
    n_obs <- length(yy)

    # Fit y = intercept + X activity without densifying sparse X.
    # The intercept is unpenalized. Centered cross-products are equivalent
    # to explicitly centering every network column, but preserve sparsity.
    x_mean <- as.numeric(Matrix::colMeans(X))
    y_mean <- mean(yy)
    raw_XtX <- as.matrix(Matrix::crossprod(X))
    centered_XtX <- raw_XtX - n_obs * tcrossprod(x_mean)
    centered_Xty <- as.numeric(Matrix::crossprod(X, yy)) -
      n_obs * x_mean * y_mean
    local_eigen <- eigen(
      centered_XtX,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    local_eigen[local_eigen < 0 & local_eigen > -1e-8] <- 0

    candidates <- if (is.null(lambda)) lambda_grid else as.numeric(lambda)
    fitted_candidates <- lapply(candidates, function(lam) {
      activity <- solve(
        centered_XtX + diag(lam, ncol(X)),
        centered_Xty
      )
      intercept <- y_mean - sum(x_mean * activity)
      fitted <- intercept + as.numeric(X %*% activity)
      rss <- sum((yy - fitted) ^ 2)
      # One additional effective degree of freedom for the intercept.
      edf <- 1 + sum(local_eigen / (local_eigen + lam))
      gcv <- (rss / n_obs) / (1 - edf / n_obs) ^ 2
      list(lambda = lam, intercept = intercept,
           activity = activity, fitted = fitted,
           rss = rss, effective_df = edf, gcv = gcv, ok = ok, y = y)
    })
    fitted_candidates[[which.min(vapply(
      fitted_candidates, `[[`, numeric(1), "gcv"
    ))]]
  })
  names(fits) <- colnames(fingerprint)

  object <- list(
    fingerprint = fingerprint,
    prior = net,
    weight_matrix = Wn,
    column_norm = norms,
    comparisons = colnames(fingerprint),
    fits = fits,
    evidence_level = evidence_level,
    call = match.call()
  )
  class(object) <- "prn"
  object
}

#' Calculate counterdirectional regulatory input
#'
#' @param object A `prn` object.
#' @return Target-by-comparison decomposition table.
#' @export
prn_counterdirection <- function(object) {
  stopifnot(inherits(object, "prn"))
  W <- object$weight_matrix

  pieces <- lapply(object$comparisons, function(comparison) {
    fit <- object$fits[[comparison]]
    edge_index <- Matrix::summary(W)
    edge_contribution <- edge_index$x * fit$activity[edge_index$j]
    positive <- as.numeric(Matrix::sparseMatrix(
      i = edge_index$i,
      j = rep.int(1L, nrow(edge_index)),
      x = pmax(edge_contribution, 0),
      dims = c(nrow(W), 1L)
    ))
    negative <- as.numeric(Matrix::sparseMatrix(
      i = edge_index$i,
      j = rep.int(1L, nrow(edge_index)),
      x = pmax(-edge_contribution, 0),
      dims = c(nrow(W), 1L)
    ))
    net <- positive - negative
    gross <- positive + negative
    counter <- gross - abs(net)
    observed <- .prn_zscore(object$fingerprint[rownames(W), comparison])

    data.frame(
      target = rownames(W),
      comparison = comparison,
      observed_response_z = observed,
      positive_input = positive,
      negative_input = negative,
      net_input = net,
      gross_input = gross,
      counterdirection = counter,
      balance_fraction = ifelse(gross > 0, counter / gross, 0),
      hidden_compensation = counter / (abs(observed) + 0.25),
      residual = observed - net,
      evidence_level = object$evidence_level,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

#' Extract comparison-specific PRN edges
#'
#' @param object A `prn` object.
#' @return Long edge table with inferred contributions.
#' @export
prn_edge_table <- function(object) {
  stopifnot(inherits(object, "prn"))
  net <- object$prior
  pieces <- lapply(object$comparisons, function(comparison) {
    activity <- object$fits[[comparison]]$activity
    names(activity) <- colnames(object$weight_matrix)
    data.frame(
      source = net$source,
      target = net$target,
      comparison = comparison,
      prior_sign = sign(net$weight),
      inferred_contribution =
        (net$weight / object$column_norm[net$source]) * activity[net$source],
      evidence_level = object$evidence_level,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

#' Summarize a PRN object
#'
#' @param object A `prn` object.
#' @return Comparison-level diagnostics.
#' @export
prn_summary <- function(object) {
  stopifnot(inherits(object, "prn"))
  pieces <- lapply(object$comparisons, function(comparison) {
    fit <- object$fits[[comparison]]
    denominator <- sum((fit$y[fit$ok] - mean(fit$y[fit$ok])) ^ 2)
    data.frame(
      comparison = comparison,
      genes = sum(fit$ok),
      regulators = ncol(object$weight_matrix),
      lambda = fit$lambda,
      effective_df = fit$effective_df,
      variance_explained = 1 - fit$rss / denominator,
      evidence_level = object$evidence_level,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

#' @export
print.prn <- function(x, ...) {
  cat("Perturbation Response Network\n")
  cat("  genes:       ", nrow(x$weight_matrix), "\n", sep = "")
  cat("  regulators:  ", ncol(x$weight_matrix), "\n", sep = "")
  cat("  comparisons: ", length(x$comparisons), "\n", sep = "")
  cat("  evidence:    ", x$evidence_level, "\n", sep = "")
  invisible(x)
}
