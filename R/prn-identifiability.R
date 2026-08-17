#' Diagnose geometric identifiability of a perturbation response network
#'
#' Counterdirection can be interpreted only when the signed network contains
#' enough distinct target patterns to separate regulator activities. This
#' function measures that network geometry. It is a diagnostic, not a p-value
#' or a guarantee of causal identification.
#'
#' @param object A fitted perturbRNet object returned by [prn_build()].
#' @param tolerance Relative singular-value tolerance used to determine rank.
#'
#' @return A list with `global`, `regulator`, `pair`, and `target` tables.
#' @export
prn_identifiability <- function(object, tolerance = sqrt(.Machine$double.eps)) {
  W <- object$weight_matrix
  if (is.null(W)) stop("`object$weight_matrix` is missing.", call. = FALSE)
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  if (!nrow(W) || !ncol(W)) stop("The weight matrix is empty.", call. = FALSE)

  regulators <- colnames(W)
  if (is.null(regulators)) regulators <- paste0("R", seq_len(ncol(W)))
  colnames(W) <- regulators

  norms <- sqrt(colSums(W^2))
  active <- is.finite(norms) & norms > 0
  Wn <- W[, active, drop = FALSE]
  Wn <- sweep(Wn, 2L, norms[active], "/")

  singular <- if (ncol(Wn)) svd(Wn, nu = 0L, nv = 0L)$d else numeric()
  cutoff <- if (length(singular)) max(singular) * tolerance else 0
  positive_singular <- singular[singular > cutoff]
  rank <- length(positive_singular)
  condition_number <- if (rank < ncol(Wn) || !rank) {
    Inf
  } else {
    max(positive_singular) / min(positive_singular)
  }
  energy <- singular^2
  probability <- if (sum(energy) > 0) energy / sum(energy) else numeric()
  effective_rank <- if (length(probability)) {
    exp(-sum(probability[probability > 0] * log(probability[probability > 0])))
  } else 0

  global <- data.frame(
    genes = nrow(W),
    regulators = ncol(W),
    nonzero_regulators = sum(active),
    rank = rank,
    rank_fraction = if (sum(active)) rank / sum(active) else NA_real_,
    effective_rank = effective_rank,
    effective_rank_fraction = if (sum(active)) effective_rank / sum(active) else NA_real_,
    condition_number = condition_number
  )

  uniqueness <- rep(NA_real_, ncol(W))
  for (j in which(active)) {
    others <- setdiff(which(active), j)
    if (!length(others)) {
      uniqueness[j] <- 1
    } else {
      fit <- stats::lm.fit(x = W[, others, drop = FALSE], y = W[, j])
      residual <- fit$residuals
      uniqueness[j] <- max(0, min(1, sum(residual^2) / sum(W[, j]^2)))
    }
  }
  regulator <- data.frame(
    regulator = regulators,
    target_count = colSums(W != 0),
    column_norm = norms,
    unique_information = uniqueness,
    geometry = cut(
      uniqueness,
      breaks = c(-Inf, 0.10, 0.40, Inf),
      labels = c("low", "moderate", "high")
    ),
    stringsAsFactors = FALSE
  )

  if (sum(active) >= 2L) {
    index <- utils::combn(which(active), 2L)
    cosine <- colSums(Wn[, match(index[1L, ], which(active)), drop = FALSE] *
                        Wn[, match(index[2L, ], which(active)), drop = FALSE])
    pair <- data.frame(
      regulator_1 = regulators[index[1L, ]],
      regulator_2 = regulators[index[2L, ]],
      signed_cosine = cosine,
      absolute_cosine = abs(cosine),
      pair_separation = 1 - abs(cosine),
      stringsAsFactors = FALSE
    )
  } else {
    pair <- data.frame(
      regulator_1 = character(), regulator_2 = character(),
      signed_cosine = numeric(), absolute_cosine = numeric(),
      pair_separation = numeric()
    )
  }

  edge <- prn_edge_table(object)
  if (!"regulator" %in% names(edge) && "source" %in% names(edge)) {
    edge$regulator <- edge$source
  }
  required <- c("comparison", "target", "regulator", "inferred_contribution")
  if (!all(required %in% names(edge))) {
    stop("`prn_edge_table(object)` lacks required columns: ",
         paste(setdiff(required, names(edge)), collapse = ", "), call. = FALSE)
  }
  pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
  separation_lookup <- stats::setNames(
    pair$pair_separation,
    pair_key(pair$regulator_1, pair$regulator_2)
  )

  split_edge <- split(edge, interaction(edge$comparison, edge$target, drop = TRUE))
  target <- do.call(rbind, lapply(split_edge, function(dat) {
    pos <- dat[is.finite(dat$inferred_contribution) & dat$inferred_contribution > 0, ]
    neg <- dat[is.finite(dat$inferred_contribution) & dat$inferred_contribution < 0, ]
    if (!nrow(pos) || !nrow(neg)) {
      return(data.frame(
        comparison = dat$comparison[1L], target = dat$target[1L],
        opposing_pairs = 0L, weighted_pair_separation = NA_real_,
        minimum_pair_separation = NA_real_, identifiability = "not_applicable"
      ))
    }
    combinations <- expand.grid(p = seq_len(nrow(pos)), n = seq_len(nrow(neg)))
    sep <- unname(separation_lookup[pair_key(
      pos$regulator[combinations$p], neg$regulator[combinations$n]
    )])
    strength <- 2 * pmin(
      pos$inferred_contribution[combinations$p],
      abs(neg$inferred_contribution[combinations$n])
    )
    keep <- is.finite(sep) & is.finite(strength) & strength > 0
    weighted <- if (any(keep)) stats::weighted.mean(sep[keep], strength[keep]) else NA_real_
    minimum <- if (any(keep)) min(sep[keep]) else NA_real_
    label <- if (!is.finite(weighted)) "unknown" else if (weighted < 0.10) {
      "low"
    } else if (weighted < 0.40) {
      "moderate"
    } else "high"
    data.frame(
      comparison = dat$comparison[1L], target = dat$target[1L],
      opposing_pairs = sum(keep), weighted_pair_separation = weighted,
      minimum_pair_separation = minimum, identifiability = label
    )
  }))
  rownames(target) <- NULL

  structure(
    list(global = global, regulator = regulator, pair = pair, target = target),
    class = "prn_identifiability"
  )
}
