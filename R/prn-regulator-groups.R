#' Group regulators whose signed footprints are not separately identifiable
#'
#' Regulators are connected when the absolute cosine similarity of their
#' normalized network columns is at least `1 - threshold`. Connected components
#' form regulator groups. Activities are oriented to the first member of each
#' group and combined into the activity direction supported by the data.
#'
#' Grouping prevents individual attribution inside a nearly collinear group;
#' it does not recover hidden activity contrasts orthogonal to the identifiable
#' group direction.
#'
#' @param object A fitted `prn` object returned by [prn_build()].
#' @param threshold Maximum pair separation for joining regulators. The default
#'   `0.10` corresponds to absolute cosine similarity of at least `0.90`.
#' @return A list with `groups`, `membership`, `activities`, and `pairs` tables.
#' @export
prn_group_regulators <- function(object, threshold = 0.10) {
  if (!inherits(object, "prn")) {
    stop("`object` must inherit from class 'prn'.", call. = FALSE)
  }
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("`threshold` must be one number between 0 and 1.", call. = FALSE)
  }

  W <- as.matrix(object$weight_matrix)
  storage.mode(W) <- "double"
  regulators <- colnames(W)
  if (is.null(regulators)) regulators <- paste0("R", seq_len(ncol(W)))
  colnames(W) <- regulators
  norms <- sqrt(colSums(W^2))
  if (any(!is.finite(norms) | norms == 0)) {
    stop("Every regulator must have a finite, nonzero footprint.", call. = FALSE)
  }
  Wn <- sweep(W, 2L, norms, "/")
  cosine <- crossprod(Wn)
  separation <- 1 - abs(cosine)
  diag(separation) <- Inf

  # Connected components of the near-collinearity graph.
  adjacency <- separation <= threshold
  visited <- rep(FALSE, length(regulators))
  component <- integer(length(regulators))
  component_n <- 0L
  for (start in seq_along(regulators)) {
    if (visited[start]) next
    component_n <- component_n + 1L
    queue <- start
    visited[start] <- TRUE
    while (length(queue)) {
      current <- queue[1L]
      queue <- queue[-1L]
      component[current] <- component_n
      neighbours <- which(adjacency[current, ] & !visited)
      if (length(neighbours)) {
        visited[neighbours] <- TRUE
        queue <- c(queue, neighbours)
      }
    }
  }

  group_ids <- paste0("RG", sprintf("%03d", component))
  orientation <- integer(length(regulators))
  anchor <- character(length(regulators))
  for (group_index in seq_len(component_n)) {
    members <- which(component == group_index)
    local_anchor <- members[1L]
    local_cosine <- cosine[local_anchor, members]
    local_orientation <- ifelse(local_cosine < 0, -1L, 1L)
    orientation[members] <- local_orientation
    anchor[members] <- regulators[local_anchor]
  }

  membership <- data.frame(
    group = group_ids,
    regulator = regulators,
    anchor = anchor,
    orientation_to_anchor = orientation,
    stringsAsFactors = FALSE
  )
  membership$group_size <- ave(
    rep.int(1L, nrow(membership)), membership$group, FUN = length
  )
  membership$individual_activity_identifiable <- membership$group_size == 1L

  groups <- do.call(rbind, lapply(unique(group_ids), function(group_id) {
    dat <- membership[membership$group == group_id, , drop = FALSE]
    members <- paste(
      ifelse(dat$orientation_to_anchor > 0, "+", "-"),
      dat$regulator,
      collapse = ";"
    )
    local <- match(dat$regulator, regulators)
    within <- if (length(local) > 1L) {
      separation[local, local, drop = FALSE][upper.tri(
        separation[local, local, drop = FALSE]
      )]
    } else numeric()
    data.frame(
      group = group_id,
      anchor = dat$anchor[1L],
      group_size = nrow(dat),
      oriented_members = members,
      maximum_pair_separation = if (length(within)) max(within) else NA_real_,
      individual_activity_identifiable = nrow(dat) == 1L,
      stringsAsFactors = FALSE
    )
  }))
  rownames(groups) <- NULL

  activities <- do.call(rbind, lapply(object$comparisons, function(comparison) {
    estimate <- object$fits[[comparison]]$activity
    if (length(estimate) != length(regulators)) {
      stop("Activity vector does not align with the weight matrix.", call. = FALSE)
    }
    if (!is.null(names(estimate))) estimate <- estimate[regulators]
    do.call(rbind, lapply(unique(group_ids), function(group_id) {
      local <- which(group_ids == group_id)
      data.frame(
        comparison = comparison,
        group = group_id,
        anchor = anchor[local[1L]],
        group_size = length(local),
        identifiable_group_activity = sum(orientation[local] * estimate[local]),
        individual_activity_identifiable = length(local) == 1L,
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(activities) <- NULL

  if (length(regulators) >= 2L) {
    index <- utils::combn(seq_along(regulators), 2L)
    pairs <- data.frame(
      regulator_1 = regulators[index[1L, ]],
      regulator_2 = regulators[index[2L, ]],
      signed_cosine = cosine[index[1L, ] + nrow(cosine) * (index[2L, ] - 1L)],
      pair_separation = separation[index[1L, ] + nrow(separation) *
                                     (index[2L, ] - 1L)],
      same_group = component[index[1L, ]] == component[index[2L, ]],
      stringsAsFactors = FALSE
    )
  } else {
    pairs <- data.frame(
      regulator_1 = character(), regulator_2 = character(),
      signed_cosine = numeric(), pair_separation = numeric(),
      same_group = logical()
    )
  }

  structure(
    list(
      threshold = threshold,
      groups = groups,
      membership = membership,
      activities = activities,
      pairs = pairs
    ),
    class = "prn_regulator_groups"
  )
}
