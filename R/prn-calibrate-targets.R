#' Degree-conditioned calibration of target counterdirection
#'
#' Calibrate target-level counterdirection against bipartite network nulls that
#' preserve every regulator's out-degree and every target's in-degree. Edge
#' signs remain attached to their source-edge records during double-edge swaps,
#' preserving each regulator's sign counts.
#'
#' @param object A fitted `prn_fit` object returned by [prn_build()].
#' @param fingerprint Numeric gene-by-comparison response matrix used to fit
#'   `object`, with gene row names and comparison column names.
#' @param prior Signed regulatory prior with source, target, and sign columns.
#' @param source_col,target_col,sign_col Column names in `prior`.
#' @param n_null Number of successful degree-preserving null fits.
#' @param swaps_per_edge Number of attempted double-edge swaps per retained
#'   prior edge for each null draw.
#' @param min_targets Minimum observed targets per regulator, as in [prn_build()].
#' @param lambda Ridge penalty. By default, reuse the unique penalty recorded in
#'   `object`; supply a value if the fit does not expose a unique penalty.
#' @param seed Random seed.
#' @param evidence_level Evidence label assigned to null fits.
#'
#' @return A list with `calibration`, `null_summary`, `draw_diagnostics`, and
#'   `settings`. `calibration` contains observed counterdirection, incoming
#'   degree, null mean and standard deviation, a degree-conditioned Z-score,
#'   null percentile, and an empirical upper-tail probability for every
#'   target-comparison pair.
#' @export
prn_calibrate_targets <- function(
    object,
    fingerprint,
    prior,
    source_col = "source",
    target_col = "target",
    sign_col = "sign",
    n_null = 100L,
    swaps_per_edge = 5,
    min_targets = 10L,
    lambda = NULL,
    seed = 1L,
    evidence_level = "degree_preserving_null") {

  fingerprint <- as.matrix(fingerprint)
  storage.mode(fingerprint) <- "double"
  if (is.null(rownames(fingerprint)) || is.null(colnames(fingerprint))) {
    stop("`fingerprint` requires gene and comparison names.", call. = FALSE)
  }
  if (any(!is.finite(fingerprint))) {
    stop("`fingerprint` must contain only finite values.", call. = FALSE)
  }
  needed <- c(source_col, target_col, sign_col)
  absent <- setdiff(needed, names(prior))
  if (length(absent)) {
    stop("Missing prior columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  n_null <- as.integer(n_null)
  min_targets <- as.integer(min_targets)
  if (!is.finite(n_null) || n_null < 2L) {
    stop("`n_null` must be at least 2.", call. = FALSE)
  }
  if (!is.finite(swaps_per_edge) || swaps_per_edge <= 0) {
    stop("`swaps_per_edge` must be positive.", call. = FALSE)
  }

  net <- data.frame(
    source = as.character(prior[[source_col]]),
    target = as.character(prior[[target_col]]),
    sign = as.numeric(prior[[sign_col]]),
    stringsAsFactors = FALSE
  )
  net <- net[
    is.finite(net$sign) & net$sign != 0 &
      net$target %in% rownames(fingerprint),
    , drop = FALSE
  ]
  if (!nrow(net)) stop("No prior targets overlap fingerprint genes.", call. = FALSE)
  if (anyDuplicated(net[c("source", "target")])) {
    stop("Retained prior has duplicated source-target pairs.", call. = FALSE)
  }
  source_degree <- table(net$source)
  retained_sources <- names(source_degree)[source_degree >= min_targets]
  net <- net[net$source %in% retained_sources, , drop = FALSE]
  if (nrow(net) < 2L || length(unique(net$source)) < 2L ||
      length(unique(net$target)) < 2L) {
    stop("Retained prior is too small for degree-preserving swaps.", call. = FALSE)
  }

  observed <- prn_counterdirection(object)
  required_counter <- c("target", "comparison", "counterdirection")
  missing_counter <- setdiff(required_counter, names(observed))
  if (length(missing_counter)) {
    stop("Observed decomposition lacks: ",
         paste(missing_counter, collapse = ", "), call. = FALSE)
  }
  observed$target <- as.character(observed$target)
  observed$comparison <- as.character(observed$comparison)
  observed$key <- paste(observed$comparison, observed$target, sep = "\r")
  if (anyDuplicated(observed$key)) {
    stop("Observed decomposition has duplicated target-comparison pairs.",
         call. = FALSE)
  }

  if (is.null(lambda)) {
    fit_summary <- prn_summary(object)
    candidates <- unique(fit_summary$lambda[is.finite(fit_summary$lambda)])
    if (length(candidates) != 1L) {
      stop("Supply `lambda`; object does not expose one unique finite penalty.",
           call. = FALSE)
    }
    lambda <- candidates
  }
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be one finite nonnegative value.", call. = FALSE)
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  null_values <- matrix(
    NA_real_, nrow = nrow(observed), ncol = n_null,
    dimnames = list(observed$key, paste0("draw_", seq_len(n_null)))
  )
  diagnostics <- vector("list", n_null)
  attempts <- 0L
  successes <- 0L
  max_attempts <- max(n_null * 3L, n_null + 10L)

  while (successes < n_null && attempts < max_attempts) {
    attempts <- attempts + 1L
    randomized <- .prn_degree_preserving_swap(
      net,
      n_attempts = as.integer(ceiling(swaps_per_edge * nrow(net)))
    )
    fit_result <- tryCatch(
      prn_build(
        fingerprint = fingerprint,
        prior = randomized$network,
        min_targets = min_targets,
        lambda = lambda,
        evidence_level = evidence_level
      ),
      error = function(e) e
    )
    if (inherits(fit_result, "error")) {
      diagnostics[[min(attempts, n_null)]] <- data.frame(
        attempt = attempts, success = FALSE,
        accepted_swaps = randomized$accepted,
        error = conditionMessage(fit_result), stringsAsFactors = FALSE
      )
      next
    }
    decomposition <- prn_counterdirection(fit_result)
    decomposition$key <- paste(
      as.character(decomposition$comparison),
      as.character(decomposition$target), sep = "\r"
    )
    matched <- match(observed$key, decomposition$key)
    if (anyNA(matched)) next
    successes <- successes + 1L
    null_values[, successes] <- decomposition$counterdirection[matched]
    diagnostics[[successes]] <- data.frame(
      attempt = attempts, success = TRUE,
      accepted_swaps = randomized$accepted,
      acceptance_fraction = randomized$accepted / randomized$attempted,
      error = NA_character_, stringsAsFactors = FALSE
    )
  }
  if (successes < n_null) {
    stop("Only ", successes, " successful null fits after ", attempts,
         " attempts.", call. = FALSE)
  }

  null_mean <- rowMeans(null_values)
  null_sd <- apply(null_values, 1L, stats::sd)
  observed_c <- observed$counterdirection
  z <- ifelse(null_sd > 0, (observed_c - null_mean) / null_sd, NA_real_)
  exceed <- rowSums(null_values >= observed_c)
  percentile <- rowMeans(null_values <= observed_c)
  degree <- table(net$target)

  calibration <- observed
  calibration$incoming_degree <- as.integer(degree[calibration$target])
  calibration$null_mean <- null_mean
  calibration$null_sd <- null_sd
  calibration$degree_conditioned_z <- z
  calibration$null_percentile <- percentile
  calibration$empirical_p_upper <- (1 + exceed) / (1 + n_null)
  calibration$key <- NULL

  null_summary <- data.frame(
    comparison = calibration$comparison,
    target = calibration$target,
    incoming_degree = calibration$incoming_degree,
    observed_counterdirection = observed_c,
    null_mean = null_mean,
    null_sd = null_sd,
    degree_conditioned_z = z,
    null_percentile = percentile,
    empirical_p_upper = (1 + exceed) / (1 + n_null),
    stringsAsFactors = FALSE
  )

  list(
    calibration = calibration,
    null_summary = null_summary,
    draw_diagnostics = do.call(rbind, diagnostics[!vapply(diagnostics, is.null, logical(1))]),
    settings = list(
      n_null = n_null,
      swaps_per_edge = swaps_per_edge,
      min_targets = min_targets,
      lambda = lambda,
      seed = seed,
      retained_edges = nrow(net),
      regulators = length(unique(net$source)),
      targets = length(unique(net$target))
    )
  )
}

.prn_degree_preserving_swap <- function(net, n_attempts) {
  out <- net
  n <- nrow(out)
  keys <- paste(out$source, out$target, sep = "\r")
  occupied <- new.env(hash = TRUE, parent = emptyenv())
  for (key in keys) assign(key, TRUE, envir = occupied)
  accepted <- 0L

  for (attempt in seq_len(n_attempts)) {
    index <- sample.int(n, 2L, replace = FALSE)
    i <- index[1L]
    j <- index[2L]
    s1 <- out$source[i]
    s2 <- out$source[j]
    t1 <- out$target[i]
    t2 <- out$target[j]
    if (s1 == s2 || t1 == t2) next
    proposed_1 <- paste(s1, t2, sep = "\r")
    proposed_2 <- paste(s2, t1, sep = "\r")
    if (exists(proposed_1, envir = occupied, inherits = FALSE) ||
        exists(proposed_2, envir = occupied, inherits = FALSE)) next
    old_1 <- paste(s1, t1, sep = "\r")
    old_2 <- paste(s2, t2, sep = "\r")
    rm(list = c(old_1, old_2), envir = occupied)
    assign(proposed_1, TRUE, envir = occupied)
    assign(proposed_2, TRUE, envir = occupied)
    out$target[i] <- t2
    out$target[j] <- t1
    accepted <- accepted + 1L
  }

  list(network = out, attempted = n_attempts, accepted = accepted)
}
