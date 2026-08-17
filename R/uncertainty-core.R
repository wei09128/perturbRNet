#' Fit an uncertainty-aware Perturbation Response Network
#'
#' Repeatedly samples plausible regulatory signs, fits one PRN per sampled
#' network, scores each sampled network by data support, and summarizes which
#' target-level counterdirection signals and regulator activities remain stable.
#'
#' @param fingerprint Gene-by-comparison response matrix.
#' @param prior Directed edge table.
#' @param source_col,target_col,sign_col Prior column names.
#' @param confidence_col Column containing the probability that the reported
#'   sign is correct. Values must lie in [0,1]. A value of 1 fixes the sign;
#'   0.5 represents an uninformative sign; 0 means the opposite sign is known.
#' @param n_draws Number of plausible signed networks to sample.
#' @param min_targets Minimum targets per regulator passed to `prn_build()`.
#' @param lambda Ridge penalty passed to `prn_build()`.
#' @param counter_threshold Counterdirection magnitude counted as present.
#' @param top_fraction Fraction of targets called top candidates in each draw.
#' @param seed Random seed.
#' @param keep_draws Retain draw-level target results for diagnostics.
#' @return A `prn_uncertain` object containing expected-sign and both
#'   unweighted and data-supported ensemble summaries.
#' @export
prn_build_uncertain <- function(
    fingerprint,
    prior,
    source_col = "source",
    target_col = "target",
    sign_col = "sign",
    confidence_col = "confidence",
    n_draws = 100L,
    min_targets = 10L,
    lambda = NULL,
    counter_threshold = 0,
    top_fraction = 0.10,
    seed = 1L,
    keep_draws = FALSE) {

  needed <- c(source_col, target_col, sign_col, confidence_col)
  absent <- setdiff(needed, names(prior))
  if (length(absent)) {
    stop("Missing prior columns: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  confidence <- as.numeric(prior[[confidence_col]])
  raw_weight <- as.numeric(prior[[sign_col]])
  if (any(!is.finite(confidence)) || any(confidence < 0 | confidence > 1)) {
    stop("Sign confidence must contain finite probabilities in [0,1].",
         call. = FALSE)
  }
  if (any(!is.finite(raw_weight)) || any(raw_weight == 0)) {
    stop("Prior signs/weights must be finite and nonzero.", call. = FALSE)
  }
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 2) {
    stop("n_draws must be at least 2.", call. = FALSE)
  }
  if (length(top_fraction) != 1L || !is.finite(top_fraction) ||
      top_fraction <= 0 || top_fraction > 1) {
    stop("top_fraction must lie in (0,1].", call. = FALSE)
  }

  n_draws <- as.integer(n_draws)
  reported_sign <- sign(raw_weight)
  magnitude <- abs(raw_weight)

  # Deterministic companion analysis using the expected sign:
  # E[S] = c*s + (1-c)*(-s) = s*(2c-1).
  # Thus c=0.5 contributes no directional information, while c<0.5 reverses
  # the reported sign because its opposite is more probable.
  weighted_prior <- prior
  weighted_prior[[sign_col]] <- magnitude * reported_sign * (2 * confidence - 1)
  weighted_prior <- weighted_prior[
    weighted_prior[[sign_col]] != 0,
    , drop = FALSE
  ]
  weighted_fit <- prn_build(
    fingerprint = fingerprint,
    prior = weighted_prior,
    source_col = source_col,
    target_col = target_col,
    sign_col = sign_col,
    min_targets = min_targets,
    lambda = lambda,
    evidence_level = "expected_sign_prior"
  )

  set.seed(as.integer(seed))
  target_draws <- vector("list", n_draws)
  activity_draws <- vector("list", n_draws)
  draw_diagnostics <- vector("list", n_draws)

  for (draw_id in seq_len(n_draws)) {
    keep_reported <- stats::runif(length(confidence)) <= confidence
    sampled_prior <- prior
    sampled_prior[[sign_col]] <- magnitude * reported_sign *
      ifelse(keep_reported, 1, -1)

    fit <- prn_build(
      fingerprint = fingerprint,
      prior = sampled_prior,
      source_col = source_col,
      target_col = target_col,
      sign_col = sign_col,
      min_targets = min_targets,
      lambda = lambda,
      evidence_level = "sign_uncertainty_ensemble"
    )

    target_result <- prn_counterdirection(fit)
    target_result$draw_id <- draw_id

    # Top-candidate membership is comparison-specific and remains meaningful
    # even when counterdirection scales differ between sampled networks.
    target_result$top_candidate <- FALSE
    for (comparison in unique(target_result$comparison)) {
      idx <- which(target_result$comparison == comparison)
      k <- max(1L, ceiling(length(idx) * top_fraction))
      ord <- idx[order(
        target_result$counterdirection[idx],
        decreasing = TRUE,
        na.last = NA
      )]
      if (length(ord)) {
        selected <- ord[seq_len(min(k, length(ord)))]
        target_result$top_candidate[selected] <- TRUE
      }
    }
    target_draws[[draw_id]] <- target_result

    activity_draws[[draw_id]] <- do.call(rbind, lapply(
      fit$comparisons,
      function(comparison) {
        activity <- fit$fits[[comparison]]$activity
        data.frame(
          source = names(activity),
          comparison = comparison,
          draw_id = draw_id,
          activity = as.numeric(activity),
          stringsAsFactors = FALSE
        )
      }
    ))

    diagnostic <- prn_summary(fit)
    diagnostic$gcv <- vapply(
      diagnostic$comparison,
      function(comparison) fit$fits[[comparison]]$gcv,
      numeric(1)
    )
    diagnostic$draw_id <- draw_id
    draw_diagnostics[[draw_id]] <- diagnostic
  }

  target_draws <- do.call(rbind, target_draws)
  activity_draws <- do.call(rbind, activity_draws)
  draw_diagnostics <- do.call(rbind, draw_diagnostics)

  # Convert relative GCV support into comparison-specific normalized weights.
  # This is an empirical model-support weight, not a literal posterior
  # probability over biochemical networks.
  support_key <- interaction(
    draw_diagnostics$comparison,
    drop = TRUE,
    lex.order = TRUE
  )
  support_split <- split(draw_diagnostics, support_key)
  draw_support <- do.call(rbind, lapply(support_split, function(dat) {
    log_support <- -0.5 * dat$genes * log(pmax(dat$gcv, .Machine$double.xmin))
    shifted <- log_support - max(log_support)
    weight <- exp(shifted)
    weight <- weight / sum(weight)
    data.frame(
      comparison = dat$comparison,
      draw_id = dat$draw_id,
      log_support = log_support,
      model_support_weight = weight,
      stringsAsFactors = FALSE
    )
  }))
  rownames(draw_support) <- NULL

  target_draws <- merge(
    target_draws,
    draw_support,
    by = c("comparison", "draw_id"),
    all.x = TRUE,
    sort = FALSE
  )
  activity_draws <- merge(
    activity_draws,
    draw_support,
    by = c("comparison", "draw_id"),
    all.x = TRUE,
    sort = FALSE
  )

  weighted_quantile <- function(x, w, probability) {
    ok <- is.finite(x) & is.finite(w) & w >= 0
    x <- x[ok]
    w <- w[ok]
    if (!length(x) || sum(w) <= 0) return(NA_real_)
    ord <- order(x)
    x <- x[ord]
    w <- w[ord] / sum(w)
    x[which(cumsum(w) >= probability)[1]]
  }

  target_key <- interaction(
    target_draws$target,
    target_draws$comparison,
    drop = TRUE,
    lex.order = TRUE
  )
  target_split <- split(target_draws, target_key)
  target_summary <- do.call(rbind, lapply(target_split, function(dat) {
    w <- dat$model_support_weight
    w <- w / sum(w)
    supported_mean <- sum(w * dat$counterdirection)
    supported_sd <- sqrt(sum(w * (dat$counterdirection - supported_mean) ^ 2))
    data.frame(
      target = dat$target[1],
      comparison = dat$comparison[1],
      observed_response_z = dat$observed_response_z[1],
      mean_counterdirection = mean(dat$counterdirection, na.rm = TRUE),
      sd_counterdirection = stats::sd(dat$counterdirection, na.rm = TRUE),
      q025_counterdirection = unname(stats::quantile(
        dat$counterdirection, 0.025, na.rm = TRUE
      )),
      median_counterdirection = stats::median(
        dat$counterdirection, na.rm = TRUE
      ),
      q975_counterdirection = unname(stats::quantile(
        dat$counterdirection, 0.975, na.rm = TRUE
      )),
      probability_counterdirection = mean(
        dat$counterdirection > counter_threshold,
        na.rm = TRUE
      ),
      top_candidate_stability = mean(dat$top_candidate, na.rm = TRUE),
      supported_mean_counterdirection = supported_mean,
      supported_sd_counterdirection = supported_sd,
      supported_q025_counterdirection = weighted_quantile(
        dat$counterdirection, w, 0.025
      ),
      supported_median_counterdirection = weighted_quantile(
        dat$counterdirection, w, 0.5
      ),
      supported_q975_counterdirection = weighted_quantile(
        dat$counterdirection, w, 0.975
      ),
      supported_probability_counterdirection = sum(
        w * (dat$counterdirection > counter_threshold)
      ),
      supported_top_candidate_stability = sum(w * dat$top_candidate),
      mean_balance_fraction = mean(dat$balance_fraction, na.rm = TRUE),
      mean_hidden_compensation = mean(
        dat$hidden_compensation, na.rm = TRUE
      ),
      n_draws_observed = sum(is.finite(dat$counterdirection)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(target_summary) <- NULL

  activity_key <- interaction(
    activity_draws$source,
    activity_draws$comparison,
    drop = TRUE,
    lex.order = TRUE
  )
  activity_split <- split(activity_draws, activity_key)
  activity_summary <- do.call(rbind, lapply(activity_split, function(dat) {
    w <- dat$model_support_weight
    w <- w / sum(w)
    supported_mean <- sum(w * dat$activity)
    data.frame(
      source = dat$source[1],
      comparison = dat$comparison[1],
      mean_activity = mean(dat$activity, na.rm = TRUE),
      sd_activity = stats::sd(dat$activity, na.rm = TRUE),
      q025_activity = unname(stats::quantile(
        dat$activity, 0.025, na.rm = TRUE
      )),
      q975_activity = unname(stats::quantile(
        dat$activity, 0.975, na.rm = TRUE
      )),
      probability_positive = mean(dat$activity > 0, na.rm = TRUE),
      supported_mean_activity = supported_mean,
      supported_sd_activity = sqrt(sum(
        w * (dat$activity - supported_mean) ^ 2
      )),
      supported_q025_activity = weighted_quantile(dat$activity, w, 0.025),
      supported_q975_activity = weighted_quantile(dat$activity, w, 0.975),
      supported_probability_positive = sum(w * (dat$activity > 0)),
      n_draws_observed = sum(is.finite(dat$activity)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(activity_summary) <- NULL

  object <- list(
    weighted_fit = weighted_fit,
    target_summary = target_summary,
    activity_summary = activity_summary,
    draw_diagnostics = draw_diagnostics,
    draw_support = draw_support,
    target_draws = if (isTRUE(keep_draws)) target_draws else NULL,
    activity_draws = if (isTRUE(keep_draws)) activity_draws else NULL,
    n_draws = n_draws,
    counter_threshold = counter_threshold,
    top_fraction = top_fraction,
    call = match.call()
  )
  class(object) <- "prn_uncertain"
  object
}

#' @export
print.prn_uncertain <- function(x, ...) {
  cat("Uncertainty-aware Perturbation Response Network\n")
  cat("  sampled networks: ", x$n_draws, "\n", sep = "")
  cat("  targets:          ", length(unique(x$target_summary$target)), "\n", sep = "")
  cat("  comparisons:      ",
      length(unique(x$target_summary$comparison)), "\n", sep = "")
  cat("  regulators:       ",
      length(unique(x$activity_summary$source)), "\n", sep = "")
  invisible(x)
}
