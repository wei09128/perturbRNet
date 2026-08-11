# Minimal simulation demonstrating the four PRN response states.
# Run after installing perturbRNet:
#   Rscript system.file("examples/simulate_prn.R", package = "perturbRNet")

library(perturbRNet)

set.seed(11)

# Four latent regulators and 400 targets. The prior is known here so the
# simulation tests decomposition separately from prior-network recovery.
regulators <- paste0("R", 1:4)
genes <- paste0("G", seq_len(400))

prior <- do.call(rbind, lapply(seq_along(regulators), function(i) {
  idx <- ((i - 1) * 100 + 1):(i * 100)
  data.frame(
    source = regulators[i],
    target = genes[idx],
    sign = sample(c(-1, 1), 100, replace = TRUE)
  )
}))

# Add 50 shared targets for an explicitly counterbalanced R1/R2 motif.
shared <- genes[1:50]
prior <- rbind(
  prior,
  data.frame(source = "R2", target = shared, sign = -prior$sign[1:50])
)

W <- matrix(0, nrow = length(genes), ncol = length(regulators),
            dimnames = list(genes, regulators))
for (i in seq_len(nrow(prior))) {
  W[prior$target[i], prior$source[i]] <- prior$sign[i]
}

activities <- cbind(
  inactive = c(0, 0, 0, 0),
  reinforcing = c(2, -2, 0, 0),
  counterbalanced = c(2, 2, 0, 0),
  context_switch = c(-2, 2, 1, -1)
)
rownames(activities) <- regulators

fingerprint <- W %*% activities +
  matrix(rnorm(length(genes) * ncol(activities), sd = 0.15),
         nrow = length(genes), dimnames = list(genes, colnames(activities)))

fit <- prn_build(
  fingerprint = fingerprint,
  prior = prior,
  min_targets = 10,
  lambda = 0.1,
  evidence_level = "simulation_known_prior"
)

print(fit)
print(prn_summary(fit))

decomposition <- prn_counterdirection(fit)
aggregate(
  cbind(counterdirection, balance_fraction, hidden_compensation) ~ comparison,
  data = decomposition,
  FUN = mean
)

