# perturbRNet

`perturbRNet` constructs **Perturbation Response Networks (PRNs)** for
scRNA-seq and snRNA-seq experiments. It is designed for arbitrary perturbations,
including genetic perturbations, drugs, radiation, infection, dose series,
time courses, and interactions with age, genotype, sex, tissue, or treatment.

The central distinction is:

> A small net molecular response may represent either inactivity or strong,
> counterbalanced regulatory inputs.

For target `j` with inferred signed incoming contributions `e_ij`, PRN reports:

```text
net              N_j = sum_i e_ij
gross            G_j = sum_i |e_ij|
counterdirection C_j = G_j - |N_j|
balance          B_j = C_j / (G_j + epsilon)
```

## Statistical unit

Cells are not treated as independent biological replicates. The default
workflow aggregates cells by donor × cell type × experimental condition and
fits donor-aware pseudobulk models. This avoids pseudoreplication while still
allowing cell-type-specific PRNs.

## Evidence levels

PRN deliberately separates:

1. `experimental_direct`: perturbation-to-response edges supported by assigned
   interventions;
2. `prior_informed`: regulator-to-target edges decomposed using a signed prior;
3. `candidate_directional`: data-driven direction requiring multiple targeted
   perturbations or temporal information;
4. `associational`: unoriented response similarity, never labeled causal.

## Minimal workflow

```r
library(perturbRNet)

pb <- prn_aggregate_pseudobulk(
  counts = counts,
  metadata = cell_metadata,
  sample_col = "donor",
  cell_type_col = "cell_type",
  condition_cols = c("perturbation", "dose", "time")
)

# Fit each cell type separately or include valid interactions in the design.
de <- prn_fit_contrasts(
  counts = pb$counts,
  metadata = pb$metadata,
  design_formula = ~ 0 + perturbation + batch,
  contrasts = list(drug_vs_control = "perturbationdrug-perturbationcontrol")
)

fingerprint <- prn_make_fingerprint(de)

prn <- prn_build(
  fingerprint = fingerprint,
  prior = signed_prior,
  source_col = "source",
  target_col = "target",
  sign_col = "sign"
)

prn_summary(prn)
head(prn_counterdirection(prn))
head(prn_edge_table(prn))
```

## What version 0.1 does—and does not claim

Version 0.1 provides donor-aware perturbation fingerprints and prior-informed
signed decomposition. It does not infer a causal gene-to-gene network from one
contrast. De novo direction will be enabled only for designs with sufficient
independent interventions or time ordering.

