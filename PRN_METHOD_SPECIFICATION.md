# Perturbation Response Network: method specification v0.1

## Objective

PRN is a perturbation-indexed, context-resolved network that preserves
direction, response magnitude, uncertainty, and opposing incoming regulatory
flow. It is intended to complement—not rename—differential expression, GSEA,
WGCNA, and conventional gene-regulatory-network inference.

## Input tensor

For gene `g`, perturbation `p`, context `c`, cell type `k`, and response feature
`q` (contrast, dose derivative, time derivative, threshold, saturation):

`F[g,p,c,k,q]`

The fingerprint tensor is estimated at the biological-replicate level. For
single-cell and single-nucleus data, cells are aggregated within donor before
inferential modeling.

## Network layers

### Layer A: experimental response

An assigned perturbation supports an edge `p -> g` whose value is the estimated
response fingerprint. This is the strongest direction available without a
mechanistic prior.

### Layer B: signed regulatory decomposition

Given signed prior matrix `W[g,r]`, estimate regulator activity `a[r,p,c,k]`:

`F[,p,c,k] = W a[,p,c,k] + epsilon`

Version 0.1 uses L2-normalized columns and ridge regularization. Alternative
Bayesian, sparse, and nonlinear estimators can share the same PRN object.

### Layer C: de novo direction

Gene-to-gene direction is estimated only when independent targeted
interventions, valid instruments, or temporal ordering make it identifiable.
Otherwise the edge is returned as associational or prior-informed.

## Counterdirection

For incoming contribution `e[r,g,p,c,k]`:

`positive = sum(max(e,0))`

`negative = sum(max(-e,0))`

`net = positive - negative`

`gross = positive + negative`

`counterdirection = gross - abs(net)`

`balance_fraction = counterdirection / (gross + epsilon)`

High counterdirection with small observed response defines a candidate hidden
compensation state. This does not prove the individual regulators are causal;
it identifies a falsifiable signed decomposition supported by the chosen prior
and perturbation fingerprint.

## Required validation

1. Simulated networks with known inactive, reinforcing, partially buffered,
   and exactly canceling motifs.
2. Null simulations preserving gene-level variance and regulon degree.
3. Stability across donor bootstrap, prior perturbation, and regularization.
4. Recovery of held-out responses in multi-perturbation datasets.
5. Comparison with DE/GSEA, WGCNA, signed regulon scoring, and established GRN
   methods using metrics appropriate to each method's claim.
6. External validation using Perturb-seq, CRISPRi/a, drug-dose, or time-course
   datasets with known intervention targets.

## First application

The NASA radiation cohorts provide response axes for young shielded, old
shielded, old direct-GCR, radiation-quality modification, and age/context
change. They are suitable for Layers A and B. They are not by themselves
sufficient for unrestricted de novo gene-to-gene causal direction.

