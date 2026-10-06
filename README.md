# HEART

**H**eterogeneity-aware **E**mpowerment via **A**daptive **R**obust **T**ransfer for target-population inference.

HEART tests variant–trait associations in a small target population (for example, an
under-represented ancestry in a biobank) and borrows strength from a larger source population.
It does not assume that genetic effects are shared between the two populations.
Source-specific signals do not inflate the error rate, and HEART is asymptotically no less
powerful than a target-only analysis.

## Installation

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/HEART", build_vignettes = TRUE)
```

## Method in brief

For each variant *j*:

1. **Estimation and decorrelation** (`heart_assoc()`). HEART computes the target-only
   estimator β̂<sup>tar</sup> and an auxiliary estimator β̂<sup>aux</sup>. The auxiliary
   estimator is pooled, doubly robust (DR, using surrogate-based synthetic outcomes), or
   source-only. It is then decorrelated from the target-only estimator using their joint
   sandwich covariance:

   β̂<sup>de</sup> = β̂<sup>aux</sup> − ρ̂ β̂<sup>tar</sup>, with ρ̂ = Ĉov(β̂<sup>aux</sup>, β̂<sup>tar</sup>) / V̂ar(β̂<sup>tar</sup>).

   This step gives asymptotically independent p-values (p<sub>tar</sub>, p<sub>de</sub>).
2. **Adaptive transfer via null decomposition** (`heart_test()`). The target null
   H<sub>0</sub>: β<sup>tar</sup> = 0 is split into H<sub>00</sub> ∪ H<sub>01</sub>. A
   four-group mixture model on (p<sub>tar</sub>, p<sub>de</sub>) gives the HEART statistic

   T = π̂<sub>00</sub>/(Ĉ<sub>1</sub>Ĉ<sub>2</sub>α̂<sub>1</sub>α̂<sub>2</sub>) · p<sub>tar</sub><sup>1−α̂<sub>1</sub></sup> p<sub>de</sub><sup>1−α̂<sub>2</sub></sup> + π̂<sub>01</sub>/(Ĉ<sub>1</sub>α̂<sub>1</sub>) · p<sub>tar</sub><sup>1−α̂<sub>1</sub></sup>.

   HEART estimates the composite null distribution of T from the data. This gives
   per-variant HEART p-values, FDR q-values and FWER-adjusted p-values.

## Quick start

```r
library(HEART)
sim <- heart_simulate(n_tar = 2000, n_src = 20000, J = 2000, seed = 1)

res <- heart(Y = sim$Y, G = sim$G, D = sim$D, X = sim$X, aux = "pooled", alpha = 0.05)
res                       # summary of estimated parameters and discoveries
head(res$results)         # p_tar, p_de, log_T, p_heart, q_heart, padj_fwer, reject_*
```

## Genome-wide analysis in chunks

`heart_test()` must see **all** variants, because the four-group parameters are estimated
genome-wide. Step 1 can be split across jobs:

```r
library(BEDMatrix)
G <- BEDMatrix("chr22.bed", simple_names = TRUE)        # rows must align with Y, D, X
idx <- match(sample_ids, rownames(G))

# --- one array job: a chunk of variants ---
fit_chunk <- heart_assoc(Y, G[idx, chunk_cols], D, X, aux = "pooled")
saveRDS(fit_chunk, sprintf("assoc_chunk%03d.rds", chunk_id))

# --- after all jobs finish ---
fit <- do.call(rbind, lapply(list.files(pattern = "^assoc_chunk"), readRDS))
res <- heart_test(fit$p_tar, fit$p_de_pooled, snp = fit$snp, alpha = 0.05)
```

For the Gaussian family the per-variant cost is linear in the number of covariates. The
computation uses the Frisch–Waugh–Lovell representation of the sandwich (co)variances and
is vectorised over blocks of variants.

## Doubly robust auxiliary estimator with surrogates

```r
# 1. Train a prediction model for Y from surrogates S on a held-out split,
#    then predict Y_hat on the analysis sample (e.g. ranger::ranger).
# 2. Density ratio c_i = P(D = 1 | S) / P(D = 1):
c_i <- heart_density_ratio(D, W = S)
# 3. Use aux = "dr":
fit <- heart_assoc(Y, G, D, X, aux = "dr", Y_hat = Y_hat, c_i = c_i)
res <- heart_test(fit$p_tar, fit$p_de_dr)
```

## Summary statistics only

When the target and source studies are independent, `heart_meta()` builds the decorrelated
inverse-variance meta-analytic estimator from (β̂, SE) of each study:

```r
m   <- heart_meta(beta_tar, se_tar, beta_src, se_src, snp = snp_ids)
res <- heart_test(m$p_tar, m$p_de_meta, snp = m$snp)
```

## Main functions

| Function | Purpose |
|---|---|
| `heart()` | End-to-end analysis (Step 1 + Step 2) |
| `heart_assoc()` | Target-only, auxiliary (`pooled` / `dr` / `source`) and decorrelated estimators; Gaussian or binomial |
| `heart_density_ratio()` | Density ratio for the DR estimator |
| `heart_meta()` | Decorrelated meta-analysis from summary statistics |
| `heart_test()` | HEART statistic, composite null distribution, FDR / FWER control |
| `heart_pi_estimate()`, `heart_tail_estimate()` | Four-group proportions and tail parameters |
| `heart_statistic()`, `heart_null_cdf()` | Statistic and null CDF (F<sub>00</sub>, Q̂, F̂<sub>01</sub>, N̂, F̂<sub>0</sub>) |
| `heart_simulate()` | Toy data under the four-group model |

## Citation

Yang H, Wang R, Lin X. *Heterogeneity-aware Empowerment via Adaptive Robust Transfer for
target-population inference.* (manuscript)

## License

MIT
