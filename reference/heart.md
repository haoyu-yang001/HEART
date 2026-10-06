# HEART: end-to-end analysis from individual-level data

Convenience wrapper that runs
[`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md)
(Step 1: target-only, auxiliary and decorrelated estimators) followed by
[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md)
(Step 2: adaptive transfer via null decomposition) using the first
auxiliary estimator in `aux`.

## Usage

``` r
heart(
  Y,
  G,
  D,
  X = NULL,
  aux = "pooled",
  family = c("gaussian", "binomial"),
  Y_hat = NULL,
  c_i = NULL,
  alpha = 0.05,
  block_size = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- Y:

  Numeric phenotype vector of length n (0/1 for `family = "binomial"`).

- G:

  Genotype matrix (n x J), a numeric vector (one variant), or any object
  supporting `G[, j, drop = FALSE]` and
  [`ncol()`](https://rdrr.io/r/base/nrow.html) such as a `BEDMatrix`.
  Missing genotypes are mean-imputed within each variant and each
  population.

- D:

  Population indicator of length n: 1 = target, 0 = source.

- X:

  Optional covariates to adjust for (matrix or data.frame, n rows); an
  intercept is always added. For multi-ancestry data include genetic
  principal components so that the pooled estimator is not confounded by
  population stratification.

- aux:

  Character vector of auxiliary estimators to compute, any of
  `"pooled"`, `"dr"`, `"source"`. The first one is used by `heart()`.

- family:

  `"gaussian"` (identity link) or `"binomial"` (logit link).

- Y_hat:

  Synthetic (predicted) outcome for every individual, required for
  `aux = "dr"`. It should be produced by a model trained on data that
  are not used here (e.g. a separate split), using surrogates only.

- c_i:

  Density ratio \\c_i = P(D=1 \| S_i, Z_i) / P(D=1)\\ for every
  individual, required for `aux = "dr"`; see
  [`heart_density_ratio()`](https://haoyu-yang001.github.io/HEART/reference/heart_density_ratio.md).

- alpha:

  Nominal FDR / FWER level used for the `reject_*` columns.

- block_size:

  Number of variants processed together (Gaussian family). Defaults to a
  value that keeps each block below roughly 150 MB.

- verbose:

  Print progress.

- ...:

  Further arguments passed to
  [`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md).

## Value

An object of class `"heart"` (see
[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md))
with an additional element `assoc` holding the output of
[`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md).

## Details

For genome-wide analyses that do not fit in memory, run
[`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md)
on chunks of variants (e.g. per chromosome or per array job), row-bind
the results, and call
[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md)
once on all variants.

## Examples

``` r
sim <- heart_simulate(n_tar = 1500, n_src = 6000, J = 1000, seed = 1)
res <- heart(sim$Y, sim$G, sim$D, sim$X, aux = "pooled")
res
#> HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer
#>   tests: 1000
#>   four-group proportions: pi00 = 0.8100, pi01 = 0.0750, pi10 = 0.1150, pi11 = 0.0000
#>   tail parameters: alpha1 = 0.401, C1 = 1.014, alpha2 = 0.095, C2 = 1.020
#>   discoveries at level 0.05: FDR = 6, FWER = 2 (target-only BH = 2, Bonferroni = 0)
# power gain over the target-only analysis on the simulated truth
table(HEART = res$results$reject_fdr, truth = sim$truth$beta_tar != 0)
#>        truth
#> HEART   FALSE TRUE
#>   FALSE   949   45
#>   TRUE      0    6
```
