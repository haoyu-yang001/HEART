# HEART test: adaptive transfer via null decomposition

Step 2 of HEART. Given target-only p-values \\p\_{tar,j}\\ and
decorrelated p-values \\p\_{de,j}\\ for \\J\\ variants (typically
genome-wide), this function

1.  estimates the four-group proportions
    ([`heart_pi_estimate()`](https://haoyu-yang001.github.io/HEART/reference/heart_pi_estimate.md)),

2.  estimates the tail parameters \\(\alpha_1, C_1)\\ from \\p\_{tar}\\
    and \\(\alpha_2, C_2)\\ from \\p\_{de}\\
    ([`heart_tail_estimate()`](https://haoyu-yang001.github.io/HEART/reference/heart_tail_estimate.md)),

3.  computes the HEART statistic \\T_j\\
    ([`heart_statistic()`](https://haoyu-yang001.github.io/HEART/reference/heart_statistic.md)),

4.  estimates its composite null distribution
    ([`heart_null_cdf()`](https://haoyu-yang001.github.io/HEART/reference/heart_null_cdf.md)),
    and

5.  returns HEART p-values, FDR q-values and FWER-adjusted p-values.

## Usage

``` r
heart_test(
  p_tar,
  p_de,
  alpha = 0.05,
  snp = NULL,
  par = NULL,
  tail_k = 50,
  pi0_cap = 0.996,
  alpha_bounds = c(1e-05, 1 - 1e-05),
  C_bounds = c(0.1, 3),
  lambdas = c(0.5, 0.6, 0.7, 0.8),
  p_min = 1e-300,
  exact = FALSE,
  n_exact = 5000L,
  n_grid = 500L
)
```

## Arguments

- p_tar, p_de:

  Target-only and decorrelated p-values (same length).

- alpha:

  Nominal FDR / FWER level used for the `reject_*` columns.

- snp:

  Optional variant identifiers.

- par:

  Optional list of model parameters (`pi00`, `pi01`, `pi10`, `alpha1`,
  `alpha2`, `C1`, `C2`, and optionally `pi11`) to bypass estimation.

- tail_k, pi0_cap, alpha_bounds, C_bounds:

  Passed to
  [`heart_tail_estimate()`](https://haoyu-yang001.github.io/HEART/reference/heart_tail_estimate.md).

- lambdas:

  Passed to
  [`heart_pi_estimate()`](https://haoyu-yang001.github.io/HEART/reference/heart_pi_estimate.md).

- p_min:

  P-values are bounded below by `p_min` to avoid `log(0)`.

- exact, n_exact, n_grid:

  Passed to
  [`heart_null_cdf()`](https://haoyu-yang001.github.io/HEART/reference/heart_null_cdf.md).

## Value

An object of class `"heart"`: a list with

- `results`:

  data.frame with `snp`, `p_tar`, `p_de`, `log_T`, `p_heart`, `q_heart`,
  `padj_fwer`, `reject_fdr`, `reject_fwer`.

- `par`:

  Estimated parameters.

- `alpha`, `n_tests`:

  Level and number of valid tests.

## Details

The FDR procedure rejects \\\\j: T_j \le \hat t\_\alpha\\\\ with \\\hat
t\_\alpha = \sup\\t: \hat N(t)/\bar F(t) \le \alpha\\\\; equivalently
`q_heart <= alpha`. The FWER procedure rejects \\\\j: \hat N(T_j) \le
\alpha/J\\\\; equivalently `padj_fwer <= alpha`. The per-variant HEART
p-value is \\\hat F_0(T_j)\\, the estimated composite-null CDF evaluated
at the observed statistic, which can be used e.g. for Manhattan plots
with a genome-wide threshold.

The parameters must be estimated from all variants jointly; when
variants are analysed in chunks with
[`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md),
combine the chunks before calling `heart_test()`.

## Examples

``` r
sim <- heart_simulate(n_tar = 1500, n_src = 6000, J = 1000, seed = 1)
fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X)
res <- heart_test(fit$p_tar, fit$p_de_pooled, snp = fit$snp)
res
#> HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer
#>   tests: 1000
#>   four-group proportions: pi00 = 0.8100, pi01 = 0.0750, pi10 = 0.1150, pi11 = 0.0000
#>   tail parameters: alpha1 = 0.401, C1 = 1.014, alpha2 = 0.095, C2 = 1.020
#>   discoveries at level 0.05: FDR = 6, FWER = 2 (target-only BH = 2, Bonferroni = 0)
```
