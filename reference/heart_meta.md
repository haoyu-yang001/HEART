# Decorrelated inverse-variance meta-analysis from summary statistics

When only summary statistics from independent target and source studies
are available, the inverse-variance meta-analytic estimator
\\\hat\beta^{meta}\_j = w^{tar}\_j\hat\beta^{tar}\_j +
w^{src}\_j\hat\beta^{src}\_j\\ can serve as the auxiliary estimator.
Because the two studies are independent, \\\hat\rho_j = w^{tar}\_j\\ and
the decorrelated estimator is \\\hat\beta^{de}\_j =
w^{src}\_j\hat\beta^{src}\_j\\ with variance \\(w^{src}\_j)^2
v^{src}\_j\\; its p-value therefore coincides with the source-study
p-value. The output can be passed to
[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md)
with `p_de = p_de_meta`.

## Usage

``` r
heart_meta(beta_tar, se_tar, beta_src, se_src, snp = NULL)
```

## Arguments

- beta_tar, se_tar:

  Target-study estimates and standard errors.

- beta_src, se_src:

  Source-study estimates and standard errors (aligned to the same effect
  allele).

- snp:

  Optional variant identifiers.

## Value

A data.frame with the target-only, meta-analytic and decorrelated
estimates, standard errors and two-sided p-values.

## Examples

``` r
heart_meta(beta_tar = c(0.1, 0), se_tar = c(0.05, 0.05),
           beta_src = c(0.08, 0.05), se_src = c(0.01, 0.01))
#>   snp beta_tar se_tar      p_tar  beta_meta     se_meta       p_meta   rho_meta
#> 1   1      0.1   0.05 0.04550026 0.08076923 0.009805807 1.767630e-16 0.03846154
#> 2   2      0.0   0.05 1.00000000 0.04807692 0.009805807 9.443044e-07 0.03846154
#>   beta_de_meta  se_de_meta    p_de_meta
#> 1   0.07692308 0.009615385 1.244192e-15
#> 2   0.04807692 0.009615385 5.733031e-07
```
