# HEART test statistic

Computes the simplified HEART statistic (with \\\pi\_{10}\\ set to zero)
\$\$T_j = \hat\pi\_{00}\hat C_1^{-1}\hat
C_2^{-1}\hat\alpha_1^{-1}\hat\alpha_2^{-1} p\_{tar,j}^{1-\hat\alpha_1}
p\_{de,j}^{1-\hat\alpha_2} + \hat\pi\_{01}\hat C_1^{-1}\hat\alpha_1^{-1}
p\_{tar,j}^{1-\hat\alpha_1}.\$\$ Small values are evidence against
\\H\_{0,j}: \beta^{tar}\_j = 0\\.

## Usage

``` r
heart_statistic(p_tar, p_de, par, log = FALSE)
```

## Arguments

- p_tar, p_de:

  Target-only and decorrelated p-values.

- par:

  List with `pi00`, `pi01`, `alpha1`, `alpha2`, `C1`, `C2` (as returned
  in `$par` by
  [`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md)).

- log:

  Return \\\log T_j\\ (recommended for very small p-values).

## Value

Numeric vector of statistics.
