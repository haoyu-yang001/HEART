# Composite null distribution of the HEART statistic

Evaluates, at the statistic values `t`, the global-null CDF \\F\_{00}(t)
= \int_0^1 G_t(q)dq\\, its empirical counterpart \\\hat Q(t) =
J^{-1}\sum_j G_t(p\_{de,j})\\, the reconstructed \\\hat F\_{01}(t) =
\\\hat Q(t) - \hat\pi\_{00}F\_{00}(t)\\/(\hat\pi\_{01}+\hat\pi\_{11})\\
(truncated to \[0, 1\]; set to 0 if \\\hat\pi\_{01}+\hat\pi\_{11} =
0\\), the null numerator \\\hat N(t) = \hat\pi\_{00}F\_{00}(t) +
\hat\pi\_{01}\hat F\_{01}(t)\\ and the composite null CDF \\\hat F_0(t)
= \hat N(t) / (\hat\pi\_{00} + \hat\pi\_{01})\\, where \\G_t(q) =
\min\\1, (t/(\kappa_1 q^{1-\alpha_2} + \kappa_2))^{1/(1-\alpha_1)}\\\\.

## Usage

``` r
heart_null_cdf(log_t, p_de, par, exact = FALSE, n_exact = 5000L, n_grid = 500L)
```

## Arguments

- log_t:

  Values of \\\log t\\.

- p_de:

  All decorrelated p-values (used for \\\hat Q\\).

- par:

  Parameter list as in
  [`heart_statistic()`](https://haoyu-yang001.github.io/HEART/reference/heart_statistic.md),
  also containing `pi11` (or `pi10`, from which \\\pi\_{11} = 1 -
  \pi\_{00} - \pi\_{01} - \pi\_{10}\\).

- exact:

  If `FALSE`, \\F\_{00}\\ is computed exactly at the `n_exact` smallest
  values and interpolated (log-log) on an `n_grid` grid elsewhere.

- n_exact, n_grid:

  See `exact`.

## Value

A data.frame with columns `log_t`, `F00`, `Q`, `F01`, `N`, `F0`, sorted
by `log_t`.
