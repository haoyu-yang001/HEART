# Per-variant target-only, auxiliary and decorrelated estimators

Step 1 of HEART. For every variant (column of `G`) this function
computes the target-only estimator \\\hat\beta^{tar}\_j\\, one or more
auxiliary estimators \\\hat\beta^{aux}\_j\\ that borrow information from
the source sample, their joint sandwich covariance, and the decorrelated
estimator \$\$\hat\beta^{de}\_j = \hat\beta^{aux}\_j - \hat\rho_j
\hat\beta^{tar}\_j, \qquad \hat\rho_j =
\widehat{Cov}(\hat\beta^{aux}\_j, \hat\beta^{tar}\_j) /
\widehat{Var}(\hat\beta^{tar}\_j),\$\$ whose Wald p-value \\p\_{de,j}\\
is asymptotically independent of the target-only p-value \\p\_{tar,j}\\.
The pair \\(p\_{tar,j}, p\_{de,j})\\ is the input of
[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md).

## Usage

``` r
heart_assoc(
  Y,
  G,
  D,
  X = NULL,
  aux = "pooled",
  family = c("gaussian", "binomial"),
  Y_hat = NULL,
  c_i = NULL,
  impute = c("zero", "mean"),
  block_size = NULL,
  verbose = FALSE
)
```

## Arguments

- Y:

  Numeric phenotype vector of length n (0/1 for `family = "binomial"`).

- G:

  Genotype matrix (n x J), a numeric vector (one variant), or any object
  supporting `G[, j, drop = FALSE]` and
  [`ncol()`](https://rdrr.io/r/base/nrow.html) such as a `BEDMatrix`.
  Missing genotypes are filled according to `impute`.

- D:

  Population indicator of length n: 1 = target, 0 = source.

- X:

  Optional covariates to adjust for (matrix or data.frame, n rows); an
  intercept is always added. For multi-ancestry data include genetic
  principal components so that the pooled estimator is not confounded by
  population stratification.

- aux:

  Character vector of auxiliary estimators to compute, any of
  `"pooled"`, `"dr"`, `"source"`. The first one is used by
  [`heart()`](https://haoyu-yang001.github.io/HEART/reference/heart.md).

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

- impute:

  How to fill missing genotypes: `"zero"` (set to 0, as in the original
  HEART analysis code) or `"mean"` (variant mean within each
  population).

- block_size:

  Number of variants processed together (Gaussian family). Defaults to a
  value that keeps each block below roughly 150 MB.

- verbose:

  Print progress.

## Value

A data.frame with one row per variant and columns

- `snp`:

  Variant identifier (column names of `G` if present).

- `beta_tar`, `se_tar`, `p_tar`:

  Target-only estimate, standard error and two-sided p-value.

- `beta_<aux>`, `se_<aux>`, `p_<aux>`:

  Auxiliary estimate(s).

- `cov_<aux>`, `rho_<aux>`:

  Estimated covariance with the target-only estimator and the projection
  coefficient.

- `beta_de_<aux>`, `se_de_<aux>`, `p_de_<aux>`:

  Decorrelated estimate.

Variants that are monomorphic in the target sample get `NA`. The
attribute `"aux"` stores the auxiliary estimators that were computed.

## Details

All estimators are M-estimators of the working model \\E(Y \| G_j, X) =
\mu(\beta_0 + G_j\beta + X^\top\beta_X)\\ with the identity link
(`family = "gaussian"`) or the logit link (`family = "binomial"`):

- `"pooled"`:

  Regression on the pooled target + source sample (the auxiliary
  estimator used for the main analyses of the HEART paper).

- `"dr"`:

  Doubly robust (AIPW) estimator that combines a surrogate-based
  synthetic outcome `Y_hat` with the density ratio `c_i` (see
  [`heart_density_ratio()`](https://haoyu-yang001.github.io/HEART/reference/heart_density_ratio.md)).
  Requires `Y_hat` and `c_i`.

- `"source"`:

  Source-only regression. It is independent of the target-only
  estimator, so the decorrelated estimator equals itself.

Variances are heteroscedasticity-robust sandwich variances; covariances
between estimators use the joint sandwich formula on the shared sample.
For the Gaussian family the computations use the Frisch-Waugh-Lovell
representation and are vectorised over blocks of variants, so the cost
per variant is linear in the number of covariates.

## See also

[`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md),
[`heart()`](https://haoyu-yang001.github.io/HEART/reference/heart.md),
[`heart_density_ratio()`](https://haoyu-yang001.github.io/HEART/reference/heart_density_ratio.md)

## Examples

``` r
sim <- heart_simulate(n_tar = 1000, n_src = 5000, J = 200, seed = 1)
fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X, aux = c("pooled", "source"))
head(fit)
#>    snp    beta_tar     se_tar      p_tar   beta_pooled  se_pooled   p_pooled
#> 1 snp1  0.09441611 0.08371737 0.25940614  0.0524877443 0.02899070 0.07021782
#> 2 snp2 -0.07608608 0.06003379 0.20501675  0.0009902377 0.02630731 0.96997379
#> 3 snp3 -0.05291813 0.06635645 0.42517064  0.0187360747 0.02373323 0.42985218
#> 4 snp4 -0.02651929 0.05056922 0.59998950 -0.0208583680 0.02129047 0.32723209
#> 5 snp5  0.08404123 0.09741733 0.38830640 -0.0566990057 0.03116421 0.06885608
#> 6 snp6 -0.09647071 0.05238216 0.06552329 -0.0369561852 0.02129053 0.08259858
#>     cov_pooled rho_pooled beta_de_pooled se_de_pooled p_de_pooled beta_source
#> 1 0.0008296028  0.1183693     0.04131178   0.02724447  0.12943432  0.04576733
#> 2 0.0006580912  0.1825974     0.01488336   0.02391462  0.53370910  0.01886940
#> 3 0.0005861838  0.1331274     0.02578093   0.02202792  0.24185006  0.02902192
#> 4 0.0004491676  0.1756450    -0.01620039   0.01934916  0.40244338 -0.01866754
#> 5 0.0009659298  0.1017825    -0.06525293   0.02954476  0.02720162 -0.07524719
#> 6 0.0004418871  0.1610440    -0.02142016   0.01954798  0.27317794 -0.02512159
#>    se_source   p_source cov_source rho_source beta_de_source se_de_source
#> 1 0.03109186 0.14101990          0          0     0.04576733   0.03109186
#> 2 0.02922456 0.51849364          0          0     0.01886940   0.02922456
#> 3 0.02562483 0.25739492          0          0     0.02902192   0.02562483
#> 4 0.02344557 0.42591140          0          0    -0.01866754   0.02344557
#> 5 0.03304810 0.02279227          0          0    -0.07524719   0.03304810
#> 6 0.02340175 0.28305046          0          0    -0.02512159   0.02340175
#>   p_de_source
#> 1  0.14101990
#> 2  0.51849364
#> 3  0.25739492
#> 4  0.42591140
#> 5  0.02279227
#> 6  0.28305046
```
