# Density ratio between the target and the pooled population

Estimates \\c_i = P(D_i = 1 \| W_i) / \hat\pi\\ with \\\hat\pi =
n^{-1}\sum_i D_i\\ by a logistic regression of the population indicator
`D` on the features `W` (typically surrogates `S` and, optionally,
covariates). Because \\P(D=1\|W)/\pi = f\_{target}(W)/f\_{pooled}(W)\\,
\\c_i\\ is the weight that aligns the pooled sample with the target
population and is used by the doubly robust auxiliary estimator
(`aux = "dr"` in
[`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md)).

## Usage

``` r
heart_density_ratio(
  D,
  W,
  quadratic = FALSE,
  ridge = 1e-08,
  maxit = 50L,
  tol = 1e-08
)
```

## Arguments

- D:

  0/1 population indicator (1 = target).

- W:

  Features (matrix or data.frame); factors are expanded by
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

- quadratic:

  If `TRUE`, squared terms of the numeric columns are added.

- ridge:

  Initial ridge penalty.

- maxit, tol:

  Maximum IRLS iterations and convergence tolerance.

## Value

Numeric vector `c_i` with attributes `coef` (fitted logistic
coefficients), `e_hat` (fitted P(D = 1 \| W)) and `pi_hat`.

## Details

The logistic model is fitted by iteratively reweighted least squares
with a small ridge penalty that is increased automatically if the
weighted Gram matrix is numerically singular (e.g. under
quasi-separation).

## Examples

``` r
sim <- heart_simulate(n_tar = 500, n_src = 2000, J = 10, seed = 1)
c_i <- heart_density_ratio(sim$D, cbind(S = sim$S, sim$X))
tapply(c_i, sim$D, mean)
#>         0         1 
#> 0.9774112 1.0903551 
```
