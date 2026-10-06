# Simulate target / source data under the four-group model

Generates a toy two-population data set for examples and tests. Variants
are independent; each variant is assigned to one of four configurations
with probabilities `prop`: `"00"` (null in both populations), `"01"`
(source-specific signal: target effect zero, source effect non-zero),
`"10"` (target-specific signal) and `"11"` (shared signal, with source
effect equal to the target effect times a random factor around 1).
Allele frequencies and covariate distributions differ between
populations; the genetic contribution is centred within each population
so that the pooled analysis is free of population stratification (in
real data, include genetic principal components in `X`). A surrogate `S`
that is informative of `Y` is also returned so that the doubly robust
estimator can be illustrated.

## Usage

``` r
heart_simulate(
  n_tar = 2000,
  n_src = 20000,
  J = 2000,
  prop = c(`00` = 0.9, `01` = 0.04, `10` = 0.01, `11` = 0.05),
  h2_snp = 0.004,
  seed = NULL
)
```

## Arguments

- n_tar, n_src:

  Target and source sample sizes.

- J:

  Number of variants.

- prop:

  Probabilities of the configurations `c("00", "01", "10", "11")`.

- h2_snp:

  Variance explained by each causal variant in its population.

- seed:

  Optional random seed.

## Value

A list with `Y`, `D` (1 = target), `X` (two covariates), `G` (n x J
genotype matrix), `S` (surrogate) and `truth`, a data.frame with the
configuration and the true target and source effects of each variant.

## Examples

``` r
sim <- heart_simulate(n_tar = 500, n_src = 2000, J = 50, seed = 1)
table(sim$truth$config)
#> 
#> 00 01 10 11 
#> 43  3  1  3 
```
