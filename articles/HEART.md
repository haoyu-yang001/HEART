# Introduction to HEART

HEART tests the association between each variant and a phenotype **in a
target population**. It borrows information from a larger source
population whose genetic effects may differ. This vignette walks through
the two steps of the method on simulated data.

``` r

library(HEART)
sim <- heart_simulate(n_tar = 3000, n_src = 20000, J = 3000, seed = 2026)
table(sim$truth$config)
#> 
#>   00   01   10   11 
#> 2732  110   30  128
```

The configurations are `00` (null in both populations), `01`
(source-specific signal, so the target null holds), `10`
(target-specific signal) and `11` (shared signal).

## Step 1: target-only, auxiliary and decorrelated estimators

``` r

fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X, aux = c("pooled", "source"))
head(fit[, c("snp", "beta_tar", "p_tar", "beta_pooled", "rho_pooled",
             "beta_de_pooled", "p_de_pooled")])
#>    snp    beta_tar      p_tar  beta_pooled rho_pooled beta_de_pooled
#> 1 snp1 -0.02314106 0.60962163 -0.012256548 0.09230360   -0.010120545
#> 2 snp2  0.05525188 0.18973649  0.002104585 0.12017348   -0.004535225
#> 3 snp3 -0.05259263 0.54332479 -0.008172641 0.06696513   -0.004650768
#> 4 snp4  0.04210665 0.37582215 -0.014303461 0.13383141   -0.019938654
#> 5 snp5 -0.08442441 0.02783308  0.001396242 0.14270131    0.013443716
#> 6 snp6 -0.08371033 0.26529879  0.031850898 0.13492930    0.043145875
#>   p_de_pooled
#> 1   0.4785980
#> 2   0.7565946
#> 3   0.8371745
#> 4   0.2565701
#> 5   0.3557711
#> 6   0.1170738
```

The pooled estimator is not valid for the target null. Its p-values are
small for the source-specific variants (`01`), whose target effect is
zero:

``` r

src_only <- sim$truth$config == "01"
c(pooled = mean(fit$p_pooled[src_only] < 1e-4),
  target = mean(fit$p_tar[src_only] < 1e-4))
#>    pooled    target 
#> 0.8909091 0.0000000
```

## Step 2: HEART test

``` r

res <- heart_test(fit$p_tar, fit$p_de_pooled, snp = fit$snp, alpha = 0.05)
res
#> HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer
#>   tests: 3000
#>   four-group proportions: pi00 = 0.8815, pi01 = 0.0807, pi10 = 0.0263, pi11 = 0.0115
#>   tail parameters: alpha1 = 0.129, C1 = 0.838, alpha2 = 0.056, C2 = 0.964
#>   discoveries at level 0.05: FDR = 84, FWER = 16 (target-only BH = 24, Bonferroni = 4)
```

Compare the discoveries with the truth:

``` r

signal <- sim$truth$beta_tar != 0
bh <- p.adjust(fit$p_tar, "BH") <= 0.05
rbind(
  HEART = c(discoveries = sum(res$results$reject_fdr),
            FDP = mean(!signal[res$results$reject_fdr]),
            power = mean(res$results$reject_fdr[signal])),
  `Target-only BH` = c(sum(bh), mean(!signal[bh]), mean(bh[signal]))
)
#>                discoveries        FDP     power
#> HEART                   84 0.03571429 0.5126582
#> Target-only BH          24 0.04166667 0.1455696
```

## Inspecting the null distribution

[`heart_null_cdf()`](https://haoyu-yang001.github.io/HEART/reference/heart_null_cdf.md)
returns the components of the estimated composite null distribution:
F₀₀, Q̂, F̂₀₁, N̂ and F̂₀.

``` r

nc <- res$null_cdf
head(nc)
#>        log_t          F00            Q          F01            N           F0
#> 1 -14.234983 1.768597e-09 1.020412e-08 9.374244e-08 9.127820e-09 9.486187e-09
#> 2 -11.313273 5.056040e-08 2.917139e-07 2.679895e-06 2.609447e-07 2.711897e-07
#> 3  -9.978783 2.338409e-07 1.349171e-06 1.239446e-05 1.206864e-06 1.254247e-06
#> 4  -9.743951 3.061689e-07 1.766476e-06 1.622813e-05 1.580153e-06 1.642191e-06
#> 5  -9.709737 3.184294e-07 1.837214e-06 1.687798e-05 1.643430e-06 1.707953e-06
#> 6  -9.572134 3.729018e-07 2.151499e-06 1.976523e-05 1.924565e-06 2.000125e-06
```
