# Estimate the four-group proportions from a pair of p-value sequences

Estimates \\\pi\_{00}, \pi\_{01}, \pi\_{10}, \pi\_{11}\\ of the
four-group model for \\(p\_{tar,j}, p\_{de,j})\\, where the first index
refers to \\\beta^{tar}\_j\\ and the second to \\\beta^{de}\_j\\ (0 =
null, 1 = non-null). The marginal and joint null proportions are
estimated by Storey-type ratios \$\$\hat\pi\_{0\cdot}(\lambda) =
\frac{\\\\p\_{tar,j} \> \lambda\\}{J(1-\lambda)},\quad
\hat\pi\_{\cdot0}(\lambda) = \frac{\\\\p\_{de,j} \>
\lambda\\}{J(1-\lambda)},\quad \hat\pi\_{00}(\lambda) =
\frac{\\\\p\_{tar,j} \> \lambda, p\_{de,j} \>
\lambda\\}{J(1-\lambda)^2},\$\$ with \\\hat\pi\_{01} =
(\hat\pi\_{0\cdot} - \hat\pi\_{00}) \vee 0\\ and \\\hat\pi\_{10} =
(\hat\pi\_{\cdot0} - \hat\pi\_{00}) \vee 0\\. A marginal proportion is
set to one when a one-sided Kolmogorov-Smirnov test does not reject
uniformity. The tuning parameter \\\lambda\\ is chosen among `lambdas`
by the QQ-slope criterion for the maximum p-value (as in the HDMT
procedure of Dai, Stanford and LeBlanc, 2022).

## Usage

``` r
heart_pi_estimate(
  p_tar,
  p_de,
  lambdas = c(0.5, 0.6, 0.7, 0.8),
  ks_level = 0.05
)
```

## Arguments

- p_tar, p_de:

  Target-only and decorrelated p-values (same length).

- lambdas:

  Candidate values of \\\lambda\\.

- ks_level:

  Level of the Kolmogorov-Smirnov uniformity check.

## Value

A list with `pi00`, `pi01`, `pi10`, `pi11`, the marginal null
proportions `pi0_tar` (\\\hat\pi\_{0\cdot}\\) and `pi0_de`
(\\\hat\pi\_{\cdot0}\\), and the selected `lambda`.
