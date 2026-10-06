# Estimate the left-tail parameters of the alternative p-value density

Under the two-group model \\f(p) = \pi_0 + (1-\pi_0) f_1(p)\\ with the
local approximation \\f_1(p) \approx C\alpha p^{\alpha-1}\\ near zero,
\\(\alpha, C)\\ solve the conditional moment equations \$\$J^{-1}\sum_j
\log(p_j) I(p_j \le r) - \pi_0(r\log r - r) - C(1-\pi_0) r^{\alpha}(\log
r - \alpha^{-1}) = 0,\$\$ \$\$J^{-1}\sum_j I(p_j \le r) = \pi_0 r +
(1-\pi_0) C r^{\alpha},\$\$ which have a closed-form solution. The
threshold \\r\\ is the \\m\\-th smallest p-value with \\m = \lceil
k\sqrt{J}\rceil\\. The solutions are truncated to `alpha_bounds` and
`C_bounds`.

## Usage

``` r
heart_tail_estimate(
  p,
  pi0,
  k = 50,
  pi0_cap = 0.996,
  alpha_bounds = c(1e-05, 1 - 1e-05),
  C_bounds = c(0.1, 3)
)
```

## Arguments

- p:

  Vector of p-values.

- pi0:

  Estimated null proportion of `p`.

- k:

  Multiplier defining the tail size \\m = k\sqrt{J}\\. `k = 1`
  corresponds to the \\J^{-1/2}\\ quantile.

- pi0_cap:

  Upper cap applied to `pi0` (avoids division by zero when the p-values
  look uniform).

- alpha_bounds, C_bounds:

  Truncation ranges for \\\alpha\\ and \\C\\.

## Value

Named numeric vector `c(alpha = , C = , r = )`.
