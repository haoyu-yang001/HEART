# HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer

HEART tests variant-trait associations in a (small) target population
while borrowing strength from a (large) source population, without
assuming that effects are shared.

## Details

The workflow has two steps:

1.  [`heart_assoc()`](https://haoyu-yang001.github.io/HEART/reference/heart_assoc.md)
    computes, for every variant, the target-only estimator, an auxiliary
    estimator (pooled, doubly robust or source-only), and the
    decorrelated auxiliary estimator, yielding asymptotically
    independent p-values \\(p\_{tar,j}, p\_{de,j})\\. With only summary
    statistics,
    [`heart_meta()`](https://haoyu-yang001.github.io/HEART/reference/heart_meta.md)
    gives the meta-analytic counterpart.

2.  [`heart_test()`](https://haoyu-yang001.github.io/HEART/reference/heart_test.md)
    combines the p-value pairs through a four-group mixture model into
    the HEART statistic and estimates its composite null distribution,
    providing FDR and FWER control.

[`heart()`](https://haoyu-yang001.github.io/HEART/reference/heart.md)
runs both steps on individual-level data.

## See also

Useful links:

- <https://github.com/haoyu-yang001/HEART>

- <https://haoyu-yang001.github.io/HEART/>

- Report bugs at <https://github.com/haoyu-yang001/HEART/issues>

## Author

**Maintainer**: Haoyu Yang <haoyuyang@hsph.harvard.edu>

Authors:

- Ruoyu Wang

- Xihong Lin <xlin@hsph.harvard.edu>
