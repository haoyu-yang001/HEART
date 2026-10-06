#' HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer
#'
#' HEART tests variant-trait associations in a (small) target population while
#' borrowing strength from a (large) source population, without assuming that
#' effects are shared.
#'
#' The workflow has two steps:
#' 1. [heart_assoc()] computes, for every variant, the target-only estimator,
#'    an auxiliary estimator (pooled, doubly robust or source-only), and the
#'    decorrelated auxiliary estimator, yielding asymptotically independent
#'    p-values \eqn{(p_{tar,j}, p_{de,j})}. With only summary statistics,
#'    [heart_meta()] gives the meta-analytic counterpart.
#' 2. [heart_test()] combines the p-value pairs through a four-group mixture
#'    model into the HEART statistic and estimates its composite null
#'    distribution, providing FDR and FWER control.
#'
#' [heart()] runs both steps on individual-level data.
#'
#' @keywords internal
"_PACKAGE"
