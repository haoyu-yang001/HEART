#' HEART: end-to-end analysis from individual-level data
#'
#' Convenience wrapper that runs [heart_assoc()] (Step 1: target-only,
#' auxiliary and decorrelated estimators) followed by [heart_test()]
#' (Step 2: adaptive transfer via null decomposition) using the first
#' auxiliary estimator in `aux`.
#'
#' For genome-wide analyses that do not fit in memory, run [heart_assoc()] on
#' chunks of variants (e.g. per chromosome or per array job), row-bind the
#' results, and call [heart_test()] once on all variants.
#'
#' @inheritParams heart_assoc
#' @inheritParams heart_test
#' @param ... Further arguments passed to [heart_test()].
#' @return An object of class `"heart"` (see [heart_test()]) with an
#'   additional element `assoc` holding the output of [heart_assoc()].
#' @examples
#' sim <- heart_simulate(n_tar = 1500, n_src = 6000, J = 1000, seed = 1)
#' res <- heart(sim$Y, sim$G, sim$D, sim$X, aux = "pooled")
#' res
#' # power gain over the target-only analysis on the simulated truth
#' table(HEART = res$results$reject_fdr, truth = sim$truth$beta_tar != 0)
#' @export
heart <- function(Y, G, D, X = NULL, aux = "pooled",
                  family = c("gaussian", "binomial"),
                  Y_hat = NULL, c_i = NULL, alpha = 0.05,
                  block_size = NULL, verbose = FALSE, ...) {
  fit <- heart_assoc(Y = Y, G = G, D = D, X = X, aux = aux, family = family,
                     Y_hat = Y_hat, c_i = c_i, block_size = block_size,
                     verbose = verbose)
  a <- attr(fit, "aux")[1]
  res <- heart_test(fit$p_tar, fit[[paste0("p_de_", a)]], alpha = alpha,
                    snp = fit$snp, ...)
  res$assoc <- fit
  res$aux <- a
  res
}
