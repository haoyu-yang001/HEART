#' Decorrelated inverse-variance meta-analysis from summary statistics
#'
#' When only summary statistics from independent target and source studies are
#' available, the inverse-variance meta-analytic estimator
#' \eqn{\hat\beta^{meta}_j = w^{tar}_j\hat\beta^{tar}_j + w^{src}_j\hat\beta^{src}_j}
#' can serve as the auxiliary estimator. Because the two studies are independent,
#' \eqn{\hat\rho_j = w^{tar}_j} and the decorrelated estimator is
#' \eqn{\hat\beta^{de}_j = w^{src}_j\hat\beta^{src}_j} with variance
#' \eqn{(w^{src}_j)^2 v^{src}_j}; its p-value therefore coincides with the
#' source-study p-value. The output can be passed to [heart_test()] with
#' `p_de = p_de_meta`.
#'
#' @param beta_tar,se_tar Target-study estimates and standard errors.
#' @param beta_src,se_src Source-study estimates and standard errors (aligned
#'   to the same effect allele).
#' @param snp Optional variant identifiers.
#' @return A data.frame with the target-only, meta-analytic and decorrelated
#'   estimates, standard errors and two-sided p-values.
#' @examples
#' heart_meta(beta_tar = c(0.1, 0), se_tar = c(0.05, 0.05),
#'            beta_src = c(0.08, 0.05), se_src = c(0.01, 0.01))
#' @export
heart_meta <- function(beta_tar, se_tar, beta_src, se_src, snp = NULL) {
  n <- length(beta_tar)
  if (!all(lengths(list(se_tar, beta_src, se_src)) == n)) {
    stop("All inputs must have the same length.", call. = FALSE)
  }
  if (is.null(snp)) snp <- seq_len(n)
  vt <- se_tar^2
  vs <- se_src^2
  ok <- is.finite(beta_tar) & is.finite(vt) & is.finite(beta_src) & is.finite(vs) &
    vt > 0 & vs > 0
  w_tar <- ifelse(ok, (1 / vt) / (1 / vt + 1 / vs), NA_real_)
  w_src <- 1 - w_tar
  beta_meta <- w_tar * beta_tar + w_src * beta_src
  var_meta <- 1 / (1 / vt + 1 / vs)
  beta_de <- beta_meta - w_tar * beta_tar
  var_de <- w_src^2 * vs
  data.frame(
    snp = snp,
    beta_tar = beta_tar, se_tar = se_tar, p_tar = .wald_p(beta_tar, vt),
    beta_meta = beta_meta, se_meta = sqrt(var_meta), p_meta = .wald_p(beta_meta, var_meta),
    rho_meta = w_tar,
    beta_de_meta = beta_de, se_de_meta = sqrt(var_de), p_de_meta = .wald_p(beta_de, var_de),
    stringsAsFactors = FALSE
  )
}
