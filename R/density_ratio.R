#' Density ratio between the target and the pooled population
#'
#' Estimates \eqn{c_i = P(D_i = 1 | W_i) / \hat\pi} with
#' \eqn{\hat\pi = n^{-1}\sum_i D_i} by a logistic regression of the population
#' indicator `D` on the features `W` (typically surrogates `S` and, optionally,
#' covariates). Because \eqn{P(D=1|W)/\pi = f_{target}(W)/f_{pooled}(W)},
#' \eqn{c_i} is the weight that aligns the pooled sample with the target
#' population and is used by the doubly robust auxiliary estimator
#' (`aux = "dr"` in [heart_assoc()]).
#'
#' The logistic model is fitted by iteratively reweighted least squares with a
#' small ridge penalty that is increased automatically if the weighted Gram
#' matrix is numerically singular (e.g. under quasi-separation).
#'
#' @param D 0/1 population indicator (1 = target).
#' @param W Features (matrix or data.frame); factors are expanded by
#'   [stats::model.matrix()].
#' @param quadratic If `TRUE`, squared terms of the numeric columns are added.
#' @param ridge Initial ridge penalty.
#' @param maxit,tol Maximum IRLS iterations and convergence tolerance.
#'
#' @return Numeric vector `c_i` with attributes `coef` (fitted logistic
#'   coefficients), `e_hat` (fitted P(D = 1 | W)) and `pi_hat`.
#' @examples
#' sim <- heart_simulate(n_tar = 500, n_src = 2000, J = 10, seed = 1)
#' c_i <- heart_density_ratio(sim$D, cbind(S = sim$S, sim$X))
#' tapply(c_i, sim$D, mean)
#' @export
heart_density_ratio <- function(D, W, quadratic = FALSE, ridge = 1e-8,
                                maxit = 50L, tol = 1e-8) {
  D <- as.numeric(D)
  if (is.data.frame(W)) {
    W <- stats::model.matrix(~ ., data = W)
    W <- W[, colnames(W) != "(Intercept)", drop = FALSE]
  }
  W <- as.matrix(W)
  if (nrow(W) != length(D)) stop("`W` and `D` must have the same number of rows.", call. = FALSE)
  if (anyNA(W) || anyNA(D)) stop("`W` and `D` must not contain NA.", call. = FALSE)
  if (is.null(colnames(W))) colnames(W) <- paste0("W", seq_len(ncol(W)))
  if (quadratic) {
    W2 <- W^2
    keep <- apply(W, 2L, function(x) length(unique(x)) > 2L)
    W2 <- W2[, keep, drop = FALSE]
    if (ncol(W2) > 0L) colnames(W2) <- paste0(colnames(W2), "^2")
    W <- cbind(W, W2)
  }
  Xd <- cbind(`(Intercept)` = 1, W)

  beta <- numeric(ncol(Xd))
  for (it in seq_len(maxit)) {
    eta <- pmax(pmin(drop(Xd %*% beta), 30), -30)
    mu <- stats::plogis(eta)
    Wt <- pmax(mu * (1 - mu), 1e-12)
    z <- eta + (D - mu) / Wt
    XtWX <- crossprod(Xd, Xd * Wt)
    XtWz <- crossprod(Xd, Wt * z)
    lam <- ridge
    repeat {
      M <- XtWX
      if (lam > 0) diag(M) <- diag(M) + lam
      beta_new <- tryCatch({
        R <- chol(M)
        backsolve(R, forwardsolve(t(R), XtWz))
      }, error = function(e) NA_real_)
      if (all(is.finite(beta_new))) break
      lam <- if (lam == 0) 1e-6 else lam * 10
      if (lam > 1e8) {
        stop("Weighted Gram matrix is singular even with strong ridge penalty; ",
             "check for separation or collinearity in `W`.", call. = FALSE)
      }
    }
    beta_new <- drop(beta_new)
    step <- max(abs(beta_new - beta))
    beta <- beta_new
    if (step < tol) break
  }
  names(beta) <- colnames(Xd)
  e_hat <- stats::plogis(drop(Xd %*% beta))
  pi_hat <- mean(D)
  c_i <- e_hat / pi_hat
  attr(c_i, "coef") <- beta
  attr(c_i, "e_hat") <- e_hat
  attr(c_i, "pi_hat") <- pi_hat
  c_i
}
