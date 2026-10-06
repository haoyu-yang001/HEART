#' Per-variant target-only, auxiliary and decorrelated estimators
#'
#' Step 1 of HEART. For every variant (column of `G`) this function computes
#' the target-only estimator \eqn{\hat\beta^{tar}_j}, one or more auxiliary
#' estimators \eqn{\hat\beta^{aux}_j} that borrow information from the source
#' sample, their joint sandwich covariance, and the decorrelated estimator
#' \deqn{\hat\beta^{de}_j = \hat\beta^{aux}_j - \hat\rho_j \hat\beta^{tar}_j,
#'   \qquad \hat\rho_j = \widehat{Cov}(\hat\beta^{aux}_j, \hat\beta^{tar}_j) /
#'   \widehat{Var}(\hat\beta^{tar}_j),}
#' whose Wald p-value \eqn{p_{de,j}} is asymptotically independent of the
#' target-only p-value \eqn{p_{tar,j}}. The pair \eqn{(p_{tar,j}, p_{de,j})}
#' is the input of [heart_test()].
#'
#' All estimators are M-estimators of the working model
#' \eqn{E(Y | G_j, X) = \mu(\beta_0 + G_j\beta + X^\top\beta_X)} with the
#' identity link (`family = "gaussian"`) or the logit link
#' (`family = "binomial"`):
#' \describe{
#'   \item{`"pooled"`}{Regression on the pooled target + source sample (the
#'     auxiliary estimator used for the main analyses of the HEART paper).}
#'   \item{`"dr"`}{Doubly robust (AIPW) estimator that combines a surrogate-based
#'     synthetic outcome `Y_hat` with the density ratio `c_i`
#'     (see [heart_density_ratio()]). Requires `Y_hat` and `c_i`.}
#'   \item{`"source"`}{Source-only regression. It is independent of the
#'     target-only estimator, so the decorrelated estimator equals itself.}
#' }
#' Variances are heteroscedasticity-robust sandwich variances; covariances
#' between estimators use the joint sandwich formula on the shared sample.
#' For the Gaussian family the computations use the Frisch-Waugh-Lovell
#' representation and are vectorised over blocks of variants, so the cost per
#' variant is linear in the number of covariates.
#'
#' @param Y Numeric phenotype vector of length n (0/1 for `family = "binomial"`).
#' @param G Genotype matrix (n x J), a numeric vector (one variant), or any
#'   object supporting `G[, j, drop = FALSE]` and `ncol()` such as a
#'   `BEDMatrix`. Missing genotypes are mean-imputed within each variant and
#'   each population.
#' @param D Population indicator of length n: 1 = target, 0 = source.
#' @param X Optional covariates to adjust for (matrix or data.frame, n rows);
#'   an intercept is always added. For multi-ancestry data include genetic
#'   principal components so that the pooled estimator is not confounded by
#'   population stratification.
#' @param aux Character vector of auxiliary estimators to compute, any of
#'   `"pooled"`, `"dr"`, `"source"`. The first one is used by [heart()].
#' @param family `"gaussian"` (identity link) or `"binomial"` (logit link).
#' @param Y_hat Synthetic (predicted) outcome for every individual, required for
#'   `aux = "dr"`. It should be produced by a model trained on data that are not
#'   used here (e.g. a separate split), using surrogates only.
#' @param c_i Density ratio \eqn{c_i = P(D=1 | S_i, Z_i) / P(D=1)} for every
#'   individual, required for `aux = "dr"`; see [heart_density_ratio()].
#' @param block_size Number of variants processed together (Gaussian family).
#'   Defaults to a value that keeps each block below roughly 150 MB.
#' @param verbose Print progress.
#'
#' @return A data.frame with one row per variant and columns
#'   \describe{
#'     \item{`snp`}{Variant identifier (column names of `G` if present).}
#'     \item{`beta_tar`, `se_tar`, `p_tar`}{Target-only estimate, standard error
#'       and two-sided p-value.}
#'     \item{`beta_<aux>`, `se_<aux>`, `p_<aux>`}{Auxiliary estimate(s).}
#'     \item{`cov_<aux>`, `rho_<aux>`}{Estimated covariance with the target-only
#'       estimator and the projection coefficient.}
#'     \item{`beta_de_<aux>`, `se_de_<aux>`, `p_de_<aux>`}{Decorrelated estimate.}
#'   }
#'   Variants that are monomorphic in the target sample get `NA`.
#'   The attribute `"aux"` stores the auxiliary estimators that were computed.
#'
#' @seealso [heart_test()], [heart()], [heart_density_ratio()]
#' @examples
#' sim <- heart_simulate(n_tar = 1000, n_src = 5000, J = 200, seed = 1)
#' fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X, aux = c("pooled", "source"))
#' head(fit)
#' @export
heart_assoc <- function(Y, G, D, X = NULL,
                        aux = "pooled",
                        family = c("gaussian", "binomial"),
                        Y_hat = NULL, c_i = NULL,
                        block_size = NULL, verbose = FALSE) {
  family <- match.arg(family)
  aux <- unique(match.arg(aux, c("pooled", "dr", "source"), several.ok = TRUE))

  Y <- as.numeric(Y)
  n <- length(Y)
  if (anyNA(Y)) stop("`Y` contains NA; remove those individuals first.", call. = FALSE)
  D <- .check_D(D, n)
  if (family == "binomial" && !all(Y %in% c(0, 1))) {
    stop("`family = \"binomial\"` requires a 0/1 phenotype.", call. = FALSE)
  }
  X1 <- .covariate_design(X, n)
  if (anyNA(X1)) stop("`X` contains NA; remove or impute first.", call. = FALSE)

  if (is.null(dim(G))) G <- matrix(G, ncol = 1L)
  if (nrow(G) != n) stop("`G` must have n = length(Y) rows.", call. = FALSE)
  J <- ncol(G)
  snp <- colnames(G)
  if (is.null(snp)) snp <- paste0("snp", seq_len(J))

  if ("dr" %in% aux) {
    if (is.null(Y_hat) || is.null(c_i)) {
      stop("`aux = \"dr\"` requires both `Y_hat` and `c_i`.", call. = FALSE)
    }
    Y_hat <- as.numeric(Y_hat)
    c_i <- as.numeric(c_i)
    if (length(Y_hat) != n || length(c_i) != n) {
      stop("`Y_hat` and `c_i` must have length n.", call. = FALSE)
    }
    if (anyNA(Y_hat) || anyNA(c_i)) stop("`Y_hat` and `c_i` must not contain NA.", call. = FALSE)
  }

  pre <- list(Y = Y, D = D, X1 = X1, n = n, aux = aux)
  if ("dr" %in% aux) {
    pi_hat <- mean(D)
    w <- D / pi_hat
    pre$pi_hat <- pi_hat
    pre$w <- w
    pre$v <- c_i - w
    pre$r <- Y - Y_hat
  }

  if (family == "gaussian") {
    pre <- .precompute_gaussian(pre)
    if (is.null(block_size)) block_size <- max(1L, min(J, floor(2e6 / n)))
    starts <- seq.int(1L, J, by = block_size)
    res <- vector("list", length(starts))
    for (k in seq_along(starts)) {
      idx <- starts[k]:min(J, starts[k] + block_size - 1L)
      Gb <- .impute_block(G[, idx, drop = FALSE], D)
      res[[k]] <- .assoc_gaussian_block(Gb, pre)
      if (verbose) message(sprintf("variants %d-%d of %d done", idx[1], max(idx), J))
    }
  } else {
    res <- vector("list", J)
    for (j in seq_len(J)) {
      g <- .impute_block(G[, j, drop = FALSE], D)[, 1L]
      res[[j]] <- .assoc_binomial_one(g, pre)
      if (verbose && j %% 100L == 0L) message(sprintf("variant %d of %d done", j, J))
    }
  }
  out <- do.call(rbind, res)
  out <- data.frame(snp = snp, out, row.names = NULL, check.names = FALSE,
                    stringsAsFactors = FALSE)
  attr(out, "aux") <- aux
  attr(out, "family") <- family
  out
}

# ---------------------------------------------------------------------------
# Gaussian family: closed forms via Frisch-Waugh-Lovell.
#
# For a design Z = (X1, G) and weights a_i, e_2' (sum a Z Z')^{-1} Z_i equals
# Gt_i / s, where Gt = G - X1 (X1' A X1)^{-1} X1' A G is the A-weighted
# residual of G on X1 and s = sum a_i Gt_i^2. This holds for every row i,
# including rows with a_i = 0, so all sandwich (co)variances of the SNP
# coefficient reduce to inner products of n-vectors.
# ---------------------------------------------------------------------------

.precompute_gaussian <- function(pre) {
  X1 <- pre$X1; Y <- pre$Y; D <- pre$D
  pre$PT <- .sym_inv(crossprod(X1, X1 * D))
  pre$Y0T <- Y - drop(X1 %*% (pre$PT %*% crossprod(X1, D * Y)))
  if ("pooled" %in% pre$aux) {
    pre$PP <- .sym_inv(crossprod(X1))
    pre$Y0P <- Y - drop(X1 %*% (pre$PP %*% crossprod(X1, Y)))
  }
  if ("source" %in% pre$aux) {
    D0 <- 1 - D
    pre$PS <- .sym_inv(crossprod(X1, X1 * D0))
    pre$Y0S <- Y - drop(X1 %*% (pre$PS %*% crossprod(X1, D0 * Y)))
  }
  if ("dr" %in% pre$aux) {
    # b_i = w_i Y_i + v_i r_i ; DR solves sum Z_i {w_i (Y_i - Z_i'beta) + v_i r_i} = 0
    pre$b <- pre$w * Y + pre$v * pre$r
    PW <- pre$PT * pre$pi_hat            # (X1' W X1)^{-1} with W = D / pi
    pre$Y0DR <- Y - drop(X1 %*% (PW %*% crossprod(X1, pre$b)))
  }
  pre
}

.assoc_gaussian_block <- function(Gb, pre) {
  n <- pre$n; X1 <- pre$X1; D <- pre$D
  B <- ncol(Gb)
  each <- function(x) rep(x, each = n)

  # ---- target-only ----
  GT <- Gb - X1 %*% (pre$PT %*% crossprod(X1, D * Gb))
  sT <- colSums(D * GT^2)
  bad <- !(sT > 1e-8 * sum(D))
  sT[bad] <- NA_real_
  bT <- colSums(D * GT * pre$Y0T) / sT
  eT <- pre$Y0T - GT * each(bT)
  IFT <- D * GT * eT / each(sT)
  vT <- colSums(IFT^2)

  out <- list(beta_tar = bT, se_tar = sqrt(vT), p_tar = .wald_p(bT, vT))

  for (a in pre$aux) {
    if (a == "pooled") {
      GP <- Gb - X1 %*% (pre$PP %*% crossprod(X1, Gb))
      sP <- colSums(GP^2)
      bA <- colSums(GP * pre$Y0P) / sP
      eP <- pre$Y0P - GP * each(bA)
      IFA <- GP * eP / each(sP)
      vA <- colSums(IFA^2)
      cA <- colSums(IFA * IFT)
    } else if (a == "source") {
      D0 <- 1 - D
      GS <- Gb - X1 %*% (pre$PS %*% crossprod(X1, D0 * Gb))
      sS <- colSums(D0 * GS^2)
      sS[!(sS > 1e-8 * sum(D0))] <- NA_real_
      bA <- colSums(D0 * GS * pre$Y0S) / sS
      eS <- pre$Y0S - GS * each(bA)
      IFA <- D0 * GS * eS / each(sS)
      vA <- colSums(IFA^2)
      cA <- rep(0, B)  # disjoint samples
    } else if (a == "dr") {
      sW <- sT / pre$pi_hat
      bA <- colSums(GT * pre$b) / sW
      eDR <- pre$Y0DR - GT * each(bA)
      IFA <- GT * (pre$w * eDR + pre$v * pre$r) / each(sW)
      vA <- colSums(IFA^2)
      cA <- colSums(IFA * IFT)
    }
    out <- c(out, .decorrelate(bA, vA, cA, bT, vT, a))
  }
  out <- as.data.frame(out, check.names = FALSE)
  out[bad, ] <- NA_real_
  out
}

# Decorrelation of an auxiliary estimator against the target-only estimator.
.decorrelate <- function(bA, vA, cA, bT, vT, name) {
  rho <- ifelse(vT > 0, cA / vT, 0)
  b_de <- bA - rho * bT
  v_de <- vA - ifelse(vT > 0, cA^2 / vT, 0)
  v_de <- pmax(v_de, 0)
  res <- list(bA, sqrt(vA), .wald_p(bA, vA), cA, rho,
              b_de, sqrt(v_de), .wald_p(b_de, v_de))
  names(res) <- c(paste0(c("beta_", "se_", "p_", "cov_", "rho_"), name),
                  paste0(c("beta_de_", "se_de_", "p_de_"), name))
  res
}

# ---------------------------------------------------------------------------
# Binomial family (logit link): per-variant IRLS / Newton.
# ---------------------------------------------------------------------------

.logit_fit <- function(Z, y, wt, start = NULL, maxit = 50L, tol = 1e-10) {
  beta <- if (is.null(start)) rep(0, ncol(Z)) else start
  for (it in seq_len(maxit)) {
    eta <- drop(Z %*% beta)
    mu <- stats::plogis(eta)
    W <- wt * pmax(mu * (1 - mu), 1e-12)
    U <- crossprod(Z, wt * (y - mu))
    A <- crossprod(Z, Z * W)
    step <- drop(.sym_inv(A) %*% U)
    beta <- beta + step
    if (max(abs(step)) < tol) break
  }
  mu <- stats::plogis(drop(Z %*% beta))
  list(beta = beta, mu = mu,
       A = crossprod(Z, Z * (wt * mu * (1 - mu))))
}

.assoc_binomial_one <- function(g, pre) {
  Y <- pre$Y; D <- pre$D
  Z <- cbind(pre$X1[, 1L, drop = FALSE], G = g, pre$X1[, -1L, drop = FALSE])
  na_row <- function() {
    cols <- c("beta_tar", "se_tar", "p_tar")
    for (a in pre$aux) {
      cols <- c(cols, paste0(c("beta_", "se_", "p_", "cov_", "rho_"), a),
                paste0(c("beta_de_", "se_de_", "p_de_"), a))
    }
    as.data.frame(as.list(stats::setNames(rep(NA_real_, length(cols)), cols)),
                  check.names = FALSE)
  }
  if (stats::var(g[D == 1]) < 1e-12) return(na_row())

  # target-only
  ft <- .logit_fit(Z, Y, D)
  hT <- drop(Z %*% .sym_inv(ft$A)[, 2L])
  IFT <- hT * D * (Y - ft$mu)
  bT <- ft$beta[2L]
  vT <- sum(IFT^2)
  out <- list(beta_tar = bT, se_tar = sqrt(vT), p_tar = .wald_p(bT, vT))

  for (a in pre$aux) {
    if (a == "pooled") {
      fa <- .logit_fit(Z, Y, rep(1, length(Y)), start = ft$beta)
      IFA <- drop(Z %*% .sym_inv(fa$A)[, 2L]) * (Y - fa$mu)
      cA <- sum(IFA * IFT)
    } else if (a == "source") {
      D0 <- 1 - D
      if (stats::var(g[D0 == 1]) < 1e-12) return(na_row())
      fa <- .logit_fit(Z, Y, D0, start = ft$beta)
      IFA <- drop(Z %*% .sym_inv(fa$A)[, 2L]) * D0 * (Y - fa$mu)
      cA <- 0
    } else if (a == "dr") {
      fa <- .dr_logit_fit(Z, Y, pre$w, pre$v, pre$r, start = ft$beta)
      IFA <- drop(Z %*% .sym_inv(fa$A)[, 2L]) *
        (pre$w * (Y - fa$mu) + pre$v * pre$r)
      cA <- sum(IFA * IFT)
    }
    bA <- fa$beta[2L]
    vA <- sum(IFA^2)
    out <- c(out, .decorrelate(bA, vA, cA, bT, vT, a))
  }
  as.data.frame(out, check.names = FALSE)
}

# Newton-Raphson for sum_i Z_i {w_i (Y_i - mu_i) + v_i r_i} = 0.
.dr_logit_fit <- function(Z, y, w, v, r, start, maxit = 50L, tol = 1e-10) {
  beta <- start
  aug <- crossprod(Z, v * r)
  for (it in seq_len(maxit)) {
    mu <- stats::plogis(drop(Z %*% beta))
    U <- crossprod(Z, w * (y - mu)) + aug
    A <- crossprod(Z, Z * (w * pmax(mu * (1 - mu), 1e-12)))
    step <- drop(.sym_inv(A) %*% U)
    beta <- beta + step
    if (max(abs(step)) < tol) break
  }
  mu <- stats::plogis(drop(Z %*% beta))
  list(beta = beta, mu = mu, A = crossprod(Z, Z * (w * mu * (1 - mu))))
}
