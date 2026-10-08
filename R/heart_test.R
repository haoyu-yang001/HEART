#' Estimate the four-group proportions from a pair of p-value sequences
#'
#' Estimates \eqn{\pi_{00}, \pi_{01}, \pi_{10}, \pi_{11}} of the four-group
#' model for \eqn{(p_{tar,j}, p_{de,j})}, where the first index refers to
#' \eqn{\beta^{tar}_j} and the second to \eqn{\beta^{de}_j}
#' (0 = null, 1 = non-null). The marginal and joint null proportions are
#' estimated by Storey-type ratios
#' \deqn{\hat\pi_{0\cdot}(\lambda) = \frac{\#\{p_{tar,j} > \lambda\}}{J(1-\lambda)},\quad
#'       \hat\pi_{\cdot0}(\lambda) = \frac{\#\{p_{de,j} > \lambda\}}{J(1-\lambda)},\quad
#'       \hat\pi_{00}(\lambda) = \frac{\#\{p_{tar,j} > \lambda, p_{de,j} > \lambda\}}{J(1-\lambda)^2},}
#' with \eqn{\hat\pi_{01} = (\hat\pi_{0\cdot} - \hat\pi_{00}) \vee 0} and
#' \eqn{\hat\pi_{10} = (\hat\pi_{\cdot0} - \hat\pi_{00}) \vee 0}. A marginal
#' proportion is set to one when a one-sided Kolmogorov-Smirnov test does not
#' reject uniformity. The tuning parameter \eqn{\lambda} is chosen among
#' `lambdas` by the QQ-slope criterion for the maximum p-value (as in the
#' HDMT procedure of Dai, Stanford and LeBlanc, 2022).
#'
#' @param p_tar,p_de Target-only and decorrelated p-values (same length).
#' @param lambdas Candidate values of \eqn{\lambda}.
#' @param ks_level Level of the Kolmogorov-Smirnov uniformity check.
#'
#' @return A list with `pi00`, `pi01`, `pi10`, `pi11`, the marginal null
#'   proportions `pi0_tar` (\eqn{\hat\pi_{0\cdot}}) and `pi0_de`
#'   (\eqn{\hat\pi_{\cdot0}}), and the selected `lambda`.
#' @export
heart_pi_estimate <- function(p_tar, p_de, lambdas = c(0.5, 0.6, 0.7, 0.8),
                              ks_level = 0.05) {
  P <- cbind(as.numeric(p_tar), as.numeric(p_de))
  ok <- stats::complete.cases(P)
  if (any(!ok)) warning("Removing ", sum(!ok), " pairs with NA p-values.", call. = FALSE)
  P <- P[ok, , drop = FALSE]
  if (nrow(P) < 10L) stop("Too few valid p-value pairs.", call. = FALSE)

  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- vapply(pcut, function(u) mean(P[, 1] >= u) / (1 - u), numeric(1))
  frac2 <- vapply(pcut, function(u) mean(P[, 2] >= u) / (1 - u), numeric(1))
  frac12 <- vapply(pcut, function(u) mean(P[, 2] >= u & P[, 1] >= u) / (1 - u)^2,
                   numeric(1))
  ks_p <- function(x) {
    suppressWarnings(stats::ks.test(x, "punif", 0, 1, alternative = "greater")$p.value)
  }
  unif1 <- ks_p(P[, 1]) > ks_level
  unif2 <- ks_p(P[, 2]) > ks_level

  pmax_sorted <- sort(pmax(P[, 1], P[, 2]))
  nmed <- nrow(P)
  nnulls <- sum(pmax_sorted > 0.8)

  out <- matrix(NA_real_, length(lambdas), 5L,
                dimnames = list(NULL, c("a10", "a01", "a00", "a1", "a2")))
  slope <- rep(Inf, length(lambdas))
  for (l in seq_along(lambdas)) {
    lambda <- lambdas[l]
    a00 <- min(frac12[pcut >= lambda][1], 1)
    a1 <- if (unif1) 1 else min(frac1[pcut >= lambda][1], 1)
    a2 <- if (unif2) 1 else min(frac2[pcut >= lambda][1], 1)
    a01 <- 0; a10 <- 0
    if (a00 == 1 || (a1 == 1 && a2 == 1)) {
      a00 <- 1
    } else if (a1 == 1 && a2 != 1) {
      a01 <- max(0, a1 - a00)
      a00 <- 1 - a01
    } else if (a1 != 1 && a2 == 1) {
      a10 <- max(0, a2 - a00)
      a00 <- 1 - a10
    } else {
      a10 <- max(0, a2 - a00)
      a01 <- max(0, a1 - a00)
      if (1 - a00 - a01 - a10 < 0) {
        a10 <- 1 - a1
        a01 <- 1 - a2
        a00 <- 1 - a10 - a01
      }
    }
    out[l, ] <- c(a10, a01, a00, a1, a2)

    # QQ slope of the maximum p-value against its null distribution
    b <- a01 + a10
    a <- 1 - b
    if (a <= 0) {
      slope[l] <- -Inf
      break
    }
    if (nnulls >= 2L) {
      i <- seq_len(nmed)
      pexp <- (-b + sqrt(b^2 + 4 * a * i / nmed)) / (2 * a)
      tail_idx <- (nmed - nnulls + 1L):nmed
      xx <- -log10(pexp[tail_idx])
      yy <- -log10(pmax_sorted[tail_idx])
      slope[l] <- unname(stats::coef(stats::lm(yy ~ xx - 1))[1])
    }
  }
  best <- which.min(slope)
  if (length(best) == 0L) best <- 1L
  r <- out[best, ]
  pi00 <- unname(r["a00"]); pi01 <- unname(r["a01"]); pi10 <- unname(r["a10"])
  list(pi00 = pi00, pi01 = pi01, pi10 = pi10,
       pi11 = max(0, 1 - pi00 - pi01 - pi10),
       pi0_tar = unname(r["a1"]), pi0_de = unname(r["a2"]),
       lambda = lambdas[best])
}

#' Estimate the left-tail parameters of the alternative p-value density
#'
#' Under the two-group model \eqn{f(p) = \pi_0 + (1-\pi_0) f_1(p)} with the
#' local approximation \eqn{f_1(p) \approx C\alpha p^{\alpha-1}} near zero,
#' \eqn{(\alpha, C)} solve the conditional moment equations
#' \deqn{J^{-1}\sum_j \log(p_j) I(p_j \le r) - \pi_0(r\log r - r) -
#'       C(1-\pi_0) r^{\alpha}(\log r - \alpha^{-1}) = 0,}
#' \deqn{J^{-1}\sum_j I(p_j \le r) = \pi_0 r + (1-\pi_0) C r^{\alpha},}
#' which have a closed-form solution. The threshold \eqn{r} is the
#' \eqn{m}-th smallest p-value with \eqn{m = \lceil k\sqrt{J}\rceil}.
#' The solutions are truncated to `alpha_bounds` and `C_bounds`.
#'
#' @param p Vector of p-values.
#' @param pi0 Estimated null proportion of `p`.
#' @param k Multiplier defining the tail size \eqn{m = k\sqrt{J}}.
#'   `k = 1` corresponds to the \eqn{J^{-1/2}} quantile.
#' @param pi0_cap Upper cap applied to `pi0` (avoids division by zero when the
#'   p-values look uniform).
#' @param alpha_bounds,C_bounds Truncation ranges for \eqn{\alpha} and \eqn{C}.
#'
#' @return Named numeric vector `c(alpha = , C = , r = )`.
#' @export
heart_tail_estimate <- function(p, pi0, k = 50, pi0_cap = 0.996,
                                alpha_bounds = c(1e-5, 1 - 1e-5),
                                C_bounds = c(0.1, 3)) {
  p <- as.numeric(p[!is.na(p)])
  J <- length(p)
  m <- min(J, max(1L, ceiling(k * sqrt(J))))
  r <- sort(p, partial = m)[m]
  r <- min(max(r, .Machine$double.xmin), 1 - 1e-12)
  pi0c <- min(pi0, pi0_cap)

  L <- sum(log(p[p <= r])) / J
  Fr <- mean(p <= r)
  excess <- Fr - pi0c * r               # = (1 - pi0) C r^alpha
  if (!(excess > 0)) {
    alpha <- alpha_bounds[2]
    C <- C_bounds[1]
  } else {
    inv_alpha <- log(r) - (L - pi0c * (r * log(r) - r)) / excess
    alpha <- if (inv_alpha <= 1) alpha_bounds[2] else 1 / inv_alpha
    alpha <- min(max(alpha, alpha_bounds[1]), alpha_bounds[2])
    C <- excess / ((1 - pi0c) * r^alpha)
  }
  C <- min(max(C, C_bounds[1]), C_bounds[2])
  c(alpha = alpha, C = C, r = r)
}

# Constants of the simplified statistic Gamma(p, q) = p^(1-a1) (k1 q^(1-a2) + k2).
.heart_consts <- function(par) {
  with(par, list(
    a = 1 / (1 - alpha1),          # exponent in G_t
    b = 1 - alpha2,
    one_m_a1 = 1 - alpha1,
    logk1 = log(pi00) - log(C1) - log(C2) - log(alpha1) - log(alpha2),
    logk2 = log(pi01) - log(C1) - log(alpha1)
  ))
}

# log(k1 q^b + k2)
.log_denom <- function(logq, cst) .log_add_exp(cst$logk1 + cst$b * logq, cst$logk2)

#' HEART test statistic
#'
#' Computes the simplified HEART statistic (with \eqn{\pi_{10}} set to zero)
#' \deqn{T_j = \hat\pi_{00}\hat C_1^{-1}\hat C_2^{-1}\hat\alpha_1^{-1}\hat\alpha_2^{-1}
#'   p_{tar,j}^{1-\hat\alpha_1} p_{de,j}^{1-\hat\alpha_2} +
#'   \hat\pi_{01}\hat C_1^{-1}\hat\alpha_1^{-1} p_{tar,j}^{1-\hat\alpha_1}.}
#' Small values are evidence against \eqn{H_{0,j}: \beta^{tar}_j = 0}.
#'
#' @param p_tar,p_de Target-only and decorrelated p-values.
#' @param par List with `pi00`, `pi01`, `alpha1`, `alpha2`, `C1`, `C2`
#'   (as returned in `$par` by [heart_test()]).
#' @param log Return \eqn{\log T_j} (recommended for very small p-values).
#' @return Numeric vector of statistics.
#' @export
heart_statistic <- function(p_tar, p_de, par, log = FALSE) {
  cst <- .heart_consts(par)
  lt <- cst$one_m_a1 * log(p_tar) + .log_denom(log(p_de), cst)
  if (log) lt else exp(lt)
}

# F00(t) = int_0^1 G_t(q) dq for a single log t (exact, adaptive quadrature).
# G_t(q) = 1 for q <= q0 and t^a (k1 q^b + k2)^(-a) beyond, so
# F00 = q0 + int_{log q0}^0 exp{s + a log t - a log(k1 e^(b s) + k2)} ds.
# Integrating over s = log q keeps the integrand a smooth combination of
# exponentials for every (alpha1, alpha2), including alpha2 close to 1.
.F00_one <- function(logt, cst) {
  a <- cst$a; b <- cst$b
  if (!is.finite(cst$logk1)) {                       # pi00 = 0: G_t constant in q
    return(min(1, exp(a * (logt - cst$logk2))))
  }
  if (is.finite(cst$logk2) && logt <= cst$logk2) {
    s0 <- -Inf
  } else {
    log_ustar <- logt + log1p(-exp(cst$logk2 - logt)) - cst$logk1   # log(q0^b)
    if (log_ustar >= 0) return(1)
    s0 <- log_ustar / b
  }
  mass0 <- if (is.finite(s0)) exp(s0) else 0
  if (s0 >= 0) return(1)
  logf <- function(s) s + a * logt - a * .log_add_exp(cst$logk1 + b * s, cst$logk2)
  # logf is concave in s, so {s: logf(s) >= M - 60} is an interval around the
  # maximiser; outside it the integrand is negligible (relative size < e^-60).
  cand <- c(0, if (is.finite(s0)) s0)
  ab <- a * b
  if (ab > 1 && is.finite(cst$logk2)) {
    sc <- (cst$logk2 - cst$logk1 - log(ab - 1)) / b
    if (sc > s0 && sc < 0) cand <- c(cand, sc)
  }
  lf <- logf(cand)
  smax <- cand[which.max(lf)]
  M <- max(lf)
  cut <- M - 60
  g <- function(s) logf(s) - cut
  lo <- if (is.finite(s0) && g(s0) >= 0) s0 else {
    left <- if (is.finite(s0)) s0 else smax - 1
    step <- 1
    while (!is.finite(s0) && g(left) > 0) { step <- step * 2; left <- smax - step }
    stats::uniroot(g, c(left, smax), tol = 1e-10)$root
  }
  hi <- if (g(0) >= 0) 0 else stats::uniroot(g, c(smax, 0), tol = 1e-10)$root
  f <- function(s) exp(logf(s) - M)
  pieces <- unique(c(lo, smax, hi))
  val <- 0
  for (i in seq_len(length(pieces) - 1L)) {
    l <- pieces[i]; h <- pieces[i + 1L]
    if (h <= l) next
    val <- val + tryCatch(
      stats::integrate(f, l, h, rel.tol = 1e-10, abs.tol = 0, subdivisions = 2000L)$value,
      error = function(e) stats::integrate(f, l, h, subdivisions = 2000L)$value)
  }
  min(1, mass0 + exp(M) * val)
}

# F00 at many (sorted, unique) log t values: exact for the smallest `n_exact`,
# log-log linear interpolation on a grid for the rest (F00 is smooth there).
.F00_vec <- function(logt, cst, exact = FALSE, n_exact = 5000L, n_grid = 500L) {
  nt <- length(logt)
  if (exact || nt <= n_exact + n_grid) {
    return(vapply(logt, .F00_one, numeric(1), cst = cst))
  }
  out <- numeric(nt)
  out[seq_len(n_exact)] <- vapply(logt[seq_len(n_exact)], .F00_one, numeric(1), cst = cst)
  grid <- seq(logt[n_exact], logt[nt], length.out = n_grid)
  Fg <- vapply(grid, .F00_one, numeric(1), cst = cst)
  rest <- (n_exact + 1L):nt
  out[rest] <- exp(stats::approx(grid, log(Fg), xout = logt[rest], rule = 2)$y)
  pmin(out, 1)
}

# Q_hat(t) = J^{-1} sum_j G_t(p_de_j), exact, O(log J) per t after sorting.
.Q_hat <- function(logt, p_de, cst) {
  a <- cst$a
  h <- sort(.log_denom(log(p_de), cst))      # G_t(q_j) = min(1, exp(a (log t - h_j)))
  J <- length(h)
  y <- rev(-a * h)                           # ascending
  cl <- numeric(J)                           # cumulative logsumexp of y
  cl[1] <- y[1]
  if (J > 1L) for (i in 2:J) {
    cl[i] <- y[i] + log1p(exp(cl[i - 1] - y[i]))
  }
  suffix <- rev(cl)                          # suffix[k] = log sum_{i >= k} exp(-a h_i)
  n1 <- findInterval(logt, h)                # #{h_j <= log t}
  tail_sum <- numeric(length(logt))
  has_tail <- n1 < J
  tail_sum[has_tail] <- exp(a * logt[has_tail] + suffix[n1[has_tail] + 1L])
  pmin(1, (n1 + tail_sum) / J)
}

#' Composite null distribution of the HEART statistic
#'
#' Evaluates, at the statistic values `t`, the global-null CDF
#' \eqn{F_{00}(t) = \int_0^1 G_t(q)dq}, its empirical counterpart
#' \eqn{\hat Q(t) = J^{-1}\sum_j G_t(p_{de,j})}, the reconstructed
#' \eqn{\hat F_{01}(t) = \{\hat Q(t) - \hat\pi_{00}F_{00}(t)\}/(\hat\pi_{01}+\hat\pi_{11})}
#' (truncated to \[0, 1\] and made nondecreasing in \eqn{t} by a running maximum;
#' set to 0 if \eqn{\hat\pi_{01}+\hat\pi_{11} = 0}), the null numerator
#' \eqn{\hat N(t) = \hat\pi_{00}F_{00}(t) + \hat\pi_{01}\hat F_{01}(t)} and the
#' composite null CDF \eqn{\hat F_0(t) = \hat N(t) / (\hat\pi_{00} + \hat\pi_{01})},
#' where \eqn{G_t(q) = \min\{1, (t/(\kappa_1 q^{1-\alpha_2} + \kappa_2))^{1/(1-\alpha_1)}\}}.
#'
#' @param log_t Values of \eqn{\log t}.
#' @param p_de All decorrelated p-values (used for \eqn{\hat Q}).
#' @param par Parameter list as in [heart_statistic()], also containing
#'   `pi11` (or `pi10`, from which \eqn{\pi_{11} = 1 - \pi_{00} - \pi_{01} - \pi_{10}}).
#' @param exact If `FALSE`, \eqn{F_{00}} is computed exactly at the `n_exact`
#'   smallest values and interpolated (log-log) on an `n_grid` grid elsewhere.
#' @param n_exact,n_grid See `exact`.
#' @return A data.frame with columns `log_t`, `F00`, `Q`, `F01`, `N`, `F0`,
#'   sorted by `log_t`.
#' @export
heart_null_cdf <- function(log_t, p_de, par, exact = FALSE,
                           n_exact = 5000L, n_grid = 500L) {
  if (par$pi00 + par$pi01 <= 0) stop("pi00 + pi01 must be positive.", call. = FALSE)
  cst <- .heart_consts(par)
  lt <- sort(unique(log_t))
  F00 <- .F00_vec(lt, cst, exact = exact, n_exact = n_exact, n_grid = n_grid)
  Q <- .Q_hat(lt, p_de, cst)
  pi11 <- if (!is.null(par$pi11)) par$pi11 else max(0, 1 - par$pi00 - par$pi01 - par$pi10)
  den <- par$pi01 + pi11
  F01 <- if (den <= 0) rep(0, length(lt)) else cummax(pmin(pmax((Q - par$pi00 * F00) / den, 0), 1))
  N <- par$pi00 * F00 + par$pi01 * F01
  F0 <- pmin(1, N / (par$pi00 + par$pi01))
  data.frame(log_t = lt, F00 = F00, Q = Q, F01 = F01, N = N, F0 = F0)
}

#' HEART test: adaptive transfer via null decomposition
#'
#' Step 2 of HEART. Given target-only p-values \eqn{p_{tar,j}} and decorrelated
#' p-values \eqn{p_{de,j}} for \eqn{J} variants (typically genome-wide), this
#' function
#' 1. estimates the four-group proportions ([heart_pi_estimate()]),
#' 2. estimates the tail parameters \eqn{(\alpha_1, C_1)} from \eqn{p_{tar}}
#'    and \eqn{(\alpha_2, C_2)} from \eqn{p_{de}} ([heart_tail_estimate()]),
#' 3. computes the HEART statistic \eqn{T_j} ([heart_statistic()]),
#' 4. estimates its composite null distribution ([heart_null_cdf()]), and
#' 5. returns HEART p-values, FDR q-values and FWER-adjusted p-values.
#'
#' The FDR procedure rejects \eqn{\{j: T_j \le \hat t_\alpha\}} with
#' \eqn{\hat t_\alpha = \sup\{t: \hat N(t)/\bar F(t) \le \alpha\}}; equivalently
#' `q_heart <= alpha`. The FWER procedure rejects
#' \eqn{\{j: \hat N(T_j) \le \alpha/J\}}; equivalently `padj_fwer <= alpha`.
#' The per-variant HEART p-value is \eqn{\hat F_0(T_j)}, the estimated
#' composite-null CDF evaluated at the observed statistic, which can be used
#' e.g. for Manhattan plots with a genome-wide threshold.
#'
#' The parameters must be estimated from all variants jointly; when variants
#' are analysed in chunks with [heart_assoc()], combine the chunks before
#' calling `heart_test()`.
#'
#' @param p_tar,p_de Target-only and decorrelated p-values (same length).
#' @param alpha Nominal FDR / FWER level used for the `reject_*` columns.
#' @param snp Optional variant identifiers.
#' @param par Optional list of model parameters (`pi00`, `pi01`, `pi10`,
#'   `alpha1`, `alpha2`, `C1`, `C2`, and optionally `pi11`) to bypass
#'   estimation.
#' @param tail_k,pi0_cap,alpha_bounds,C_bounds Passed to
#'   [heart_tail_estimate()].
#' @param lambdas Passed to [heart_pi_estimate()].
#' @param p_min P-values are bounded below by `p_min` to avoid `log(0)`.
#' @param exact,n_exact,n_grid Passed to [heart_null_cdf()].
#'
#' @return An object of class `"heart"`: a list with
#'   \describe{
#'     \item{`results`}{data.frame with `snp`, `p_tar`, `p_de`, `log_T`,
#'       `p_heart`, `q_heart`, `padj_fwer`, `reject_fdr`, `reject_fwer`.}
#'     \item{`par`}{Estimated parameters.}
#'     \item{`alpha`, `n_tests`}{Level and number of valid tests.}
#'   }
#' @examples
#' sim <- heart_simulate(n_tar = 1500, n_src = 6000, J = 1000, seed = 1)
#' fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X)
#' res <- heart_test(fit$p_tar, fit$p_de_pooled, snp = fit$snp)
#' res
#' @export
heart_test <- function(p_tar, p_de, alpha = 0.05, snp = NULL, par = NULL,
                       tail_k = 50, pi0_cap = 0.996,
                       alpha_bounds = c(1e-5, 1 - 1e-5), C_bounds = c(0.1, 3),
                       lambdas = c(0.5, 0.6, 0.7, 0.8),
                       p_min = 1e-300, exact = FALSE,
                       n_exact = 5000L, n_grid = 500L) {
  p_tar <- as.numeric(p_tar); p_de <- as.numeric(p_de)
  if (length(p_tar) != length(p_de)) stop("`p_tar` and `p_de` must have the same length.", call. = FALSE)
  J_all <- length(p_tar)
  if (is.null(snp)) snp <- seq_len(J_all)
  ok <- is.finite(p_tar) & is.finite(p_de)
  if (any(!ok)) warning(sum(!ok), " variants with missing p-values are excluded.", call. = FALSE)
  pt <- pmin(pmax(p_tar[ok], p_min), 1)
  pd <- pmin(pmax(p_de[ok], p_min), 1)
  J <- length(pt)

  if (is.null(par)) {
    pis <- heart_pi_estimate(pt, pd, lambdas = lambdas)
    tl1 <- heart_tail_estimate(pt, pis$pi0_tar, k = tail_k, pi0_cap = pi0_cap,
                               alpha_bounds = alpha_bounds, C_bounds = C_bounds)
    tl2 <- heart_tail_estimate(pd, pis$pi0_de, k = tail_k, pi0_cap = pi0_cap,
                               alpha_bounds = alpha_bounds, C_bounds = C_bounds)
    par <- list(pi00 = pis$pi00, pi01 = pis$pi01, pi10 = pis$pi10, pi11 = pis$pi11,
                alpha1 = unname(tl1["alpha"]), alpha2 = unname(tl2["alpha"]),
                C1 = unname(tl1["C"]), C2 = unname(tl2["C"]),
                pi0_tar = pis$pi0_tar, pi0_de = pis$pi0_de, lambda = pis$lambda)
  } else {
    need <- c("pi00", "pi01", "pi10", "alpha1", "alpha2", "C1", "C2")
    if (!all(need %in% names(par))) {
      stop("`par` must contain: ", paste(need, collapse = ", "), call. = FALSE)
    }
    if (is.null(par$pi11)) par$pi11 <- max(0, 1 - par$pi00 - par$pi01 - par$pi10)
  }

  log_T <- heart_statistic(pt, pd, par, log = TRUE)
  nc <- heart_null_cdf(log_T, pd, par, exact = exact, n_exact = n_exact, n_grid = n_grid)
  m <- match(log_T, nc$log_t)
  N_j <- nc$N[m]
  p_heart <- nc$F0[m]
  padj_fwer <- pmin(1, J * N_j)

  # q-values: min over l >= rank of N(T_(l)) J / l, with l = #{T <= T_(l)}
  k_j <- rank(log_T, ties.method = "max")
  mfdr <- pmin(1, N_j * J / k_j)
  o <- order(log_T, decreasing = TRUE)
  q_heart <- numeric(J)
  q_heart[o] <- cummin(mfdr[o])

  res <- data.frame(snp = snp, p_tar = p_tar, p_de = p_de,
                    log_T = NA_real_, p_heart = NA_real_, q_heart = NA_real_,
                    padj_fwer = NA_real_, reject_fdr = NA, reject_fwer = NA,
                    stringsAsFactors = FALSE)
  res$log_T[ok] <- log_T
  res$p_heart[ok] <- p_heart
  res$q_heart[ok] <- q_heart
  res$padj_fwer[ok] <- padj_fwer
  res$reject_fdr[ok] <- q_heart <= alpha
  res$reject_fwer[ok] <- padj_fwer <= alpha

  structure(list(results = res, par = par, alpha = alpha, n_tests = J,
                 null_cdf = nc),
            class = "heart")
}

#' @export
print.heart <- function(x, ...) {
  r <- x$results
  p <- x$par
  cat("HEART: Heterogeneity-aware Empowerment via Adaptive Robust Transfer\n")
  cat(sprintf("  tests: %d\n", x$n_tests))
  cat(sprintf("  four-group proportions: pi00 = %.4f, pi01 = %.4f, pi10 = %.4f, pi11 = %.4f\n",
              p$pi00, p$pi01, p$pi10, p$pi11))
  cat(sprintf("  tail parameters: alpha1 = %.3f, C1 = %.3f, alpha2 = %.3f, C2 = %.3f\n",
              p$alpha1, p$C1, p$alpha2, p$C2))
  cat(sprintf("  discoveries at level %g: FDR = %d, FWER = %d (target-only BH = %d, Bonferroni = %d)\n",
              x$alpha, sum(r$reject_fdr, na.rm = TRUE), sum(r$reject_fwer, na.rm = TRUE),
              sum(stats::p.adjust(r$p_tar, "BH") <= x$alpha, na.rm = TRUE),
              sum(stats::p.adjust(r$p_tar, "bonferroni") <= x$alpha, na.rm = TRUE)))
  invisible(x)
}

#' @export
summary.heart <- function(object, n = 10L, ...) {
  r <- object$results
  r <- r[order(r$log_T), , drop = FALSE]
  print(object)
  cat("\nTop variants:\n")
  print(utils::head(r, n), row.names = FALSE)
  invisible(object)
}
