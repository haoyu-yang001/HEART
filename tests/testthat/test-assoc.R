# Reference implementation with explicit sandwich matrices (identity link).
ref_sandwich <- function(Y, g, D, X, c_i = NULL, Y_hat = NULL) {
  Z <- cbind(1, g, X)
  At <- crossprod(Z, Z * D)
  bt <- solve(At, crossprod(Z, D * Y))
  phi <- Z * c(D * (Y - Z %*% bt))
  Ap <- crossprod(Z)
  bp <- solve(Ap, crossprod(Z, Y))
  psi <- Z * c(Y - Z %*% bp)
  It <- phi %*% solve(At); Ip <- psi %*% solve(Ap)
  out <- c(beta_tar = bt[2], var_tar = sum(It[, 2]^2),
           beta_pooled = bp[2], var_pooled = sum(Ip[, 2]^2),
           cov_pooled = sum(It[, 2] * Ip[, 2]))
  if (!is.null(c_i)) {
    pi <- mean(D); w <- D / pi; v <- c_i - w; r <- Y - Y_hat
    Aw <- crossprod(Z, Z * w)
    bd <- solve(Aw, crossprod(Z, w * Y + v * r))
    psd <- Z * c(w * (Y - Z %*% bd) + v * r)
    Id <- psd %*% solve(Aw)
    out <- c(out, beta_dr = bd[2], var_dr = sum(Id[, 2]^2), cov_dr = sum(It[, 2] * Id[, 2]))
  }
  out
}

test_that("Gaussian estimators match explicit sandwich formulas", {
  sim <- heart_simulate(n_tar = 300, n_src = 1200, J = 6, seed = 42)
  Y_hat <- stats::fitted(stats::lm(sim$Y ~ sim$S + sim$X))
  c_i <- heart_density_ratio(sim$D, cbind(S = sim$S))
  fit <- heart_assoc(sim$Y, sim$G, sim$D, sim$X, aux = c("pooled", "dr", "source"),
                     Y_hat = Y_hat, c_i = c_i, block_size = 4)
  for (j in seq_len(ncol(sim$G))) {
    ref <- ref_sandwich(sim$Y, sim$G[, j], sim$D, sim$X, c_i, Y_hat)
    expect_equal(fit$beta_tar[j], unname(ref["beta_tar"]), tolerance = 1e-8)
    expect_equal(fit$se_tar[j]^2, unname(ref["var_tar"]), tolerance = 1e-8)
    expect_equal(fit$beta_pooled[j], unname(ref["beta_pooled"]), tolerance = 1e-8)
    expect_equal(fit$se_pooled[j]^2, unname(ref["var_pooled"]), tolerance = 1e-8)
    expect_equal(fit$cov_pooled[j], unname(ref["cov_pooled"]), tolerance = 1e-8)
    expect_equal(fit$beta_dr[j], unname(ref["beta_dr"]), tolerance = 1e-8)
    expect_equal(fit$se_dr[j]^2, unname(ref["var_dr"]), tolerance = 1e-8)
    expect_equal(fit$cov_dr[j], unname(ref["cov_dr"]), tolerance = 1e-8)
    rho <- ref["cov_pooled"] / ref["var_tar"]
    expect_equal(fit$beta_de_pooled[j], unname(ref["beta_pooled"] - rho * ref["beta_tar"]),
                 tolerance = 1e-8)
    expect_equal(fit$se_de_pooled[j]^2,
                 unname(ref["var_pooled"] - ref["cov_pooled"]^2 / ref["var_tar"]),
                 tolerance = 1e-8)
  }
  # source-only estimator is independent of the target-only one
  expect_true(all(fit$cov_source == 0))
  expect_equal(fit$p_de_source, fit$p_source)
})

test_that("block size does not change results and monomorphic variants give NA", {
  sim <- heart_simulate(n_tar = 200, n_src = 600, J = 5, seed = 1)
  G <- sim$G
  G[sim$D == 1, 3] <- 1
  a <- heart_assoc(sim$Y, G, sim$D, sim$X, block_size = 1)
  b <- heart_assoc(sim$Y, G, sim$D, sim$X, block_size = 5)
  expect_equal(a, b)
  expect_true(all(is.na(a[3, -1])))
  expect_false(anyNA(a[-3, -1]))
})

test_that("missing genotypes are mean-imputed", {
  sim <- heart_simulate(n_tar = 200, n_src = 600, J = 2, seed = 2)
  G <- sim$G
  G[c(1, 500), 1] <- NA
  Gi <- G
  Gi[c(1, 500), 1] <- mean(G[, 1], na.rm = TRUE)
  expect_equal(heart_assoc(sim$Y, G, sim$D, sim$X), heart_assoc(sim$Y, Gi, sim$D, sim$X))
})

test_that("binomial estimators match glm", {
  sim <- heart_simulate(n_tar = 400, n_src = 1500, J = 3, seed = 7)
  set.seed(1)
  yb <- stats::rbinom(length(sim$Y), 1, stats::plogis(-0.3 + 0.4 * sim$X[, 1]))
  fit <- heart_assoc(yb, sim$G, sim$D, sim$X, family = "binomial")
  for (j in 1:3) {
    gt <- stats::glm(yb ~ sim$G[, j] + sim$X, family = stats::binomial, subset = sim$D == 1)
    gp <- stats::glm(yb ~ sim$G[, j] + sim$X, family = stats::binomial)
    expect_equal(fit$beta_tar[j], unname(stats::coef(gt)[2]), tolerance = 1e-6)
    expect_equal(fit$beta_pooled[j], unname(stats::coef(gp)[2]), tolerance = 1e-6)
  }
})

test_that("input validation", {
  sim <- heart_simulate(n_tar = 50, n_src = 100, J = 2, seed = 1)
  expect_error(heart_assoc(sim$Y, sim$G, sim$D + 1, sim$X), "0/1")
  expect_error(heart_assoc(sim$Y, sim$G, sim$D, sim$X, aux = "dr"), "Y_hat")
  expect_error(heart_assoc(sim$Y, sim$G, sim$D, cbind(sim$X, sim$X[, 1])), "rank deficient")
})
