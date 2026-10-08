make_par <- function() {
  list(pi00 = 0.9, pi01 = 0.05, pi10 = 0.01, pi11 = 0.04,
       alpha1 = 0.2, alpha2 = 0.3, C1 = 1.1, C2 = 0.8)
}

test_that("F00 agrees with direct integration of G_t over q", {
  par <- make_par()
  cst <- HEART:::.heart_consts(par)
  Gt <- function(q, t) {
    pmin(1, (t / (exp(cst$logk1) * q^cst$b + exp(cst$logk2)))^cst$a)
  }
  for (t in c(1e-3, 0.05, 0.3, 1)) {
    ref <- stats::integrate(Gt, 0, 1, t = t, rel.tol = 1e-12)$value
    expect_equal(HEART:::.F00_one(log(t), cst), ref, tolerance = 1e-7)
  }
  # pi01 = 0 and pi00 = 0 edge cases
  par0 <- modifyList(par, list(pi01 = 0)); c0 <- HEART:::.heart_consts(par0)
  expect_true(is.finite(HEART:::.F00_one(log(1e-4), c0)))
  par1 <- modifyList(par, list(pi00 = 0)); c1 <- HEART:::.heart_consts(par1)
  expect_equal(HEART:::.F00_one(log(1e-3), c1),
               min(1, (1e-3 / exp(c1$logk2))^c1$a))
})

test_that("Q_hat equals the brute-force average of G_t(p_de)", {
  par <- make_par()
  cst <- HEART:::.heart_consts(par)
  set.seed(3)
  pd <- c(stats::runif(500), stats::rbeta(50, 0.1, 1))
  lt <- log(c(1e-4, 1e-2, 0.1, 0.5, 2))
  brute <- sapply(lt, function(l) mean(pmin(1, exp(cst$a * (l - HEART:::.log_denom(log(pd), cst))))))
  expect_equal(HEART:::.Q_hat(lt, pd, cst), brute, tolerance = 1e-12)
})

test_that("tail estimates solve the moment equations", {
  set.seed(9)
  p <- c(stats::runif(9000), stats::rbeta(1000, 0.1, 1))
  pi0 <- 0.9
  est <- heart_tail_estimate(p, pi0, k = 5, C_bounds = c(1e-5, Inf))
  r <- est[["r"]]; a <- est[["alpha"]]; C <- est[["C"]]
  J <- length(p)
  expect_equal(mean(p <= r), pi0 * r + (1 - pi0) * C * r^a, tolerance = 1e-10)
  expect_equal(sum(log(p[p <= r])) / J - pi0 * (r * log(r) - r) -
                 C * (1 - pi0) * r^a * (log(r) - 1 / a), 0, tolerance = 1e-10)
})

test_that("heart_test output is coherent", {
  set.seed(10)
  J <- 5000
  cfg <- sample(c("00", "01", "11"), J, TRUE, prob = c(0.9, 0.05, 0.05))
  z1 <- stats::rnorm(J, ifelse(cfg == "11", 5, 0))
  z2 <- stats::rnorm(J, ifelse(cfg != "00", 5, 0))
  p1 <- 2 * stats::pnorm(-abs(z1)); p2 <- 2 * stats::pnorm(-abs(z2))
  res <- heart_test(p1, p2, alpha = 0.1)
  r <- res$results
  expect_s3_class(res, "heart")
  expect_true(all(r$p_heart >= 0 & r$p_heart <= 1))
  expect_true(all(r$q_heart >= r$p_heart * 0 & r$q_heart <= 1))
  o <- order(r$log_T)
  expect_false(is.unsorted(r$q_heart[o]))
  expect_false(is.unsorted(r$p_heart[o]))
  expect_true(all(r$reject_fwer <= r$reject_fdr))
  # source-only signals ("01") should rarely be rejected
  fdp <- sum(r$reject_fdr & cfg != "11") / max(1, sum(r$reject_fdr))
  expect_lt(fdp, 0.2)
  expect_output(print(res), "HEART")
})

test_that("global null: no FWER discoveries and roughly uniform p_heart", {
  set.seed(11)
  J <- 4000
  res <- heart_test(stats::runif(J), stats::runif(J), alpha = 0.05)
  expect_lte(sum(res$results$reject_fwer), 1)
  expect_gt(mean(res$results$p_heart), 0.3)
})

test_that("NA p-values are excluded with a warning", {
  set.seed(12)
  p1 <- stats::runif(500); p2 <- stats::runif(500)
  p1[3] <- NA
  expect_warning(res <- heart_test(p1, p2), "excluded")
  expect_true(is.na(res$results$p_heart[3]))
  expect_equal(res$n_tests, 499)
})

test_that("heart_meta decorrelated p-value equals the source p-value", {
  m <- heart_meta(c(0.1, -0.02), c(0.05, 0.04), c(0.08, 0.03), c(0.01, 0.02))
  expect_equal(m$p_de_meta, 2 * stats::pnorm(-abs(c(0.08, 0.03) / c(0.01, 0.02))))
  expect_equal(m$rho_meta, (1 / c(0.05, 0.04)^2) / (1 / c(0.05, 0.04)^2 + 1 / c(0.01, 0.02)^2))
})

test_that("end-to-end wrapper runs", {
  sim <- heart_simulate(n_tar = 600, n_src = 3000, J = 300, seed = 5)
  res <- heart(sim$Y, sim$G, sim$D, sim$X)
  expect_s3_class(res, "heart")
  expect_equal(nrow(res$results), 300)
  expect_equal(res$aux, "pooled")
})

test_that("F01 follows (Q - pi00 F00) / (pi01 + pi11), made nondecreasing", {
  par <- make_par()
  set.seed(4)
  pd <- c(stats::runif(900), stats::rbeta(100, 0.1, 1))
  lt <- log(c(1e-3, 0.02, 0.2))
  nc <- heart_null_cdf(lt, pd, par, exact = TRUE)
  ref <- cummax(pmin(pmax((nc$Q - par$pi00 * nc$F00) / (par$pi01 + par$pi11), 0), 1))
  expect_equal(nc$F01, ref)
  expect_equal(nc$N, par$pi00 * nc$F00 + par$pi01 * nc$F01)
  expect_false(is.unsorted(nc$F01))
  # pi11 is derived from pi10 when not supplied
  par2 <- par[setdiff(names(par), "pi11")]
  expect_equal(heart_null_cdf(lt, pd, par2, exact = TRUE), nc)
})

test_that("F00 is accurate when alpha2 is at its upper bound (uninformative p_de)", {
  par <- list(pi00 = 1, pi01 = 0, pi10 = 0, pi11 = 0,
              alpha1 = 0.3, alpha2 = 1 - 1e-5, C1 = 0.5, C2 = 0.1)
  cst <- HEART:::.heart_consts(par)
  # with alpha2 -> 1, k1 q^(1 - alpha2) ~ k1 for all but astronomically small q,
  # so F00(t) ~ (t / k1)^(1 / (1 - alpha1))
  for (lt in c(-30, -5, 0)) {
    approx <- min(1, exp(cst$a * (lt - cst$logk1)))
    expect_equal(HEART:::.F00_one(lt, cst), approx, tolerance = 1e-3)
  }
  set.seed(1)
  res <- heart_test(stats::runif(3000), stats::runif(3000), par = par)
  expect_gt(min(res$results$p_heart), 1e-6)
})
