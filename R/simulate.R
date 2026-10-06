#' Simulate target / source data under the four-group model
#'
#' Generates a toy two-population data set for examples and tests. Variants
#' are independent; each variant is assigned to one of four configurations
#' with probabilities `prop`:
#' `"00"` (null in both populations), `"01"` (source-specific signal: target
#' effect zero, source effect non-zero), `"10"` (target-specific signal) and
#' `"11"` (shared signal, with source effect equal to the target effect times
#' a random factor around 1). Allele frequencies and covariate distributions
#' differ between populations; the genetic contribution is centred within each
#' population so that the pooled analysis is free of population stratification
#' (in real data, include genetic principal components in `X`). A surrogate
#' `S` that is informative of `Y` is also returned so that the doubly robust
#' estimator can be illustrated.
#'
#' @param n_tar,n_src Target and source sample sizes.
#' @param J Number of variants.
#' @param prop Probabilities of the configurations `c("00", "01", "10", "11")`.
#' @param h2_snp Variance explained by each causal variant in its population.
#' @param seed Optional random seed.
#'
#' @return A list with `Y`, `D` (1 = target), `X` (two covariates), `G`
#'   (n x J genotype matrix), `S` (surrogate) and `truth`, a data.frame with
#'   the configuration and the true target and source effects of each variant.
#' @examples
#' sim <- heart_simulate(n_tar = 500, n_src = 2000, J = 50, seed = 1)
#' table(sim$truth$config)
#' @export
heart_simulate <- function(n_tar = 2000, n_src = 20000, J = 2000,
                           prop = c(`00` = 0.90, `01` = 0.04, `10` = 0.01, `11` = 0.05),
                           h2_snp = 0.004, seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  prop <- prop / sum(prop)
  n <- n_tar + n_src
  D <- c(rep(1, n_tar), rep(0, n_src))

  maf_src <- stats::runif(J, 0.05, 0.5)
  maf_tar <- pmin(pmax(maf_src + stats::rnorm(J, 0, 0.1), 0.05), 0.5)
  G <- matrix(0, n, J, dimnames = list(NULL, paste0("snp", seq_len(J))))
  G[D == 1, ] <- matrix(stats::rbinom(n_tar * J, 2, rep(maf_tar, each = n_tar)), n_tar, J)
  G[D == 0, ] <- matrix(stats::rbinom(n_src * J, 2, rep(maf_src, each = n_src)), n_src, J)

  config <- sample(c("00", "01", "10", "11"), J, replace = TRUE, prob = prop)
  sgn <- sample(c(-1, 1), J, replace = TRUE)
  eff_tar <- sgn * sqrt(h2_snp / (2 * maf_tar * (1 - maf_tar)))
  eff_src <- sgn * sqrt(h2_snp / (2 * maf_src * (1 - maf_src)))
  beta_tar <- ifelse(config %in% c("10", "11"), eff_tar, 0)
  beta_src <- ifelse(config == "01", eff_src,
                     ifelse(config == "11", eff_tar * stats::rnorm(J, 1, 0.2), 0))

  X <- cbind(X1 = stats::rnorm(n, mean = ifelse(D == 1, 0.3, 0), sd = ifelse(D == 1, 1.2, 1)),
             X2 = stats::rbinom(n, 1, ifelse(D == 1, 0.6, 0.5)))
  S <- 0.5 * X[, "X1"] + stats::rnorm(n)
  # genetic contribution, centred within each population so that allele
  # frequency differences do not induce population stratification
  gb <- numeric(n)
  gb[D == 1] <- G[D == 1, , drop = FALSE] %*% beta_tar
  gb[D == 0] <- G[D == 0, , drop = FALSE] %*% beta_src
  gb[D == 1] <- gb[D == 1] - mean(gb[D == 1])
  gb[D == 0] <- gb[D == 0] - mean(gb[D == 0])
  Y <- gb + 0.3 * X[, "X1"] + 0.2 * X[, "X2"] + 0.6 * S + stats::rnorm(n)

  list(Y = Y, D = D, X = X, G = G, S = S,
       truth = data.frame(snp = colnames(G), config = config,
                          beta_tar = beta_tar, beta_src = beta_src,
                          maf_tar = maf_tar, maf_src = maf_src,
                          stringsAsFactors = FALSE))
}
