# Internal helpers -----------------------------------------------------------

# Inverse of a symmetric positive (semi-)definite matrix, with fallbacks.
.sym_inv <- function(A) {
  out <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (is.null(out)) {
    out <- tryCatch(solve(A), error = function(e) NULL)
  }
  if (is.null(out)) {
    out <- qr.solve(A, diag(nrow(A)))
  }
  out
}

# Build the covariate design (intercept + covariates) as a numeric matrix.
.covariate_design <- function(X, n) {
  if (is.null(X)) {
    X1 <- matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    return(X1)
  }
  if (is.data.frame(X)) {
    X <- stats::model.matrix(~ ., data = X)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be numeric or a data.frame.", call. = FALSE)
  if (nrow(X) != n) stop("`X` must have the same number of rows as `Y`.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  X1 <- cbind(`(Intercept)` = 1, X)
  qrX <- qr(X1)
  if (qrX$rank < ncol(X1)) {
    stop("The covariate design (intercept + X) is rank deficient; ",
         "remove collinear or constant columns of `X`.", call. = FALSE)
  }
  X1
}

# Coerce a genotype block to a numeric matrix and fill missing genotypes:
# "zero" sets them to 0 (as in the original analysis code); "mean" uses the
# variant mean within each population (D = 1 / D = 0).
.impute_block <- function(Gb, D, impute = "zero") {
  Gb <- as.matrix(Gb)
  storage.mode(Gb) <- "double"
  if (!anyNA(Gb)) return(Gb)
  if (impute == "zero") {
    Gb[is.na(Gb)] <- 0
    return(Gb)
  }
  for (grp in list(D == 1, D == 0)) {
    sub <- Gb[grp, , drop = FALSE]
    if (!anyNA(sub)) next
    mu <- colMeans(sub, na.rm = TRUE)
    mu[!is.finite(mu)] <- 0
    idx <- which(is.na(sub), arr.ind = TRUE)
    sub[idx] <- mu[idx[, 2L]]
    Gb[grp, ] <- sub
  }
  Gb
}

# Two-sided Wald p-value from estimate and variance.
.wald_p <- function(beta, var) {
  z <- beta / sqrt(var)
  2 * stats::pnorm(-abs(z))
}

# Numerically stable log(exp(a) + exp(b)), vectorised.
.log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

.check_D <- function(D, n) {
  if (length(D) != n) stop("`D` must have length equal to length(Y).", call. = FALSE)
  D <- as.numeric(D)
  if (anyNA(D) || !all(D %in% c(0, 1))) {
    stop("`D` must be a 0/1 vector (1 = target, 0 = source) without NA.", call. = FALSE)
  }
  if (sum(D) < 2 || sum(1 - D) < 1) {
    stop("Both target (D = 1) and source (D = 0) samples are required.", call. = FALSE)
  }
  D
}
