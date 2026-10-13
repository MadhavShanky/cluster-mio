suppressMessages(library(Rcpp)); suppressMessages(library(RcppEigen))
sourceCpp("src_scs.cpp")

prof_obj <- function(Xfull, Y, cols, mu) {
  Xs <- Xfull[, cols, drop = FALSE]
  G  <- crossprod(Xs) + mu * diag(ncol(Xs))
  a  <- Y - Xs %*% solve(G, crossprod(Xs, Y))
  as.numeric(crossprod(Y, a)) / (2 * length(Y))
}

gen <- function(K, p_fix, nk, lambda, x_corr) {
  n <- K * nk; label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)
  X <- matrix(rnorm(n * p_fix), n, p_fix)
  if (x_corr) X <- X + 0.8 * matrix(rnorm(K * p_fix), K, p_fix)[label, ]
  X <- cbind(1, X); pf <- ncol(X)
  gamma <- rep(0, K); active <- sample(K, lambda); gamma[active] <- rnorm(lambda, 0, 3)
  Y <- as.numeric(X %*% rnorm(pf) + A %*% gamma + rnorm(n))
  list(Xfull = cbind(X, A), Y = Y, pf = pf, K = K, clust = (pf + 1):(pf + K))
}

one <- function(K, p_fix, nk, lambda, mu, x_corr) {
  d <- gen(K, p_fix, nk, lambda, x_corr)
  # brute force ground truth
  combs <- combn(K, lambda); best <- Inf; bestset <- NULL
  for (j in 1:ncol(combs)) {
    o <- prof_obj(d$Xfull, d$Y, c(1:d$pf, d$clust[combs[, j]]), mu)
    if (o < best) { best <- o; bestset <- combs[, j] }
  }
  r <- scs_solve_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lambda), mu, n_restart = 4L, seed = 1L)
  sel <- sort(r$selected[[1]])
  list(match = setequal(bestset, sel), gap = (r$obj - best) / abs(best),
       objchk = abs(r$obj - prof_obj(d$Xfull, d$Y, c(1:d$pf, d$clust[sel]), mu)))
}

bench <- function(K, p_fix, nk, lambda, mu, x_corr, reps = 50) {
  set.seed(100)
  r <- replicate(reps, one(K, p_fix, nk, lambda, mu, x_corr), simplify = FALSE)
  cat(sprintf("K=%d pf=%d nk=%d lam=%d mu=%g corr=%-5s -> match %.2f  maxGap %.2e  objErr %.1e\n",
              K, p_fix, nk, lambda, mu, x_corr,
              mean(sapply(r, `[[`, "match")), max(sapply(r, `[[`, "gap")), max(sapply(r, `[[`, "objchk"))))
}

bench(10, 3, 40, 3, 1.0, FALSE)
bench(10, 3, 40, 3, 1.0, TRUE)
bench(12, 5, 20, 4, 1.0, TRUE)
bench(14, 8, 50, 5, 1.0, TRUE)
bench(12, 3, 8,  4, 0.01, TRUE)
