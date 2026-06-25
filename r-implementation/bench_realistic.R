# bench_realistic.R — turbo timing at REALISTIC budgets (flag a handful of clusters), large K.
suppressMessages(library(Rcpp)); suppressMessages(library(RcppEigen))
sourceCpp("src_scs.cpp")
gen <- function(K, p_fix, nk, lambda) {
  n <- K * nk; label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)
  X <- cbind(1, matrix(rnorm(n * p_fix), n, p_fix)); pf <- ncol(X)
  gamma <- rep(0, K); gamma[sample(K, lambda)] <- rnorm(lambda, 0, 3)
  Y <- as.numeric(X %*% rnorm(pf) + A %*% gamma + rnorm(n))
  list(Xfull = cbind(X, A), Y = Y, pf = pf, K = K)
}
timeit <- function(K, p_fix, nk, lambda, reps = 10) {
  set.seed(1); ts <- numeric(reps)
  for (i in 1:reps) { d <- gen(K, p_fix, nk, lambda)
    t0 <- Sys.time(); scs_solve_turbo_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lambda), 1.0, 4L, 1L)
    ts[i] <- as.numeric(Sys.time()-t0, units="secs")*1000 }
  cat(sprintf("K=%-4d nk=%-3d lam=%-3d -> %8.2f ms (median %.2f)\n", K, nk, lambda, mean(ts), median(ts)))
}
cat("=== REALISTIC budgets (flag <=~5% of clusters), p_fix=4 ===\n")
timeit(170, 4, 30, 10)
timeit(170, 4, 30, 20)
timeit(500, 4, 30, 15)
timeit(500, 4, 30, 25)
timeit(1000, 4, 20, 20)
timeit(1000, 4, 20, 30)
timeit(1000, 4, 20, 50)
