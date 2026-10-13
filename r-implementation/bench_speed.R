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
timeit <- function(K, p_fix, nk, lambda, reps = 20) {
  set.seed(1); ts <- numeric(reps)
  for (i in 1:reps) {
    d <- gen(K, p_fix, nk, lambda)
    t0 <- Sys.time()
    scs_solve_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lambda), 1.0, n_restart = 4L, seed = 1L)
    ts[i] <- as.numeric(Sys.time() - t0, units = "secs") * 1000
  }
  cat(sprintf("K=%-4d p_fix=%-3d nk=%-3d lam=%-3d  -> %7.1f ms (median %.1f)\n",
              K, p_fix, nk, lambda, mean(ts), median(ts)))
}

cat("== paper sim settings (single fit; paper Table reports ms) ==\n")
timeit(4,  10, 50, 2)    # Low
timeit(10, 25, 50, 4)    # Medium
timeit(14, 35, 50, 5)    # High
cat("== hospital scale (K >> sim) ==\n")
timeit(100, 5, 50, 30, reps = 10)
timeit(170, 3, 30, 50, reps = 10)
timeit(300, 3, 30, 90, reps = 6)
