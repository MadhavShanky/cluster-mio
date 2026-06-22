suppressMessages(library(Rcpp)); suppressMessages(library(RcppEigen))
sourceCpp("src_scs.cpp")

prof_obj <- function(Xfull, Y, cols, mu) {
  Xs <- Xfull[, cols, drop = FALSE]
  G  <- crossprod(Xs) + mu * diag(ncol(Xs))
  a  <- Y - Xs %*% solve(G, crossprod(Xs, Y))
  as.numeric(crossprod(Y, a)) / (2 * length(Y))
}
gen <- function(K, p_fix, nk, lambda, x_corr = TRUE) {
  n <- K * nk; label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)
  X <- matrix(rnorm(n * p_fix), n, p_fix)
  if (x_corr) X <- X + 0.8 * matrix(rnorm(K * p_fix), K, p_fix)[label, ]
  X <- cbind(1, X); pf <- ncol(X)
  gamma <- rep(0, K); gamma[sample(K, lambda)] <- rnorm(lambda, 0, 3)
  Y <- as.numeric(X %*% rnorm(pf) + A %*% gamma + rnorm(n))
  list(Xfull = cbind(X, A), Y = Y, pf = pf, K = K, clust = (pf + 1):(pf + K))
}

cat("=== TURBO vs brute vs fast (mu=1) ===\n")
set.seed(21)
for (cfg in list(c(10,3,40,3), c(12,5,20,4), c(14,8,50,5), c(12,3,8,4), c(16,2,30,6))) {
  K<-cfg[1]; pf0<-cfg[2]; nk<-cfg[3]; lam<-cfg[4]
  okb <- okf <- okobj <- TRUE; mg <- 0
  for (rep in 1:30) {
    d <- gen(K, pf0, nk, lam)
    rt <- scs_solve_turbo_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lam), 1.0, 4L, 1L)
    rf <- scs_solve_fast_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lam), 1.0, 4L, 1L)
    combs <- combn(K, lam); best <- Inf; bestset <- NULL
    for (j in 1:ncol(combs)) { o <- prof_obj(d$Xfull, d$Y, c(1:d$pf, d$clust[combs[,j]]), 1.0)
      if (o < best) { best <- o; bestset <- combs[,j] } }
    selt <- sort(rt$selected)
    okb <- okb && setequal(selt, bestset)
    okf <- okf && setequal(selt, sort(rf$selected))
    okobj <- okobj && abs(rt$obj - rt$obj_from_gain) < 1e-7   # gain identity self-consistency
    mg <- max(mg, abs(rt$obj - best)/abs(best))
  }
  cat(sprintf("K=%d pf=%d nk=%d lam=%d: turbo==brute %s | turbo==fast %s | gainID %s | maxGap %.1e\n",
              K, pf0, nk, lam, okb, okf, okobj, mg))
}

cat("\n=== SPEED: turbo (mean ms) ===\n")
timeit <- function(K, p_fix, nk, lambda, reps = 10, fn = scs_solve_turbo_cpp) {
  set.seed(1); ts <- numeric(reps)
  for (i in 1:reps) { d <- gen(K, p_fix, nk, lambda)
    t0 <- Sys.time(); fn(d$Xfull, d$Y, d$pf, K, as.integer(lambda), 1.0, 4L, 1L)
    ts[i] <- as.numeric(Sys.time()-t0, units="secs")*1000 }
  cat(sprintf("K=%-4d p_fix=%-3d nk=%-3d lam=%-3d -> %8.2f ms (median %.2f)\n",
              K, p_fix, nk, lambda, mean(ts), median(ts)))
}
timeit(14, 35, 50, 5, 20)
timeit(100, 5, 50, 30)
timeit(170, 4, 30, 50)
timeit(300, 3, 30, 90)
timeit(500, 3, 30, 150, 5)
timeit(1000, 3, 20, 200, 3)
