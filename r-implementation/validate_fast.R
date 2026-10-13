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

cat("=== fast vs slow vs brute, + beta check (mu=1) ===\n")
set.seed(11)
for (cfg in list(c(10,3,40,3), c(12,5,20,4), c(14,8,50,5), c(12,3,8,4))) {
  K<-cfg[1]; pf0<-cfg[2]; nk<-cfg[3]; lam<-cfg[4]
  agree_fs <- agree_bf <- betamatch <- TRUE; mg <- 0
  for (rep in 1:30) {
    d <- gen(K, pf0, nk, lam, x_corr = TRUE)
    rf <- scs_solve_fast_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lam), 1.0, 4L, 1L)
    rs <- scs_solve_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lam), 1.0, 4L, 1L)
    # brute
    combs <- combn(K, lam); best <- Inf; bestset <- NULL
    for (j in 1:ncol(combs)) { o <- prof_obj(d$Xfull, d$Y, c(1:d$pf, d$clust[combs[,j]]), 1.0)
      if (o < best) { best <- o; bestset <- combs[,j] } }
    self <- sort(rf$selected); sels <- sort(rs$selected[[1]])
    agree_fs <- agree_fs && setequal(self, sels) && abs(rf$obj - rs$obj) < 1e-9
    agree_bf <- agree_bf && setequal(self, bestset)
    # beta parity vs direct ridge solve on the chosen support
    cols <- c(1:d$pf, d$clust[self]); Xs <- d$Xfull[,cols]
    bsup <- solve(crossprod(Xs) + 1.0*diag(ncol(Xs)), crossprod(Xs, d$Y))
    bfull <- rep(0, d$pf+K); bfull[cols] <- bsup
    betamatch <- betamatch && max(abs(bfull - rf$beta)) < 1e-7
    mg <- max(mg, abs(rf$obj - best)/abs(best))
  }
  cat(sprintf("K=%d pf=%d nk=%d lam=%d: fast==slow %s | fast==brute %s | beta ok %s | maxGap %.1e\n",
              K, pf0, nk, lam, agree_fs, agree_bf, betamatch, mg))
}

cat("\n=== SPEED: fast path (mean ms over reps) ===\n")
timeit <- function(K, p_fix, nk, lambda, reps = 20) {
  set.seed(1); ts <- numeric(reps)
  for (i in 1:reps) { d <- gen(K, p_fix, nk, lambda)
    t0 <- Sys.time(); scs_solve_fast_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lambda), 1.0, 4L, 1L)
    ts[i] <- as.numeric(Sys.time()-t0, units="secs")*1000 }
  cat(sprintf("K=%-4d p_fix=%-3d nk=%-3d lam=%-3d -> %7.2f ms (median %.2f)\n",
              K, p_fix, nk, lambda, mean(ts), median(ts)))
}
timeit(4, 10, 50, 2); timeit(10, 25, 50, 4); timeit(14, 35, 50, 5)
timeit(100, 5, 50, 30, 10); timeit(170, 4, 30, 50, 10); timeit(300, 3, 30, 90, 6)
