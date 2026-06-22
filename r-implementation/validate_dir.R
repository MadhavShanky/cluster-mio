suppressMessages(library(Rcpp)); suppressMessages(library(RcppEigen))
sourceCpp("src_scs.cpp")

inner <- function(Xfull, Y, cols, mu) {  # returns obj and fitted gamma on cluster cols
  Xs <- Xfull[, cols, drop = FALSE]
  G  <- crossprod(Xs) + mu * diag(ncol(Xs))
  coef <- solve(G, crossprod(Xs, Y))
  a <- Y - Xs %*% coef
  list(obj = as.numeric(crossprod(Y, a)) / (2 * length(Y)), coef = coef)
}
gen <- function(K, p_fix, nk, lambda) {
  n <- K * nk; label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)
  X <- cbind(1, matrix(rnorm(n * p_fix), n, p_fix) + 0.8 * matrix(rnorm(K * p_fix), K, p_fix)[label, ]); pf <- ncol(X)
  gamma <- rep(0, K); idx <- sample(K, lambda); gamma[idx] <- rnorm(lambda, 0, 3)
  Y <- as.numeric(X %*% rnorm(pf) + A %*% gamma + rnorm(n))
  list(Xfull = cbind(X, A), Y = Y, pf = pf, K = K, clust = (pf + 1):(pf + K))
}

# brute force asymmetric: min obj over supports (size 0..cap) whose fitted gamma has <=lpos pos & <=lneg neg
brute_dir <- function(d, lpos, lneg, mu, tol = 1e-9) {
  K <- d$K; cap <- min(K, lpos + lneg); best <- Inf; bp <- bn <- integer(0)
  fixed_only <- inner(d$Xfull, d$Y, 1:d$pf, mu)$obj
  cand_best <- fixed_only; bestpos <- integer(0); bestneg <- integer(0)
  for (sz in 0:cap) {
    if (sz == 0) { if (fixed_only < cand_best) cand_best <- fixed_only; next }
    cb <- combn(K, sz)
    for (j in 1:ncol(cb)) {
      cols <- c(1:d$pf, d$clust[cb[, j]])
      r <- inner(d$Xfull, d$Y, cols, mu); g <- r$coef[(d$pf + 1):length(r$coef)]
      np <- sum(g > tol); nn <- sum(g < -tol)
      if (np <= lpos && nn <= lneg && r$obj < cand_best - 1e-12) {
        cand_best <- r$obj; bestpos <- sort(cb[, j][g > tol]); bestneg <- sort(cb[, j][g < -tol])
      }
    }
  }
  list(obj = cand_best, pos = bestpos, neg = bestneg)
}

set.seed(33)
for (cfg in list(c(8,2,30,3,2,1), c(9,3,25,4,2,2), c(10,2,40,5,3,1), c(8,3,20,3,1,2))) {
  K<-cfg[1]; pf0<-cfg[2]; nk<-cfg[3]; lam<-cfg[4]; lp<-cfg[5]; ln<-cfg[6]
  okobj <- okset <- TRUE; mg <- 0
  for (rep in 1:25) {
    d <- gen(K, pf0, nk, lam)
    r <- scs_solve_dir_cpp(d$Xfull, d$Y, d$pf, K, as.integer(lp), as.integer(ln), 1.0, 5L, 1L)
    b <- brute_dir(d, lp, ln, 1.0)
    okobj <- okobj && abs(r$obj - b$obj) < 1e-7
    okset <- okset && setequal(r$selected_pos, b$pos) && setequal(r$selected_neg, b$neg)
    mg <- max(mg, abs(r$obj - b$obj)/abs(b$obj))
  }
  cat(sprintf("K=%d pf=%d lam=%d Lpos=%d Lneg=%d: obj==brute %s | set==brute %s | maxRelGap %.1e\n",
              K, pf0, lam, lp, ln, okobj, okset, mg))
}
