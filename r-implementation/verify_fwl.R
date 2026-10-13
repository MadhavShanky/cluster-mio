# Verify whether Q=0 cardinality-constrained random-intercept selection (the cluster-mio
# inner/outer problem) reduces to an EXACT one-shot ranking of clusters (architecture C),
# vs. requiring true combinatorial search. Brute force over all size-lambda subsets is ground truth.

set.seed(1)

# profiled objective c(S) exactly as inner_opt computes it:
#   Xs = [X | A_S],  alpha = Y - Xs (muI + Xs'Xs)^{-1} Xs'Y,  obj = <Y,alpha>/(2n)
prof_obj <- function(Xfull, Y, cols, mu) {
  Xs <- Xfull[, cols, drop = FALSE]
  G  <- t(Xs) %*% Xs + mu * diag(ncol(Xs))
  a  <- Y - Xs %*% solve(G, t(Xs) %*% Y)
  as.numeric(crossprod(Y, a)) / (2 * length(Y))
}

run_one <- function(K, p_fix, nk, lambda, mu, x_corr_cluster) {
  n <- K * nk
  label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)            # n x K one-hot
  # fixed design: optionally correlate X with cluster identity (breaks X ⟂ clusters)
  X <- matrix(rnorm(n * p_fix), n, p_fix)
  if (x_corr_cluster) {
    cluster_shift <- matrix(rnorm(K * p_fix), K, p_fix)[label, ]
    X <- X + 0.8 * cluster_shift                     # X now correlated with cluster
  }
  X <- cbind(1, X)                                   # intercept
  pf <- ncol(X)
  gamma <- rep(0, K); active <- sample(K, lambda); gamma[active] <- rnorm(lambda, 0, 3)
  Y <- X %*% rnorm(pf) + A %*% gamma + rnorm(n, 0, 1)
  Xfull <- cbind(X, A)
  fixed_cols <- 1:pf
  clust_cols <- (pf + 1):(pf + K)

  # ground truth: brute force best subset of size lambda
  combs <- combn(K, lambda)
  best <- Inf; best_set <- NULL
  for (j in 1:ncol(combs)) {
    cols <- c(fixed_cols, clust_cols[combs[, j]])
    o <- prof_obj(Xfull, Y, cols, mu)
    if (o < best) { best <- o; best_set <- combs[, j] }
  }

  # architecture C candidate: rank clusters by marginal gradient at the fixed-only fit
  s_fix <- rep(0, pf + K); s_fix[fixed_cols] <- 1
  Xs <- Xfull[, fixed_cols, drop = FALSE]
  G  <- t(Xs) %*% Xs + mu * diag(pf)
  a  <- Y - Xs %*% solve(G, t(Xs) %*% Y)            # residual after fixed-only fit
  score <- ( as.numeric(t(Xfull[, clust_cols]) %*% a) )^2   # (a_k' alpha)^2 ~ -grad
  rank_set <- order(score, decreasing = TRUE)[1:lambda]
  rank_obj <- prof_obj(Xfull, Y, c(fixed_cols, clust_cols[rank_set]), mu)

  list(bf_obj = best, rank_obj = rank_obj,
       match = setequal(best_set, rank_set),
       gap = (rank_obj - best) / abs(best))
}

cat("=== X ORTHOGONAL to clusters (x_corr_cluster=FALSE) ===\n")
res1 <- replicate(40, run_one(10, 3, 40, 3, 1.0, FALSE), simplify = FALSE)
cat(sprintf("exact-match rate: %.2f   mean rel gap: %.2e   max rel gap: %.2e\n",
            mean(sapply(res1, `[[`, "match")),
            mean(sapply(res1, `[[`, "gap")), max(sapply(res1, `[[`, "gap"))))

cat("=== X CORRELATED with clusters (x_corr_cluster=TRUE) ===\n")
res2 <- replicate(40, run_one(10, 3, 40, 3, 1.0, TRUE), simplify = FALSE)
cat(sprintf("exact-match rate: %.2f   mean rel gap: %.2e   max rel gap: %.2e\n",
            mean(sapply(res2, `[[`, "match")),
            mean(sapply(res2, `[[`, "gap")), max(sapply(res2, `[[`, "gap"))))

cat("=== smaller nk=8 (less info per cluster, stronger coupling) ===\n")
res3 <- replicate(40, run_one(10, 3, 8, 3, 1.0, TRUE), simplify = FALSE)
cat(sprintf("exact-match rate: %.2f   mean rel gap: %.2e   max rel gap: %.2e\n",
            mean(sapply(res3, `[[`, "match")),
            mean(sapply(res3, `[[`, "gap")), max(sapply(res3, `[[`, "gap"))))
