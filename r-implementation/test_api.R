setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
source("scs.R")
set.seed(1)
K <- 30; nk <- 40; p <- 4; n <- K * nk
cluster <- factor(rep(1:K, each = nk))
X <- matrix(rnorm(n * p), n, p)
gamma_true <- rep(0, K); active <- sample(K, 6); gamma_true[active] <- rnorm(6, 0, 3)
y <- as.numeric(cbind(1, X) %*% rnorm(p + 1) + gamma_true[cluster] + rnorm(n))

cat("==== symmetric fit (lambda=6) ====\n")
fit <- scs_fit(X, y, cluster, lambda = 6)
print(fit)
cat("true active:", sort(active), "\n")
cat("recovered  :", sort(as.integer(fit$selected)), "\n\n")

cat("==== asymmetric fit (worst=4, best=2) ====\n")
fit2 <- scs_fit(X, y, cluster, lambda = c(worst = 4, best = 2))
print(fit2)
cat("  pos:", fit2$selected_pos, " neg:", fit2$selected_neg, "\n\n")

cat("==== predict (new cluster via CART soft map) ====\n")
Xnew <- matrix(rnorm(5 * p), 5, p)
print(round(predict(fit, newdata = Xnew), 3))
cat("predict (known cluster):", round(predict(fit, Xnew[1:2,,drop=FALSE], newcluster = active[1]), 3), "\n\n")

cat("==== stability selection (B=60) ====\n")
st <- stability(fit, B = 60L, pi_thr = 0.6)
print(st)
cat("\nhead selection frequencies:\n"); print(round(head(st$freq, 10), 3))
