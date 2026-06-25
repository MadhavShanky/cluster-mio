# test_deliverables.R — exercises D1 (relative mu + unshrunk refit), C (audit_budget), D (cluster_ranking).
setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
source("scs.R")
set.seed(1)
K <- 30; nk <- 40; p <- 4; n <- K * nk
cluster <- factor(rep(1:K, each = nk))
X <- matrix(rnorm(n * p), n, p)
gamma_true <- rep(0, K); active <- sample(K, 6); gamma_true[active] <- rnorm(6, 0, 3)
y <- as.numeric(cbind(1, X) %*% rnorm(p + 1) + gamma_true[cluster] + rnorm(n))

cat("==== D1: relative-mu default + unshrunk refit ====\n")
fit <- scs_fit(X, y, cluster, lambda = 6)
cat(sprintf("relative mu = %.4g (eps=%g)\n", fit$mu, fit$eps))
cat("true active:", sort(active), "\n")
cat("recovered  :", sort(as.integer(fit$selected)), "\n")
# unshrunk check: refit gamma on selected support == plain lm on [1|X|A_sel]
sel_i <- as.integer(fit$selected)
A <- model.matrix(~ cluster - 1)
lmfit <- lm(y ~ X + A[, sel_i] - 1 + 1)  # intercept + X + selected dummies
cat(sprintf("max|refit gamma - lm gamma| = %.2e (should be ~0 -> 'unshrunk' literally holds)\n",
            max(abs(fit$gamma[fit$selected] - tail(coef(lmfit), length(sel_i))))))

cat("\n==== D1: epsilon-sensitivity ====\n")
print(scs_eps_sensitivity(X, y, cluster, lambda = 6))

cat("\n==== C: audit budget (value of each review slot) ====\n")
pth <- scs_path(X, y, cluster, lambda_max = 12)
ab <- audit_budget(pth, target = 0.90)
print(ab)
stopifnot(all(diff(ab$table$detected) >= -1e-8))  # detected deviation nondecreasing in budget

cat("\n==== D: lambda-path cluster ranking ====\n")
rk <- cluster_ranking(pth)
print(rk)
cat("top-6 ranked vs true active:\n")
cat("  ranked:", head(rk$ranking$cluster, 6), "\n")
cat("  active:", sort(active), "\n")

cat("\nALL DELIVERABLE CHECKS RAN\n")
