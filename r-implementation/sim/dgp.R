# dgp.R — data-generating process for the SCS simulations.
#   y_{ik} = x_{ik}'beta + gamma_k + eps_{ik},  eps ~ N(0, sigma^2)
#   x_{ik} = m_k + e_{ik},  m_k ~ N(0, rho_x I_p) [between],  e_{ik} ~ N(0,(1-rho_x) I_p) [within]
#   (S) sparse:   round(pi*K) clusters deviate with gamma_k = +/- delta*sigma; rest = 0
#   (G) gaussian: gamma_k ~ N(0, tau^2) with tau^2 = icc/(1-icc)*sigma^2  (=> ICC = tau^2/(tau^2+sigma^2))

# Regimes:
#  sparse   : gamma_k = +/- delta on round(pi*K) random clusters (gamma _|_ X) -- selection/estimation.
#  sparse_x : gamma_k linked to a PREDICTION-ONLY cluster covariate Z (NOT in the estimation X, so the LMM
#             stays exogenous, E[gamma|X]=0). active = clusters with the most extreme standardized Z-score;
#             gamma graded (a'Z), RMS scaled to delta (NOT sign(.) -- avoids being trivially CART-aligned).
#             predict.scs / a Z->gamma_hat rule can then recover new-cluster effects (deliverable D gate).
#  gaussian : gamma_k ~ N(0, tau^2), tau^2 = icc/(1-icc)*sigma^2 (CONDITIONAL random-intercept ICC) --
#             L0 target MIS-specified here -> misspecification check (L0 should not dominate; the L0 target is misspecified here).
# Robustness knobs: err in {normal,t(nu)} (t scaled to var sigma^2); unequal_nk -> n_k ~ Gamma(shape 1.5).
gen_scs <- function(K, nk, p, regime = c("sparse", "sparse_x", "gaussian"),
                    pi_active = 0.2, delta = 3, icc = 0.3, rho_x = 0.5,
                    sigma = 1, q = 3, err = c("normal", "t"), nu = 3,
                    unequal_nk = FALSE, seed = NULL) {
  regime <- match.arg(regime); err <- match.arg(err)
  if (!is.null(seed)) set.seed(seed)
  nk_vec <- if (unequal_nk) pmax(5L, as.integer(round(rgamma(K, shape = 1.5, scale = nk / 1.5)))) else rep(as.integer(nk), K)
  n  <- sum(nk_vec)
  cl <- rep(seq_len(K), times = nk_vec)
  mk <- matrix(rnorm(K * p, 0, sqrt(rho_x)),       K, p)   # between-cluster covariate means
  ew <- matrix(rnorm(n * p, 0, sqrt(1 - rho_x)),   n, p)   # within-cluster covariate noise
  X  <- mk[cl, , drop = FALSE] + ew
  beta <- rnorm(p)
  gamma <- numeric(K); active <- integer(0); Z <- NULL
  if (regime == "sparse") {
    Ka <- max(1L, round(pi_active * K))
    active <- sort(sample.int(K, Ka))
    gamma[active] <- delta * sigma * sample(c(-1, 1), Ka, replace = TRUE)
  } else if (regime == "sparse_x") {
    Z <- matrix(rnorm(K * q), K, q)                          # cluster-level PREDICTION-ONLY covariates
    a <- rnorm(q); s <- as.numeric(scale(Z %*% a))           # standardized linear score
    Ka <- max(1L, round(pi_active * K))
    active <- sort(order(abs(s), decreasing = TRUE)[seq_len(Ka)])
    gamma[active] <- delta * sigma * s[active] / sqrt(mean(s[active]^2))  # graded, Z-linked, RMS = delta
  } else {
    tau2  <- icc / (1 - icc) * sigma^2
    gamma <- rnorm(K, 0, sqrt(tau2))
    active <- seq_len(K)              # no exact zeros under the Gaussian regime
  }
  eps <- if (err == "t") rt(n, nu) * sigma * sqrt((nu - 2) / nu) else rnorm(n, 0, sigma)
  y <- as.numeric(X %*% beta) + gamma[cl] + eps
  list(X = X, y = y, cluster = factor(cl), beta = beta, gamma = gamma, active = active, Z = Z,
       regime = regime, K = K, nk = nk, nk_vec = nk_vec, p = p, q = q,
       pi_active = pi_active, delta = delta, icc = icc, rho_x = rho_x, sigma = sigma,
       err = err, nu = nu, unequal_nk = unequal_nk)
}

# selection scores of a recovered active set vs truth (regime "sparse")
sel_scores <- function(selected, active, K) {
  sel <- as.integer(selected); act <- as.integer(active)
  tp <- length(intersect(sel, act)); fp <- length(setdiff(sel, act)); fn <- length(setdiff(act, sel))
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  rec  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  fdr  <- if (tp + fp > 0) fp / (tp + fp) else 0
  list(precision = prec, recall = rec, f1 = f1, fdr = fdr, tp = tp, fp = fp, fn = fn)
}
