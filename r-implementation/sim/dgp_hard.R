# dgp_hard.R -- harder-regime DGPs: regimes that do NOT favor the SCS
# point-mass estimand. Mirrors gen_scs() in dgp.R (same x_{ik}=m_k+e_{ik} design, same eps), adding:
#   null_sd  > 0 : the NON-active "null" clusters carry small nonzero gamma ~ N(0, null_sd^2) (so "true
#                  zero" is not clean -- tests spurious flagging of near-zero clusters).
#   regime "dense_t": ALL gamma ~ scaled t_nu, dense (no zeros), scaled to target ICC -- heavy-tailed
#                  dense effects (L0 target misspecified; misspecification stress test beyond the Gaussian regime).
#   size_effect TRUE: deviating clusters are systematically SMALLER (n_k correlated with |gamma_k|),
#                  an unfavorable size/effect coupling.
# TRUE active set for selection metrics = the +/- delta clusters only (null_sd clusters are NOT active).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("sim/dgp.R")   # for sel_scores() and the reference design
}))

gen_hard <- function(K, nk, p, regime = c("sparse", "dense_t"),
                     pi_active = 0.1, delta = 1, icc = 0.3, rho_x = 0.5,
                     sigma = 1, err = c("normal", "t"), nu = 3,
                     null_sd = 0, size_effect = FALSE, seed = NULL) {
  regime <- match.arg(regime); err <- match.arg(err)
  if (!is.null(seed)) set.seed(seed)

  gamma <- numeric(K); active <- integer(0)
  if (regime == "sparse") {
    Ka <- max(1L, round(pi_active * K))
    active <- sort(sample.int(K, Ka))
    gamma[active] <- delta * sigma * sample(c(-1, 1), Ka, replace = TRUE)
    if (null_sd > 0) {                                   # small nonzero "nulls" (true zero not clean)
      nulls <- setdiff(seq_len(K), active)
      gamma[nulls] <- rnorm(length(nulls), 0, null_sd * sigma)
    }
  } else {                                               # dense_t: heavy-tailed dense effects
    tau2 <- icc / (1 - icc) * sigma^2
    g <- rt(K, nu) * sqrt((nu - 2) / nu)                 # var 1
    gamma <- g * sqrt(tau2)
    active <- seq_len(K)                                 # all nonzero
  }

  # cluster sizes: optionally couple size to |gamma| (deviating clusters smaller -> harder to detect)
  if (size_effect) {
    base <- pmax(5L, as.integer(round(rgamma(K, shape = 1.5, scale = nk / 1.5))))
    ord  <- order(abs(gamma), decreasing = TRUE)         # largest |gamma| -> smallest n_k
    nk_vec <- integer(K); nk_vec[ord] <- sort(base)      # ascending sizes assigned to descending |gamma|
  } else {
    nk_vec <- rep(as.integer(nk), K)
  }
  n  <- sum(nk_vec); cl <- rep(seq_len(K), times = nk_vec)
  mk <- matrix(rnorm(K * p, 0, sqrt(rho_x)),     K, p)
  ew <- matrix(rnorm(n * p, 0, sqrt(1 - rho_x)), n, p)
  X  <- mk[cl, , drop = FALSE] + ew
  beta <- rnorm(p)
  eps <- if (err == "t") rt(n, nu) * sigma * sqrt((nu - 2) / nu) else rnorm(n, 0, sigma)
  y <- as.numeric(X %*% beta) + gamma[cl] + eps
  list(X = X, y = y, cluster = factor(cl), beta = beta, gamma = gamma, active = active,
       regime = regime, K = K, nk = nk, nk_vec = nk_vec, p = p,
       pi_active = pi_active, delta = delta, icc = icc, rho_x = rho_x, sigma = sigma,
       err = err, nu = nu, null_sd = null_sd, size_effect = size_effect,
       # mk_full is the between-cluster covariate mean -- the driver drops a column of X to induce
       # omitted-variable (misspecified-FE) confounding that is cluster-correlated through mk.
       mk = mk)
}
