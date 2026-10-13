# planted_pfer_calibration.R -- design-calibrated calibration of the CPSS false-flag reference WITH signal.
# Complements null_pfer_calibration.R (all gamma_k = 0), which is the easiest case for the condition
# N subset L_{q/K}. Here the published watch-list hospitals carry their refit effects and every other
# hospital is null:
#   y* = X beta_hat + gamma_planted[cluster] + N(0, sigma_hat^2),
# with (beta_hat, gamma_planted, sigma_hat) from the unpenalized SCS refit on the real data at the
# published budget, and gamma_planted the refit effects of the stable set (B = 400, pi_thr = 0.8).
# Each replicate reruns the exact published procedure (scs_fit at lambda, stability B = 400
# complementary pairs) and records V = #flags outside the planted set, TP, q, and the reference.
# Env: STATE = ny|tx, SIM_REPS, SIM_SEED0, SIM_SHARD, SIM_B (default 400), SIM_PITHR (default 0.8).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R")
}))
state <- tolower(Sys.getenv("STATE", "ny"))
R     <- as.integer(Sys.getenv("SIM_REPS", "20"))
seed0 <- as.integer(Sys.getenv("SIM_SEED0", "5000"))
shard <- Sys.getenv("SIM_SHARD", "local")
B     <- as.integer(Sys.getenv("SIM_B", "400"))
pi_thr <- as.numeric(Sys.getenv("SIM_PITHR", "0.8"))

# ---- design: identical cleaning to null_pfer_calibration.R (published cohorts) ----
src <- readLines("sim/null_pfer_calibration.R")
i0 <- grep("^if \\(state == \"ny\"\\)", src); i1 <- grep("^\\} else stop\\(\"STATE must be ny or tx\"\\)", src)
eval(parse(text = src[i0:i1]))
y  <- log(d$los); n <- length(y); K <- nlevels(fac)
X  <- model.matrix(form, data = d)[, -1, drop = FALSE]

# ---- planted truth from the real-data fit ----
fit0 <- scs_fit(X, y, fac, lambda = lam, n_restart = 6L, seed = 1L)
st0  <- stability(fit0, B = B, pi_thr = pi_thr, seed = 1L, pairing = "complementary")
planted <- as.character(st0$flagged)
gam <- setNames(numeric(K), levels(fac)); gam[planted] <- fit0$gamma[planted]
mu_hat <- as.numeric(cbind(1, X) %*% fit0$fixef) + gam[as.character(fac)]
res0 <- y - (as.numeric(cbind(1, X) %*% fit0$fixef) + fit0$gamma[as.character(fac)])
sig_hat <- sqrt(sum(res0^2) / (n - ncol(X) - 1 - length(fit0$selected)))
cat(sprintf("[%s] N=%d K=%d lambda=%d B=%d | planted %d hospitals (|gamma| %.3f-%.3f) sigma_hat=%.3f\n",
            toupper(state), n, K, lam, B, length(planted), min(abs(gam[planted])), max(abs(gam[planted])), sig_hat))

Vv <- TPv <- qv <- refv <- numeric(R)
for (r in seq_len(R)) {
  set.seed(seed0 + r)
  ys  <- mu_hat + rnorm(n, 0, sig_hat)
  fit <- scs_fit(X, ys, fac, lambda = lam, n_restart = 4L, seed = 1L)
  st  <- stability(fit, B = B, pi_thr = pi_thr, seed = r, pairing = "complementary")
  fl  <- as.character(st$flagged)
  Vv[r] <- length(setdiff(fl, planted)); TPv[r] <- length(intersect(fl, planted))
  qv[r] <- st$q; refv[r] <- st$expected_false_flags_bound
  cat(sprintf("  rep %d/%d: V=%d TP=%d q=%.2f ref=%.2f\n", r, R, Vv[r], TPv[r], qv[r], refv[r]))
}
dir.create(file.path("data", state), showWarnings = FALSE, recursive = TRUE)
outf <- sprintf("data/%s/planted_pfer_%s.csv", state, shard)
write.csv(data.frame(state = state, rep = seed0 + seq_len(R), V = Vv, TP = TPv, n_planted = length(planted),
                     q = qv, ref = refv), outf, row.names = FALSE)
cat(sprintf("mean V = %.3f (ref mean %.3f); mean TP = %.2f of %d; saved %s\n",
            mean(Vv), mean(refv), mean(TPv), length(planted), outf))
