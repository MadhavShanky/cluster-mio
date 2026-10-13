# null_pfer_calibration.R -- C2 (referee round 3): design-calibrated NULL calibration of the
# stability watch-list's false-flag count under the EXACT implemented procedure (observation-level
# subsampling). Tests empirically whether the reported Meinshausen-Buhlmann PFER value tracks the
# realized expected number of false flags when all cluster effects are zero and clusters are unequal /
# non-exchangeable -- the validity question Codex flagged (obs-level subsampling does not satisfy the MB
# exchangeability assumption, so the bound is not theorem-backed; this provides simulation-calibrated
# evidence instead).
#
# DGP (null): Y0 = X beta_hat + N(0, sigma_hat^2), all gamma_k = 0, with X, cluster, n_k taken from the
# REAL design (SPARCS NY or PUDF TX). beta_hat, sigma_hat from the fixed-effects-only (null) OLS fit.
# Procedure: for each replicate, run stability(scs_fit(X, Y0, fac, lambda), B, frac=0.5, pi_thr) exactly
# as in the paper; record V = #flagged (all are FALSE flags under the null), q, and which clusters flag.
#
# Usage (local):   Rscript sim/null_pfer_calibration.R           # STATE=ny default, 200 reps
#   env STATE=ny|tx ; SIM_NULL_REPS=<R> ; SIM_NULL_SEED0=<int> ; SIM_NULL_SHARD=<tag for output file>
# Cluster rep-sharding: launch array tasks with disjoint SIM_NULL_SEED0 and a per-task SIM_NULL_SHARD,
#   then pool the per-rep CSVs (no MC loss). Mirrors sim/run_sims.R rep-shard pattern.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R")
}))

state <- tolower(Sys.getenv("STATE", "ny"))
R     <- as.integer(Sys.getenv("SIM_NULL_REPS", "200"))
seed0 <- as.integer(Sys.getenv("SIM_NULL_SEED0", "1000"))
shard <- Sys.getenv("SIM_NULL_SHARD", "local")
B     <- as.integer(Sys.getenv("SIM_NULL_B", "100"))
pi_thr <- as.numeric(Sys.getenv("SIM_NULL_PITHR", "0.8"))

# ---- load the real design (full case-mix spec), matching realdata_sensitivity.R / tx_sensitivity.R ----
if (state == "ny") {
  d <- read.csv("data/sparcs/drg194_2024.csv", stringsAsFactors = FALSE)
  d$los <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$length_of_stay)))
  sev_lv <- c("Minor","Moderate","Major","Extreme")
  d$sev <- factor(d$apr_severity_of_illness, levels = sev_lv)
  d$rom <- factor(d$apr_risk_of_mortality,  levels = sev_lv)
  d$age <- factor(d$age_group); d$sex <- factor(d$gender); d$adm <- factor(d$type_of_admission)
  d$fac <- d$permanent_facility_id
  keep <- !is.na(d$los) & d$los > 0 & d$sex %in% c("M","F") &
          complete.cases(d[, c("sev","rom","age","sex","adm","fac")])
  d <- d[keep, ]; tb <- table(d$fac); d <- d[d$fac %in% names(tb)[tb >= 25L], ]
  fac <- droplevels(factor(d$fac)); for (v in c("age","sex","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
  form <- ~ age + sex + sev + rom + adm
  lam  <- as.integer(Sys.getenv("SIM_NULL_LAMBDA", "21"))   # NY watch-list budget (Supplement S4)
} else if (state == "tx") {
  # cleaning MUST mirror tx_hf_analysis.R / tx_sensitivity.R exactly so the null cohort matches the
  # published K=208 analysis (sev/rom validity via complete.cases + MIN_NK>=25 low-volume cutoff).
  MIN_NK <- 25L
  d <- read.csv("data/tx/tx_hf_1q2019.csv", stringsAsFactors = FALSE,
                colClasses = c(thcic_id = "character", age = "character",
                               rom = "character", sev = "character", adm = "character"))
  d$sex <- trimws(d$sex)
  d <- d[d$sex %in% c("M","F") & !is.na(d$los) & d$los > 0, ]
  d$sev <- factor(d$sev, levels = c("1","2","3","4"))        # 0/blank ("no class") -> NA, dropped
  d$rom <- factor(d$rom, levels = c("1","2","3","4"))
  d$age <- factor(d$age); d$sexf <- factor(d$sex); d$adm <- factor(d$adm)
  d <- d[complete.cases(d[, c("sev","rom","age","adm","sexf")]), ]
  tb <- table(d$thcic_id); d <- d[d$thcic_id %in% names(tb)[tb >= MIN_NK], ]
  fac <- droplevels(factor(d$thcic_id))
  for (v in c("age","sexf","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
  form <- ~ age + sexf + sev + rom + adm
  lam  <- as.integer(Sys.getenv("SIM_NULL_LAMBDA", "20"))   # TX watch-list budget (Supplement S4)
} else stop("STATE must be ny or tx")

y  <- log(d$los); n <- length(y); K <- nlevels(fac)
X  <- model.matrix(form, data = d)[, -1, drop = FALSE]
nk <- as.integer(table(fac))

# ---- null fit: fixed effects only (no cluster effect) -> beta_hat, sigma_hat ----
fe  <- lm.fit(cbind(1, X), y)
mu_hat  <- as.numeric(cbind(1, X) %*% fe$coefficients)
sig_hat <- sqrt(sum(fe$residuals^2) / (n - length(fe$coefficients)))
nominal_mb <- function(q) q^2 / ((2 * pi_thr - 1) * K)   # MB formula the paper reports
cat(sprintf("[%s] N=%d K=%d lambda=%d B=%d pi_thr=%.2f sigma_hat=%.3f | NULL DGP (all gamma_k=0)\n",
            toupper(state), n, K, lam, B, pi_thr, sig_hat))

# ---- replicate loop ----
Vv <- numeric(R); qv <- numeric(R); boundv <- numeric(R)
flagcount <- setNames(numeric(K), levels(fac))
for (r in 1:R) {
  set.seed(seed0 + r)
  y0  <- mu_hat + rnorm(n, 0, sig_hat)
  fit <- scs_fit(X, y0, fac, lambda = lam, n_restart = 4L, seed = 1L)
  st  <- stability(fit, B = B, frac = 0.5, pi_thr = pi_thr, seed = r, pairing = "complementary")
  Vv[r]     <- length(st$flagged)
  qv[r]     <- st$q
  boundv[r] <- st$expected_false_flags_bound
  if (length(st$flagged)) flagcount[st$flagged] <- flagcount[st$flagged] + 1
  if (r %% 25 == 0) cat(sprintf("  rep %d/%d : meanV=%.3f meanQ=%.2f nomMB=%.2f\n",
                                r, R, mean(Vv[1:r]), mean(qv[1:r]), mean(boundv[1:r])))
}

# ---- per-rep output (poolable across shards) ----
dir.create(file.path("data", state), showWarnings = FALSE, recursive = TRUE)
perrep <- data.frame(state = state, rep = seed0 + (1:R), V = Vv, q = qv, mb_nominal = boundv)
outf <- sprintf("data/%s/null_pfer_%s.csv", state, shard)
write.csv(perrep, outf, row.names = FALSE)

# ---- per-cluster false-flag frequency vs size (does false flagging concentrate in large clusters?) ----
sz_terc <- cut(nk, quantile(nk, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("small","mid","large"))
byterc  <- tapply(flagcount / R, sz_terc, mean)
write.csv(data.frame(cluster = levels(fac), nk = nk, false_flag_freq = flagcount / R,
                     size_tercile = sz_terc),
          sprintf("data/%s/null_pfer_bycluster_%s.csv", state, shard), row.names = FALSE)

# ---- summary ----
cat(sprintf("\n==== NULL PFER CALIBRATION [%s], R=%d reps ====\n", toupper(state), R))
cat(sprintf("realized E[V] (mean false flags)  = %.3f  (MC SE %.3f)\n", mean(Vv), sd(Vv)/sqrt(R)))
cat(sprintf("  V quantiles  50/90/95/max        = %.0f / %.0f / %.0f / %.0f\n",
            quantile(Vv,.5), quantile(Vv,.9), quantile(Vv,.95), max(Vv)))
cat(sprintf("mean q (avg selected)              = %.3f\n", mean(qv)))
cat(sprintf("nominal MB formula value q^2/((2pi-1)K) = %.3f\n", mean(boundv)))
cat(sprintf("P(realized V > nominal MB value)   = %.3f\n", mean(Vv > boundv)))
cat(sprintf("false-flag freq by size tercile    : small=%.3f mid=%.3f large=%.3f\n",
            byterc["small"], byterc["mid"], byterc["large"]))
cat(sprintf("saved %s\n", outf))
