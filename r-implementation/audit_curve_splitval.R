# audit_curve_splitval.R -- quantify the in-sample optimism of the audit-budget capture curve (ref2 item B).
# The figure R3 / audit_budget() curve is computed in sample: detected(lambda)=obj(0)-obj(lambda), where the
# support at each lambda is chosen to MINIMIZE that same in-sample objective. The manuscript now states this
# is optimistic and that "a sample-split or cross-validated version removes the in-sample optimism". This
# script produces that held-out curve and the gap, for NY SPARCS and TX PUDF.
#
# Method: repeated 50/50 observation splits. On the SELECT half run scs_path -> support S_lambda. On the
# EVAL half, for each lambda, refit y ~ [1 : X : D_{S_lambda}] by unpenalized LS and record the held-out
# residual SS; detected_eval(lambda)=RSS_eval(0)-RSS_eval(lambda), normalized by its value at lambda_max.
# Compare lambda*(90%) in sample vs held out. Held-out needs MORE units -> demonstrates the optimism.
# Outputs data/audit_splitval.csv (per state x lambda: frac_insample, frac_heldout) + a printed summary.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R")
}))
RNGkind("Mersenne-Twister")
LMAX   <- 25L
NSPLIT <- as.integer(Sys.getenv("SPLITVAL_NSPLIT", "20"))
TARGET <- 0.90

## ---- state loaders: reproduce the exact analysis cohorts -----------------------------------------
load_ny <- function() {
  d <- read.csv("data/sparcs/drg194_2024.csv", stringsAsFactors = FALSE)
  los <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$length_of_stay)))
  lv <- c("Minor","Moderate","Major","Extreme")
  d$sev <- factor(d$apr_severity_of_illness, levels = lv); d$rom <- factor(d$apr_risk_of_mortality, levels = lv)
  d$age <- factor(d$age_group); d$sex <- factor(d$gender); d$adm <- factor(d$type_of_admission)
  d$fac <- d$permanent_facility_id; d$los <- los
  keep <- !is.na(d$los) & d$los > 0 & d$sex %in% c("M","F") &
          complete.cases(d[, c("sev","rom","age","sex","adm","fac")])
  d <- d[keep, ]; tb <- table(d$fac); d <- d[d$fac %in% names(tb)[tb >= 25L], ]
  fac <- droplevels(factor(d$fac)); for (v in c("age","sex","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
  list(X = model.matrix(~ age + sex + sev + rom + adm, data = d)[, -1, drop = FALSE], y = log(d$los), fac = fac)
}
load_tx <- function() {
  d <- read.csv("data/tx/tx_hf_1q2019.csv", stringsAsFactors = FALSE,
                colClasses = c(thcic_id = "character", age = "character",
                               rom = "character", sev = "character", adm = "character"))
  d$sex <- trimws(d$sex); d <- d[d$sex %in% c("M","F") & !is.na(d$los) & d$los > 0, ]
  d$sev <- factor(d$sev, levels = c("1","2","3","4")); d$rom <- factor(d$rom, levels = c("1","2","3","4"))
  d$age <- factor(d$age); d$sexf <- factor(d$sex); d$adm <- factor(d$adm)
  d <- d[complete.cases(d[, c("sev","rom","age","adm","sexf")]), ]
  tb <- table(d$thcic_id); d <- d[d$thcic_id %in% names(tb)[tb >= 25L], ]
  fac <- droplevels(factor(d$thcic_id)); for (v in c("age","sexf","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
  list(X = model.matrix(~ age + sexf + sev + rom + adm, data = d)[, -1, drop = FALSE], y = log(d$los), fac = fac)
}

# unpenalized held-out residual SS for a fixed support S (cluster-level names) on rows `idx`
rss_support <- function(X, y, fac, idx, S) {
  Xe <- cbind(1, X[idx, , drop = FALSE]); ye <- y[idx]; fe <- droplevels(fac[idx])
  if (length(S)) {
    Suse <- intersect(S, levels(fe))
    if (length(Suse)) { D <- model.matrix(~ fe - 1); colnames(D) <- levels(fe); Xe <- cbind(Xe, D[, Suse, drop = FALSE]) }
  }
  r <- stats::lm.fit(Xe, ye)$residuals
  sum(r^2)
}

# first lambda whose fraction-of-max >= target
lam_star <- function(lambda, frac, target) { hit <- which(frac >= target); if (length(hit)) lambda[min(hit)] else max(lambda) }

run_state <- function(name, dat) {
  X <- dat$X; y <- dat$y; fac <- dat$fac; n <- length(y); lev <- levels(fac)
  lmax <- min(LMAX, nlevels(fac) - 1L); budgets <- 0:lmax
  ## in-sample curve (matches audit_budget on the full data)
  pin <- scs_path(X, y, fac, lambda_max = lmax, n_restart = 4L, seed = 1L)
  det_in <- pin$obj[1] - pin$obj; frac_in <- if (tail(det_in,1) > 0) det_in / tail(det_in,1) else det_in * 0
  ## held-out curve, averaged over NSPLIT 50/50 observation splits
  acc <- matrix(0, nrow = length(budgets), ncol = NSPLIT)
  set.seed(1L)
  for (s in seq_len(NSPLIT)) {
    sel <- sample.int(n, floor(n / 2)); ev <- setdiff(seq_len(n), sel)
    csel <- droplevels(fac[sel])
    ptr <- tryCatch(scs_path(X[sel, , drop = FALSE], y[sel], csel,
                             lambda_max = min(lmax, nlevels(csel) - 1L), n_restart = 4L, seed = s),
                    error = function(e) NULL)
    if (is.null(ptr)) { acc[, s] <- NA; next }
    rss0 <- rss_support(X, y, fac, ev, character(0))
    for (i in seq_along(budgets)) {
      b <- budgets[i]; col <- paste0("l", b)
      S <- if (col %in% colnames(ptr$member)) rownames(ptr$member)[ptr$member[, col] == 1L] else character(0)
      acc[i, s] <- rss0 - rss_support(X, y, fac, ev, S)        # detected on held-out
    }
  }
  det_ev <- rowMeans(acc, na.rm = TRUE); det_ev[det_ev < 0] <- 0
  frac_ev <- if (max(det_ev) > 0) det_ev / max(det_ev) else det_ev
  ls_in <- lam_star(budgets, frac_in, TARGET); ls_ev <- lam_star(budgets, frac_ev, TARGET)
  cat(sprintf("%-9s  lambda*(90%%) in-sample = %d   held-out = %d   (optimism gap = %d units)\n",
              name, ls_in, ls_ev, ls_ev - ls_in))
  data.frame(state = name, lambda = budgets, frac_insample = round(frac_in, 4), frac_heldout = round(frac_ev, 4))
}

out <- rbind(run_state("NY_SPARCS", load_ny()), run_state("TX_PUDF", load_tx()))
write.csv(out, "data/audit_splitval.csv", row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), paste("RNGkind:", paste(RNGkind(), collapse = ","))),
           "data/_audit_splitval_sessionInfo.txt")
cat(sprintf("saved data/audit_splitval.csv (NSPLIT=%d)\n", NSPLIT))
