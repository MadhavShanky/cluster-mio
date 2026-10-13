# tx_sensitivity.R -- risk-adjustment sensitivity of the TX PUDF watch-list.
# MIRROR of realdata_sensitivity.R (NY SPARCS) for the Texas cohort, so the manuscript can report the
# cross-state replication's sensitivity instead of leaving Texas "unchecked" (ref2 item O). Cleaning is
# copied verbatim from tx_hf_analysis.R so the cohort matches the published K=208 analysis.
# Outputs data/tx/sensitivity.csv : per spec, the stable watch-list size + Jaccard vs the full spec.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R")
}))
MIN_NK <- 25L
d <- read.csv("data/tx/tx_hf_1q2019.csv", stringsAsFactors = FALSE,
              colClasses = c(thcic_id = "character", age = "character",
                             rom = "character", sev = "character", adm = "character"))
## clean -- identical to tx_hf_analysis.R ----------------------------------------------------------
d$sex <- trimws(d$sex)
d <- d[d$sex %in% c("M","F") & !is.na(d$los) & d$los > 0, ]
d$sev  <- factor(d$sev, levels = c("1","2","3","4"))        # 0/blank ("no class") -> NA, dropped
d$rom  <- factor(d$rom, levels = c("1","2","3","4"))
d$age  <- factor(d$age); d$sexf <- factor(d$sex); d$adm <- factor(d$adm)
d <- d[complete.cases(d[, c("sev","rom","age","adm","sexf")]), ]
tb <- table(d$thcic_id); d <- d[d$thcic_id %in% names(tb)[tb >= MIN_NK], ]
fac <- droplevels(factor(d$thcic_id))
for (v in c("age","sexf","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
y <- log(d$los); K <- nlevels(fac)
cat(sprintf("TX cohort N=%d K=%d\n", nrow(d), K))

# FE specifications: full + leave-one-domain-out + minimal (mirror NY specs in realdata_sensitivity.R)
specs <- list(
  full       = ~ age + sexf + sev + rom + adm,
  drop_sev   = ~ age + sexf + rom + adm,
  drop_rom   = ~ age + sexf + sev + adm,
  drop_age   = ~ sexf + sev + rom + adm,
  minimal    = ~ age + sexf + adm
)
jacc <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

# fix the budget from the full spec so specs are compared at the same audit size
Xf <- model.matrix(specs$full, data = d)[, -1, drop = FALSE]
lam <- max(1L, audit_budget(scs_path(Xf, y, fac, lambda_max = 25L, n_restart = 4L, seed = 1L), target = 0.90)$suggested_budget)
cat(sprintf("budget lambda=%d (from full spec)\n", lam))

flag <- list(); rows <- list()
for (nm in names(specs)) {
  X <- model.matrix(specs[[nm]], data = d)[, -1, drop = FALSE]
  fit <- scs_fit(X, y, fac, lambda = lam, n_restart = 6L, seed = 1L)
  st  <- stability(fit, B = 400L, pi_thr = 0.8, seed = 1L, pairing = "complementary")
  flag[[nm]] <- as.character(st$flagged)
  rows[[nm]] <- data.frame(spec = nm, p = ncol(X), n_flagged = length(st$flagged),
                           mb_bound = round(st$expected_false_flags_bound, 2),
                           jacc_vs_full = round(jacc(flag[[nm]], flag$full), 3),
                           n_shared_full = length(intersect(flag[[nm]], flag$full)))
  cat(sprintf("  %-9s p=%d  stable=%d  Jaccard_vs_full=%.3f\n", nm, ncol(X), length(st$flagged), rows[[nm]]$jacc_vs_full))
}
res <- do.call(rbind, rows)
write.csv(res, "data/tx/sensitivity.csv", row.names = FALSE)
saveRDS(flag, "data/tx/sensitivity_flags.rds")  # per-spec flagged ids, for flagged_table.R core marking
core <- Reduce(intersect, flag)
cat(sprintf("\ncore watch-list (flagged under ALL %d specs): %d hospitals\n", length(specs), length(core)))
cat(sprintf("full-spec watch-list: %d; mean Jaccard of reduced specs vs full: %.3f\n",
            length(flag$full), mean(sapply(setdiff(names(specs), "full"), function(s) rows[[s]]$jacc_vs_full))))
cat("saved data/tx/sensitivity.csv\n")
