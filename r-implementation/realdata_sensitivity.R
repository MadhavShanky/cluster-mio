# realdata_sensitivity.R -- risk-adjustment sensitivity of the SPARCS watch-list.
# Vary the fixed-effect (case-mix) specification and measure how the stable flagged set moves.
# Outputs data/sparcs/sensitivity.csv : per spec, the stable watch-list size + Jaccard vs the full spec.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R")
}))
d <- read.csv("data/sparcs/drg194_2024.csv", stringsAsFactors = FALSE)
los <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$length_of_stay)))
sev_lv <- c("Minor", "Moderate", "Major", "Extreme")
d$sev <- factor(d$apr_severity_of_illness, levels = sev_lv)
d$rom <- factor(d$apr_risk_of_mortality,  levels = sev_lv)
d$age <- factor(d$age_group); d$sex <- factor(d$gender); d$adm <- factor(d$type_of_admission)
d$fac <- d$permanent_facility_id; d$los <- los
keep <- !is.na(d$los) & d$los > 0 & d$sex %in% c("M","F") &
        complete.cases(d[, c("sev","rom","age","sex","adm","fac")])
d <- d[keep, ]
tb <- table(d$fac); d <- d[d$fac %in% names(tb)[tb >= 25L], ]
fac <- droplevels(factor(d$fac)); for (v in c("age","sex","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
y <- log(d$los); K <- nlevels(fac)
cat(sprintf("cohort N=%d K=%d\n", nrow(d), K))

# FE specifications: full + leave-one-domain-out + minimal
specs <- list(
  full       = ~ age + sex + sev + rom + adm,
  drop_sev   = ~ age + sex + rom + adm,
  drop_rom   = ~ age + sex + sev + adm,
  drop_age   = ~ sex + sev + rom + adm,
  minimal    = ~ age + sex + adm
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
write.csv(res, "data/sparcs/sensitivity.csv", row.names = FALSE)
saveRDS(flag, "data/sparcs/sensitivity_flags.rds")  # per-spec flagged ids, for flagged_table.R core marking
# core = hospitals flagged under EVERY spec
core <- Reduce(intersect, flag)
cat(sprintf("\ncore watch-list (flagged under ALL %d specs): %d hospitals\n", length(specs), length(core)))
cat(sprintf("full-spec watch-list: %d; mean Jaccard of reduced specs vs full: %.3f\n",
            length(flag$full), mean(sapply(setdiff(names(specs), "full"), function(s) rows[[s]]$jacc_vs_full))))
cat("saved data/sparcs/sensitivity.csv\n")
