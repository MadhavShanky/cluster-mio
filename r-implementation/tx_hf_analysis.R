# tx_hf_analysis.R -- TX PUDF hospital-profiling, parallel to pull_sparcs_hf.R (NY SPARCS).
# Cohort: APR-DRG 194 (Heart Failure), 1Q2019, parsed from the fixed-width PUDF base1 file
# (positions from InpatientDataDictionary1Q2019.pdf). Outcome = log(LOS), risk-adjusted,
# hospital = THCIC_ID. Runs SCS audit-budget (C), cluster ranking (D), symmetric flag, asymmetric tails (B).
suppressWarnings(suppressMessages({
  setwd("C:/Users/tomch/AIProjects/research/MIO LMM"); source("scs.R")
}))
MIN_NK <- 25L
d <- read.csv("data/tx/tx_hf_1q2019.csv", stringsAsFactors = FALSE,
              colClasses = c(thcic_id = "character", age = "character",
                             rom = "character", sev = "character", adm = "character"))
cat(sprintf("raw HF rows: %d, hospitals: %d\n", nrow(d), length(unique(d$thcic_id))))

## clean -------------------------------------------------------------------------------------------
d$sex <- trimws(d$sex)
d <- d[d$sex %in% c("M","F") & !is.na(d$los) & d$los > 0, ]
d$sev  <- factor(d$sev, levels = c("1","2","3","4"))        # 0/blank ("no class") -> NA, dropped
d$rom  <- factor(d$rom, levels = c("1","2","3","4"))
d$age  <- factor(d$age); d$sexf <- factor(d$sex); d$adm <- factor(d$adm)
d <- d[complete.cases(d[, c("sev","rom","age","adm","sexf")]), ]
## drop low-volume hospitals so each effect is estimable
tb <- table(d$thcic_id); d <- d[d$thcic_id %in% names(tb)[tb >= MIN_NK], ]
fac <- droplevels(factor(d$thcic_id))
for (v in c("age","sexf","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
y <- log(d$los)
X <- model.matrix(~ age + sexf + sev + rom + adm, data = d)[, -1, drop = FALSE]
K <- nlevels(fac)
cat(sprintf("analysis cohort: N=%d, K=%d hospitals (>= %d HF discharges), p=%d covariates\n",
            nrow(d), K, MIN_NK, ncol(X)))
cat(sprintf("log-LOS: mean %.2f sd %.2f | raw LOS median %.0f\n", mean(y), sd(y), exp(median(y))))

## SCS: audit-budget path (C) + cluster ranking (D) ------------------------------------------------
pth <- scs_path(X, y, fac, lambda_max = 25L, n_restart = 4L, seed = 1L)
ab  <- audit_budget(pth, target = 0.90)
cat("\n--- audit-budget (C): hospitals to flag ---\n")
cat(sprintf("suggested budget (90%% of detectable deviation): lambda = %d | elbow: lambda = %d\n",
            ab$suggested_budget, ab$suggested_budget_elbow))
rk <- cluster_ranking(pth)
cat("\n--- top-12 most atypical hospitals (D: lambda-path entry order) ---\n")
print(utils::head(rk$ranking, 12), row.names = FALSE)

## symmetric flag at the suggested budget ----------------------------------------------------------
lam <- max(1L, ab$suggested_budget)
fit <- scs_fit(X, y, fac, lambda = lam, n_restart = 6L, seed = 1L)
g <- sort(fit$gamma[fit$selected], decreasing = TRUE)
cat(sprintf("\n--- flagged %d hospitals at lambda=%d (gamma = adjusted log-LOS deviation) ---\n", length(g), lam))
print(round(g, 3))
cat("  (exp(gamma) = multiplicative effect on risk-adjusted LOS; >1 longer, <1 shorter than expected)\n")

## asymmetric worst/best tails (B) -----------------------------------------------------------------
fitA <- scs_fit(X, y, fac, lambda = c(worst = 5, best = 5), n_restart = 6L, seed = 1L)
cat("\n--- asymmetric audit list (B): 5 worst (longest) + 5 best (shortest) adjusted LOS ---\n")
cat("WORST:\n"); print(round(sort(fitA$gamma[fitA$selected_pos], decreasing = TRUE), 3))
cat("BEST:\n");  print(round(sort(fitA$gamma[fitA$selected_neg]), 3))

saveRDS(list(path = pth, audit = ab, ranking = rk, fit = fit, fitA = fitA,
             N = nrow(d), K = K, state = "TX", drg = "194", quarter = "1Q2019"),
        "data/tx/tx_hf_scs_result.rds")
cat("\nsaved data/tx/tx_hf_scs_result.rds\n")
