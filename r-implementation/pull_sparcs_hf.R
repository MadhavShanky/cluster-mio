# pull_sparcs_hf.R -- real-data hospital-profiling demo on NY SPARCS (open API).
# Cohort: APR-DRG 194 (Heart Failure), discharge_year 2024. Outcome = log(LOS), risk-adjusted,
# hospital = permanent_facility_id. Runs SCS: audit-budget path (C), cluster ranking (D),
# symmetric flag, and asymmetric worst/best tails (B). Default cohort -- change DRG via APRDRG below.
suppressWarnings(suppressMessages({
  setwd("C:/Users/tomch/AIProjects/research/MIO LMM"); source("scs.R")
}))
APRDRG <- "194"; YEAR <- "2024"; MIN_NK <- 25L   # min discharges/hospital for a stable effect

dir.create("data/sparcs", showWarnings = FALSE)
base  <- "https://health.data.ny.gov/resource/sf4k-39ay.csv"
where <- sprintf("apr_drg_code='%s' AND discharge_year='%s'", APRDRG, YEAR)
url   <- paste0(base, "?$where=", gsub(" ", "%20", where), "&$limit=200000")
dest  <- sprintf("data/sparcs/drg%s_%s.csv", APRDRG, YEAR)
download.file(url, dest, mode = "wb", quiet = TRUE)
d <- read.csv(dest, stringsAsFactors = FALSE)
cat(sprintf("pulled %d rows for APR-DRG %s (%s)\n", nrow(d), APRDRG, YEAR))

## clean -------------------------------------------------------------------------------------------
los <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d$length_of_stay)))   # "120 +" -> 120
sev_lv <- c("Minor", "Moderate", "Major", "Extreme")
d$sev <- factor(d$apr_severity_of_illness, levels = sev_lv)
d$rom <- factor(d$apr_risk_of_mortality,  levels = sev_lv)
d$age <- factor(d$age_group); d$sex <- factor(d$gender); d$adm <- factor(d$type_of_admission)
d$fac <- d$permanent_facility_id
d$los <- los
keep <- !is.na(d$los) & d$los > 0 & d$sex %in% c("M","F") &
        complete.cases(d[, c("sev","rom","age","sex","adm","fac")])
d <- d[keep, ]
## drop low-volume hospitals so each hospital effect is estimable
tb <- table(d$fac); d <- d[d$fac %in% names(tb)[tb >= MIN_NK], ]
fac <- droplevels(factor(d$fac))
for (v in c("age","sex","sev","rom","adm")) d[[v]] <- droplevels(d[[v]])
y <- log(d$los)
X <- model.matrix(~ age + sex + sev + rom + adm, data = d)[, -1, drop = FALSE]
K <- nlevels(fac)
cat(sprintf("analysis cohort: N=%d, K=%d hospitals (>= %d discharges), p=%d covariates\n",
            nrow(d), K, MIN_NK, ncol(X)))
cat(sprintf("log-LOS: mean %.2f sd %.2f | raw LOS median %.0f\n", mean(y), sd(y), exp(median(y))))

## SCS: audit-budget path (C) + cluster ranking (D) ------------------------------------------------
pth <- scs_path(X, y, fac, lambda_max = 25L, n_restart = 4L, seed = 1L)
ab  <- audit_budget(pth, target = 0.90)
cat("\n--- audit-budget (deliverable C): how many hospitals to flag ---\n")
cat(sprintf("suggested budget (90%% of detectable deviation): lambda = %d\n", ab$suggested_budget))
cat(sprintf("marginal-gain elbow: lambda = %d\n", ab$suggested_budget_elbow))
rk <- cluster_ranking(pth)
cat("\n--- top-12 most atypical hospitals (deliverable D: lambda-path entry order) ---\n")
print(utils::head(rk$ranking, 12), row.names = FALSE)

## symmetric flag at the suggested budget ----------------------------------------------------------
lam <- max(1L, ab$suggested_budget)
fit <- scs_fit(X, y, fac, lambda = lam, n_restart = 6L, seed = 1L)
g <- sort(fit$gamma[fit$selected], decreasing = TRUE)
cat(sprintf("\n--- flagged %d hospitals at lambda=%d (gamma = adjusted log-LOS deviation) ---\n",
            length(g), lam))
print(round(g, 3))
cat("  (exp(gamma) = multiplicative effect on risk-adjusted LOS; >1 longer, <1 shorter than expected)\n")

## asymmetric worst/best tails (deliverable B) -----------------------------------------------------
fitA <- scs_fit(X, y, fac, lambda = c(worst = 5, best = 5), n_restart = 6L, seed = 1L)
cat("\n--- asymmetric audit list (deliverable B): 5 worst (longest) + 5 best (shortest) adjusted LOS ---\n")
cat("WORST (longer-than-expected LOS):\n"); print(round(sort(fitA$gamma[fitA$selected_pos], decreasing = TRUE), 3))
cat("BEST (shorter-than-expected LOS):\n");  print(round(sort(fitA$gamma[fitA$selected_neg]), 3))

saveRDS(list(path = pth, audit = ab, ranking = rk, fit = fit, fitA = fitA,
             N = nrow(d), K = K, drg = APRDRG, year = YEAR),
        "data/sparcs/hf_scs_result.rds")
cat("\nsaved data/sparcs/hf_scs_result.rds\n")
