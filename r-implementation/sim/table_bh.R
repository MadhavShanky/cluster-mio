# table_bh.R -- compact Supplement table from sim/agg/agg_bh.csv (run sim/aggregate_bh.R first).
# Rows: design group x delta; columns: mean TPR / FDR per method, averaged over the cells in the group
# (K in {100, 300}, pi in the group's grid). Base design split by balanced / unequal cluster sizes.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); library(dplyr); library(tidyr)
}))
a <- read.csv("sim/agg/agg_bh.csv")
a$group <- ifelse(a$design == "base", ifelse(a$unequal, "independent, unequal", "independent, balanced"),
                  ifelse(a$design == "clustercv", "cluster-level covariate", "between-cluster covariates"))
g <- a %>% filter(method %in% c("BH05", "Bonf05", "SCS", "SCS_CPSS"), delta %in% c(1, 2)) %>%
  group_by(group, delta, method) %>% summarise(tpr = mean(tpr), fdr = mean(fdr), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = c(tpr, fdr))
ord <- c("independent, balanced", "independent, unequal", "cluster-level covariate", "between-cluster covariates")
g <- g %>% mutate(group = factor(group, ord)) %>% arrange(group, delta)
print(as.data.frame(g %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
k100 <- a %>% filter(design != "base", K == 100, method == "BH05") %>% summarise(lo = min(fdr), hi = max(fdr))
cat(sprintf("\nBH05 realized FDR in the two hard designs at K=100: %.3f-%.3f\n", k100$lo, k100$hi))
fmt <- function(t, f) sprintf("%.2f / %.3f", t, f)
lines <- with(g, sprintf("%s & %g & %s & %s & %s & %s\\\\", group, delta, fmt(tpr_BH05, fdr_BH05),
                         fmt(tpr_Bonf05, fdr_Bonf05), fmt(tpr_SCS, fdr_SCS), fmt(tpr_SCS_CPSS, fdr_SCS_CPSS)))
writeLines(lines, "sim/agg/table_bh_rows.tex"); cat(lines, sep = "\n")
