# aggregate_bh.R -- summarize sim/results_bh/cell_*.csv (bh_comparator_sim.R).
# Per (K, pi, delta, unequal, method): mean TPR = TP / n_active, mean FDR = FP / max(n_sel, 1),
# mean number of false flags V = FP, mean selected count; for SCS_CPSS also the mean CPSS reference and
# P(V > reference). Writes sim/agg/agg_bh.csv and prints a compact table.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); library(dplyr)
}))
dirs <- c(base = "sim/results_bh", clustercv = "sim/results_bh_clustercv", highcorr = "sim/results_bh_highcorr")
fs <- unlist(lapply(dirs, list.files, pattern = "^cell_\\d+\\.csv$", full.names = TRUE))
d <- bind_rows(lapply(fs, function(f) {
  x <- read.csv(f); if (!"design" %in% names(x)) x$design <- "base"; x }))
cat(sprintf("%d files, %d rows, reps per cell: %s\n", length(fs), nrow(d),
            paste(range(table(d$cell[d$method == "SCS"])), collapse = "-")))
a <- d %>% mutate(tpr = TP / n_active, fdr = ifelse(n_sel > 0, FP / n_sel, 0)) %>%
  group_by(design, K, pi, delta, unequal, method) %>%
  summarise(R = n(), tpr = mean(tpr), fdr = mean(fdr), V = mean(FP), n_sel = mean(n_sel),
            ref = mean(ref), p_exceed = mean(FP > ref), se_fdr = sd(FP / pmax(n_sel, 1)) / sqrt(n()),
            .groups = "drop")
dir.create("sim/agg", showWarnings = FALSE)
write.csv(a, "sim/agg/agg_bh.csv", row.names = FALSE)
w <- a %>% select(design, K, pi, delta, unequal, method, tpr, fdr) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(tpr, fdr))
options(width = 220)
print(as.data.frame(w %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
cat("\nCPSS: mean V vs mean reference, P(V > ref)\n")
print(as.data.frame(a %>% filter(method == "SCS_CPSS") %>% select(design, K, pi, delta, unequal, V, ref, p_exceed) %>%
                      mutate(across(c(V, ref, p_exceed), ~ round(.x, 3)))), row.names = FALSE)
