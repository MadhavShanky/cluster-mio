# tables.R -- "summary of all simulation results" main-text table from sim/results/cell_*.csv.
# Rows = representative scenarios; cols = methods x {selection F1, gamma-error, within-cluster MSE, runtime}.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); library(dplyr)
}))
files <- list.files("sim/results", "^cell_[0-9]+\\.csv$", full.names = TRUE)
raw <- bind_rows(lapply(files, function(f) suppressWarnings(read.csv(f, stringsAsFactors = FALSE))))
keyc <- c("K","nk","p","rho_x","regime","pi","delta","icc","err","unequal")
agg <- raw %>% group_by(across(all_of(keyc))) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), reps = dplyr::n(), .groups = "drop")
ref <- function(d) dplyr::filter(d, err == "normal", !unequal, nk == 30, p == 5, rho_x == 0.5)

scen <- bind_rows(
  ref(agg) %>% filter(regime == "sparse",   K == 100, delta == 2, pi %in% c(0.05, 0.2, 0.6)),
  ref(agg) %>% filter(regime == "sparse_x", K == 100, delta == 2, pi == 0.2),
  ref(agg) %>% filter(regime == "gaussian", K == 100, icc == 0.5))

g <- function(col) if (col %in% names(scen)) round(scen[[col]], 3) else NA
tab <- data.frame(
  scenario = with(scen, ifelse(regime == "gaussian", sprintf("Gaussian ICC=%.1f", icc),
                        sprintf("%s pi=%.2f delta=%g", regime, pi, delta))),
  K = scen$K, reps = scen$reps,
  F1_L0 = g("f1_l0"), F1_stab = g("f1_stab"), F1_glasso = g("f1_glasso"), F1_LMM = g("f1_lmm"),
  FDR_L0 = g("fdr_l0"), FDR_stab = g("fdr_stab"),
  gerr_L0 = g("gerr_l0"), gerr_LMM = g("gerr_lmm"), gerr_glasso = g("gerr_glasso"),
  wMSE_L0 = g("wmse_l0"), wMSE_LMM = g("wmse_lmm"),
  cbMSE_L0 = g("cbmse_l0"), cbMSE_LMM = g("cbmse_lmm"),
  runtime_L0_s = g("runtime_l0"), stringsAsFactors = FALSE)

write.csv(tab, "sim/summary_table.csv", row.names = FALSE)
# markdown for quick paste into notes / manuscript draft
md <- c(paste0("| ", paste(names(tab), collapse = " | "), " |"),
        paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"),
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, "sim/summary_table.md")
cat("=== summary table (representative scenarios) ===\n"); print(tab, row.names = FALSE)
cat("\nwrote sim/summary_table.csv + .md\n")
