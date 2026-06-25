# figures.R -- aggregate sim/results/cell_*.csv into the paper figures (F1-F4 + S-figures).
# Robust to partial results (guards empty subsets); fills in as cells land. Uses sim/theme_scs.R.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("sim/theme_scs.R"); library(dplyr)
}))
dir.create("sim/figs", showWarnings = FALSE)
REF <- function(d) dplyr::filter(d, err == "normal", !unequal, nk == 30, p == 5, rho_x == 0.5)  # reference slice

files <- list.files("sim/results", "^cell_[0-9]+\\.csv$", full.names = TRUE)
raw <- bind_rows(lapply(files, function(f) suppressWarnings(read.csv(f, stringsAsFactors = FALSE))))
keyc <- c("K","nk","p","rho_x","regime","pi","delta","icc","err","unequal")
agg <- raw %>% group_by(across(all_of(keyc))) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), reps = dplyr::n(), .groups = "drop")
cat(sprintf("aggregated %d reps over %d cells (regimes: %s)\n", nrow(raw), nrow(agg),
            paste(sort(unique(agg$regime)), collapse = ", ")))

save_if <- function(p, file, w = 8, h = 4.4) if (!is.null(p)) scs_save(p, file, w = w, h = h)

## ---- F1: phase map -- L0 estimation advantage over LMM (headline figure: sparsity x K x signal) -------------
f1 <- REF(agg) %>% filter(regime == "sparse") %>% mutate(adv = log2(gerr_lmm / gerr_l0))
pF1 <- if (nrow(f1)) {
  lim <- c(-1, 1) * max(abs(f1$adv), na.rm = TRUE)
  scs_phase_tile(f1, x = "pi", y = "K", value = "adv", diverging = TRUE, facet = "delta", limit = lim,
                 title = "Where L0 cluster-selection beats shrinkage",
                 subtitle = "gamma-recovery advantage over LMM-Gaussian, sparse regime (faceted by signal delta)",
                 xlab = "fraction of clusters deviating (pi)", value_name = "L0 advantage\nlog2(err_LMM/err_L0)")
} else NULL
save_if(pF1, "sim/figs/F1_phase_map.png")

## ---- F2: selection FDR/power -- TPR & FDR vs delta, methods, faceted by K ------------------------
mk_long <- function(d, mets) bind_rows(lapply(names(mets), function(m) {
  cols <- mets[[m]]; data.frame(K = d$K, delta = d$delta, method = m,
    tpr = if (!is.na(cols[1])) d[[cols[1]]] else NA, fdr = if (!is.na(cols[2])) d[[cols[2]]] else NA) }))
f2src <- REF(agg) %>% filter(regime == "sparse", pi == 0.1)
pF2 <- if (nrow(f2src)) {
  L <- mk_long(f2src, list("L0-MIO" = c("tpr_l0","fdr_l0"), "L0-MIO+stab" = c(NA,"fdr_stab"),
                           "LMM-Gaussian" = c("tpr_lmm","fdr_lmm")))
  ggplot(L, aes(delta, colour = method)) +
    geom_line(aes(y = tpr, linetype = "TPR")) + geom_point(aes(y = tpr)) +
    geom_line(aes(y = fdr, linetype = "FDR")) +
    facet_wrap(~K, nrow = 1) + scale_color_scs(name = NULL) +
    scale_linetype_manual(values = c(TPR = "solid", FDR = "dashed"), name = NULL) +
    labs(title = "Selection power and false-discovery vs signal strength",
         subtitle = "sparse regime, pi = 0.1; TPR (solid) and FDR (dashed) by method, faceted by K",
         x = expression("signal "*delta), y = "rate") + theme_scs()
} else NULL
save_if(pF2, "sim/figs/F2_selection.png")

## ---- F3: computational scaling -- L0 runtime vs K (log-log) --------------------------------------
f3 <- REF(agg) %>% filter(regime == "sparse", pi == 0.1, delta == 2)
pF3 <- if (nrow(f3) > 1) ggplot(f3, aes(K, runtime_l0)) +
  geom_line(colour = scs_cols["L0-MIO"], linewidth = 0.8) + geom_point(colour = scs_cols["L0-MIO"], size = 1.8) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Computational scaling of the L0 solver",
       subtitle = "wall-clock per fit incl. CV, vs number of clusters K (log-log); sparse pi=0.1, delta=2",
       x = "clusters K", y = "seconds (CV + final fit)",
       caption = "Solver matches brute force on all checked instances.") +
  theme_scs() else NULL
save_if(pF3, "sim/figs/F3_scaling.png", w = 6.5)

## ---- F4: capabilities -- (a) Sx prediction gate, (b) Gaussian-regime check ----------------------------
sxd <- REF(agg) %>% filter(regime == "sparse_x")
pF4a <- if (nrow(sxd)) {
  D <- bind_rows(data.frame(K = sxd$K, pi = sxd$pi, method = "L0-MIO", mse = sxd$cbmse_l0),
                 data.frame(K = sxd$K, pi = sxd$pi, method = "LMM-Gaussian", mse = sxd$cbmse_lmm))
  ggplot(D, aes(factor(pi), mse, fill = method)) + geom_col(position = "dodge") + facet_wrap(~K, nrow = 1) +
    scale_fill_scs(name = NULL) +
    labs(title = "New-cluster prediction gate (regime Sx: gamma linked to a prediction-only covariate)",
         subtitle = "cluster-blocked test-MSE via a Z->gamma rule; L0's unshrunk sparse effects predict new clusters better",
         x = "fraction deviating (pi)", y = "cluster-blocked test MSE") + theme_scs()
} else NULL
save_if(pF4a, "sim/figs/F4a_prediction_gate.png")

gd <- agg %>% filter(regime == "gaussian", K != 10)  # K=10 gaussian is a single stray ICC cell, not part of the {30,100,300} sweep
pF4b <- if (nrow(gd)) {
  D <- bind_rows(data.frame(icc = gd$icc_true, K = gd$K, method = "L0-MIO", gerr = gd$gerr_l0),
                 data.frame(icc = gd$icc_true, K = gd$K, method = "LMM-Gaussian", gerr = gd$gerr_lmm))
  ggplot(D, aes(icc, gerr, colour = method)) + geom_line() + geom_point() + facet_wrap(~K, nrow = 1) +
    scale_color_scs(name = NULL) +
    labs(title = "Gaussian random-effect DGP: gamma-recovery error vs ICC",
         subtitle = "gamma-recovery error vs true ICC, faceted by K; L0 target is misspecified under this DGP",
         x = "true ICC", y = expression("||"*hat(gamma)-gamma*"||")) + theme_scs()
} else NULL
save_if(pF4b, "sim/figs/F4b_gaussian_honesty.png")

## ---- S: ICC recovery + robustness sweeps ---------------------------------------------------------
pS_icc <- if (nrow(gd)) ggplot(gd, aes(icc_true, icc_lmm)) + geom_abline(linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = factor(K)), size = 2) + geom_line(aes(colour = factor(K))) +
  labs(title = "ICC recovery (LMM-Gaussian)", x = "true ICC", y = "estimated ICC", colour = "K") + theme_scs() else NULL
save_if(pS_icc, "sim/figs/S_icc_recovery.png", w = 6)

# robustness: reference sparse cell vs its sweeps (unequal n_k, t-errors, etc.)
sw <- agg %>% filter(regime == "sparse", K == 100, pi == 0.1, delta == 2) %>%
  mutate(setting = case_when(unequal ~ "unequal n_k", err == "t" ~ "t errors", p == 25 ~ "p=25",
                             rho_x == 0 ~ "rho_x=0", nk == 20 ~ "n_k=20", nk == 50 ~ "n_k=50", TRUE ~ "reference"))
pS_rob <- if (nrow(sw) > 1) ggplot(sw, aes(reorder(setting, gerr_l0))) +
  geom_point(aes(y = gerr_l0, colour = "L0-MIO"), size = 2.5) + geom_point(aes(y = gerr_lmm, colour = "LMM-Gaussian"), size = 2.5) +
  scale_color_scs(name = NULL) + coord_flip() +
  labs(title = "Robustness sweeps (K=100, pi=0.1, delta=2)", x = NULL, y = expression("||"*hat(gamma)-gamma*"||")) +
  theme_scs() else NULL
save_if(pS_rob, "sim/figs/S_robustness.png", w = 6.5)

cat("figures written to sim/figs/ (those with data)\n")
