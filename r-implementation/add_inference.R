# add_inference.R -- deliverable A (stability-selection FDR) + BLUP comparison, on BOTH real cohorts.
# Reuses the saved symmetric SCS fits (fit$data holds X/y/cluster). Shows: (1) the L0 flagged hospitals
# survive a Meinshausen-Buhlmann stability bound on expected false flags; (2) the gate L0 opens that BLUP
# cannot -- a sparse, FDR-controlled audit list vs. K shrunken-toward-zero estimates with no native cutoff.
suppressWarnings(suppressMessages({
  setwd("C:/Users/tomch/AIProjects/research/MIO LMM"); source("scs.R"); library(lme4)
}))

analyze <- function(rds, label) {
  R <- readRDS(rds); fit <- R$fit
  X <- fit$data$X; y <- fit$data$y; cl <- fit$data$cluster
  K <- length(fit$levels)
  cat(sprintf("\n================ %s : N=%d, K=%d, point-flagged=%d ================\n",
              label, length(y), K, length(fit$selected)))

  ## --- deliverable A: stability selection + MB expected-false-flag bound ---
  st <- stability(fit, B = 100L, frac = 0.5, pi_thr = 0.6, seed = 1L)
  cat(sprintf("Stability selection (B_used=%d, avg q=%.2f selected/subsample):\n", st$B_used, st$q))
  cat(sprintf("  MB bound on E[# false flags] at pi_thr=%.2f: %.2f\n", st$pi_thr, st$expected_false_flags_bound))
  conf <- intersect(fit$selected, st$flagged)
  cat(sprintf("  stable flagged set (freq>=%.2f): %d hospitals; %d/%d point-flagged survive\n",
              st$pi_thr, length(st$flagged), length(conf), length(fit$selected)))
  cat("  top stable hospitals (selection frequency):\n")
  print(round(utils::head(st$freq[st$flagged], 10), 2))

  ## --- BLUP comparison: the capability gap ---
  df <- data.frame(y = y, X, cluster = cl)
  fm <- suppressWarnings(suppressMessages(lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  blup <- ranef(fm)$cluster[, 1]; names(blup) <- rownames(ranef(fm)$cluster)
  blup <- blup[as.character(fit$levels)]; blup[is.na(blup)] <- 0
  g0 <- fit$gamma[as.character(fit$levels)]
  lam <- length(fit$selected)
  top_blup <- names(sort(abs(blup), decreasing = TRUE))[seq_len(max(1L, lam))]
  ov <- length(intersect(fit$selected, top_blup))
  cat(sprintf("\nBLUP vs L0:\n"))
  cat(sprintf("  overlap of L0 flagged with top-%d |BLUP| hospitals: %d/%d\n", lam, ov, lam))
  cat(sprintf("  corr(L0 gamma, BLUP) over all K: %.3f\n", cor(g0, blup)))
  cat(sprintf("  exact zeros: BLUP %d/%d vs L0 %d/%d  -> only L0 yields a sparse audit list\n",
              sum(blup == 0), K, sum(g0 == 0), K))
  cat(sprintf("  max |effect|: BLUP %.2f vs L0 %.2f  (BLUP shrinks the extremes toward 0)\n",
              max(abs(blup)), max(abs(g0))))
  invisible(list(stability = st, blup = blup, g0 = g0, K = K, flagged = fit$selected, conf = conf))
}

s1 <- analyze("data/sparcs/hf_scs_result.rds", "NY SPARCS HF")
s2 <- analyze("data/tx/tx_hf_scs_result.rds",  "TX PUDF HF")
saveRDS(list(sparcs = s1, tx = s2), "data/inference_both.rds")
cat("\n\n=== cross-state summary ===\n")
cat(sprintf("NY: %d/%d point-flagged survive stability; MB false-flag bound %.2f\n",
            length(s1$conf), length(s1$flagged), s1$stability$expected_false_flags_bound))
cat(sprintf("TX: %d/%d point-flagged survive stability; MB false-flag bound %.2f\n",
            length(s2$conf), length(s2$flagged), s2$stability$expected_false_flags_bound))
cat("saved data/inference_both.rds\n")
