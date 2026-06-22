# figures_realdata.R -- paper figures for the real-data hospital-profiling section (NY SPARCS + TX PUDF).
# Reuses saved fits + inference (no refitting). Three figures, faceted by state; stability at pi_thr = 0.8.
suppressWarnings(suppressMessages({
  setwd("C:/Users/tomch/AIProjects/research/MIO LMM"); source("sim/theme_scs.R")
}))
PI_THR <- 0.8
dir.create("figs", showWarnings = FALSE)
I  <- readRDS("data/inference_both.rds")
RS <- readRDS("data/sparcs/hf_scs_result.rds")
RT <- readRDS("data/tx/tx_hf_scs_result.rds")
NY <- "NY SPARCS  (HF, K=161)"; TX <- "TX PUDF  (HF, K=208)"

## re-threshold stability at pi_thr = 0.8 + recompute the MB bound (q, K from the saved object) ------
stab_at <- function(s, thr) {
  fr <- s$stability$freq; q <- s$stability$q; K <- s$stability$K
  list(freq = fr, flagged = names(fr)[fr >= thr],
       bound = q^2 / ((2 * thr - 1) * K), q = q, K = K)
}
b1 <- stab_at(I$sparcs, PI_THR); b2 <- stab_at(I$tx, PI_THR)
cat(sprintf("pi_thr=%.2f | NY: %d stable, MB bound %.2f  ||  TX: %d stable, MB bound %.2f\n",
            PI_THR, length(b1$flagged), b1$bound, length(b2$flagged), b2$bound))

## ---- FIG R1: caterpillar -- L0 (sparse, un-shrunk) vs BLUP (dense, shrunk) -----------------------
mk_cat <- function(s, state) {
  d <- data.frame(hosp = names(s$blup), blup = as.numeric(s$blup),
                  l0 = as.numeric(s$g0[names(s$blup)]),
                  flagged = names(s$blup) %in% s$flagged, state = state)
  d <- d[order(-d$blup), ]; d$rank <- seq_len(nrow(d)); d
}
catdf <- rbind(mk_cat(I$sparcs, NY), mk_cat(I$tx, TX))
fl <- catdf[catdf$flagged, ]
pR1 <- ggplot(catdf, aes(rank)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(aes(y = blup, colour = "LMM-Gaussian"), size = 0.9, alpha = 0.55) +
  geom_segment(data = fl, aes(x = rank, xend = rank, y = blup, yend = l0), colour = "grey70", linewidth = 0.3) +
  geom_point(data = fl, aes(y = l0, colour = "L0-MIO"), size = 1.9) +
  facet_wrap(~state, scales = "free_x") +
  scale_color_scs(name = NULL) +
  labs(title = "L0 flags a sparse set and reports it un-shrunk; BLUP spreads and shrinks",
       subtitle = "Risk-adjusted hospital effect on log length-of-stay; hospitals ranked by BLUP",
       x = "hospital (ranked by BLUP random intercept)",
       y = expression(hat(gamma)[k]~"  (log-LOS deviation)"),
       caption = "Forest = L0-MIO flagged hospitals (unshrunk LS refit); blue = LMM-Gaussian BLUP (all hospitals, none exactly 0). Segment = shrinkage gap.") +
  theme_scs()
scs_save(pR1, "figs/R1_caterpillar_L0_vs_BLUP.png", w = 8.4, h = 4.4)

## ---- FIG R2: stability selection (pi_thr = 0.8) --------------------------------------------------
mk_stab <- function(s, state) {
  fr <- s$stability$freq
  data.frame(hosp = names(fr), freq = as.numeric(fr), rank = seq_along(fr), state = state)
}
stabdf <- rbind(mk_stab(I$sparcs, NY), mk_stab(I$tx, TX))
stabdf$grp <- ifelse(stabdf$freq >= PI_THR, "L0-MIO+stab", "OLS")  # forest vs grey via palette
lab <- data.frame(state = c(NY, TX),
                  txt = c(sprintf("%d stable (freq>=%.1f)\nMB E[false]<=%.1f", length(b1$flagged), PI_THR, b1$bound),
                          sprintf("%d stable (freq>=%.1f)\nMB E[false]<=%.1f", length(b2$flagged), PI_THR, b2$bound)))
pR2 <- ggplot(stabdf, aes(rank, freq, colour = grp)) +
  geom_hline(yintercept = PI_THR, linetype = "dashed", colour = scs_emphasis, linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text(data = lab, aes(x = Inf, y = 0.05, label = txt), inherit.aes = FALSE,
            hjust = 1.02, vjust = 0, size = 3, colour = "grey30") +
  facet_wrap(~state, scales = "free_x") +
  scale_color_scs(name = NULL, labels = c("L0-MIO+stab" = "stable flag", "OLS" = "below threshold")) +
  labs(title = "Stability selection: bootstrap selection frequency per hospital",
       subtitle = sprintf("B = 100 subsamples; dashed line = pi_thr = %.1f (Meinshausen-Buhlmann threshold)", PI_THR),
       x = "hospital (ranked by selection frequency)", y = "selection frequency",
       caption = "Hospitals above the dashed line form the FDR-aware stable audit list. MB bound = expected number of false flags.") +
  theme_scs()
scs_save(pR2, "figs/R2_stability.png", w = 8.4, h = 4.4)

## ---- FIG R3: audit-budget curve -----------------------------------------------------------------
mk_aud <- function(res, state) {
  tb <- res$audit$table
  data.frame(lambda = tb$lambda, frac = tb$frac_of_max, state = state, sugg = res$audit$suggested_budget)
}
auddf <- rbind(mk_aud(RS, NY), mk_aud(RT, TX))
sg <- do.call(rbind, by(auddf, auddf$state, function(d) {
  s <- d$sugg[1]; data.frame(state = d$state[1], sugg = s, frac = approx(d$lambda, d$frac, s)$y) }))
pR3 <- ggplot(auddf, aes(lambda, frac)) +
  geom_hline(yintercept = 0.9, linetype = "dotted", colour = "grey55") +
  geom_line(colour = scs_cols["L0-MIO"], linewidth = 0.8) +
  geom_point(colour = scs_cols["L0-MIO"], size = 1) +
  geom_segment(data = sg, aes(x = sugg, xend = sugg, y = 0, yend = frac), colour = scs_emphasis, linewidth = 0.5) +
  geom_point(data = sg, aes(x = sugg, y = frac), colour = scs_emphasis, size = 2.4) +
  geom_text(data = sg, aes(x = sugg, y = frac, label = sprintf("lambda*=%d", sugg)),
            colour = scs_emphasis, vjust = -0.9, hjust = -0.05, size = 3) +
  facet_wrap(~state, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Audit budget: deviation captured vs number of hospitals reviewed",
       subtitle = "Dotted line = 90% of detectable risk-adjusted deviation; coral = suggested budget",
       x = expression("audit budget"~lambda~"(hospitals flagged)"),
       y = "fraction of max detectable deviation",
       caption = "The L0 solve returns, for each budget lambda, the hospitals that maximize detected risk-adjusted deviation (deliverable C).") +
  theme_scs()
scs_save(pR3, "figs/R3_audit_budget.png", w = 8.4, h = 4.4)

# low-res view copies (<2000px) for inline inspection; paper versions stay at 300 dpi above
ggsave("figs/R1_view.png", pR1, width = 8.4, height = 4.4, dpi = 110)
ggsave("figs/R2_view.png", pR2, width = 8.4, height = 4.4, dpi = 110)
ggsave("figs/R3_view.png", pR3, width = 8.4, height = 4.4, dpi = 110)
cat("saved: figs/R1_caterpillar_L0_vs_BLUP.png, figs/R2_stability.png, figs/R3_audit_budget.png (+ _view copies)\n")
