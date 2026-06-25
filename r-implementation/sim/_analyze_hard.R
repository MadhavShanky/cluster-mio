setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
suppressMessages(library(dplyr))
fs <- list.files("sim/results_hard", "^cell_[0-9]+\\.csv$", full.names = TRUE)
d <- as.data.frame(bind_rows(lapply(fs, read.csv)))
for (c0 in c("gerr_l0_eff","gerr_lmm_eff","f1_stab","fdr_stab","f1_l0","fdr_l0","tpr_l0","icc"))
  if (is.null(d[[c0]])) d[[c0]] <- NA
cat("total reps:", nrow(d), " cells:", length(fs), " mean fail_lmm:", round(mean(d$fail_lmm, na.rm=TRUE),3), "\n\n")
ag <- function(sub, by, cols) {
  z <- aggregate(sub[, cols], by = sub[, by, drop=FALSE], FUN = function(v) mean(v, na.rm=TRUE))
  z[, cols] <- round(z[, cols], 3); z
}
cat("== BOUNDARY (weak signal): selection vs delta ==\n")
b <- d[d$regime=="sparse" & d$null_sd==0 & !d$size_effect & !d$drop_x & d$delta<=1, ]
print(ag(b, c("K","delta"), c("tpr_l0","f1_l0","fdr_l0","f1_stab","fdr_stab")), row.names=FALSE)

cat("\n== NULLSD (small nonzero nulls): spurious flagging ==\n")
ns <- d[d$null_sd>0, ]
print(ag(ns, c("K","null_sd"), c("tpr_l0","fdr_l0","f1_stab","fdr_stab")), row.names=FALSE)

cat("\n== DENSE_T (heavy-tailed dense): misspecification check, gerr L0 vs LMM ==\n")
dt <- d[d$regime=="dense_t", ]
print(ag(dt, c("K","icc"), c("gerr_l0","gerr_lmm")), row.names=FALSE)

cat("\n== SIZE/EFFECT COUPLING: detection ==\n")
sz <- d[d$size_effect==TRUE, ]
print(ag(sz, c("K","delta"), c("tpr_l0","f1_l0","fdr_l0","f1_stab")), row.names=FALSE)

cat("\n== MISSPEC (omitted cluster-correlated covariate) ==\n")
ms <- d[d$drop_x==TRUE, ]
print(ag(ms, c("K"), c("tpr_l0","fdr_l0","f1_stab","fdr_stab","gerr_l0","gerr_l0_eff")), row.names=FALSE)
