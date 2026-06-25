setwd("C:/Users/tomch/AIProjects/research/MIO LMM/sim/results")
s <- read.csv("cell_0008.csv"); cat("SPARSE K=10:\n")
print(round(colMeans(s[, c("lambda","f1_l0","fdr_l0","f1_lmm","f1_glasso","f1_stab","mb_bound",
                           "gerr_l0","gerr_lmm","gerr_glasso","wmse_l0","wmse_lmm","cbmse_l0","cbmse_lmm")], na.rm = TRUE), 3))
x <- read.csv("cell_0112.csv"); cat("\nSPARSE_X K=30 (prediction win = cbmse_l0 < cbmse_lmm):\n")
print(round(colMeans(x[, c("lambda","f1_l0","f1_glasso","gerr_l0","gerr_lmm","wmse_l0","wmse_lmm","cbmse_l0","cbmse_lmm")], na.rm = TRUE), 3))
g <- read.csv("cell_0097.csv"); cat("\nGAUSSIAN K=30 (misspecification check: L0 should not beat LMM on gerr; the L0 target is misspecified here):\n")
print(round(colMeans(g[, c("icc_true","icc_lmm","gerr_l0","gerr_lmm","wmse_l0","wmse_lmm")], na.rm = TRUE), 3))
