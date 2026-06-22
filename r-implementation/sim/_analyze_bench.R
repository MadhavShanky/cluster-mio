setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
d <- read.csv("sim/results_bench/bench_l0learn.csv")
cat("rows:", nrow(d), " L0Learn non-NA obj:", sum(is.finite(d$obj_l0learn)), "/", nrow(d), "\n\n")

en <- d[is.finite(d$obj_brute), ]
cat("== Enumerable instances (brute available):", nrow(en), "rows ==\n")
cat("SCS objective gap to brute: max =", format(max(en$gap_scs), digits=3),
    " mean =", format(mean(en$gap_scs), digits=3), "\n")
cat("SCS support == brute (Jaccard=1):", round(100*mean(en$jacc_scs_brute==1),1), "%\n")
l0 <- en[is.finite(en$gap_l0learn), ]
cat("L0Learn rows w/ exact-size soln:", nrow(l0), "of", nrow(en), "\n")
if(nrow(l0)) {
  cat("L0Learn objective gap to brute: max =", format(max(l0$gap_l0learn), digits=3),
      " mean =", format(mean(l0$gap_l0learn), digits=3), "\n")
  cat("L0Learn support == brute:", round(100*mean(l0$jacc_l0_brute==1),1), "%\n")
  cat("SCS strictly better than L0Learn (lower obj):", sum(l0$obj_scs < l0$obj_l0learn - 1e-9),
      " | tie:", sum(abs(l0$obj_scs-l0$obj_l0learn)<1e-9), " | L0Learn better:", sum(l0$obj_l0learn < l0$obj_scs - 1e-9), "\n")
}
cat("\n== Timing (mean seconds, single fit) ==\n")
agg <- aggregate(cbind(t_scs, t_l0learn) ~ K, data=d, FUN=function(x) mean(x, na.rm=TRUE))
print(agg, digits=3)
cat("\nSCS mean time by K=300:", format(mean(d$t_scs[d$K==300]),digits=3),
    "K=1000:", format(mean(d$t_scs[d$K==1000]),digits=3), "s\n")
