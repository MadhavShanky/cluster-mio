setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
d <- read.csv("sim/results_bench/bench_l0learn.csv")
cat("rows:", nrow(d), "\n")
en <- d[is.finite(d$gap_scs), ]
cat("enumerable rows:", nrow(en), "\n")
cat("max gap_scs (SCS vs brute, want ~0):", format(max(en$gap_scs), digits = 3), "\n")
cat("jacc_scs_brute min/max (want 1):", min(en$jacc_scs_brute), "/", max(en$jacc_scs_brute), "\n")
cat("SCS time K=300:", format(mean(d$t_scs[d$K == 300]), digits = 3),
    "s; K=1000:", format(mean(d$t_scs[d$K == 1000]), digits = 3), "s\n")
