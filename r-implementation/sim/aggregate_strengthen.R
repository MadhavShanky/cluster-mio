#!/usr/bin/env Rscript
# Aggregate the 4 strengthen-campaign arrays (consistency / scheme / restart / nullpfer).
# Run on dengly: Rscript sim/aggregate_strengthen.R
# Writes tidy summaries to reviews/agg/ and prints headline numbers.
suppressWarnings(suppressMessages({}))
root <- "/home/dengly/mio-sim"
setwd(root)
outdir <- file.path(root, "reviews", "agg")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rbind_dir <- function(dir, pat = "\\.csv$") {
  fs <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(fs)) stop("no files in ", dir)
  do.call(rbind, lapply(fs, function(f) read.csv(f, stringsAsFactors = FALSE)))
}
# safe mean
m <- function(x) mean(x, na.rm = TRUE)

cat("################ SIM #1  CONSISTENCY / SCREENING ################\n")
con <- rbind_dir("sim/results_consistency")
cat(sprintf("rows=%d  reps/cell range: %s\n", nrow(con),
            paste(range(table(con$cid)), collapse="-")))
grpcols <- c("K","unequal","delta","s","mmin","snr","eps","n_restart","lambda_kind","lambda")
con$g <- apply(con[grpcols], 1, paste, collapse="|")
agg1 <- do.call(rbind, by(con, con$g, function(d) {
  data.frame(d[1, grpcols],
    nrep   = nrow(d),
    p_exact  = m(d$exact),
    p_screen = m(d$screen),
    f1     = m(d$f1),  fdr = m(d$fdr), tpr = m(d$tpr),
    n_false = m(d$n_false), n_sel = m(d$n_sel),
    jacc_eps = m(d$jacc_vs_epsref), jacc_nr = m(d$jacc_vs_nrref),
    signal = d$mmin[1]*d$delta[1]^2/log(d$K[1]),
    row.names = NULL)
}))
agg1 <- agg1[order(agg1$lambda_kind, agg1$K, agg1$unequal, agg1$delta, agg1$eps, agg1$n_restart), ]
write.csv(agg1, file.path(outdir, "agg_consistency.csv"), row.names = FALSE)
# headline: screening vs exact at lambda=s, eps=0, default restart
base1 <- subset(agg1, lambda_kind=="exact" & eps==min(agg1$eps) & n_restart==min(agg1$n_restart))
cat("-- lambda=s (exact), eps=min, restart=min : screening should approach 1 before exact --\n")
print(base1[, c("K","unequal","delta","signal","p_exact","p_screen","fdr","n_false")], digits=3)
cat("\n-- ridge robustness: Jaccard vs eps-ref (1 = ridge does not move support) --\n")
print(aggregate(jacc_eps ~ eps + lambda_kind, data=agg1, FUN=function(x) round(m(x),3)))

cat("\n################ SIM #2  OBS vs CLUSTER_STRAT SUBSAMPLING ################\n")
sc <- rbind_dir("sim/results_scheme")
gcol2 <- c("K","unequal","delta","s","scheme","eps")
sc$g <- apply(sc[gcol2], 1, paste, collapse="|")
agg2 <- do.call(rbind, by(sc, sc$g, function(d) {
  data.frame(d[1, gcol2], nrep=nrow(d),
    V_false=m(d$V_false), mb_ref=m(d$mb_ref), n_watch=m(d$n_watch),
    tpr=m(d$tpr), tpr_small=m(d$tpr_small), tpr_mid=m(d$tpr_mid), tpr_large=m(d$tpr_large),
    n_active_small=m(d$n_active_small), jacc_vs_obs=m(d$jacc_vs_obs), row.names=NULL)
}))
agg2 <- agg2[order(agg2$K, agg2$unequal, agg2$delta, agg2$scheme), ]
write.csv(agg2, file.path(outdir, "agg_scheme.csv"), row.names = FALSE)
cat("-- TPR by cluster-size band: does cluster_strat lift the SMALL-cluster cap? --\n")
print(agg2[, c("K","unequal","delta","scheme","tpr_small","tpr_mid","tpr_large","V_false","mb_ref")], digits=3)

cat("\n################ SIM #4  RESTART -> ROBUST CORE ################\n")
rs <- rbind_dir("sim/results_restart")
gcol4 <- c("K","unequal","delta","s","n_restart")
rs$g <- apply(rs[gcol4], 1, paste, collapse="|")
agg4 <- do.call(rbind, by(rs, rs$g, function(d) {
  data.frame(d[1, gcol4], nrep=nrow(d),
    n_core=m(d$n_core), n_watch=m(d$n_watch), tpr=m(d$tpr), V_false=m(d$V_false),
    max_null_freq=m(d$max_null_freq), jacc_watch=m(d$jacc_watch_vs_nrref), row.names=NULL)
}))
agg4 <- agg4[order(agg4$K, agg4$unequal, agg4$delta, agg4$n_restart), ]
write.csv(agg4, file.path(outdir, "agg_restart.csv"), row.names = FALSE)
cat("-- robust core (freq>=0.9) + watch-list stability vs n_restart --\n")
print(agg4[, c("K","unequal","delta","n_restart","n_core","n_watch","jacc_watch","max_null_freq")], digits=3)

cat("\n################ SIM #3  NULL PFER CALIBRATION ################\n")
nullsum <- list(); byclust <- list()
for (st in c("ny","tx")) {
  fs <- list.files(file.path("data", st), pattern = sprintf("^null_pfer_%s[0-9]+\\.csv$", st), full.names = TRUE)
  d <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors=FALSE))
  R <- nrow(d); seV <- sd(d$V, na.rm=TRUE)/sqrt(R)
  nullsum[[st]] <- data.frame(state=st, R=R, meanV=m(d$V), seV=seV,
    q=d$q[1], mb_nominal=d$mb_nominal[1],
    P_V_gt_q = m(d$V > d$q), P_V_ge_1 = m(d$V >= 1),
    P_V_gt_mb = m(d$V > d$mb_nominal[1]))
  # by-cluster tercile pooling
  bf <- list.files(file.path("data", st), pattern = sprintf("^null_pfer_bycluster_%s[0-9]+\\.csv$", st), full.names = TRUE)
  bc <- do.call(rbind, lapply(bf, read.csv, stringsAsFactors=FALSE))
  ter <- aggregate(false_flag_freq ~ size_tercile, data=bc, FUN=function(x) round(m(x),4))
  ter$state <- st
  byclust[[st]] <- ter
}
nullsum <- do.call(rbind, nullsum); byclust <- do.call(rbind, byclust)
write.csv(nullsum, file.path(outdir, "agg_nullpfer.csv"), row.names = FALSE)
write.csv(byclust, file.path(outdir, "agg_nullpfer_bytercile.csv"), row.names = FALSE)
cat("-- pooled R per state, realized E[V] vs nominal MB bound --\n")
print(nullsum, digits=4, row.names=FALSE)
cat("-- false-flag freq by cluster-size tercile --\n")
print(byclust, digits=4, row.names=FALSE)

cat("\nAll summaries written to ", outdir, "\n")
