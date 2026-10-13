# aggregate_planted.R -- pool data/{ny,tx}/planted_pfer_<state><shard>.csv (planted_pfer_calibration.R).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
}))
for (st in c("ny", "tx")) {
  fs <- list.files(file.path("data", st), pattern = sprintf("^planted_pfer_%s\\d+\\.csv$", st), full.names = TRUE)
  if (!length(fs)) next
  d <- do.call(rbind, lapply(fs, read.csv)); R <- nrow(d)
  cat(sprintf("%s: %d shards, R=%d | planted=%d | mean V=%.3f (MC SE %.3f), max V=%d | P(V>ref)=%.3f | ref=%.2f | mean TP=%.2f (%.1f%% of planted), P(all planted flagged)=%.3f | mean q=%.2f\n",
              toupper(st), length(fs), R, d$n_planted[1], mean(d$V), sd(d$V) / sqrt(R), max(d$V),
              mean(d$V > d$ref), mean(d$ref), mean(d$TP), 100 * mean(d$TP) / d$n_planted[1],
              mean(d$TP == d$n_planted), mean(d$q)))
}
