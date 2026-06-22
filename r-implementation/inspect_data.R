# inspect_data.R -- structural audit of the 3 open datasets for SCS/LMM fitness.
# For each: cluster unit K, within-cluster sizes, outcome type (continuous vs count), p covariates.
setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
line <- function(s) cat("\n==== ", s, " ====\n", sep = "")

## 1. NHSR LOS_model -- candidate hospital-profiling example -----------------------------------------
line("NHSR LOS_model")
los <- readRDS("data/NHSR_LOS_model.rds")
cat("dim:", paste(dim(los), collapse = " x "), "\n"); cat("cols:", paste(names(los), collapse = ", "), "\n")
print(utils::head(los, 3))
for (nm in names(los)) {
  v <- los[[nm]]
  if (is.factor(v) || is.character(v) || (is.numeric(v) && length(unique(v)) < 30))
    cat(sprintf("  %-14s: %d unique%s\n", nm, length(unique(v)),
                if (length(unique(v)) <= 12) paste0(" {", paste(sort(unique(v)), collapse=","), "}") else ""))
  else cat(sprintf("  %-14s: numeric range [%.2f, %.2f], mean %.2f\n", nm, min(v,na.rm=T), max(v,na.rm=T), mean(v,na.rm=T)))
}
# guess the org/cluster column
orgcol <- names(los)[which.max(sapply(los, function(c) if (is.numeric(c)) 0 else length(unique(c))))]
cat("likely cluster col:", orgcol, "\n")
tb <- table(los[[orgcol]]); cat("K =", length(tb), " | n_k range [", min(tb), ",", max(tb), "] median", median(tb), "\n")

## 2. SpatialEpi pennLC -- spatial county example ----------------------------------------------------
line("SpatialEpi pennLC")
pl <- readRDS("data/SpatialEpi_pennLC.rds")
cat("components:", paste(names(pl), collapse = ", "), "\n")
d2 <- pl$data; cat("data dim:", paste(dim(d2), collapse=" x "), " cols:", paste(names(d2), collapse=", "), "\n")
print(utils::head(d2, 3))
cat("K counties =", length(unique(d2$county)), " | strata per county =", nrow(d2)/length(unique(d2$county)), "\n")
cat("outcome 'cases': range [", min(d2$cases), ",", max(d2$cases), "] -- COUNT data\n")

## 3. TLC trial -- candidate longitudinal child example ----------------------------------------------
line("TLC trial")
raw <- readLines("data/tlc.txt", n = 5)
cat("first lines:\n"); cat(paste0("  ", raw), sep = "\n"); cat("\n")
tlc <- tryCatch(read.table("data/tlc.txt", header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(tlc)) {
  cat("dim:", paste(dim(tlc), collapse = " x "), "\n"); print(utils::head(tlc, 3))
  cat("interpretation: wide format, 1 row/child; lead at weeks 0,1,4,6 typical. K children =", nrow(tlc), "\n")
}
cat("\n==== done ====\n")
