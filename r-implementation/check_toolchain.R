cat("R", as.character(getRversion()), "\n")
for (p in c("Rcpp","RcppEigen","RcppArmadillo","pkgbuild")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-14s %s\n", p, ifelse(ok, as.character(packageVersion(p)), "MISSING")))
}
tb <- tryCatch(pkgbuild::has_build_tools(debug = FALSE), error = function(e) paste("ERR:", conditionMessage(e)))
cat("has_build_tools:", as.character(tb), "\n")
cat("make:", Sys.which("make"), "\n")
cat("g++:", Sys.which("g++"), "\n")
