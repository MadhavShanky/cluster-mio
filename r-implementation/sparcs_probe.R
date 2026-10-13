# sparcs_probe.R -- verify NY SPARCS De-Identified Inpatient Discharges via the open Socrata API.
# Pulls a sample (no login/DUA), confirms: K hospitals, continuous LOS/charges, case-mix covariates.
setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
dir.create("data/sparcs", showWarnings = FALSE)
url  <- "https://health.data.ny.gov/resource/sf4k-39ay.csv?$limit=40000"
dest <- "data/sparcs/sparcs_sample.csv"
ok <- tryCatch({ download.file(url, dest, quiet = TRUE, mode = "wb"); TRUE },
               error = function(e) { cat("DOWNLOAD FAILED:", conditionMessage(e), "\n"); FALSE })
if (ok) {
  d <- read.csv(dest, stringsAsFactors = FALSE, check.names = TRUE)
  cat("sample dim:", paste(dim(d), collapse = " x "), "\n\ncolumns:\n")
  print(names(d))
  # locate facility + LOS + charges columns by name pattern
  fc <- grep("facility", names(d), ignore.case = TRUE, value = TRUE)
  lc <- grep("length|los",  names(d), ignore.case = TRUE, value = TRUE)
  cc <- grep("charge|cost",  names(d), ignore.case = TRUE, value = TRUE)
  cat("\nfacility cols:", paste(fc, collapse = ", "),
      "\nLOS cols:",      paste(lc, collapse = ", "),
      "\ncharge cols:",   paste(cc, collapse = ", "), "\n")
  idcol <- grep("permanent_facility_id|facility_id", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(idcol) && length(fc)) idcol <- fc[1]
  if (!is.na(idcol)) {
    K <- length(unique(d[[idcol]]))
    cat(sprintf("\nIN SAMPLE: K = %d distinct '%s' | rows/facility median %.0f\n",
                K, idcol, median(table(d[[idcol]]))))
  }
  losc <- lc[1]
  if (!is.na(losc)) {
    los <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", d[[losc]])))
    cat(sprintf("LOS ('%s'): range [%.0f, %.0f], median %.0f, mean %.1f -- continuous\n",
                losc, min(los, na.rm=T), max(los, na.rm=T), median(los, na.rm=T), mean(los, na.rm=T)))
  }
  cat("\nfirst row:\n"); print(t(d[1, , drop = FALSE]))
}
