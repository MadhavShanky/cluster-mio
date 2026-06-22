# dl_data.R -- attempt to fetch the OPENLY-available example datasets to ./data.
# Credentialed sets (eICU, CAMP) cannot be auto-downloaded; this script only grabs the
# open ones and reports a status table. Run: Rscript dl_data.R  (see run notes in chat).
setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
dir.create("data", showWarnings = FALSE)
repo <- "https://cloud.r-project.org"
status <- list()
note   <- function(name, ok, msg = "") { status[[name]] <<- list(ok = ok, msg = msg) }

# 1. NHSRdatasets::LOS_model -- open, individual-level hospital length-of-stay (reproducible illustration).
tryCatch({
  if (!requireNamespace("NHSRdatasets", quietly = TRUE))
    install.packages("NHSRdatasets", repos = repo)
  e <- new.env(); data("LOS_model", package = "NHSRdatasets", envir = e)
  d <- e$LOS_model
  saveRDS(d, "data/NHSR_LOS_model.rds"); write.csv(d, "data/NHSR_LOS_model.csv", row.names = FALSE)
  note("NHSR_LOS_model", TRUE, sprintf("%d rows x %d cols", nrow(d), ncol(d)))
}, error = function(err) note("NHSR_LOS_model", FALSE, conditionMessage(err)))

# 2. SpatialEpi::pennLC -- open, Pennsylvania lung-cancer counts by county (spatial hot/cold-spot scenario).
tryCatch({
  if (!requireNamespace("SpatialEpi", quietly = TRUE))
    install.packages("SpatialEpi", repos = repo)
  e <- new.env(); data("pennLC", package = "SpatialEpi", envir = e)
  saveRDS(e$pennLC, "data/SpatialEpi_pennLC.rds")
  note("SpatialEpi_pennLC", TRUE, sprintf("%d county-strata rows", nrow(e$pennLC$data)))
}, error = function(err) note("SpatialEpi_pennLC", FALSE, conditionMessage(err)))

# 3. TLC trial -- open longitudinal classic (blood-lead trajectories), Fitzmaurice ALA data portal.
#    URL may move; tolerant of failure.
tryCatch({
  url <- "https://content.sph.harvard.edu/fitzmaur/ala2e/tlc-data.txt"
  dest <- "data/tlc.txt"
  download.file(url, dest, quiet = TRUE, mode = "wb")
  note("TLC", file.exists(dest) && file.info(dest)$size > 0, url)
}, error = function(err) note("TLC", FALSE, conditionMessage(err)))

cat("\n=== open-dataset download status ===\n")
for (nm in names(status)) cat(sprintf("  %-18s %s  %s\n", nm,
  if (status[[nm]]$ok) "OK " else "FAIL", status[[nm]]$msg))
cat("\nCredentialed (NOT auto-downloadable): eICU (PhysioNet), CAMP (NHLBI BioLINCC) -- see chat.\n")
