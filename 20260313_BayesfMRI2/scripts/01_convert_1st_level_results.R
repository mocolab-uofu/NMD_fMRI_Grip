# Convert supplied first-level bglm.Rdata files to bglm.rds for group analysis.
# Run in R or with Rscript after setting data_path to the external analysis directory.
# Read and write within 20251123_1stlevel; existing bglm.rds files are skipped.

# Locate the supplied first-level models
data_path <- "/Volumes/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20260313_BayesfMRI2"
result_path <- file.path(data_path, "20251123_1stlevel")

rdata_files <- list.files(result_path, pattern = "bglm\\.Rdata$", recursive = TRUE, full.names = TRUE)
rdata_files <- rdata_files[grepl("sub-NMD[0-9]+/02_BayesfMRI/bglm\\.Rdata$", rdata_files)]
stopifnot(length(rdata_files) > 0)

# Save each bglm object in the format read by BayesGLM2
# Keep subjects separate in memory and retain the original .Rdata files.
for (rdata_file in rdata_files) {
  rds_file <- file.path(dirname(rdata_file), "bglm.rds")
  if (file.exists(rds_file)) {
    cat("Skipping (already exists):", rds_file, "\n")
    next
  }

  cat("Converting:", rdata_file, "\n")
  tryCatch({
    model_env <- new.env()
    load(rdata_file, envir = model_env)
    if (!exists("bglm", envir = model_env)) {
      warning("No object named 'bglm' in: ", rdata_file)
    } else {
      saveRDS(model_env$bglm, rds_file)
      cat("Saved:", rds_file, "\n")
    }
    rm(model_env)
    gc()
  }, error = function(err) {
    cat("Skipping (bad file):", rdata_file, "\n")
    cat("Reason:", conditionMessage(err), "\n")
    gc()
  })
}

cat("All done.\n")
