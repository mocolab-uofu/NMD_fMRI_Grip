# Convert supplied first-level bglm.Rdata files to bglm.rds for group analysis.
# Run in R or with Rscript after setting data_path to the external analysis directory.
# Read and write within 20251123_1stlevel; existing bglm.rds files are skipped.

# Step1. Locate the supplied first-level models
data_path <- "/Volumes/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20260313_BayesfMRI2"
result_path <- file.path(data_path, "20251123_1stlevel")

files <- list.files(result_path, pattern = "bglm\\.Rdata$", recursive = TRUE, full.names = TRUE)
files <- files[grepl("sub-NMD[0-9]+/02_BayesfMRI/bglm\\.Rdata$", files)]
stopifnot(length(files) > 0)

# Step2. Save each bglm object in the format read by BayesGLM2
# Keep subjects separate in memory and retain the original .Rdata files.
for (f in files) {
  out_file <- file.path(dirname(f), "bglm.rds")
  if (file.exists(out_file)) {
    cat("Skipping (already exists):", out_file, "\n")
    next
  }

  cat("Converting:", f, "\n")
  tryCatch({
    e <- new.env()
    load(f, envir = e)
    if (!exists("bglm", envir = e)) {
      warning("No object named 'bglm' in: ", f)
    } else {
      saveRDS(e$bglm, out_file)
      cat("Saved:", out_file, "\n")
    }
    rm(e)
    gc()
  }, error = function(err) {
    cat("Skipping (bad file):", f, "\n")
    cat("Reason:", conditionMessage(err), "\n")
    gc()
  })
}

cat("All done.\n")
