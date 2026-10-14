# Fit the NMD group model with covariates and five contrasts using BayesGLM2.
# Run in R or with Rscript after setting the external data and Workbench paths.
# Read first-level models and CSV matrices; save the group model and run information externally.

# Step1. Set the environment and external paths
library(ciftiTools)
library(INLA)
library(BayesfMRI)
library(peakRAM)

COMPUTER <- Sys.info()[["nodename"]]
wb_path <- "/path/to/workbench"  # Replace with the local Connectome Workbench installation.
ciftiTools.setOption("wb_path", wb_path)

data_path <- "/Volumes/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20260313_BayesfMRI2"
result_path <- file.path(data_path, "20251123_1stlevel")
output_path <- file.path(data_path, "results")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Step2. Match first-level models to the design matrix
rds_files <- list.files(result_path, pattern = "bglm\\.rds$", recursive = TRUE, full.names = TRUE)
rds_files <- rds_files[grepl("sub-NMD[0-9]+/02_BayesfMRI/bglm\\.rds$", rds_files)]
id_num <- as.numeric(sub(".*sub-NMD([0-9]+).*", "\\1", rds_files))
ord <- order(id_num)
rds_files <- rds_files[ord]
id_num <- id_num[ord]
names(rds_files) <- sprintf("sub-NMD%03d", id_num)

design_matrix <- read.csv(file.path(data_path, "design_matrix.csv"))
# Check IDs before removing them; never silently reorder or drop design rows.
stopifnot(length(rds_files) > 0, !anyDuplicated(id_num))
stopifnot(identical(names(design_matrix), c("record_id", "group_CTR", "group_DEP", "age", "sex", "upgrade")))
stopifnot(identical(names(rds_files), paste0("sub-", design_matrix$record_id)))
design_matrix <- design_matrix[, -1]
contrast_matrix <- read.csv(file.path(data_path, "contrast_matrix.csv"), header = FALSE)
stopifnot(nrow(contrast_matrix) == 5, ncol(contrast_matrix) == ncol(design_matrix))

cat("All information loaded successfully.\n")

# Step3. Fit and save the group model
# Preserve the original model settings. alpha_nested selects nested activation mode;
# the package ignores alpha in this mode. Sampling counts retain the package defaults.
bglm2_profile <- peakRAM({
  bglm2 <- BayesGLM2(
    results = rds_files, design_matrix = design_matrix, contrast_matrix = contrast_matrix,
    excursion_type = c(">", ">", "!=", "!=", "!="),
    contrast_names = c("CTR > DEP", "DEP > CTR", "CTR average", "DEP average", "CTR and DEP average"),
    gamma = 0, alpha = 0.05, alpha_nested = c(0.05, 0.01)
  )
})

# Checkpoint: this fit previously took about 13.9 hours; plotting reads this saved object.
saveRDS(bglm2, file.path(output_path, "bglm2.rds"))

# Step4. Save software versions and resource measurements
# Retain the historical test_subjects/test_files keys; here they describe the full sample.
run_info <- list(
  BayesfMRI_version = as.character(packageVersion("BayesfMRI")),
  INLA_version = as.character(packageVersion("INLA")),
  ciftiTools_version = as.character(packageVersion("ciftiTools")),
  R_version = R.version.string,
  COMPUTER = COMPUTER,
  result_path = result_path,
  test_subjects = names(rds_files),
  test_files = unname(rds_files),
  BayesGLM2_time_sec = bglm2_profile$Elapsed_Time_sec[1],
  BayesGLM2_peak_RAM_MiB = bglm2_profile$Peak_RAM_Used_MiB[1],
  BayesGLM2_total_RAM_MiB = bglm2_profile$Total_RAM_Used_MiB[1],
  BayesGLM2_profile = bglm2_profile
)
saveRDS(run_info, file.path(output_path, "bglm2_runinfo.rds"))

cat("Done.\n")
