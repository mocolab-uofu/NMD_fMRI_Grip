# Fit all 40 NMD participants on Hazel with BayesfMRI 0.15.0.
# Run through the paired Slurm script; read existing first-level models and CSV matrices.
# Save one posterior sample set and area-filtered group results in NMD_RUN_DIR.

# Select the pinned runtime and persistent output directory
library(BayesfMRI)
library(INLA)
library(ciftiTools)
library(peakRAM)
INLA::inla.setOption(num.threads = "1:1")
ciftiTools.setOption("wb_path", Sys.getenv("NMD_WB_PATH"))

upstream_commit <- "6f5b08fa498fe87af1f77a72acc63d28645ab198"
env_path <- "/rs1/researchers/a/afmejia/Software/ywu226/envs/nmd-bayesfmri-0.15.0-20261025"
stopifnot(getRversion() == "4.4.1", packageVersion("BayesfMRI") == "0.15.0")
stopifnot(packageVersion("INLA") == "24.12.11", packageVersion("ciftiTools") == "0.19.0")
stopifnot(normalizePath(find.package("BayesfMRI")) == file.path(env_path, "r-library/BayesfMRI"))
stopifnot(identical(readLines(file.path(env_path, "logs/BayesfMRI_upstream_commit.txt")), upstream_commit))

run_path <- Sys.getenv("NMD_RUN_DIR")
code_commit <- Sys.getenv("NMD_CODE_COMMIT")
job_id <- Sys.getenv("SLURM_JOB_ID")
stopifnot(nzchar(run_path), dir.exists(run_path), nzchar(code_commit), nzchar(job_id))
run_path <- normalizePath(run_path)
for (folder in c("results", "posterior_samples", "analysis_records", "tables")) dir.create(file.path(run_path, folder), showWarnings = FALSE)
output_path <- file.path(run_path, "results")
check_path <- file.path(run_path, "analysis_records")
record_path <- file.path(run_path, "run_record.org")
data_path <- "/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20260313_BayesfMRI2"
result_path <- file.path(data_path, "20251123_1stlevel")

# Keep the full sample and exact design-row order
rds_files <- list.files(result_path, pattern = "bglm\\.rds$", recursive = TRUE, full.names = TRUE)
rds_files <- rds_files[grepl("sub-NMD[0-9]+/02_BayesfMRI/bglm\\.rds$", rds_files)]
id_num <- as.integer(sub(".*sub-NMD([0-9]+).*", "\\1", rds_files))
ord <- order(id_num)
rds_files <- rds_files[ord]
id_num <- id_num[ord]
names(rds_files) <- sprintf("sub-NMD%03d", id_num)
expected_ids <- c(1, 2, 3, 4, 5, 6, 7, 9, 14, 19, 22, 23, 25, 26, 32, 33, 35, 38, 39, 40,
                  43, 45, 49, 50, 56, 57, 60, 62, 63, 64, 65, 67, 68, 71, 73, 76, 83, 84, 86, 87)
stopifnot(length(rds_files) == 40L, !anyDuplicated(id_num), all(id_num == expected_ids))
stopifnot(all(file.access(rds_files, 4) == 0), all(file.info(rds_files)$size > 0))

design_file <- file.path(data_path, "design_matrix.csv")
contrast_file <- file.path(data_path, "contrast_matrix.csv")
design_with_ids <- read.csv(design_file)
stopifnot(identical(names(design_with_ids), c("record_id", "group_CTR", "group_DEP", "age", "sex", "upgrade")))
stopifnot(identical(names(rds_files), paste0("sub-", design_with_ids$record_id)))
design_matrix <- design_with_ids[, -1]
stopifnot(all(is.finite(as.matrix(design_matrix))), qr(as.matrix(design_matrix))$rank == 5L)
stopifnot(all(design_matrix$group_CTR %in% 0:1), all(design_matrix$group_DEP %in% 0:1))
stopifnot(all(design_matrix$group_CTR + design_matrix$group_DEP == 1))
stopifnot(sum(design_matrix$group_CTR) == 22, sum(design_matrix$group_DEP) == 18)
contrast_matrix <- read.csv(contrast_file, header = FALSE)
expected_contrasts <- rbind(c(1, -1, 0, 0, 0), c(-1, 1, 0, 0, 0), c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0.5, 0.5, 0, 0, 0))
stopifnot(identical(dim(contrast_matrix), c(5L, 5L)), all(as.matrix(contrast_matrix) == expected_contrasts))
contrast_names <- c("CTR > DEP", "DEP > CTR", "CTR average", "DEP average", "CTR and DEP average")
excursion_type <- c(">", ">", "!=", "!=", "!=")
matrix_md5 <- tools::md5sum(c(design_file, contrast_file))
input_manifest <- data.frame(subject = names(rds_files), path = unname(rds_files), bytes = file.info(rds_files)$size,
                             modified = format(file.info(rds_files)$mtime, tz = "UTC", usetz = TRUE), row.names = NULL)
write.csv(input_manifest, file.path(check_path, "input_manifest.csv"), row.names = FALSE)
print(dim(design_matrix))
print(colSums(design_matrix[, c("group_CTR", "group_DEP")]))
print(contrast_matrix)

# Fit once at area 20; preserve the completed fit before deriving other thresholds
area_limits <- c(20, 30, 40) # Squared model-coordinate units, not assumed anatomical mm^2.
seed <- 20261025L
num_cores <- 4L
stopifnot(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) >= num_cores)
fit_path <- file.path(output_path, "bglm2_area_20.rds")
info_path <- file.path(output_path, "bglm2_runinfo.rds")
if (file.exists(fit_path)) {
  base <- readRDS(fit_path)
  run_info <- readRDS(info_path)
  stopifnot(identical(run_info$test_files, unname(rds_files)), identical(run_info$matrix_md5, matrix_md5))
  stopifnot(identical(run_info$input_manifest, input_manifest), identical(run_info$design_matrix, design_matrix))
  stopifnot(run_info$seed == seed, run_info$BayesfMRI_upstream_commit == upstream_commit)
  cat("Resuming filtering from the saved area-20 checkpoint.\n")
} else {
  # A partial sample run is evidence of an interrupted fit, not a resumable group object.
  stopifnot(length(list.files(file.path(run_path, "posterior_samples"), all.files = FALSE)) == 0L)
  run_info <- list(status = "fitting", code_commit = code_commit, analysis_job_id = job_id,
    BayesfMRI_version = as.character(packageVersion("BayesfMRI")), BayesfMRI_upstream_commit = upstream_commit,
    INLA_version = as.character(packageVersion("INLA")), ciftiTools_version = as.character(packageVersion("ciftiTools")),
    R_version = R.version.string, COMPUTER = Sys.info()[["nodename"]], result_path = result_path,
    test_subjects = names(rds_files), test_files = unname(rds_files), input_manifest = input_manifest,
    design_matrix = design_matrix, contrast_matrix = contrast_matrix, matrix_md5 = matrix_md5,
    contrast_names = contrast_names, excursion_type = excursion_type, gamma = 0, alpha = 0.05, alpha_nested = NULL,
    area_limits = area_limits, area_units = "squared model-coordinate units", seed = seed, rng_kind = RNGkind(),
    nsamp_theta = 50L, nsamp_beta = 100L, num_cores = num_cores, worker_threads = 1L,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), run_path = run_path)
  saveRDS(run_info, info_path)
  capture.output(sessionInfo(), file = file.path(check_path, "session_info.txt"))
  set.seed(seed)
  bglm2_profile <- peakRAM({
    base <- BayesGLM2(results = rds_files, design_matrix = design_matrix, contrast_matrix = contrast_matrix,
      excursion_type = excursion_type, contrast_names = contrast_names, gamma = 0, alpha = 0.05, alpha_nested = NULL,
      nsamp_theta = 50, nsamp_beta = 100, num_cores = num_cores, area.limit = 20, sample_dir = file.path(run_path, "posterior_samples"))
  })
  saveRDS(base, paste0(fit_path, ".tmp"))
  stopifnot(file.rename(paste0(fit_path, ".tmp"), fit_path))
  run_info$BayesGLM2_time_sec <- bglm2_profile$Elapsed_Time_sec[1]
  run_info$BayesGLM2_peak_RAM_MiB <- bglm2_profile$Peak_RAM_Used_MiB[1]
  run_info$BayesGLM2_total_RAM_MiB <- bglm2_profile$Total_RAM_Used_MiB[1]
  run_info$BayesGLM2_profile <- bglm2_profile
  run_info$status <- "fit_complete"
  saveRDS(run_info, info_path)
}

stopifnot(inherits(base, "BGLM2"), base$BayesGLM2_results$area.limit == 20)
stopifnot(base$BayesGLM2_results$nsamp_theta == 50, base$BayesGLM2_results$nsamp_beta == 100)
stopifnot(all(base$BayesGLM2_results$alpha == 0.05), is.null(base$BayesGLM2_results$alpha_nested))
stopifnot(identical(base$contrast_estimate_xii$meta$cifti$names, contrast_names))
manifest_path <- base$BayesGLM2_results$sample_manifest
manifest <- readRDS(manifest_path)
stopifnot(manifest$status == "complete", all(c("cortexL", "cortexR") %in% names(manifest$models)))
sample_files <- unique(unlist(lapply(manifest$models, `[[`, "sample.files"), use.names = FALSE))
stopifnot(length(sample_files) == 500L, all(file.exists(sample_files)))
sample_md5 <- tools::md5sum(sample_files)
if (!is.null(run_info$sample_md5)) stopifnot(identical(sample_md5, run_info$sample_md5))
write.csv(data.frame(path = names(sample_md5), md5 = unname(sample_md5)), file.path(check_path, "sample_checksums.csv"), row.names = FALSE)
mesh_rows <- list()
for (hem in c("cortexL", "cortexR")) {
  saved <- manifest$models[[hem]]
  stopifnot(saved$status == "complete", identical(saved$valid, base$masks$Mdat[[hem]]))
  stopifnot(nrow(saved$reference_mesh$loc) == length(saved$valid), all(is.finite(saved$reference_mesh$loc)))
  geometry <- BayesfMRI:::surface_cluster_geometry(saved$reference_mesh, saved$valid)
  mesh_rows[[hem]] <- data.frame(hemisphere = hem, reference_subject = names(rds_files)[1],
    vertices = length(saved$valid), valid_vertices = sum(saved$valid), faces = nrow(saved$reference_mesh$graph$tv),
    total_mesh_area = sum(geometry$areas), valid_mesh_area = sum(geometry$areas[saved$valid]),
    min_vertex_area = min(geometry$areas[saved$valid]), max_vertex_area = max(geometry$areas[saved$valid]),
    area_units = "squared model-coordinate units", source_manifest = manifest_path)
}
write.csv(do.call(rbind, mesh_rows), file.path(check_path, "mesh_summary.csv"), row.names = FALSE)

# Recompute each threshold from the same samples and original valid mask
for (area_limit in area_limits[-1]) {
  area_path <- file.path(output_path, paste0("bglm2_area_", area_limit, ".rds"))
  if (file.exists(area_path)) next
  res <- base
  derived <- manifest
  derived$status <- "filtering"
  derived$area.limit <- area_limit
  derived$source_manifest <- manifest_path
  derived$refilter <- list(function_name = "BayesfMRI:::excursions.mc.no.spurious", code_commit = code_commit,
    job_id = job_id, area.limit = area_limit, gamma = 0, alpha = 0.05, source_call = manifest$call,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  derived_path <- file.path(check_path, paste0("refilter_manifest_area_", area_limit, ".rds"))
  saveRDS(derived, derived_path)
  for (hem in c("cortexL", "cortexR")) {
    saved <- manifest$models[[hem]]
    geometry <- BayesfMRI:::surface_cluster_geometry(saved$reference_mesh, saved$valid)
    for (cc in seq_along(contrast_names)) {
      cat("Filtering", hem, contrast_names[cc], "at area", area_limit, "\n")
      files_cc <- vapply(saved$sample.files, `[`, character(1), cc)
      excur <- BayesfMRI:::excursions.mc.no.spurious(sample.files = files_cc, weights = saved$weights,
        mesh = saved$reference_mesh, area.limit = area_limit, u = 0, alpha = 0.05,
        type = excursion_type[cc], geometry = geometry)
      res$BayesGLM2_results$model_results[[hem]]$ppm[, cc] <- excur$F
      res$BayesGLM2_results$model_results[[hem]]$active[, cc] <- excur$E
      res$BayesGLM2_results$model_results[[hem]]$area_filter[[cc]] <- excur$area_filter
      res$BayesGLM2_results$model_results[[hem]]$area_filter[[cc]]$direction <- excur$direction
    }
    # Match BayesGLM2's mesh-to-CIFTI masking, including missing data inside maskIn.
    active_xii <- res$BayesGLM2_results$model_results[[hem]]$active
    active_xii[res$masks$In[[hem]] & !res$masks$Mdat[[hem]], ] <- NA
    active_xii <- active_xii[res$masks$In[[hem]], , drop = FALSE]
    part <- if (hem == "cortexL") "cortex_left" else "cortex_right"
    stopifnot(identical(dim(active_xii), dim(res$activations_xii$data[[part]])))
    res$activations_xii$data[[part]][, ] <- active_xii
    derived$models[[hem]]$area_filter <- res$BayesGLM2_results$model_results[[hem]]$area_filter
  }
  res$BayesGLM2_results$area.limit <- area_limit
  res$BayesGLM2_results$sample_manifest <- derived_path
  derived$status <- "complete"
  derived$refilter$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  saveRDS(derived, derived_path)
  saveRDS(res, paste0(area_path, ".tmp"))
  stopifnot(file.rename(paste0(area_path, ".tmp"), area_path))
}

# Check the completed formal results and summarize signed connected components
area_rows <- cluster_rows <- list()
for (area_limit in area_limits) {
  res <- readRDS(file.path(output_path, paste0("bglm2_area_", area_limit, ".rds")))
  stopifnot(res$BayesGLM2_results$area.limit == area_limit)
  stopifnot(identical(res$contrast_estimate_xii, base$contrast_estimate_xii))
  stopifnot(identical(res$BayesGLM2_results$model_results$subcort, base$BayesGLM2_results$model_results$subcort))
  stopifnot(identical(res$activations_xii$data$subcort, base$activations_xii$data$subcort))
  stopifnot(identical(res$BayesGLM2_results$sample_dir, base$BayesGLM2_results$sample_dir))
  for (hem in c("cortexL", "cortexR")) {
    saved <- manifest$models[[hem]]
    geometry <- BayesfMRI:::surface_cluster_geometry(saved$reference_mesh, saved$valid)
    model <- res$BayesGLM2_results$model_results[[hem]]
    stopifnot(identical(model$estimates, base$BayesGLM2_results$model_results[[hem]]$estimates))
    expected_xii <- model$active
    expected_xii[res$masks$In[[hem]] & !res$masks$Mdat[[hem]], ] <- NA
    expected_xii <- expected_xii[res$masks$In[[hem]], , drop = FALSE]
    part <- if (hem == "cortexL") "cortex_left" else "cortex_right"
    stopifnot(isTRUE(all.equal(unname(expected_xii), unname(res$activations_xii$data[[part]]))))
    for (cc in seq_along(contrast_names)) {
      report <- model$area_filter[[cc]]
      active <- model$active[, cc] == 1L
      stopifnot(report$area.limit == area_limit, !any(active & !saved$valid))
      stopifnot(identical(unname(report$sample.files), unname(vapply(saved$sample.files, `[`, character(1), cc))))
      stopifnot(isTRUE(all.equal(report$weights, saved$weights)), all(report$sample_counts == 100L))
      stopifnot(report$joint_probability > 0.95, all(model$ppm[active, cc] > 0.95))
      clusters <- BayesfMRI:::activation_clusters(active, geometry, report$direction)
      areas <- vapply(clusters, `[[`, numeric(1), "area")
      stopifnot(all(areas >= area_limit), identical(clusters, report$clusters))
      area_rows[[length(area_rows) + 1L]] <- data.frame(area_limit = area_limit, contrast = contrast_names[cc], hemisphere = hem,
        active_vertices = sum(active), clusters = length(clusters), total_active_area = sum(areas),
        min_cluster_area = if (length(areas)) min(areas) else NA_real_, joint_probability = report$joint_probability,
        area_units = "squared model-coordinate units")
      for (kk in seq_along(clusters)) {
        cluster_rows[[length(cluster_rows) + 1L]] <- data.frame(area_limit = area_limit, contrast = contrast_names[cc], hemisphere = hem,
          cluster = kk, direction = clusters[[kk]]$direction, vertices = length(clusters[[kk]]$vertices), area = areas[kk],
          area_units = "squared model-coordinate units")
      }
    }
  }
}
area_summary <- do.call(rbind, area_rows)
cluster_summary <- data.frame(area_limit = numeric(), contrast = character(), hemisphere = character(), cluster = integer(),
                              direction = integer(), vertices = integer(), area = numeric(), area_units = character())
if (length(cluster_rows)) cluster_summary <- do.call(rbind, cluster_rows)
write.csv(area_summary, file.path(run_path, "tables/area_summary.csv"), row.names = FALSE)
write.csv(cluster_summary, file.path(run_path, "tables/cluster_summary.csv"), row.names = FALSE)
stopifnot(identical(tools::md5sum(sample_files), sample_md5))
run_info$status <- "analysis_complete"
run_info$sample_manifest <- manifest_path
run_info$sample_md5 <- sample_md5
run_info$last_analysis_code_commit <- code_commit
run_info$last_analysis_job_id <- job_id
run_info$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
saveRDS(run_info, info_path)
writeLines(c("Formal analysis complete: 40 participants, 5 contrasts, areas 20/30/40.",
  "Shared posterior checksums unchanged; coefficients and voxel results identical across thresholds.",
  "CIFTI masks match mesh activations; all retained signed clusters meet the area limit.",
  paste("Completed:", run_info$completed_utc)), file.path(check_path, "run_checks.txt"))
cat("\n* Analysis complete\n", "Job: ", job_id, "\nCode: ", code_commit, "\nCompleted: ", run_info$completed_utc,
    "\nSamples: ", manifest_path, "\n", sep = "", file = record_path, append = TRUE)
print(area_summary)
cat("Formal area-filtered analysis complete.\n")
