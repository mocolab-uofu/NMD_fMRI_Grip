# Organize the completed area-20/30/40 outputs once, on Hazel's shared Projects storage.
# Activate the existing R environment and set NMD_CODE_COMMIT to this checkout's revision.
# Update path metadata only; retain numerical results and the completed job-1571150 figures.

# Establish the exact source, destination and historical evidence
root <- "/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering"
old_run <- file.path(root, "full/job-1558105")
new_run <- file.path(root, "full/group_area_20_30_40")
replot <- "replots/title-fix-20261027"
sample_manifest <- "samples/BayesGLM2-3770226c94905a/manifest.rds"
code_commit <- Sys.getenv("NMD_CODE_COMMIT")
stopifnot(nzchar(code_commit), dir.exists(old_run), !file.exists(new_run), file.access(old_run, 2) == 0)
queued <- system2("squeue", c("-h", "-u", "ywu226", "-o", "%i"), stdout = TRUE)
stopifnot(is.null(attr(queued, "status")), length(queued) == 0L)
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
stage <- file.path(old_run, ".layout_migration")
stopifnot(!file.exists(stage)) # Inspect a previous partial migration before retrying.

old_rel <- c(paste0("results/bglm2_area_", c(20, 30, 40), ".rds"), "results/bglm2_runinfo.rds",
             sample_manifest, paste0("results/refilter_manifest_area_", c(30, 40), ".rds"))
new_rel <- c(old_rel[1:4], sub("^samples/", "posterior_samples/", old_rel[5]),
             sub("^results/", "analysis_records/", old_rel[6:7]))
objects <- lapply(file.path(old_run, old_rel), readRDS)
before_md5 <- unname(tools::md5sum(file.path(old_run, old_rel)))
plot_model_checks <- read.csv(file.path(old_run, replot, "checks/model_checksums.csv"), stringsAsFactors = FALSE)
stopifnot(identical(before_md5[1:3], plot_model_checks$md5))
samples <- read.csv(file.path(old_run, "checks/sample_checksums.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(samples) == 500L, !anyDuplicated(samples$path), all(file.exists(samples$path)))
stopifnot(identical(samples$path, names(objects[[4]]$sample_md5)), identical(samples$md5, unname(objects[[4]]$sample_md5)))
batch_files <- file.path(unique(dirname(samples$path)), "batch.rds")
stopifnot(length(batch_files) == 100L, all(file.exists(batch_files)))
batch_md5 <- unname(tools::md5sum(batch_files))
table_files <- list.files(file.path(old_run, "tables"), full.names = TRUE)
table_md5 <- unname(tools::md5sum(table_files))
figures <- read.csv(file.path(old_run, replot, "plots/figure_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(figures) == 125L, sum(figures$kind != "legend") == 95L, sum(figures$kind == "legend") == 30L)
stopifnot(all(figures$plot_job_id == 1571150), all(file.exists(figures$path)))
stopifnot(identical(unname(tools::md5sum(figures$path)), figures$md5))
actual_png <- list.files(file.path(old_run, replot, "plots"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
stopifnot(setequal(actual_png, figures$path))
old_png <- list.files(file.path(old_run, "plots"), pattern = "\\.png$", recursive = TRUE)
new_png <- list.files(file.path(old_run, replot, "plots"), pattern = "\\.png$", recursive = TRUE)
stopifnot(setequal(old_png, new_png))
main_record <- readLines(file.path(old_run, "run_record.org"))
plot_record <- readLines(file.path(old_run, replot, "run_record.org"))

# Replace active path strings and path-valued names, leaving original calls untouched
path_from <- c(file.path(old_run, old_rel[6:7]), file.path(old_run, "samples"), file.path(old_run, "checks"), old_run)
path_to <- c(file.path(new_run, new_rel[6:7]), file.path(new_run, "posterior_samples"), file.path(new_run, "analysis_records"), new_run)
replace_paths <- function(x, from = path_from, to = path_to) {
  for (i in seq_along(from)) {
    hit <- !is.na(x) & (x == from[i] | startsWith(x, paste0(from[i], "/")))
    x[hit] <- paste0(to[i], substring(x[hit], nchar(from[i]) + 1L))
  }
  x
}
relocate <- function(x, from = path_from, to = path_to) {
  if (is.language(x)) return(x) # Historical source calls describe the original execution.
  if (is.character(x)) x <- replace_paths(x, from, to)
  if (is.list(x)) for (i in seq_along(x)) x[i] <- list(relocate(x[[i]], from, to))
  if (!is.null(names(x))) names(x) <- replace_paths(names(x), from, to)
  x
}
updated <- lapply(objects, relocate)
for (i in seq_along(objects)) stopifnot(identical(relocate(updated[[i]], path_to, path_from), objects[[i]]))
for (i in 1:3) {
  stopifnot(identical(objects[[i]][1:5], updated[[i]][1:5]))
  for (hem in names(objects[[i]]$BayesGLM2_results$model_results)) {
    before <- objects[[i]]$BayesGLM2_results$model_results[[hem]]
    after <- updated[[i]]$BayesGLM2_results$model_results[[hem]]
    before$area_filter <- after$area_filter <- NULL
    stopifnot(identical(before, after))
  }
}

# Stage readable RDS replacements and retain temporary recovery copies on Projects
dir.create(stage)
for (i in seq_along(objects)) {
  original_file <- file.path(stage, "original", old_rel[i])
  staged_file <- file.path(stage, "updated", new_rel[i])
  dir.create(dirname(original_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(staged_file), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(old_run, old_rel[i]), original_file))
  saveRDS(updated[[i]], staged_file)
  stopifnot(identical(readRDS(staged_file), updated[[i]]))
}
writeLines(c(paste("Started:", started), paste("Code:", code_commit), paste("Source:", old_run), paste("Destination:", new_run)), file.path(stage, "migration.txt"))

# Rename directories on the same filesystem, then install the checked metadata
stopifnot(file.rename(old_run, new_run))
stage <- file.path(new_run, ".layout_migration")
stopifnot(file.rename(file.path(new_run, "checks"), file.path(new_run, "analysis_records")))
stopifnot(file.rename(file.path(new_run, "samples"), file.path(new_run, "posterior_samples")))
records <- file.path(new_run, "analysis_records")
for (i in seq_along(objects)) {
  if (i %in% 6:7) stopifnot(file.rename(file.path(new_run, old_rel[i]), file.path(new_run, new_rel[i])))
  stopifnot(file.rename(file.path(stage, "updated", new_rel[i]), file.path(new_run, new_rel[i])))
  stopifnot(identical(readRDS(file.path(new_run, new_rel[i])), updated[[i]]))
}

# Keep historic plotting records with their job IDs and expose the latest records
old_job <- file.path(records, "plots-job-1558106")
for (name in c("plot_checks.txt", "plot_session_info.txt")) {
  stopifnot(!file.exists(file.path(old_job, name)))
  stopifnot(file.rename(file.path(records, name), file.path(old_job, name)))
}
stopifnot(file.copy(file.path(new_run, "plots/figure_manifest.csv"), file.path(old_job, "figure_manifest.csv")))
latest_checks <- file.path(new_run, replot, "checks")
latest_job <- file.path(records, "plots-job-1571150")
stopifnot(!file.exists(latest_job), file.rename(file.path(latest_checks, "plots-job-1571150"), latest_job))
stopifnot(file.copy(file.path(new_run, replot, "plots/figure_manifest.csv"), file.path(latest_job, "figure_manifest.csv")))
for (name in c("model_checksums.csv", "models_before_submission.md5", "plot_checks.txt", "plot_session_info.txt")) {
  stopifnot(file.copy(file.path(latest_checks, name), file.path(latest_job, name)))
}
for (name in c("plot_checks.txt", "plot_session_info.txt")) stopifnot(file.copy(file.path(latest_checks, name), file.path(records, name)))

# Promote the corrected images; hold the superseded set until all final checks pass
stopifnot(file.rename(file.path(new_run, "plots"), file.path(stage, "superseded_plots")))
stopifnot(file.rename(file.path(new_run, replot, "plots"), file.path(new_run, "plots")))
figures$path <- replace_paths(figures$path, file.path(old_run, replot, "plots"), file.path(new_run, "plots"))
write.csv(figures, file.path(new_run, "plots/figure_manifest.csv"), row.names = FALSE)
samples$path <- replace_paths(samples$path)
write.csv(samples, file.path(records, "sample_checksums.csv"), row.names = FALSE)
mesh <- read.csv(file.path(records, "mesh_summary.csv"), stringsAsFactors = FALSE)
mesh$source_manifest <- replace_paths(mesh$source_manifest)
write.csv(mesh, file.path(records, "mesh_summary.csv"), row.names = FALSE)

# Check every posterior file and all corrected figures at their final locations
stopifnot(identical(unname(tools::md5sum(samples$path)), samples$md5))
stopifnot(identical(unname(tools::md5sum(replace_paths(batch_files))), batch_md5))
stopifnot(identical(unname(tools::md5sum(replace_paths(table_files))), table_md5))
stopifnot(identical(unname(tools::md5sum(figures$path)), figures$md5))
stopifnot(setequal(list.files(file.path(new_run, "plots"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE), figures$path))
for (i in seq_len(nrow(figures))) {
  size <- dim(png::readPNG(figures$path[i]))[1:2]
  stopifnot(identical(as.integer(size), as.integer(c(figures$height[i], figures$width[i]))))
}
for (i in 1:3) {
  fit <- readRDS(file.path(new_run, new_rel[i]))$BayesGLM2_results
  stopifnot(dir.exists(fit$sample_dir), file.exists(fit$sample_manifest))
  manifest <- readRDS(fit$sample_manifest)
  stopifnot(manifest$status == "complete", manifest$area.limit == fit$area.limit)
  if (i > 1L) stopifnot(file.exists(manifest$source_manifest))
  for (hem in c("cortexL", "cortexR")) {
    stopifnot(dir.exists(manifest$models[[hem]]$directory))
    stopifnot(all(file.exists(unlist(manifest$models[[hem]]$sample.files, use.names = FALSE))))
    for (cc in seq_along(fit$model_results[[hem]]$area_filter)) {
      stopifnot(identical(unname(fit$model_results[[hem]]$area_filter[[cc]]$sample.files),
                          unname(vapply(manifest$models[[hem]]$sample.files, `[`, character(1), cc))))
    }
  }
}
after_md5 <- unname(tools::md5sum(file.path(new_run, new_rel)))
write.csv(data.frame(old_path = file.path(old_run, old_rel), path = file.path(new_run, new_rel),
                    before_md5 = before_md5, after_md5 = after_md5), file.path(records, "migration_checksums.csv"), row.names = FALSE)
write.csv(data.frame(path = file.path(new_run, new_rel[1:3]), md5 = after_md5[1:3]), file.path(records, "model_checksums.csv"), row.names = FALSE)
write.csv(data.frame(old_path = c(path_from, file.path(old_run, replot, "plots")),
                    path = c(path_to, file.path(new_run, "plots"))), file.path(records, "directory_mapping.csv"), row.names = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
stopifnot(length(script_path) == 1L, file.copy(script_path, file.path(records, "06_organize_area_outputs.R")))

# Consolidate the current entry point and the unchanged historical execution record
header <- c("#+TITLE: NMD Area-Filtered Group Analysis", "* Current result location",
  paste("Run directory:", new_run), paste("Mounted directory:", sub("/rs1/researchers/a/afmejia", "/Volumes/afmejia", new_run, fixed = TRUE)),
  "Analysis job: 1558105; current plot job: 1571150 (corrected two-line titles).",
  "Review plots/comparisons/ and tables/; provenance is in analysis_records/.",
  "The sections below preserve historical commands, paths, checksums and deliveries.",
  "Use the current directory above and the migration mapping instead of those old paths.", "")
combined <- c(header, "* Historical analysis and initial plotting record", sub("^(\\*+)", "*\\1", main_record[-1]),
  "", "* Historical title-correction record", sub("^(\\*+)", "*\\1", plot_record[-1]),
  "", "* Directory organization", paste("Started UTC:", started), paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Organization code:", code_commit), paste("Previous run directory:", old_run), paste("Current run directory:", new_run),
  "checks/ -> analysis_records/; samples/ -> posterior_samples/.",
  "results/refilter_manifest_area_{30,40}.rds -> analysis_records/ with matching model references.",
  "replots/title-fix-20261027/plots/ -> plots/; only job-1571150 figures retained.",
  "All numerical model contents, masks and filtering diagnostics are unchanged except path strings.",
  "All 500 posterior matrices, 100 batch records, two summary tables and 125 corrected PNGs passed byte-checksum checks.",
  "RDS metadata checksums changed only through path relocation; see analysis_records/migration_checksums.csv.",
  "Historical plot-job checksum files retain their original paths and values; current model checksums are in analysis_records/model_checksums.csv.",
  "No model fitting, sampling, area filtering or rendering was run. No results were downloaded or pushed to GitHub.",
  "Future plot edits update this plots/ directory unless a separate destination is explicitly requested.")
writeLines(combined, file.path(new_run, "run_record.org"))
writeLines(c("Directory migration complete.", "Three models: numerical contents unchanged; all active sample and manifest references resolve.",
  "500 posterior files and 100 batch records: original byte checksums retained.",
  "125 corrected PNGs (95 main, 30 legends): checksums, dimensions and manifest verified.",
  "Two summary tables unchanged; original analysis and plot provenance retained.", paste("Code:", code_commit),
  paste("Completed:", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(records, "migration_checks.txt"))

# Remove superseded figures and temporary recovery files only after the final checks
stopifnot(unlink(file.path(new_run, "replots"), recursive = TRUE) == 0L)
stopifnot(unlink(stage, recursive = TRUE) == 0L)
stopifnot(!dir.exists(old_run), !dir.exists(stage), !dir.exists(file.path(new_run, "replots")))
cat("Migration complete:", new_run, "\n")
print(table(figures$kind))
