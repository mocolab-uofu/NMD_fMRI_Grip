# Plot the completed 40-participant area-filtered analysis on Hazel.
# Run through the paired Slurm script and private display; read NMD_RUN_DIR/results.
# Update maps and comparisons in place; NMD_PLOT_RUN_DIR defaults to NMD_RUN_DIR.

# Read the completed formal run and fixed display surfaces
library(BayesfMRI)
library(ciftiTools)
library(png)
library(grid)
ciftiTools.setOption("wb_path", Sys.getenv("NMD_WB_PATH"))
stopifnot(packageVersion("BayesfMRI") == "0.15.0", packageVersion("ciftiTools") == "0.19.0", !rgl::rgl.useNULL())
run_path <- Sys.getenv("NMD_RUN_DIR")
stopifnot(nzchar(run_path), dir.exists(run_path), nzchar(Sys.getenv("SLURM_JOB_ID")))
run_path <- normalizePath(run_path)
run_info <- readRDS(file.path(run_path, "results/bglm2_runinfo.rds"))
stopifnot(run_info$status == "analysis_complete")
area_limits <- c(20, 30, 40)
model_files <- file.path(run_path, "results", paste0("bglm2_area_", area_limits, ".rds"))
model_md5 <- tools::md5sum(model_files)
models <- lapply(model_files, readRDS)
names(models) <- as.character(area_limits)
base <- models[["20"]]
contrast_names <- c("CTR > DEP", "DEP > CTR", "CTR average", "DEP average", "CTR and DEP average")
contrast_slugs <- c("CTR_gt_DEP", "DEP_gt_CTR", "CTR_average", "DEP_average", "CTR_DEP_average")
stopifnot(identical(base$contrast_estimate_xii$meta$cifti$names, contrast_names))
surface_path <- "/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/surface_models"
templates <- c("inflated", "very_inflated")
plot_run_path <- Sys.getenv("NMD_PLOT_RUN_DIR", run_path)
stopifnot(dir.exists(plot_run_path))
plot_run_path <- normalizePath(plot_run_path)
plot_path <- file.path(plot_run_path, "plots")
dir.create(plot_path, showWarnings = FALSE)
dir.create(file.path(plot_run_path, "analysis_records"), showWarnings = FALSE)
figures <- list()

# Use the same view, size, and contrast-specific color limits for all thresholds
for (cc in seq_along(contrast_names)) {
  zlim <- if (grepl("average", contrast_names[cc])) c(-1, 1) else c(-0.5, 0.5)
  for (template in templates) {
    surfL <- file.path(surface_path, paste0("S1200.L.", template, "_MSMAll.10k_fs_LR.surf.gii"))
    surfR <- file.path(surface_path, paste0("S1200.R.", template, "_MSMAll.10k_fs_LR.surf.gii"))
    stopifnot(file.exists(surfL), file.exists(surfR))
    left <- read_surf(surfL)
    right <- read_surf(surfR)
    stopifnot(nrow(left$vertices) == length(base$masks$In$cortexL), nrow(right$vertices) == length(base$masks$In$cortexR))

    # The unmasked coefficient estimate is shared by all area thresholds.
    folder <- file.path(plot_path, "surface/contrasts", contrast_slugs[cc])
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    fname <- file.path(folder, paste0(template, "_contrasts.png"))
    plot(base, stat = "contrasts", idx = cc, what = "surface", surfL = left, surfR = right, view = "both",
         title = "", zlim = zlim,
         width = 1200, height = 800, legend_embed = TRUE, fname = fname)
    if (rgl::cur3d() != 0) rgl::close3d()
    # Grid text supports line breaks that the OpenGL title renderer displays as boxes.
    view_comp(fname, ncol = 1, title = paste0(contrast_names[cc], "\nUnmasked coefficient estimate"),
              title_height = 0.1, title_fsize = 2, fname = fname)
    figures[[length(figures) + 1L]] <- data.frame(path = fname, kind = "surface", contrast = contrast_names[cc],
      template = template, area_limit = NA_real_, stat = "contrasts", area_filtered = FALSE)

    for (area_limit in area_limits) {
      res <- models[[as.character(area_limit)]]
      stopifnot(res$BayesGLM2_results$area.limit == area_limit)
      active <- c(res$activations_xii$data$cortex_left[, cc], res$activations_xii$data$cortex_right[, cc])
      empty <- !any(active == 1, na.rm = TRUE)
      title <- paste0(contrast_names[cc], "\nArea >= ", area_limit, " model units^2; alpha = 0.05",
                      if (empty) "; no surviving clusters" else "")
      folder <- file.path(plot_path, "surface", paste0("area_", area_limit), contrast_slugs[cc])
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      for (stat in c("activations", "activated_contrast")) {
        fname <- file.path(folder, paste0(template, "_", stat, ".png"))
        plot(res, stat = stat, idx = cc, what = "surface", surfL = left, surfR = right, view = "both",
             title = "", zlim = if (stat == "activations") c(0, 1) else zlim,
             width = 1200, height = 800, legend_embed = stat != "activations", fname = fname)
        if (rgl::cur3d() != 0) rgl::close3d()
        view_comp(fname, ncol = 1, title = title, title_height = 0.1, title_fsize = 2, fname = fname)
        figures[[length(figures) + 1L]] <- data.frame(path = fname, kind = "surface", contrast = contrast_names[cc],
          template = template, area_limit = area_limit, stat = stat, area_filtered = TRUE)
      }
    }

    # Assemble the three masked-effect maps without changing their colors or geometry.
    panel_files <- file.path(plot_path, "surface", paste0("area_", area_limits), contrast_slugs[cc], paste0(template, "_activated_contrast.png"))
    panels <- lapply(panel_files, png::readPNG)
    panel_dims <- vapply(panels, function(x) dim(x)[1:2], integer(2))
    stopifnot(length(unique(panel_dims[2, ])) == 1L)
    panel_height <- max(panel_dims[1, ]) # Empty activation maps may omit their color bar.
    folder <- file.path(plot_path, "comparisons", contrast_slugs[cc])
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    fname <- file.path(folder, paste0(template, "_area_20_30_40.png"))
    png(fname, width = 3 * panel_dims[2, 1], height = panel_height, bg = "white")
    grid.newpage()
    for (j in seq_along(panels)) {
      relative_height <- panel_dims[1, j] / panel_height
      grid.raster(panels[[j]], x = (j - 0.5) / 3, y = 1 - relative_height / 2, width = 1 / 3, height = relative_height, interpolate = FALSE)
    }
    dev.off()
    figures[[length(figures) + 1L]] <- data.frame(path = fname, kind = "comparison", contrast = contrast_names[cc],
      template = template, area_limit = NA_real_, stat = "activated_contrast_20_30_40", area_filtered = TRUE)
  }

  # Voxel inference is unchanged by cortical area filtering; export this supplement once.
  folder <- file.path(plot_path, "volume", contrast_slugs[cc])
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  for (stat in c("contrasts", "activations", "activated_contrast")) {
    empty <- stat != "contrasts" && !any(base$activations_xii$data$subcort[, cc] == 1, na.rm = TRUE)
    title <- paste0(contrast_names[cc], "\nVoxel supplement: no area filtering",
                    if (stat != "contrasts") "; alpha = 0.05" else "", if (empty) "; no activation" else "")
    fname <- file.path(folder, paste0(stat, ".png"))
    plot(base, stat = stat, idx = cc, what = "volume", title = title,
         zlim = if (stat == "activations") c(0, 1) else zlim, width = 1200, height = 800,
         legend_embed = stat != "activations", fname = fname)
    figures[[length(figures) + 1L]] <- data.frame(path = fname, kind = "volume", contrast = contrast_names[cc],
      template = NA_character_, area_limit = NA_real_, stat = stat, area_filtered = FALSE)
  }
}

# Record every main figure and any separate legends after checking PNG readability
figure_manifest <- do.call(rbind, figures)
stopifnot(nrow(figure_manifest) == 95L, !anyDuplicated(figure_manifest$path), all(file.exists(figure_manifest$path)))
extra_png <- setdiff(list.files(plot_path, pattern = "\\.png$", recursive = TRUE, full.names = TRUE), figure_manifest$path)
if (length(extra_png)) {
  figure_manifest <- rbind(figure_manifest, data.frame(path = extra_png, kind = "legend", contrast = NA_character_,
    template = NA_character_, area_limit = NA_real_, stat = "legend", area_filtered = NA))
}
dimensions <- matrix(NA_integer_, nrow(figure_manifest), 2, dimnames = list(NULL, c("height", "width")))
for (i in seq_len(nrow(figure_manifest))) {
  img <- png::readPNG(figure_manifest$path[i])
  stopifnot(all(dim(img)[1:2] > 0L), all(is.finite(img)))
  dimensions[i, ] <- dim(img)[1:2]
}
figure_manifest <- cbind(figure_manifest, dimensions)
figure_manifest$bytes <- file.info(figure_manifest$path)$size
figure_manifest$md5 <- unname(tools::md5sum(figure_manifest$path))
figure_manifest$plot_job_id <- Sys.getenv("SLURM_JOB_ID")
figure_manifest$code_commit <- Sys.getenv("NMD_CODE_COMMIT")
figure_manifest$analysis_job_id <- run_info$analysis_job_id
figure_manifest$analysis_code_commit <- run_info$code_commit
write.csv(figure_manifest, file.path(plot_path, "figure_manifest.csv"), row.names = FALSE)
stopifnot(identical(model_md5, tools::md5sum(model_files)))
write.csv(data.frame(path = model_files, md5 = unname(model_md5)), file.path(plot_run_path, "analysis_records/model_checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(plot_run_path, "analysis_records/plot_session_info.txt"))
writeLines(c("Plot export complete: 60 filtered surface maps, 10 coefficient references, 10 comparisons, 15 voxel supplements.",
  "Saved model checksums unchanged; no fitting, sampling or refiltering performed.",
  paste("Separate legends:", length(extra_png)), paste("Completed:", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(plot_run_path, "analysis_records/plot_checks.txt"))
# Keep each job's completion evidence while updating the same public figure paths.
job_record <- file.path(plot_run_path, "analysis_records", paste0("plots-job-", Sys.getenv("SLURM_JOB_ID")))
dir.create(job_record, showWarnings = FALSE)
record_files <- c(file.path(plot_path, "figure_manifest.csv"), file.path(plot_run_path, "analysis_records", c("model_checksums.csv", "plot_session_info.txt", "plot_checks.txt")))
stopifnot(all(file.copy(record_files, job_record, overwrite = TRUE)))
cat("\n* Plot export complete\nJob: ", Sys.getenv("SLURM_JOB_ID"), "\nCode: ", Sys.getenv("NMD_CODE_COMMIT"),
    "\nCompleted: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\nMain PNG files: 95\n", sep = "",
    file = file.path(plot_run_path, "run_record.org"), append = TRUE)
print(table(figure_manifest$kind))
