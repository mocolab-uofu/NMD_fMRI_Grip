# Export interactive cortical effects and subcortical volumes from the completed Hazel models.
# Run through the paired Slurm script. Inspect CTR > DEP first, then submit the remaining contrasts.
# The first job stages its maps under analysis_records; the second completes the exports/ directory.

library(BayesfMRI)
library(ciftiTools)
ciftiTools.setOption("wb_path", Sys.getenv("NMD_WB_PATH"))
stopifnot(packageVersion("BayesfMRI") == "0.15.0", packageVersion("ciftiTools") == "0.19.0", rmarkdown::pandoc_available())

# Select the same maps, names and geometry as the signed-effect PNGs and voxel supplements
project <- "/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah"
run <- file.path(project, "20261025_BayesfMRI_AreaFiltering/full/group_area_20_30_40")
job <- Sys.getenv("SLURM_JOB_ID")
phase <- Sys.getenv("NMD_EXPORT_PHASE")
first_job <- Sys.getenv("NMD_FIRST_EXPORT_JOB")
record <- file.path(run, "analysis_records", paste0("exports-job-", job))
first_record <- file.path(run, "analysis_records", paste0("exports-job-", first_job))
stage <- file.path(first_record, "maps")
stopifnot(phase %in% c("first", "remaining"), dir.exists(record), !file.exists(file.path(run, "exports")))
if (phase == "first") stopifnot(dir.create(stage))

contrasts <- c("CTR > DEP", "DEP > CTR", "CTR average", "DEP average", "CTR and DEP average")
areas <- c(20, 30, 40)
templates <- c("inflated", "very_inflated")
figures <- read.csv(file.path(run, "plots/figure_manifest.csv"))
keep <- (figures$kind == "surface" & figures$stat %in% c("contrasts", "activated_contrast")) |
        (figures$kind == "volume" & figures$stat %in% c("contrasts", "activations"))
maps <- figures[keep, c("path", "kind", "contrast", "template", "area_limit", "stat")]
maps$path <- sub("^(.*/)?plots/", "", maps$path)
maps$path <- paste0(sub("\\.png$", "", maps$path), ifelse(maps$kind == "volume", ".nii", ".html"))
stopifnot(!anyDuplicated(maps$path), setequal(maps$contrast, contrasts), all(table(maps$contrast, maps$kind)[, "surface"] == 8L),
          all(table(maps$contrast, maps$kind)[, "volume"] == 2L))
selected <- if (phase == "first") maps$contrast == contrasts[1] else maps$contrast != contrasts[1]
if (phase == "remaining") {
  previous <- read.csv(file.path(first_record, "first_maps.csv"))
  stopifnot(identical(previous$path, maps$path[!selected]), identical(previous$md5, unname(tools::md5sum(file.path(stage, previous$path)))),
            identical(readLines(file.path(first_record, "code_commit.txt")), Sys.getenv("NMD_CODE_COMMIT")))
}

# Match the three saved fits to the original checksums and common reference coefficients
files <- file.path(run, "results", paste0("bglm2_area_", areas, ".rds"))
checksums <- read.csv(file.path(run, "analysis_records/model_checksums.csv"))
stopifnot(identical(unname(tools::md5sum(files)), checksums$md5[match(files, checksums$path)]))
write.csv(checksums, file.path(record, "input_checksums.csv"), row.names = FALSE)
models <- setNames(lapply(files, readRDS), areas)
reference <- models[["20"]]
beta <- reference$contrast_estimate_xii
active <- reference$activations_xii
stopifnot(identical(beta$meta$cifti$names, contrasts))
for (limit in areas) {
  model <- models[[as.character(limit)]]
  stopifnot(model$BayesGLM2_results$area.limit == limit, identical(model$contrast_estimate_xii, beta),
            identical(model$activations_xii$meta$cifti$names, contrasts), identical(model$activations_xii$data$subcort, active$data$subcort),
            identical(model$activations_xii$meta$subcort, active$meta$subcort))
}

# Restore the voxel-to-world transform only after matching the full mask and ordered labels
inputs <- read.csv(file.path(run, "analysis_records/input_manifest.csv"))
geometry_file <- file.path(dirname(dirname(inputs$path[1])), "03_Activation/activation.Rdata")
first_level <- new.env()
load(geometry_file, envir = first_level)
geometry <- first_level$act$activations_xii$single_sess$meta$subcort
for (spatial in list(beta$meta$subcort, active$meta$subcort)) {
  stopifnot(identical(geometry$mask, spatial$mask), identical(geometry$labels, spatial$labels))
}
affine <- geometry$trans_mat
spacing <- sqrt(colSums(affine[1:3, 1:3]^2))
stopifnot(identical(dim(geometry$mask), c(91L, 109L, 91L)), all(spacing == 2), identical(geometry$trans_units, "mm"), det(affine) != 0)
writeLines(c(geometry_file, capture.output(print(affine)), "Spatial units: mm"), file.path(record, "volume_geometry.txt"))

for (i in which(selected & maps$kind == "volume")) {
  cc <- match(maps$contrast[i], contrasts)
  values <- if (maps$stat[i] == "contrasts") beta$data$subcort[, cc] else active$data$subcort[, cc]
  if (maps$stat[i] == "activations") stopifnot(all(na.omit(values) %in% 0:1))
  values[is.na(values)] <- NaN
  volume <- array(NaN, dim = dim(geometry$mask))
  volume[geometry$mask] <- values
  volume <- RNifti::`pixdim<-`(volume, spacing)
  volume <- RNifti::`sform<-`(volume, structure(affine, code = 2L))
  volume <- RNifti::`pixunits<-`(volume, "mm")
  path <- file.path(stage, maps$path[i])
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(volume, path, datatype = "double")
  saved <- RNifti::readNifti(path)
  transform <- RNifti::xform(saved, useQuaternionFirst = FALSE)
  stopifnot(identical(dim(saved), dim(geometry$mask)), identical(as.numeric(transform), as.numeric(affine)),
            all(RNifti::pixdim(saved)[1:3] == spacing), as.character(RNifti::pixunits(saved))[1] == "mm",
            all(is.na(saved[!geometry$mask])), identical(as.numeric(saved[geometry$mask]), values))
}

# Validate each displayed vertex against the full model arrays before passing it to the script 05 renderer
surfaces <- list()
for (template in templates) {
  surface_files <- file.path(project, "surface_models", paste0("S1200.", c("L", "R"), ".", template, "_MSMAll.10k_fs_LR.surf.gii"))
  surfaces[[template]] <- lapply(surface_files, read_surf)
  stopifnot(identical(vapply(surfaces[[template]], function(s) nrow(s$vertices), integer(1)), unname(lengths(beta$meta$cortex$medial_wall_mask))))
}
for (i in which(selected & maps$kind == "surface")) {
  cc <- match(maps$contrast[i], contrasts)
  filtered <- maps$stat[i] == "activated_contrast"
  limit <- maps$area_limit[i]
  model <- if (filtered) models[[as.character(limit)]] else reference
  surface <- surfaces[[maps$template[i]]]
  display <- beta
  for (hem in c("cortexL", "cortexR")) {
    side <- if (hem == "cortexL") "left" else "right"
    part <- paste0("cortex_", side)
    inside <- model$masks$In[[hem]]
    valid <- model$masks$Mdat[[hem]]
    fit <- model$BayesGLM2_results$model_results[[hem]]
    stopifnot(identical(inside, beta$meta$cortex$medial_wall_mask[[side]]))
    expected <- fit$estimates[, cc]
    expected[!valid] <- NA_real_
    if (filtered) {
      mask <- model$activations_xii$data[[part]][, cc]
      retained <- !is.na(mask) & mask == 1
      fitted_set <- valid & fit$active[, cc] %in% 1
      stopifnot(identical(which(retained), which(fitted_set[inside])))
      display$data[[part]][!retained, cc] <- NA_real_
      expected[!fitted_set] <- NA_real_
    }
    stopifnot(identical(as.numeric(display$data[[part]][, cc]), as.numeric(expected[inside])))
  }
  label <- if (filtered) paste0("Masked effect | area >= ", limit, " model units^2 | alpha = 0.05") else "Unmasked effect"
  title <- paste(contrasts[cc], gsub("_", " ", maps$template[i]), label, sep = " | ")
  scale <- if (cc <= 2) c(-0.5, 0.5) else c(-1, 1)
  path <- file.path(stage, maps$path[i])
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Pandoc may retain a dependency directory after embedding its contents. Keep it in job scratch.
  widget_file <- tempfile(pattern = "surface-", fileext = ".html")
  view_xifti(display, idx = cc, what = "surface", surfL = surface[[1]], surfR = surface[[2]], view = "both",
             zlim = scale, title = title, width = 1200, height = 800, widget = TRUE, legend_embed = TRUE, fname = widget_file)
  if (rgl::cur3d() != 0) rgl::close3d()
  stopifnot(file.copy(widget_file, path, overwrite = TRUE))
  page <- xml2::read_html(path, options = c("HUGE", "NONET", "RECOVER"))
  widgets <- xml2::xml_find_all(page, "//*[contains(@class, 'rglWebGL') and @id]")
  scripts <- xml2::xml_find_all(page, "//script[@type='application/json' and @data-for]")
  resources <- xml2::xml_find_all(page, "//script[@src]/@src | //link[@rel='stylesheet']/@href | //img[@src]/@src")
  stopifnot(length(widgets) > 0L, all(xml2::xml_attr(widgets, "id") %in% xml2::xml_attr(scripts, "data-for")),
            all(startsWith(xml2::xml_text(resources), "data:")))
  cat(maps$path[i], "\n")
}

# Retain the first maps for visual inspection, then deliver the full validated set on Projects
capture.output(sessionInfo(), file = file.path(record, "session_info.txt"))
manifest <- if (phase == "first") maps[selected, ] else maps
stopifnot(setequal(list.files(stage, recursive = TRUE), manifest$path))
manifest$bytes <- file.info(file.path(stage, manifest$path))$size
manifest$md5 <- unname(tools::md5sum(file.path(stage, manifest$path)))
stopifnot(all(manifest$bytes > 0), !anyNA(manifest$md5))
if (phase == "first") {
  write.csv(manifest, file.path(record, "first_maps.csv"), row.names = FALSE)
} else {
  write.csv(manifest, file.path(stage, "map_manifest.csv"), row.names = FALSE)
  stopifnot(file.rename(stage, file.path(run, "exports")))
  cat("\n* Cortical HTML and subcortical NIfTI exports\nJobs: ", first_job, " and ", job, "\nCode: ", Sys.getenv("NMD_CODE_COMMIT"),
      "\nManifest: exports/map_manifest.csv\nSource checksums and geometry: analysis_records/exports-job-", job,
      "\nCompleted UTC: ", format(Sys.time(), tz = "UTC"), "\n", sep = "", file = file.path(run, "run_record.org"), append = TRUE)
}
print(table(manifest$stat, manifest$kind))
