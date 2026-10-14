# Plot group contrasts, nested activations, and activation-masked effects.
# Run in R or with Rscript after setting the external data and Workbench paths.
# Read results/bglm2.rds and save PNG files under the external plots/bglm2_plots directory.

# Step1. Set the environment and read the saved group model
library(ciftiTools)
library(INLA)
library(BayesfMRI)

wb_path <- "/path/to/workbench"  # Replace with the local Connectome Workbench installation.
ciftiTools.setOption("wb_path", wb_path)

data_path <- "/Volumes/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20260313_BayesfMRI2"
output_path <- file.path(data_path, "results")
plot_path <- file.path(data_path, "plots")
bglm2_plot_path <- file.path(plot_path, "bglm2_plots")
dir.create(bglm2_plot_path, showWarnings = FALSE, recursive = TRUE)
res <- readRDS(file.path(output_path, "bglm2.rds"))

# Step2. Export the original plot set for each contrast
# Retain both surface templates and surface/volume exports for comparability.
# Activation and masked-effect plots use level 1 (alpha = 0.05 in the main fit).
for (i in seq_along(res$activation_levels_xii)) {
  contrast_name <- names(res$activation_levels_xii)[i]
  contrast_plot_path <- file.path(bglm2_plot_path, contrast_name)
  dir.create(contrast_plot_path, showWarnings = FALSE, recursive = TRUE)

  if (grepl("average", contrast_name)) {
    zlim <- c(-1, 1)
  } else {
    zlim <- c(-0.5, 0.5)
  }

  for (template in c("inflated", "very inflated")) {
    template_plot_path <- file.path(contrast_plot_path, template)
    dir.create(template_plot_path, showWarnings = FALSE, recursive = TRUE)

    for (what in c("surface", "volume")) {
      plot(res, stat = "contrast", idx = i, what = what,
           surfL = template, surfR = template, zlim = zlim,
           fname = file.path(template_plot_path, paste0(template, "_", what, "_", contrast_name, ".png")))

      plot(res, stat = "nested_activations", idx = i, level = 1, what = what,
           surfL = template, surfR = template,
           fname = file.path(template_plot_path, paste0(template, "_", what, "_", contrast_name, "_activation.png")))

      plot(res, stat = "activated_level_contrast", idx = i, level = 1, what = what,
           surfL = template, surfR = template, zlim = zlim,
           fname = file.path(template_plot_path, paste0(template, "_", what, "_", contrast_name, "_masked.png")))
    }
  }
}
