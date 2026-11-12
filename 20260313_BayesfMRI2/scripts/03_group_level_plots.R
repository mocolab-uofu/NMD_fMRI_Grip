# Plot group contrasts, nested activations, and activation-masked effects.
# Run in R or with Rscript after setting the external data and Workbench paths.
# Read results/bglm2.rds and save PNG files under the external plots/bglm2_plots directory.

# Set the environment and read the saved group model
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
group_model <- readRDS(file.path(output_path, "bglm2.rds"))

# Export the original plot set for each contrast
# Retain both surface templates and surface/volume exports for comparability.
# Activation and masked-effect plots use level 1 (alpha = 0.05 in the main fit).
for (i in seq_along(group_model$activation_levels_xii)) {
  contrast_name <- names(group_model$activation_levels_xii)[i]
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
      plot_prefix <- file.path(template_plot_path, paste(template, what, contrast_name, sep = "_"))
      plot(group_model, stat = "contrast", idx = i, what = what,
           surfL = template, surfR = template, zlim = zlim,
           fname = paste0(plot_prefix, ".png"))

      plot(group_model, stat = "nested_activations", idx = i, level = 1, what = what,
           surfL = template, surfR = template,
           fname = paste0(plot_prefix, "_activation.png"))

      plot(group_model, stat = "activated_level_contrast", idx = i, level = 1, what = what,
           surfL = template, surfR = template, zlim = zlim,
           fname = paste0(plot_prefix, "_masked.png"))
    }
  }
}
