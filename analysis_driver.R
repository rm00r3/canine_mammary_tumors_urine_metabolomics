# Mammary Urine Analysis Driver
# Required input files:
# - newjunemsms_allmodes_msms_data_new.csv (samples in columns)
# - Meta_all_fixed_test.csv (metadata with sample IDs in rownames-compatible first column)
#
# Expected project structure:
# - analysis_driver.R
# - R/
# - outputs/
#   - model1/
#   - model2/
#   - model3/
#   - model4/
#
# How to choose model(s):
# - options(mammary.models = "model1"); source("analysis_driver.R")
# - options(mammary.models = "model2"); source("analysis_driver.R")
# - options(mammary.models = "model3"); source("analysis_driver.R")
# - options(mammary.models = "model4"); source("analysis_driver.R")
# - options(mammary.models = "all");    source("analysis_driver.R")
#
# Optional input/output overrides:
# - options(mammary.data_file = "/abs/path/to/data.csv")
# - options(mammary.meta_file = "/abs/path/to/meta.csv")
# - options(mammary.output_root = "/abs/path/to/outputs")
#
# Outputs are written under outputs/<model_tag>/ with original filenames preserved.

.driver_file <- NULL
for (i in rev(seq_along(sys.frames()))) {
  fi <- sys.frames()[[i]]
  if (!is.null(fi$ofile)) {
    .driver_file <- normalizePath(as.character(fi$ofile), winslash = "/", mustWork = FALSE)
    break
  }
}
if (is.null(.driver_file)) {
  .driver_file <- normalizePath("analysis_driver.R", winslash = "/", mustWork = FALSE)
}
.project_root <- dirname(.driver_file)

source(file.path(.project_root, "R", "config.R"), local = FALSE)
source(file.path(.project_root, "R", "io.R"), local = FALSE)
source(file.path(.project_root, "R", "alignment_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "metaboanalyst_pipeline.R"), local = FALSE)
source(file.path(.project_root, "R", "hc3_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "pairwise_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "logistic_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "robust_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "permutation_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "diagnostics_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "residual_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "plotting_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "heatmap_helpers.R"), local = FALSE)
source(file.path(.project_root, "R", "model_runner.R"), local = FALSE)

driver_cfg <- get_driver_config(project_root = .project_root)
all_model_cfgs <- get_model_configs()
models_to_run <- parse_models_to_run(getOption("mammary.models", "all"))

if (!isTRUE(getOption("mammary.auto_run", TRUE))) {
  message("Loaded analysis modules; auto-run disabled by options(mammary.auto_run = FALSE).")
} else {
  validate_required_packages(models_to_run)
  ensure_dir(driver_cfg$output_root)
  for (mid in models_to_run) {
    ensure_dir(file.path(driver_cfg$output_root, mid))
  }
  for (mid in models_to_run) {
    run_single_model(model_cfg = all_model_cfgs[[mid]], driver_cfg = driver_cfg)
  }
}
