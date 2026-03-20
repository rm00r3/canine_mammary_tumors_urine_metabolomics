ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

with_working_dir <- function(path, code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  force(code)
}

model_output_dir <- function(driver_cfg, model_cfg) {
  file.path(driver_cfg$output_root, model_cfg$model_tag)
}

write_csv_compat <- function(x, file, use_readr = FALSE, row.names = FALSE) {
  if (use_readr) {
    readr::write_csv(x, file)
  } else {
    utils::write.csv(x, file = file, row.names = row.names)
  }
  invisible(file)
}
