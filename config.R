`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

get_project_root <- function(default = getwd()) {
  frames <- sys.frames()
  ofiles <- vapply(
    frames,
    function(x) {
      if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_
    },
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles)) {
    return(dirname(normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(default, winslash = "/", mustWork = FALSE)
}

validate_required_packages <- function(model_ids = names(get_model_configs())) {
  base_pkgs <- c(
    "MetaboAnalystR", "ggplot2", "dplyr", "tidyr", "sandwich", "lmtest",
    "broom", "limma", "readr", "tibble", "purrr", "MASS", "ggpubr",
    "pheatmap", "ggrepel"
  )
  model_pkgs <- c()
  if (any(model_ids %in% c("model1", "model2", "model3", "model4"))) {
    model_pkgs <- c(model_pkgs, "logistf", "robustbase")
  }
  needed <- unique(c(base_pkgs, model_pkgs))
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      paste0(
        "Missing required packages: ",
        paste(sort(missing), collapse = ", "),
        ". Install them before running analysis_driver.R."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

get_model_configs <- function() {
  list(
    model1 = list(
      model_tag = "model1",
      OUT_PREFIX = "M1",
      group_col = "tumour_malignancy",
      control_label = "1_Control",
      pairwise_disease_levels = c("2_Benign", "3_Malignant and Benign", "4_Malignant"),
      pairwise_case_label = NULL,
      covariates = c("Less_exact_diet_groups", "Age_group", "Sterilization_status", "Dog_size_class_kg"),
      covariatecheck_file = "covariatecheck_M1.csv",
      covariatecheck_pairwise_file = "covariatecheck_M1_pairwise.csv",
      residual_covariate_file = "model1_residuals_covariates_only_ALL.csv",
      per_contrast_file = "M1_per_contrast_pvals_hc3.csv",
      per_contrast_top_file = NULL,
      pairwise_file = "M1_pairwise_pvals_hc3.csv",
      pairwise_top_file = "M1_pairwise_top50_adjSig.csv",
      pairwise_plot_file = "M1_pairwise_unadj_vs_adj.pdf",
      pairwise_plot_labeled_file = "M1_pairwise_unadj_vs_adj_labeled.pdf",
      pairwise_bar_dir = NULL,
      heatmap_file = "M1_heatmap_adjSig.pdf",
      or_mode = "per_contrast_joined",
      or_raw_file = "M1_or_table_raw_bycontrast.csv",
      or_fmt_file = "M1_or_table_formatted_bycontrast.csv",
      or_diag_file = NULL,
      run_pairwise = TRUE,
      run_pairwise_bars = FALSE,
      run_covariate_barplots = TRUE,
      heatmap_stop_if_empty = TRUE
    ),
    model2 = list(
      model_tag = "model2",
      OUT_PREFIX = "M2",
      group_col = "Mammary_tumour",
      control_label = "Control",
      pairwise_disease_levels = NULL,
      pairwise_case_label = "Case",
      covariates = c("Less_exact_diet_groups", "Age_group", "Sterilization_status", "Dog_size_class_kg"),
      covariatecheck_file = "covariatecheck_M2.csv",
      covariatecheck_pairwise_file = NULL,
      residual_covariate_file = NULL,
      per_contrast_file = "M2_per_contrast_pvals_hc3.csv",
      per_contrast_top_file = "M2_top50_per_contrast_adjSig.csv",
      pairwise_file = NULL,
      pairwise_top_file = NULL,
      pairwise_plot_file = NULL,
      pairwise_plot_labeled_file = NULL,
      pairwise_bar_dir = NULL,
      heatmap_file = "M2_heatmap_adjSig.pdf",
      or_mode = "case_control",
      or_raw_file = "M2_or_table_raw.csv",
      or_fmt_file = "M2_or_table_formatted.csv",
      or_diag_file = NULL,
      run_pairwise = FALSE,
      run_pairwise_bars = FALSE,
      run_covariate_barplots = FALSE,
      heatmap_stop_if_empty = FALSE
    ),
    model3 = list(
      model_tag = "model3",
      OUT_PREFIX = "M3",
      group_col = "Tumour_size",
      control_label = "1. No tumour",
      pairwise_disease_levels = c("2. Small tumour", "3. Medium tumour", "4. Large tumour"),
      pairwise_case_label = NULL,
      covariates = c("Less_exact_diet_groups", "Age_group", "Sterilization_status", "Dog_size_class_kg"),
      covariatecheck_file = "covariatecheck_M3.csv",
      covariatecheck_pairwise_file = NULL,
      residual_covariate_file = "model3_residuals_covariates_only_ALL.csv",
      per_contrast_file = "M3_per_contrast_pvals_hc3.csv",
      per_contrast_top_file = "M3_top50_per_contrast_adjSig.csv",
      pairwise_file = "M3_pairwise_pvals_hc3.csv",
      pairwise_top_file = NULL,
      pairwise_plot_file = "M3_pairwise_unadj_vs_adj.pdf",
      pairwise_plot_labeled_file = NULL,
      pairwise_bar_dir = "M3_pairwise_bars",
      heatmap_file = "M3_heatmap_adjSig.pdf",
      or_mode = "per_contrast",
      or_raw_file = "M3_or_table_raw.csv",
      or_fmt_file = "M3_or_table_formatted.csv",
      or_diag_file = "M3_or_table_with_diagnostics.csv",
      run_pairwise = TRUE,
      run_pairwise_bars = TRUE,
      run_covariate_barplots = TRUE,
      heatmap_stop_if_empty = TRUE
    ),
    model4 = list(
      model_tag = "model4",
      OUT_PREFIX = "M4",
      group_col = "Tumours",
      control_label = "1. No tumour",
      pairwise_disease_levels = c("2. Single", "3. Multiple"),
      pairwise_case_label = NULL,
      covariates = c("Less_exact_diet_groups", "Age_group", "Sterilization_status", "Dog_size_class_kg"),
      covariatecheck_file = "covariatecheck_M4.csv",
      covariatecheck_pairwise_file = NULL,
      residual_covariate_file = "model4_residuals_covariates_only_ALL.csv",
      per_contrast_file = "M4_per_contrast_pvals_hc3.csv",
      per_contrast_top_file = "M4_top50_per_contrast_adjSig.csv",
      pairwise_file = NULL,
      pairwise_top_file = NULL,
      pairwise_plot_file = NULL,
      pairwise_plot_labeled_file = NULL,
      pairwise_bar_dir = NULL,
      heatmap_file = "M4_heatmap_adjSig.pdf",
      or_mode = "per_contrast",
      or_raw_file = "M4_or_table_raw.csv",
      or_fmt_file = "M4_or_table_formatted.csv",
      or_diag_file = NULL,
      run_pairwise = TRUE,
      run_pairwise_bars = FALSE,
      run_covariate_barplots = FALSE,
      heatmap_stop_if_empty = TRUE
    )
  )
}

get_driver_config <- function(project_root = get_project_root()) {
  default_data <- "~/Desktop/Mammary_urine/Mammary_urine/data_raw/newjunemsms_allmodes_msms_data_new.csv"
  default_meta <- "~/Desktop/Mammary_urine/Mammary_urine/data_raw/Meta_all_fixed_test.csv"
  list(
    project_root = project_root,
    data_file = path.expand(getOption("mammary.data_file", default_data)),
    meta_file = path.expand(getOption("mammary.meta_file", default_meta)),
    output_root = getOption("mammary.output_root", file.path(project_root, "outputs")),
    alpha_sig = getOption("mammary.alpha_sig", 0.05),
    label_size = getOption("mammary.label_size", 2.5),
    show_only_sig = getOption("mammary.show_only_sig", FALSE),
    axes_use_fdr = getOption("mammary.axes_use_fdr", TRUE),
    sig_use_fdr = getOption("mammary.sig_use_fdr", TRUE),
    label_mode = getOption("mammary.label_mode", "adj")
  )
}

parse_models_to_run <- function(x = getOption("mammary.models", "all")) {
  model_ids <- names(get_model_configs())
  if (is.null(x)) return(model_ids)
  if (is.character(x) && length(x) == 1L && identical(tolower(x), "all")) return(model_ids)
  if (is.character(x) && length(x) == 1L && grepl(",", x, fixed = TRUE)) {
    x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  }
  x <- unique(as.character(x))
  bad <- setdiff(x, model_ids)
  if (length(bad)) {
    stop("Unknown model id(s): ", paste(bad, collapse = ", "),
         ". Allowed: ", paste(model_ids, collapse = ", "), ", all", call. = FALSE)
  }
  x
}
