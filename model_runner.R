empty_perm_table <- function() {
  tibble::tibble(
    feature = character(),
    Contrast = character(),
    t_obs = numeric(),
    p_perm = numeric(),
    p_perm_fdr = numeric()
  )
}

run_single_model <- function(model_cfg, driver_cfg) {
  out_dir <- model_output_dir(driver_cfg, model_cfg)
  ensure_dir(out_dir)
  message("Running ", model_cfg$model_tag, " -> ", out_dir)

  with_working_dir(out_dir, {
    OUT_PREFIX <- model_cfg$OUT_PREFIX
    alpha_sig <- driver_cfg$alpha_sig
    label_size <- driver_cfg$label_size
    show_only_sig <- driver_cfg$show_only_sig
    axes_use_fdr <- driver_cfg$axes_use_fdr
    sig_use_fdr <- driver_cfg$sig_use_fdr
    label_mode <- driver_cfg$label_mode
    group_col <- model_cfg$group_col
    control_label <- model_cfg$control_label
    covariates_adj <- model_cfg$covariates

    message("OUT_PREFIX = ", OUT_PREFIX)
    met <- run_metaboanalyst_pipeline(model_cfg, driver_cfg)
    aligned <- align_expression_and_meta(met$X, met$meta, group_col = group_col)
    X <- aligned$X
    meta <- aligned$meta
    plot_metabo_results(met$metabo_result)

    ov <- check_overlap(meta, group_col, covariates_adj)
    utils::write.csv(ov, file = model_cfg$covariatecheck_file, row.names = FALSE)

    design <- stats::model.matrix(
      stats::as.formula(paste0("~ ", group_col, " + ", paste(covariates_adj, collapse = " + "))),
      data = droplevels(meta)
    )
    prot_idx <- c(
      which(colnames(design) == "(Intercept)"),
      grep(paste0("^", group_col), colnames(design))
    )
    design <- ensure_full_rank(design, protect_idx = prot_idx, verbose = TRUE)

    pw <- NULL
    if (isTRUE(model_cfg$run_pairwise)) {
      lvl_order <- model_cfg$pairwise_disease_levels
      if (!is.null(lvl_order)) {
        available_levels <- levels(factor(meta[[group_col]]))
        missing_levels <- setdiff(lvl_order, available_levels)
        if (length(missing_levels) > 0L) {
          stop(
            paste0(
              "Configured pairwise_disease_levels not found in aligned metadata for ",
              model_cfg$model_tag,
              " (group_col = '", group_col, "'): ",
              paste(shQuote(missing_levels), collapse = ", "),
              ". Available levels: ",
              paste(shQuote(available_levels), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }
      pw_unadj <- pairwise_pvals_fast_ols(
        X, meta, group_col = group_col,
        covariates = NULL,
        disease_levels = lvl_order,
        control_label = control_label
      ) |>
        dplyr::rename(p_unadj_raw_hc3 = p_raw_hc3, p_unadj_fdr_hc3 = p_fdr_hc3)
      pw_adj <- pairwise_pvals_fast_ols(
        X, meta, group_col = group_col,
        covariates = covariates_adj,
        disease_levels = lvl_order,
        control_label = control_label
      ) |>
        dplyr::rename(p_adj_raw_hc3 = p_raw_hc3, p_adj_fdr_hc3 = p_fdr_hc3)

      by_cols <- intersect(names(pw_unadj), names(pw_adj))
      by_cols <- intersect(by_cols, c("feature", "Contrast", "L1", "L2"))
      pw <- dplyr::left_join(pw_unadj, pw_adj, by = by_cols, suffix = c("_unadj", "_adj")) |>
        dplyr::mutate(
          x = -log10(pmax(p_unadj_fdr_hc3, .Machine$double.xmin)),
          y = -log10(pmax(p_adj_fdr_hc3, .Machine$double.xmin)),
          is_sig_unadj = p_unadj_fdr_hc3 < alpha_sig,
          is_sig_adj = p_adj_fdr_hc3 < alpha_sig,
          sig_status = dplyr::case_when(
            is_sig_unadj & is_sig_adj ~ "Sig in both",
            is_sig_unadj & !is_sig_adj ~ "Sig unadjusted only",
            !is_sig_unadj & is_sig_adj ~ "Sig adjusted only",
            TRUE ~ "Not sig"
          )
        )

      if (!is.null(model_cfg$pairwise_file)) {
        if (!is.null(model_cfg$covariatecheck_pairwise_file) && identical(model_cfg$model_tag, "model1")) {
          ov_m1 <- check_overlap(meta, group_col = "tumour_malignancy", covars = covariates_adj)
          utils::write.csv(ov_m1, model_cfg$covariatecheck_pairwise_file, row.names = FALSE)
        }
        if (identical(model_cfg$model_tag, "model1")) {
          pw_out <- dplyr::select(
            pw,
            dplyr::any_of(c(
              "feature", "L1", "L2", "Contrast",
              "p_unadj_raw_hc3", "p_unadj_fdr_hc3",
              "p_adj_raw_hc3", "p_adj_fdr_hc3",
              "is_sig_unadj", "is_sig_adj", "sig_status"
            ))
          )
        } else {
          pw_out <- dplyr::select(
            pw,
            dplyr::any_of(c(
              "feature", "Contrast",
              "p_unadj_raw_hc3", "p_unadj_fdr_hc3",
              "p_adj_raw_hc3", "p_adj_fdr_hc3",
              "is_sig_unadj", "is_sig_adj", "sig_status"
            ))
          )
        }
        readr::write_csv(pw_out, model_cfg$pairwise_file)
      }
      if (!is.null(model_cfg$pairwise_top_file)) {
        pw_top <- pw |>
          dplyr::filter(is_sig_adj) |>
          dplyr::group_by(Contrast) |>
          dplyr::arrange(p_adj_fdr_hc3, .by_group = TRUE) |>
          dplyr::mutate(rank_adj_fdr = dplyr::row_number()) |>
          dplyr::slice_head(n = 50) |>
          dplyr::ungroup()
        readr::write_csv(pw_top, model_cfg$pairwise_top_file)
      }
      if (!is.null(model_cfg$pairwise_plot_file)) {
        title_text <- switch(
          model_cfg$model_tag,
          model1 = "Disease-only pairwise (tumour malignancy): unadjusted vs covariate-adjusted",
          model3 = "Disease-only pairwise (size): unadjusted vs covariate-adjusted",
          model4 = "Disease-only pairwise (tumour number): unadjusted vs covariate-adjusted",
          "Disease-only pairwise: unadjusted vs covariate-adjusted"
        )
        plot_pairwise_scatter(
          pw = pw,
          alpha_sig = alpha_sig,
          label_size = label_size,
          title_text = title_text,
          output_pdf = model_cfg$pairwise_plot_file,
          output_pdf_labeled = model_cfg$pairwise_plot_labeled_file
        )
      }
      if (isTRUE(model_cfg$run_pairwise_bars) && !is.null(model_cfg$pairwise_bar_dir)) {
        plot_model3_pairwise_bars(
          X = X,
          meta = meta,
          pw = pw,
          group_col = group_col,
          covariates_adj = covariates_adj,
          alpha_sig = alpha_sig,
          out_dir = model_cfg$pairwise_bar_dir
        )
      }
    }

    pc <- build_per_contrast(
      X = X,
      meta = meta,
      group_col = group_col,
      control_label = control_label,
      covariates_adj = covariates_adj,
      alpha_sig = alpha_sig,
      axes_use_fdr = axes_use_fdr,
      sig_use_fdr = sig_use_fdr
    )

    if (identical(model_cfg$model_tag, "model1")) {
      pc <- pc |>
        dplyr::mutate(
          Contrast_raw = Contrast,
          Contrast_plot = factor(
            pretty_contrast_label(Contrast_raw),
            levels = pretty_contrast_label(unique(Contrast_raw))
          )
        )
      plot_per_contrast_scatter(
        pc = pc,
        alpha_sig = alpha_sig,
        axes_use_fdr = axes_use_fdr,
        label_mode = label_mode,
        label_size = label_size,
        show_only_sig = show_only_sig,
        output_pdf = "per_contrast_unadj_vs_adj.pdf",
        output_pdf_labeled = "per_contrast_unadj_vs_adj_labeled.pdf",
        facet_col = "Contrast_plot"
      )
      pc_out <- pc |>
        dplyr::select(
          feature, Contrast_raw, Contrast_plot,
          p_unadj_raw_hc3, p_unadj_fdr_hc3,
          p_adj_raw_hc3, p_adj_fdr_hc3,
          is_sig_unadj, is_sig_adj, sig_status
        )
      utils::write.csv(pc_out, model_cfg$per_contrast_file, row.names = FALSE)
    } else {
      pc <- pc |>
        dplyr::mutate(
          Contrast = factor(
            Contrast,
            levels = unique(Contrast),
            labels = pretty_contrast_label(unique(Contrast))
          )
        )
      plot_per_contrast_scatter(
        pc = pc,
        alpha_sig = alpha_sig,
        axes_use_fdr = axes_use_fdr,
        label_mode = label_mode,
        label_size = label_size,
        show_only_sig = show_only_sig,
        output_pdf = "per_contrast_unadj_vs_adj.pdf",
        output_pdf_labeled = "per_contrast_unadj_vs_adj_labeled.pdf",
        facet_col = "Contrast"
      )
      pc_out <- pc |>
        dplyr::select(
          feature, Contrast,
          p_unadj_raw_hc3, p_unadj_fdr_hc3,
          p_adj_raw_hc3, p_adj_fdr_hc3,
          is_sig_unadj, is_sig_adj, sig_status
        )
      utils::write.csv(pc_out, model_cfg$per_contrast_file, row.names = FALSE)
      if (!is.null(model_cfg$per_contrast_top_file)) {
        pc_top <- pc |>
          dplyr::filter(is_sig_adj) |>
          dplyr::group_by(Contrast) |>
          dplyr::arrange(p_adj_fdr_hc3, .by_group = TRUE) |>
          dplyr::mutate(rank_adj_fdr = dplyr::row_number()) |>
          dplyr::slice_head(n = 50) |>
          dplyr::ungroup()
        utils::write.csv(pc_top, model_cfg$per_contrast_top_file, row.names = FALSE)
      }
    }

    if (identical(model_cfg$or_mode, "per_contrast_joined")) {
      feat_by_contrast <- pc |>
        dplyr::filter(is_sig_adj) |>
        dplyr::group_by(Contrast_raw) |>
        dplyr::summarise(features = list(unique(feature)), .groups = "drop")

      if (nrow(feat_by_contrast)) {
        or_tabs <- lapply(seq_len(nrow(feat_by_contrast)), function(i) {
          this_raw <- feat_by_contrast$Contrast_raw[i]
          feats <- feat_by_contrast$features[[i]]
          tmp <- per_contrast_or_hc3(
            X = X, meta = meta, features = feats,
            group_col = group_col, control_label = control_label, covariates = covariates_adj
          )
          dplyr::filter(tmp, Contrast == this_raw)
        })
        or_tab <- dplyr::bind_rows(or_tabs)
      } else {
        or_tab <- tibble::tibble()
      }

      if (nrow(or_tab)) {
        or_tab_joined <- or_tab |>
          dplyr::left_join(
            pc |>
              dplyr::select(feature, Contrast_raw, q_lm = p_adj_fdr_hc3),
            by = c("feature" = "feature", "Contrast" = "Contrast_raw")
          ) |>
          dplyr::group_by(Contrast) |>
          dplyr::mutate(q_or = stats::p.adjust(p, "fdr")) |>
          dplyr::ungroup() |>
          dplyr::mutate(Contrast_plot = pretty_contrast_label(Contrast))
      } else {
        or_tab_joined <- tibble::tibble(
          feature = character(), Contrast = character(), OR = numeric(), OR_low = numeric(),
          OR_high = numeric(), p = numeric(), q = numeric(), se_source = character(),
          n = numeric(), n_case = numeric(), n_ctrl = numeric(), sd_all = numeric(),
          sd_case = numeric(), sd_ctrl = numeric(), overlap_flag = logical(), reason = character(),
          q_lm = numeric(), q_or = numeric(), Contrast_plot = character()
        )
      }
      utils::write.csv(or_tab_joined, model_cfg$or_raw_file, row.names = FALSE)
      or_tab_fmt <- or_tab_joined |>
        dplyr::mutate(
          `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high),
          `p (logistic)` = formatC(p, format = "e", digits = 2),
          `q_OR (FDR)` = formatC(q_or, format = "e", digits = 2),
          `q_LM (FDR)` = formatC(q_lm, format = "e", digits = 2)
        ) |>
        dplyr::select(
          Feature = feature, Contrast = Contrast_plot,
          `OR (95% CI)`, `p (logistic)`, `q_OR (FDR)`, `q_LM (FDR)`
        )
      utils::write.csv(or_tab_fmt, model_cfg$or_fmt_file, row.names = FALSE)
    }

    if (identical(model_cfg$or_mode, "case_control")) {
      stop_if_empty(meta, group_col, control_label, covariates_adj, min_per_class = 3)
      feat_for_or <- pc |>
        dplyr::filter(p_adj_fdr_hc3 < alpha_sig) |>
        dplyr::pull(feature) |>
        unique()
      if (length(feat_for_or)) {
        or_cc <- case_control_or_hc3(
          X = X, meta = meta, features = feat_for_or,
          group_col = group_col, case_label = model_cfg$pairwise_case_label,
          control_label = control_label, covariates = covariates_adj
        )
        or_cc_out <- or_cc |>
          dplyr::mutate(
            `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high),
            `p (HC3/Firth)` = formatC(p, format = "e", digits = 2),
            `q (FDR)` = formatC(q, format = "e", digits = 2)
          ) |>
          dplyr::select(Feature = feature, Contrast, `OR (95% CI)`, `p (HC3/Firth)`, `q (FDR)`, Note = se_source)
        utils::write.csv(or_cc, model_cfg$or_raw_file, row.names = FALSE)
        utils::write.csv(or_cc_out, model_cfg$or_fmt_file, row.names = FALSE)
      } else {
        message("No HC3-FDR significant features - skipping OR table for ", model_cfg$model_tag)
      }
    }

    if (identical(model_cfg$or_mode, "per_contrast")) {
      feat_for_or <- pc |>
        dplyr::filter(p_adj_fdr_hc3 < alpha_sig) |>
        dplyr::pull(feature) |>
        unique()
      if (length(feat_for_or)) {
        or_tab <- per_contrast_or_hc3(
          X = X, meta = meta, features = feat_for_or,
          group_col = group_col, control_label = control_label, covariates = covariates_adj
        )
      } else {
        message(model_cfg$model_tag, ": no HC3-FDR significant features - skipping OR table.")
        or_tab <- tibble::tibble(
          feature = character(), Contrast = character(), OR = numeric(), OR_low = numeric(),
          OR_high = numeric(), p = numeric(), q = numeric(), se_source = character(),
          n = numeric(), n_case = numeric(), n_ctrl = numeric(), sd_all = numeric(),
          sd_case = numeric(), sd_ctrl = numeric(), overlap_flag = logical(), reason = character()
        )
      }
      utils::write.csv(or_tab, model_cfg$or_raw_file, row.names = FALSE)
      if (!is.null(model_cfg$or_diag_file)) utils::write.csv(or_tab, model_cfg$or_diag_file, row.names = FALSE)
      or_tab_out <- or_tab |>
        dplyr::mutate(
          `OR (95% CI)` = ifelse(is.finite(OR), sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high), "NA (NA-NA)"),
          `p (HC3/Firth)` = ifelse(is.finite(p), formatC(p, format = "e", digits = 2), "NaN"),
          `q (FDR)` = ifelse(is.finite(q), formatC(q, format = "e", digits = 2), "NaN")
        )
      if ("se_source" %in% names(or_tab_out)) {
        or_tab_out <- or_tab_out |>
          dplyr::select(feature, Contrast, `OR (95% CI)`, `p (HC3/Firth)`, `q (FDR)`, se_source) |>
          dplyr::rename(Feature = feature, Note = se_source)
      } else {
        or_tab_out <- or_tab_out |>
          dplyr::select(feature, Contrast, `OR (95% CI)`, `p (HC3/Firth)`, `q (FDR)`) |>
          dplyr::mutate(Note = NA_character_) |>
          dplyr::rename(Feature = feature)
      }
      utils::write.csv(or_tab_out, model_cfg$or_fmt_file, row.names = FALSE)
    }

    if (!is.null(model_cfg$residual_covariate_file)) {
      export_covariate_only_residuals(
        X = X,
        meta = meta,
        covariates = covariates_adj,
        outfile = model_cfg$residual_covariate_file
      )
    }

    if (isTRUE(model_cfg$run_covariate_barplots)) {
      plot_covariate_adjusted_feature_bars(
        X = X, meta = meta, pc = pc,
        group_col = group_col, control_label = control_label,
        covariates_adj = covariates_adj, alpha_sig = alpha_sig, prefix = OUT_PREFIX
      )
    }

    diag_adj <- influence_diag_strict(X, meta, covariates_adj, group_col, control_label)
    utils::write.csv(diag_adj, "influence_diagnostics_adjusted.csv", row.names = FALSE)
    lev_tbl <- sample_leverage_table(meta, covariates_adj, group_col, control_label)
    utils::write.csv(lev_tbl, "sample_leverage_strict.csv", row.names = FALSE)
    cook_tbl <- sample_cook_counts(X, meta, covariates_adj, group_col, control_label, cook_thresh_strict = 1)
    utils::write.csv(cook_tbl, "sample_cook_repeat_offenders.csv", row.names = FALSE)

    n <- nrow(stats::model.matrix(
      stats::as.formula(paste("~", group_col, "+", paste(covariates_adj, collapse = "+"))),
      data = meta[stats::complete.cases(meta[, c(group_col, covariates_adj), drop = FALSE]), ]
    ))
    p <- ncol(stats::model.matrix(
      stats::as.formula(paste("~", group_col, "+", paste(covariates_adj, collapse = "+"))),
      data = meta[stats::complete.cases(meta[, c(group_col, covariates_adj), drop = FALSE]), ]
    ))
    lev_thr <- 3 * p / n
    cat(sprintf("Strict thresholds: leverage > %.3f (3p/n with n=%d, p=%d); Cook's D > 1\n", lev_thr, n, p))
    cat(sprintf("Samples with leverage > %.3f: %d/%d\n", lev_thr, sum(lev_tbl$flag_leverage), nrow(lev_tbl)))
    cat(sprintf("Top Cook repeat-offenders (prop of features > 1): %s\n",
      paste(sprintf("%s=%.3f", head(cook_tbl$sample, 5), head(cook_tbl$cook_prop, 5)), collapse = ", ")))

    top_feats <- pc |>
      dplyr::filter(p_adj_fdr_hc3 < alpha_sig) |>
      dplyr::pull(feature) |>
      unique() |>
      head(30)

    if (identical(model_cfg$model_tag, "model2")) {
      has_adj_sig <- any(pc$p_adj_fdr_hc3 < alpha_sig, na.rm = TRUE)
      if (isTRUE(has_adj_sig) && length(top_feats)) {
        rob_tab <- robust_fit_features(X, meta, top_feats, covariates_adj, group_col, control_label)
        if (nrow(rob_tab)) utils::write.csv(rob_tab, "robust_regression_subset.csv", row.names = FALSE)
        perm_tab <- permute_FL_pvals(
          X, meta, top_feats, covariates_adj, group_col, control_label,
          B = 1000, seed = 1
        )
        if (nrow(perm_tab)) utils::write.csv(perm_tab, "permutation_FL_subset.csv", row.names = FALSE)
      } else {
        message("No HC3-FDR significant features - skipping robust/permutation in Model 2.")
      }
    } else if (identical(model_cfg$model_tag, "model4")) {
      if (length(top_feats)) {
        rob_tab <- robust_fit_features(X, meta, top_feats, covariates_adj, group_col, control_label)
        utils::write.csv(rob_tab, "robust_regression_subset.csv", row.names = FALSE)
        perm_tab <- permute_FL_pvals(
          X, meta, top_feats, covariates_adj, group_col, control_label,
          B = 1000, seed = 1, robust_projection = TRUE
        )
        utils::write.csv(perm_tab, "permutation_FL_subset.csv", row.names = FALSE)
      } else {
        message("Model 4: no HC3-FDR significant features; skipping robust/permutation blocks.")
      }
    } else {
      rob_tab <- robust_fit_features(X, meta, top_feats, covariates_adj, group_col, control_label)
      utils::write.csv(rob_tab, "robust_regression_subset.csv", row.names = FALSE)
      if (length(top_feats)) {
        perm_tab <- permute_FL_pvals(
          X, meta, top_feats, covariates_adj, group_col, control_label,
          B = 1000, seed = 1
        )
      } else {
        perm_tab <- empty_perm_table()
      }
      utils::write.csv(perm_tab, "permutation_FL_subset.csv", row.names = FALSE)
    }

    res_adj_norm <- get_residuals(
      X, meta, covars = covariates_adj,
      group_col = group_col, control_label = control_label
    )
    shap_p <- apply(res_adj_norm$resid, 1, function(r) {
      if (sum(is.finite(r)) >= 3) stats::shapiro.test(r)$p.value else NA_real_
    })
    shap_fdr <- stats::p.adjust(shap_p, "fdr")
    normality_plots(res_adj_norm, shap_p, shap_fdr)

    if (identical(model_cfg$model_tag, "model1")) {
      run_na_robust_cook_summary(
        X, meta, covariates_adj, group_col, OUT_PREFIX,
        default_diag_summary_file = paste0(OUT_PREFIX, "_diag_summary.csv"),
        also_write_diagnostics_summary = TRUE
      )
      write_diag_readme(OUT_PREFIX, summary_file = paste0(OUT_PREFIX, "_diag_summary.csv"))
    } else {
      run_na_robust_cook_summary(
        X, meta, covariates_adj, group_col, OUT_PREFIX,
        default_diag_summary_file = paste0(OUT_PREFIX, "_diag_summary.csv"),
        also_write_diagnostics_summary = FALSE
      )
      write_diag_readme(OUT_PREFIX, summary_file = paste0(OUT_PREFIX, "_diag_summary.csv"))
    }

    tryCatch(
      {
        export_adjusted_significant_heatmap(
          X = X,
          meta = meta,
          pc = pc,
          group_col = group_col,
          alpha_sig = alpha_sig,
          filename = model_cfg$heatmap_file,
          stop_if_empty = model_cfg$heatmap_stop_if_empty,
          symmetric_palette = model_cfg$model_tag %in% c("model1", "model3")
        )
      },
      error = function(e) {
        message("Heatmap step skipped: ", conditionMessage(e))
        NULL
      }
    )

    invisible(TRUE)
  })
}
