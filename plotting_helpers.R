plot_metabo_results <- function(res) {
  ln <- tolower(names(res))
  pcol <- names(res)[which.max(c(
    match("p.value", ln), match("pvalue", ln), match("p_val", ln),
    as.numeric(grepl("^p(\\.|_)?value$", ln))
  ))]
  if (is.na(pcol)) stop("Couldn't find a P-value column.")
  adj_candidates <- c("adj.p.val", "adj.pvalue", "fdr", "q.value", "qvalue")
  adj_idx <- match(adj_candidates, ln)
  adjcol <- names(res)[adj_idx[!is.na(adj_idx)][1]]
  effect_candidates <- grep(
    "(^coef$)|(^logfc$)|contrast|vs|beta|estimate|^t$|^b$|^F$",
    names(res), ignore.case = TRUE, value = TRUE
  )
  effectcol <- NULL
  if (length(effect_candidates)) {
    pick <- setdiff(effect_candidates, c("F", "f"))
    effectcol <- if (length(pick)) pick[1] else effect_candidates[1]
  }

  feat_lab_col <- NULL
  for (cand in c("Name", "name", "Metabolite", "metabolite", "Feature", "feature", "Compound", "compound")) {
    if (cand %in% names(res)) {
      feat_lab_col <- cand
      break
    }
  }
  if (is.null(feat_lab_col) && !is.null(rownames(res))) {
    res$.__feat__ <- rownames(res)
    feat_lab_col <- ".__feat__"
  }

  res$neglog10p <- -log10(res[[pcol]])
  if (!is.null(adjcol)) {
    res$signif <- res[[adjcol]] < 0.05
  } else {
    res$signif <- res[[pcol]] < 0.05
  }

  if (!is.null(effectcol) && !identical(effectcol, "F")) {
    p_volcano <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[effectcol]], y = neglog10p, color = signif)) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::labs(
        x = effectcol,
        y = "-log10(p)",
        color = if (!is.null(adjcol)) paste0(adjcol, " < 0.05") else paste0(pcol, " < 0.05")
      ) +
      ggplot2::theme_minimal(base_size = 12)
    ord <- order(res[[pcol]], na.last = NA)
    lab_df <- res[head(ord, 12), , drop = FALSE]
    if (!is.null(feat_lab_col) && requireNamespace("ggrepel", quietly = TRUE)) {
      p_volcano <- p_volcano + ggrepel::geom_text_repel(
        data = lab_df,
        ggplot2::aes(label = .data[[feat_lab_col]]),
        size = 3
      )
    }
    print(p_volcano)
    ggplot2::ggsave("lm_volcano.pdf", p_volcano, width = 7, height = 5)
  } else if ("F" %in% names(res)) {
    p_anova <- ggplot2::ggplot(res, ggplot2::aes(x = .data[["F"]], y = neglog10p)) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::labs(x = "F-statistic", y = "-log10(p)") +
      ggplot2::theme_minimal(base_size = 12)
    print(p_anova)
    ggplot2::ggsave("anova_scatter.pdf", p_anova, width = 7, height = 5)
  } else {
    stop("Couldn't determine an effect/coef column in MetaboAnalyst result table.")
  }

  out <- res
  names(out)[names(out) == pcol] <- "P.Value"
  if (!is.null(adjcol)) names(out)[names(out) == adjcol] <- "adj.P.Val"
  if (!is.null(effectcol)) names(out)[names(out) == effectcol] <- "Effect"
  utils::write.csv(out, "lm_results_tidy.csv", row.names = FALSE)
}

build_per_contrast <- function(
  X, meta, group_col, control_label, covariates_adj,
  alpha_sig, axes_use_fdr, sig_use_fdr
) {
  pc_unadj_hc3 <- per_contrast_pvals_hc3(
    X, meta, covars = NULL, group_col = group_col, control_label = control_label
  ) |>
    dplyr::rename(p_unadj_raw_hc3 = p_raw_hc3, p_unadj_fdr_hc3 = p_fdr_hc3)
  pc_adj_hc3 <- per_contrast_pvals_hc3(
    X, meta, covars = covariates_adj, group_col = group_col, control_label = control_label
  ) |>
    dplyr::rename(p_adj_raw_hc3 = p_raw_hc3, p_adj_fdr_hc3 = p_fdr_hc3)

  pc <- dplyr::inner_join(pc_unadj_hc3, pc_adj_hc3, by = c("feature", "Contrast"))
  pc <- pc |>
    dplyr::mutate(
      x = if (axes_use_fdr) safe_nlog10(p_unadj_fdr_hc3) else safe_nlog10(p_unadj_raw_hc3),
      y = if (axes_use_fdr) safe_nlog10(p_adj_fdr_hc3) else safe_nlog10(p_adj_raw_hc3),
      is_sig_unadj = if (sig_use_fdr) p_unadj_fdr_hc3 < alpha_sig else p_unadj_raw_hc3 < alpha_sig,
      is_sig_adj = if (sig_use_fdr) p_adj_fdr_hc3 < alpha_sig else p_adj_raw_hc3 < alpha_sig,
      sig_status = dplyr::case_when(
        is_sig_unadj & is_sig_adj ~ "Sig in both",
        is_sig_unadj & !is_sig_adj ~ "Sig unadjusted only",
        !is_sig_unadj & is_sig_adj ~ "Sig adjusted only",
        TRUE ~ "Not sig"
      )
    )
  pc
}

plot_per_contrast_scatter <- function(
  pc, alpha_sig, axes_use_fdr, label_mode, label_size, show_only_sig,
  output_pdf = "per_contrast_unadj_vs_adj.pdf",
  output_pdf_labeled = "per_contrast_unadj_vs_adj_labeled.pdf",
  facet_col = "Contrast"
) {
  pc_plot <- if (show_only_sig) dplyr::filter(pc, sig_status != "Not sig") else pc
  p <- ggplot2::ggplot(pc_plot, ggplot2::aes(x = x, y = y, color = sig_status)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), nrow = 1) +
    ggplot2::labs(
      title = "Per-contrast: unadjusted vs covariate-adjusted",
      subtitle = paste0(
        "Axes & significance use ", if (axes_use_fdr) "FDR" else "raw p",
        " (alpha = ", alpha_sig, ")"
      ),
      x = expression("-log"[10] * "(p) unadjusted"),
      y = expression("-log"[10] * "(p) adjusted"),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  ggplot2::ggsave(output_pdf, p, width = 12, height = 4.8)

  label_filter <- switch(
    label_mode,
    "adj" = quote(is_sig_adj),
    "changed" = quote(sig_status %in% c("Sig adjusted only", "Sig unadjusted only")),
    "any" = quote(sig_status != "Not sig"),
    quote(is_sig_adj)
  )
  lab_df <- pc |>
    dplyr::filter(!!label_filter) |>
    dplyr::select(dplyr::all_of(c(facet_col, "feature", "x", "y", "sig_status"))) |>
    dplyr::distinct()
  message("Labels to draw: ", nrow(lab_df))
  if (nrow(lab_df) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p_lab <- p + ggrepel::geom_text_repel(
      data = lab_df,
      mapping = ggplot2::aes(x = x, y = y, label = feature, color = sig_status),
      inherit.aes = FALSE,
      size = label_size,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      seed = 42
    )
    print(p_lab)
    ggplot2::ggsave(output_pdf_labeled, p_lab, width = 12, height = 4.8)
  }
  invisible(p)
}

plot_pairwise_scatter <- function(
  pw, alpha_sig, label_size,
  title_text, output_pdf, output_pdf_labeled = NULL
) {
  p <- ggplot2::ggplot(pw, ggplot2::aes(x = x, y = y, color = sig_status)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::facet_wrap(~ Contrast, nrow = 1, scales = "free_x") +
    ggplot2::labs(
      title = title_text,
      subtitle = paste0("Axes & significance use FDR (alpha = ", alpha_sig, ")"),
      x = expression("-log"[10] * "(p) unadjusted"),
      y = expression("-log"[10] * "(p) adjusted"),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(output_pdf, p, width = 12, height = 4.8)

  if (!is.null(output_pdf_labeled) && requireNamespace("ggrepel", quietly = TRUE)) {
    lab_df <- pw |>
      dplyr::filter(sig_status == "Sig in both") |>
      dplyr::distinct(Contrast, feature, x, y, sig_status)
    if (nrow(lab_df)) {
      p_lab <- p + ggrepel::geom_text_repel(
        data = lab_df,
        mapping = ggplot2::aes(x = x, y = y, label = feature, color = sig_status),
        inherit.aes = FALSE,
        size = label_size,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        seed = 42
      )
      ggplot2::ggsave(output_pdf_labeled, p_lab, width = 12, height = 4.8)
    }
  }
  invisible(p)
}

plot_covariate_adjusted_feature_bars <- function(
  X, meta, pc, group_col, control_label, covariates_adj, alpha_sig, prefix
) {
  annotate_use <- "fdr"
  pc_map <- pc |>
    dplyr::mutate(
      group2 = sub("_vs_.*$", "", as.character(Contrast)),
      p_use = dplyr::case_when(
        annotate_use == "fdr" ~ p_adj_fdr_hc3,
        TRUE ~ p_adj_raw_hc3
      )
    ) |>
    dplyr::select(feature, group2, p_use)

  sig_pc <- pc |>
    dplyr::filter(is.finite(p_adj_fdr_hc3), p_adj_fdr_hc3 < alpha_sig) |>
    dplyr::select(feature) |>
    dplyr::distinct()
  feat_use <- sig_pc$feature
  if (!length(feat_use)) {
    message("No adjusted-significant features found; skipping residual bar plots.")
    return(invisible(NULL))
  }

  m2 <- meta
  for (v in covariates_adj) if (is.character(m2[[v]])) m2[[v]] <- factor(m2[[v]])
  Z <- stats::model.matrix(stats::as.formula(paste("~", paste(covariates_adj, collapse = " + "))), data = m2)
  fit_cov <- limma::lmFit(X, Z)
  Res_cov <- limma::residuals.MArrayLM(fit_cov, y = X)

  lvl_order <- switch(
    prefix,
    "M1" = c("1_Control", "2_Benign", "3_Malignant and Benign", "4_Malignant"),
    "M3" = c("1. No tumour", "2. Small tumour", "3. Medium tumour", "4. Large tumour"),
    unique(meta[[group_col]])
  )

  keep_feats <- intersect(feat_use, rownames(Res_cov))
  if (!length(keep_feats)) return(invisible(NULL))
  df_long <- Res_cov[keep_feats, , drop = FALSE] |>
    tibble::as_tibble(rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "SampleID", values_to = "resid") |>
    dplyr::mutate(group = meta[SampleID, group_col, drop = TRUE]) |>
    dplyr::filter(!is.na(group)) |>
    dplyr::mutate(group = factor(group, levels = lvl_order))

  safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  for (f in unique(df_long$feature)) {
    dfi <- dplyr::filter(df_long, feature == f)
    summ <- dfi |>
      dplyr::group_by(group) |>
      dplyr::summarise(
        mean = mean(resid, na.rm = TRUE),
        se = stats::sd(resid, na.rm = TRUE) / sqrt(sum(is.finite(resid))),
        .groups = "drop"
      ) |>
      dplyr::mutate(group = factor(group, levels = lvl_order))

    n_lab <- max(1L, floor(0.30 * nrow(dfi)))
    dfi_low <- dplyr::slice_min(dfi, resid, n = n_lab, with_ties = FALSE)
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data = summ, ggplot2::aes(x = group, y = mean), width = 0.6, alpha = 0.85) +
      ggplot2::geom_errorbar(
        data = summ,
        ggplot2::aes(x = group, ymin = mean - se, ymax = mean + se),
        width = 0.2
      ) +
      ggplot2::geom_point(
        data = dfi,
        ggplot2::aes(x = group, y = resid),
        position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 42),
        alpha = 0.7,
        size = 1.6
      ) +
      ggplot2::labs(
        title = f,
        x = group_col,
        y = "Covariate-adjusted residual (log10 intensity)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::scale_x_discrete(limits = lvl_order, drop = FALSE)

    levs <- setdiff(levels(dfi$group), control_label)
    y_base <- max(summ$mean + summ$se, na.rm = TRUE)
    y_span <- diff(range(dfi$resid, na.rm = TRUE))
    if (!is.finite(y_span) || y_span == 0) y_span <- 1
    y_step <- 0.10 * y_span
    ann_list <- lapply(seq_along(levs), function(i) {
      g2 <- levs[i]
      pv <- pc_map |>
        dplyr::filter(feature == f, group2 == as.character(g2)) |>
        dplyr::pull(p_use)
      if (!length(pv) || !is.finite(pv)) return(NULL)
      lab <- as.character(stars_from_p(pv))
      if (lab == "") return(NULL)
      tibble::tibble(group1 = control_label, group2 = g2, y.position = y_base + i * y_step, label = lab)
    })
    ann <- dplyr::bind_rows(ann_list)
    if (nrow(ann) && requireNamespace("ggpubr", quietly = TRUE)) {
      p <- p +
        ggplot2::expand_limits(y = max(ann$y.position, na.rm = TRUE) + 0.05 * y_span) +
        ggpubr::stat_pvalue_manual(
          ann,
          label = "label",
          xmin = "group1",
          xmax = "group2",
          y.position = "y.position",
          tip.length = 0.01,
          size = 5
        )
    }

    p <- p + ggplot2::geom_point(
      data = dfi_low,
      ggplot2::aes(x = group, y = resid),
      shape = 21,
      stroke = 0.6,
      size = 2.6
    )
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = dfi_low,
        ggplot2::aes(x = group, y = resid, label = SampleID),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.2,
        seed = 42
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = dfi_low,
        ggplot2::aes(x = group, y = resid, label = SampleID),
        vjust = 1.2,
        size = 3
      )
    }
    ggplot2::ggsave(paste0(prefix, "_bar_", safe(f), ".pdf"), p, width = 6.5, height = 4.2)
  }
  invisible(NULL)
}
