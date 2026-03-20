pairwise_pvals_fast_ols <- function(
  X, meta, group_col, covariates = NULL,
  disease_levels = NULL, control_label = NULL, verbose = TRUE
) {
  mk_empty <- function() {
    tibble::tibble(
      feature = character(),
      L1 = character(),
      L2 = character(),
      Contrast = character(),
      p_raw_hc3 = numeric(),
      se_source = character(),
      p_fdr_hc3 = numeric()
    )
  }

  stopifnot(is.matrix(X))
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("feature_", seq_len(nrow(X)))
  }
  g_all <- factor(meta[[group_col]])

  if (is.null(disease_levels)) {
    levs <- levels(g_all)
    if (is.null(control_label)) {
      ctrl_like <- grep("control|no tumour|no tumor", levs, ignore.case = TRUE, value = TRUE)
      control_label <- ctrl_like[1]
    }
    disease_levels <- setdiff(levs, control_label)
  }
  disease_levels <- intersect(disease_levels, levels(g_all))
  if (length(disease_levels) < 2L) return(mk_empty())

  make_na_rows <- function(L1, L2) {
    tibble::tibble(
      feature = rownames(X),
      L1 = L1,
      L2 = L2,
      Contrast = paste0(L2, "_vs_", L1),
      p_raw_hc3 = NA_real_,
      se_source = "HC3"
    )
  }

  one_pair_hc3 <- function(L1, L2) {
    keep_pair <- meta[[group_col]] %in% c(L1, L2)
    if (!any(keep_pair)) return(make_na_rows(L1, L2))

    X_pair <- X[, keep_pair, drop = FALSE]
    meta_pair <- meta[keep_pair, , drop = FALSE]

    need <- c(group_col, if (is.null(covariates)) character() else covariates)
    cc <- stats::complete.cases(meta_pair[, need, drop = FALSE])
    X_pair <- X_pair[, cc, drop = FALSE]
    meta_pair <- meta_pair[cc, , drop = FALSE]

    grp <- factor(meta_pair[[group_col]], levels = c(L1, L2))
    if (sum(grp == L1, na.rm = TRUE) < 1L || sum(grp == L2, na.rm = TRUE) < 1L) {
      if (verbose) {
        message("pairwise HC3 [", L2, "_vs_", L1, "]: non-estimable after complete-case filtering.")
      }
      return(make_na_rows(L1, L2))
    }

    meta_pair$.__pw_group <- stats::relevel(grp, ref = L1)
    if (!is.null(covariates) && length(covariates)) {
      for (v in covariates) {
        if (is.character(meta_pair[[v]])) meta_pair[[v]] <- factor(meta_pair[[v]])
      }
    }

    frm <- stats::as.formula(
      paste0(
        "~ .__pw_group",
        if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""
      )
    )
    design_raw <- stats::model.matrix(frm, data = meta_pair)
    group_idx_raw <- grep("^\\.__pw_group", colnames(design_raw))
    if (length(group_idx_raw) != 1L) {
      if (verbose) {
        message("pairwise HC3 [", L2, "_vs_", L1, "]: non-estimable group effect in design.")
      }
      return(make_na_rows(L1, L2))
    }

    protect_idx <- c(which(colnames(design_raw) == "(Intercept)"), group_idx_raw)
    design <- ensure_full_rank(design_raw, protect_idx = protect_idx, verbose = FALSE)
    dropped <- setdiff(colnames(design_raw), colnames(design))
    if (verbose && length(dropped)) {
      message("pairwise HC3 [", L2, "_vs_", L1, "]: dropped aliased cols: ", paste(dropped, collapse = ", "))
    }

    group_colname <- colnames(design_raw)[group_idx_raw]
    g_idx <- which(colnames(design) == group_colname)
    if (length(g_idx) != 1L || qr(design)$rank < ncol(design)) {
      if (verbose) {
        message("pairwise HC3 [", L2, "_vs_", L1, "]: non-estimable pairwise design.")
      }
      return(make_na_rows(L1, L2))
    }

    df_t <- nrow(design) - ncol(design)
    if (!is.finite(df_t) || df_t <= 0) {
      if (verbose) {
        message("pairwise HC3 [", L2, "_vs_", L1, "]: non-positive residual df.")
      }
      return(make_na_rows(L1, L2))
    }

    pvals <- apply(X_pair, 1, function(y) {
      fit_lm <- tryCatch(stats::lm(y ~ design - 1), error = function(e) NULL)
      if (is.null(fit_lm)) return(NA_real_)

      cf <- tryCatch(stats::coef(fit_lm), error = function(e) NULL)
      if (is.null(cf) || is.na(cf[g_idx])) return(NA_real_)

      V <- tryCatch(sandwich::vcovHC(fit_lm, type = "HC3"), error = function(e) NULL)
      if (is.null(V) || !is.matrix(V)) return(NA_real_)

      se <- suppressWarnings(sqrt(V[g_idx, g_idx]))
      if (!is.finite(se) || se <= 0) return(NA_real_)

      2 * stats::pt(abs(cf[g_idx] / se), df = df_t, lower.tail = FALSE)
    })

    tibble::tibble(
      feature = rownames(X_pair),
      L1 = L1,
      L2 = L2,
      Contrast = paste0(L2, "_vs_", L1),
      p_raw_hc3 = as.numeric(pvals),
      se_source = "HC3"
    )
  }

  pairs <- utils::combn(disease_levels, 2, simplify = FALSE)
  out <- dplyr::bind_rows(lapply(pairs, function(pair_levels) {
    one_pair_hc3(pair_levels[[1]], pair_levels[[2]])
  }))
  if (!nrow(out)) return(mk_empty())

  out |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(p_fdr_hc3 = stats::p.adjust(p_raw_hc3, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::select(feature, L1, L2, Contrast, p_raw_hc3, se_source, p_fdr_hc3)
}

stars_from_p <- function(p) {
  cut(
    p,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("****", "***", "**", "*", ""),
    right = TRUE
  )
}

plot_model3_pairwise_bars <- function(
  X, meta, pw, group_col, covariates_adj, alpha_sig, out_dir
) {
  if (!nrow(pw)) {
    message("No pairwise table rows for Model 3; skipping pairwise bars.")
    return(invisible(NULL))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  m2 <- meta
  for (v in covariates_adj) if (is.character(m2[[v]])) m2[[v]] <- factor(m2[[v]])
  Z <- stats::model.matrix(stats::as.formula(paste("~", paste(covariates_adj, collapse = " + "))), data = m2)
  fit_cov <- limma::lmFit(X, Z)
  Res_cov <- limma::residuals.MArrayLM(fit_cov, y = X)

  pw_sig <- pw |>
    dplyr::filter(is.finite(p_adj_fdr_hc3), p_adj_fdr_hc3 < alpha_sig) |>
    tidyr::separate_wider_delim("Contrast", delim = "_vs_", names = c("G2", "G1"))

  if (!nrow(pw_sig)) {
    message("No pairwise adjusted-significant features; skipping pairwise bar plots.")
    return(invisible(NULL))
  }

  safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  for (i in seq_len(nrow(pw_sig))) {
    rowi <- pw_sig[i, , drop = FALSE]
    g1 <- as.character(rowi$G1)
    g2 <- as.character(rowi$G2)
    f <- as.character(rowi$feature)
    p_fdr <- as.numeric(rowi$p_adj_fdr_hc3)

    keep_smp <- rownames(meta)[meta[[group_col]] %in% c(g1, g2)]
    if (!length(keep_smp) || !(f %in% rownames(Res_cov))) next

    dfi <- tibble::tibble(
      SampleID = keep_smp,
      group = factor(meta[keep_smp, group_col, drop = TRUE], levels = c(g1, g2)),
      resid = as.numeric(Res_cov[f, keep_smp])
    ) |>
      dplyr::filter(is.finite(resid))

    if (nrow(dfi) < 4 || nlevels(dfi$group) < 2) next

    summ <- dfi |>
      dplyr::group_by(group) |>
      dplyr::summarise(
        mean = mean(resid, na.rm = TRUE),
        se = stats::sd(resid, na.rm = TRUE) / sqrt(sum(is.finite(resid))),
        .groups = "drop"
      )

    n_lab <- max(1L, floor(0.15 * nrow(dfi)))
    dfi_low <- dplyr::slice_min(dfi, resid, n = n_lab, with_ties = FALSE)
    rng <- range(dfi$resid, na.rm = TRUE)
    yspan <- if (is.finite(diff(rng))) diff(rng) else 1
    ybase <- max(summ$mean + summ$se, na.rm = TRUE)
    ypos <- ybase + 0.10 * yspan
    lab <- as.character(stars_from_p(p_fdr))

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
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::labs(
        title = paste0(f, "  (", g2, " vs ", g1, ")"),
        x = group_col,
        y = "Covariate-adjusted residual (log10 intensity)"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    if (requireNamespace("ggpubr", quietly = TRUE) && nchar(lab) > 0) {
      ann <- tibble::tibble(group1 = g1, group2 = g2, y.position = ypos, label = lab)
      p <- p +
        ggplot2::expand_limits(y = ypos + 0.05 * yspan) +
        ggpubr::stat_pvalue_manual(
          ann,
          label = "label",
          y.position = "y.position",
          tip.length = 0.01,
          size = 5
        )
    }

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

    ggplot2::ggsave(
      file.path(out_dir, paste0("M3_pairwise_bar_", safe(g2), "_vs_", safe(g1), "_", safe(f), ".pdf")),
      p,
      width = 6.5,
      height = 4.2
    )
  }
  invisible(NULL)
}
