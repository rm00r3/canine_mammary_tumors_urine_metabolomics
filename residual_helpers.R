get_residuals <- function(X, meta, covars = NULL, group_col, control_label) {
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  if (!is.null(covars) && length(covars)) {
    for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])
  }
  frm <- stats::as.formula(
    paste0("~ ", group_col, if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else "")
  )
  design <- stats::model.matrix(frm, data = meta1)
  fit <- limma::lmFit(X1, design)
  fitted <- fit$coefficients %*% t(design)
  resid <- X1 - fitted
  list(resid = resid, samples = colnames(X1))
}

export_covariate_only_residuals <- function(X, meta, covariates, outfile) {
  stopifnot(is.matrix(X), ncol(X) == nrow(meta))
  if (!identical(colnames(X), rownames(meta))) {
    common <- intersect(colnames(X), rownames(meta))
    stopifnot(length(common) >= 3L)
    X <- X[, common, drop = FALSE]
    meta <- meta[common, , drop = FALSE]
  }

  for (v in covariates) {
    if (is.character(meta[[v]])) meta[[v]] <- factor(meta[[v]])
    if (is.factor(meta[[v]])) meta[[v]] <- droplevels(meta[[v]])
  }

  varying_covs <- covariates[vapply(covariates, function(v) {
    x <- meta[[v]]
    if (is.factor(x)) nlevels(x) > 1 else stats::sd(as.numeric(x), na.rm = TRUE) > 0
  }, logical(1L))]

  cc <- if (length(varying_covs)) stats::complete.cases(meta[, varying_covs, drop = FALSE]) else rep(TRUE, nrow(meta))
  X0 <- X[, cc, drop = FALSE]
  meta0 <- meta[cc, , drop = FALSE]

  design <- if (length(varying_covs)) {
    stats::model.matrix(stats::as.formula(paste0("~ ", paste(varying_covs, collapse = " + "))), data = meta0)
  } else {
    stats::model.matrix(~ 1, data = meta0)
  }

  if (ncol(design) > 1L) {
    nzv <- vapply(seq_len(ncol(design)), function(j) stats::sd(design[, j]) > 0, logical(1L))
    nzv[1] <- TRUE
    design <- design[, nzv, drop = FALSE]
  }
  qrD <- qr(design)
  if (qrD$rank < ncol(design)) design <- design[, qrD$pivot[1:qrD$rank], drop = FALSE]
  if (ncol(design) >= ncol(X0)) {
    warning("Design has >= number of samples; falling back to intercept-only.")
    design <- stats::model.matrix(~ 1, data = meta0)
  }

  fit_all <- limma::lmFit(X0, design)
  Res_all_cov <- limma::residuals.MArrayLM(fit_all, y = X0)
  utils::write.csv(Res_all_cov, outfile, row.names = TRUE)
  message("Wrote: ", outfile, "  (", nrow(Res_all_cov), " features x ", ncol(Res_all_cov), " samples)")
  Res_all_cov
}

normality_plots <- function(res_adj_norm, shap_p, shap_fdr) {
  df_norm <- tibble::tibble(
    feature = names(shap_p) %||% seq_along(shap_p),
    p = as.numeric(shap_p),
    q = as.numeric(shap_fdr)
  ) |>
    dplyr::filter(is.finite(p), is.finite(q))

  n_tested <- nrow(df_norm)
  n_p05 <- sum(df_norm$p < 0.05)
  n_q05 <- sum(df_norm$q < 0.05)

  p_hist <- ggplot2::ggplot(df_norm, ggplot2::aes(p)) +
    ggplot2::geom_histogram(bins = 40, boundary = 0, closed = "left") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "Shapiro-Wilk residual normality: p-value distribution",
      subtitle = sprintf("n=%d features; p<0.05: %d (%.1f%%)", n_tested, n_p05, 100 * n_p05 / n_tested),
      x = "p-value",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  q_hist <- ggplot2::ggplot(df_norm, ggplot2::aes(q)) +
    ggplot2::geom_histogram(bins = 40, boundary = 0, closed = "left") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "Benjamini-Hochberg FDR: q-value distribution",
      subtitle = sprintf("q<0.05: %d (%.1f%%)", n_q05, 100 * n_q05 / n_tested),
      x = "q-value",
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  p_ecdf <- ggplot2::ggplot(df_norm, ggplot2::aes(p)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(title = "ECDF of Shapiro p-values", x = "p-value", y = "F(p)") +
    ggplot2::theme_minimal(base_size = 12)

  q_ecdf <- ggplot2::ggplot(df_norm, ggplot2::aes(q)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(title = "ECDF of BH q-values", x = "q-value", y = "F(q)") +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave("supp_normality_p_hist.pdf", p_hist, width = 6.5, height = 4.2)
  ggplot2::ggsave("supp_normality_q_hist.pdf", q_hist, width = 6.5, height = 4.2)
  ggplot2::ggsave("supp_normality_p_ecdf.pdf", p_ecdf, width = 4.5, height = 4.2)
  ggplot2::ggsave("supp_normality_q_ecdf.pdf", q_ecdf, width = 4.5, height = 4.2)

  ord_all <- names(sort(shap_p, na.last = NA))
  n_feat <- sum(!is.na(shap_p))
  k_each <- if (n_feat < 200) 2 else if (n_feat < 1000) 3 else 4
  k_each <- min(k_each, 5)
  feat_worst <- head(ord_all, k_each)
  feat_best <- tail(ord_all, k_each)
  k_mid <- if (n_feat >= 1000) 2 else 1
  mid_idx <- round(seq(0.5, 0.5, length.out = k_mid) * length(ord_all))
  feat_mid <- ord_all[pmax(1, pmin(length(ord_all), mid_idx))]
  rep_feats <- unique(na.omit(c(feat_worst, feat_mid, feat_best)))
  rep_feats <- rep_feats[rep_feats %in% rownames(res_adj_norm$resid)]
  message("Representative features: ", paste(rep_feats, collapse = ", "))

  df_res <- dplyr::bind_rows(lapply(rep_feats, function(f) {
    tibble::tibble(feature = f, resid = as.numeric(res_adj_norm$resid[f, ]))
  }))
  if (!nrow(df_res)) return(invisible(NULL))
  df_res <- df_res |>
    dplyr::group_by(feature) |>
    dplyr::mutate(z = {
      s <- stats::sd(resid, na.rm = TRUE)
      m <- mean(resid, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      (resid - m) / s
    }) |>
    dplyr::ungroup()

  p_dens <- ggplot2::ggplot(df_res, ggplot2::aes(x = z)) +
    ggplot2::geom_density() +
    ggplot2::stat_function(fun = stats::dnorm, linetype = "dashed") +
    ggplot2::facet_wrap(~ feature, scales = "free_y") +
    ggplot2::labs(
      title = "Residual density (z-scaled) with N(0,1) overlay",
      x = "Residual (z)",
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave("residual_density_panels.pdf", p_dens, width = 10, height = 6)

  p_qq <- ggplot2::ggplot(df_res, ggplot2::aes(sample = z)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::facet_wrap(~ feature) +
    ggplot2::labs(
      title = "QQ plots of residuals (z-scaled)",
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave("residual_qq_panels.pdf", p_qq, width = 10, height = 6)
  invisible(NULL)
}
