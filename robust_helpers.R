robust_fit_features <- function(X, meta, features, covars, group_col, control_label) {
  if (!length(features)) {
    return(tibble::tibble(
      feature = character(),
      Contrast = character(),
      estimate = numeric(),
      std.error = numeric(),
      statistic = numeric(),
      p.value = numeric(),
      p_fdr = numeric()
    ))
  }

  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[features, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]
  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])

  frm <- stats::as.formula(
    paste0("y ~ ", group_col, if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else "")
  )
  res <- lapply(rownames(X1), function(f) {
    y <- as.numeric(X1[f, ])
    fit <- robustbase::lmrob(frm, data = cbind(meta1, y = y))
    broom::tidy(fit) |>
      dplyr::filter(grepl(paste0("^", group_col), term)) |>
      dplyr::mutate(
        feature = f,
        Contrast = paste0(sub(paste0("^", group_col), "", term), "_vs_", control_label)
      ) |>
      dplyr::select(feature, Contrast, estimate, std.error, statistic, p.value)
  })
  out <- dplyr::bind_rows(res)
  if (!nrow(out)) {
    return(tibble::tibble(
      feature = character(),
      Contrast = character(),
      estimate = numeric(),
      std.error = numeric(),
      statistic = numeric(),
      p.value = numeric(),
      p_fdr = numeric()
    ))
  }
  out |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(p_fdr = stats::p.adjust(p.value, "fdr")) |>
    dplyr::ungroup()
}
