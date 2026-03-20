safe_nlog10 <- function(p) {
  -log10(pmax(p, .Machine$double.xmin))
}

ensure_full_rank <- function(Z, protect_idx = integer(0L), verbose = TRUE) {
  keep <- unique(protect_idx[protect_idx %in% seq_len(ncol(Z))])
  for (j in setdiff(seq_len(ncol(Z)), keep)) {
    Zcand <- cbind(Z[, keep, drop = FALSE], Z[, j, drop = FALSE])
    if (qr(Zcand)$rank > ncol(Z[, keep, drop = FALSE])) keep <- c(keep, j)
  }
  if (verbose && length(keep) < ncol(Z)) {
    message("Dropped aliased columns: ",
            paste(setdiff(colnames(Z), colnames(Z)[keep]), collapse = ", "))
  }
  Z[, keep, drop = FALSE]
}

check_overlap <- function(meta, group_col, covars) {
  purrr::map_dfr(covars, function(v) {
    tab <- table(meta[[group_col]], meta[[v]], useNA = "no")
    as.data.frame.matrix(tab) |>
      tibble::rownames_to_column("group") |>
      tidyr::pivot_longer(-group, names_to = "level", values_to = "n") |>
      dplyr::mutate(
        covariate = v,
        zero_in_group = n == 0,
        case_only = (n > 0) & ave(n, level, FUN = function(x) sum(x > 0)) == 1
      )
  })
}

pretty_contrast_label <- function(txt) {
  sub("^X", "", gsub("_vs_.*$", "", txt))
}

per_contrast_pvals_hc3 <- function(X, meta, covars = NULL, group_col, control_label) {
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  if (!is.null(covars) && length(covars)) {
    for (v in covars) {
      if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])
    }
  }

  frm <- stats::as.formula(
    paste0(
      "~ ", group_col,
      if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""
    )
  )
  design <- stats::model.matrix(frm, data = meta1)

  levs <- levels(grp)
  others <- setdiff(levs, control_label)
  if (!length(others)) stop("Only control samples available.")

  coef_names <- colnames(design)
  coef_idx <- sapply(
    others,
    function(lv) which(coef_names == paste0(group_col, lv)),
    USE.NAMES = TRUE
  )

  do_one <- function(y) {
    fit_lm <- stats::lm(y ~ design - 1)
    V <- sandwich::vcovHC(fit_lm, type = "HC3")
    cf <- stats::coef(fit_lm)
    tibble::tibble(
      term = names(coef_idx),
      p_raw_hc3 = vapply(coef_idx, function(j) {
        if (is.na(cf[j])) return(NA_real_)
        se <- sqrt(V[j, j])
        if (!is.finite(se) || se == 0) return(NA_real_)
        2 * stats::pt(abs(cf[j] / se), df = nrow(design) - ncol(design), lower.tail = FALSE)
      }, numeric(1L))
    )
  }

  lst <- apply(X1, 1, do_one)
  out <- dplyr::bind_rows(lapply(seq_along(lst), function(i) {
    tibble::tibble(
      feature = rownames(X1)[i],
      Contrast = paste0(lst[[i]]$term, "_vs_", control_label),
      p_raw_hc3 = lst[[i]]$p_raw_hc3
    )
  }))
  out |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(p_fdr_hc3 = stats::p.adjust(p_raw_hc3, "fdr")) |>
    dplyr::ungroup()
}
