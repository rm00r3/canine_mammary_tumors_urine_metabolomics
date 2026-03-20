permute_FL_pvals <- function(
  X, meta, features, covars, group_col, control_label,
  B = 1000, seed = 1, robust_projection = FALSE
) {
  set.seed(seed)
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[features, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])

  Fm <- stats::as.formula(
    paste0("~ ", group_col, if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else "")
  )
  Rm <- if (length(covars)) {
    stats::as.formula(paste0("~ ", paste(covars, collapse = " + ")))
  } else {
    ~ 1
  }
  Xf <- stats::model.matrix(Fm, meta1)
  Xr <- stats::model.matrix(Rm, meta1)

  grp_cols <- setdiff(colnames(Xf), colnames(Xr))
  grp_idx <- match(grp_cols, colnames(Xf))
  levs <- setdiff(levels(grp), control_label)

  n <- nrow(Xf)
  pF <- ncol(Xf)
  coef_tstat <- function(y, j) {
    fit <- stats::lm.fit(Xf, y)
    res <- y - Xf %*% fit$coefficients
    sigma2 <- sum(res^2) / max(1, n - pF)
    XtXinv <- tryCatch(solve(crossprod(Xf)), error = function(e) MASS::ginv(crossprod(Xf)))
    se <- sqrt(sigma2 * XtXinv[j, j])
    beta <- fit$coefficients[j]
    if (!is.finite(se) || se == 0) return(NA_real_)
    beta / se
  }

  if (robust_projection) {
    XtX <- crossprod(Xr)
    XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
    P_r <- Xr %*% XtX_inv %*% t(Xr)
    yhatR <- X1 %*% P_r
    residR <- X1 - yhatR
  } else {
    fitR_coef <- stats::lm.fit(Xr, t(X1))$coefficients
    yhatR <- t(Xr %*% fitR_coef)
    residR <- X1 - yhatR
  }

  t_obs <- lapply(seq_along(grp_idx), function(k) {
    j <- grp_idx[k]
    apply(X1, 1, coef_tstat, j = j)
  })
  ge_counts <- lapply(t_obs, function(x) integer(length(x)))

  for (b in seq_len(B)) {
    perm <- sample.int(n)
    Xb <- yhatR + residR[, perm, drop = FALSE]
    tb <- lapply(seq_along(grp_idx), function(k) {
      j <- grp_idx[k]
      apply(Xb, 1, coef_tstat, j = j)
    })
    for (k in seq_along(grp_idx)) {
      ge_counts[[k]] <- ge_counts[[k]] + as.integer(abs(tb[[k]]) >= abs(t_obs[[k]]))
    }
  }

  out <- lapply(seq_along(grp_idx), function(k) {
    tibble::tibble(
      feature = rownames(X1),
      Contrast = paste0(levs[k], "_vs_", control_label),
      t_obs = t_obs[[k]],
      p_perm = (ge_counts[[k]] + 1) / (B + 1)
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(p_perm_fdr = stats::p.adjust(p_perm, "fdr")) |>
    dplyr::ungroup()
  out
}
