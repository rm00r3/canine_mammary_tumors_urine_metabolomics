influence_diag_strict <- function(X, meta, covars, group_col, control_label) {
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])

  frm <- stats::as.formula(paste0("~ ", group_col, " + ", paste(covars, collapse = " + ")))
  design <- stats::model.matrix(frm, data = meta1)
  n <- nrow(design)
  p <- ncol(design)
  lev_thresh_strict <- 3 * p / n
  cook_thresh_strict <- 1

  H <- design %*% solve(crossprod(design)) %*% t(design)
  max_hat_all <- max(diag(H), na.rm = TRUE)
  do_one <- function(y) {
    fit <- stats::lm(y ~ design - 1)
    c(max_hat = max_hat_all, max_cook = max(stats::cooks.distance(fit), na.rm = TRUE))
  }
  M <- t(apply(X1, 1, do_one))
  tibble::as_tibble(M, rownames = "feature") |>
    dplyr::mutate(
      n = n,
      p = p,
      lev_thresh_strict = lev_thresh_strict,
      cook_thresh_strict = cook_thresh_strict,
      flag_hat_strict = max_hat > lev_thresh_strict,
      flag_cook_strict = max_cook > cook_thresh_strict
    )
}

sample_leverage_table <- function(meta, covars, group_col, control_label) {
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])
  frm <- stats::as.formula(paste0("~ ", group_col, " + ", paste(covars, collapse = " + ")))
  design <- stats::model.matrix(frm, data = meta1)
  n <- nrow(design)
  p <- ncol(design)

  H <- design %*% solve(crossprod(design)) %*% t(design)
  lev <- diag(H)
  lev_thresh_strict <- 3 * p / n
  tibble::tibble(
    sample = rownames(meta1),
    leverage = lev,
    flag_leverage = lev > lev_thresh_strict,
    lev_thresh_strict = lev_thresh_strict
  ) |>
    dplyr::arrange(dplyr::desc(leverage))
}

sample_cook_counts <- function(X, meta, covars, group_col, control_label, cook_thresh_strict = 1) {
  need <- c(group_col, covars %||% character())
  cc <- stats::complete.cases(meta[, need, drop = FALSE])
  X1 <- X[, cc, drop = FALSE]
  meta1 <- meta[cc, , drop = FALSE]

  grp <- stats::relevel(factor(meta1[[group_col]]), ref = control_label)
  meta1[[group_col]] <- grp
  for (v in covars) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])
  frm <- stats::as.formula(paste0("~ ", group_col, " + ", paste(covars, collapse = " + ")))
  design <- stats::model.matrix(frm, data = meta1)
  n_feat <- nrow(X1)
  counts <- setNames(integer(nrow(design)), rownames(meta1))

  for (i in seq_len(n_feat)) {
    y <- as.numeric(X1[i, ])
    fit <- stats::lm(y ~ design - 1)
    cd <- stats::cooks.distance(fit)
    counts <- counts + as.integer(cd > cook_thresh_strict)
  }

  tibble::tibble(
    sample = names(counts),
    cook_count = as.integer(counts),
    cook_prop = counts / n_feat
  ) |>
    dplyr::arrange(dplyr::desc(cook_prop))
}

run_na_robust_cook_summary <- function(
  X, meta, covariates_adj, group_col, out_prefix,
  default_diag_summary_file = paste0(out_prefix, "_diag_summary.csv"),
  also_write_diagnostics_summary = FALSE
) {
  OUT_PREFIX <- out_prefix
  COOK_THRESH <- NA_real_
  SUBJ_PROP_FLAG <- 0.01

  m2 <- meta
  for (v in covariates_adj) if (is.character(m2[[v]])) m2[[v]] <- factor(m2[[v]])
  frm <- stats::as.formula(
    paste0("~ ", group_col,
           if (length(covariates_adj)) paste0(" + ", paste(covariates_adj, collapse = " + ")) else "")
  )
  Z_full <- stats::model.matrix(frm, data = m2)

  qrD <- qr(Z_full)
  if (qrD$rank < ncol(Z_full)) {
    Z_full <- Z_full[, qrD$pivot[1:qrD$rank], drop = FALSE]
  }
  p <- ncol(Z_full)
  if (!is.finite(COOK_THRESH) || is.na(COOK_THRESH)) {
    COOK_THRESH <- 4 / max(1, nrow(Z_full) - p)
  }

  ids <- rownames(m2)
  stopifnot(identical(colnames(X), ids))
  subj_hits <- setNames(integer(length(ids)), ids)
  subj_tested <- setNames(integer(length(ids)), ids)
  feat_hits <- integer(nrow(X))
  names(feat_hits) <- rownames(X)

  for (f in rownames(X)) {
    y <- as.numeric(X[f, ])
    cc <- is.finite(y) & stats::complete.cases(Z_full)
    n_cc <- sum(cc)
    if (n_cc <= p + 1) next
    Z <- Z_full[cc, , drop = TRUE]
    yy <- y[cc]
    fit <- try(stats::lm(yy ~ Z - 1), silent = TRUE)
    if (inherits(fit, "try-error")) next
    cd <- try(stats::cooks.distance(fit), silent = TRUE)
    if (inherits(cd, "try-error") || !length(cd)) next
    s_ids <- ids[cc]
    subj_tested[s_ids] <- subj_tested[s_ids] + 1L
    flagged <- cd > COOK_THRESH & is.finite(cd)
    if (any(flagged)) {
      feat_hits[f] <- 1L
      subj_hits[s_ids[flagged]] <- subj_hits[s_ids[flagged]] + 1L
    }
  }

  cooks_subj_df <- tibble::tibble(
    SampleID = names(subj_hits),
    cook_count = as.integer(subj_hits),
    tested_in = as.integer(subj_tested),
    cook_prop = ifelse(tested_in > 0, cook_count / tested_in, NA_real_),
    flag_high = is.finite(cook_prop) & (cook_prop >= SUBJ_PROP_FLAG)
  )
  cooks_feat_df <- tibble::tibble(
    feature = names(feat_hits),
    any_flagged = feat_hits > 0L,
    n_flagged = as.integer(feat_hits > 0L),
    prop_flagged = NA_real_
  )

  sum_cook_subj <- cooks_subj_df |>
    dplyr::summarise(
      units_tested = sum(tested_in > 0),
      n_flagged = sum(flag_high, na.rm = TRUE),
      pct_flagged = if (units_tested > 0) 100 * n_flagged / units_tested else 0
    ) |>
    dplyr::mutate(
      unit = "subjects",
      threshold_rule = paste0("Cook's > ", format(COOK_THRESH, digits = 3),
                              "; subject flagged if prop >= ", format(SUBJ_PROP_FLAG, digits = 3)),
      action_taken = "None (no exclusions); sensitivity unchanged"
    )
  sum_cook_feat <- cooks_feat_df |>
    dplyr::summarise(
      units_tested = nrow(cooks_feat_df),
      n_flagged = sum(any_flagged, na.rm = TRUE),
      pct_flagged = if (units_tested > 0) 100 * n_flagged / units_tested else 0
    ) |>
    dplyr::mutate(
      unit = "features",
      threshold_rule = paste0("Cook's > ", format(COOK_THRESH, digits = 3), " in >=1 subject"),
      action_taken = "None (no exclusions); sensitivity unchanged"
    )

  diag_summary <- dplyr::bind_rows(sum_cook_feat, sum_cook_subj) |>
    dplyr::mutate(Model = OUT_PREFIX, Test = "Cook's influence") |>
    dplyr::select(
      Model, Test, unit, units_tested, n_flagged, pct_flagged,
      threshold_rule, action_taken
    )

  readr::write_csv(cooks_subj_df, paste0(OUT_PREFIX, "_cooks_subject_counts.csv"))
  readr::write_csv(cooks_feat_df, paste0(OUT_PREFIX, "_cooks_feature_summary.csv"))
  readr::write_csv(diag_summary, default_diag_summary_file)
  if (also_write_diagnostics_summary) {
    readr::write_csv(diag_summary, paste0(OUT_PREFIX, "_diagnostics_summary.csv"))
  }
  invisible(list(
    cooks_subject = cooks_subj_df,
    cooks_feature = cooks_feat_df,
    diag_summary = diag_summary
  ))
}

write_diag_readme <- function(prefix, summary_file = paste0(prefix, "_diag_summary.csv")) {
  fn_sum <- summary_file
  fn_subj <- paste0(prefix, "_cooks_subject_counts.csv")
  fn_feat <- paste0(prefix, "_cooks_feature_summary.csv")
  df <- try(readr::read_csv(fn_sum, show_col_types = FALSE), silent = TRUE)
  getrow <- function(u, col) {
    if (inherits(df, "try-error") || !nrow(df)) return(NA)
    v <- df[df$unit == u, col, drop = TRUE]
    if (!length(v)) NA else v[1]
  }
  lines <- c(
    paste0("# Diagnostics summary for ", prefix),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    "",
    "Files produced:",
    paste0(" - ", fn_subj),
    paste0(" - ", fn_feat),
    paste0(" - ", fn_sum),
    "",
    "Cook's influence - subjects:",
    sprintf(" - Units tested: %s; Flagged: %s (%.2f%%)",
            getrow("subjects", "units_tested") %||% 0,
            getrow("subjects", "n_flagged") %||% 0,
            getrow("subjects", "pct_flagged") %||% 0),
    "Cook's influence - features:",
    sprintf(" - Units tested: %s; Flagged: %s (%.2f%%)",
            getrow("features", "units_tested") %||% 0,
            getrow("features", "n_flagged") %||% 0,
            getrow("features", "pct_flagged") %||% 0),
    "",
    "Thresholds:",
    paste0(" - Subject rule: ", getrow("subjects", "threshold_rule") %||% "n/a"),
    paste0(" - Feature rule: ", getrow("features", "threshold_rule") %||% "n/a"),
    "",
    "Action taken:",
    paste0(" - ", if (!inherits(df, "try-error") && nrow(df)) unique(df$action_taken)[1] else "None (no exclusions)")
  )
  writeLines(lines, con = paste0(prefix, "_diagnostics_README.txt"))
  invisible(paste0(prefix, "_diagnostics_README.txt"))
}
