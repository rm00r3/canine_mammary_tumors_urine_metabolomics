stop_if_empty <- function(meta, group_col, control_label, covariates, min_per_class = 3) {
  grp <- droplevels(factor(meta[[group_col]]))
  tab <- table(grp, useNA = "ifany")
  message("Group counts (all rows):")
  print(tab)
  if (length(covariates)) {
    cc_cov <- stats::complete.cases(meta[, covariates, drop = FALSE])
    tab_cc <- table(grp[cc_cov], useNA = "ifany")
    message("Group counts after complete-cases over covariates:")
    print(tab_cc)
    if (any(tab_cc[names(tab_cc) != control_label] < min_per_class) ||
        (tab_cc[control_label] %||% 0) < min_per_class) {
      message("Some contrasts may have too few complete cases for stable ORs.")
    }
  }
  invisible(NULL)
}

get_coef_with_fallback <- function(df, frm, coef_name = "zexp") {
  fit_mle <- try(stats::glm(frm, data = df, family = stats::binomial(), control = stats::glm.control(maxit = 50)), silent = TRUE)
  if (!inherits(fit_mle, "try-error")) {
    est <- NA_real_
    se <- NA_real_
    p <- NA_real_
    src <- NA_character_

    ct_hc3 <- try({
      V <- sandwich::vcovHC(fit_mle, type = "HC3")
      lmtest::coeftest(fit_mle, vcov. = V)
    }, silent = TRUE)
    if (!inherits(ct_hc3, "try-error") && coef_name %in% rownames(ct_hc3)) {
      est <- unname(ct_hc3[coef_name, "Estimate"])
      se <- unname(ct_hc3[coef_name, "Std. Error"])
      p <- as.numeric(ct_hc3[coef_name, grep("^Pr", colnames(ct_hc3))[1]])
      src <- "HC3"
    } else {
      sm <- try(summary(fit_mle)$coefficients, silent = TRUE)
      if (!inherits(sm, "try-error") && coef_name %in% rownames(sm)) {
        est <- unname(sm[coef_name, "Estimate"])
        se <- unname(sm[coef_name, "Std. Error"])
        p <- as.numeric(sm[coef_name, grep("^Pr", colnames(sm))[1]])
        src <- "MLE"
      }
    }
    if (is.finite(est) && is.finite(se) && se > 0) {
      return(list(beta = est, se = se, p = p, source = src))
    }
  }

  if (requireNamespace("logistf", quietly = TRUE)) {
    fit_f <- try(logistf::logistf(frm, data = df), silent = TRUE)
    if (!inherits(fit_f, "try-error") && !is.null(fit_f$coefmat) &&
        coef_name %in% rownames(fit_f$coefmat)) {
      row <- fit_f$coefmat[coef_name, , drop = FALSE]
      return(list(
        beta = as.numeric(row[, "coef"]),
        se = as.numeric(row[, "se(coef)"]),
        p = as.numeric(row[, "prob"]),
        source = "Firth"
      ))
    }
  }

  if (requireNamespace("brglm2", quietly = TRUE)) {
    fit_b <- try(stats::glm(
      frm,
      data = df,
      family = stats::binomial(),
      method = brglm2::brglmFit,
      type = "AS_mean",
      control = stats::glm.control(maxit = 50)
    ), silent = TRUE)
    if (!inherits(fit_b, "try-error")) {
      sm <- try(summary(fit_b)$coefficients, silent = TRUE)
      if (!inherits(sm, "try-error") && coef_name %in% rownames(sm)) {
        est <- unname(sm[coef_name, "Estimate"])
        se <- unname(sm[coef_name, "Std. Error"])
        z <- est / se
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        return(list(beta = est, se = se, p = p, source = "brglm2"))
      }
    }
  }
  list(beta = NA_real_, se = NA_real_, p = NA_real_, source = NA_character_)
}

per_contrast_or_hc3 <- function(
  X, meta, features, group_col, control_label, covariates, min_per_class = 3
) {
  grp_all <- droplevels(factor(meta[[group_col]]))
  others <- setdiff(levels(grp_all), control_label)
  if (!length(others)) stop("No non-control levels in ", group_col)
  out <- list()

  for (lv in others) {
    keep_rows <- meta[[group_col]] %in% c(control_label, lv)
    metas <- droplevels(meta[keep_rows, , drop = FALSE])
    y_full <- as.integer(metas[[group_col]] == lv)

    for (v in covariates) if (is.character(metas[[v]])) metas[[v]] <- factor(metas[[v]])
    rhs <- if (length(covariates)) paste(c("zexp", covariates), collapse = " + ") else "zexp"
    frm <- stats::as.formula(paste0("y ~ ", rhs))

    res_lv <- lapply(features, function(f) {
      x <- as.numeric(X[f, keep_rows])
      cc <- if (length(covariates)) {
        stats::complete.cases(metas[, covariates, drop = FALSE]) & is.finite(x)
      } else {
        is.finite(x)
      }
      y <- y_full[cc]

      if (length(unique(y)) != 2L) {
        return(tibble::tibble(
          feature = f,
          Contrast = paste0(lv, "_vs_", control_label),
          OR = NA_real_,
          OR_low = NA_real_,
          OR_high = NA_real_,
          p = NA_real_,
          q = NA_real_,
          se_source = NA_character_,
          n = length(y),
          n_case = sum(y == 1L),
          n_ctrl = sum(y == 0L),
          sd_all = NA_real_,
          sd_case = NA_real_,
          sd_ctrl = NA_real_,
          overlap_flag = NA,
          reason = "monoclass_after_CC"
        ))
      }

      z <- scale(x[cc])[, 1]
      sd_all <- stats::sd(z)
      sd_case <- stats::sd(z[y == 1L])
      sd_ctrl <- stats::sd(z[y == 0L])
      n_case <- sum(y == 1L)
      n_ctrl <- sum(y == 0L)

      if (n_case < min_per_class || n_ctrl < min_per_class) {
        return(tibble::tibble(
          feature = f,
          Contrast = paste0(lv, "_vs_", control_label),
          OR = NA_real_,
          OR_low = NA_real_,
          OR_high = NA_real_,
          p = NA_real_,
          q = NA_real_,
          se_source = NA_character_,
          n = length(y),
          n_case = n_case,
          n_ctrl = n_ctrl,
          sd_all = sd_all,
          sd_case = sd_case,
          sd_ctrl = sd_ctrl,
          overlap_flag = NA,
          reason = "too_few_per_class"
        ))
      }
      if (!is.finite(sd_all) || sd_all == 0) {
        return(tibble::tibble(
          feature = f,
          Contrast = paste0(lv, "_vs_", control_label),
          OR = NA_real_,
          OR_low = NA_real_,
          OR_high = NA_real_,
          p = NA_real_,
          q = NA_real_,
          se_source = NA_character_,
          n = length(y),
          n_case = n_case,
          n_ctrl = n_ctrl,
          sd_all = sd_all,
          sd_case = sd_case,
          sd_ctrl = sd_ctrl,
          overlap_flag = NA,
          reason = "zero_variance"
        ))
      }

      zc <- z[y == 0L]
      za <- z[y == 1L]
      overlap_flag <- !(max(zc, na.rm = TRUE) < min(za, na.rm = TRUE) ||
                        max(za, na.rm = TRUE) < min(zc, na.rm = TRUE))
      df <- cbind(metas[cc, , drop = FALSE], y = y, zexp = as.numeric(z))
      cf <- get_coef_with_fallback(df, frm, coef_name = "zexp")

      if (!is.finite(cf$beta) || !is.finite(cf$se) || cf$se <= 0) {
        return(tibble::tibble(
          feature = f,
          Contrast = paste0(lv, "_vs_", control_label),
          OR = NA_real_,
          OR_low = NA_real_,
          OR_high = NA_real_,
          p = NA_real_,
          q = NA_real_,
          se_source = cf$source,
          n = length(y),
          n_case = n_case,
          n_ctrl = n_ctrl,
          sd_all = sd_all,
          sd_case = sd_case,
          sd_ctrl = sd_ctrl,
          overlap_flag = overlap_flag,
          reason = "coef_or_se_na"
        ))
      }

      tibble::tibble(
        feature = f,
        Contrast = paste0(lv, "_vs_", control_label),
        OR = exp(cf$beta),
        OR_low = exp(cf$beta - 1.96 * cf$se),
        OR_high = exp(cf$beta + 1.96 * cf$se),
        p = as.numeric(cf$p),
        q = NA_real_,
        se_source = cf$source,
        n = length(y),
        n_case = n_case,
        n_ctrl = n_ctrl,
        sd_all = sd_all,
        sd_case = sd_case,
        sd_ctrl = sd_ctrl,
        overlap_flag = overlap_flag,
        reason = NA_character_
      )
    })
    out[[lv]] <- dplyr::bind_rows(res_lv)
  }

  ans <- dplyr::bind_rows(out)
  if (!nrow(ans)) return(ans)
  ans |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(q = stats::p.adjust(p, method = "fdr")) |>
    dplyr::ungroup()
}

case_control_or_hc3 <- function(
  X, meta, features, group_col, case_label = "Case", control_label = "Control", covariates
) {
  keep <- meta[[group_col]] %in% c(control_label, case_label)
  X1 <- X[features, keep, drop = FALSE]
  meta1 <- meta[keep, , drop = FALSE]
  y <- as.integer(meta1[[group_col]] == case_label)
  for (v in covariates) if (is.character(meta1[[v]])) meta1[[v]] <- factor(meta1[[v]])

  frm <- stats::as.formula(
    paste0("y ~ zexp", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else "")
  )
  rows <- lapply(rownames(X1), function(f) {
    x <- as.numeric(X1[f, ])
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      return(tibble::tibble(
        feature = f, OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_,
        p = NA_real_, se_source = NA_character_
      ))
    }
    df <- cbind(meta1, y = y, zexp = (x - m) / s)
    cf <- get_coef_with_fallback(df, frm, "zexp")
    if (!is.finite(cf$beta) || !is.finite(cf$se) || cf$se <= 0) {
      return(tibble::tibble(
        feature = f, OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_,
        p = as.numeric(cf$p), se_source = cf$source
      ))
    }
    tibble::tibble(
      feature = f,
      OR = exp(cf$beta),
      OR_low = exp(cf$beta - 1.96 * cf$se),
      OR_high = exp(cf$beta + 1.96 * cf$se),
      p = as.numeric(cf$p),
      se_source = cf$source
    )
  })
  out <- dplyr::bind_rows(rows)
  out$p <- as.numeric(out$p)
  out |>
    dplyr::mutate(
      q = stats::p.adjust(p, method = "fdr"),
      Contrast = paste0(case_label, "_vs_", control_label)
    ) |>
    dplyr::select(feature, Contrast, OR, OR_low, OR_high, p, q, se_source)
}
