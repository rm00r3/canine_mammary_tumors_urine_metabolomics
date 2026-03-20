export_adjusted_significant_heatmap <- function(
  X, meta, pc, group_col, alpha_sig, filename, stop_if_empty = TRUE, symmetric_palette = FALSE
) {
  stopifnot(is.data.frame(pc), is.matrix(X), ncol(X) == nrow(meta), group_col %in% colnames(meta))
  sig_feats <- pc |>
    dplyr::group_by(feature) |>
    dplyr::summarise(best_q = suppressWarnings(min(p_adj_fdr_hc3, na.rm = TRUE)), .groups = "drop") |>
    dplyr::filter(is.finite(best_q) & best_q < alpha_sig) |>
    dplyr::arrange(best_q) |>
    dplyr::pull(feature)

  if (!length(sig_feats)) {
    msg <- paste0("No adjusted-significant features at FDR < ", alpha_sig)
    if (stop_if_empty) stop(msg, call. = FALSE)
    message(msg, " - skipping heatmap.")
    return(invisible(NULL))
  }

  Xsig <- X[intersect(sig_feats, rownames(X)), , drop = FALSE]
  if (!nrow(Xsig)) {
    if (stop_if_empty) stop("Selected features not found in X rownames.", call. = FALSE)
    message("Selected features not found in X rownames; skipping heatmap.")
    return(invisible(NULL))
  }

  Xsig[!is.finite(Xsig)] <- NA
  mu <- rowMeans(Xsig, na.rm = TRUE)
  sdv <- apply(Xsig, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xz <- sweep(sweep(Xsig, 1, mu, "-"), 1, sdv, "/")
  Xz[!is.finite(Xz)] <- 0

  grp_vec <- as.character(meta[colnames(Xz), group_col, drop = TRUE])
  grp_vec[is.na(grp_vec) | grp_vec == ""] <- "Missing"
  ord <- order(grp_vec)
  Xz <- Xz[, ord, drop = FALSE]
  grp_vec <- grp_vec[ord]
  ann_df <- data.frame(setNames(list(factor(grp_vec)), group_col), row.names = colnames(Xz))

  if (symmetric_palette) {
    rng <- range(Xz, finite = TRUE)
    v <- max(abs(rng))
    if (!is.finite(v) || v == 0) v <- 1e-6
    brks <- seq(-v, v, length.out = 101)
    cols <- grDevices::colorRampPalette(c("#5858f4", "#dcdcff", "#ff0000"))(length(brks) - 1)
  } else {
    rng <- range(Xz, finite = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
      eps <- ifelse(is.finite(rng[1]) && rng[1] != 0, abs(rng[1]) * 1e-6, 1e-6)
      rng <- c(-eps, eps)
    }
    brks <- seq(rng[1], rng[2], length.out = 101)
    cols <- grDevices::colorRampPalette(c("#214478", "#ffffff", "#781f19"))(length(brks) - 1)
  }

  pheatmap::pheatmap(
    Xz,
    scale = "none",
    cluster_rows = nrow(Xz) > 1,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = ann_df,
    color = cols,
    breaks = brks,
    border_color = NA,
    filename = filename,
    width = 7.5,
    height = if (symmetric_palette) 1 + 0.12 * nrow(Xz) else 4 + 0.12 * nrow(Xz)
  )
  message("Heatmap drawn for ", nrow(Xz), " adjusted-significant features (FDR < ", alpha_sig, ").")
  invisible(filename)
}
