align_expression_and_meta <- function(X, meta, group_col) {
  X <- as.matrix(X)
  ids_meta <- rownames(meta)
  if (is.null(ids_meta)) {
    stop("Metadata rownames are required for sample alignment.")
  }

  ov_col <- intersect(colnames(X), ids_meta)
  ov_row <- intersect(rownames(X), ids_meta)

  if (length(ov_col) >= 2L) {
    message(sprintf("Detected samples in COLUMNS (%d overlap).", length(ov_col)))
    X <- X[, ov_col, drop = FALSE]
    meta <- meta[ov_col, , drop = FALSE]
  } else if (length(ov_row) >= 2L) {
    message(sprintf("Detected samples in ROWS (%d overlap) - transposing.", length(ov_row)))
    X <- t(X)
    X <- X[, ov_row, drop = FALSE]
    meta <- meta[ov_row, , drop = FALSE]
  } else {
    stop("No sample ID overlap between expression matrix and metadata.")
  }

  keep_feat <- rowSums(is.finite(X)) >= 2L
  X <- X[keep_feat, , drop = FALSE]

  message(sprintf("Features: %d | Samples: %d", nrow(X), ncol(X)))
  if (!group_col %in% colnames(meta)) {
    stop("Missing group column in metadata: ", group_col)
  }
  message("Group counts AFTER alignment:")
  print(table(meta[[group_col]], useNA = "ifany"))

  list(X = X, meta = meta)
}
