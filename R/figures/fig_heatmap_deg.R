# =============================================================================
#  fig_heatmap_deg.R — команда `fig_heatmap_deg`: тепловые карты значимых DEG
#  (из Pictures.R ч.2).
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, evs_dir])
box::use(R/core/celltypes[classify_cell_type, load_degs_classified, GROUP_ORDER])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "pheatmap", "dplyr", "tidyr", "RColorBrewer"))
  suppressMessages({ library(openxlsx); library(pheatmap); library(dplyr); library(tidyr); library(grDevices) })
  
  degs <- load_degs_classified(tissue, cc)
  heatmap_outdir <- file.path(evs_dir(tissue, cc), paste0(TAG[[cc]], " Heatmaps"))
  dir.create(heatmap_outdir, showWarnings = FALSE, recursive = TRUE)
  cat("[fig_heatmap_deg]", tissue, TAG[[cc]], "- DEG:", nrow(degs), "\n")
  
  filter_genes_for_heatmap <- function(data) {
    n_ct <- n_distinct(data$Cell_Type); min_ct <- min(HM_MIN_CT, n_ct)
    gene_ct_count <- data %>% dplyr::group_by(gene) %>%
      dplyr::summarise(n_cell_types = n_distinct(Cell_Type), min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(n_cell_types >= min_ct) %>% dplyr::arrange(min_padj) %>% head(HM_TOP_N)
    data %>% dplyr::filter(gene %in% gene_ct_count$gene)
  }
  build_heatmap_matrix <- function(data) {
    mat_long <- data %>% dplyr::select(gene, log2FoldChange, Cell_Type) %>%
      dplyr::group_by(gene, Cell_Type) %>%
      dplyr::summarise(log2FoldChange = mean(log2FoldChange), .groups = "drop") %>%
      pivot_wider(names_from = Cell_Type, values_from = log2FoldChange, values_fill = 0)
    mat <- as.matrix(mat_long[, -1]); rownames(mat) <- mat_long$gene; mat
  }
  order_columns_by_group <- function(mat) {
    ct_df <- data.frame(ct = colnames(mat), grp = sapply(colnames(mat), classify_cell_type), stringsAsFactors = FALSE)
    ct_df$grp <- factor(ct_df$grp, levels = GROUP_ORDER); ct_df <- ct_df[order(ct_df$grp, ct_df$ct), ]
    mat[, ct_df$ct, drop = FALSE]
  }
  make_annotation <- function(mat) {
    annotation_col <- data.frame(Group = sapply(colnames(mat), classify_cell_type), row.names = colnames(mat))
    annotation_col$Group <- factor(annotation_col$Group, levels = GROUP_ORDER); annotation_col
  }
  draw_heatmap <- function(data, filename_base, title_text) {
    if (nrow(data) == 0) { cat("  нет данных:", title_text, "\n"); return(invisible(NULL)) }
    if (n_distinct(data$Cell_Type) < 2) { cat("  <2 типов клеток, пропуск:", title_text, "\n"); return(invisible(NULL)) }
    filtered_data <- filter_genes_for_heatmap(data)
    if (nrow(filtered_data) == 0) { cat("  нет генов после фильтра:", title_text, "\n"); return(invisible(NULL)) }
    mat <- build_heatmap_matrix(filtered_data); mat <- order_columns_by_group(mat)
    if (nrow(mat) > 1 && ncol(mat) > 1) { row_vars <- apply(mat, 1, var); mat <- mat[row_vars > 0, , drop = FALSE] }
    if (nrow(mat) < 2) { cat("  мало генов:", title_text, "\n"); return(invisible(NULL)) }
    annotation_col <- make_annotation(mat)
    heatmap_colors <- colorRampPalette(c(HM_COL_LOW, HM_COL_MID, HM_COL_HIGH))(100)
    ann_colors <- list(Group = c(Glutamatergic = "#E41A1C", GABAergic = "#377EB8", Glia = "#4DAF4A", Microglia = "#984EA3"))
    h_adj <- max(HM_HEIGHT, nrow(mat) * 0.15 + 2)
    png(file.path(heatmap_outdir, paste0(filename_base, ".png")), width = HM_WIDTH, height = h_adj, units = "in", res = 300)
    pheatmap(mat, scale = HM_SCALE, color = heatmap_colors,
             clustering_distance_rows = HM_DISTANCE, clustering_distance_cols = HM_DISTANCE,
             clustering_method = HM_CLUSTER_METHOD, cluster_rows = nrow(mat) > 1, cluster_cols = ncol(mat) > 1,
             show_rownames = TRUE, show_colnames = TRUE, fontsize_row = HM_FONT_ROW, fontsize_col = HM_FONT_COL,
             border_color = NA, annotation_col = annotation_col, annotation_colors = ann_colors,
             main = title_text, display_numbers = FALSE)
    dev.off(); cat("  сохранена:", filename_base, ".png\n")
  }
  
  group_titles <- c(Glutamatergic = "Glutamatergic neurons", GABAergic = "GABAergic neurons",
                    Glia = "Glia", Microglia = "Microglia")
  draw_heatmap(degs, "Heatmap_All_CellTypes", paste(tissue, TAG[[cc]]))
  for (grp in GROUP_ORDER)
    draw_heatmap(degs %>% dplyr::filter(Group == grp), paste0("Heatmap_", grp), group_titles[grp])
  cat("  Heatmaps ->", heatmap_outdir, "\n")
}
