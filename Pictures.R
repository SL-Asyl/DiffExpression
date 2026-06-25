# --- Пути ---
degs_file <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx"
heatmap_outdir <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13 Heatmaps"
volcano_outdir <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13 VolcanoPlots"
heatmap_txt_dir <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13 Heatmaps"

# --- Heatmap настройки ---
heatmap_width <- 6
heatmap_height <- 12
heatmap_scale <- "row"
heatmap_show_rownames <- TRUE
heatmap_fontsize_row <- 8
heatmap_fontsize_col <- 8
heatmap_border <- FALSE
heatmap_cluster_method <- "ward.D"  # "complete", "ward.D", "ward.D2"
heatmap_distance <- "euclidean"
heatmap_color_low <- "#0000FF"
heatmap_color_mid <- "#FFFFFF"
heatmap_color_high <- "#FF0000"
heatmap_remove_identical <- TRUE
heatmap_top_n <- 50
heatmap_min_cell_types <- 2

# --- Volcano plot настройки ---
comparison_label <- "K18-ACE2-KI + FS vs K18-ACE2-KI + FS + par"
figure_width <- 8
figure_height <- 8
p_min <- 1e-100
fdr_cutoff <- 0.05
fc_cutoff <- 2
point_size <- 80
point_shape <- 16
point_alpha <- 0.7
show_grid <- TRUE
ticks_fontsize <- 18
axis_label_fontsize <- 24
legend_fontsize <- 14
label_fontsize <- 14
label_color <- "#FF0000"
n_label_genes <- 20
xmin <- NULL
xmax <- NULL
ymin <- -0.1
ymax <- NULL

color_up <- "#F39B7F"
color_down <- "#4DBBD5"
color_ns <- "#A6A6A6"

# --- УСТАНОВКА И ЗАГРУЗКА ПАКЕТОВ ---

required_packages <- c("openxlsx", "pheatmap", "ggplot2", "ggrepel",
                        "dplyr", "tidyr", "RColorBrewer", "grDevices")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(openxlsx)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(grDevices)

# --- ГРУППИРОВКА КЛЕТОЧНЫХ ТИПОВ ---

classify_cell_type <- function(ct) {
  if (grepl("[_-]Glut$", ct)) return("Glutamatergic")
  if (grepl("[_-]Gaba$", ct)) return("GABAergic")
  if (grepl("[_-]Inh", ct))   return("GABAergic")
  
  # кортикальные глутаматергические проекционные типы (любой суффикс после слоя)
  if (grepl("^L[0-9]", ct) || grepl("[_-]CTX$", ct)) return("Glutamatergic")
  if (grepl("[_-](IT|ET|PT|CT|NP)([_-]CTX)?$", ct)) return("Glutamatergic")
  
  # известные GABAergic подклассы без суффикса _Gaba
  if (ct %in% c("Pvalb", "Lamp5", "Vip", "Sncg", "Sst")) return("GABAergic")
  
  if (ct == "Microglia" || ct == "Microglia_NN") return("Microglia")
  if (ct %in% c("Astro", "Oligo", "OPC", "Endo", "Peri") ||
      grepl("[_-]NN$", ct)) return("Glia")
  
  return("Other")
}

group_order <- c("Glutamatergic", "GABAergic", "Glia", "Microglia")

# --- ЧТЕНИЕ ДАННЫХ ---
cell_type_clean_pattern <- "^10_vs_13_"
degs <- read.xlsx(degs_file)
degs$Cell_Type <- gsub(cell_type_clean_pattern, "", degs$Cell_Type)
degs$Group <- sapply(degs$Cell_Type, classify_cell_type)

cat("Загружено", nrow(degs), "DEGs из", length(unique(degs$Cell_Type)),
    "клеточных типов\n")

# --- ТЕПЛОВЫЕ КАРТЫ (HEATMAPS) ---

dir.create(heatmap_outdir, showWarnings = FALSE, recursive = TRUE)

filter_genes_for_heatmap <- function(data) {
  n_ct <- n_distinct(data$Cell_Type)
  min_ct <- min(heatmap_min_cell_types, n_ct)

  gene_ct_count <- data %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(n_cell_types = n_distinct(Cell_Type),
              min_padj = min(padj, na.rm = TRUE),
              .groups = "drop") %>%
    dplyr::filter(n_cell_types >= min_ct) %>%
    dplyr::arrange(min_padj) %>%
    head(heatmap_top_n)

  data %>% dplyr::filter(gene %in% gene_ct_count$gene)
}

build_heatmap_matrix <- function(data) {
  mat_long <- data %>%
    dplyr::select(gene, log2FoldChange, Cell_Type) %>%
    dplyr::group_by(gene, Cell_Type) %>%
    dplyr::summarise(log2FoldChange = mean(log2FoldChange), .groups = "drop") %>%
    pivot_wider(names_from = Cell_Type, values_from = log2FoldChange, values_fill = 0)

  mat <- as.matrix(mat_long[, -1])
  rownames(mat) <- mat_long$gene
  mat
}

order_columns_by_group <- function(mat) {
  cell_types <- colnames(mat)
  ct_groups <- sapply(cell_types, classify_cell_type)
  ct_df <- data.frame(ct = cell_types, grp = ct_groups, stringsAsFactors = FALSE)
  ct_df$grp <- factor(ct_df$grp, levels = group_order)
  ct_df <- ct_df[order(ct_df$grp, ct_df$ct), ]
  mat[, ct_df$ct, drop = FALSE]
}

make_annotation <- function(mat) {
  ct_groups <- sapply(colnames(mat), classify_cell_type)
  annotation_col <- data.frame(Group = ct_groups, row.names = colnames(mat))
  annotation_col$Group <- factor(annotation_col$Group, levels = group_order)
  annotation_col
}

save_heatmap_txt <- function(mat, filename) {
  ct_groups <- sapply(colnames(mat), classify_cell_type)

  header_group <- paste(c("group", ct_groups), collapse = "\t")
  header_sample <- paste(c("sample", colnames(mat)), collapse = "\t")

  lines <- c(header_group, header_sample)
  for (i in seq_len(nrow(mat))) {
    line <- paste(c(rownames(mat)[i], mat[i, ]), collapse = "\t")
    lines <- c(lines, line)
  }

  fpath <- file.path(heatmap_txt_dir, filename)
  writeLines(lines, fpath)
  cat("  TXT файл сохранён:", fpath, "\n")
}

draw_heatmap <- function(data, filename_base, title_text) {
  if (nrow(data) == 0) {
    cat("  Нет данных для", title_text, "- пропускаем\n")
    return(invisible(NULL))
  }
  
  # --- НОВОЕ: пропускаем группы с одним типом клеток ---
  n_ct_in_group <- n_distinct(data$Cell_Type)
  if (n_ct_in_group < 2) {
    cat("  Пропускаем", title_text, "- только", n_ct_in_group, "тип(а) клеток (heatmap не строится)\n")
    return(invisible(NULL))
  }

  filtered_data <- filter_genes_for_heatmap(data)
  if (nrow(filtered_data) == 0) {
    cat("  После фильтрации не осталось генов для", title_text, "\n")
    return(invisible(NULL))
  }

  cat("  Отобрано генов:", n_distinct(filtered_data$gene),
      "(встречаются в >=", heatmap_min_cell_types, "типах клеток, топ-",
      heatmap_top_n, "по padj)\n")

  mat <- build_heatmap_matrix(filtered_data)
  mat <- order_columns_by_group(mat)

  if (heatmap_remove_identical && nrow(mat) > 1 && ncol(mat) > 1) {
    row_vars <- apply(mat, 1, var)
    mat <- mat[row_vars > 0, , drop = FALSE]
  }

  if (nrow(mat) < 2) {
    cat("  Недостаточно генов для heatmap:", title_text, "\n")
    return(invisible(NULL))
  }

  save_heatmap_txt(mat, paste0(filename_base, ".txt"))

  annotation_col <- make_annotation(mat)

  n_colors <- 100
  heatmap_colors <- colorRampPalette(c(heatmap_color_low, heatmap_color_mid,
                                       heatmap_color_high))(n_colors)

  ann_colors <- list(Group = c(
    Glutamatergic = "#E41A1C",
    GABAergic     = "#377EB8",
    Glia          = "#4DAF4A",
    Microglia     = "#984EA3"
  ))

  do_cluster_cols <- ncol(mat) > 1
  do_cluster_rows <- nrow(mat) > 1

  h_adj <- max(heatmap_height, nrow(mat) * 0.15 + 2)

  fpath <- file.path(heatmap_outdir, paste0(filename_base, ".png"))
  png(fpath, width = heatmap_width, height = h_adj, units = "in", res = 300)

  pheatmap(mat,
           scale = heatmap_scale,
           color = heatmap_colors,
           clustering_distance_rows = heatmap_distance,
           clustering_distance_cols = heatmap_distance,
           clustering_method = heatmap_cluster_method,
           cluster_rows = do_cluster_rows,
           cluster_cols = do_cluster_cols,
           show_rownames = heatmap_show_rownames,
           show_colnames = TRUE,
           fontsize_row = heatmap_fontsize_row,
           fontsize_col = heatmap_fontsize_col,
           border_color = ifelse(heatmap_border, "grey60", NA),
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           main = title_text,
           display_numbers = FALSE)

  dev.off()
  cat("  Heatmap сохранена:", filename_base, ".png\n")
}

# 1. Общая heatmap
group_titles <- c(
  "Glutamatergic" = "Glutamatergic neurons",
  "GABAergic"     = "GABAergic neurons",
  "Glia"          = "Glia",
  "Microglia"     = "Microglia"
)

draw_heatmap(degs, "Heatmap_All_CellTypes", "K18-ACE2-KI + FS vs K18-ACE2-KI + FS + par")

for (grp in group_order) {
  
  sub <- degs %>% dplyr::filter(Group == grp)
  fname <- paste0("Heatmap_", grp)
  
  draw_heatmap(
    sub,
    fname,
    group_titles[grp]
  )
}

cat("\n=== Все тепловые карты построены ===\n\n")

# --- ЗАДАЧА 2: VOLCANO PLOTS ---

dir.create(volcano_outdir, showWarnings = FALSE, recursive = TRUE)

draw_volcano <- function(data, filename_base, title_text) {
  if (nrow(data) == 0) {
    cat("  Нет данных для volcano:", title_text, "\n")
    return(invisible(NULL))
  }

  df <- data %>%
    dplyr::mutate(
      padj = ifelse(padj < p_min, p_min, padj),
      neg_log10_padj = -log10(padj),
      Category = dplyr::case_when(
        log2FoldChange >= log2(fc_cutoff) & padj <= fdr_cutoff ~ "Upregulated",
        log2FoldChange <= -log2(fc_cutoff) & padj <= fdr_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )

  df$Category <- factor(df$Category,
                         levels = c("Upregulated", "Downregulated", "Not significant"))

  color_map <- c("Upregulated" = color_up,
                 "Downregulated" = color_down,
                 "Not significant" = color_ns)

  x_lo <- ifelse(is.null(xmin), min(df$log2FoldChange, na.rm = TRUE) * 1.1, xmin)
  x_hi <- ifelse(is.null(xmax), max(df$log2FoldChange, na.rm = TRUE) * 1.1, xmax)

  y_hi <- ifelse(is.null(ymax),
                 max(df$neg_log10_padj, na.rm = TRUE) * 1.05,
                 ymax)

  sig_df <- df %>% dplyr::filter(Category != "Not significant")
  label_df <- data.frame()
  if (n_label_genes > 0 && nrow(sig_df) > 0) {
    top_by_padj <- sig_df %>% dplyr::arrange(padj) %>% head(n_label_genes)
    top_by_fc <- sig_df %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% head(n_label_genes)
    label_df <- dplyr::bind_rows(top_by_padj, top_by_fc) %>% dplyr::distinct(gene, .keep_all = TRUE)
  }

  n_up <- sum(df$Category == "Upregulated")
  n_down <- sum(df$Category == "Downregulated")
  n_ns <- sum(df$Category == "Not significant")

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = Category)) +
    geom_point(size = point_size / 20, shape = point_shape, alpha = point_alpha) +
    scale_color_manual(
      values = color_map,
      labels = c(
        paste0("Up regulated (", n_up, ")"),
        paste0("Down regulated (", n_down, ")"),
        paste0("Not sig (", n_ns, ")")
      )
    ) +
    geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
               linetype = "dashed", color = "#000000", linewidth = 0.5) +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "#000000", linewidth = 0.5) +
    labs(x = comparison_label,
         y = expression(-log[10](FDR)),
         title = title_text,
         color = NULL) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(ymin, y_hi)) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_text(size = ticks_fontsize),
      axis.title = element_text(size = axis_label_fontsize),
      legend.text = element_text(size = legend_fontsize),
      legend.position = "top",
      panel.grid = if (show_grid) element_line(color = "grey90") else element_blank(),
      plot.title = element_text(hjust = 0.5, size = axis_label_fontsize)
    )

  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = label_fontsize / 3,
      color = label_color,
      max.overlaps = 30,
      segment.color = "grey50",
      segment.size = 0.3,
      box.padding = 0.5,
      point.padding = 0.3
    )
  }

  fpath <- file.path(volcano_outdir, paste0(filename_base, ".png"))
  ggsave(fpath, plot = p, width = figure_width, height = figure_height, dpi = 300)
  cat("  Volcano plot сохранён:", filename_base, ".png\n")
}

# 1. Все DEGs
draw_volcano(degs, "10 vs 13 Volcano plot", "K18-ACE2-KI + FS vs K18-ACE2-KI + FS + par")

# 2-5. По группам
volcano_groups <- list(
  list(group = "Glutamatergic", title = "Glutamatergic neurons",
       file = "Volcano_Glutamatergic"),
  list(group = "GABAergic",     title = "GABAergic neurons",
       file = "Volcano_GABAergic"),
  list(group = "Glia",          title = "Glia",
       file = "Volcano_Glia"),
  list(group = "Microglia",     title = "Microglia",
       file = "Volcano_Microglia")
)

for (vg in volcano_groups) {
  sub <- degs %>% dplyr::filter(Group == vg$group)
  
  draw_volcano(
    sub,
    vg$file,
    vg$title
  )
}

cat("\n=== Все volcano plots построены ===\n")
cat("Heatmaps:", heatmap_outdir, "\n")
cat("Volcano plots:", volcano_outdir, "\n")
