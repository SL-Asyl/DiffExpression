# =============================================================================
#  fig_volcano.R — команда `fig_volcano`: volcano по датасету (из Pictures.R ч.1)
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, evs_dir])
box::use(R/core/celltypes[load_degs_classified, GROUP_ORDER])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "ggplot2", "ggrepel", "dplyr"))
  suppressMessages({ library(openxlsx); library(ggplot2); library(ggrepel); library(dplyr) })
  
  degs <- load_degs_classified(tissue, cc)
  volcano_outdir <- file.path(evs_dir(tissue, cc), paste0(TAG[[cc]], " VolcanoPlots"))
  dir.create(volcano_outdir, showWarnings = FALSE, recursive = TRUE)
  comparison_label <- paste(tissue, TAG[[cc]])
  cat("[fig_volcano]", tissue, TAG[[cc]], "- DEG:", nrow(degs), "\n")
  
  draw_volcano <- function(data, filename_base, title_text) {
    if (nrow(data) == 0) { cat("  нет данных:", title_text, "\n"); return(invisible(NULL)) }
    df <- data %>% dplyr::mutate(
      padj = ifelse(padj < VOLC_PMIN, VOLC_PMIN, padj),
      neg_log10_padj = -log10(padj),
      Category = dplyr::case_when(
        log2FoldChange >= log2(VOLC_FC) & padj <= VOLC_FDR ~ "Upregulated",
        log2FoldChange <= -log2(VOLC_FC) & padj <= VOLC_FDR ~ "Downregulated",
        TRUE ~ "Not significant"))
    df$Category <- factor(df$Category, levels = c("Upregulated", "Downregulated", "Not significant"))
    color_map <- c("Upregulated" = VOLC_COL_UP, "Downregulated" = VOLC_COL_DOWN, "Not significant" = VOLC_COL_NS)
    x_lo <- min(df$log2FoldChange, na.rm = TRUE) * 1.1; x_hi <- max(df$log2FoldChange, na.rm = TRUE) * 1.1
    y_hi <- max(df$neg_log10_padj, na.rm = TRUE) * 1.05
    sig_df <- df %>% dplyr::filter(Category != "Not significant"); label_df <- data.frame()
    if (VOLC_N_LABEL > 0 && nrow(sig_df) > 0) {
      top_by_padj <- sig_df %>% dplyr::arrange(padj) %>% head(VOLC_N_LABEL)
      top_by_fc <- sig_df %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% head(VOLC_N_LABEL)
      label_df <- dplyr::bind_rows(top_by_padj, top_by_fc) %>% dplyr::distinct(gene, .keep_all = TRUE)
    }
    n_up <- sum(df$Category == "Upregulated"); n_down <- sum(df$Category == "Downregulated")
    n_ns <- sum(df$Category == "Not significant")
    p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = Category)) +
      geom_point(size = VOLC_PT_SIZE / 20, shape = 16, alpha = 0.7) +
      scale_color_manual(values = color_map,
                         labels = c(paste0("Up regulated (", n_up, ")"), paste0("Down regulated (", n_down, ")"),
                                    paste0("Not sig (", n_ns, ")"))) +
      geom_vline(xintercept = c(-log2(VOLC_FC), log2(VOLC_FC)), linetype = "dashed", color = "#000000", linewidth = 0.5) +
      geom_hline(yintercept = -log10(VOLC_FDR), linetype = "dashed", color = "#000000", linewidth = 0.5) +
      labs(x = comparison_label, y = expression(-log[10](FDR)), title = title_text, color = NULL) +
      coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(-0.1, y_hi)) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(size = 18), axis.title = element_text(size = 24),
            legend.text = element_text(size = 14), legend.position = "top",
            panel.grid = element_line(color = "grey90"), plot.title = element_text(hjust = 0.5, size = 24))
    if (nrow(label_df) > 0)
      p <- p + geom_text_repel(data = label_df, aes(label = gene), size = 14 / 3, color = "#FF0000",
                               max.overlaps = 30, segment.color = "grey50", segment.size = 0.3, box.padding = 0.5, point.padding = 0.3)
    ggsave(file.path(volcano_outdir, paste0(filename_base, ".png")), plot = p, width = VOLC_W, height = VOLC_H, dpi = 300)
    cat("  сохранён:", filename_base, ".png\n")
  }
  
  draw_volcano(degs, paste0(TAG[[cc]], " Volcano plot"), paste(tissue, TAG[[cc]]))
  volcano_groups <- list(
    list(group = "Glutamatergic", title = "Glutamatergic neurons", file = "Volcano_Glutamatergic"),
    list(group = "GABAergic",     title = "GABAergic neurons",     file = "Volcano_GABAergic"),
    list(group = "Glia",          title = "Glia",                  file = "Volcano_Glia"),
    list(group = "Microglia",     title = "Microglia",             file = "Volcano_Microglia"))
  for (vg in volcano_groups)
    draw_volcano(degs %>% dplyr::filter(Group == vg$group), vg$file, vg$title)
  cat("  Volcano ->", volcano_outdir, "\n")
}
