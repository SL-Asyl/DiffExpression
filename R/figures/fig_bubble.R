# =============================================================================
#  fig_bubble.R — команда `fig_bubble`: bubble-plots GO/KEGG из ГОТОВЫХ
#  мастер-таблиц (объединение Filtration.R ч.2 и BubblePlots.R). Таблицы не
#  пересчитывает.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, evs_dir])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "dplyr", "ggplot2"))
  suppressMessages({ library(openxlsx); library(dplyr); library(ggplot2) })
  
  folder <- evs_dir(tissue, cc)
  plots_dir <- file.path(folder, "EnrichmentPlots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  go_ontologies <- c("BP", "MF", "CC")
  
  wrap_labels <- function(x, width = BUBBLE_WRAP)
    vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
  parse_gene_ratio <- function(gr_col) {
    if (all(is.na(gr_col))) return(rep(NA_real_, length(gr_col)))
    result <- rep(NA_real_, length(gr_col))
    is_pct <- grepl("%$", gr_col) & !is.na(gr_col)
    if (any(is_pct)) result[is_pct] <- as.numeric(sub("%$", "", gr_col[is_pct])) / 100
    is_frac <- grepl("/", gr_col) & !is_pct & !is.na(gr_col)
    if (any(is_frac)) { parts <- strsplit(gr_col[is_frac], "/")
    result[is_frac] <- sapply(parts, function(x) as.numeric(x[1]) / as.numeric(x[2])) }
    is_num <- !is_pct & !is_frac & !is.na(gr_col)
    if (any(is_num)) result[is_num] <- suppressWarnings(as.numeric(gr_col[is_num]))
    result
  }
  draw_bubble_plot <- function(data, filename_base, title_text) {
    if (is.null(data) || nrow(data) == 0) { cat("  нет данных:", title_text, "\n"); return(invisible(NULL)) }
    df <- data
    if ("GeneRatio_num" %in% colnames(df) && any(!is.na(df$GeneRatio_num))) df$GR <- as.numeric(df$GeneRatio_num)
    else if ("GeneRatio" %in% colnames(df)) df$GR <- parse_gene_ratio(df$GeneRatio)
    else { cat("  нет GeneRatio:", title_text, "\n"); return(invisible(NULL)) }
    if (all(is.na(df$GR)) && "Count" %in% colnames(df) && "Genes" %in% colnames(df)) {
      df$Count <- as.numeric(df$Count)
      total_genes <- length(unique(unlist(strsplit(paste(df$Genes, collapse = ", "), ",\\s*"))))
      df$GR <- df$Count / max(total_genes, 1)
    }
    if (!"Count" %in% colnames(df))
      df$Count <- if ("Genes" %in% colnames(df)) sapply(strsplit(as.character(df$Genes), ",\\s*"), length) else 1
    df$Count <- as.numeric(df$Count)
    padj_col <- intersect(c("p.adjust", "Adjusted.P.value", "padj"), colnames(df))[1]
    if (is.na(padj_col)) { cat("  нет p.adjust:", title_text, "\n"); return(invisible(NULL)) }
    df$Padj <- as.numeric(df[[padj_col]])
    df <- df %>% filter(!is.na(GR), !is.na(Padj), !is.na(Description), GR > 0) %>% arrange(Padj) %>% head(BUBBLE_TOP_N)
    if (nrow(df) == 0) { cat("  пусто после фильтра:", title_text, "\n"); return(invisible(NULL)) }
    if (any(duplicated(df$Description)))
      df <- df %>% group_by(Description) %>% slice_min(order_by = Padj, n = 1, with_ties = FALSE) %>% ungroup()
    df$Description <- factor(df$Description, levels = rev(df$Description))
    p <- ggplot(df, aes(x = GR, y = Description, size = Count, color = Padj)) +
      geom_point() +
      scale_size_continuous(range = BUBBLE_PT_RANGE, name = "Count") +
      scale_color_gradient(low = BUBBLE_COL_LOW, high = BUBBLE_COL_HIGH, name = "Adjusted\nP value", trans = "log10") +
      scale_y_discrete(labels = wrap_labels) +
      labs(x = "GeneRatio", y = NULL, title = wrap_labels(title_text, 48)) +
      theme_bw(base_size = BUBBLE_FONT) +
      theme(aspect.ratio = 1,
            axis.text.y  = element_text(size = BUBBLE_FONT - 1, lineheight = 0.9),
            axis.text.x  = element_text(size = BUBBLE_FONT - 1),
            axis.title.x = element_text(size = BUBBLE_FONT + 1),
            plot.title   = element_text(hjust = 0.5, size = BUBBLE_FONT + 3, face = "bold"),
            legend.title = element_text(size = BUBBLE_FONT),
            legend.text  = element_text(size = BUBBLE_FONT - 2), legend.position = "right")
    ggsave(file.path(plots_dir, paste0(filename_base, ".png")), plot = p, width = BUBBLE_SIZE, height = BUBBLE_SIZE, dpi = 300)
    cat("  сохранён:", filename_base, ".png\n")
  }
  find_one <- function(folder, pattern) {
    f <- list.files(folder, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(f) == 0) NA_character_ else f[1]
  }
  
  cat("[fig_bubble]", tissue, TAG[[cc]], "\n")
  go_all_f  <- find_one(folder, "MASTER_GO_ALL_CELL_TYPES\\.xlsx$")
  go_filt_f <- find_one(folder, "MASTER_GO_FILTERED\\.xlsx$")
  kegg_f    <- find_one(folder, "MASTER_KEGG_ALL_CELL_TYPES\\.xlsx$")
  
  if (!is.na(go_all_f)) {
    go_all <- read_xlsx_safe(go_all_f)
    for (ont in go_ontologies)
      draw_bubble_plot(go_all[go_all$ONTOLOGY == ont, , drop = FALSE],
                       paste0("Bubble_GO_", ont, "_unfiltered"), paste("GO", ont, "(unfiltered)"))
    if ("MP" %in% go_all$ONTOLOGY)
      draw_bubble_plot(go_all[go_all$ONTOLOGY == "MP", , drop = FALSE], "Bubble_MP_unfiltered", "Mammalian Phenotype (unfiltered)")
  } else cat("  нет MASTER_GO_ALL — пропуск unfiltered GO\n")
  if (!is.na(kegg_f)) draw_bubble_plot(read_xlsx_safe(kegg_f), "Bubble_KEGG_unfiltered", "KEGG (unfiltered)")
  if (!is.na(go_filt_f)) {
    go_filt <- read_xlsx_safe(go_filt_f)
    for (ont in go_ontologies)
      draw_bubble_plot(go_filt[go_filt$ONTOLOGY == ont, , drop = FALSE],
                       paste0("Bubble_GO_", ont, "_filtered"), paste("GO", ont, "(filtered)"))
  } else cat("  нет MASTER_GO_FILTERED — пропуск filtered GO\n")
}
