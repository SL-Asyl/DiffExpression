# --- Пути ---
go_master_file   <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13_MASTER_GO_ALL_CELL_TYPES.xlsx"
kegg_master_file <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/10 vs 13_MASTER_KEGG_ALL_CELL_TYPES.xlsx"
output_dir       <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS"
plots_dir        <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS/EnrichmentPlots"

# --- Параметры семантической фильтрации ---
# rrvgo: работает офлайн через GO.db, вычисляет попарное семантическое сходство
# терминов через Information Content и объединяет избыточные термины в кластеры,
# оставляя наиболее значимый представитель.
# Применяется только к GO (BP, MF, CC) — для MP и KEGG нет иерархической
# онтологии, поэтому семантическая фильтрация к ним неприменима.
similarity_threshold <- 0.7
similarity_method    <- "Rel"  # "Lin", "Rel", "Jiang", "Resnik", "Wang"
organism_db          <- "org.Mm.eg.db"

# --- Параметры bubble plot ---
bubble_width  <- 10
bubble_height <- 8
top_n_terms   <- 20
bubble_point_size_range <- c(3, 10)
color_low  <- "#FF0000"
color_high <- "#0000FF"

# --- УСТАНОВКА И ЗАГРУЗКА ПАКЕТОВ ---

cran_pkgs <- c("openxlsx", "dplyr", "ggplot2", "stringr", "tidyr")
bioc_pkgs <- c("rrvgo", "org.Mm.eg.db", "GO.db")

missing_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

library(openxlsx)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(rrvgo)

# --- ЧТЕНИЕ ДАННЫХ ---

go_master <- read.xlsx(go_master_file)
kegg_master <- read.xlsx(kegg_master_file)

cat("GO/MP мастер-таблица:", nrow(go_master), "строк\n")
cat("KEGG мастер-таблица:", nrow(kegg_master), "строк\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# --- СЕМАНТИЧЕСКАЯ ФИЛЬТРАЦИЯ (GO BP, MF, CC)

go_ontologies <- c("BP", "MF", "CC")

split_go_by_ontology <- function(go_data) {
  result <- list()
  for (ont in unique(go_data$ONTOLOGY)) {
    result[[ont]] <- go_data %>% filter(ONTOLOGY == ont)
  }
  result
}

filter_go_semantic <- function(go_subset, ontology_type) {
  if (nrow(go_subset) == 0) return(go_subset)

  go_ids <- go_subset$ID
  valid_go <- grepl("^GO:", go_ids)

  if (sum(valid_go) < 2) {
    cat(" (менее 2 GO-терминов, фильтрация невозможна)")
    return(go_subset)
  }

  go_terms <- go_ids[valid_go]
  pvals <- go_subset$p.adjust[valid_go]
  names(pvals) <- go_terms

  sim_matrix <- tryCatch({
    calculateSimMatrix(go_terms,
                       orgdb = organism_db,
                       ont = ontology_type,
                       method = similarity_method)
  }, error = function(e) {
    cat("\n  ОШИБКА расчёта сходства для", ontology_type, ":", e$message)
    return(NULL)
  })

  if (is.null(sim_matrix) || nrow(sim_matrix) < 2) {
    cat(" (матрица сходства не построена)")
    return(go_subset)
  }

  pvals_matched <- pvals[rownames(sim_matrix)]
  pvals_matched <- pvals_matched[!is.na(pvals_matched)]

  reduced <- reduceSimMatrix(sim_matrix, pvals_matched,
                             threshold = similarity_threshold,
                             orgdb = organism_db)

  best_per_cluster <- reduced %>%
    group_by(parentTerm) %>%
    slice_min(order_by = score, n = 1, with_ties = FALSE) %>%
    pull(go)

  keep_ids <- unique(c(best_per_cluster, go_ids[!valid_go]))
  filtered <- go_subset %>% filter(ID %in% keep_ids)
  filtered
}

go_split <- split_go_by_ontology(go_master)

cat("\nОнтологии в мастер-таблице:", paste(names(go_split), collapse = ", "), "\n\n")

go_filtered_list <- list()

for (ont_name in go_ontologies) {
  if (!ont_name %in% names(go_split)) {
    cat("Онтология", ont_name, "отсутствует в данных\n")
    next
  }
  cat("Фильтрация", ont_name, "- исходно:", nrow(go_split[[ont_name]]), "терминов...")
  go_filtered_list[[ont_name]] <- filter_go_semantic(go_split[[ont_name]], ont_name)
  cat(" после:", nrow(go_filtered_list[[ont_name]]), "\n")
}

# --- Сборка фильтрованной мастер-таблицы ---
filtered_go_combined <- bind_rows(go_filtered_list)

# MP не фильтруем, но включаем в итоговую таблицу как есть
if ("MP" %in% names(go_split)) {
  cat("MP: семантическая фильтрация неприменима (нет иерархической онтологии),",
      nrow(go_split[["MP"]]), "терминов оставлены без изменений\n")
  filtered_go_combined <- bind_rows(filtered_go_combined, go_split[["MP"]])
}

cat("\nИтого фильтрованная GO мастер-таблица:", nrow(filtered_go_combined), "строк",
    "(было", nrow(go_master), ")\n")

# Сохранение единой фильтрованной мастер-таблицы
wb <- createWorkbook()
addWorksheet(wb, "GO_MP_Filtered")
writeData(wb, "GO_MP_Filtered", filtered_go_combined)
setColWidths(wb, "GO_MP_Filtered", 1:ncol(filtered_go_combined), widths = "auto")
addFilter(wb, "GO_MP_Filtered", row = 1, cols = 1:ncol(filtered_go_combined))
freezePane(wb, "GO_MP_Filtered", firstRow = TRUE)
fpath <- file.path(output_dir, "MASTER_GO_FILTERED.xlsx")
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("Фильтрованная мастер-таблица GO сохранена:", fpath, "\n")

# KEGG не фильтруем — нет онтологической иерархии
cat("KEGG: семантическая фильтрация неприменима (нет иерархической онтологии)\n")

cat("\n=== Семантическая фильтрация завершена ===\n\n")

# --- ВИЗУАЛИЗАЦИЯ ENRICHMENT BUBBLE PLOTS ---

parse_gene_ratio <- function(gr_col) {
  if (all(is.na(gr_col))) return(rep(NA_real_, length(gr_col)))

  result <- rep(NA_real_, length(gr_col))

  is_pct <- grepl("%$", gr_col) & !is.na(gr_col)
  if (any(is_pct)) {
    result[is_pct] <- as.numeric(sub("%$", "", gr_col[is_pct])) / 100
  }

  is_frac <- grepl("/", gr_col) & !is_pct & !is.na(gr_col)
  if (any(is_frac)) {
    parts <- strsplit(gr_col[is_frac], "/")
    result[is_frac] <- sapply(parts, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  }

  is_num <- !is_pct & !is_frac & !is.na(gr_col)
  if (any(is_num)) {
    result[is_num] <- suppressWarnings(as.numeric(gr_col[is_num]))
  }

  result
}

draw_bubble_plot <- function(data, filename_base, title_text, require_success = FALSE) {
  if (is.null(data) || nrow(data) == 0) {
    msg <- paste("Нет данных для bubble plot:", title_text)
    if (require_success) warning(msg) else cat("  ", msg, "\n")
    return(invisible(NULL))
  }

  df <- data

  # Вычислить GeneRatio
  if ("GeneRatio_num" %in% colnames(df) && any(!is.na(df$GeneRatio_num))) {
    df$GR <- as.numeric(df$GeneRatio_num)
  } else if ("GeneRatio" %in% colnames(df)) {
    df$GR <- parse_gene_ratio(df$GeneRatio)
  } else {
    msg <- paste("Нет столбца GeneRatio для:", title_text)
    if (require_success) warning(msg) else cat("  ", msg, "\n")
    return(invisible(NULL))
  }

  # Для MP: GeneRatio может быть NA — вычислить из Count / общего числа генов
  if (all(is.na(df$GR)) && "Count" %in% colnames(df) && "Genes" %in% colnames(df)) {
    df$Count <- as.numeric(df$Count)
    total_genes <- length(unique(unlist(strsplit(paste(df$Genes, collapse = ", "), ",\\s*"))))
    df$GR <- df$Count / max(total_genes, 1)
  }

  if (!"Count" %in% colnames(df)) {
    if ("Genes" %in% colnames(df)) {
      df$Count <- sapply(strsplit(as.character(df$Genes), ",\\s*"), length)
    } else {
      df$Count <- 1
    }
  }
  df$Count <- as.numeric(df$Count)

  padj_col <- intersect(c("p.adjust", "Adjusted.P.value", "padj"), colnames(df))[1]
  if (is.na(padj_col)) {
    msg <- paste("Нет столбца p.adjust для:", title_text)
    if (require_success) warning(msg) else cat("  ", msg, "\n")
    return(invisible(NULL))
  }
  df$Padj <- as.numeric(df[[padj_col]])

  df <- df %>%
    filter(!is.na(GR), !is.na(Padj), !is.na(Description), GR > 0) %>%
    arrange(Padj) %>%
    head(top_n_terms)

  if (nrow(df) == 0) {
    msg <- paste("Пустой набор после фильтрации:", title_text,
                 "(возможно все GeneRatio = NA или 0)")
    if (require_success) warning(msg) else cat("  ", msg, "\n")
    return(invisible(NULL))
  }

  # Убрать дубликаты описаний — оставить строку с наименьшим p.adjust
  if (any(duplicated(df$Description))) {
    df <- df %>%
      group_by(Description) %>%
      slice_min(order_by = Padj, n = 1, with_ties = FALSE) %>%
      ungroup()
  }

  df$Description <- factor(df$Description, levels = rev(df$Description))

  p <- ggplot(df, aes(x = GR, y = Description, size = Count, color = Padj)) +
    geom_point() +
    scale_size_continuous(range = bubble_point_size_range, name = "Count") +
    scale_color_gradient(low = color_low, high = color_high,
                         name = "Adjusted\nP value",
                         trans = "log10") +
    labs(x = "GeneRatio", y = NULL, title = title_text) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )

  h_adj <- max(bubble_height, nrow(df) * 0.35 + 2)

  fpath <- file.path(plots_dir, paste0(filename_base, ".png"))
  ggsave(fpath, plot = p, width = bubble_width, height = h_adj, dpi = 300)
  cat("  Bubble plot сохранён:", filename_base, ".png\n")
  return(TRUE)
}



# --- Без фильтрации ---
cat("--- Bubble plots БЕЗ фильтрации ---\n")

for (ont_name in go_ontologies) {
  if (!ont_name %in% names(go_split)) {
    cat("  ВНИМАНИЕ: онтология", ont_name, "отсутствует в данных\n")
    next
  }
  draw_bubble_plot(go_split[[ont_name]],
                   paste0("Bubble_GO_", ont_name, "_unfiltered"),
                   paste("GO", ont_name, "(unfiltered)"))
}

if ("MP" %in% names(go_split)) {
  draw_bubble_plot(go_split[["MP"]],
                   "Bubble_MP_unfiltered",
                   "Mammalian Phenotype (unfiltered)")
} else {
  cat("  ВНИМАНИЕ: MP отсутствует в мастер-таблице GO\n")
}

draw_bubble_plot(kegg_master, "Bubble_KEGG_unfiltered", "KEGG (unfiltered)")

# --- После фильтрации (только GO BP, MF, CC) ---
cat("\n--- Bubble plots ПОСЛЕ фильтрации ---\n")

for (ont_name in go_ontologies) {
  if (!ont_name %in% names(go_filtered_list)) {
    cat("  ОШИБКА: фильтрованные данные для", ont_name, "отсутствуют!\n")
    next
  }
  result <- draw_bubble_plot(
    go_filtered_list[[ont_name]],
    paste0("Bubble_GO_", ont_name, "_filtered"),
    paste("GO", ont_name, "(filtered)"),
    require_success = TRUE
  )
  if (is.null(result)) {
    cat("  ОШИБКА: bubble plot для GO", ont_name,
        "filtered НЕ создан. Строк данных:",
        nrow(go_filtered_list[[ont_name]]),
        ", из них с GeneRatio:",
        sum(!is.na(go_filtered_list[[ont_name]]$GeneRatio_num)), "\n")
  }
}

cat("\n=== Все enrichment bubble plots построены ===\n")
cat("Фильтрованная мастер-таблица:", file.path(output_dir, "MASTER_GO_FILTERED.xlsx"), "\n")
cat("Графики:", plots_dir, "\n")
