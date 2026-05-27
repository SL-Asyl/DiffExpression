# ==============================================================================
# AUTOMATED ORA BATCH ANALYSIS (clusterProfiler, KEGG, GO, MP)
# ==============================================================================

# ---------- Пакеты ------------------------------------------------------------
# Список пакетов из CRAN
cran_packages <- c("openxlsx", "enrichR", "dplyr", "tidyr", "ggplot2", "BiocManager")

# Список пакетов из Bioconductor
bioc_packages <- c("clusterProfiler", "AnnotationDbi", "org.Mm.eg.db", 
                   "ReactomePA", "rWikiPathways")

# 1. Установка пропущенных пакетов из CRAN
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  cat("Установка пропущенных пакетов CRAN:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, dependencies = TRUE)
}

# 2. Установка пропущенных пакетов из Bioconductor
library(BiocManager)
missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(missing_bioc) > 0) {
  cat("Установка пропущенных пакетов Bioconductor:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# 3. Загрузка всех библиотек в сессию R
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(openxlsx)
library(enrichR)
library(ReactomePA)
library(rWikiPathways)

# ---------- НАСТРОЙКИ ПОЛЬЗОВАТЕЛЯ (USER INPUT) -------------------------------
input_dir   <- "raw_data"     # Папка, где лежат исходные Excel
outdir      <- "ORA_results"  # Папка, куда сохранятся результаты

logFC_cutoff <- 1
FDR_cutoff   <- 0.05

dir.create(outdir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Вспомогательные функции
# ------------------------------------------------------------------------------
run_kegg_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for KEGG ORA: ", label); return(NULL) }
  enrichKEGG(gene = genes, organism = "mmu", keyType = "ncbi-geneid",
             universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             pAdjustMethod = "BH")
}

run_go_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for GO ORA: ", label); return(NULL) }
  enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL",
           universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
           pAdjustMethod = "BH", readable = TRUE)
}

mp_dbs <- c("MGI_Mammalian_Phenotype_Level_4_2021", "MGI_Mammalian_Phenotype_Level_3_2021")
run_mp_enrichr <- function(symbols, label) {
  if (length(symbols) < 5) { message("  Too few genes for MP ORA: ", label); return(NULL) }
  enrichr(symbols, mp_dbs)
}

format_ora_df <- function(res) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  df <- as.data.frame(res@result)
  
  if ("geneID" %in% colnames(df)) {
    df$Genes <- vapply(df$geneID, function(x) {
      ids <- unlist(strsplit(x, "/"))
      if (all(grepl("^[0-9]+$", ids))) {
        paste(mapIds(org.Mm.eg.db, keys = ids, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"), collapse = ", ")
      } else paste(ids, collapse = ", ")
    }, character(1))
    df <- select(df, -geneID)
  }
  
  for (col in intersect(c("GeneRatio", "BgRatio"), colnames(df))) {
    df[[col]] <- sprintf("%.1f%%", vapply(strsplit(df[[col]], "/"), \(x) as.numeric(x[1])/as.numeric(x[2])*100, numeric(1)))
  }
  for (col in intersect(c("FoldEnrichment","zScore","pvalue","p.adjust","qvalue"), colnames(df))) {
    df[[col]] <- signif(as.numeric(df[[col]]), 3)
  }
  
  df <- df %>%
    separate_rows(Genes, sep = ",\\s*") %>% 
    filter(Genes != "" & !is.na(Genes)) %>%
    select(Genes, everything()) %>%
    arrange(Genes, pvalue)
  
  df
}

format_mp_df <- function(mp_res, level = "MGI_Mammalian_Phenotype_Level_4_2021") {
  if (is.null(mp_res) || !level %in% names(mp_res)) return(NULL)
  df <- mp_res[[level]] %>%
    rename(Description = Term, pvalue = P.value, p.adjust = Adjusted.P.value, Genes = Genes) %>%
    mutate(pvalue = signif(pvalue, 3), p.adjust = signif(p.adjust, 3)) %>%
    separate_rows(Genes, sep = ";\\s*") %>% 
    filter(Genes != "" & !is.na(Genes)) %>%
    select(Genes, everything()) %>%
    arrange(Genes, pvalue)
  df
}

write_block <- function(wb, sheet, df, start_row, title = NULL) {
  if (is.null(df) || nrow(df) == 0) return(start_row)
  if (!is.null(title)) { writeData(wb, sheet, title, startRow = start_row, colNames = FALSE); start_row <- start_row + 1 }
  writeData(wb, sheet, df, startRow = start_row, colNames = TRUE)
  start_row + nrow(df) + 2
}

# ------------------------------------------------------------------------------
# ПОИСК СЫРЫХ ФАЙЛОВ В ПАПКЕ
# ------------------------------------------------------------------------------
raw_files <- list.files(path = input_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(raw_files) == 0) {
  stop("Критическая ошибка: не найдено файлов .xlsx в папке '", input_dir, "'. Проверьте путь!")
}

cat("Найдено файлов для анализа:", length(raw_files), "\n\n")

# ==============================================================================
# ОСНОВНОЙ ЦИКЛ АВТОМАТИЗАЦИИ
# ==============================================================================
for (file_path in raw_files) {
  
  # Достаем имя типа клеток из названия файла (например, "5_vs_6_Astro_TE_NN")
  cell_type <- gsub(".xlsx", "", basename(file_path))
  
  cat("============================================================\n")
  cat("ЗАПУСК АНАЛИЗА ДЛЯ:", cell_type, "\n")
  cat("============================================================\n")
  
  # Читаем файл данных
  tb_raw <- tryCatch({
    read.xlsx(file_path)
  }, error = function(e) {
    message("  Ошибка при чтении файла ", file_path, " Скипаем."); next
  })
  
  # 1. PREPARE DATA
  tb <- tb_raw %>%
    rename(gene = names, log2FoldChange = logfoldchanges, padj = pvals_adj) %>%
    select(gene, log2FoldChange, padj) %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj))
  
  # Проверка, есть ли вообще строчки в файле
  if (nrow(tb) == 0) { message("  Файл пустой после очистки NA. Скипаем."); next }
  
  suppressMessages({
    tb$ENTREZID <- mapIds(org.Mm.eg.db,
                          keys = tb$gene,
                          keytype = "SYMBOL",
                          column = "ENTREZID",
                          multiVals = "first")
  })
  
  tb <- tb %>% filter(!is.na(ENTREZID))
  if (nrow(tb) == 0) { message("  Не удалось сопоставить ENTREZ ID для генов. Скипаем."); next }
  
  # 2. DEFINE GENE SETS
  universe_genes <- unique(tb$ENTREZID)
  
  up_genes     <- tb %>% filter(padj < FDR_cutoff, log2FoldChange >= logFC_cutoff) %>% pull(ENTREZID) %>% unique()
  down_genes   <- tb %>% filter(padj < FDR_cutoff, log2FoldChange <= -logFC_cutoff) %>% pull(ENTREZID) %>% unique()
  
  up_symbols   <- tb %>% filter(padj < FDR_cutoff, log2FoldChange >= logFC_cutoff) %>% pull(gene) %>% unique()
  down_symbols <- tb %>% filter(padj < FDR_cutoff, log2FoldChange <= -logFC_cutoff) %>% pull(gene) %>% unique()
  
  cat("  Up genes:", length(up_genes), "| Down genes:", length(down_genes), "| Universe:", length(universe_genes), "\n")
  
  # 4. RUN ORA
  kegg_up   <- run_kegg_ora(up_genes, universe_genes, "UP")
  kegg_down <- run_kegg_ora(down_genes, universe_genes, "DOWN")
  go_up     <- run_go_ora(up_genes, universe_genes, "UP")
  go_down   <- run_go_ora(down_genes, universe_genes, "DOWN")
  mp_up     <- run_mp_enrichr(up_symbols, "UP")
  mp_down   <- run_mp_enrichr(down_symbols, "DOWN")
  
  # 5. SAVE RESULTS TO EXCEL С ИНДИВИДУАЛЬНЫМ ИМЕНЕМ
  wb <- createWorkbook()
  
  # KEGG
  if (!is.null(kegg_up)) { addWorksheet(wb, "KEGG_UP"); writeData(wb, "KEGG_UP", format_ora_df(kegg_up)) }
  if (!is.null(kegg_down)) { addWorksheet(wb, "KEGG_DOWN"); writeData(wb, "KEGG_DOWN", format_ora_df(kegg_down)) }
  
  # GO + MP
  for (type in c("UP","DOWN")) {
    go_res <- if(type=="UP") go_up else go_down
    mp_res <- if(type=="UP") mp_up else mp_down
    sheet_name <- paste0("GO_MP_", type)
    if (!is.null(go_res) || !is.null(mp_res)) {
      addWorksheet(wb, sheet_name); row <- 1
      if (!is.null(go_res)) row <- write_block(wb, sheet_name, format_ora_df(go_res), row, paste("GO enrichment (",type,")",sep=""))
      if (!is.null(mp_res)) row <- write_block(wb, sheet_name, format_mp_df(mp_res, "MGI_Mammalian_Phenotype_Level_4_2021"), row, paste("Mammalian Phenotype enrichment (",type,")",sep=""))
    }
  }
  
  # Генерируем уникальное имя файла для типа клеток
  output_excel_name <- file.path(outdir, paste0("ORA_KEGG_GO_", cell_type, ".xlsx"))
  
  tryCatch({
    saveWorkbook(wb, output_excel_name, overwrite = TRUE)
    cat("  Таблица результатов сохранена в:", output_excel_name, "\n")
  }, error = function(e) {
    message("  ОШИБКА сохранения Excel для ", cell_type, ". Возможно, файл открыт в другой программе.")
  })
  
  cat("  Анализ для", cell_type, "успешно завершен.\n\n")
}

cat("============================================================\n")
cat("ПАКЕТНЫЙ АНАЛИЗ ВСЕХ ТАБЛИЦ ПОЛНОСТЬЮ ЗАВЕРШЕН.\n")
cat("Все файлы Excel сохранены в папку:", outdir, "\n")
cat("============================================================\n")
