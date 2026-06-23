cran_packages <- c("openxlsx", "enrichR", "dplyr", "tidyr", "stringr", "BiocManager",
                   "future", "future.apply", "progressr")
bioc_packages <- c("clusterProfiler", "AnnotationDbi", "org.Mm.eg.db",
                   "ReactomePA", "rWikiPathways")

missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  cat("Установка пропущенных пакетов CRAN:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, dependencies = TRUE)
}

library(BiocManager)
missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(missing_bioc) > 0) {
  cat("Установка пропущенных пакетов Bioconductor:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(enrichR)
library(ReactomePA)
library(rWikiPathways)
library(future)
library(future.apply)
library(progressr)

handlers(handler_txtprogressbar())

input_dir  <- "исходный путь"
output_dir <- "путь сохранения" 

logFC_cutoff <- 1
FDR_cutoff   <- 0.05

mp_db <- "MGI_Mammalian_Phenotype_Level_4_2021"


cell_type_clean_pattern <- "5_vs_6_" # если нужно убрать приставку для чистых названий

n_workers          <- 2          # количество параллельной обработки: не рекомендуется больше 3
retry_max_tries    <- 3          # сколько раз повторить сетевой вызов при ошибке
retry_wait_sec     <- 5          # базовая пауза перед повтором (растёт с каждой попыткой)
jitter_range_sec   <- c(0.5, 2)  # случайная пауза перед каждым сетевым вызовом

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

.symbol_cache <- new.env(parent = emptyenv())
.entrez_cache <- new.env(parent = emptyenv())

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  uncached <- symbols[!vapply(symbols, exists, logical(1), envir = .entrez_cache, inherits = FALSE)]
  
  if (length(uncached) > 0) {
    ids <- suppressMessages(mapIds(org.Mm.eg.db,
                                   keys = uncached,
                                   keytype = "SYMBOL",
                                   column = "ENTREZID",
                                   multiVals = "first"))
    for (sym in uncached) {
      val <- ids[[sym]]
      assign(sym, if (is.null(val) || is.na(val)) NA_character_ else val, envir = .entrez_cache)
    }
  }
  
  vapply(symbols, function(sym) get(sym, envir = .entrez_cache), character(1), USE.NAMES = TRUE)
}

get_symbols_for_ids <- function(slash_separated_ids) {
  all_ids <- unique(unlist(strsplit(slash_separated_ids, "/", fixed = TRUE)))
  is_entrez <- grepl("^[0-9]+$", all_ids)
  entrez_ids <- all_ids[is_entrez]
  
  uncached <- entrez_ids[!vapply(entrez_ids, exists, logical(1), envir = .symbol_cache, inherits = FALSE)]
  
  if (length(uncached) > 0) {
    syms <- suppressMessages(mapIds(org.Mm.eg.db,
                                    keys = uncached,
                                    keytype = "ENTREZID",
                                    column = "SYMBOL",
                                    multiVals = "first"))
    for (id in uncached) {
      val <- syms[[id]]
      assign(id, if (is.null(val) || is.na(val)) id else val, envir = .symbol_cache)
    }
  }
  
  non_entrez <- all_ids[!is_entrez]
  for (id in non_entrez) {
    if (!exists(id, envir = .symbol_cache, inherits = FALSE)) assign(id, id, envir = .symbol_cache)
  }
  
  vapply(slash_separated_ids, function(x) {
    ids <- strsplit(x, "/", fixed = TRUE)[[1]]
    paste(vapply(ids, function(id) get(id, envir = .symbol_cache), character(1)), collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

jitter_pause <- function() Sys.sleep(runif(1, jitter_range_sec[1], jitter_range_sec[2]))

with_retry <- function(expr_fun, label) {
  for (attempt in seq_len(retry_max_tries)) {
    res <- tryCatch(expr_fun(), error = function(e) {
      message("  [", label, "] попытка ", attempt, "/", retry_max_tries,
              " не удалась: ", conditionMessage(e))
      NULL
    })
    if (!is.null(res)) return(res)
    if (attempt < retry_max_tries) Sys.sleep(retry_wait_sec * attempt)
  }
  message("  [", label, "] все попытки исчерпаны, пропускаем.")
  NULL
}

run_kegg_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for KEGG ORA: ", label); return(NULL) }
  jitter_pause()
  with_retry(function() {
    enrichKEGG(gene = genes, organism = "mmu", keyType = "ncbi-geneid",
               universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               pAdjustMethod = "BH")
  }, label = paste("KEGG", label))
}

run_go_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for GO ORA: ", label); return(NULL) }
  jitter_pause()
  with_retry(function() {
    enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL",
             universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             pAdjustMethod = "BH", readable = TRUE)
  }, label = paste("GO", label))
}

run_mp_enrichr <- function(symbols, label) {
  if (length(symbols) < 5) { message("  Too few genes for MP ORA: ", label); return(NULL) }
  jitter_pause()
  with_retry(function() enrichr(symbols, mp_db), label = paste("Enrichr", label))
}

format_kegg_master <- function(res, cell_type, direction) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  df <- as.data.frame(res@result)
  
  df$Genes <- get_symbols_for_ids(df$geneID)
  
  gr <- do.call(rbind, strsplit(df$GeneRatio, "/", fixed = TRUE))
  df$GeneRatio_num <- as.numeric(gr[, 1]) / as.numeric(gr[, 2])
  df$GeneRatio <- sprintf("%.1f%%", df$GeneRatio_num * 100)
  
  bg <- do.call(rbind, strsplit(df$BgRatio, "/", fixed = TRUE))
  df$BgRatio <- sprintf("%.1f%%", as.numeric(bg[, 1]) / as.numeric(bg[, 2]) * 100)
  
  for (col in c("pvalue", "p.adjust", "qvalue")) df[[col]] <- signif(as.numeric(df[[col]]), 3)
  
  df$Cell_Type <- cell_type
  df$Direction <- direction
  
  df %>% select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
                Count, Genes, Cell_Type, Direction, GeneRatio_num)
}

format_go_master <- function(res, cell_type, direction) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  df <- as.data.frame(res@result)
  
  df$Genes <- get_symbols_for_ids(df$geneID)
  
  gr <- do.call(rbind, strsplit(df$GeneRatio, "/", fixed = TRUE))
  df$GeneRatio_num <- as.numeric(gr[, 1]) / as.numeric(gr[, 2])
  df$GeneRatio <- sprintf("%.1f%%", df$GeneRatio_num * 100)
  
  bg <- do.call(rbind, strsplit(df$BgRatio, "/", fixed = TRUE))
  df$BgRatio <- sprintf("%.1f%%", as.numeric(bg[, 1]) / as.numeric(bg[, 2]) * 100)
  
  for (col in intersect(c("pvalue", "p.adjust", "qvalue"), colnames(df))) {
    df[[col]] <- signif(as.numeric(df[[col]]), 3)
  }
  
  df$Cell_Type <- cell_type
  df$Direction <- direction
  
  df %>% select(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
                Count, Genes, Cell_Type, Direction, GeneRatio_num)
}

format_mp_master <- function(mp_res, cell_type, direction) {
  if (is.null(mp_res) || !mp_db %in% names(mp_res)) return(NULL)
  df <- mp_res[[mp_db]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  id   <- str_extract(df$Term, "MP:\\d+")
  desc <- str_trim(gsub("\\(?MP:\\d+\\)?", "", df$Term))
  genes <- gsub(";\\s*", ", ", df$Genes)
  cnt   <- vapply(strsplit(df$Genes, ";\\s*"), length, integer(1))
  
  out <- data.frame(
    ID            = id,
    Description   = desc,
    ONTOLOGY      = "MP",
    GeneRatio     = NA_character_,
    BgRatio       = NA_character_,
    pvalue        = signif(as.numeric(df$P.value), 3),
    p.adjust      = signif(as.numeric(df$Adjusted.P.value), 3),
    qvalue        = NA_real_,
    Count         = cnt,
    Genes         = genes,
    Cell_Type     = cell_type,
    Direction     = direction,
    GeneRatio_num = NA_real_,
    stringsAsFactors = FALSE
  )
  
  out %>% filter(!is.na(ID))
}

process_file <- function(file_path, p = NULL) {
  suppressMessages({
    library(clusterProfiler); library(AnnotationDbi); library(org.Mm.eg.db)
    library(dplyr); library(tidyr); library(stringr); library(openxlsx); library(enrichR)
  })
  
  cell_type <- gsub("\\.xlsx$", "", basename(file_path))
  
  if (!is.null(p)) on.exit(p(sprintf("%s", cell_type)), add = TRUE)
  
  cat("ЗАПУСК АНАЛИЗА ДЛЯ:", cell_type, "\n")
  
  tb_raw <- tryCatch(read.xlsx(file_path), error = function(e) NULL)
  if (is.null(tb_raw)) { message("  Ошибка при чтении файла ", file_path, ". Пропускаем."); return(NULL) }
  
  tb <- tb_raw %>%
    rename(gene = names, log2FoldChange = logfoldchanges, padj = pvals_adj) %>%
    select(gene, log2FoldChange, padj) %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj))
  
  if (nrow(tb) == 0) { message("  Файл пустой после очистки NA. Пропускаем."); return(NULL) }
  
  tb$ENTREZID <- map_symbols_to_entrez(tb$gene)[tb$gene]
  
  tb <- tb %>% filter(!is.na(ENTREZID))
  if (nrow(tb) == 0) { message("  Не удалось сопоставить ENTREZ ID для генов. Пропускаем."); return(NULL) }
  
  universe_genes <- unique(tb$ENTREZID)
  
  up   <- tb %>% filter(padj < FDR_cutoff, log2FoldChange >=  logFC_cutoff)
  down <- tb %>% filter(padj < FDR_cutoff, log2FoldChange <= -logFC_cutoff)
  
  up_genes     <- unique(up$ENTREZID);  down_genes   <- unique(down$ENTREZID)
  up_symbols   <- unique(up$gene);      down_symbols <- unique(down$gene)
  
  cat("  [", cell_type, "] Up:", length(up_genes), "| Down:", length(down_genes),
      "| Universe:", length(universe_genes), "\n")
  
  kegg_up   <- run_kegg_ora(up_genes, universe_genes, paste(cell_type, "UP"))
  kegg_down <- run_kegg_ora(down_genes, universe_genes, paste(cell_type, "DOWN"))
  go_up     <- run_go_ora(up_genes, universe_genes, paste(cell_type, "UP"))
  go_down   <- run_go_ora(down_genes, universe_genes, paste(cell_type, "DOWN"))
  mp_up     <- run_mp_enrichr(up_symbols, paste(cell_type, "UP"))
  mp_down   <- run_mp_enrichr(down_symbols, paste(cell_type, "DOWN"))
  
  kegg_rows <- list(format_kegg_master(kegg_up, cell_type, "UP"),
                    format_kegg_master(kegg_down, cell_type, "DOWN"))
  go_rows   <- list(format_go_master(go_up, cell_type, "UP"),
                    format_go_master(go_down, cell_type, "DOWN"),
                    format_mp_master(mp_up, cell_type, "UP"),
                    format_mp_master(mp_down, cell_type, "DOWN"))
  
  kegg_rows <- kegg_rows[!vapply(kegg_rows, is.null, logical(1))]
  go_rows   <- go_rows[!vapply(go_rows, is.null, logical(1))]
  
  cat("  Анализ для", cell_type, "успешно завершен.\n\n")
  list(kegg = kegg_rows, go = go_rows)
}

raw_files <- list.files(path = input_dir, pattern = "\\.xlsx$", full.names = TRUE)
if (length(raw_files) == 0) {
  stop("Не найдено файлов .xlsx в папке '", input_dir, "'. Проверьте путь!")
}

cat("Найдено файлов для анализа:", length(raw_files), "\n")
cat("Параллельных воркеров:", n_workers, "\n\n")

plan(multisession, workers = n_workers)

t0 <- Sys.time()
with_progress({
  p <- progressor(along = raw_files)
  results <- future_lapply(raw_files, function(fp) process_file(fp, p), future.seed = TRUE)
})

plan(sequential)

master_kegg_list <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$kegg))
master_go_list   <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$go))

if (length(master_kegg_list) > 0) {
  master_kegg <- bind_rows(master_kegg_list) %>%
    mutate(Cell_Type_Clean = gsub(cell_type_clean_pattern, "", Cell_Type)) %>%
    select(Cell_Type, Cell_Type_Clean, ID, Description, GeneRatio, GeneRatio_num,
           BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
  
  wb_kegg <- createWorkbook()
  addWorksheet(wb_kegg, "KEGG_Master")
  writeData(wb_kegg, "KEGG_Master", master_kegg)
  addFilter(wb_kegg, "KEGG_Master", row = 1, cols = 1:ncol(master_kegg))
  freezePane(wb_kegg, "KEGG_Master", firstRow = TRUE)
  
  saveWorkbook(wb_kegg, file = file.path(output_dir, "MASTER_KEGG_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Мастер-таблица KEGG сохранена (", nrow(master_kegg), "строк).\n\n")
} else {
  cat("Значимые данные KEGG не обнаружены.\n\n")
}

if (length(master_go_list) > 0) {
  master_go <- bind_rows(master_go_list) %>%
    mutate(Cell_Type_Clean = gsub(cell_type_clean_pattern, "", Cell_Type)) %>%
    select(Cell_Type, Cell_Type_Clean, ID, Description, ONTOLOGY, GeneRatio, GeneRatio_num,
           BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
  
  wb_go <- createWorkbook()
  addWorksheet(wb_go, "GO_MP_Master")
  writeData(wb_go, "GO_MP_Master", master_go)
  addFilter(wb_go, "GO_MP_Master", row = 1, cols = 1:ncol(master_go))
  freezePane(wb_go, "GO_MP_Master", firstRow = TRUE)
  
  saveWorkbook(wb_go, file = file.path(output_dir, "MASTER_GO_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Мастер-таблица GO сохранена (", nrow(master_go), "строк).\n\n")
} else {
  stop("Данные GO/MP не были собраны.")
}
cat("Время обработки всех файлов:", round(difftime(Sys.time(), t0, units = "secs"), 1), "сек.\n\n")
cat("Готово. Мастер-таблицы сохранены в папку:", output_dir, "\n")
