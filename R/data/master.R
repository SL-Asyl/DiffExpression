# =============================================================================
#  master.R — команда `master`: ORA (GO/KEGG/MP) -> *_MASTER_GO/KEGG [интернет]
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, resolve_raw_dir, evs_dir, kegg_file, go_all_file])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "enrichR", "dplyr", "tidyr", "stringr",
                "future", "future.apply", "progressr"),
              bioc = c("clusterProfiler", "AnnotationDbi", "org.Mm.eg.db"))
  suppressMessages({
    library(clusterProfiler); library(AnnotationDbi); library(org.Mm.eg.db)
    library(dplyr); library(tidyr); library(stringr); library(openxlsx)
    library(enrichR); library(future); library(future.apply); library(progressr)
  })
  handlers(handler_txtprogressbar())
  
  input_dir <- resolve_raw_dir(tissue, cc)
  outdir    <- evs_dir(tissue, cc)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  mp_db <- MP_DB
  
  .symbol_cache <- new.env(parent = emptyenv())
  .entrez_cache <- new.env(parent = emptyenv())
  
  map_symbols_to_entrez <- function(symbols) {
    symbols <- unique(symbols)
    uncached <- symbols[!vapply(symbols, exists, logical(1), envir = .entrez_cache, inherits = FALSE)]
    if (length(uncached) > 0) {
      ids <- suppressMessages(mapIds(org.Mm.eg.db, keys = uncached, keytype = "SYMBOL",
                                     column = "ENTREZID", multiVals = "first"))
      for (sym in uncached) { val <- ids[[sym]]
      assign(sym, if (is.null(val) || is.na(val)) NA_character_ else val, envir = .entrez_cache) }
    }
    vapply(symbols, function(sym) get(sym, envir = .entrez_cache), character(1), USE.NAMES = TRUE)
  }
  get_symbols_for_ids <- function(slash_separated_ids) {
    all_ids <- unique(unlist(strsplit(slash_separated_ids, "/", fixed = TRUE)))
    is_entrez <- grepl("^[0-9]+$", all_ids); entrez_ids <- all_ids[is_entrez]
    uncached <- entrez_ids[!vapply(entrez_ids, exists, logical(1), envir = .symbol_cache, inherits = FALSE)]
    if (length(uncached) > 0) {
      syms <- suppressMessages(mapIds(org.Mm.eg.db, keys = uncached, keytype = "ENTREZID",
                                      column = "SYMBOL", multiVals = "first"))
      for (id in uncached) { val <- syms[[id]]; assign(id, if (is.null(val) || is.na(val)) id else val, envir = .symbol_cache) }
    }
    non_entrez <- all_ids[!is_entrez]
    for (id in non_entrez) if (!exists(id, envir = .symbol_cache, inherits = FALSE)) assign(id, id, envir = .symbol_cache)
    vapply(slash_separated_ids, function(x) {
      ids <- strsplit(x, "/", fixed = TRUE)[[1]]
      paste(vapply(ids, function(id) get(id, envir = .symbol_cache), character(1)), collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
  }
  jitter_pause <- function() Sys.sleep(runif(1, JITTER_RANGE_SEC[1], JITTER_RANGE_SEC[2]))
  with_retry <- function(expr_fun, label) {
    for (attempt in seq_len(RETRY_MAX_TRIES)) {
      res <- tryCatch(expr_fun(), error = function(e) {
        message("  [", label, "] попытка ", attempt, "/", RETRY_MAX_TRIES, " не удалась: ", conditionMessage(e)); NULL })
      if (!is.null(res)) return(res)
      if (attempt < RETRY_MAX_TRIES) Sys.sleep(RETRY_WAIT_SEC * attempt)
    }
    message("  [", label, "] все попытки исчерпаны."); NULL
  }
  run_kegg_ora <- function(genes, universe_genes, label) {
    if (length(genes) < 5) { message("  Мало генов для KEGG ORA: ", label); return(NULL) }
    jitter_pause()
    with_retry(function() enrichKEGG(gene = genes, organism = "mmu", keyType = "ncbi-geneid",
                                     universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH"),
               label = paste("KEGG", label))
  }
  run_go_ora <- function(genes, universe_genes, label) {
    if (length(genes) < 5) { message("  Мало генов для GO ORA: ", label); return(NULL) }
    jitter_pause()
    with_retry(function() enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL",
                                   universe = universe_genes, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                                   pAdjustMethod = "BH", readable = TRUE), label = paste("GO", label))
  }
  run_mp_enrichr <- function(symbols, label) {
    if (length(symbols) < 5) { message("  Мало генов для MP ORA: ", label); return(NULL) }
    jitter_pause()
    with_retry(function() enrichr(symbols, mp_db), label = paste("Enrichr", label))
  }
  format_kegg_master <- function(res, cell_type, direction) {
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    df <- as.data.frame(res@result); df$Genes <- get_symbols_for_ids(df$geneID)
    gr <- do.call(rbind, strsplit(df$GeneRatio, "/", fixed = TRUE))
    df$GeneRatio_num <- as.numeric(gr[, 1]) / as.numeric(gr[, 2]); df$GeneRatio <- sprintf("%.1f%%", df$GeneRatio_num * 100)
    bg <- do.call(rbind, strsplit(df$BgRatio, "/", fixed = TRUE))
    df$BgRatio <- sprintf("%.1f%%", as.numeric(bg[, 1]) / as.numeric(bg[, 2]) * 100)
    for (col in c("pvalue", "p.adjust", "qvalue")) df[[col]] <- signif(as.numeric(df[[col]]), 3)
    df$Cell_Type <- cell_type; df$Direction <- direction
    df %>% dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
                         Count, Genes, Cell_Type, Direction, GeneRatio_num)
  }
  format_go_master <- function(res, cell_type, direction) {
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    df <- as.data.frame(res@result); df$Genes <- get_symbols_for_ids(df$geneID)
    gr <- do.call(rbind, strsplit(df$GeneRatio, "/", fixed = TRUE))
    df$GeneRatio_num <- as.numeric(gr[, 1]) / as.numeric(gr[, 2]); df$GeneRatio <- sprintf("%.1f%%", df$GeneRatio_num * 100)
    bg <- do.call(rbind, strsplit(df$BgRatio, "/", fixed = TRUE))
    df$BgRatio <- sprintf("%.1f%%", as.numeric(bg[, 1]) / as.numeric(bg[, 2]) * 100)
    for (col in intersect(c("pvalue", "p.adjust", "qvalue"), colnames(df))) df[[col]] <- signif(as.numeric(df[[col]]), 3)
    df$Cell_Type <- cell_type; df$Direction <- direction
    df %>% dplyr::select(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
                         Count, Genes, Cell_Type, Direction, GeneRatio_num)
  }
  format_mp_master <- function(mp_res, cell_type, direction) {
    if (is.null(mp_res) || !mp_db %in% names(mp_res)) return(NULL)
    df <- mp_res[[mp_db]]; if (is.null(df) || nrow(df) == 0) return(NULL)
    id <- str_extract(df$Term, "MP:\\d+"); desc <- str_trim(gsub("\\(?MP:\\d+\\)?", "", df$Term))
    genes <- gsub(";\\s*", ", ", df$Genes); cnt <- vapply(strsplit(df$Genes, ";\\s*"), length, integer(1))
    out <- data.frame(ID = id, Description = desc, ONTOLOGY = "MP", GeneRatio = NA_character_, BgRatio = NA_character_,
                      pvalue = signif(as.numeric(df$P.value), 3), p.adjust = signif(as.numeric(df$Adjusted.P.value), 3),
                      qvalue = NA_real_, Count = cnt, Genes = genes, Cell_Type = cell_type, Direction = direction,
                      GeneRatio_num = NA_real_, stringsAsFactors = FALSE)
    out %>% dplyr::filter(!is.na(ID))
  }
  
  process_file <- function(file_path, p = NULL) {
    suppressMessages({
      library(clusterProfiler); library(AnnotationDbi); library(org.Mm.eg.db)
      library(dplyr); library(tidyr); library(stringr); library(openxlsx); library(enrichR)
    })
    cell_type <- gsub("\\.xlsx$", "", basename(file_path))
    if (!is.null(p)) on.exit(p(sprintf("%s", cell_type)), add = TRUE)
    cat("ЗАПУСК ORA ДЛЯ:", cell_type, "\n")
    tb_raw <- tryCatch(read.xlsx(file_path), error = function(e) NULL)
    if (is.null(tb_raw)) { message("  Ошибка чтения ", file_path); return(NULL) }
    tb <- tb_raw %>% dplyr::rename(gene = names, log2FoldChange = logfoldchanges, padj = pvals_adj) %>%
      dplyr::select(gene, log2FoldChange, padj) %>%
      dplyr::filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj))
    if (nrow(tb) == 0) { message("  Пусто после очистки NA"); return(NULL) }
    tb$ENTREZID <- map_symbols_to_entrez(tb$gene)[tb$gene]
    tb <- tb %>% dplyr::filter(!is.na(ENTREZID))
    if (nrow(tb) == 0) { message("  Нет ENTREZ ID"); return(NULL) }
    universe_genes <- unique(tb$ENTREZID)
    up   <- tb %>% dplyr::filter(padj < FDR_CUTOFF, log2FoldChange >=  LOGFC_CUTOFF)
    down <- tb %>% dplyr::filter(padj < FDR_CUTOFF, log2FoldChange <= -LOGFC_CUTOFF)
    up_genes <- unique(up$ENTREZID); down_genes <- unique(down$ENTREZID)
    up_symbols <- unique(up$gene);   down_symbols <- unique(down$gene)
    cat("  [", cell_type, "] Up:", length(up_genes), "| Down:", length(down_genes),
        "| Universe:", length(universe_genes), "\n")
    kegg_up   <- run_kegg_ora(up_genes, universe_genes, paste(cell_type, "UP"))
    kegg_down <- run_kegg_ora(down_genes, universe_genes, paste(cell_type, "DOWN"))
    go_up     <- run_go_ora(up_genes, universe_genes, paste(cell_type, "UP"))
    go_down   <- run_go_ora(down_genes, universe_genes, paste(cell_type, "DOWN"))
    mp_up     <- run_mp_enrichr(up_symbols, paste(cell_type, "UP"))
    mp_down   <- run_mp_enrichr(down_symbols, paste(cell_type, "DOWN"))
    kegg_rows <- list(format_kegg_master(kegg_up, cell_type, "UP"), format_kegg_master(kegg_down, cell_type, "DOWN"))
    go_rows   <- list(format_go_master(go_up, cell_type, "UP"), format_go_master(go_down, cell_type, "DOWN"),
                      format_mp_master(mp_up, cell_type, "UP"), format_mp_master(mp_down, cell_type, "DOWN"))
    kegg_rows <- kegg_rows[!vapply(kegg_rows, is.null, logical(1))]
    go_rows   <- go_rows[!vapply(go_rows, is.null, logical(1))]
    cat("  Готово:", cell_type, "\n\n")
    list(kegg = kegg_rows, go = go_rows)
  }
  
  raw_files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(raw_files) == 0) stop("Нет .xlsx в ", input_dir)
  cat("[master]", tissue, TAG[[cc]], "- файлов:", length(raw_files), "\n")
  
  plan(multisession, workers = N_WORKERS)
  with_progress({
    p <- progressor(along = raw_files)
    results <- future_lapply(raw_files, function(fp) process_file(fp, p), future.seed = TRUE)
  })
  plan(sequential)
  
  master_kegg_list <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$kegg))
  master_go_list   <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$go))
  
  if (length(master_kegg_list) > 0) {
    master_kegg <- bind_rows(master_kegg_list) %>%
      dplyr::mutate(Cell_Type_Clean = sub("^\\d+_vs_\\d+_", "", Cell_Type)) %>%
      dplyr::select(Cell_Type, Cell_Type_Clean, ID, Description, GeneRatio, GeneRatio_num,
                    BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
    wb_kegg <- createWorkbook(); addWorksheet(wb_kegg, "KEGG_Master")
    writeData(wb_kegg, "KEGG_Master", master_kegg)
    addFilter(wb_kegg, "KEGG_Master", row = 1, cols = 1:ncol(master_kegg)); freezePane(wb_kegg, "KEGG_Master", firstRow = TRUE)
    saveWorkbook(wb_kegg, file = kegg_file(tissue, cc), overwrite = TRUE)
    cat("  KEGG сохранена (", nrow(master_kegg), "строк).\n")
  } else cat("  KEGG не обнаружены.\n")
  
  if (length(master_go_list) > 0) {
    master_go <- bind_rows(master_go_list) %>%
      dplyr::mutate(Cell_Type_Clean = sub("^\\d+_vs_\\d+_", "", Cell_Type)) %>%
      dplyr::select(Cell_Type, Cell_Type_Clean, ID, Description, ONTOLOGY, GeneRatio, GeneRatio_num,
                    BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
    wb_go <- createWorkbook(); addWorksheet(wb_go, "GO_MP_Master")
    writeData(wb_go, "GO_MP_Master", master_go)
    addFilter(wb_go, "GO_MP_Master", row = 1, cols = 1:ncol(master_go)); freezePane(wb_go, "GO_MP_Master", firstRow = TRUE)
    saveWorkbook(wb_go, file = go_all_file(tissue, cc), overwrite = TRUE)
    cat("  GO сохранена (", nrow(master_go), "строк).\n")
  } else stop("Данные GO/MP не собраны.")
}
