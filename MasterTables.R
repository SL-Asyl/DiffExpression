cran_packages <- c("openxlsx", "enrichR", "dplyr", "tidyr", "stringr", "BiocManager",
                   "future", "ps")
bioc_packages <- c("clusterProfiler", "AnnotationDbi", "org.Mm.eg.db", "GO.db", "BiocParallel")

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
library(GO.db)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(BiocParallel)
library(future)
# enrichR подключается ниже отдельно, с повторами (при подключении он проверяет сервер)

# fgsea в главной сессии — последовательно (в воркерах выставляется отдельно)
BiocParallel::register(BiocParallel::SerialParam())

# Готовые таблицы GO передаются воркерам; поднимаем лимит на размер передаваемых объектов
options(future.globals.maxSize = 3 * 1024^3)

# =============================================================================
# ВЫБОР ПАПОК
# =============================================================================
# Диалог выбора папки: на Windows — системный (utils::choose.dir),
# на других ОС — через tcltk, если нет графики — ввод пути в консоли.
pick_dir <- function(caption, default = "") {
  path <- NA_character_
  if (.Platform$OS.type == "windows") {
    path <- utils::choose.dir(default = default, caption = caption)
  } else if (capabilities("tcltk") && requireNamespace("tcltk", quietly = TRUE)) {
    path <- tcltk::tk_choose.dir(default = default, caption = caption)
  } else {
    cat(caption, "\nВведите путь: ")
    path <- readLines(con = "stdin", n = 1)
  }
  if (length(path) == 0 || is.na(path) || !nzchar(path)) return(NA_character_)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

input_dir <- pick_dir("Выберите папку с .xlsx-файлами для анализа")
if (is.na(input_dir) || !dir.exists(input_dir)) {
  stop("Папка с входными файлами не выбрана. Анализ отменён.")
}
cat("Входная папка:", input_dir, "\n")

# Папка для результатов: если в диалоге нажать «Отмена», создаётся
# соседняя папка "<имя входной папки>_enrichment"
default_outdir <- file.path(dirname(input_dir), paste0(basename(input_dir), "_enrichment"))
outdir <- pick_dir("Выберите папку для результатов (Отмена — создать рядом с входной)",
                   default = dirname(input_dir))
if (is.na(outdir)) outdir <- default_outdir
if (normalizePath(outdir, winslash = "/", mustWork = FALSE) == input_dir) {
  # результаты в той же папке перепутались бы с входными .xlsx при повторном запуске
  outdir <- default_outdir
  cat("Папка результатов совпала с входной — результаты будут сохранены в:", outdir, "\n")
}
cat("Папка результатов:", outdir, "\n\n")

# =============================================================================
# ПАРАМЕТРЫ
# =============================================================================
# --- ORA (по списку значимых генов) ---
logFC_cutoff <- 1
FDR_cutoff   <- 0.05

mp_db <- "MGI_Mammalian_Phenotype_Level_4_2021"

# --- GSEA (по ВСЕМ генам, без порогов) ---
run_gsea            <- TRUE       # выключатель GSEA-блока
gsea_rank_metric    <- "lfc"      # "lfc"      : ранжирование по log2FC
# "signed_p" : sign(log2FC) * -log10(padj)
gsea_go_ont         <- "BP"       # "BP", "CC", "MF" или "ALL"
gsea_min_genes      <- 200        # минимум генов в ранжированном списке
gsea_min_set        <- 10         # минимальный размер набора генов
gsea_max_set        <- 500        # максимальный размер набора генов
gsea_pvalue_cutoff  <- 0.25       # мягкий порог для сохранения (значимыми считаем p.adjust < 0.05)
gsea_sig_cutoff     <- 0.05       # порог для колонки Significant

# --- Чистые названия типов клеток (колонка Cell_Type_Clean) ---
# NULL  -> общая приставка/окончание у всех файлов определяется и убирается автоматически
#          (например, "6_vs_10_Astro.xlsx", "6_vs_10_Micro.xlsx" -> "Astro", "Micro").
# "..." -> убрать вручную заданный текст (как раньше).
cell_type_clean_pattern <- NULL

# --- Параллельность ---
# Число воркеров ограничено ядрами и СВОБОДНОЙ памятью (считается прямо перед запуском).
# Перед запуском закройте браузер и другие тяжёлые программы.
n_workers_max      <- 6      # жёсткий верхний предел
mem_per_worker_gb  <- 1.5    # оценка памяти на один воркер; если память уходит за 90% — поставьте 2
mem_reserve_gb     <- 3      # сколько памяти оставить системе и главной сессии R

# --- Сеть (только Enrichr) ---
# Принцип «всё или ничего»: если Enrichr не ответил хотя бы на один запрос после всех
# повторов, анализ останавливается целиком (фоновые воркеры тоже) и мастер-таблицы не создаются.
enrichr_max_tries  <- 5          # попыток на каждый запрос к Enrichr
retry_wait_sec     <- 5          # базовая пауза перед повтором (5, 10, 15, 20 с ...)
pause_range_sec    <- c(0.3, 1)  # пауза между последовательными запросами к Enrichr

# --- Кэш аннотаций (KEGG, GO, перевод символов) ---
# Хранится в пользовательской папке кэша R, не зависит от выбранной папки результатов
# (на Windows обычно C:/Users/<имя>/AppData/Local/R/cache/R/enrichment_script)
cache_dir          <- tools::R_user_dir("enrichment_script", which = "cache")
kegg_cache_file    <- file.path(cache_dir, "kegg_mmu_cache.rds")
kegg_cache_days    <- 30
go_cache_file      <- file.path(cache_dir, "go_mmu_cache.rds")   # пересобирается при обновлении org.Mm.eg.db / GO.db

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# =============================================================================
short_pause <- function() Sys.sleep(runif(1, pause_range_sec[1], pause_range_sec[2]))

with_retry <- function(expr_fun, label, tries = 3) {
  for (attempt in seq_len(tries)) {
    res <- tryCatch(expr_fun(), error = function(e) {
      message("  [", label, "] попытка ", attempt, "/", tries,
              " не удалась: ", conditionMessage(e))
      NULL
    })
    if (!is.null(res)) return(res)
    if (attempt < tries) Sys.sleep(retry_wait_sec * attempt)
  }
  NULL
}

# Число из текста с возможной десятичной запятой (как в GSEA-скрипте)
to_num <- function(x) as.numeric(sub(",", ".", trimws(as.character(x)), fixed = TRUE))

# Автоматическая очистка названий: убирает общую для ВСЕХ файлов приставку и окончание,
# обрезая их по разделителю (_ пробел - .), чтобы не откусить кусок самого названия.
clean_cell_type_names <- function(x, manual_pattern = NULL) {
  if (!is.null(manual_pattern) && nzchar(manual_pattern)) {
    return(gsub(manual_pattern, "", x, fixed = TRUE))
  }
  if (length(x) < 2) return(x)
  
  common_prefix <- function(s) {
    ch <- strsplit(s, "")
    n  <- min(lengths(ch))
    i  <- 0
    while (i < n && length(unique(vapply(ch, `[`, character(1), i + 1))) == 1) i <- i + 1
    substr(s[1], 1, i)
  }
  rev_str <- function(s) vapply(strsplit(s, ""), function(z) paste(rev(z), collapse = ""), character(1))
  
  pre <- sub("[^_ .-]*$", "", common_prefix(x))                    # до последнего разделителя включительно
  suf <- sub("^[^_ .-]*", "", rev_str(common_prefix(rev_str(x))))  # с первого разделителя
  
  out <- substr(x, nchar(pre) + 1, nchar(x) - nchar(suf))
  bad <- !nzchar(trimws(out))
  out[bad] <- x[bad]
  out
}

# =============================================================================
# ПОДКЛЮЧЕНИЕ enrichR (один раз, только в главной сессии)
# =============================================================================
enrichr_attached <- isTRUE(with_retry(function() {
  suppressPackageStartupMessages(library(enrichR))
  TRUE
}, label = "подключение enrichR", tries = enrichr_max_tries))

if (!enrichr_attached || !isTRUE(getOption("enrichR.live", TRUE))) {
  stop("Сервер Enrichr (maayanlab.cloud) недоступен. Без MP-анализа обработка не выполняется.\n",
       "Проверьте, открывается ли https://maayanlab.cloud/Enrichr/ в браузере, и запустите скрипт ещё раз.")
}

# =============================================================================
# АННОТАЦИИ: готовятся ОДИН РАЗ в главной сессии и раздаются воркерам
# =============================================================================
# --- KEGG ---
load_kegg_mmu <- function() {
  if (file.exists(kegg_cache_file)) {
    age_days <- as.numeric(difftime(Sys.time(), file.mtime(kegg_cache_file), units = "days"))
    if (age_days < kegg_cache_days) {
      cat("KEGG: используется кэш (", round(age_days, 1), "дн.)\n")
      return(readRDS(kegg_cache_file))
    }
  }
  
  cat("KEGG: скачивание данных для mmu...\n")
  dl <- getFromNamespace("download_KEGG", "clusterProfiler")
  raw <- with_retry(function() dl("mmu"), label = "download KEGG")
  if (is.null(raw)) stop("Не удалось скачать KEGG. Проверьте подключение к интернету.")
  
  t2g <- as.data.frame(raw$KEGGPATHID2EXTID, stringsAsFactors = FALSE)
  t2n <- as.data.frame(raw$KEGGPATHID2NAME,  stringsAsFactors = FALSE)
  colnames(t2g) <- c("term", "gene")
  colnames(t2n) <- c("term", "name")
  t2g$gene <- sub("^mmu:", "", t2g$gene)                      # для mmu ID гена KEGG = Entrez
  t2n$name <- sub(" - Mus musculus.*$", "", t2n$name)         # убрать хвост с названием вида
  
  kegg <- list(t2g = t2g, t2n = t2n)
  saveRDS(kegg, kegg_cache_file)
  cat("KEGG: загружено путей:", length(unique(t2g$term)), "\n")
  kegg
}

# --- GO + перевод SYMBOL <-> ENTREZ ---
# Тот же запрос (GOALL + ONTOLOGYALL по всем Entrez), что clusterProfiler делает внутри
# enrichGO/gseGO, но один раз, а не в каждом воркере.
load_go_mmu <- function() {
  ver <- paste(packageVersion("org.Mm.eg.db"), packageVersion("GO.db"))
  if (file.exists(go_cache_file)) {
    obj <- readRDS(go_cache_file)
    if (identical(obj$version, ver)) { cat("GO: используется кэш\n"); return(obj) }
  }
  
  cat("GO: построение аннотации (один раз, 1-3 мин)...\n")
  kk <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "ENTREZID")
  anno <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = kk,
                                                 columns = c("GOALL", "ONTOLOGYALL"),
                                                 keytype = "ENTREZID"))
  anno <- unique(anno[!is.na(anno$GOALL), c("GOALL", "ENTREZID", "ONTOLOGYALL")])
  
  t2g_by_ont <- lapply(c(BP = "BP", CC = "CC", MF = "MF"), function(o) {
    a <- anno[anno$ONTOLOGYALL == o, ]
    data.frame(term = a$GOALL, gene = a$ENTREZID, stringsAsFactors = FALSE)
  })
  
  term_ont <- tapply(anno$ONTOLOGYALL, anno$GOALL, `[`, 1)
  term_ont <- setNames(as.character(term_ont), names(term_ont))
  
  tn <- suppressMessages(AnnotationDbi::select(GO.db, keys = names(term_ont),
                                               columns = "TERM", keytype = "GOID"))
  t2n <- data.frame(term = tn$GOID, name = tn$TERM, stringsAsFactors = FALSE)
  
  # перевод символов: та же логика, что mapIds(..., multiVals = "first")
  sym2entrez <- suppressMessages(mapIds(org.Mm.eg.db,
                                        keys = AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL"),
                                        keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
  entrez2sym <- suppressMessages(mapIds(org.Mm.eg.db, keys = kk,
                                        keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))
  
  obj <- list(version = ver, t2g_by_ont = t2g_by_ont, t2n = t2n, term_ont = term_ont,
              sym2entrez = sym2entrez, entrez2sym = entrez2sym)
  saveRDS(obj, go_cache_file)
  rm(anno); invisible(gc(verbose = FALSE))
  obj
}

kegg_data  <- load_kegg_mmu()
kegg_t2g   <- kegg_data$t2g
kegg_t2n   <- kegg_data$t2n

go_data    <- load_go_mmu()
go_t2g     <- go_data$t2g_by_ont   # список: BP, CC, MF
go_t2n     <- go_data$t2n
go_ont     <- go_data$term_ont     # GO ID -> онтология
sym2entrez <- go_data$sym2entrez
entrez2sym <- go_data$entrez2sym
rm(kegg_data, go_data); invisible(gc(verbose = FALSE))

map_symbols_to_entrez <- function(symbols) unname(sym2entrez[symbols])

get_symbols_for_ids <- function(slash_separated_ids) {
  vapply(strsplit(slash_separated_ids, "/", fixed = TRUE), function(ids) {
    s <- unname(entrez2sym[ids])
    s[is.na(s)] <- ids[is.na(s)]                              # не нашлось — оставить как есть
    paste(s, collapse = ", ")
  }, character(1))
}

# =============================================================================
# ЧТЕНИЕ ФАЙЛА (одна функция и для Enrichr-этапа, и для воркеров —
# гарантирует, что списки генов везде одинаковые)
# =============================================================================
read_deg_table <- function(file_path) {
  tb_raw <- tryCatch(openxlsx::read.xlsx(file_path), error = function(e) NULL)
  if (is.null(tb_raw)) { message("  Ошибка при чтении файла ", file_path, ". Пропускаем."); return(NULL) }
  
  tb <- tb_raw %>%
    dplyr::rename(gene = names, log2FoldChange = logfoldchanges, padj = pvals_adj) %>%
    dplyr::select(gene, log2FoldChange, padj) %>%
    dplyr::mutate(log2FoldChange = to_num(log2FoldChange),
                  padj           = to_num(padj)) %>%
    dplyr::filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj))
  
  if (nrow(tb) == 0) { message("  Файл пустой после очистки NA: ", basename(file_path)); return(NULL) }
  
  tb$ENTREZID <- map_symbols_to_entrez(tb$gene)
  tb <- tb %>% dplyr::filter(!is.na(ENTREZID))
  if (nrow(tb) == 0) { message("  Не удалось сопоставить ENTREZ ID: ", basename(file_path)); return(NULL) }
  
  tb
}

up_subset   <- function(tb) tb[tb$padj < FDR_cutoff & tb$log2FoldChange >=  logFC_cutoff, ]
down_subset <- function(tb) tb[tb$padj < FDR_cutoff & tb$log2FoldChange <= -logFC_cutoff, ]

# =============================================================================
# ORA
# =============================================================================
run_kegg_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for KEGG ORA: ", label); return(NULL) }
  tryCatch(
    enricher(gene = genes, universe = universe_genes,
             TERM2GENE = kegg_t2g, TERM2NAME = kegg_t2n,
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH"),
    error = function(e) { message("  [KEGG ", label, "] ошибка: ", conditionMessage(e)); NULL })
}

# Как enrichGO(ont = "ALL"): BP, CC, MF считаются по отдельности (своя поправка BH
# и свой universe для каждой), затем результаты объединяются
run_go_ora <- function(genes, universe_genes, label) {
  if (length(genes) < 5) { message("  Too few genes for GO ORA: ", label); return(NULL) }
  lres <- lapply(names(go_t2g), function(o) {
    r <- tryCatch(
      enricher(gene = genes, universe = universe_genes,
               TERM2GENE = go_t2g[[o]], TERM2NAME = go_t2n,
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH"),
      error = function(e) { message("  [GO ", o, " ", label, "] ошибка: ", conditionMessage(e)); NULL })
    if (is.null(r) || nrow(r@result) == 0) return(NULL)
    r@result$ONTOLOGY <- o
    r
  })
  lres <- lres[!vapply(lres, is.null, logical(1))]
  if (length(lres) == 0) return(NULL)
  res <- lres[[1]]
  res@result <- do.call(rbind, lapply(lres, function(r) r@result))
  res
}

# Enrichr: при сбое после всех повторов — ОСТАНОВКА всего анализа (всё или ничего).
# Пропуск возможен только по тому же правилу, что и для KEGG/GO: меньше 5 генов в списке.
run_mp_enrichr <- function(symbols, label) {
  if (length(symbols) < 5) { message("  Too few genes for MP ORA: ", label); return(NULL) }
  short_pause()
  res <- with_retry(function() {
    # enrichR печатает «Uploading data… Querying… Parsing…» через cat — перехватываем
    utils::capture.output(r <- suppressMessages(enrichR::enrichr(symbols, mp_db)))
    if (is.null(r) || !mp_db %in% names(r) || !is.data.frame(r[[mp_db]])) {
      stop("Enrichr вернул пустой или некорректный ответ")
    }
    r
  }, label = paste("Enrichr", label), tries = enrichr_max_tries)
  
  if (is.null(res)) {
    stop("Enrichr не ответил для '", label, "' после ", enrichr_max_tries, " попыток. ",
         "Анализ остановлен, мастер-таблицы НЕ созданы. Запустите скрипт позже.", call. = FALSE)
  }
  res
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
  
  for (col in intersect(c("pvalue", "p.adjust", "qvalue"), colnames(df))) {
    df[[col]] <- signif(as.numeric(df[[col]]), 3)
  }
  
  df$Cell_Type <- cell_type
  df$Direction <- direction
  
  df %>% dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
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
  
  df %>% dplyr::select(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
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
  
  out %>% dplyr::filter(!is.na(ID))
}

# =============================================================================
# GSEA
# =============================================================================
# Ранжированный вектор (имена = Entrez, значения = метрика, по убыванию).
# Берутся ВСЕ гены, порогов нет. Дубликаты Entrez: остаётся ген с макс. |log2FC|.
build_ranks <- function(tb) {
  tb <- tb %>%
    dplyr::filter(is.finite(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  score <- if (gsea_rank_metric == "signed_p") {
    sign(tb$log2FoldChange) * -log10(pmax(tb$padj, 1e-300))
  } else {
    tb$log2FoldChange
  }
  
  sort(setNames(score, tb$ENTREZID), decreasing = TRUE)
}

run_kegg_gsea <- function(ranks, label) {
  tryCatch(
    GSEA(geneList = ranks, TERM2GENE = kegg_t2g, TERM2NAME = kegg_t2n,
         minGSSize = gsea_min_set, maxGSSize = gsea_max_set,
         pvalueCutoff = gsea_pvalue_cutoff, pAdjustMethod = "BH",
         verbose = FALSE, seed = TRUE),
    error = function(e) { message("  [GSEA KEGG ", label, "] ошибка: ", conditionMessage(e)); NULL })
}

# Как gseGO: те же наборы генов (GOALL выбранной онтологии), тот же алгоритм fgsea
run_go_gsea <- function(ranks, label) {
  t2g <- if (gsea_go_ont == "ALL") do.call(rbind, go_t2g) else go_t2g[[gsea_go_ont]]
  tryCatch(
    GSEA(geneList = ranks, TERM2GENE = t2g, TERM2NAME = go_t2n,
         minGSSize = gsea_min_set, maxGSSize = gsea_max_set,
         pvalueCutoff = gsea_pvalue_cutoff, pAdjustMethod = "BH",
         verbose = FALSE, seed = TRUE),
    error = function(e) { message("  [GSEA GO ", label, "] ошибка: ", conditionMessage(e)); NULL })
}

format_gsea_master <- function(res, cell_type, db) {
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  df <- as.data.frame(res@result)
  
  # гены, давшие вклад в сигнал (leading edge) -> символы
  df$Core_Genes <- get_symbols_for_ids(df$core_enrichment)
  df$Core_Count <- lengths(strsplit(df$core_enrichment, "/", fixed = TRUE))
  
  for (col in intersect(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), colnames(df))) {
    df[[col]] <- signif(as.numeric(df[[col]]), 3)
  }
  
  df$ONTOLOGY <- if (db == "GO") unname(go_ont[df$ID]) else NA_character_
  
  df$DB          <- db
  df$Cell_Type   <- cell_type
  df$Direction   <- ifelse(df$NES > 0, "UP", "DOWN")     # NES>0: путь повышен, NES<0: снижен
  df$Significant <- df$p.adjust < gsea_sig_cutoff
  
  df %>% dplyr::select(dplyr::any_of(c("Cell_Type", "DB", "ONTOLOGY", "ID", "Description",
                                       "Direction", "NES", "enrichmentScore", "setSize",
                                       "pvalue", "p.adjust", "qvalue", "Significant",
                                       "Core_Count", "Core_Genes")))
}

# =============================================================================
# ОБРАБОТКА ОДНОГО ФАЙЛА В ВОРКЕРЕ (KEGG, GO, GSEA; MP параллельно считает главная сессия)
# =============================================================================
process_file <- function(file_path) {
  # воркеру не нужны org.Mm.eg.db, GO.db и enrichR: всё приходит готовым
  suppressMessages({
    library(clusterProfiler); library(dplyr); library(tidyr)
    library(stringr); library(openxlsx)
  })
  
  # fgsea внутри воркера — последовательно, без запуска собственного snow-кластера
  BiocParallel::register(BiocParallel::SerialParam())
  
  cell_type <- gsub("\\.xlsx$", "", basename(file_path))
  
  cat("ЗАПУСК АНАЛИЗА ДЛЯ:", cell_type, "\n")
  
  tb <- read_deg_table(file_path)
  if (is.null(tb)) return(NULL)
  
  universe_genes <- unique(tb$ENTREZID)
  
  # ------------------------- ORA -------------------------
  up_genes   <- unique(up_subset(tb)$ENTREZID)
  down_genes <- unique(down_subset(tb)$ENTREZID)
  
  cat("  [", cell_type, "] Up:", length(up_genes), "| Down:", length(down_genes),
      "| Universe:", length(universe_genes), "\n")
  
  kegg_rows <- list(format_kegg_master(run_kegg_ora(up_genes,   universe_genes, paste(cell_type, "UP")),   cell_type, "UP"),
                    format_kegg_master(run_kegg_ora(down_genes, universe_genes, paste(cell_type, "DOWN")), cell_type, "DOWN"))
  go_rows   <- list(format_go_master(run_go_ora(up_genes,   universe_genes, paste(cell_type, "UP")),   cell_type, "UP"),
                    format_go_master(run_go_ora(down_genes, universe_genes, paste(cell_type, "DOWN")), cell_type, "DOWN"))
  
  kegg_rows <- kegg_rows[!vapply(kegg_rows, is.null, logical(1))]
  go_rows   <- go_rows[!vapply(go_rows, is.null, logical(1))]
  
  # ------------------------- GSEA ------------------------
  gsea_rows <- list()
  if (run_gsea) {
    ranks <- build_ranks(tb)
    if (length(ranks) < gsea_min_genes) {
      message("  [", cell_type, "] Слишком мало генов для GSEA (", length(ranks), "). Пропускаем GSEA.")
    } else {
      cat("  [", cell_type, "] GSEA: генов в ранжированном списке:", length(ranks), "\n")
      gsea_rows <- list(format_gsea_master(run_kegg_gsea(ranks, cell_type), cell_type, "KEGG"),
                        format_gsea_master(run_go_gsea(ranks, cell_type),   cell_type, "GO"))
      gsea_rows <- gsea_rows[!vapply(gsea_rows, is.null, logical(1))]
    }
  }
  
  cat("  Анализ для", cell_type, "успешно завершен.\n\n")
  out <- list(kegg = kegg_rows, go = go_rows, gsea = gsea_rows)
  rm(tb); invisible(gc(verbose = FALSE))   # освободить память воркера перед следующим файлом
  out
}

# =============================================================================
# ЗАПУСК
# =============================================================================
raw_files <- list.files(path = input_dir, pattern = "\\.xlsx$", full.names = TRUE)
raw_files <- raw_files[!grepl("^~\\$", basename(raw_files))]   # временные файлы открытого Excel
if (length(raw_files) == 0) {
  stop("Не найдено файлов .xlsx в папке '", input_dir, "'. Проверьте путь!")
}

cell_types <- gsub("\\.xlsx$", "", basename(raw_files))
clean_map  <- setNames(clean_cell_type_names(cell_types, cell_type_clean_pattern), cell_types)
cat("Найдено файлов для анализа:", length(raw_files), "\n")
cat("Пример чистых названий:", paste(head(unique(clean_map), 5), collapse = ", "), "\n\n")

t0 <- Sys.time()

# Число воркеров — по свободной памяти ПОСЛЕ загрузки аннотаций в главную сессию
avail_gb <- tryCatch(ps::ps_system_memory()$avail / 1024^3, error = function(e) NA_real_)
n_by_mem <- if (is.na(avail_gb)) 2 else floor((avail_gb - mem_reserve_gb) / mem_per_worker_gb)
n_workers <- max(1, min(n_workers_max, future::availableCores(omit = 1), n_by_mem, length(raw_files)))
cat("Свободно памяти:", if (is.na(avail_gb)) "не удалось определить" else paste(round(avail_gb, 1), "ГБ"),
    "-> воркеров:", n_workers, "\n\n")

# =============================================================================
# Две части анализа идут ОДНОВРЕМЕННО:
#  - фоновые воркеры считают KEGG, GO и GSEA (тяжёлая часть);
#  - главная сессия тем временем по очереди отправляет запросы в Enrichr (MP).
# Прогресс обеих частей показывается в одной обновляемой строке.
# Всё или ничего: если Enrichr окончательно не ответил, фоновые воркеры
# останавливаются и мастер-таблицы не создаются.
# =============================================================================
n_files      <- length(raw_files)
progress_dir <- tempfile("enrichment_progress_")   # воркеры кладут сюда отметку о каждом готовом файле
dir.create(progress_dir)

text_bar <- function(done, total, width = 12) {
  k <- if (total > 0) round(width * done / total) else width
  paste0("[", strrep("#", k), strrep("-", width - k), "]")
}

ui <- new.env()          # состояние строки индикатора
ui$mp   <- 0
ui$line <- ""

draw_progress <- function(mp_done = ui$mp, enr_done = enr_done_count()) {
  ui$mp <- mp_done
  el <- as.integer(difftime(Sys.time(), t0, units = "secs"))
  line <- sprintf("MP-фенотипы %s %d/%d | KEGG/GO/GSEA %s %d/%d | %02d:%02d",
                  text_bar(mp_done, n_files), mp_done, n_files,
                  text_bar(enr_done, n_files), enr_done, n_files, el %/% 60, el %% 60)
  pad <- max(0, nchar(ui$line) - nchar(line))          # затереть хвост прошлой строки
  cat("\r", line, strrep(" ", pad), sep = "")
  ui$line <- line
  utils::flush.console()
}

clear_progress <- function() {
  if (nzchar(ui$line)) cat("\r", strrep(" ", nchar(ui$line)), "\r", sep = "")
  ui$line <- ""
}

with_progress_line <- function(expr) {
  withCallingHandlers(expr,
                      message = function(m) {
                        clear_progress(); cat(conditionMessage(m)); draw_progress()
                        invokeRestart("muffleMessage")
                      },
                      warning = function(w) {
                        clear_progress(); cat("Предупреждение:", conditionMessage(w), "\n"); draw_progress()
                        invokeRestart("muffleWarning")
                      })
}

enr_done_count <- function() length(list.files(progress_dir, pattern = "\\.done$"))

stop_all <- function(msg) {
  if (nzchar(ui$line)) cat("\n")
  plan(sequential)                          # останавливает фоновые воркеры
  unlink(progress_dir, recursive = TRUE)
  stop(msg, call. = FALSE)
}

# ------------------------- Запуск фоновых воркеров (KEGG, GO, GSEA) -----------
plan(multisession, workers = n_workers)

# файлы раздаются вперемешку (1,4,7… / 2,5,8… / …) для равномерной нагрузки;
# по одной задаче на воркер, чтобы тяжёлые таблицы GO передавались каждому один раз
chunks <- split(seq_len(n_files), rep_len(seq_len(n_workers), n_files))
futs <- lapply(chunks, function(idx) {
  future({
    lapply(idx, function(i) {
      r <- process_file(raw_files[i])
      file.create(file.path(progress_dir, paste0(i, ".done")))   # отметка для индикатора
      r
    })
  }, seed = TRUE, stdout = FALSE)
})

# ------------------------- Параллельно: MP через Enrichr ----------------------
mp_rows_by_file <- vector("list", n_files)
draw_progress(0, enr_done_count())
for (i in seq_len(n_files)) {
  ct <- cell_types[i]
  rows_i <- with_progress_line(tryCatch({
    tb <- read_deg_table(raw_files[i])
    if (is.null(tb)) {
      list()                                # нечитаемый файл пропускается целиком и воркером
    } else {
      rows <- list(format_mp_master(run_mp_enrichr(unique(up_subset(tb)$gene),   paste(ct, "UP")),   ct, "UP"),
                   format_mp_master(run_mp_enrichr(unique(down_subset(tb)$gene), paste(ct, "DOWN")), ct, "DOWN"))
      rows[!vapply(rows, is.null, logical(1))]
    }
  }, error = function(e) stop_all(conditionMessage(e))))
  mp_rows_by_file[i] <- list(rows_i)
  draw_progress(i, enr_done_count())
}

# ------------------------- Ожидание фоновых воркеров ---------------------------
repeat {
  draw_progress(n_files, enr_done_count())
  if (all(vapply(futs, resolved, logical(1)))) break
  Sys.sleep(1)
}
draw_progress(n_files, enr_done_count())
cat("\n\n")

results <- vector("list", n_files)
for (k in seq_along(futs)) {
  vals <- tryCatch(value(futs[[k]]),
                   error = function(e) stop_all(paste("Ошибка в фоновом воркере:", conditionMessage(e))))
  results[chunks[[k]]] <- vals
}

plan(sequential)
unlink(progress_dir, recursive = TRUE)

master_kegg_list <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$kegg))
master_gsea_list <- do.call(c, lapply(results, function(r) if (is.null(r)) list() else r$gsea))
# GO и MP объединяются по файлам в том же порядке, что и раньше: GO UP, GO DOWN, MP UP, MP DOWN
master_go_list   <- do.call(c, lapply(seq_along(results), function(i) {
  if (is.null(results[[i]])) list() else c(results[[i]]$go, mp_rows_by_file[[i]])
}))

# ------------------------- MASTER KEGG (ORA) -------------------------
if (length(master_kegg_list) > 0) {
  master_kegg <- bind_rows(master_kegg_list) %>%
    dplyr::mutate(Cell_Type_Clean = unname(clean_map[Cell_Type])) %>%
    dplyr::select(Cell_Type, Cell_Type_Clean, ID, Description, GeneRatio, GeneRatio_num,
                  BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
  
  wb_kegg <- createWorkbook()
  addWorksheet(wb_kegg, "KEGG_Master")
  writeData(wb_kegg, "KEGG_Master", master_kegg)
  addFilter(wb_kegg, "KEGG_Master", row = 1, cols = 1:ncol(master_kegg))
  freezePane(wb_kegg, "KEGG_Master", firstRow = TRUE)
  
  saveWorkbook(wb_kegg, file = file.path(outdir, "MASTER_KEGG_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Мастер-таблица KEGG сохранена (", nrow(master_kegg), "строк).\n\n")
} else {
  cat("Значимые данные KEGG не обнаружены.\n\n")
}

# ------------------------- MASTER GO/MP (ORA) ------------------------
if (length(master_go_list) > 0) {
  master_go <- bind_rows(master_go_list) %>%
    dplyr::mutate(Cell_Type_Clean = unname(clean_map[Cell_Type])) %>%
    dplyr::select(Cell_Type, Cell_Type_Clean, ID, Description, ONTOLOGY, GeneRatio, GeneRatio_num,
                  BgRatio, pvalue, p.adjust, qvalue, Count, Direction, Genes)
  
  wb_go <- createWorkbook()
  addWorksheet(wb_go, "GO_MP_Master")
  writeData(wb_go, "GO_MP_Master", master_go)
  addFilter(wb_go, "GO_MP_Master", row = 1, cols = 1:ncol(master_go))
  freezePane(wb_go, "GO_MP_Master", firstRow = TRUE)
  
  saveWorkbook(wb_go, file = file.path(outdir, "MASTER_GO_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Мастер-таблица GO сохранена (", nrow(master_go), "строк).\n\n")
} else {
  cat("ВНИМАНИЕ: данные GO/MP (ORA) не были собраны.\n\n")
}

# ------------------------- MASTER GSEA -------------------------------
if (run_gsea) {
  if (length(master_gsea_list) > 0) {
    master_gsea <- bind_rows(master_gsea_list) %>%
      dplyr::mutate(Cell_Type_Clean = unname(clean_map[Cell_Type])) %>%
      dplyr::select(Cell_Type, Cell_Type_Clean, dplyr::everything()) %>%
      dplyr::arrange(Cell_Type, DB, p.adjust)
    
    wb_gsea <- createWorkbook()
    addWorksheet(wb_gsea, "GSEA_Master")
    writeData(wb_gsea, "GSEA_Master", master_gsea)
    addFilter(wb_gsea, "GSEA_Master", row = 1, cols = 1:ncol(master_gsea))
    freezePane(wb_gsea, "GSEA_Master", firstRow = TRUE)
    
    saveWorkbook(wb_gsea, file = file.path(outdir, "MASTER_GSEA_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
    cat("-> Мастер-таблица GSEA сохранена (", nrow(master_gsea), "строк, из них p.adjust <",
        gsea_sig_cutoff, ":", sum(master_gsea$Significant, na.rm = TRUE), ").\n\n")
  } else {
    cat("Наборы GSEA не найдены.\n\n")
  }
}

cat("Время обработки всех файлов:", round(difftime(Sys.time(), t0, units = "secs"), 1), "сек.\n\n")
cat("Готово. Мастер-таблицы сохранены в папку:", outdir, "\n")