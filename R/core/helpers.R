# =============================================================================
#  helpers.R — общие хелперы: установка пакетов, устойчивое чтение .xlsx,
#  пути датасетов, разбор ключей, выбор папки сырых данных.
# =============================================================================

box::use(R/config)
# перенести все настройки из config в это пространство имён (чтобы обращаться
# к ним без префикса config$)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())

#' @export
ensure_pkgs <- function(cran = character(0), bioc = character(0)) {
  miss_cran <- cran[!(cran %in% rownames(installed.packages()))]
  if (length(miss_cran)) install.packages(miss_cran)
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    miss_bioc <- bioc[!(bioc %in% rownames(installed.packages()))]
    if (length(miss_bioc)) BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }
}

#' @export
read_xlsx_safe <- function(path, sheet = NULL) {
  tryCatch({
    if (is.null(sheet)) as.data.frame(readxl::read_excel(path, guess_max = 1e6))
    else                as.data.frame(readxl::read_excel(path, sheet = sheet, guess_max = 1e6))
  }, error = function(e) {
    if (is.null(sheet)) openxlsx::read.xlsx(path) else openxlsx::read.xlsx(path, sheet = sheet)
  })
}

#' @export
gene_median_file <- function(path) {
  d <- read_xlsx_safe(path)
  lfc <- suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(d$log2FoldChange)))))
  g   <- as.character(d$gene)
  ok  <- !is.na(lfc) & !is.na(g) & nzchar(g)
  tapply(lfc[ok], g[ok], median)
}

#' @export
evs_dir <- function(tissue, cc) file.path(BASE_DIR, if (tissue == "HPC") HD[[cc]] else PD[[cc]])
#' @export
deg_file <- function(tissue, cc) file.path(evs_dir(tissue, cc),
                                           paste0(TAG[[cc]], "_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx"))
#' @export
go_all_file <- function(tissue, cc) file.path(evs_dir(tissue, cc),
                                              paste0(TAG[[cc]], "_MASTER_GO_ALL_CELL_TYPES.xlsx"))
#' @export
kegg_file <- function(tissue, cc) file.path(evs_dir(tissue, cc),
                                            paste0(TAG[[cc]], "_MASTER_KEGG_ALL_CELL_TYPES.xlsx"))
#' @export
go_filt_file <- function(tissue, cc) file.path(evs_dir(tissue, cc),
                                               paste0(TAG[[cc]], "_MASTER_GO_FILTERED.xlsx"))
#' @export
gsea_file <- function(tissue, cc) file.path(evs_dir(tissue, cc),
                                            paste0(TAG[[cc]], "_MASTER_GSEA.xlsx"))

#' @export
parse_dataset <- function(key) list(tissue = sub("_.*$", "", key), cc = sub("^[A-Za-z]+_", "", key))

#' @export
resolve_raw_dir <- function(tissue, cc) {
  key <- paste0(tissue, "_", cc)
  d <- if (key %in% names(RAW_DIRS)) RAW_DIRS[[key]] else NA_character_
  if ((is.na(d) || !nzchar(d)) && isTRUE(RAW_CHOOSE_INTERACTIVE)) {
    msg <- sprintf("Выберите папку с сырыми таблицами ДЭ для %s %s", tissue, TAG[[cc]])
    message(msg)
    d <- tryCatch(utils::choose.dir(caption = msg), error = function(e) NA_character_)
    if (is.na(d) || !nzchar(d)) d <- trimws(readline(paste0(msg, "\nВведите путь: ")))
  }
  if (is.na(d) || !nzchar(d))
    stop("Не задана папка сырых данных для ", key,
         ". Заполните RAW_DIRS[['", key, "']] в config.R или включите RAW_CHOOSE_INTERACTIVE.")
  if (!dir.exists(d)) stop("Папка сырых данных не существует: ", d)
  d
}
