# =============================================================================
#  degs.R — команда `degs`: сырые scanpy-таблицы -> *_SIGNIFICANT_DEGs.xlsx
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, resolve_raw_dir, evs_dir, deg_file])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "future.apply", "dplyr"))
  suppressMessages({ library(openxlsx); library(future.apply); library(dplyr) })
  
  input_dir <- resolve_raw_dir(tissue, cc)
  outdir    <- evs_dir(tissue, cc)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) == 0) stop("Нет файлов .xlsx в ", input_dir)
  cat("[degs]", tissue, TAG[[cc]], "- файлов:", length(files), "\n")
  
  process_file <- function(file_path) {
    cell_type <- sub("\\.xlsx$", "", basename(file_path))
    cell_type <- sub("^\\d+_vs_\\d+_", "", cell_type)
    tb_raw <- tryCatch(read.xlsx(file_path), error = function(e) NULL)
    if (is.null(tb_raw)) return(NULL)
    tb <- data.frame(gene = tb_raw$names, log2FoldChange = tb_raw$logfoldchanges,
                     padj = tb_raw$pvals_adj, stringsAsFactors = FALSE)
    tb <- tb[complete.cases(tb), ]
    if (nrow(tb) == 0) return(NULL)
    tb$log2FoldChange <- as.numeric(sub(",", ".", trimws(tb$log2FoldChange)))
    idx <- tb$padj < FDR_CUTOFF & abs(tb$log2FoldChange) >= LOGFC_CUTOFF
    tb <- tb[idx, c("gene", "log2FoldChange", "padj")]
    if (nrow(tb) == 0) return(NULL)
    tb <- tb[order(-tb$log2FoldChange), ]
    tb$Cell_Type <- cell_type
    tb
  }
  
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
  res <- future_lapply(files, process_file)
  plan(sequential)
  res <- Filter(Negate(is.null), res)
  
  if (length(res) == 0) { cat("  Значимых генов не обнаружено.\n"); return(invisible()) }
  master_table <- dplyr::bind_rows(res)
  wb <- createWorkbook(); addWorksheet(wb, "Significant_DEGs")
  writeData(wb, "Significant_DEGs", master_table)
  setColWidths(wb, "Significant_DEGs", 1:ncol(master_table), widths = "auto")
  addFilter(wb, "Significant_DEGs", row = 1, cols = 1:ncol(master_table))
  out <- deg_file(tissue, cc)
  saveWorkbook(wb, out, overwrite = TRUE)
  cat("  сохранено:", out, "(", nrow(master_table), "строк )\n")
}
