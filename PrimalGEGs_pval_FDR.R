required_packages <- c("openxlsx", "future.apply", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(openxlsx)
library(future.apply)
library(dplyr)

input_dir <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/10 vs 13"
outdir    <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 10 vs 13 EVS"

logFC_cutoff <- 1
FDR_cutoff   <- 0.05

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
if (length(files) == 0) stop("Нет файлов .xlsx")

cat("Файлов для обработки:", length(files), "\n")

plan(multisession, workers = max(1, parallel::detectCores() - 1))

process_file <- function(file_path) {
  
  cell_type <- basename(file_path)
  cell_type <- sub(".xlsx$", "", cell_type)
  cell_type <- sub("5_vs_6_", "", cell_type)
  
  tb_raw <- tryCatch(read.xlsx(file_path), error = function(e) NULL)
  if (is.null(tb_raw)) return(NULL)
  
  tb <- data.frame(
    gene = tb_raw$names,
    log2FoldChange = tb_raw$logfoldchanges,
    padj = tb_raw$pvals_adj,
    stringsAsFactors = FALSE
  )
  
  tb <- tb[complete.cases(tb), ]
  if (nrow(tb) == 0) return(NULL)
  
  tb$log2FoldChange <- as.numeric(
    sub(",", ".", trimws(tb$log2FoldChange))
  )
  
  idx <- tb$padj < FDR_cutoff & abs(tb$log2FoldChange) >= logFC_cutoff
  tb <- tb[idx, c("gene", "log2FoldChange", "padj")]
  
  if (nrow(tb) == 0) return(NULL)
  
  tb <- tb[order(-tb$log2FoldChange), ]
  
  tb$Cell_Type <- cell_type
  
  return(tb)
}

res <- future_lapply(files, process_file)
res <- Filter(Negate(is.null), res)

if (length(res) > 0) {
  master_table <- dplyr::bind_rows(res)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Significant_DEGs")
  writeData(wb, "Significant_DEGs", master_table)
  setColWidths(wb, "Significant_DEGs", 1:ncol(master_table), widths = "auto")
  addFilter(wb, "Significant_DEGs", row = 1, cols = 1:ncol(master_table))
  
  output_file <- file.path(outdir, "ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("Файл успешно сохранен:", output_file, "\n")
} else {
  cat("Значимых генов не обнаружено ни в одном файле.\n")
}
