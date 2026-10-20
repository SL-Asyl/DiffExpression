# =============================================================================
#  gsea.R — команда `gsea`: GSEA по полным scanpy-таблицам -> *_MASTER_GSEA
#  Читает ВСЕ гены (не фильтрует по порогу). [нужен интернет]
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, resolve_raw_dir, gsea_file])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "dplyr"), bioc = c("clusterProfiler", "org.Mm.eg.db"))
  suppressMessages({ library(clusterProfiler); library(org.Mm.eg.db); library(openxlsx); library(dplyr) })
  
  input_dir <- resolve_raw_dir(tissue, cc)
  outfile   <- gsea_file(tissue, cc)
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  
  run_gsea_one <- function(file_path) {
    cell_type <- sub("\\.xlsx$", "", basename(file_path)); cell_type <- sub("^\\d+_vs_\\d+_", "", cell_type)
    tb <- tryCatch(read.xlsx(file_path), error = function(e) NULL); if (is.null(tb)) return(NULL)
    tb <- data.frame(gene = tb$names, lfc = as.numeric(sub(",", ".", trimws(tb$logfoldchanges))), stringsAsFactors = FALSE)
    tb <- tb[complete.cases(tb), ]; tb <- tb[is.finite(tb$lfc), ]
    if (nrow(tb) < 200) return(NULL)
    map <- suppressWarnings(mapIds(org.Mm.eg.db, keys = tb$gene, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
    tb$entrez <- map[tb$gene]; tb <- tb[!is.na(tb$entrez), ]; tb <- tb[!duplicated(tb$entrez), ]
    ranks <- sort(setNames(tb$lfc, tb$entrez), decreasing = TRUE)
    res <- list()
    kk <- tryCatch(gseKEGG(geneList = ranks, organism = "mmu", pvalueCutoff = 0.25, verbose = FALSE, seed = TRUE),
                   error = function(e) NULL)
    if (!is.null(kk) && nrow(kk@result) > 0) {
      r <- kk@result; r$Cell_Type <- cell_type; r$DB <- "KEGG"
      r$Direction <- ifelse(r$NES > 0, "UP (активирован)", "DOWN (подавлен)"); res$KEGG <- r
    }
    gg <- tryCatch(gseGO(geneList = ranks, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.25, verbose = FALSE, seed = TRUE),
                   error = function(e) NULL)
    if (!is.null(gg) && nrow(gg@result) > 0) {
      r <- setReadable(gg, org.Mm.eg.db, keyType = "ENTREZID")@result
      r$Cell_Type <- cell_type; r$DB <- "GO_BP"
      r$Direction <- ifelse(r$NES > 0, "UP (активирован)", "DOWN (подавлен)"); res$GO <- r
    }
    bind_rows(res)
  }
  
  files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
  cat("[gsea]", tissue, TAG[[cc]], "- файлов:", length(files), "\n")
  master <- bind_rows(lapply(files, function(f) {
    cat("  GSEA:", basename(f), "\n")
    tryCatch(run_gsea_one(f), error = function(e) { cat("    пропуск:", e$message, "\n"); NULL })
  }))
  if (nrow(master) > 0) {
    keep <- c("Cell_Type", "DB", "ID", "Description", "Direction", "NES",
              "setSize", "pvalue", "p.adjust", "core_enrichment")
    master <- master[, intersect(keep, colnames(master))]
    master <- master[order(master$Cell_Type, master$p.adjust), ]
    write.xlsx(master, outfile)
    cat("  сохранено:", outfile, "(", nrow(master), "путей )\n")
  } else cat("  Значимых наборов GSEA не найдено.\n")
}
