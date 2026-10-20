# =============================================================================
#  graphpad.R — команда `graphpad`: биномиальный + Спирмен по HPC -> .csv Prism
#  Использует PAIR.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, deg_file])

#' @export
run <- function() {
  ensure_pkgs(c("readxl", "openxlsx")); suppressMessages({ library(readxl); library(openxlsx) })
  ccA <- PAIR[1]; ccB <- PAIR[2]
  out_dir <- file.path(BASE_DIR, "Prism_interaction"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("[graphpad] HPC:", TAG[[ccA]], "(C57) vs", TAG[[ccB]], "(k18)\n")
  
  read_deg <- function(path) {
    d <- read_xlsx_safe(path)
    d$log2FoldChange <- suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(d$log2FoldChange)))))
    d$ct <- sub("^[0-9]+_vs_[0-9]+_", "", as.character(d$Cell_Type))
    d[!is.na(d$log2FoldChange) & !is.na(d$gene), c("gene", "ct", "log2FoldChange")]
  }
  A <- read_deg(deg_file("HPC", ccA)); B <- read_deg(deg_file("HPC", ccB))
  
  m <- merge(A, B, by = c("gene", "ct"), suffixes = c("_C57", "_hACE2")); m <- m[order(m$gene, m$ct), ]
  n_pairs <- nrow(m); disc <- sign(m$log2FoldChange_C57) != sign(m$log2FoldChange_hACE2); n_disc <- sum(disc)
  bt <- binom.test(n_disc, n_pairs, p = 0.5)
  
  write.csv(data.frame(category = c("Разнонаправленные (opposite sign)", "Однонаправленные (same sign)"),
                       count = c(n_disc, n_pairs - n_disc)),
            file.path(out_dir, "HPC_FS_binomial_discordance.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(data.frame(gene = m$gene, cell_type = m$ct,
                       `log2FC_C57Bl6N (no FS vs FS)`  = round(m$log2FoldChange_C57, 4),
                       `log2FC_k18hACE2 (no FS vs FS)` = round(m$log2FoldChange_hACE2, 4), check.names = FALSE),
            file.path(out_dir, "HPC_FS_log2FC_paired_pairs.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  sp_pairs <- suppressWarnings(cor.test(m$log2FoldChange_C57, m$log2FoldChange_hACE2, method = "spearman"))
  
  medA <- tapply(A$log2FoldChange, A$gene, median); medB <- tapply(B$log2FoldChange, B$gene, median)
  g <- sort(intersect(names(medA), names(medB)))
  write.csv(data.frame(gene = g,
                       `median_log2FC_C57Bl6N (no FS vs FS)`  = round(as.numeric(medA[g]), 4),
                       `median_log2FC_k18hACE2 (no FS vs FS)` = round(as.numeric(medB[g]), 4), check.names = FALSE),
            file.path(out_dir, "HPC_FS_log2FC_paired_genes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  sp_genes <- suppressWarnings(cor.test(medA[g], medB[g], method = "spearman"))
  
  readme <- c(
    "ДАННЫЕ ДЛЯ GraphPad Prism — взаимодействие «FS x генотип» в HPC (биномиальный + Спирмен)",
    "", "Все общие значимые пары/гены по ВСЕМ типам клеток (без пан-клеточного фильтра).",
    "Источник: *_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx.", "",
    "1) HPC_FS_binomial_discordance.csv — Parts of whole -> Binomial test, доля 0,5.",
    sprintf("   Из %d общих пар %d разнонаправлены; p = %.2g.", n_pairs, n_disc, bt$p.value), "",
    "2) HPC_FS_log2FC_paired_pairs.csv — XY, Correlation -> Spearman.",
    sprintf("   rho = %.2f (p = %.2g, n = %d).", sp_pairs$estimate, sp_pairs$p.value, n_pairs),
    "   HPC_FS_log2FC_paired_genes.csv (медианы по генам):",
    sprintf("   rho = %.2f (p = %.2g, n = %d).", sp_genes$estimate, sp_genes$p.value, length(g)))
  writeLines(readme, file.path(out_dir, "README_HPC_binomial_spearman.txt"), useBytes = TRUE)
  
  cat(sprintf("  ПАРЫ: n=%d; разнонапр.=%d (%.0f%%); биномиальный p=%.3g\n", n_pairs, n_disc, 100 * n_disc / n_pairs, bt$p.value))
  cat(sprintf("  Спирмен (пары): rho=%.3f; p=%.3g\n", sp_pairs$estimate, sp_pairs$p.value))
  cat(sprintf("  Спирмен (гены): rho=%.3f; p=%.3g; n=%d\n", sp_genes$estimate, sp_genes$p.value, length(g)))
  cat("  файлы ->", out_dir, "\n")
}
