# =============================================================================
#  common_genes.R — команда `common_genes`: общие гены двух контрастов,
#  широкая таблица log2FC + Спирмен -> common_genes_wide_log2FC.xlsx (PAIR).
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, deg_file])

#' @export
run <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "dplyr", "tidyr"))
  suppressMessages({ library(readxl); library(openxlsx); library(dplyr); library(tidyr) })
  ccA <- PAIR[1]; ccB <- PAIR[2]; tag_A <- ccA; tag_B <- ccB
  out_xlsx <- file.path(BASE_DIR, "common_genes_wide_log2FC.xlsx")
  cat("[common_genes] HPC:", TAG[[ccA]], "vs", TAG[[ccB]], "| scope =", CG_SCOPE, "\n")
  
  read_degs <- function(path) {
    df <- read_xlsx_safe(path, sheet = CG_SHEET)
    req <- c("gene", "log2FoldChange", "Cell_Type")
    if (!all(req %in% names(df))) stop("Нет колонок ", paste(req, collapse = ", "), "; есть: ", paste(names(df), collapse = ", "))
    data.frame(gene = as.character(df$gene),
               cell = sub("^\\d+_vs_\\d+_", "", as.character(df$Cell_Type)),
               lfc  = suppressWarnings(as.numeric(df$log2FoldChange)), stringsAsFactors = FALSE) %>%
      filter(!is.na(lfc), nzchar(gene)) %>% group_by(gene, cell) %>% summarise(lfc = median(lfc), .groups = "drop")
  }
  A <- read_degs(deg_file("HPC", ccA)); B <- read_degs(deg_file("HPC", ccB))
  if (CG_SCOPE == "neuronal") { A <- A %>% filter(grepl(CG_NEURONAL_RX, cell)); B <- B %>% filter(grepl(CG_NEURONAL_RX, cell)) }
  
  ovA <- A %>% group_by(gene) %>% summarise(overall = median(lfc), .groups = "drop")
  ovB <- B %>% group_by(gene) %>% summarise(overall = median(lfc), .groups = "drop")
  common <- sort(intersect(ovA$gene, ovB$gene))
  wideA <- A %>% filter(gene %in% common) %>% pivot_wider(names_from = cell, values_from = lfc)
  wideB <- B %>% filter(gene %in% common) %>% pivot_wider(names_from = cell, values_from = lfc)
  shared_ct <- intersect(setdiff(names(wideA), "gene"), setdiff(names(wideB), "gene"))
  pair_n <- sapply(shared_ct, function(ct)
    sum(!is.na(wideA[[ct]][match(common, wideA$gene)]) & !is.na(wideB[[ct]][match(common, wideB$gene)])))
  shared_ct <- names(sort(pair_n[pair_n >= CG_MIN_PAIRS], decreasing = TRUE))
  
  out <- data.frame(gene = common, stringsAsFactors = FALSE)
  out[[paste0("overall__", tag_A)]] <- ovA$overall[match(common, ovA$gene)]
  out[[paste0("overall__", tag_B)]] <- ovB$overall[match(common, ovB$gene)]
  for (ct in shared_ct) {
    out[[paste0(ct, "__", tag_A)]] <- wideA[[ct]][match(common, wideA$gene)]
    out[[paste0(ct, "__", tag_B)]] <- wideB[[ct]][match(common, wideB$gene)]
  }
  sp <- function(a, b) {
    ok <- !is.na(a) & !is.na(b)
    if (sum(ok) < 3) return(c(rho = NA, p = NA, n = sum(ok)))
    r <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman"))
    c(rho = round(unname(r$estimate), 3), p = signif(r$p.value, 3), n = sum(ok))
  }
  corr <- rbind(data.frame(set = "overall", t(sp(out[[2]], out[[3]]))),
                do.call(rbind, lapply(shared_ct, function(ct)
                  data.frame(set = ct, t(sp(out[[paste0(ct, "__", tag_A)]], out[[paste0(ct, "__", tag_B)]]))))))
  corr <- corr[order(-corr$n), ]
  openxlsx::write.xlsx(list(common_genes = out, correlations = corr), file = out_xlsx, overwrite = TRUE, keepNA = FALSE)
  cat("  сохранено ->", out_xlsx, "| общих генов:", length(common), "| типов:", length(shared_ct), "\n")
  cat("  overall Spearman: rho =", corr$rho[corr$set == "overall"],
      "p =", corr$p[corr$set == "overall"], "n =", corr$n[corr$set == "overall"], "\n")
}
