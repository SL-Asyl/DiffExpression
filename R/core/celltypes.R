# =============================================================================
#  celltypes.R — классификация типов клеток и загрузка DEG с группировкой
#  (используется в fig_volcano и fig_heatmap_deg).
# =============================================================================

box::use(R/core/helpers[read_xlsx_safe, deg_file])

#' @export
GROUP_ORDER <- c("Glutamatergic", "GABAergic", "Glia", "Microglia")

#' @export
classify_cell_type <- function(ct) {
  if (grepl("[_-]Glut$", ct)) return("Glutamatergic")
  if (grepl("[_-]Gaba$", ct)) return("GABAergic")
  if (grepl("[_-]Inh", ct))   return("GABAergic")
  if (grepl("^L[0-9]", ct) || grepl("[_-]CTX$", ct)) return("Glutamatergic")
  if (grepl("[_-](IT|ET|PT|CT|NP)([_-]CTX)?$", ct)) return("Glutamatergic")
  if (ct %in% c("Pvalb", "Lamp5", "Vip", "Sncg", "Sst")) return("GABAergic")
  if (ct == "Microglia" || ct == "Microglia_NN") return("Microglia")
  if (ct %in% c("Astro", "Oligo", "OPC", "Endo", "Peri") || grepl("[_-]NN$", ct)) return("Glia")
  "Other"
}

#' @export
load_degs_classified <- function(tissue, cc) {
  degs <- read_xlsx_safe(deg_file(tissue, cc))
  degs$Cell_Type <- sub("^\\d+_vs_\\d+_", "", degs$Cell_Type)
  degs$log2FoldChange <- suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(degs$log2FoldChange)))))
  degs$padj <- suppressWarnings(as.numeric(degs$padj))
  degs$Group <- sapply(degs$Cell_Type, classify_cell_type)
  degs
}
