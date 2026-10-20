# =============================================================================
#  rrvgo.R — команда `rrvgo`: семантическая свёртка GO -> *_MASTER_GO_FILTERED
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, go_all_file, go_filt_file])

#' @export
run <- function(tissue, cc) {
  ensure_pkgs(c("openxlsx", "dplyr"), bioc = c("rrvgo", "org.Mm.eg.db", "GO.db"))
  suppressMessages({ library(openxlsx); library(dplyr); library(rrvgo) })
  
  go_master <- read_xlsx_safe(go_all_file(tissue, cc))
  cat("[rrvgo]", tissue, TAG[[cc]], "- GO строк:", nrow(go_master), "\n")
  
  go_ontologies <- c("BP", "MF", "CC")
  split_go_by_ontology <- function(go_data) {
    result <- list()
    for (ont in unique(go_data$ONTOLOGY)) result[[ont]] <- go_data %>% filter(ONTOLOGY == ont)
    result
  }
  filter_go_semantic <- function(go_subset, ontology_type) {
    if (nrow(go_subset) == 0) return(go_subset)
    go_ids <- go_subset$ID; valid_go <- grepl("^GO:", go_ids)
    if (sum(valid_go) < 2) { cat(" (<2 GO, пропуск)"); return(go_subset) }
    go_terms <- go_ids[valid_go]; pvals <- go_subset$p.adjust[valid_go]; names(pvals) <- go_terms
    sim_matrix <- tryCatch(calculateSimMatrix(go_terms, orgdb = ORGANISM_DB, ont = ontology_type, method = SIM_METHOD),
                           error = function(e) { cat("\n  ОШИБКА сходства", ontology_type, ":", e$message); NULL })
    if (is.null(sim_matrix) || nrow(sim_matrix) < 2) { cat(" (матрица не построена)"); return(go_subset) }
    pvals_matched <- pvals[rownames(sim_matrix)]; pvals_matched <- pvals_matched[!is.na(pvals_matched)]
    reduced <- reduceSimMatrix(sim_matrix, pvals_matched, threshold = SIM_THRESHOLD, orgdb = ORGANISM_DB)
    best_per_cluster <- reduced %>% group_by(parentTerm) %>%
      slice_min(order_by = score, n = 1, with_ties = FALSE) %>% pull(go)
    keep_ids <- unique(c(best_per_cluster, go_ids[!valid_go]))
    go_subset %>% filter(ID %in% keep_ids)
  }
  
  go_split <- split_go_by_ontology(go_master)
  cat("\n  онтологии:", paste(names(go_split), collapse = ", "), "\n")
  go_filtered_list <- list()
  for (ont_name in go_ontologies) {
    if (!ont_name %in% names(go_split)) { cat("  нет", ont_name, "\n"); next }
    cat("  фильтрация", ont_name, "- было:", nrow(go_split[[ont_name]]))
    go_filtered_list[[ont_name]] <- filter_go_semantic(go_split[[ont_name]], ont_name)
    cat(" -> стало:", nrow(go_filtered_list[[ont_name]]), "\n")
  }
  filtered_go_combined <- bind_rows(go_filtered_list)
  if ("MP" %in% names(go_split)) {
    cat("  MP: без фильтрации,", nrow(go_split[["MP"]]), "терминов\n")
    filtered_go_combined <- bind_rows(filtered_go_combined, go_split[["MP"]])
  }
  wb <- createWorkbook(); addWorksheet(wb, "GO_MP_Filtered")
  writeData(wb, "GO_MP_Filtered", filtered_go_combined)
  setColWidths(wb, "GO_MP_Filtered", 1:ncol(filtered_go_combined), widths = "auto")
  addFilter(wb, "GO_MP_Filtered", row = 1, cols = 1:ncol(filtered_go_combined)); freezePane(wb, "GO_MP_Filtered", firstRow = TRUE)
  saveWorkbook(wb, go_filt_file(tissue, cc), overwrite = TRUE)
  cat("  сохранено:", go_filt_file(tissue, cc), "(", nrow(filtered_go_combined), "строк )\n")
}
