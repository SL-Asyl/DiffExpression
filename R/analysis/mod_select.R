# =============================================================================
#  mod_select.R — АВТОМАТИЧЕСКОЕ выделение функциональных модулей через rrvgo.
#  Полностью заменяет прежний ручной отбор генов по ключевым словам.
#
#  Логика: значимые / reversal гены -> GO BP-термины, их содержащие ->
#  семантическое сходство терминов (rrvgo calculateSimMatrix) ->
#  иерархическая кластеризация (hclust + cutree) -> присвоение каждого гена
#  кластеру по большинству его терминов -> тепловые карты и списки генов.
#
#  Команды:
#    hpc()  — HPC, reversal-набор (значим в обоих FS-сравнениях PAIR, с разворотом
#             знака медианного log2FC между линиями), k = MOD_K_HPC.
#    pfc()  — PFC (реверсии нет): отдельно ответ FS у C57Bl/6N (PAIR[1], k = MOD_K_PFC_C57)
#             и у k18-hACE2 (PAIR[2], k = MOD_K_PFC_K18).
#
#  Параметры кластеризации и стиля карт — в config.R (блок MOD_*).
#  Выход в BASE_DIR: Fig_mod_<tag>_<name>.png и module_<tag>_<name>_genes.txt.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, gene_median_file, deg_file, go_all_file])

# медианы log2FC по контрастам тепловой карты (MOD_CONTRASTS)
.med <- function(tissue)
  setNames(lapply(MOD_CONTRASTS, function(cc) gene_median_file(deg_file(tissue, cc))), MOD_CONTRASTS)

# GO BP-термины, содержащие заданные гены (объединение по нескольким контрастам)
.term_genes <- function(tissue, ccs, sig) {
  tg <- list(); tp <- list(); td <- c()
  for (cc in ccs) {
    go <- read_xlsx_safe(go_all_file(tissue, cc)); go <- go[!is.na(go$ONTOLOGY) & go$ONTOLOGY == "BP", ]
    for (i in seq_len(nrow(go))) {
      gid <- go$ID[i]; g <- intersect(unlist(strsplit(as.character(go$Genes[i]), "[/,; ]+")), sig)
      if (!length(g)) next
      tg[[gid]] <- union(tg[[gid]], g)
      tp[[gid]] <- min(c(tp[[gid]], as.numeric(go$p.adjust[i]))); td[gid] <- go$Description[i]
    }
  }
  list(tg = tg, tp = tp, td = td)
}

# человекочитаемое имя кластера по ведущим (самым значимым) GO-терминам
.name <- function(td_sorted) {
  s <- tolower(paste(head(td_sorted, 5), collapse = " | ")); h <- function(k) grepl(k, s)
  if (h("oxidative phosphorylation|respirat|atp synthesis|electron transport|glycoly|tricarboxylic|mitochond")) return("Energy metabolism")
  if (h("synaptic vesicle|vesicle-mediated|exocytosis|neurotransmitter secretion")) return("Synaptic vesicle cycle")
  if (h("cyclic nucleotide|cyclic purine|camp|cgmp")) return("cAMP and nucleotide signaling")
  if (h("calcium|sequestered calcium")) return("Calcium signaling and ion channels")
  if (h("postsynaptic density|postsynaptic specialization|glutamate receptor|ionotropic")) return("Postsynaptic organization and excitability")
  if (h("potassium ion|action potential|ion channel|cation channel|membrane potential|localization to synapse")) return("Ion channels and excitability")
  if (h("cell-cell adhesion|adhesion|axon|projection guidance")) return("Cell adhesion and axon guidance")
  tools::toTitleCase(td_sorted[1])
}

# семантическая кластеризация терминов + присвоение генов кластерам
.cluster <- function(tg, tp, sig, k) {
  terms <- names(tg)[sapply(tg, length) >= 2]
  sm <- rrvgo::calculateSimMatrix(terms, orgdb = ORGANISM_DB, ont = "BP", method = SIM_METHOD); terms <- rownames(sm)
  hc <- hclust(as.dist(1 - sm), method = "average"); cl <- cutree(hc, k = k)
  gm <- sapply(sig, function(ge) {
    tv <- terms[sapply(terms, function(t) ge %in% tg[[t]])]
    if (!length(tv)) return(NA_integer_)
    as.integer(names(sort(table(cl[tv]), decreasing = TRUE))[1])
  })
  list(cl = cl, gm = gm)
}

# тепловая карта модуля (англ. подписи, MOD_CONTRASTS по столбцам)
.draw <- function(genes, med, title, out, rank_cc) {
  if (!length(genes)) return(invisible())
  genes <- genes[order(-abs(med[[rank_cc]][genes]))]; if (MOD_TOP_N > 0) genes <- head(genes, MOD_TOP_N)
  df <- expand.grid(gene = genes, ci = seq_along(MOD_CONTRASTS), stringsAsFactors = FALSE)
  df$lfc <- mapply(function(g, ci) { d <- med[[MOD_CONTRASTS[ci]]]; if (g %in% names(d)) as.numeric(d[[g]]) else NA_real_ }, df$gene, df$ci)
  df$contrast <- factor(MOD_COL_ORDER[df$ci], levels = MOD_COL_ORDER); df$gene <- factor(df$gene, levels = rev(genes))
  df$lfc_cap <- pmax(pmin(df$lfc, MOD_CAP), -MOD_CAP); df$lab <- ifelse(is.na(df$lfc), "–", sprintf("%+.1f", df$lfc))
  p <- ggplot2::ggplot(df, ggplot2::aes(contrast, gene, fill = lfc_cap)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5, width = 0.98, height = 0.98) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.4, colour = "grey10") +
    ggplot2::scale_fill_gradient2(low = MOD_COL_LOW, mid = MOD_COL_MID, high = MOD_COL_HIGH, midpoint = 0,
                                  limits = c(-MOD_CAP, MOD_CAP), na.value = "grey95", name = "median\nlog2FC") +
    ggplot2::scale_x_discrete(position = "top") + ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold.italic", size = 7.5, colour = "grey5"),
      axis.text.x.top = ggplot2::element_text(size = 8.5, colour = "grey15"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5, colour = "grey5"),
      legend.title = ggplot2::element_text(size = 9.5, colour = "grey10"),
      legend.text = ggplot2::element_text(size = 9, colour = "grey10"),
      legend.key.width = ggplot2::unit(0.6, "cm"), legend.key.height = ggplot2::unit(1.25, "cm"),
      legend.box.spacing = ggplot2::unit(0.1, "cm"), legend.margin = ggplot2::margin(0, 0, 0, 2))
  ggplot2::ggsave(out, p, width = MOD_FIG_W, height = max(MOD_ROW_H * length(genes) + 1.3, 4.5),
                  dpi = 300, bg = "white", limitsize = FALSE)
  cat("  ->", basename(out), "(", length(genes), "генов )\n")
}

# общий обработчик: кластеризует набор генов, рисует карты и пишет списки
.run_set <- function(tissue, sig, ccs, k, rank_cc, med, tag) {
  tgo <- .term_genes(tissue, ccs, sig)
  res <- .cluster(tgo$tg, tgo$tp, sig, k)
  for (m in sort(unique(res$cl))) {
    genes <- sort(sig[which(res$gm == m)]); if (length(genes) < 5) next
    mem <- names(res$cl)[res$cl == m]; nm <- .name(tgo$td[mem[order(unlist(tgo$tp[mem]))]])
    slug <- gsub("[^A-Za-z]+", "_", nm)
    .draw(genes, med, sprintf("%s, %s", nm, tissue), file.path(BASE_DIR, sprintf("Fig_mod_%s_%s.png", tag, slug)), rank_cc)
    writeLines(genes, file.path(BASE_DIR, sprintf("module_%s_%s_genes.txt", tag, slug)))
    cat(sprintf("  [%s] %d генов -> %s\n", nm, length(genes), slug))
  }
}

#' @export
hpc <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2"), bioc = c("rrvgo", "org.Mm.eg.db", "GO.db"))
  suppressMessages({ library(readxl); library(openxlsx); library(ggplot2); library(rrvgo) })
  tissue <- "HPC"; med <- .med(tissue)
  a <- med[[PAIR[1]]]; b <- med[[PAIR[2]]]; common <- intersect(names(a), names(b))
  sig <- sort(common[a[common] > 0 & b[common] < 0])   # reversal: вверх в PAIR[1], вниз в PAIR[2]
  cat("[mod_hpc] reversal генов:", length(sig), "\n")
  .run_set(tissue, sig, MOD_CONTRASTS, MOD_K_HPC, PAIR[1], med, "hpc")
}

#' @export
pfc <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2"), bioc = c("rrvgo", "org.Mm.eg.db", "GO.db"))
  suppressMessages({ library(readxl); library(openxlsx); library(ggplot2); library(rrvgo) })
  tissue <- "PFC"; med <- .med(tissue)
  d1 <- read_xlsx_safe(deg_file(tissue, PAIR[1])); s1 <- unique(as.character(d1$gene)); s1 <- s1[!is.na(s1) & nzchar(s1)]
  d2 <- read_xlsx_safe(deg_file(tissue, PAIR[2])); s2 <- unique(as.character(d2$gene)); s2 <- s2[!is.na(s2) & nzchar(s2)]
  cat("[mod_pfc] значимых у C57 (", TAG[[PAIR[1]]], "):", length(s1), "| у k18 (", TAG[[PAIR[2]]], "):", length(s2), "\n")
  .run_set(tissue, s1, PAIR[1], MOD_K_PFC_C57, PAIR[1], med, "pfc_c57")
  .run_set(tissue, s2, PAIR[2], MOD_K_PFC_K18, PAIR[2], med, "pfc_k18")
}
