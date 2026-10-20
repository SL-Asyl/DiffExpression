# =============================================================================
#  engine.R — общий модульный движок: загрузка медиан, узкие тепловые карты,
#  GO-аннотация ген->термины, сохранение/чтение списков генов.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[read_xlsx_safe, gene_median_file, deg_file, go_all_file])

#' @export
load_med <- function(tissue)
  setNames(lapply(CONTRASTS4, function(cc) gene_median_file(deg_file(tissue, cc))), CONTRASTS4)

#' @export
draw_narrow_panel <- function(genes, med, out_path, title_text) {
  df <- expand.grid(gene = genes, ci = seq_along(CONTRASTS4), stringsAsFactors = FALSE)
  df$lfc <- mapply(function(g, ci) { d <- med[[CONTRASTS4[ci]]]
  if (g %in% names(d)) as.numeric(d[[g]]) else NA_real_ }, df$gene, df$ci)
  df$contrast <- factor(MOD_COL_ORDER[df$ci], levels = MOD_COL_ORDER)
  df$gene     <- factor(df$gene, levels = rev(genes))
  df$lfc_cap  <- pmax(pmin(df$lfc, MOD_CAP), -MOD_CAP)
  df$lab      <- ifelse(is.na(df$lfc), "–", sprintf("%+.1f", df$lfc))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = contrast, y = gene, fill = lfc_cap)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4, width = 0.98, height = 0.98) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 1.9, colour = "grey10") +
    ggplot2::scale_fill_gradient2(low = "#2A78D6", mid = "#F0EFEC", high = "#E34948",
                                  midpoint = 0, limits = c(-MOD_CAP, MOD_CAP), na.value = "grey95", name = "median\nlog2FC") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "bold.italic", size = 6.5),
                   axis.text.x.top = ggplot2::element_text(size = 5.5),
                   plot.title  = ggplot2::element_text(face = "bold", size = 7, hjust = 0.5),
                   legend.title = ggplot2::element_text(size = 5.5), legend.text = ggplot2::element_text(size = 5),
                   legend.key.width = ggplot2::unit(0.18, "cm"), legend.key.height = ggplot2::unit(0.35, "cm"))
  ggplot2::ggsave(out_path, p, width = MOD_FIG_W, height = MOD_ROW_H * length(genes) + 1.3,
                  dpi = 300, bg = "white", limitsize = FALSE)
  cat("  сохранено:", basename(out_path), "(", length(genes), "генов )\n")
}

#' @export
heatmap_module <- function(genes, med, title, out_file, rank_cc) {
  if (length(genes) == 0) { cat("  пустой модуль, пропуск:", out_file, "\n"); return(invisible()) }
  rk  <- med[[rank_cc]]
  val <- vapply(genes, function(g) if (g %in% names(rk)) abs(rk[[g]]) else 0, numeric(1))
  genes <- genes[order(-val)]
  n <- length(genes); n_panels <- max(1, ceiling(n / MOD_GENES_PER_PANEL))
  base_name <- sub("\\.png$", "", out_file); sfx <- c("a","b","c","d","e","f"); ru <- c("а","б","в","г","д","е")
  for (p_i in seq_len(n_panels)) {
    from <- (p_i - 1) * MOD_GENES_PER_PANEL + 1; to <- min(p_i * MOD_GENES_PER_PANEL, n)
    gg <- genes[from:to]
    if (n_panels == 1) { out <- out_file; ttl <- title }
    else { out <- paste0(base_name, "_", sfx[p_i], ".png"); ttl <- sprintf("%s (%s): гены %d–%d", title, ru[p_i], from, to) }
    draw_narrow_panel(gg, med, file.path(BASE_DIR, out), ttl)
  }
}

#' @export
build_gene2terms <- function(tissue) {
  gene2terms <- new.env(hash = TRUE)
  for (cc in CONTRASTS4) {
    go <- read_xlsx_safe(go_all_file(tissue, cc))
    go <- go[!is.na(go$ONTOLOGY) & go$ONTOLOGY == "BP", ]
    for (i in seq_len(nrow(go))) {
      desc  <- tolower(as.character(go$Description[i]))
      genes <- unlist(strsplit(as.character(go$Genes[i]), "[/,;[:space:]]+")); genes <- genes[nzchar(genes)]
      for (g in genes) {
        cur <- if (exists(g, envir = gene2terms, inherits = FALSE)) get(g, envir = gene2terms) else character(0)
        assign(g, union(cur, desc), envir = gene2terms)
      }
    }
  }
  gene2terms
}

#' @export
save_gene_list <- function(genes, name) {
  path <- file.path(BASE_DIR, paste0(name, "_genes.txt"))
  writeLines(genes, path); cat("  список генов:", path, "(", length(genes), ")\n")
}

#' @export
load_gene_list <- function(name) {
  path <- file.path(BASE_DIR, paste0(name, "_genes.txt"))
  if (!file.exists(path)) stop("Нет файла списка генов: ", path, " — сначала запустите команду mod_*")
  readLines(path)
}
