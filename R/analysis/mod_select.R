# =============================================================================
#  mod_select.R — отбор генов модулей (команды mod_hpc / mod_pfc_energy /
#  mod_pfc_excit). Используют PAIR. Сохраняют списки генов в *_genes.txt.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs, read_xlsx_safe, gene_median_file, deg_file, go_all_file])
box::use(R/core/engine[build_gene2terms, save_gene_list])

#' @export
hpc <- function() {
  ensure_pkgs(c("readxl", "openxlsx")); suppressMessages({ library(readxl); library(openxlsx) })
  ccA <- PAIR[1]; ccB <- PAIR[2]
  cat("[mod_hpc] разворот HPC:", TAG[[ccA]], "vs", TAG[[ccB]], "\n")
  a <- gene_median_file(deg_file("HPC", ccA)); b <- gene_median_file(deg_file("HPC", ccB))
  common <- intersect(names(a), names(b)); reversal <- common[a[common] > 0 & b[common] < 0]
  
  gene2terms <- build_gene2terms("HPC")
  terms_of <- function(gene) if (exists(gene, envir = gene2terms, inherits = FALSE)) get(gene, envir = gene2terms) else character(0)
  GROUPS1 <- list(
    "наведение аксонов"    = c("axon guidance","neuron projection guidance","axonogenesis","axon development","ephrin receptor signaling"),
    "адгезия"              = c("cell adhesion","cell-cell adhesion","homophilic cell adhesion","heterophilic cell adhesion","synaptic membrane adhesion"),
    "организация синапсов" = c("synapse organization","synapse assembly","regulation of synapse structure","postsynapse organization","presynapse","dendritic spine","neuron projection morphogenesis"))
  GROUPS2 <- list(
    "ионные каналы"        = c("ion channel","potassium ion","sodium ion","cation transmembrane","monoatomic cation","ligand-gated","voltage-gated","action potential","membrane potential"),
    "Ca-сигнализация"      = c("calcium ion"),
    "глутаматный рецептор" = c("glutamate receptor","ionotropic glutamate","metabotropic glutamate","neurotransmitter receptor"))
  KW1 <- unlist(GROUPS1, use.names = FALSE); KW2 <- unlist(GROUPS2, use.names = FALSE)
  nhits <- function(gene, kws) {
    tl <- terms_of(gene); if (!length(tl)) return(0L); n <- 0L
    for (t in tl) for (k in kws) if (grepl(k, t, fixed = TRUE)) n <- n + 1L
    n
  }
  mod1 <- character(0); mod2 <- character(0)
  for (g in sort(reversal)) {
    h1 <- nhits(g, KW1); h2 <- nhits(g, KW2)
    if (h1 == 0 && h2 == 0) next
    if (h1 >= h2) mod1 <- c(mod1, g) else mod2 <- c(mod2, g)
  }
  cat("  модуль адгезии:", length(mod1), "| ионный:", length(mod2), "\n")
  save_gene_list(mod1, "module_hpc_adhesion"); save_gene_list(mod2, "module_hpc_ion")
}

#' @export
pfc_energy <- function() {
  ensure_pkgs(c("readxl", "openxlsx")); suppressMessages({ library(readxl); library(openxlsx) })
  cc <- PAIR[1]; cat("[mod_pfc_energy] PFC, значимость в", TAG[[cc]], "\n")
  med <- gene_median_file(deg_file("PFC", cc))
  CORE <- c("oxidative phosphorylation","respiratory electron transport","electron transport chain",
            "respiratory chain","nadh dehydrogenase","cytochrome complex",
            "atp synthesis coupled electron transport","aerobic respiration","cellular respiration","oxidoreductase complex")
  go <- read_xlsx_safe(go_all_file("PFC", cc)); core_genes <- character(0)
  for (i in seq_len(nrow(go))) {
    d <- tolower(as.character(go$Description[i]))
    if (any(vapply(CORE, function(k) grepl(k, d, fixed = TRUE), logical(1)))) {
      gs <- unlist(strsplit(as.character(go$Genes[i]), "[/,;[:space:]]+")); core_genes <- union(core_genes, gs[nzchar(gs)])
    }
  }
  module <- sort(intersect(core_genes, names(med)))
  cat("  энергетический модуль PFC:", length(module), "генов\n")
  save_gene_list(module, "module_pfc_energy")
}

#' @export
pfc_excit <- function() {
  ensure_pkgs(c("readxl", "openxlsx")); suppressMessages({ library(readxl); library(openxlsx) })
  cc <- PAIR[2]; cat("[mod_pfc_excit] PFC, значимость в", TAG[[cc]], "\n")
  med <- gene_median_file(deg_file("PFC", cc))
  CORE <- c("action potential","membrane potential","potassium ion transmembrane","sodium ion transmembrane",
            "calcium ion transmembrane","cation channel","anion channel","ion channel","voltage-gated","gated channel",
            "glutamate receptor","ionotropic glutamate","glutamatergic synaptic","postsynaptic density",
            "postsynaptic membrane","regulation of postsynaptic","neurotransmitter receptor","synaptic transmission",
            "regulation of membrane potential","monoatomic cation","monoatomic ion channel")
  go <- read_xlsx_safe(go_all_file("PFC", cc)); core_genes <- character(0)
  for (i in seq_len(nrow(go))) {
    d <- tolower(as.character(go$Description[i]))
    if (any(vapply(CORE, function(k) grepl(k, d, fixed = TRUE), logical(1)))) {
      gs <- unlist(strsplit(as.character(go$Genes[i]), "[/,;[:space:]]+")); core_genes <- union(core_genes, gs[nzchar(gs)])
    }
  }
  module <- sort(intersect(core_genes, names(med)))
  cat("  общий модуль возбудимости PFC:", length(module), "генов\n")
  GLU  <- "^(Gria|Grik|Grin|Grm|Dlg|Dlgap|Shank|Neto|Shisa|Gphn|Nlgn|Homer|Lin7)"
  CHAN <- "^(Kcn|Cacna|Scn|Hcn|Gabr|Ano|Asic|Dpp6|Dpp10|Fgf1[234]|Trp|Clcn|Cnga|Cngb)"
  CHAN_EXTRA <- c("Atp1a2","Atp2b1","Slc8a1","Slc9a9","Slc4a4","Atp13a5","Stk39")
  classify <- function(g) {
    if (grepl(GLU, g, ignore.case = TRUE)) "M2"
    else if (grepl(CHAN, g, ignore.case = TRUE) || g %in% CHAN_EXTRA) "M1" else "M3"
  }
  cls <- vapply(module, classify, character(1))
  cat("  каналы:", sum(cls == "M1"), "| глутамат/PSD:", sum(cls == "M2"), "| адгезия/сигнал:", sum(cls == "M3"), "\n")
  save_gene_list(module[cls == "M1"], "module_pfc_excit_channels")
  save_gene_list(module[cls == "M2"], "module_pfc_excit_glutamate")
  save_gene_list(module[cls == "M3"], "module_pfc_excit_adhesion")
}
