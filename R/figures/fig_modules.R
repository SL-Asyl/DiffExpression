# =============================================================================
#  fig_modules.R — рисунки модулей (fig_mod_hpc_adhesion / _ion /
#  pfc_energy / pfc_excit). Читают списки генов из команд mod_*.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[ensure_pkgs])
box::use(R/core/engine[load_med, load_gene_list, heatmap_module, draw_narrow_panel])

#' @export
hpc_adhesion <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2")); suppressMessages({ library(readxl); library(openxlsx); library(ggplot2) })
  genes <- load_gene_list("module_hpc_adhesion"); medH <- load_med("HPC"); medP <- load_med("PFC")
  cat("[fig_mod_hpc_adhesion]", length(genes), "генов\n")
  heatmap_module(genes, medH, "Модуль адгезии и наведения аксонов (HPC)", "Fig_module_adhesion_HPC.png", PAIR[1])
  heatmap_module(genes, medP, "Тот же модуль в PFC", "Fig_module_adhesion_PFC.png", PAIR[1])
}

#' @export
hpc_ion <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2")); suppressMessages({ library(readxl); library(openxlsx); library(ggplot2) })
  genes <- load_gene_list("module_hpc_ion"); medH <- load_med("HPC"); medP <- load_med("PFC")
  cat("[fig_mod_hpc_ion]", length(genes), "генов\n")
  heatmap_module(genes, medH, "Модуль ионных каналов / Ca / глутаматных рецепторов (HPC)", "Fig_module_ion_HPC.png", PAIR[1])
  heatmap_module(genes, medP, "Тот же модуль в PFC", "Fig_module_ion_PFC.png", PAIR[1])
}

#' @export
pfc_energy <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2")); suppressMessages({ library(readxl); library(openxlsx); library(ggplot2) })
  genes <- load_gene_list("module_pfc_energy"); medP <- load_med("PFC")
  cat("[fig_mod_pfc_energy]", length(genes), "генов\n")
  heatmap_module(genes, medP, "Энергетический модуль PFC", "Fig_module_energy_PFC.png", PAIR[1])
}

#' @export
pfc_excit <- function() {
  ensure_pkgs(c("readxl", "openxlsx", "ggplot2")); suppressMessages({ library(readxl); library(openxlsx); library(ggplot2) })
  medP <- load_med("PFC"); rc <- PAIR[2]
  m1 <- load_gene_list("module_pfc_excit_channels")
  m2 <- load_gene_list("module_pfc_excit_glutamate")
  m3 <- load_gene_list("module_pfc_excit_adhesion")
  cat("[fig_mod_pfc_excit] каналы:", length(m1), "глутамат:", length(m2), "адгезия:", length(m3), "\n")
  draw_narrow_panel(m1[order(-abs(medP[[rc]][m1]))], medP, file.path(BASE_DIR, "Fig_submodule_PFC_a.png"),
                    sprintf("Модуль ионных каналов и возбудимости PFC (%d генов)", length(m1)))
  draw_narrow_panel(m2[order(-abs(medP[[rc]][m2]))], medP, file.path(BASE_DIR, "Fig_submodule_PFC_b.png"),
                    sprintf("Модуль глутаматных рецепторов и PSD PFC (%d генов)", length(m2)))
  ord3 <- m3[order(-abs(medP[[rc]][m3]))]; n3 <- length(ord3); k <- min(MOD_GENES_PER_PANEL, n3)
  draw_narrow_panel(ord3[seq_len(k)], medP, file.path(BASE_DIR, "Fig_submodule_PFC_c.png"),
                    sprintf("Модуль синаптической адгезии и сигнализации PFC (топ-%d из %d)", k, n3))
  if (n3 > k)
    draw_narrow_panel(ord3[(k + 1):n3], medP, file.path(BASE_DIR, "Fig_submodule_PFC_c2.png"),
                      sprintf("Модуль синаптической адгезии и сигнализации PFC (остальные %d из %d)", n3 - k, n3))
}
