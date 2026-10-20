# =============================================================================
#  meta.R — составные команды: prep / figures / modules / all.
# =============================================================================

box::use(R/config)
for (.n in ls(config)) assign(.n, get(.n, envir = config), environment())
box::use(R/core/helpers[parse_dataset])
box::use(
  R/data/degs, R/data/master, R/data/rrvgo,
  R/figures/fig_bubble, R/figures/fig_volcano, R/figures/fig_heatmap_deg,
  R/analysis/mod_select, R/figures/fig_modules,
  R/analysis/graphpad, R/analysis/common_genes
)

#' @export
prep <- function() {
  for (key in DATASETS) {
    d <- parse_dataset(key); cat("\n===== prep:", key, "=====\n")
    degs$run(d$tissue, d$cc); master$run(d$tissue, d$cc); rrvgo$run(d$tissue, d$cc)
  }
}

#' @export
figures <- function() {
  for (key in DATASETS) {
    d <- parse_dataset(key); cat("\n===== figures:", key, "=====\n")
    fig_bubble$run(d$tissue, d$cc); fig_volcano$run(d$tissue, d$cc); fig_heatmap_deg$run(d$tissue, d$cc)
  }
  fig_modules$hpc_adhesion(); fig_modules$hpc_ion(); fig_modules$pfc_energy(); fig_modules$pfc_excit()
}

#' @export
modules <- function() {
  mod_select$hpc(); mod_select$pfc_energy(); mod_select$pfc_excit()
  fig_modules$hpc_adhesion(); fig_modules$hpc_ion(); fig_modules$pfc_energy(); fig_modules$pfc_excit()
}

#' @export
all <- function() {
  prep(); graphpad$run(); common_genes$run(); modules(); figures()
}
