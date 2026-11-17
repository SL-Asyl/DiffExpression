#!/usr/bin/env Rscript
# =============================================================================
#  run.R — ТОЧКА ВХОДА. Запускает одну команду конвейера.
#
#  КАК ЗАПУСКАТЬ
#  ------------
#  (1) Командная строка:   Rscript run.R <команда>        напр.  Rscript run.R degs
#  (2) В RStudio:          задать SCENARIO <- "degs" ниже и нажать Source.
#
#  Перед запуском проверить настройки в R/config.R (пути, TISSUE, CONTRAST, PAIR,
#  RAW_DIRS). Список всех команд — в R/config.R нет; см. README.md или switch ниже.
#
#  ВАЖНО (RStudio): box кэширует модули. После правки config.R или модулей
#  выполнить box::reload(...) или перезапустить R, иначе подхватится старое.
#  При запуске через Rscript каждый вызов — свежий процесс, кэш не мешает.
# =============================================================================

if (!requireNamespace("box", quietly = TRUE)) install.packages("box")

# --- определить папку скрипта (корень проекта) ---
.get_script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
.root <- .get_script_dir()
setwd(.root)
# корень проекта на box.path -> модули адресуются как R/<папка>/<модуль>
options(box.path = unique(c(.root, getOption("box.path"))))

box::use(
  R/config,
  R/data/degs, R/data/master, R/data/rrvgo, R/data/gsea,
  R/figures/fig_bubble, R/figures/fig_volcano, R/figures/fig_heatmap_deg,
  R/analysis/graphpad, R/analysis/common_genes,
  R/analysis/mod_select, R/meta
)

dispatch <- function(scenario) {
  tis <- config$TISSUE; cc <- config$CONTRAST
  switch(scenario,
         degs                 = degs$run(tis, cc),
         master               = master$run(tis, cc),
         rrvgo                = rrvgo$run(tis, cc),
         gsea                 = gsea$run(tis, cc),
         fig_bubble           = fig_bubble$run(tis, cc),
         fig_volcano          = fig_volcano$run(tis, cc),
         fig_heatmap_deg      = fig_heatmap_deg$run(tis, cc),
         graphpad             = graphpad$run(),
         common_genes         = common_genes$run(),
         mod_hpc              = mod_select$hpc(),
         mod_pfc              = mod_select$pfc(),
         prep                 = meta$prep(),
         figures              = meta$figures(),
         modules              = meta$modules(),
         all                  = meta$all(),
         stop("Неизвестная команда: '", scenario, "'. См. README.md.")
  )
  cat("\n=== команда '", scenario, "' выполнена ===\n", sep = "")
}

# --- выбор команды: аргумент командной строки имеет приоритет над SCENARIO ---
SCENARIO <- "degs"   # <-- задать здесь при запуске внутри RStudio

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && nzchar(.args[1])) SCENARIO <- .args[1]

dispatch(SCENARIO)
