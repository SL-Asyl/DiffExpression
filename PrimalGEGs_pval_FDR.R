# ==============================================================================
# СКРИПТ ВЫГРУЗКИ СВОДНОЙ ТАБЛИЦЫ ДИФФЕРЕНЦИАЛЬНО ЭКСПРЕССИРОВАННЫХ ГЕНОВ (DEGs)
# Модуль: Сбор значимых генов по индивидуальным вкладкам для каждого типа клеток
# ==============================================================================

# Автоматическая проверка и установка необходимых пакетов
required_packages <- c("openxlsx", "dplyr", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(openxlsx)
library(dplyr)
library(stringr)

# ---------- НАСТРОЙКИ ПОЛЬЗОВАТЕЛЯ (USER INPUT) -------------------------------
input_dir    <- "C:/Users/ES/Desktop/ФМБА/PTSD/PTSD_ PFC difexp_Egorova/PFC 9 vs 10_ Egorova"     # Папка с исходными Excel файлами от секвенирования
outdir       <- "C:/Users/ES/Desktop/ФМБА/PTSD/Criteria_ORA/9 vs 10 EVS/ORA_results"  # Папка для сохранения итогового сводного файла

logFC_cutoff <- 1              # Порог изменения экспрессии (в 2 раза)
FDR_cutoff   <- 0.05           # Порог значимости (FDR)

# Создаем папку для результатов, если её ещё нет
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Находим все исходные файлы
raw_files <- list.files(path = input_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(raw_files) == 0) {
  stop("Критическая ошибка: не найдено файлов .xlsx в папке '", input_dir, "'. Проверьте путь!")
}

cat("=== ЗАПУСК ЭКСПОРТА СПИСКОВ ДИФФЕРЕНЦИАЛЬНЫХ ГЕНОВ ===\n")
cat("Найдено файлов для обработки: ", length(raw_files), "\n\n")

# Создаем единый воркбук Excel, куда будем складывать листы
wb_all_degs <- createWorkbook()

# ==============================================================================
# ОСНОВНОЙ ЦИКЛ ОБРАБОТКИ ФАЙЛОВ
# ==============================================================================
for (file_path in raw_files) {
  
  # Очищаем имя файла, чтобы получить чистое название типа клеток
  cell_type_name <- gsub(".xlsx", "", basename(file_path))
  cell_type_name <- gsub("5_vs_6_", "", cell_type_name) # убираем префикс сравнения, если есть
  
  # Excel ограничивает длину имени вкладки в 31 символ. Обрезаем, если имя слишком длинное
  sheet_name <- substr(cell_type_name, 1, 31)
  
  cat("Обработка кластера:", cell_type_name, "... ")
  
  # Чтение сырых данных с защитой от ошибок
  tb_raw <- tryCatch({
    read.xlsx(file_path)
  }, error = function(e) {
    cat("[ОШИБКА ЧТЕНИЯ] Пропуск.\n")
    next
  })
  
  # Стандартизация колонок и базовая очистка от NA
  tb <- tb_raw %>%
    rename(gene = names, log2FoldChange = logfoldchanges, padj = pvals_adj) %>%
    select(gene, log2FoldChange, padj) %>%
    filter(!is.na(gene), !is.na(log2FoldChange), !is.na(padj))
  
  if (nrow(tb) == 0) {
    cat("[ТАБЛИЦА ПУСТАЯ после фильтрации NA] Пропуск.\n")
    next
  }
  
  # Фильтрация генов по порогам FDR и Log2FC
  filtered_degs <- tb %>%
    filter(padj < FDR_cutoff, abs(log2FoldChange) >= logFC_cutoff) %>%
    # Добавляем текстовый маркер направления для удобства фильтрации в Excel
    mutate(Direction = ifelse(log2FoldChange >= logFC_cutoff, "UP", "DOWN")) %>%
    # Сортируем гены по силе изменения (сначала самые мощные UP, затем самые мощные DOWN)
    arrange(desc(log2FoldChange)) %>%
    # Красиво переупорядочиваем колонки
    select(gene, Direction, log2FoldChange, padj)
  
  # Если значимых генов в этом типе клеток нет, вкладку не создаем
  if (nrow(filtered_degs) == 0) {
    cat("0 значимых генов. Пропуск.\n")
    next
  }
  
  # Запись данных на индивидуальный лист
  addWorksheet(wb_all_degs, sheetName = sheet_name)
  writeData(wb_all_degs, sheet = sheet_name, x = filtered_degs)
  
  # Программное включение кнопок-автофильтров на шапку таблицы
  addFilter(wb_all_degs, sheet = sheet_name, row = 1, cols = 1:ncol(filtered_degs))
  
  # Делаем ширину колонок автоматической по размеру контента, чтобы текст не обрезался
  setColWidths(wb_all_degs, sheet = sheet_name, cols = 1:ncol(filtered_degs), widths = "auto")
  
  cat("Успешно добавлено генов:", nrow(filtered_degs), "\n")
}

# ==============================================================================
# СОХРАНЕНИЕ ИТОГОВОГО СВОДНОГО ФАЙЛА
# ==============================================================================
output_file <- file.path(outdir, "ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx")

tryCatch({
  saveWorkbook(wb_all_degs, output_file, overwrite = TRUE)
  cat("\n============================================================\n")
  cat("СБОРКА ЗАВЕРШЕНА! Сводный файл сохранен в:\n")
  cat(output_file, "\n")
  cat("============================================================\n")
}, error = function(e) {
  message("\nКРИТИЧЕСКАЯ ОШИБКА: Не удалось сохранить итоговый файл. Закройте его, если он открыт в Excel!")
})

