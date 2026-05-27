# ==============================================================================
# КОМПЛЕКСНЫЙ КОНВЕЙЕР СБОРКИ МАСТЕР-ТАБЛИЦ
# Модуль: Автоматический парсинг KEGG, GO (BP, CC, MF) и Mammalian Phenotype (MP)
# ==============================================================================

# 1. Автоматическая проверка, установка и загрузка необходимых пакетов
required_packages <- c("openxlsx", "dplyr", "tidyr", "ggplot2", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# ---------- НАСТРОЙКИ ПОЛЬЗОВАТЕЛЯ --------------------------------------------
base_dir     <- "C:/Users/ES/Desktop/ФМБА/PTSD/Criteria_ORA/5 vs 6 EVS"
ora_dir      <- file.path(base_dir, "ORA_results") # Папка для сохранения результатов ORA
p_adj_cutoff <- 0.05                               # Порог значимости (FDR)

# Создаем подпапку результатов, если её ещё нет
dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)

# Находим все файлы результатов в исходной папке
ora_files <- list.files(path = base_dir, pattern = "^ORA_KEGG_GO_.*\\.xlsx$", full.names = TRUE)

if (length(ora_files) == 0) {
  stop("Критическая ошибка: Файлы ORA_KEGG_GO_*.xlsx не найдены в папке: ", base_dir)
}

cat("=== ЗАПУСК ПАЙПЛАЙНА СБОРКИ ДАННЫХ ===\n")
cat("Найдено файлов для анализа: ", length(ora_files), "\n\n")

# ==============================================================================
# БЛОК 1. СБОР ДАННЫХ KEGG И ПОСТРОЕНИЕ МАСШТАБНОЙ ТАБЛИЦЫ KEGG
# ==============================================================================
cat("--- [ЭТАП 1/2] Парсинг путей KEGG ---\n")
master_kegg_list <- list()

for (file_path in ora_files) {
  cell_type <- gsub("ORA_KEGG_GO_", "", basename(file_path))
  cell_type <- gsub("\\.xlsx$", "", cell_type)
  
  wb <- loadWorkbook(file_path)
  sheets <- names(wb)
  
  for (direction in c("UP", "DOWN")) {
    sheet_name <- paste0("KEGG_", direction)
    
    if (sheet_name %in% sheets) {
      df <- read.xlsx(file_path, sheet = sheet_name)
      
      if (!is.null(df) && nrow(df) > 0) {
        df_collapsed <- df %>%
          group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
          summarise(Genes = paste(unique(Genes), collapse = ", "), .groups = 'drop') %>%
          mutate(
            Cell_Type = cell_type,
            Direction = direction,
            GeneRatio_num = as.numeric(gsub("%", "", GeneRatio)) / 100
          )
        
        master_kegg_list[[length(master_kegg_list) + 1]] <- df_collapsed
      }
    }
  }
}

if (length(master_kegg_list) > 0) {
  master_kegg <- bind_rows(master_kegg_list) %>%
    mutate(Cell_Type_Clean = gsub("5_vs_6_", "", Cell_Type))
  
  # --- ДОБАВЛЕНИЕ АВТОФИЛЬТРОВ ДЛЯ KEGG ---
  wb_kegg <- createWorkbook()
  addWorksheet(wb_kegg, "KEGG_Master")
  writeData(wb_kegg, "KEGG_Master", master_kegg)
  
  # Включаем фильтры: указываем лист, строки и колонки (от 1 до последней)
  addFilter(wb_kegg, "KEGG_Master", row = 1, cols = 1:ncol(master_kegg))
  
  saveWorkbook(wb_kegg, file = file.path(ora_dir, "MASTER_KEGG_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Успех: Мастер-таблица KEGG сохранена с АВТОФИЛЬТРАМИ.\n\n")
} else {
  cat("[Предупреждение]: Значимые данные KEGG не обнаружены.\n\n")
}

# ==============================================================================
# БЛОК 2. СБОР ДАННЫХ GO И MP ИЗ СЛОЖНЫХ ДВУХБЛОЧНЫХ ЛИСТОВ
# ==============================================================================
cat("--- [ЭТАП 2/2] Парсинг онтологий GO и фенотипов MP ---\n")
master_go_list <- list()

for (file_path in ora_files) {
  cell_type <- gsub("ORA_KEGG_GO_", "", basename(file_path))
  cell_type <- gsub("\\.xlsx$", "", cell_type)
  
  wb <- loadWorkbook(file_path)
  sheets <- names(wb)
  
  go_sheets <- sheets[str_detect(sheets, "GO|MP")]
  
  for (sheet_name in go_sheets) {
    direction <- case_when(
      str_detect(sheet_name, "UP") ~ "UP",
      str_detect(sheet_name, "DOWN") ~ "DOWN",
      TRUE ~ "Unknown"
    )
    
    raw_sheet <- read.xlsx(file_path, sheet = sheet_name, colNames = FALSE, skipEmptyRows = FALSE)
    if (is.null(raw_sheet) || nrow(raw_sheet) == 0) next
    
    go_start_idx <- which(apply(raw_sheet, 1, function(x) any(str_detect(na.omit(as.character(x)), "(?i)GO enrichment"))))
    mp_start_idx <- which(apply(raw_sheet, 1, function(x) any(str_detect(na.omit(as.character(x)), "(?i)Mammalian Phenotype"))))
    
    # --- СУББЛОК 2А. ОБРАБОТКА ТАБЛИЦЫ GO ENRICHMENT ---
    if (length(go_start_idx) > 0) {
      go_header_row <- go_start_idx + 1
      go_end_row <- if(length(mp_start_idx) > 0) mp_start_idx - 1 else nrow(raw_sheet)
      
      if (go_end_row > go_header_row) {
        go_df <- read.xlsx(file_path, sheet = sheet_name, startRow = go_header_row + 1, 
                           rows = (go_header_row + 1):go_end_row, colNames = FALSE)
        
        if (!is.null(go_df) && nrow(go_df) > 0 && ncol(go_df) >= 13) {
          colnames(go_df)[1:13] <- c("Genes", "ONTOLOGY", "ID", "Description", "GeneRatio", 
                                     "BgRatio", "RichFactor", "FoldEnrichment", "zScore", 
                                     "pvalue", "p.adjust", "qvalue", "Count")
          
          go_df <- go_df %>%
            mutate(across(c(pvalue, p.adjust, qvalue), ~as.numeric(gsub(",", ".", as.character(.)))))
          
          go_collapsed <- go_df %>%
            filter(!is.na(ID)) %>%
            group_by(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
            summarise(Genes = paste(unique(Genes), collapse = ", "), .groups = 'drop') %>%
            mutate(
              Cell_Type = cell_type,
              Direction = direction,
              GeneRatio_num = as.numeric(gsub("%", "", GeneRatio)) / 100
            ) %>%
            dplyr::select(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, Genes, Cell_Type, Direction, GeneRatio_num)
          
          master_go_list[[length(master_go_list) + 1]] <- go_collapsed
        }
      }
    }
    
    # --- СУББЛОК 2Б. ОБРАБОТКА ТАБЛИЦЫ MAMMALIAN PHENOTYPE (MP) ---
    if (length(mp_start_idx) > 0) {
      mp_header_row <- mp_start_idx + 1
      mp_df <- read.xlsx(file_path, sheet = sheet_name, startRow = mp_header_row + 1, colNames = FALSE)
      
      if (!is.null(mp_df) && nrow(mp_df) > 0 && ncol(mp_df) >= 5) {
        colnames(mp_df)[1:5] <- c("Genes", "Raw_Description", "Overlap", "pvalue", "p.adjust")
        
        mp_df <- mp_df %>%
          filter(!is.na(Raw_Description)) %>%
          mutate(
            ID = str_extract(Raw_Description, "MP:\\d+"),
            Description = str_trim(gsub("MP:\\d+", "", Raw_Description)),
            pvalue = as.numeric(gsub(",", ".", as.character(pvalue))),
            p.adjust = as.numeric(gsub(",", ".", as.character(p.adjust)))
          )
        
        mp_collapsed <- mp_df %>%
          filter(!is.na(ID)) %>%
          group_by(ID, Description, pvalue, p.adjust) %>%
          summarise(
            Genes = paste(unique(Genes), collapse = ", "),
            Count = n_distinct(Genes),
            .groups = 'drop'
          ) %>%
          mutate(
            ONTOLOGY = "MP",
            GeneRatio = NA_character_,
            BgRatio = NA_character_,
            qvalue = NA_real_,
            GeneRatio_num = NA_real_,
            Cell_Type = cell_type,
            Direction = direction
          ) %>%
          dplyr::select(ID, Description, ONTOLOGY, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, Genes, Cell_Type, Direction, GeneRatio_num)
        
        master_go_list[[length(master_go_list) + 1]] <- mp_collapsed
      }
    }
    
  }
}

# Финальное объединение таблиц GO и MP с генерацией автофильтра
if (length(master_go_list) > 0) {
  master_go <- bind_rows(master_go_list) %>%
    mutate(Cell_Type_Clean = gsub("5_vs_6_", "", Cell_Type))
  
  # --- ДОБАВЛЕНИЕ АВТОФИЛЬТРОВ ДЛЯ GO/MP ---
  wb_go <- createWorkbook()
  addWorksheet(wb_go, "GO_MP_Master")
  writeData(wb_go, "GO_MP_Master", master_go)
  
  # Включаем фильтры на первую строку для всех колонок таблицы
  addFilter(wb_go, "GO_MP_Master", row = 1, cols = 1:ncol(master_go))
  
  saveWorkbook(wb_go, file = file.path(ora_dir, "MASTER_GO_ALL_CELL_TYPES.xlsx"), overwrite = TRUE)
  cat("-> Успех: Мастер-таблица GO/MP сохранена с АВТОФИЛЬТРАМИ.\n\n")
} else {
  stop("Критическая ошибка: Данные GO/MP не были собраны.")
}

cat("=== СБОРКА ОБОИХ МОДУЛЕЙ ЗАВЕРШЕНА УСПЕШНО ===\n")

