# =============================================================================
#  config.R — ЕДИНЫЙ файл настроек. ПРАВИТСЯ ВРУЧНУЮ перед запуском.
#  Это box-модуль: все переменные ниже видны остальным модулям через
#  box::use(R/config). Больше нигде пути задавать не нужно.
# =============================================================================

# --- корень с папками-выходами (EVS / HPC_*) ---
# TODO: укажите путь к корневой папке с выходными данными пайплайна.
BASE_DIR <- ""   # напр. "D:/DiffExpression_output"

# --- активный датасет (для по-датасетных команд) ---
TISSUE   <- "PFC"          # "HPC" | "PFC"
CONTRAST <- "5v6"          # "5v6" | "5v9" | "9v10" | "6v10"

# --- активная пара (для парных команд: модули / graphpad / common_genes) ---
#   PAIR[1] = FS у C57Bl/6N (обычно 5v6); PAIR[2] = FS у k18-hACE2 (9v10)
PAIR <- c("5v6", "9v10")

# --- пороги значимости (degs / master) ---
LOGFC_CUTOFF <- 1
FDR_CUTOFF   <- 0.05

# --- соответствие кодов контрастов именам папок/файлов ---
TAG <- c('5v6' = '5 vs 6', '5v9' = '5 vs 9', '9v10' = '9 vs 10', '6v10' = '6 vs 10')
HD  <- c('5v6' = 'HPC_5 vs 6',   '5v9' = 'HPC_5 vs 9',
         '9v10' = 'HPC_9 vs 10', '6v10' = 'HPC_6 vs 10')
PD  <- c('5v6' = 'PFC 5 vs 6 EVS',  '5v9' = 'PFC 5 vs 9 EVS',
         '9v10' = 'PFC 9 vs 10 EVS', '6v10' = 'PFC 6 vs 10 EVS')

# --- датасеты для meta-команд (цикл) ---
CONTRASTS4 <- c('5v6', '5v9', '9v10', '6v10')
DATASETS   <- c(paste0("HPC_", CONTRASTS4), paste0("PFC_", CONTRASTS4))

# --- ПАПКИ С СЫРЫМИ ТАБЛИЦАМИ ДЭ ПО ТИПАМ КЛЕТОК (вход для degs / master / gsea) ---
#   Ключ = "<TISSUE>_<contrast>"; значение = путь к папке с .xlsx по типам клеток
#   (колонки names / logfoldchanges / pvals_adj). Заполнить вручную ИЛИ оставить
#   RAW_CHOOSE_INTERACTIVE <- TRUE, чтобы выбирать папку окном при запуске.
RAW_CHOOSE_INTERACTIVE <- FALSE
RAW_DIRS <- c(
  # укажите пути к папкам сырых таблиц, напр.:
  # "PFC_5v6"  = "path/to/raw/PFC 5 vs 6",
  # "HPC_5v6"  = "path/to/raw/HPC 5 vs 6",
  # "PFC_9v10" = "path/to/raw/PFC 9 vs 10",
  # "HPC_9v10" = "path/to/raw/HPC 9 vs 10"
)

# --- параметры master (ORA) ---
MP_DB            <- "MGI_Mammalian_Phenotype_Level_4_2021"
N_WORKERS        <- 2
RETRY_MAX_TRIES  <- 3
RETRY_WAIT_SEC   <- 5
JITTER_RANGE_SEC <- c(0.5, 2)

# --- параметры rrvgo ---
SIM_THRESHOLD <- 0.7
SIM_METHOD    <- "Rel"
ORGANISM_DB   <- "org.Mm.eg.db"

# --- параметры fig_bubble ---
BUBBLE_TOP_N    <- 10
BUBBLE_SIZE     <- 9
BUBBLE_PT_RANGE <- c(4, 12)
BUBBLE_WRAP     <- 34
BUBBLE_FONT     <- 15
BUBBLE_COL_LOW  <- "#FF0000"
BUBBLE_COL_HIGH <- "#0000FF"

# --- параметры fig_heatmap_deg ---
HM_WIDTH <- 6; HM_HEIGHT <- 12; HM_SCALE <- "row"
HM_FONT_ROW <- 8; HM_FONT_COL <- 8
HM_CLUSTER_METHOD <- "ward.D"; HM_DISTANCE <- "euclidean"
HM_COL_LOW <- "#0000FF"; HM_COL_MID <- "#FFFFFF"; HM_COL_HIGH <- "#FF0000"
HM_TOP_N <- 50; HM_MIN_CT <- 2

# --- параметры fig_volcano ---
VOLC_W <- 8; VOLC_H <- 8
VOLC_PMIN <- 1e-100; VOLC_FDR <- 0.05; VOLC_FC <- 2
VOLC_PT_SIZE <- 80; VOLC_N_LABEL <- 20
VOLC_COL_UP <- "#F39B7F"; VOLC_COL_DOWN <- "#4DBBD5"; VOLC_COL_NS <- "#A6A6A6"

# --- параметры автоматических функциональных модулей (mod_select, rrvgo) ---
MOD_K_HPC     <- 2       # число семантических кластеров для HPC (reversal-набор)
MOD_K_PFC_C57 <- 2       # PFC, ответ FS у C57Bl/6N (PAIR[1])
MOD_K_PFC_K18 <- 5       # PFC, ответ FS у k18-hACE2 (PAIR[2])
MOD_TOP_N     <- 35      # максимум генов на тепловую карту (0 = все)
MOD_CAP       <- 3.2     # предел |log2FC| для цветовой шкалы
MOD_FIG_W     <- 5.2     # ширина карты, дюймы
MOD_ROW_H     <- 0.185   # высота строки (гена), дюймы
MOD_COL_LOW   <- "#2A78D6"; MOD_COL_MID <- "#F0EFEC"; MOD_COL_HIGH <- "#D42A2A"
# три контраста тепловой карты и их подписи (порядок столбцов согласован):
#   базовый генотип (5v9) первым, затем FS у C57Bl/6N и у k18-hACE2
MOD_CONTRASTS <- c("5v9", "5v6", "9v10")
MOD_COL_ORDER <- c("no FS:\nC57Bl/6N vs k18-hACE2",
                   "C57Bl/6N:\nno FS vs FS",
                   "k18-hACE2:\nno FS vs FS")

# --- параметры common_genes ---
CG_SHEET       <- "Significant_DEGs"
CG_SCOPE       <- "all"            # "all" | "neuronal"
CG_NEURONAL_RX <- "_(Glut|Gaba)$"
CG_MIN_PAIRS   <- 1
