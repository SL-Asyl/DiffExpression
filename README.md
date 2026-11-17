------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# PTSD_pipeline — подробное руководство

Модульный командный конвейер анализа дифференциальной экспрессии (k18-hACE2 / C57Bl/6N, отделы PFC и HPC). Собран на пакете **`box`**: каждая команда — отдельный модуль; запускается одной командой через `run.R`.

------------------------------------------------------------------------

## 1. Содержание

1.  [Быстрый старт](#2-быстрый-старт)
2.  [Требования](#3-требования)
3.  [Соглашения о путях и именах файлов](#4-соглашения-о-путях-и-именах-файлов)
4.  [Настройки (config.R)](#5-настройки-configr)
5.  [Команды — подробно](#6-команды--подробно)
6.  [Схема зависимостей](#7-схема-зависимостей-что-от-чего-зависит)
7.  [Рекомендуемый порядок прогона](#8-рекомендуемый-порядок-прогона-с-нуля)
8.  [Как запускать](#9-как-запускать)
9.  [Частые проблемы](#10-частые-проблемы)

------------------------------------------------------------------------

## 2. Быстрый старт

``` r
install.packages("box")          # один раз
```

1.  Открыть `R/config.R`, указать `BASE_DIR` и заполнить `RAW_DIRS` (папки с сырыми таблицами ДЭ по типам клеток).

2.  Запустить по стадиям (см. §8). Коротко — с нуля:

    ```         
    Rscript run.R prep       # degs + master + rrvgo по всем датасетам
    Rscript run.R modules    # отбор генов модулей + их рисунки
    Rscript run.R figures    # bubble / volcano / heatmap + модульные рисунки
    Rscript run.R graphpad
    Rscript run.R common_genes
    ```

    либо почти всё одной командой: `Rscript run.R all` (без `gsea`).

------------------------------------------------------------------------

## 3. Требования

- **R** (рекомендуется ≥ 4.1) и пакет **`box`**.
- Пакеты ставятся автоматически при первом запуске нужной команды (функция `ensure_pkgs`): CRAN (`openxlsx`, `readxl`, `dplyr`, `tidyr`, `ggplot2`, `ggrepel`, `pheatmap`, `future`, `future.apply`, `progressr`, `enrichR`, `RColorBrewer`) и Bioconductor (`clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `rrvgo`, `GO.db`).
- **Интернет нужен только двум командам:** `master` (запросы KEGG + Enrichr MP) и `gsea` (KEGG/GO через clusterProfiler). Остальные работают офлайн.

| Команда          | Интернет |
|------------------|:--------:|
| `master`, `gsea` |  **да**  |
| все остальные    |   нет    |

------------------------------------------------------------------------

## 4. Соглашения о путях и именах файлов

Всё строится от `BASE_DIR` (задаётся в `config.R`). Код контраста (`cc`) — один из `5v6`, `5v9`, `9v10`, `6v10`; человекочитаемый тег берётся из `TAG` (`5v6` → `"5 vs 6"`).

**Папки-выходы (EVS) по датасету `TISSUE + cc`:**

| Отдел | Папка                                             |
|-------|---------------------------------------------------|
| HPC   | `BASE_DIR/HPC_<tag>` (напр. `HPC_5 vs 6`)         |
| PFC   | `BASE_DIR/PFC <tag> EVS` (напр. `PFC 5 vs 6 EVS`) |

**Файлы в папке-выходе** (префикс = `<tag>`):

| Файл | Создаёт | Потребляют |
|------------------------|------------------------|------------------------|
| `<tag>_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx` | `degs` | fig_volcano, fig_heatmap_deg, mod\_*, graphpad, common_genes |
| `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx` | `master` | rrvgo, fig_bubble, mod\_\* |
| `<tag>_MASTER_KEGG_ALL_CELL_TYPES.xlsx` | `master` | fig_bubble |
| `<tag>_MASTER_GO_FILTERED.xlsx` | `rrvgo` | fig_bubble |
| `<tag>_MASTER_GSEA.xlsx` | `gsea` | — (самостоятельный результат) |

**Сырые таблицы ДЭ по типам клеток (вход для `degs`, `master`, `gsea`)** задаются вручную в `RAW_DIRS` по ключу `"<TISSUE>_<cc>"`, напр. `RAW_DIRS["PFC_5v6"]`. Это папка с .xlsx по типам клеток, в каждом — колонки `names`, `logfoldchanges`, `pvals_adj`.

**Файлы в корне `BASE_DIR`** (общие результаты модулей и GraphPad): `module_*_genes.txt`, `Fig_mod_*.png`, `common_genes_wide_log2FC.xlsx`, папка `Prism_interaction/`.

------------------------------------------------------------------------

## 5. Настройки (config.R)

| Параметр | Назначение |
|------------------------------------|------------------------------------|
| `BASE_DIR` | корень проекта (все EVS-папки внутри) |
| `TISSUE`, `CONTRAST` | активный датасет для по-датасетных команд |
| `PAIR` | пара контрастов для парных команд: `PAIR[1]`=FS у C57 (5v6), `PAIR[2]`=FS у k18 (9v10) |
| `LOGFC_CUTOFF`, `FDR_CUTOFF` | пороги значимости (degs, master) |
| `TAG`, `HD`, `PD` | карты кодов контрастов → теги и имена папок |
| `CONTRASTS4`, `DATASETS` | список контрастов и датасетов для meta-команд |
| `RAW_DIRS`, `RAW_CHOOSE_INTERACTIVE` | папки сырых таблиц ДЭ (ручной список / выбор окном) |
| `MP_DB`, `N_WORKERS`, `RETRY_*`, `JITTER_*` | параметры ORA (`master`) |
| `SIM_THRESHOLD`, `SIM_METHOD`, `ORGANISM_DB` | параметры свёртки GO (`rrvgo`) |
| `BUBBLE_*` | параметры bubble-plot |
| `HM_*` | параметры heatmap DEG |
| `VOLC_*` | параметры volcano |
| `MOD_*` | параметры авто-модулей: число кластеров (`MOD_K_*`), топ-N генов, контрасты и стиль карт |
| `CG_*` | параметры common_genes (лист, scope all/neuronal) |

**Гибридная модель контраста:** по-датасетные команды смотрят на `TISSUE`+`CONTRAST`; парные (модули, graphpad, common_genes) — на `PAIR`; meta-команды идут циклом по `DATASETS`.

------------------------------------------------------------------------

## 6. Команды — подробно

Обозначения: **Вход** — что должно уже лежать на диске; **Выход** — что создаётся; **Настройки** — ключевые параметры из config.R.

### Подготовка данных

#### `degs`

- **Что делает:** фильтрует сырые таблицы ДЭ до значимых DEG (`padj < FDR_CUTOFF` и `|log2FC| ≥ LOGFC_CUTOFF`), собирает в одну таблицу по всем типам клеток. Параллельно (future).
- **Вход:** `RAW_DIRS["<TISSUE>_<CONTRAST>"]` — папка сырых .xlsx (`names`, `logfoldchanges`, `pvals_adj`).
- **Выход:** `<tag>_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx` в папке-выходе.
- **Настройки:** `FDR_CUTOFF`, `LOGFC_CUTOFF`.

#### `master` *(интернет)*

- **Что делает:** ORA (over-representation) по GO (BP/MF/CC), KEGG и MP (Enrichr) отдельно для up- и down-генов каждого типа клеток; собирает мастер-таблицы.
- **Вход:** те же сырые таблицы ДЭ `RAW_DIRS["<TISSUE>_<CONTRAST>"]`.
- **Выход:** `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx` и `<tag>_MASTER_KEGG_ALL_CELL_TYPES.xlsx`.
- **Настройки:** `MP_DB`, `N_WORKERS`, `RETRY_*`, `JITTER_*`, пороги.

#### `rrvgo`

- **Что делает:** семантически сворачивает избыточные GO-термины (rrvgo, офлайн через GO.db), оставляя по одному представителю на кластер. KEGG и MP не трогает.
- **Вход:** `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx` (результат `master`).
- **Выход:** `<tag>_MASTER_GO_FILTERED.xlsx`.
- **Настройки:** `SIM_THRESHOLD`, `SIM_METHOD`, `ORGANISM_DB`.

#### `gsea` *(интернет)*

- **Что делает:** GSEA по ПОЛНЫМ таблицам (все гены, без порога); ранжирует по log2FC, ищет согласованно сдвинутые наборы (KEGG + GO BP), даёт направление (NES) и `core_enrichment` — короткий список генов пути.
- **Вход:** сырые таблицы ДЭ `RAW_DIRS["<TISSUE>_<CONTRAST>"]`.
- **Выход:** `<tag>_MASTER_GSEA.xlsx`.
- **Самостоятельна:** ни от чего в конвейере не зависит, кроме сырых данных.

### Рисунки

#### `fig_bubble`

- **Что делает:** bubble-plots обогащения из ГОТОВЫХ мастер-таблиц (ничего не пересчитывает). Рисует unfiltered (GO BP/MF/CC, MP, KEGG) и filtered (GO).
- **Вход:** `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx`, `<tag>_MASTER_GO_FILTERED.xlsx`, `<tag>_MASTER_KEGG_ALL_CELL_TYPES.xlsx` (каждый по наличию).
- **Выход:** PNG в `<папка-выхода>/EnrichmentPlots`.
- **Настройки:** `BUBBLE_*`.

#### `fig_volcano`

- **Что делает:** volcano-графики — общий и по группам клеток (Glutamatergic / GABAergic / Glia / Microglia).
- **Вход:** `<tag>_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx`.
- **Выход:** PNG в `<папка-выхода>/<tag> VolcanoPlots`.
- **Настройки:** `VOLC_*`.

#### `fig_heatmap_deg`

- **Что делает:** тепловые карты значимых DEG — общая и по группам клеток (топ-`HM_TOP_N` генов, встречающихся в ≥ `HM_MIN_CT` типах).
- **Вход:** `<tag>_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx`.
- **Выход:** PNG (+ .txt матрицы) в `<папка-выхода>/<tag> Heatmaps`.
- **Настройки:** `HM_*`.

### Данные под GraphPad (парные, используют `PAIR`)

#### `graphpad`

- **Что делает:** взаимодействие «FS × генотип» в HPC: биномиальный тест разнонаправленности общих пар «ген × тип клетки» и корреляция Спирмена log2FC между линиями (по парам и по медианам генов).
- **Вход:** `<tagA>_ALL_CELL_TYPES_SIGNIFICANT_DEGs.xlsx` (HPC, `PAIR[1]`) и `<tagB>_..._SIGNIFICANT_DEGs.xlsx` (HPC, `PAIR[2]`).
- **Выход:** в `BASE_DIR/Prism_interaction/`: `HPC_FS_binomial_discordance.csv`, `HPC_FS_log2FC_paired_pairs.csv`, `HPC_FS_log2FC_paired_genes.csv`, `README_HPC_binomial_spearman.txt`.
- **Контроль (ожидаемо):** 79 пар, 70 разнонаправлены, биномиальный p≈7,8e-13; ρ(пары)≈−0,64; ρ(гены)≈−0,39, n=153.

#### `common_genes`

- **Что делает:** общие значимые гены двух контрастов → широкая таблица log2FC (overall-медиана + по типам клеток) и корреляция Спирмена (overall + по типам).
- **Вход:** `<tagA>` и `<tagB>` HPC `..._SIGNIFICANT_DEGs.xlsx` (лист `CG_SHEET`).
- **Выход:** `BASE_DIR/common_genes_wide_log2FC.xlsx` (листы `common_genes`, `correlations`).
- **Настройки:** `CG_SCOPE` (`all` \| `neuronal`), `CG_MIN_PAIRS`.

### Функциональные модули (парные; авто-кластеризация rrvgo)

#### `mod_hpc`

- **Что делает:** автоматически выделяет функциональные модули HPC. Берёт reversal-набор (ген значим в обоих FS-сравнениях `PAIR` и разворачивает знак медианного log2FC между линиями), собирает GO BP-термины этих генов, считает семантическое сходство терминов (`rrvgo`), кластеризует их (`hclust` + `cutree`, `k = MOD_K_HPC`) и присваивает каждый ген кластеру по большинству его терминов. Имя кластера — по ведущим GO-терминам. Для каждого кластера сразу рисует тепловую карту и пишет список генов.
- **Вход:** DEG HPC для `PAIR[1]` и `PAIR[2]` (определение reversal) и для контрастов `MOD_CONTRASTS` (медианы для карт) **плюс** `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx` HPC для контрастов `MOD_CONTRASTS` (аннотация ген→GO).
- **Выход:** `Fig_mod_hpc_<name>.png` и `module_hpc_<name>_genes.txt` в `BASE_DIR` (число файлов = число кластеров).
- **Настройки:** `MOD_K_HPC`, `MOD_TOP_N`, стиль карт `MOD_*`, `SIM_METHOD`, `ORGANISM_DB`.

#### `mod_pfc`

- **Что делает:** то же для PFC, но без reversal (разнонаправленности между линиями в PFC нет): отдельно кластеризует значимые гены ответа FS у C57Bl/6N (`PAIR[1]`, `k = MOD_K_PFC_C57`) и у k18-hACE2 (`PAIR[2]`, `k = MOD_K_PFC_K18`).
- **Вход:** DEG PFC для `PAIR[1]` и `PAIR[2]` (значимые гены) и для `MOD_CONTRASTS` (медианы для карт) **плюс** `<tag>_MASTER_GO_ALL_CELL_TYPES.xlsx` PFC для `PAIR[1]` и `PAIR[2]`.
- **Выход:** `Fig_mod_pfc_c57_<name>.png` / `module_pfc_c57_<name>_genes.txt` и `Fig_mod_pfc_k18_<name>.png` / `module_pfc_k18_<name>_genes.txt`.
- **Настройки:** `MOD_K_PFC_C57`, `MOD_K_PFC_K18`, `MOD_TOP_N`, `MOD_*`, `SIM_METHOD`, `ORGANISM_DB`.

### Составные (meta)

| Команда | Что делает |
|------------------------------------|------------------------------------|
| `prep` | по каждому датасету из `DATASETS`: `degs` → `master` → `rrvgo` |
| `figures` | по каждому датасету: `fig_bubble` → `fig_volcano` → `fig_heatmap_deg` |
| `modules` | `mod_hpc` + `mod_pfc` (кластеризация + тепловые карты модулей) |
| `all` | `prep` → `graphpad` + `common_genes` → `modules` → `figures` (без `gsea`) |

------------------------------------------------------------------------

## 7. Схема зависимостей (что от чего зависит)

```         
RAW_DIRS (сырые таблицы ДЭ)
   │
   ├─ degs ─────────────► *_SIGNIFICANT_DEGs.xlsx ──────────────────────────┐
   │                                                                        │
   ├─ master (интернет) ─► *_MASTER_GO_ALL, *_MASTER_KEGG                    │
   │                            │                                           │
   │                            └─ rrvgo ─► *_MASTER_GO_FILTERED            │
   │                                                                        │
   └─ gsea (интернет) ───► *_MASTER_GSEA   (ни на что больше не влияет)      │
                                                                            │
   *_MASTER_GO_ALL / _FILTERED / _KEGG ──► fig_bubble                       │
                                                                            │
   *_SIGNIFICANT_DEGs ──► fig_volcano, fig_heatmap_deg                      │
                      ──► graphpad, common_genes (HPC, пара PAIR) ◄─────────┤
                      ──► mod_hpc / mod_pfc (+ нужны *_MASTER_GO_ALL)        │
                             └─► Fig_mod_*.png + module_*_genes.txt  ◄───────┘
```

Ключевые требования, которые легко упустить: - **`rrvgo` требует `master`** (читает `MASTER_GO_ALL`). - **`fig_bubble` требует `master`** (и `rrvgo` — для filtered-картинок). - **`mod_hpc` / `mod_pfc` требуют и `degs`, и `master`** (DEG + аннотация GO) для контрастов `MOD_CONTRASTS` соответствующего отдела; они сами и кластеризуют, и рисуют карты. Практически: перед командой `modules` должен быть сделан `degs` для контрастов `MOD_CONTRASTS` (проще всего — `prep`).

------------------------------------------------------------------------

## 8. Рекомендуемый порядок прогона (с нуля)

0.  **Настроить `R/config.R`:** `BASE_DIR`, `RAW_DIRS` (минимум для нужных датасетов), при необходимости `TISSUE`/`CONTRAST`/`PAIR`.
1.  **`prep`** — соберёт `degs` + `master` + `rrvgo` по всем `DATASETS`. *(нужен интернет из-за `master`. Если хочется по одному датасету — выставить `TISSUE`/`CONTRAST` и гонять `degs`, `master`, `rrvgo` по отдельности.)*
2.  **`modules`** — `mod_hpc` + `mod_pfc`: семантическая кластеризация, отбор генов и тепловые карты модулей.
3.  **`figures`** — `fig_bubble`, `fig_volcano`, `fig_heatmap_deg` по датасетам
    - модульные карты.
4.  **`graphpad`**, **`common_genes`** — данные для Prism и корреляций.
5.  **`gsea`** — по желанию, по нужным датасетам (интернет).

Почти всё сразу: **`all`** = шаги 1→4 (без `gsea`). `gsea` всегда отдельно.

> Для шагов 2–4, если не гонишь `prep`/`all`, убедись, что `degs` сделан для контрастов `MOD_CONTRASTS` отдела (модульные карты читают их), а для `fig_bubble` — что сделаны `master` и `rrvgo`.

------------------------------------------------------------------------

## 9. Как запускать

**Командная строка** (каждый запуск — свежий процесс R):

```         
Rscript run.R <команда>
```

`run.R` сам определяет корень проекта, ставит его на `options(box.path=...)` и находит модули по именам `R/<папка>/<модуль>`.

**RStudio:** открыть `run.R`, задать внизу `SCENARIO <- "<команда>"`, нажать Source. \> После правки `config.R` или любого модуля в той же сессии нужно \> `box::reload(...)` или перезапуск R — `box` кэширует модули. Через `Rscript` \> этого не требуется.

Полный список команд: `degs master rrvgo gsea fig_bubble fig_volcano fig_heatmap_deg graphpad common_genes mod_hpc mod_pfc prep figures modules all`.

------------------------------------------------------------------------

## 10. Частые проблемы

| Симптом | Причина / решение |
|------------------------------------|------------------------------------|
| `Не задана папка сырых данных для PFC_5v6` | не заполнен `RAW_DIRS["PFC_5v6"]` в config.R (или включи `RAW_CHOOSE_INTERACTIVE`) |
| `fig_bubble` пропускает filtered-картинки | не сделан `rrvgo` (нет `MASTER_GO_FILTERED`) |
| `master`/`gsea` падают с сетевой ошибкой | нужен интернет (KEGG/Enrichr); при обрывах помогает повтор — в `master` встроен retry |
| карты модулей с пустыми ячейками «–» | у гена нет значимого log2FC в этом контрасте — это норма, не ошибка |
| в RStudio не подхватились правки config.R | `box::reload(...)` или перезапуск R |
| `could not find function`/модуль не найден | запускай через `run.R` (он ставит `box.path`); отдельные модули сами по себе не грузятся |

------------------------------------------------------------------------
